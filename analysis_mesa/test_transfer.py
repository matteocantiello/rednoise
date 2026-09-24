#!/usr/bin/env python3
"""
Test simplified FeCZ -> surface transfer models (transfer_models.py) star by star.

For every candidate predictor X (velocity, amplitude, or frequency) the MESA main-sequence
models of one sub-grid (default MW, omega/omega_c = 0) are interpolated to the (log Teff,
log L_spec) position of each observed main-sequence star. Where the models have no FeCZ
(present in < half of the surrounding models) X = 0. Then

  velocity, amplitude:  log y = log10(10^f + 10^a X^beta)     (floor f, normalisation a, exponent beta)
     'free'  : f, a, beta fitted            'linear': beta = 1 (the model's own scaling)
  frequency:            log y = a + beta log X                 ('free' and 'linear': beta = 1)

and compared with the empirical benchmark log y = c0 + c1 log L_spec + c2 log Teff (3 params)
by BIC. Residual slopes against tau and log M (partial, at fixed position) show what a
predictor misses. beta intervals: 16-84% from 300 bootstrap resamples of the stars.

Writes data/transfer_tests_<grid>.csv and prints the ranking.
Usage: python3 test_transfer.py [Z w]      (default MW 0.0)
"""
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import transfer_models as tm

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
D = f'{HERE}/data'
Z, W = (sys.argv[1], float(sys.argv[2])) if len(sys.argv) > 2 else ('MW', 0.0)
rng = np.random.default_rng(2)
NBOOT = 300
SX, SY = 0.03, 0.08          # scaling of (log Teff, log L_spec) for the coverage test

# ------------------------------------------------------------------ observations (MS stars)
EXT = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
VM = pd.read_csv(f'{OBS}/macroturbulence_evol.csv')
prim = EXT[EXT['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic']) & EXT.a0_ok]
prim = prim[(prim.ev_tau < 1) & (prim.ev_dist < 0.5)]
vm = VM[(VM.ev_tau < 1) & (VM.ev_dist < 0.5) & np.isfinite(VM.vmac) & (VM.vmac > 0)]
OBSSETS = {
    'vmac':   dict(x=vm.logTeff_sp.values, y=vm.logL_sp.values, v=np.log10(vm.vmac.values),
                   tau=vm.ev_tau.values, logM=vm.ev_logM.values, noise=0.08),
    'vmac_hot': dict(x=vm.logTeff_sp.values[vm.logTeff_sp > 4.3], y=vm.logL_sp.values[vm.logTeff_sp > 4.3],
                     v=np.log10(vm.vmac.values[vm.logTeff_sp > 4.3]), tau=vm.ev_tau.values[vm.logTeff_sp > 4.3],
                     logM=vm.ev_logM.values[vm.logTeff_sp > 4.3], noise=0.08),
    'alpha0': dict(x=prim.lT.values, y=prim.lL.values, v=prim.log_alpha0.values,
                   tau=prim.ev_tau.values, logM=prim.ev_logM.values, noise=0.19),
    'nuchar': dict(x=prim.lT.values, y=prim.lL.values, v=prim.log_nuchar.values,
                   tau=prim.ev_tau.values, logM=prim.ev_logM.values, noise=0.19),
}

# ------------------------------------------------------------------ models
MSFILE = os.environ.get('MESA_MS', f'{D}/mesa_ms.csv')
ms = pd.read_csv(MSFILE)
sub = ms[(ms.Z == Z) & (ms.w == W)].reset_index(drop=True)
q, V, A, N = tm.all_models(sub)
pts = np.c_[sub.logTeff.values / SX, sub.logLspec.values / SY]
tree = cKDTree(pts)
present = sub.FeCZ_present.values.astype(float)
pres_int = LinearNDInterpolator(pts, present)


def predict(X, needs_fecz, ox, oy):
    """Model X at the observed positions; 0 where the FeCZ is absent; NaN off the grid."""
    P = np.c_[ox / SX, oy / SY]
    dist, _ = tree.query(P)
    ok = np.isfinite(X) & (X > 0)
    out = np.full(len(ox), np.nan)
    if ok.sum() > 50:
        out = 10**LinearNDInterpolator(pts[ok], np.log10(X[ok]))(P)
    if needs_fecz:
        pf = pres_int(P)
        out = np.where(np.isfinite(pf) & (pf < 0.5), 0.0, out)
    out[dist > 1.5] = np.nan
    return out


def bic(res, k):
    n = len(res)
    return n * np.log(np.sum(res**2) / n) + k * np.log(n)


def fit_floor(X, y, beta=None, quad=False):
    lx = np.log10(np.where(X > 0, X, np.nan))
    pos = X > 0

    def model(p):
        f, a = p[0], p[1]
        b = beta if beta is not None else p[2]
        term = np.where(pos, 10**(a + b * np.nan_to_num(lx, nan=0.0)), 0.0)
        if quad:   # independent fields add in quadrature
            return 0.5 * np.log10(10**(2 * f) + term**2)
        return np.log10(10**f + term)

    f0 = np.percentile(y, 10)
    if pos.sum() > 5:
        b0 = beta if beta is not None else np.polyfit(lx[pos], y[pos], 1)[0]
        a0 = np.median(y[pos] - b0 * lx[pos])
    else:
        b0, a0 = 1.0, 0.0
    p0 = [f0, a0] + ([] if beta is not None else [b0])
    r = least_squares(lambda p: y - model(p), p0, loss='linear')
    return r.x, y - model(r.x)


def fit_lin(X, y, beta=None):
    lx = np.log10(X)
    if beta is None:
        c = np.polyfit(lx, y, 1)
        return np.array([c[1], c[0]]), y - np.polyval(c, lx)
    a = np.mean(y - beta * lx)
    return np.array([a, beta]), y - a - beta * lx


def partial_slopes(res, o, sel):
    """Slopes of the residual on tau and log M (joint regression) for the fitted stars."""
    Xm = np.c_[np.ones(sel.sum()), o['tau'][sel], o['logM'][sel]]
    c = np.linalg.lstsq(Xm, res, rcond=None)[0]
    return c[1], c[2]




def main():
    rows = []
    for obs_name, o in OBSSETS.items():
        cands = N if obs_name == 'nuchar' else (V if obs_name.startswith('vmac') else A)
        # benchmark on the same stars as each candidate is computed below
        for name, (X, desc) in cands.items():
            needs = not name.startswith('core') and name not in ('dynamical', 'nu_max scaling')
            Xs = predict(X, needs, o['x'], o['y'])
            sel = np.isfinite(Xs) & np.isfinite(o['v'])
            if obs_name == 'nuchar':
                sel &= Xs > 0
            if sel.sum() < 30:
                continue
            y, Xv = o['v'][sel], Xs[sel]
            # benchmark plane on the same stars
            Bm = np.c_[np.ones(sel.sum()), o['y'][sel], o['x'][sel]]
            bres = y - Bm @ np.linalg.lstsq(Bm, y, rcond=None)[0]
            row = dict(observable=obs_name, predictor=name, description=desc, n=int(sel.sum()),
                       n_below_onset=int((Xv == 0).sum()), bic_plane=bic(bres, 3), rms_plane=np.std(bres))
            pos = Xv > 0
            row['spearman'] = spearmanr(Xv[pos], y[pos]).correlation if pos.sum() > 10 else np.nan
            for mode in ('free', 'linear', 'quad'):
                if mode == 'quad' and obs_name == 'nuchar':
                    continue
                beta = None if mode == 'free' else 1.0
                try:
                    if obs_name == 'nuchar':
                        p, res = fit_lin(Xv, y, beta); k = 2 if beta is None else 1
                        row[f'{mode}_a'] = p[0]; row[f'{mode}_beta'] = p[1]
                    else:
                        p, res = fit_floor(Xv, y, beta, quad=(mode == 'quad')); k = 3 if beta is None else 2
                        row[f'{mode}_floor'] = 10**p[0]; row[f'{mode}_a'] = p[1]
                        row[f'{mode}_beta'] = p[2] if beta is None else 1.0
                except Exception as e:
                    continue
                row[f'{mode}_rms'] = np.std(res)
                row[f'{mode}_dBIC'] = bic(res, k) - row['bic_plane']      # < 0: better than the plane
                row[f'{mode}_res_tau'], row[f'{mode}_res_logM'] = partial_slopes(res, o, sel)
                if mode == 'free':
                    bs = []
                    idx_all = np.arange(sel.sum())
                    for _ in range(NBOOT):
                        ii = rng.choice(idx_all, len(idx_all))
                        try:
                            pb = (fit_lin if obs_name == 'nuchar' else fit_floor)(Xv[ii], y[ii])[0]
                            bs.append(pb[1] if obs_name == 'nuchar' else pb[2])
                        except Exception:
                            pass
                    row['free_beta_lo'], row['free_beta_hi'] = np.percentile(bs, [16, 84])
            rows.append(row)

    res = pd.DataFrame(rows)
    out = f'{D}/transfer_tests_{Z}_w{W}.csv'
    res.to_csv(out, index=False, float_format='%.4g')
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 28)
    for obs_name in OBSSETS:
        r = res[res.observable == obs_name].sort_values('linear_dBIC')
        print(f'\n=== {obs_name}  (N ~ {r.n.max()}, per-star noise ~{OBSSETS[obs_name]["noise"]} dex; '
              f'benchmark plane rms {r.rms_plane.iloc[0]:.3f})')
        cols = ['predictor', 'n', 'n_below_onset', 'spearman', 'linear_a', 'linear_rms', 'linear_dBIC', 'quad_a', 'quad_dBIC',
                'linear_res_tau', 'linear_res_logM', 'free_beta', 'free_beta_lo', 'free_beta_hi', 'free_rms', 'free_dBIC']
        if obs_name != 'nuchar':
            cols.insert(4, 'linear_floor')
        else:
            cols = [c for c in cols if not c.startswith('quad')]
        print(r[cols].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    print('\nwritten', out)


if __name__ == '__main__':
    main()
