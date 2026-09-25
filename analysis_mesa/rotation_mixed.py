#!/usr/bin/env python3
"""
Star-by-star rotation for the transfer-model scenarios (v2 grid).

Each observed star gets its own omega/omega_crit: its equatorial velocity is estimated as
v_eq = v sin i x 4/pi (mean over random inclinations), and omega is the value at which the model
surface rotation velocity at the star's (log Teff, log L_spec) position equals v_eq, interpolated between
the omega = 0, 0.2, 0.4, 0.6 sub-grids and clipped to [0, 0.6]. Stars without v sin i get the sample median.
Every predictor is then interpolated in omega between the four maps (in log where all four are > 0,
linearly otherwise). The scenarios of scenarios.py are scored with beta = 1 (normalisation + floor), against
the empirical L_spec + Teff plane, on the SAME stars for: per-star omega, and each fixed omega.

Usage (from analysis_mesa/):  RN_DATA=data_v2 python3 rotation_mixed.py
Writes $RN_DATA/rotation_mixed.csv and $RN_DATA/rotation_mixed_omega.csv (per-star omega).
"""
import os
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.environ.setdefault('RN_DATA', f'{HERE}/data_v2')
os.environ['MESA_MS'] = f'{D}/mesa_ms_MWw0.0.csv'
import transfer_models as tm                      # noqa: E402
import test_transfer as tt                        # noqa: E402  (observations, fit helpers)
from scenarios import SCEN, NOFECZ                # noqa: E402

OMEGAS = [0.0, 0.2, 0.4, 0.6]
SX, SY = tt.SX, tt.SY


class OmegaMap:
    """Interpolation of model quantities of one rotation sub-grid onto (log Teff, log L_spec)."""

    def __init__(self, w):
        df = pd.read_csv(f'{D}/mesa_ms_MWw{w}.csv')
        self.df = df
        self.q, self.V, self.A, self.N = tm.all_models(df)
        self.pts = np.c_[df.logTeff.values / SX, df.logLspec.values / SY]
        self.tree = cKDTree(self.pts)
        self.pres = LinearNDInterpolator(self.pts, df.FeCZ_present.values.astype(float))
        self._cache = {}

    def predict(self, X, needs_fecz, ox, oy):
        P = np.c_[ox / SX, oy / SY]
        dist, _ = self.tree.query(P)
        ok = np.isfinite(X) & (X > 0)
        out = np.full(len(ox), np.nan)
        if ok.sum() > 50:
            out = 10**LinearNDInterpolator(self.pts[ok], np.log10(X[ok]))(P)
        if needs_fecz:
            pf = self.pres(P)
            out = np.where(np.isfinite(pf) & (pf < 0.5), 0.0, out)
        out[dist > 1.5] = np.nan
        return out

    def get(self, kind, key, ox, oy, tag):
        ck = (kind, key, tag)
        if ck not in self._cache:
            X = {'V': self.V, 'A': self.A, 'N': self.N}[kind][key][0]
            needs = not key.startswith(NOFECZ)
            self._cache[ck] = self.predict(X, needs, ox, oy)
        return self._cache[ck]


maps = {w: OmegaMap(w) for w in OMEGAS}

# observed stars (same selections as test_transfer), with v sin i
STARS = {
    'vmac': dict(x=tt.vm.logTeff_sp.values, y=tt.vm.logL_sp.values, v=np.log10(tt.vm.vmac.values),
                 vsini=tt.vm.vsini.values),
    'alpha0': dict(x=tt.prim.lT.values, y=tt.prim.lL.values, v=tt.prim.log_alpha0.values, vsini=tt.prim.vsini.values),
    'nuchar': dict(x=tt.prim.lT.values, y=tt.prim.lL.values, v=tt.prim.log_nuchar.values, vsini=tt.prim.vsini.values),
}
KIND = {'vmac': 'V', 'alpha0': 'A', 'nuchar': 'N'}


def star_omega(o, tag):
    """Per-star omega from v sin i and the model surface velocity at the star's position."""
    veq = np.vstack([np.zeros(len(o['x']))] +
                    [maps[w].predict(maps[w].df.vsurf_kms.values, False, o['x'], o['y']) for w in OMEGAS[1:]]).T
    vest = o['vsini'] * 4 / np.pi
    om = np.full(len(vest), np.nan)
    for i in range(len(vest)):
        if np.isfinite(vest[i]) and np.all(np.isfinite(veq[i])) and np.all(np.diff(veq[i]) > 0):
            om[i] = np.interp(vest[i], veq[i], OMEGAS)        # clips to [0, 0.6]
    med = np.nanmedian(om)
    return np.where(np.isfinite(om), om, med), om


def interp_omega(vals, om):
    """vals: (n_stars, 4) predictor at the four omegas; om: per-star omega."""
    out = np.full(len(om), np.nan)
    for i in range(len(om)):
        v = vals[i]
        if not np.all(np.isfinite(v)):
            continue
        if np.all(v > 0):
            out[i] = 10**np.interp(om[i], OMEGAS, np.log10(v))
        else:
            out[i] = np.interp(om[i], OMEGAS, v)
    return out


def score(obs, X, y, ox, oy):
    """beta = 1 fit (normalisation + floor; frequencies: offset only) and dBIC against the plane."""
    Bm = np.c_[np.ones(len(y)), oy, ox]
    bres = y - Bm @ np.linalg.lstsq(Bm, y, rcond=None)[0]
    if obs == 'nuchar':
        p, res = tt.fit_lin(X, y, 1.0); k = 1
        pf, resf = tt.fit_lin(X, y)
        return dict(dBIC=tt.bic(res, k) - tt.bic(bres, 3), a=p[0], rms=np.std(res), free_beta=pf[1],
                    free_dBIC=tt.bic(resf, 2) - tt.bic(bres, 3))
    p, res = tt.fit_floor(X, y, 1.0); k = 2
    return dict(dBIC=tt.bic(res, k) - tt.bic(bres, 3), a=p[1], floor=10**p[0], rms=np.std(res))


rows, om_rows = [], []
for obs, o in STARS.items():
    om, om_raw = star_omega(o, obs)
    om_rows.append(pd.DataFrame(dict(observable=obs, logTeff=o['x'], logLspec=o['y'], vsini=o['vsini'],
                                     omega=om, omega_from_vsini=np.isfinite(om_raw))))
    print(f'{obs}: {np.isfinite(om_raw).sum()} of {len(om)} stars with an omega from v sin i; '
          f'omega median {np.median(om):.2f}, 16-84% {np.percentile(om, 16):.2f}-{np.percentile(om, 84):.2f}, '
          f'at the 0.6 cap: {(om >= 0.599).mean():.2f}')
    for sname, keys in SCEN.items():
        key = keys[{'vmac': 0, 'alpha0': 1, 'nuchar': 2}[obs]]
        vals = np.vstack([maps[w].get(KIND[obs], key, o['x'], o['y'], obs) for w in OMEGAS]).T
        cands = {'per-star omega': interp_omega(vals, om)}
        for j, w in enumerate(OMEGAS):
            cands[f'fixed omega {w}'] = vals[:, j]
        sel = np.isfinite(o['v']) & np.all(np.isfinite(vals), axis=1)
        if obs == 'nuchar':
            sel &= np.all(vals > 0, axis=1)
        for mode, X in cands.items():
            s = score(obs, X[sel], o['v'][sel], o['x'][sel], o['y'][sel])
            rows.append(dict(scenario=sname, observable=obs, mode=mode, n=int(sel.sum()), **s))

res = pd.DataFrame(rows)
res.to_csv(f'{D}/rotation_mixed.csv', index=False, float_format='%.4g')
pd.concat(om_rows).to_csv(f'{D}/rotation_mixed_omega.csv', index=False, float_format='%.4g')
joint = res.pivot_table(index=['scenario', 'mode'], columns='observable', values='dBIC')
joint['joint'] = joint.sum(axis=1)
pd.set_option('display.width', 200)
print('\ndBIC against the empirical plane (beta = 1; < 0 is better), identical stars in every row of a scenario')
print(joint.round(1).to_string())
print('\nnu_char free beta by mode:')
print(res[res.observable == 'nuchar'].pivot_table(index='scenario', columns='mode', values='free_beta').round(2).to_string())
