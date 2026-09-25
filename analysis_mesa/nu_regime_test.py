"""Which model frequency tracks the observed nu_char? Tests on the regime table (regime_params.py).

Observed sample: primary red-noise MS stars (Bowman+2020 and Shen+2024 Galactic, a0_ok, ev_tau < 1,
ev_dist < 0.5), as in test_transfer.py. Model quantities are interpolated linearly in
(log Teff / 0.03, log L_spec / 0.08) over the MS rows of the regime table (NaN off the grid).
For each candidate X: log nu = a + beta log X, with beta = 1 (fixed) and free; plus the two-timescale model
log nu = a + b1 log nu_turn + b2 log nu_th. dBIC is against the empirical plane
log nu = c0 + c1 log L_spec + c2 log Teff on the same stars.

Usage (from analysis_mesa/): python3 nu_regime_test.py [data_v2/regime_MWw0.0.csv]
"""
import sys
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
SX, SY = 0.03, 0.08
f = sys.argv[1] if len(sys.argv) > 1 else 'data_v2/regime_MWw0.0.csv'

EXT = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
prim = EXT[EXT['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic']) & EXT.a0_ok]
prim = prim[(prim.ev_tau < 1) & (prim.ev_dist < 0.5)]
ox, oy, ov = prim.lT.values, prim.lL.values, prim.log_nuchar.values

d = pd.read_csv(f)
ms = d[(d.Xc > 1e-3) & (d.Xc < 0.7)].reset_index(drop=True)
pts = np.c_[ms.logTeff / SX, ms.logLspec / SY]
tree = cKDTree(pts)
P = np.c_[ox / SX, oy / SY]
dist, _ = tree.query(P)


def at_obs(col):
    x = ms[col].values
    ok = np.isfinite(x) & (x > 0)
    out = 10 ** LinearNDInterpolator(pts[ok], np.log10(x[ok]))(P)
    out[dist > 1.5] = np.nan
    return out


def bic(res, k):
    return len(res) * np.log(np.mean(res ** 2)) + k * np.log(len(res))


def ols(Xcols, y):
    X = np.c_[np.ones(len(y)), *Xcols]
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return c, y - X @ c


cands = ['nu_turn', 'nu_th_peak', 'nu_th_top', 'nu_th_zone', 'nu_th_above']
vals = {c: at_obs(c) for c in cands}
for c in ('Gamma_Fe', 'tau_ratio', 'Y_Fe', 'mach_max'):
    vals[c] = at_obs(c)
ok = np.all([np.isfinite(vals[c]) & (vals[c] > 0) for c in cands], axis=0)
y = ov[ok]
print(f'{ok.sum()} of {len(ov)} stars covered by FeCZ-bearing MS models')
_, rp = ols([oy[ok], ox[ok]], y)
bp = bic(rp, 3)
print(f'plane: rms {np.std(rp):.3f}')
print(f"{'predictor':30s} {'a(beta=1)':>9s} {'dBIC':>7s} | {'beta_free':>9s} {'dBIC':>7s} {'rms':>6s} {'spearman':>8s}")
from scipy.stats import spearmanr
for c in cands:
    lx = np.log10(vals[c][ok])
    a1 = np.mean(y - lx); r1 = y - lx - a1
    cf, rf = ols([lx], y)
    print(f'{c:30s} {a1:9.3f} {bic(r1, 1) - bp:7.1f} | {cf[1]:9.3f} {bic(rf, 2) - bp:7.1f} {np.std(rf):6.3f} {spearmanr(lx, y)[0]:8.3f}')
for th in ('nu_th_peak', 'nu_th_top', 'nu_th_zone'):
    cf, rf = ols([np.log10(vals['nu_turn'][ok]), np.log10(vals[th][ok])], y)
    print(f'{"turn^b1 x " + th + "^b2":30s} {"":>9s} {"":>7s} | b1={cf[1]:.2f} b2={cf[2]:.2f} {bic(rf, 3) - bp:7.1f} {np.std(rf):6.3f}')
# does the turnover offset depend on the regime?
off = y - np.log10(vals['nu_turn'][ok])
print('\nresidual log(nu_char / nu_turn) against regime parameters (Spearman rho):')
for c in ('Gamma_Fe', 'tau_ratio', 'Y_Fe', 'mach_max'):
    v = vals[c][ok]
    s = np.isfinite(v) & (v > 0)
    print(f'  {c:10s} rho = {spearmanr(np.log10(v[s]), off[s])[0]:+.3f} (n = {s.sum()})')
