#!/usr/bin/env python3
"""Where do the FeCZ regime parameters cross their thresholds along log L_spec?

Hot main-sequence models (log Teff > 4.3, X_c > 1e-3) of one regime table (regime_params.py): median of each
parameter in 0.05-dex bins of log L_spec, the first log L_spec where the median crosses each threshold (linear
interpolation between bin centres), and a 16-84% interval from 500 bootstrap resamples of whole tracks.

Usage (from analysis_mesa/):
    python3 regime_crossings.py data_v2/regime_MWw0.0.csv [data_v2/regime_LMCw0.0.csv ...]
Writes <table>_crossings.csv next to each input table.
"""
import sys
import numpy as np
import pandas as pd

THRESH = [('Gamma_Fe', 0.4), ('Gamma_Fe', 0.8), ('Gamma_Fe', 0.9), ('Gamma_Fe', 1.0),
          ('mach_max', 0.2), ('v_c', 3.0), ('Y_Fe', 1.0)]
EDGES = np.arange(2.0, 4.8, 0.05)
NBOOT = 500
rng = np.random.default_rng(5)


def crossing(df, col, thr):
    b = pd.cut(df.logLspec, EDGES)
    med = df.groupby(b, observed=True)[col].median()
    x = np.array([iv.mid for iv in med.index])
    y = med.values
    ok = np.isfinite(y)
    x, y = x[ok], y[ok]
    above = np.nonzero(y >= thr)[0]
    if len(above) == 0 or above[0] == 0:
        return np.nan
    i = above[0]
    return x[i - 1] + (thr - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])


for fn in sys.argv[1:]:
    df = pd.read_csv(fn)
    hot = df[(df.logTeff > 4.3) & (df.Xc > 1e-3)].copy()
    if 'v_c' in hot:
        hot['v_c'] = hot['v_c'].fillna(0.0)
    tracks = hot.Minit.unique()
    rows = []
    for col, thr in THRESH:
        if col not in hot:
            continue
        c0 = crossing(hot, col, thr)
        bs = []
        for _ in range(NBOOT):
            pick = rng.choice(tracks, len(tracks))
            sub = pd.concat([hot[hot.Minit == m] for m in pick])
            bs.append(crossing(sub, col, thr))
        lo, hi = np.nanpercentile(bs, [16, 84]) if np.isfinite(bs).any() else (np.nan, np.nan)
        rows.append(dict(parameter=col, threshold=thr, logLspec=c0, lo=lo, hi=hi,
                         n_models=len(hot), n_tracks=len(tracks)))
    out = pd.DataFrame(rows)
    out.to_csv(fn.replace('.csv', '_crossings.csv'), index=False, float_format='%.3f')
    print(f'\n{fn}  ({len(hot)} hot MS profiles, {len(tracks)} tracks)')
    print(out.to_string(index=False, float_format=lambda v: f'{v:.2f}'))
