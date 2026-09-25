"""Restrict v1 main-sequence extracts to the coverage of the (possibly incomplete) v2 grid.

For each rotation rate, keep only v1 tracks whose initial mass exists in v2, and cut each track at the
lowest X_c that v2 has reached for that mass. Writes data_v1matched/mesa_ms_MWw<w>.csv, used to compare v1
and v2 physics on the same observed stars (REPORT_transfer_v2.md). Then run test_transfer.py / scenarios.py
with RN_DATA=data_v1matched and MESA_MS=data_v1matched/mesa_ms_MWw<w>.csv.

Usage (from analysis_mesa/): python3 make_v1matched.py
"""
import os
import pandas as pd

os.makedirs('data_v1matched', exist_ok=True)
for w in ('0.0', '0.2', '0.4'):
    v1 = pd.read_csv(f'data/mesa_ms_MWw{w}.csv')
    v2 = pd.read_csv(f'data_v2/mesa_ms_MWw{w}.csv')
    keep = []
    for M, t in v1.groupby('Minit'):
        t2 = v2[v2.Minit == M]
        if len(t2):
            keep.append(t[t.Xc >= t2.Xc.min() - 1e-4])
    m = pd.concat(keep)
    m.to_csv(f'data_v1matched/mesa_ms_MWw{w}.csv', index=False, float_format='%.6g')
    print(f'w{w}: v1 {len(v1)} rows -> matched {len(m)} ({m.Minit.nunique()} masses); v2 {len(v2)} rows')
