"""Side-by-side scenario results for two grid versions (default: data/ = v1, data_v2/ = v2).

Usage: python3 compare_v1_v2.py [DIR_A] [DIR_B]
Prints, per rotation rate, each scenario's per-observable dBIC (against the empirical planes), the joint
dBIC, the implied sensitivity ratios, and the number of observed stars the models covered.
Observed ratios (REPORT_transfer.md): 5.3 [4.9, 5.8] / 6.4 [4.9, 8.6] / 2.2 [1.8, 2.7] along L / tau / M.
"""
import os
import sys
import pandas as pd

A = sys.argv[1] if len(sys.argv) > 1 else 'data'
B = sys.argv[2] if len(sys.argv) > 2 else 'data_v2'
cols = ['vmac_n', 'vmac_dBIC', 'alpha0_n', 'alpha0_dBIC', 'nuchar_n', 'nuchar_dBIC', 'joint_dBIC',
        'ratio_logL', 'ratio_tau', 'ratio_logM']
short = {'vmac_n': 'n_v', 'vmac_dBIC': 'dBIC_v', 'alpha0_n': 'n_a', 'alpha0_dBIC': 'dBIC_a', 'nuchar_n': 'n_nu',
         'nuchar_dBIC': 'dBIC_nu', 'joint_dBIC': 'joint', 'ratio_logL': 'r_L', 'ratio_tau': 'r_tau',
         'ratio_logM': 'r_M'}
pd.set_option('display.width', 200)
for w in ('0.0', '0.2', '0.4', '0.6'):
    frames = {}
    for tag, d in (('A', A), ('B', B)):
        f = f'{d}/scenarios_MW_w{w}.csv'
        if os.path.exists(f):
            frames[tag] = pd.read_csv(f).set_index('scenario')[cols].rename(columns=short)
    if not frames:
        continue
    print(f'\n=== MW w{w}   (A = {A}, B = {B})')
    out = pd.concat(frames, axis=1).swaplevel(axis=1).sort_index(axis=1, level=0, sort_remaining=False)
    order = [c for c in short.values() for t in ('A', 'B') if (c, t) in out.columns]
    out = out[[(c, t) for c in short.values() for t in ('A', 'B') if (c, t) in out.columns]]
    print(out.round(2).to_string())
