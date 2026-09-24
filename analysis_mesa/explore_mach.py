"""Exploration: v_macro / c_s(model) against FeCZ quantities; alpha0 against observed surface Mach."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import transfer_models as tm
from test_transfer import OBSSETS, predict, V, A, N, sub, prim, VM, fit_floor, bic
q = tm.inputs(sub)
cs = q['c_s'] / 1e5
pred = lambda X, needs=True, o=None: predict(X, needs, o['x'], o['y'])

print('--- v_macro / c_s  (MS stars; c_s = model photospheric sound speed at the star position)')
for obs in ('vmac', 'vmac_hot'):
    o = OBSSETS[obs]
    cs_s = pred(cs, False, o)
    y = o['v'] - np.log10(cs_s)
    ok = np.isfinite(y)
    print(f'{obs}: N={ok.sum()}, median v_macro/c_s = {10**np.median(y[ok]):.2f} [{10**np.percentile(y[ok],16):.2f}, {10**np.percentile(y[ok],84):.2f}]')
    # does v_macro/c_s scale with c_s itself? (log v = a + b log c_s)
    c = np.polyfit(np.log10(cs_s[ok]), o['v'][ok], 1); r = o['v'][ok] - np.polyval(c, np.log10(cs_s[ok]))
    Bm = np.c_[np.ones(ok.sum()), o['y'][ok], o['x'][ok]]; br = o['v'][ok] - Bm @ np.linalg.lstsq(Bm, o['v'][ok], rcond=None)[0]
    print(f'   log v_macro = a + {c[0]:.2f} log c_s : rms {np.std(r):.3f}, dBIC vs plane {bic(r,2)-bic(br,3):+.1f}')
    tests = {'M_x': q['M_x'], 'M_t': q['M_t'], 'F_c/F': q['FcF'], 'v_x/c_s': q['v_x'] / q['c_s'],
             'Gamma_surf': sub.Gamma_surf.values, 'nHP': q['nHP'], 'H_c/R': q['H_c'] / q['R'],
             'M_t F_c/F': q['M_t'] * q['FcF'], 'depth/R': q['depth'] / q['R']}
    for k, X in tests.items():
        Xs = pred(X, k not in ('Gamma_surf',), o)
        s = ok & np.isfinite(Xs) & (Xs > 0)
        if s.sum() < 20:
            print(f"   v/c_s vs {k:11s}: n/a"); continue
        rho = spearmanr(Xs[s], y[s]).correlation
        cc = np.polyfit(np.log10(Xs[s]), y[s], 1)
        print(f'   v/c_s vs {k:11s}: rho {rho:+.2f}  slope {cc[0]:+.3f}  (N {s.sum()})')

print('\n--- alpha0 vs observed v_macro and model c_s (stars with both)')
vmk = VM.set_index('key')
EXT = pd.read_csv('/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data/rednoise_sHRD_extended_evol.csv')
p = prim.copy()
p['key'] = p.star.str.replace(' ', '', regex=False).str.upper()
vk = VM.copy(); vk['key2'] = vk.key.str.replace(' ', '', regex=False).str.upper()
m = p.merge(vk[['key2', 'vmac']], left_on='key', right_on='key2', how='inner')
m = m[np.isfinite(m.vmac) & (m.vmac > 0)]
cs_m = predict(cs, False, m.lT.values, m.lL.values)
ok = np.isfinite(cs_m)
la, lv, lc = m.log_alpha0.values[ok], np.log10(m.vmac.values[ok]), np.log10(cs_m[ok])
print(f'N = {ok.sum()} MS stars with alpha0, v_macro and model c_s')
X1 = np.c_[np.ones(ok.sum()), lv]; c1 = np.linalg.lstsq(X1, la, rcond=None)[0]; r1 = la - X1 @ c1
X2 = np.c_[np.ones(ok.sum()), lv, lc]; c2 = np.linalg.lstsq(X2, la, rcond=None)[0]; r2 = la - X2 @ c2
X3 = np.c_[np.ones(ok.sum()), lv - lc]; c3 = np.linalg.lstsq(X3, la, rcond=None)[0]; r3 = la - X3 @ c3
print(f'  log a0 = a + {c1[1]:.2f} log v_macro                  rms {np.std(r1):.3f}  BIC {bic(r1,2):.1f}')
print(f'  log a0 = a + {c2[1]:.2f} log v_macro {c2[2]:+.2f} log c_s   rms {np.std(r2):.3f}  BIC {bic(r2,3):.1f}')
print(f'  log a0 = a + {c3[1]:.2f} log (v_macro/c_s)            rms {np.std(r3):.3f}  BIC {bic(r3,2):.1f}')
for b in (1, 2):
    r = la - np.mean(la - b * (lv - lc)) - b * (lv - lc)
    print(f'  fixed: a0 ~ (v_macro/c_s)^{b}: rms {np.std(r):.3f}  BIC {bic(r,1):.1f}')
