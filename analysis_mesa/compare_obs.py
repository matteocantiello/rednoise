#!/usr/bin/env python3
"""
MESA grid vs observations (handoff §8.3). Reads data/mesa_ms.csv, data/mesa_post.csv (from
extract_mesa.py) and the observational tables of the handoff package. Writes:

  data/mesa_onset.csv         (1) onset: model thresholds in log L_spec per sub-grid
  data/mesa_slopes.csv        (2)(3)(5) model regressions (L_spec+Teff basis, tau+logM basis)
  data/mesa_ratios.csv        (3) model vs observed sensitivity ratios
  data/mesa_nuc_offset.csv    (4) log nu_char(GP) - log nu_c(model) over the sHRD
  data/mesa_coolside.csv      (6) zone properties in the cool-side regions
  data/obs_binned.csv         observed binned medians used in the figure
  figures/fig8_mesa.pdf/.png  proposed Fig. 8

Model conventions: MS rows only (unless stated), FeCZ present, primary model sub-grid
MW w0.0 (alpha_MLT = 1.6, no MLT++); w0.2/w0.4 as rotation sensitivity; each track weighted
equally and sampled uniformly in MS age (rows of mesa_ms.csv).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
D = os.environ.get('RN_DATA', f'{HERE}/data')        # RN_DATA/RN_FIG: e.g. data_v2, figures_v2
FIG = os.environ.get('RN_FIG', f'{HERE}/figures')
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(1)
NBOOT = 1000

ms = pd.read_csv(f'{D}/mesa_ms.csv')
post = pd.read_csv(f'{D}/mesa_post.csv')
EXT = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
VM = pd.read_csv(f'{OBS}/macroturbulence_evol.csv')
GP = np.load(f'{OBS}/hrd_interpolated_grids_v2.npz')

# ------------------------------------------------------------------ observed reference values
primary = EXT[(EXT['sample'] == 'Bowman2020_Galactic') | (EXT['sample'] == 'Shen2024_Galactic')]
primary = primary[primary.a0_ok & np.isfinite(primary.log_alpha0)]
a0_ob = EXT[EXT.a0_ok & (EXT.lT >= 4.0) & np.isfinite(EXT.alpha0)]            # 'log alpha0 (OB)'
vm_hot = VM[(VM.logTeff_sp > 4.3) & np.isfinite(VM.vmac) & np.isfinite(VM.logL_sp)]
print(f'observed: primary {len(primary)}, alpha0 OB {len(a0_ob)}, vmac hot {len(vm_hot)}')

OBS_SLOPES = {  # from hrd_gradients.csv / evolution_test_full.csv (paper Table fits)
    ('vmac', 'logL'): 0.291, ('alpha0', 'logL'): 1.545,
    ('vmac', 'tau'): 0.272, ('alpha0', 'tau'): 1.732,
    ('vmac', 'logM'): 0.603, ('alpha0', 'logM'): 1.321,
    ('nuchar', 'logTeff'): 1.531, ('nuchar', 'logL'): -0.233, ('nuchar', 'tau'): -0.326,
}
OBS_CI = {  # 16-84% of the coefficients
    ('vmac', 'tau'): (0.222, 0.321), ('alpha0', 'tau'): (1.567, 1.913),
    ('vmac', 'logM'): (0.572, 0.633), ('alpha0', 'logM'): (1.111, 1.523),
    ('alpha0', 'logL'): (1.430, 1.680),
}
OBS_HINGE = dict(onset=3.16, onset_lo=2.94, onset_hi=3.44, sat=3.70, sat_lo=3.70, sat_hi=3.84)
# domain of the primary sample (5-95%)
DOM_L = tuple(np.percentile(primary.lL, [5, 95]))
DOM_T = tuple(np.percentile(primary.lT, [5, 95]))
DOM_M = tuple(np.percentile(EXT.loc[EXT.ev_tau < 1, 'ev_logM'].dropna(), [5, 95]))
print(f'primary domain: logL {DOM_L[0]:.2f}-{DOM_L[1]:.2f}, logT {DOM_T[0]:.2f}-{DOM_T[1]:.2f}, '
      f'MS log M {DOM_M[0]:.2f}-{DOM_M[1]:.2f}')


def binned(x, y, edges, nmin=5):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (x >= a) & (x < b) & np.isfinite(y)
        if k.sum() >= nmin:
            out.append((0.5 * (a + b), np.median(y[k]), *np.percentile(y[k], [16, 84]), k.sum()))
    return np.array(out)


def wls(X, y, w=None):
    X = np.c_[np.ones(len(y)), X]
    w = np.ones(len(y)) if w is None else w
    sw = np.sqrt(w)
    return np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]


def fit_boot(df, cols, ycol, groups='track'):
    """Coefficients of ycol on cols, with a bootstrap over tracks (not rows) for the interval."""
    ok = np.isfinite(df[ycol]) & np.all(np.isfinite(df[cols]), axis=1)
    df = df[ok]
    X, y = df[cols].values, df[ycol].values
    coef = wls(X, y)[1:]
    keys = df[groups].unique()
    members = [np.where(df[groups].values == k)[0] for k in keys]
    boots = []
    for _ in range(NBOOT):
        idx = np.concatenate([members[i] for i in rng.integers(0, len(keys), len(keys))])
        boots.append(wls(X[idx], y[idx])[1:])
    lo, hi = np.percentile(boots, [16, 84], axis=0)
    return coef, lo, hi, len(df), len(keys)


def prep(df):
    df = df.copy()
    df['track'] = df.Z + '/' + df.w.astype(str) + '/' + df.Minit.astype(str)
    df['logM'] = np.log10(df.Minit)
    df['log_v'] = np.log10(df.FeCZ_vmax_kms)
    df['log_F'] = np.log10(df.FeCZ_FcF)
    df['log_rho'] = np.log10(df.FeCZ_rho)
    df['log_nuc'] = np.log10(df.FeCZ_nuc_d)
    df['log_mach'] = np.log10(df.FeCZ_mach_max)
    return df


ms = prep(ms)
post = prep(post)
MW0 = ms[(ms.Z == 'MW') & (ms.w == 0.0)]

# ------------------------------------------------------------------ (1) onset
onset_rows = []
for (z, w), sub in ms.groupby(['Z', 'w']):
    hot = sub[sub.logTeff > 4.3]
    edges = np.arange(1.8, 4.6, 0.05)
    x = hot.logLspec.values
    pres = binned(x, hot.FeCZ_present.astype(float).values, edges, nmin=3)
    v = binned(x, np.nan_to_num(hot.FeCZ_vmax_kms.values, nan=0.0), edges, nmin=3)
    F = binned(x, np.nan_to_num(hot.FeCZ_FcF.values, nan=0.0), edges, nmin=3)

    def cross(tab, thr):
        k = np.where(tab[:, 1] >= thr)[0]
        return tab[k[0], 0] if len(k) else np.nan

    onset_rows.append(dict(Z=z, w=w, L_first_FeCZ=cross(pres, 0.5),
                           L_v1=cross(v, 1.0), L_v3=cross(v, 3.0), L_v10=cross(v, 10.0),
                           L_F1e3=cross(F, 1e-3), L_F1e2=cross(F, 1e-2), L_F1e1=cross(F, 1e-1)))
onset = pd.DataFrame(onset_rows)
onset.to_csv(f'{D}/mesa_onset.csv', index=False, float_format='%.3f')
print('\n(1) ONSET: median over MS models with log Teff > 4.3, log L_spec where the median first exceeds')
print(onset.to_string(index=False))

# ------------------------------------------------------------------ (2)(3)(5) regressions
slope_rows = []


def add(label, sub, cols, ycol):
    coef, lo, hi, n, nt = fit_boot(sub, cols, ycol)
    for c, a, l, h in zip(cols, coef, lo, hi):
        slope_rows.append(dict(sample=label, y=ycol, x=c, coef=a, lo=l, hi=h, n_rows=n, n_tracks=nt))
    return dict(zip(cols, coef))


results = {}
for w in (0.0, 0.2, 0.4):
    sub = ms[(ms.Z == 'MW') & (ms.w == w) & ms.FeCZ_present]
    dom = sub[sub.logLspec.between(*DOM_L) & sub.logTeff.between(*DOM_T)]
    hot = sub[sub.logTeff > 4.3]
    lab = f'MW w{w}'
    r = {}
    for y in ('log_v', 'log_F', 'log_rho', 'log_nuc', 'log_mach'):
        r[('L', y)] = add(f'{lab} | sHRD domain', dom, ['logLspec', 'logTeff'], y)
    # saturation: two regimes along L_spec for hot stars
    for tag, (a, b) in (('below', (OBS_HINGE['onset'], OBS_HINGE['sat'])), ('above', (OBS_HINGE['sat'], 4.4))):
        reg = hot[hot.logLspec.between(a, b)]
        for y in ('log_v', 'log_F', 'log_mach'):
            r[(tag, y)] = add(f'{lab} | hot, {a:.2f}<logL<{b:.2f}', reg, ['logLspec', 'logTeff'], y)
    # evolution at fixed mass (complete MS tracks inside the observed mass range)
    ev = sub[sub.ms_complete & sub.logM.between(*DOM_M) & np.isfinite(sub.tau)]
    for y in ('log_v', 'log_F', 'log_rho', 'log_nuc'):
        r[('tauM', y)] = add(f'{lab} | tau+logM', ev, ['tau', 'logM'], y)
    results[w] = r

slopes = pd.DataFrame(slope_rows)
slopes.to_csv(f'{D}/mesa_slopes.csv', index=False, float_format='%.4f')


def S(w, key, y, x):
    s = slopes[(slopes['sample'].str.startswith(f'MW w{w} |')) & (slopes.y == y) & (slopes.x == x)]
    tagmap = {'L': 'sHRD domain', 'below': f'{OBS_HINGE["onset"]:.2f}<', 'above': f'{OBS_HINGE["sat"]:.2f}<',
              'tauM': 'tau+logM'}
    s = s[s['sample'].str.contains(tagmap[key], regex=False)]
    return s.iloc[0]


print('\n(2) SATURATION (hot MS models, L_spec+Teff basis): d log Q / d log L_spec')
for w in (0.0, 0.2, 0.4):
    for y in ('log_v', 'log_F', 'log_mach'):
        b, a = S(w, 'below', y, 'logLspec'), S(w, 'above', y, 'logLspec')
        print(f'  w={w} {y:9s} below {b.coef:+.2f} [{b.lo:+.2f},{b.hi:+.2f}]  above {a.coef:+.2f} [{a.lo:+.2f},{a.hi:+.2f}]')
print(f'  observed v_macro: below +0.81, above +0.29 (hot IACOB, two-hinge fit)')

# sensitivity ratios
ratio_rows = []
for axis, key, x in (('logL_spec', 'L', 'logLspec'), ('tau', 'tauM', 'tau'), ('log M', 'tauM', 'logM')):
    obs_x = {'logL_spec': 'logL', 'tau': 'tau', 'log M': 'logM'}[axis]
    oa, ov = OBS_SLOPES[('alpha0', obs_x)], OBS_SLOPES[('vmac', obs_x)]
    ci_a = OBS_CI.get(('alpha0', obs_x), (oa, oa)); ci_v = OBS_CI.get(('vmac', obs_x), (ov, ov))
    row = dict(axis=axis, obs_dalpha0=oa, obs_dvmac=ov, obs_ratio=oa / ov,
               obs_ratio_lo=ci_a[0] / ci_v[1], obs_ratio_hi=ci_a[1] / ci_v[0])
    for w in (0.0, 0.2, 0.4):
        f, v, rho = S(w, key, 'log_F', x), S(w, key, 'log_v', x), S(w, key, 'log_rho', x)
        row[f'w{w}_dlogF'] = f.coef; row[f'w{w}_dlogv'] = v.coef; row[f'w{w}_dlogrho'] = rho.coef
        row[f'w{w}_ratio'] = f.coef / v.coef
        # interval from the bootstrap ends (conservative, as for the observed ratio)
        row[f'w{w}_ratio_lo'] = min(f.lo / v.hi, f.hi / v.lo) if v.lo > 0 else np.nan
        row[f'w{w}_ratio_hi'] = max(f.lo / v.hi, f.hi / v.lo) if v.lo > 0 else np.nan
        row[f'w{w}_naive_ratio'] = (3 * v.coef + rho.coef) / v.coef   # if F_c/F ~ rho v^3 at fixed Teff
        # response exponents mapping model onto observed quantities along this axis:
        # alpha0 ~ (F_c/F)^beta_alpha, v_macro ~ v_c^beta_v
        row[f'w{w}_beta_alpha'] = oa / f.coef
        row[f'w{w}_beta_v'] = ov / v.coef
    ratio_rows.append(row)
ratios = pd.DataFrame(ratio_rows)
ratios.to_csv(f'{D}/mesa_ratios.csv', index=False, float_format='%.3f')
print('\n(3) SENSITIVITY RATIOS  d log(F_c/F) / d log v_c  (model)  vs  d log alpha0 / d log v_macro (observed)')
print(ratios[['axis', 'obs_ratio', 'obs_ratio_lo', 'obs_ratio_hi', 'w0.0_dlogF', 'w0.0_dlogv', 'w0.0_dlogrho',
              'w0.0_ratio', 'w0.0_ratio_lo', 'w0.0_ratio_hi', 'w0.0_naive_ratio', 'w0.0_beta_alpha', 'w0.0_beta_v', 'w0.2_ratio', 'w0.4_ratio']].to_string(index=False))

print('\n(5) EVOLUTION at fixed mass (MW w0.0, complete MS, log M in observed range):')
for y, o in (('log_v', 'vmac'), ('log_F', 'alpha0'), ('log_nuc', 'nuchar')):
    t = S(0.0, 'tauM', y, 'tau'); m = S(0.0, 'tauM', y, 'logM')
    print(f'  {y:8s} d/dtau {t.coef:+.2f} [{t.lo:+.2f},{t.hi:+.2f}]  d/dlogM {m.coef:+.2f}   '
          f'(obs {o}: tau {OBS_SLOPES[(o, "tau")]:+.2f}' + (f', logM {OBS_SLOPES[(o, "logM")]:+.2f})' if (o, 'logM') in OBS_SLOPES else ')'))

# ------------------------------------------------------------------ (4) nu_c vs nu_char
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

GT, GL = GP['logTeff_ext'], GP['logL_ext']
gp_nu = GP['nuchar_mean']            # log10 nu_char [1/d]
gp_ok = GP['nuchar_ok'].astype(bool)


def model_map(df, col, dx=0.03, dy=0.08):
    k = np.isfinite(df[col])
    x, y, z = df.logTeff[k].values, df.logLspec[k].values, df[col][k].values
    zi = griddata((x, y), z, (GT, GL), method='linear')
    dist, _ = cKDTree(np.c_[x / dx, y / dy]).query(np.c_[GT.ravel() / dx, GL.ravel() / dy])
    zi[dist.reshape(GT.shape) > 1] = np.nan
    return zi


nuc_rows = []
maps = {}
for w in (0.0, 0.2, 0.4):
    df = pd.concat([ms[(ms.Z == 'MW') & (ms.w == w)], post[(post.Z == 'MW') & (post.w == w)]])
    df = df[df.FeCZ_present]
    mm = model_map(df, 'log_nuc')
    maps[w] = mm
    dlt = gp_nu - mm
    for region, sel in (('OB (logT>4.3)', GT > 4.3), ('all', GT > 0)):
        k = gp_ok & np.isfinite(dlt) & sel
        c = wls(np.c_[GL[k], GT[k]], dlt[k])
        cm = wls(np.c_[GL[k], GT[k]], mm[k])
        cg = wls(np.c_[GL[k], GT[k]], gp_nu[k])
        nuc_rows.append(dict(w=w, region=region, n_cells=int(k.sum()), median_offset_dex=np.median(dlt[k]),
                             p16=np.percentile(dlt[k], 16), p84=np.percentile(dlt[k], 84),
                             offset_dlogL=c[1], offset_dlogT=c[2],
                             model_dlogL=cm[1], model_dlogT=cm[2], gp_dlogL=cg[1], gp_dlogT=cg[2]))
nuc = pd.DataFrame(nuc_rows)
nuc.to_csv(f'{D}/mesa_nuc_offset.csv', index=False, float_format='%.3f')
print('\n(4) nu_char(GP) vs nu_c(model FeCZ) on the trustworthy GP cells covered by models')
print(nuc.to_string(index=False))
t = S(0.0, 'L', 'log_nuc', 'logTeff'); l = S(0.0, 'L', 'log_nuc', 'logLspec')
print(f'  model rows (MS, sHRD domain): d log nu_c/d log Teff {t.coef:+.2f} [{t.lo:+.2f},{t.hi:+.2f}], '
      f'd/d log L {l.coef:+.2f} [{l.lo:+.2f},{l.hi:+.2f}]   (obs +1.53, -0.23)')

# ------------------------------------------------------------------ (6) cool side
cool_rows = []
allmw = pd.concat([ms[(ms.Z == 'MW') & (ms.w == 0.0)].assign(phase='MS'),
                   post[(post.Z == 'MW') & (post.w == 0.0)].assign(phase='postMS')])
REG = {'hot OB (logT>=4.3)': (4.3, 5.0, -9, 9), 'evolved OB/BSG (4.0-4.3, logL>3.6)': (4.0, 4.3, 3.6, 9),
       'YSG (3.75-4.0)': (3.75, 4.0, -9, 9)}
for reg, (t0, t1, l0, l1) in REG.items():
    sub = allmw[allmw.logTeff.between(t0, t1) & allmw.logLspec.between(l0, l1)]
    wt = sub.w_time.values
    for zn in ('FeCZ', 'HeII', 'HeI', 'HI'):
        pres = sub[f'{zn}_present'].values.astype(bool)
        row = dict(region=reg, zone=zn, n_rows=len(sub), time_frac_present=np.sum(wt[pres]) / np.sum(wt))
        for q in ('nuc_d', 'vmax_kms', 'FcF'):
            v = sub[f'{zn}_{q}'].values[pres]
            if len(v):
                o = np.argsort(v); cw = np.cumsum(wt[pres][o]) / np.sum(wt[pres])
                row[f'{q}_med'] = v[o][np.searchsorted(cw, 0.5)]
                row[f'{q}_p16'] = v[o][np.searchsorted(cw, 0.16)]
                row[f'{q}_p84'] = v[o][min(np.searchsorted(cw, 0.84), len(v) - 1)]
        cool_rows.append(row)
cool = pd.DataFrame(cool_rows)
cool.to_csv(f'{D}/mesa_coolside.csv', index=False, float_format='%.4g')
print('\n(6) COOL SIDE (MW w0.0, MS+post-MS, time-weighted): observed nu_char hot 1.47, BSG 0.27, YSG 0.30 1/d')
print(cool[['region', 'zone', 'time_frac_present', 'nuc_d_med', 'nuc_d_p16', 'nuc_d_p84', 'vmax_kms_med', 'FcF_med']]
      .to_string(index=False))

# ------------------------------------------------------------------ observed binned medians
eb = np.arange(2.5, 4.61, 0.15)
ob_v = binned(vm_hot.logL_sp.values, np.log10(vm_hot.vmac.values), eb)
ob_a = binned(a0_ob.lL.values, a0_ob.log_alpha0.values, eb, nmin=4)
pd.concat([pd.DataFrame(ob_v, columns=['logL', 'med', 'p16', 'p84', 'n']).assign(q='log vmac (logT>4.3)'),
           pd.DataFrame(ob_a, columns=['logL', 'med', 'p16', 'p84', 'n']).assign(q='log alpha0 (OB)')]
          ).to_csv(f'{D}/obs_binned.csv', index=False, float_format='%.3f')

# save what the figure needs
np.savez(f'{D}/fig8_inputs.npz', ob_v=ob_v, ob_a=ob_a, GT=GT, GL=GL, gp_nu=gp_nu, gp_ok=gp_ok,
         model_nu_w0=maps[0.0], model_nu_w2=maps[0.2])
print('\nwritten:', ', '.join(sorted(os.listdir(D))))
