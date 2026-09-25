#!/usr/bin/env python3
"""
Magnetic suppression test: do magnetic hot stars (Shen+2023, J/ApJ/955/123; same SLF method and units as Shen 2024) have lower
red-noise amplitude / different nu_char than non-magnetic stars at the same (log L_spec, log Teff)?

- Shen 2023: per-star medians over sectors (as for Shen 2024 in the primary sample); log L_spec = 4 log Teff - log g - log ell_sun.
- Dipole strengths B_d from Shultz+2022 (MOBSTER VI, J/MNRAS/513/1429, table A1).
- Baseline: the primary sample (Bowman 2020 + Shen 2024, non-magnetic), fitted with log y = c0 + c1 log L_spec + c2 log Teff;
  magnetic-star residuals from that plane, within the primary sample's domain.
- FeCZ shut-off field at each star's position from the MW v2 grid (history FeCZ_B_shutoff_conv, MS models): the FeCZ picture
  predicts suppression where the star has an FeCZ and B_d > B_shutoff.
Writes data_obs/magnetic_test.csv and data_obs/magnetic_stars.csv.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
D2 = f'{HERE}/data_v2'
sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/skill')
from kernel import vz_fetch, star_key          # noqa: E402
RAW = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vizier'
LSUN = np.log10(5777.0 ** 4 / 27400.0)


def main():
    t1 = vz_fetch('J/ApJ/955/123/table1', 'shen23_t1', outdir=RAW)
    t2 = vz_fetch('J/ApJ/955/123/table2', 'shen23_t2', outdir=RAW)
    sh = vz_fetch('J/MNRAS/513/1429/tablea1', 'shultz22_a1', outdir=RAW)
    print('Shen23 t1', len(t1), list(t1.columns)[:12]); print('Shen23 t2', len(t2), list(t2.columns)[:12])
    key2 = 'Name' if 'Name' in t2 else t2.columns[0]
    agg = t2.groupby(key2).agg(alpha0=('alpha0', 'median'), nuchar=('nuChar', 'median'), gamma=('gamma', 'median'),
                               nsec=('alpha0', 'size')).reset_index().rename(columns={key2: 'Name'})
    m = t1.merge(agg, on='Name', how='inner')
    m['lT'] = pd.to_numeric(m.logT, errors='coerce')
    m['lg'] = pd.to_numeric(m.logg, errors='coerce')
    m['lL'] = 4 * m.lT - m.lg - LSUN
    m['key'] = m.Name.map(star_key)
    sh['key'] = sh.Star.map(star_key)
    sh['Bd'] = pd.to_numeric(sh.Bd, errors='coerce')
    sh['l_Bd'] = sh.l_Bd.fillna('').astype(str).str.strip()
    m = m.merge(sh[['key', 'Bd', 'l_Bd']], on='key', how='left')
    m['Bd_src'] = np.where(np.isfinite(m.Bd), 'Shultz2022', '')
    # fall back on Petit+2013 (MNRAS 429, 398) Table 1 polar field B_p (= B_d for a dipole; parsed from the arXiv source,
    # handoff/data_survey_raw/petit2013); '>' entries are lower limits
    pt = pd.read_csv(f'{RAW}/../petit2013/petit2013_table1_Bp.csv')
    pt['key'] = pt.Star.str.replace('NGC 1624-2', 'NGC 1624 2').map(star_key)
    m = m.merge(pt[['key', 'Bp_kG', 'l_Bp']], on='key', how='left')
    fb = ~np.isfinite(m.Bd) & np.isfinite(m.Bp_kG)
    m.loc[fb, 'Bd'] = m.Bp_kG[fb]
    m.loc[fb, 'l_Bd'] = np.where(m.l_Bp[fb].fillna('') == '>', 'LL', '')
    m.loc[fb, 'Bd_src'] = 'Petit2013'
    m['l_Bd'] = m.l_Bd.fillna('')
    print(f'Shen23 stars with SLF medians: {len(m)}; with log g: {np.isfinite(m.lg).sum()}; with Shultz B_d: {np.isfinite(m.Bd).sum()}')

    E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
    prim = E[E['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic']) & E.a0_ok]
    prim_all = prim
    prim_keys = set(prim.star.map(star_key))
    m['in_primary'] = m.key.isin(prim_keys)
    # stars fitted by both Shen 2023 and the primary sample: same-method consistency check
    ov = m[m.in_primary].merge(prim_all.assign(key=prim.star.map(star_key))[['key', 'alpha0', 'nuchar']], on='key',
                               suffixes=('', '_prim'))
    if len(ov):
        d_a = np.log10(ov.alpha0) - np.log10(ov.alpha0_prim); d_n = np.log10(ov.nuchar) - np.log10(ov.nuchar_prim)
        print(f'overlap with primary sample: N = {len(ov)}; log alpha0 (Shen23 - primary) median {np.median(d_a):+.2f}, '
              f'MAD {np.median(abs(d_a - np.median(d_a))):.2f}; log nu_char {np.nanmedian(d_n):+.2f}')
        print(ov[['Name', 'alpha0', 'alpha0_prim', 'nuchar', 'nuchar_prim']].to_string(index=False))
    # the baseline must be non-magnetic: drop every known magnetic star (Shen 2023 sample, Petit 2013, Shultz 2022)
    mag_keys = set(m.key) | set(pt.key) | set(sh.key)
    n0 = len(prim)
    prim = prim[~prim.star.map(star_key).isin(mag_keys)]
    print(f'primary baseline: {n0} -> {len(prim)} after removing known magnetic stars')
    lo_L, hi_L = np.percentile(prim.lL, [2, 98]); lo_T, hi_T = np.percentile(prim.lT, [2, 98])
    dom = m[m.lL.between(lo_L, hi_L) & m.lT.between(lo_T, hi_T) & np.isfinite(m.alpha0)].copy()
    print(f'magnetic stars in the primary domain (log L {lo_L:.2f}-{hi_L:.2f}, log Teff {lo_T:.2f}-{hi_T:.2f}), '
          f'incl. those also in the primary sample: {len(dom)}')

    # FeCZ shut-off field and presence at each star's position (MW v2 grid, non-rotating, main sequence)
    os.environ.setdefault('RN_GRID', '/mnt/home/mcantiello/work/rednoise/models/grids_v2')
    sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/models/grids')
    import grid_io as g
    X, Y, BS, PR = [], [], [], []
    for mass, d in g.load_subgrid('MW', 'w0.0', verbose=False).items():
        h = d['center_h1']; i0 = g.find_zams(h); it = g.find_tams(h, 1e-3)
        sl = slice(i0, it if it is not None else len(h))
        lt = d['log_Teff'][sl]; ll = g.spec_ell(lt, d['log_g'][sl])
        bs = np.asarray(d['FeCZ_B_shutoff_conv'][sl], float)
        X.append(lt); Y.append(ll); BS.append(bs); PR.append(np.isfinite(bs) & (bs > 0))
    X, Y, BS, PR = map(np.concatenate, (X, Y, BS, PR))
    step = max(1, len(X) // 40000)
    X, Y, BS, PR = X[::step], Y[::step], BS[::step], PR[::step]
    pts = np.c_[X / 0.03, Y / 0.08]
    pres = LinearNDInterpolator(pts, PR.astype(float))
    ok = PR
    lbs = LinearNDInterpolator(pts[ok], np.log10(BS[ok]))
    P = np.c_[dom.lT / 0.03, dom.lL / 0.08]
    dom['FeCZ_present_model'] = pres(P)
    dom['B_shutoff_kG'] = 10 ** lbs(P) / 1e3
    dom.loc[dom.FeCZ_present_model < 0.5, 'B_shutoff_kG'] = np.nan
    dom['Bd_over_Bshutoff'] = dom.Bd / dom.B_shutoff_kG

    # baseline planes from the (non-magnetic) primary sample, and magnetic residuals
    rows = []
    for q, col_p, col_m in (('log alpha0', 'alpha0', 'alpha0'), ('log nu_char', 'nuchar', 'nuchar')):
        pp = prim[np.isfinite(prim[col_p]) & (prim[col_p] > 0)]
        A = np.c_[np.ones(len(pp)), pp.lL, pp.lT]
        c = np.linalg.lstsq(A, np.log10(pp[col_p]), rcond=None)[0]
        rp = np.log10(pp[col_p]) - A @ c
        mm = dom[np.isfinite(dom[col_m]) & (dom[col_m] > 0)]
        rm = np.log10(mm[col_m]) - np.c_[np.ones(len(mm)), mm.lL, mm.lT] @ c
        dom.loc[mm.index, f'res_{col_m}'] = rm
        mm = dom.loc[mm.index]
        base = dict(quantity=q, n_primary=len(pp), n_magnetic=len(mm), median_res_magnetic=np.median(rm),
                    iqr_res_magnetic=np.subtract(*np.percentile(rm, [75, 25])), rms_res_primary=np.std(rp),
                    mwu_p=mannwhitneyu(rm, rp).pvalue)
        for lab, sel in (('with model FeCZ', mm.FeCZ_present_model >= 0.5), ('without model FeCZ', mm.FeCZ_present_model < 0.5)):
            base[f'median_res_{lab}'] = np.median(rm[sel]) if sel.sum() else np.nan
            base[f'n_{lab}'] = int(sel.sum())
        rot = mm.Var.astype(str).str.startswith('rot')
        base['median_res_rot_modulated'] = np.median(rm[rot.values]); base['n_rot_modulated'] = int(rot.sum())
        base['median_res_not_rot'] = np.median(rm[~rot.values]); base['n_not_rot'] = int((~rot).sum())
        ha = mm[np.isfinite(mm.Bd)]
        if len(ha) >= 6:
            base['spearman_res_vs_logBd_incl_limits'] = spearmanr(np.log10(ha.Bd), ha[f'res_{col_m}']).correlation
            base['spearman_p_incl_limits'] = spearmanr(np.log10(ha.Bd), ha[f'res_{col_m}']).pvalue
            base['n_Bd_incl_limits'] = len(ha)
            hi = ha.Bd_over_Bshutoff > 1
            base['median_res_Bd>Bshutoff_incl_limits'] = np.median(ha[f'res_{col_m}'][hi]) if hi.any() else np.nan
            base['median_res_Bd<Bshutoff_incl_limits'] = np.median(ha[f'res_{col_m}'][~hi]) if (~hi).any() else np.nan
            base['n_Bd>Bshutoff_incl_limits'] = int(hi.sum())
        hb = mm[np.isfinite(mm.Bd) & (mm.l_Bd == '')]
        if len(hb) >= 6:
            base['spearman_res_vs_logBd'] = spearmanr(np.log10(hb.Bd.clip(lower=1e-3)), hb[f'res_{col_m}']).correlation
            base['n_Bd'] = len(hb)
            hf = hb[np.isfinite(hb.Bd_over_Bshutoff)]
            if len(hf) >= 5:
                base['spearman_res_vs_Bd_over_Bshutoff'] = spearmanr(hf.Bd_over_Bshutoff, hf[f'res_{col_m}']).correlation
                base['n_Bd_FeCZ'] = len(hf)
                base['median_res_Bd>Bshutoff'] = np.median(hf[f'res_{col_m}'][hf.Bd_over_Bshutoff > 1]) if (hf.Bd_over_Bshutoff > 1).any() else np.nan
                base['n_Bd>Bshutoff'] = int((hf.Bd_over_Bshutoff > 1).sum())
        rows.append(base)
    res = pd.DataFrame(rows)
    out = f'{HERE}/data_obs'
    res.to_csv(f'{out}/magnetic_test.csv', index=False, float_format='%.3f')
    dom.to_csv(f'{out}/magnetic_stars.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    print(res.T.to_string())
    print(dom[['Name', 'SpT', 'lT', 'lL', 'alpha0', 'nuchar', 'Bd', 'l_Bd', 'Bd_src', 'Var', 'B_shutoff_kG', 'res_alpha0', 'res_nuchar']]
          .sort_values('lL').to_string(index=False, float_format=lambda v: f'{v:.2f}'))


if __name__ == '__main__':
    main()
