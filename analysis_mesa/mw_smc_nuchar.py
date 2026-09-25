#!/usr/bin/env python3
"""
Homogeneous Galactic vs SMC nu_char: the same pipeline (vandaele_fit.py: 0.1-40 d-1 amplitude spectrum, least-squares
semi-Lorentzian, automatic prewhitening, nu_char >= 0.15 and Delta BIC > 0) on
  - SMC: Van Daele+2026 PSF light curves (data_obs/vandaele_fits.csv), placed with Bestenlehner+2025 T_eff and log g;
  - MW: SPOC 2-min light curves of the Bowman 2020 + Shen 2024 stars (tess_mw_download.py; data_obs/vandaele_fits_mw*.csv),
    placed with our catalogue (rednoise_sHRD_extended_evol.csv).
(1) validation of our MW fits against the published values: Shen 2024 per sector (same sector), Bowman 2020 per star;
(2) per-star median log nu_char; plane log nu = c0 + c1 log L_spec + c2 log Teff fitted on MW; SMC residuals inside the MW
    domain; bootstrap over stars. Variants: fit metric (lin/log), MW noise-degraded to the SMC alpha0/C_w (tag _mwdeg).
Writes data_obs/mw_smc_nuchar.csv and data_obs/mw_smc_validation.csv.
"""
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/skill')
from kernel import vz_read_tsv, star_key          # noqa: E402
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
RAW = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vizier'
LSUN = np.log10(5777.0 ** 4 / 27400.0)
rng = np.random.default_rng(5)


def safe(n):
    return re.sub(r'[^A-Za-z0-9+-]', '_', n.strip())


def good(f, metric):
    return f[(f.metric == metric) & f.prewhitened & (f.nuchar >= 0.15) & (f.dBIC_vs_white > 0)]


def names():
    b = pd.read_csv(f'{OBS}/bowman2020_A36_tablea2.tsv', sep='\t', comment='#')
    b = b[pd.to_numeric(b.TIC, errors='coerce').notna()]
    s = vz_read_tsv(f'{RAW}/shen24_t1.tsv')
    nm = list(b.Name.str.strip()) + list(s.Name.str.strip())
    return {safe(n): n for n in nm}, b


def validation(fmw, metric):
    m, b = names()
    f = good(fmw, metric).copy()
    f['Name'] = f.id.map(m)
    s2 = vz_read_tsv(f'{RAW}/shen24_t2.tsv')
    s2['Name'] = s2.Name.str.strip(); s2['sector'] = pd.to_numeric(s2.Sec, errors='coerce')
    s2['nu_pub'] = pd.to_numeric(s2['nu-char'], errors='coerce')
    x = f.merge(s2[['Name', 'sector', 'nu_pub']], on=['Name', 'sector'])
    rows = []
    d = np.log10(x.nuchar) - np.log10(x.nu_pub)
    rows.append(dict(ref='Shen2024 same sector', metric=metric, n=len(x), median=d.median(),
                     mad=np.median(abs(d - d.median())), within_0p1=(abs(d) < 0.1).mean(),
                     spearman=spearmanr(x.nuchar, x.nu_pub).correlation))
    b['Name'] = b.Name.str.strip(); b['nu_pub'] = b.nuchar.astype(float)
    y = f[f.sector <= 13].groupby('Name').nuchar.apply(lambda v: 10 ** np.median(np.log10(v))).reset_index()
    y = y.merge(b[['Name', 'nu_pub']], on='Name')
    d = np.log10(y.nuchar) - np.log10(y.nu_pub)
    rows.append(dict(ref='Bowman2020 per star (sectors <= 13)', metric=metric, n=len(y), median=d.median(),
                     mad=np.median(abs(d - d.median())), within_0p1=(abs(d) < 0.1).mean(),
                     spearman=spearmanr(y.nuchar, y.nu_pub).correlation))
    return rows


def per_star(f, metric):
    g = good(f, metric)
    return g.groupby('id').agg(lnu=('nuchar', lambda v: np.median(np.log10(v))), nsec=('nuchar', 'size'),
                               la=('alpha0', lambda v: np.nan), a0cw=('alpha0', 'size')).reset_index().drop(columns=['la'])


def place_mw(st):
    m, _ = names()
    E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
    E = E[E['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic'])].copy()
    E['key'] = E.star.map(star_key)
    st['key'] = st.id.map(m).map(lambda n: star_key(n) if isinstance(n, str) else None)
    return st.merge(E.drop_duplicates('key')[['key', 'lT', 'lL']], on='key', how='inner')


def place_smc(st):
    sys.path.insert(0, HERE)
    import vandaele_compare as vc
    g = vc.placed()
    g['id'] = 'GAIA DR3 ' + g.gaia.astype(str)
    g['lL_spec'] = 4 * g.lT - pd.to_numeric(g.logg, errors='coerce') - LSUN
    x = st.merge(g[['id', 'lT', 'lL_spec', 'lL']].rename(columns={'lL': 'lL_class'}), on='id')
    return x.rename(columns={'lL_spec': 'lL'})


def compare(mw, smc, label):
    mw = mw[np.isfinite(mw.lT) & np.isfinite(mw.lL)]
    smc = smc[np.isfinite(smc.lT) & np.isfinite(smc.lL)]
    A = np.c_[np.ones(len(mw)), mw.lL, mw.lT]
    c = np.linalg.lstsq(A, mw.lnu, rcond=None)[0]
    rmw = mw.lnu - A @ c
    lo_L, hi_L = np.percentile(mw.lL, [2, 98]); lo_T, hi_T = np.percentile(mw.lT, [2, 98])
    s = smc[smc.lL.between(lo_L, hi_L) & smc.lT.between(lo_T, hi_T)]
    B = np.c_[np.ones(len(s)), s.lL, s.lT]
    rs = s.lnu.values - B @ c
    boot = []
    for _ in range(2000):
        i = rng.integers(0, len(mw), len(mw)); j = rng.integers(0, len(s), len(s))
        cc = np.linalg.lstsq(A[i], mw.lnu.values[i], rcond=None)[0]
        boot.append(np.median(s.lnu.values[j] - B[j] @ cc))
    hot = lambda d: d[d.lT.between(4.4, 4.65)]
    return dict(variant=label, n_mw=len(mw), n_smc_domain=len(s), n_smc_all=len(smc),
                plane_c1_lL=c[1], plane_c2_lT=c[2], rms_mw=np.std(rmw),
                smc_median_residual=np.median(rs), lo95=np.percentile(boot, 2.5), hi95=np.percentile(boot, 97.5),
                mwu_p=mannwhitneyu(rs, rmw).pvalue,
                median_lnu_mw_hot=hot(mw).lnu.median(), median_lnu_smc_hot=hot(smc).lnu.median(),
                n_hot_mw=len(hot(mw)), n_hot_smc=len(hot(smc)),
                smc_lL_range=f'{s.lL.min():.2f}-{s.lL.max():.2f}' if len(s) else '')


def main():
    fsmc = pd.read_csv(f'{HERE}/data_obs/vandaele_fits.csv')
    out, val = [], []
    for tag in ('_mw', '_mwdeg', '_mwsap'):
        p = f'{HERE}/data_obs/vandaele_fits{tag}.csv'
        if not os.path.exists(p):
            continue
        fmw = pd.read_csv(p)
        for metric in ('lin', 'log'):
            if tag == '_mw':
                val += validation(fmw, metric)
            mw = place_mw(per_star(fmw, metric))
            smc = place_smc(per_star(fsmc, metric))
            out.append(compare(mw, smc, f'{tag[1:]} / {metric}'))
    pd.set_option('display.width', 250)
    v = pd.DataFrame(val); r = pd.DataFrame(out)
    print(v.to_string(index=False, float_format=lambda x: f'{x:+.3f}'))
    print(r.T.to_string())
    v.to_csv(f'{HERE}/data_obs/mw_smc_validation.csv', index=False, float_format='%.3f')
    r.to_csv(f'{HERE}/data_obs/mw_smc_nuchar.csv', index=False, float_format='%.3f')


if __name__ == '__main__':
    main()
