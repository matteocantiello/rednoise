#!/usr/bin/env python3
"""
Analysis of smc_floor_model.py outputs.
(1) injection-recovery into the real SMC noise floor: detection fraction above the floor, and the bias and scatter of the
    recovered log nu_char for the naive fit ('slf') and the floor-aware fit ('slf+floor'), vs nu_in and R = a0_in / A_floor;
(2) SMC star-sectors: detection above the floor by Van Daele class, and the estimated R;
(3) MW vs SMC nu_char at matched (L, Teff), all with Whittle fits: (a) naive, all SMC SLF-class sectors; (b) floor-aware,
    only sectors with SLF detected above the floor; (c) as (b) restricted to the R range where the injections show < 0.1 dex bias.
Writes data_obs/floor_injection_summary.csv and data_obs/floor_mw_smc.csv.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mw_smc_nuchar as ms          # noqa: E402
import vandaele_compare as vc       # noqa: E402
DET = 10.0
rng = np.random.default_rng(9)


def injection():
    d = pd.concat([pd.read_csv(p) for p in sorted(__import__('glob').glob(f'{HERE}/data_obs/floor_injection*.csv'))
                   if 'summary' not in p], ignore_index=True)
    d['det'] = d.dbic_above_floor > DET
    d['b_naive'] = np.log10(d['nu_slf'] / d.nu_in)
    d['b_floor'] = np.log10(d['nu_slf+floor'] / d.nu_in)
    mad = lambda v: np.median(abs(v - np.median(v))) if len(v) else np.nan
    s = d.groupby(['nu_in', 'R']).apply(lambda x: pd.Series(dict(
        n=len(x), det_frac=x.det.mean(), bias_naive=x.b_naive.median(), mad_naive=mad(x.b_naive),
        bias_floor_det=x.b_floor[x.det].median(), mad_floor_det=mad(x.b_floor[x.det].values),
        frac_within_0p1_floor_det=(abs(x.b_floor[x.det]) < 0.1).mean() if x.det.any() else np.nan))).reset_index()
    return d, s


def per_star(df, nucol, sel):
    g = df[sel]
    return g.groupby('id').agg(lnu=(nucol, lambda v: np.median(np.log10(v))), nsec=(nucol, 'size')).reset_index()


def compare(mw, smc, label):
    mw = mw[np.isfinite(mw.lT) & np.isfinite(mw.lL)]
    smc = smc[np.isfinite(smc.lT) & np.isfinite(smc.lL)]
    A = np.c_[np.ones(len(mw)), mw.lL, mw.lT]
    c = np.linalg.lstsq(A, mw.lnu, rcond=None)[0]
    rmw = mw.lnu - A @ c
    lo_L, hi_L = np.percentile(mw.lL, [2, 98]); lo_T, hi_T = np.percentile(mw.lT, [2, 98])
    s = smc[smc.lL.between(lo_L, hi_L) & smc.lT.between(lo_T, hi_T)]
    if len(s) < 3:
        return dict(variant=label, n_mw=len(mw), n_smc_domain=len(s))
    B = np.c_[np.ones(len(s)), s.lL, s.lT]
    rs = s.lnu.values - B @ c
    boot = []
    for _ in range(2000):
        i = rng.integers(0, len(mw), len(mw)); j = rng.integers(0, len(s), len(s))
        cc = np.linalg.lstsq(A[i], mw.lnu.values[i], rcond=None)[0]
        boot.append(np.median(s.lnu.values[j] - B[j] @ cc))
    return dict(variant=label, n_mw=len(mw), n_smc_domain=len(s), smc_median_residual=np.median(rs),
                lo95=np.percentile(boot, 2.5), hi95=np.percentile(boot, 97.5), mwu_p=mannwhitneyu(rs, rmw).pvalue,
                rms_mw=np.std(rmw))


def main():
    inj, s = injection()
    s.to_csv(f'{HERE}/data_obs/floor_injection_summary.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 220)
    print(s.to_string(index=False, float_format=lambda v: f'{v:+.2f}'))
    # R range where the floor-aware fit is unbiased (< 0.1 dex) and detected in > 50% at every nu_in
    good_R = [R for R, g in s.groupby('R') if (g.det_frac > 0.5).all() and (abs(g.bias_floor_det) < 0.1).all()]
    Rmin = min(good_R) if good_R else np.inf
    print(f'R threshold (all nu_in detected > 50% and |bias| < 0.1 dex): R >= {Rmin}')

    smc = pd.read_csv(f'{HERE}/data_obs/floor_fits_smc.csv')
    smc['id'] = 'GAIA DR3 ' + smc.gaia.astype(str)
    smc['det'] = smc.dbic_above_floor > DET
    smc['R_est'] = smc['a0_slf+floor'] / (smc.Af * 10 ** smc['logk_slf+floor'])
    print(smc.groupby('cls').agg(n=('det', 'size'), det=('det', 'mean'), R_med_det=('R_est', lambda v: v[smc.loc[v.index, 'det']].median()),
                                 lnu_naive=('nu_slf', lambda v: np.log10(v).median()),
                                 lnu_floor_det=('nu_slf+floor', lambda v: np.log10(v[smc.loc[v.index, 'det']]).median())).round(2).to_string())
    mwf = pd.read_csv(f'{HERE}/data_obs/floor_fits_mw.csv')
    mw = ms.place_mw(per_star(mwf, 'nu_slf', (mwf.dbic_slf_vs_white > DET) & (mwf.nu_slf >= 0.15)))
    slf = smc.cls == 'SLF'
    variants = {
        'naive Whittle, SLF-class sectors (nu >= 0.15)': per_star(smc, 'nu_slf', slf & (smc.nu_slf >= 0.15)),
        'naive Whittle, all classes (nu >= 0.15)': per_star(smc, 'nu_slf', smc.nu_slf >= 0.15),
        'floor-aware, detected above floor (any class)': per_star(smc, 'nu_slf+floor', smc.det),
        'floor-aware, detected above floor, SLF class': per_star(smc, 'nu_slf+floor', smc.det & slf),
        f'floor-aware, detected and R >= {Rmin}': per_star(smc, 'nu_slf+floor', smc.det & (smc.R_est >= Rmin)),
    }
    rows = []
    for lab, st in variants.items():
        rows.append(compare(mw, ms.place_smc(st), lab))
    r = pd.DataFrame(rows)
    r.to_csv(f'{HERE}/data_obs/floor_mw_smc.csv', index=False, float_format='%.3f')
    print(r.to_string(index=False, float_format=lambda v: f'{v:+.3f}'))


if __name__ == '__main__':
    main()


def expected_R():
    """R a Galactic-like SLF would have in each SMC star-sector: MW plane in log(a0/flux) (Whittle fits, PDCSAP median flux)
    evaluated at the SMC star, divided by the SMC floor prediction for that sector and Tmag."""
    import glob
    import smc_floor_model as fm
    mwf = pd.read_csv(f'{HERE}/data_obs/floor_fits_mw.csv')
    fl = {}
    for p in glob.glob(f'{fm.LC_MW}/*_SPOC.txt'):
        b = os.path.basename(p)
        fl[b] = np.nanmedian(np.loadtxt(p, usecols=1))
    mwf['flux'] = [fl.get(f'{i}_sector{s:02d}_SPOC.txt', np.nan) for i, s in zip(mwf.id, mwf.sector)]
    g = mwf[(mwf.dbic_slf_vs_white > DET) & (mwf.nu_slf >= 0.15) & (mwf.flux > 0)].copy()
    g['la_rel'] = np.log10(g.a0_slf / g.flux)
    st = g.groupby('id').agg(lnu=('la_rel', 'median')).reset_index()      # 'lnu' column reused by place_mw
    mw = ms.place_mw(st); mw = mw[np.isfinite(mw.lT) & np.isfinite(mw.lL)]
    A = np.c_[np.ones(len(mw)), mw.lL, mw.lT]; c = np.linalg.lstsq(A, mw.lnu, rcond=None)[0]
    sd = np.std(mw.lnu - A @ c)
    cal = pd.read_csv(f'{HERE}/data_obs/floor_calib.csv')
    F, sig, beta, shape = fm.floor_model(cal)
    smc = pd.read_csv(f'{HERE}/data_obs/floor_fits_smc.csv'); smc['id'] = 'GAIA DR3 ' + smc.gaia.astype(str)
    pl = ms.place_smc(smc[['id']].drop_duplicates().assign(lnu=0.0, nsec=1, a0cw=1))[['id', 'lT', 'lL']]
    x = smc.merge(pl, on='id')
    x = x[np.isfinite(x.lT) & np.isfinite(x.lL)]
    lo_L, hi_L = np.percentile(mw.lL, [2, 98]); lo_T, hi_T = np.percentile(mw.lT, [2, 98])
    x = x[x.lL.between(lo_L, hi_L) & x.lT.between(lo_T, hi_T)]
    flux_rel_floor = np.array([F(s, t, 1.0)[0] for s, t in zip(x.sector, x.Tmag)])
    x['logR_exp'] = (np.c_[np.ones(len(x)), x.lL, x.lT] @ c) - np.log10(flux_rel_floor)
    print(f'MW plane in log(a0/flux): N = {len(mw)}, rms {sd:.2f} dex; median fractional a0 at SMC positions '
          f'{np.median(np.c_[np.ones(len(x)), x.lL, x.lT] @ c):+.2f} (log); SMC floor fractional median {np.log10(np.median(flux_rel_floor)):+.2f}')
    print(f'expected log R for a Galactic-like SLF in the {len(x)} SMC star-sectors in the domain: median {x.logR_exp.median():+.2f}, '
          f'16-84% [{x.logR_exp.quantile(0.16):+.2f}, {x.logR_exp.quantile(0.84):+.2f}]; fraction with R_exp >= 1: {(x.logR_exp >= 0).mean():.2f}')
    print(f'observed detection above floor in these sectors: {(x.dbic_above_floor > DET).mean():.2f} '
          f'(SLF class {(x[x.cls == "SLF"].dbic_above_floor > DET).mean():.2f}, N = {(x.cls == "SLF").sum()})')
    x.to_csv(f'{HERE}/data_obs/floor_expected_R.csv', index=False, float_format='%.4g')
    return x


if __name__ == '__main__' and os.environ.get('EXPECTED_R'):
    expected_R()


def expected_detection():
    """Fold each SMC sector's expected (nu, R) for Galactic-like SLF through the injection detection map, and compare with the
    observed detection rate above the floor. Also the same for SLF weaker or slower by given factors."""
    from scipy.interpolate import RegularGridInterpolator
    x = pd.read_csv(f'{HERE}/data_obs/floor_expected_R.csv')
    s = pd.read_csv(f'{HERE}/data_obs/floor_injection_summary.csv')
    nus, Rs = sorted(s.nu_in.unique()), sorted(s.R.unique())
    D = s.pivot(index='nu_in', columns='R', values='det_frac').loc[nus, Rs].values
    Bn = s.pivot(index='nu_in', columns='R', values='bias_naive').loc[nus, Rs].values
    f = RegularGridInterpolator((np.log10(nus), np.log10(Rs)), D, bounds_error=False, fill_value=None)
    fb = RegularGridInterpolator((np.log10(nus), np.log10(Rs)), Bn, bounds_error=False, fill_value=None)
    mwf = pd.read_csv(f'{HERE}/data_obs/floor_fits_mw.csv')
    mw = ms.place_mw(per_star(mwf, 'nu_slf', (mwf.dbic_slf_vs_white > DET) & (mwf.nu_slf >= 0.15)))
    mw = mw[np.isfinite(mw.lT) & np.isfinite(mw.lL)]
    A = np.c_[np.ones(len(mw)), mw.lL, mw.lT]; c = np.linalg.lstsq(A, mw.lnu, rcond=None)[0]
    x['lnu_exp'] = np.c_[np.ones(len(x)), x.lL, x.lT] @ c
    obs = (x.dbic_above_floor > DET).mean()
    rows = []
    for dnu in (0.0, -0.2, -0.4, -0.6):
        for dR in (0.0, -0.3, -0.6, -1.0):
            ln = np.clip(x.lnu_exp + dnu, np.log10(nus[0]), np.log10(nus[-1]))
            lr = np.clip(x.logR_exp + dR, np.log10(Rs[0]) - 0.3, np.log10(Rs[-1]))
            p = np.clip(f(np.c_[ln, lr]), 0, 1)
            # naive-fit median log nu the pipeline would report
            rep = x.lnu_exp + dnu + fb(np.c_[ln, np.clip(lr, np.log10(Rs[0]), np.log10(Rs[-1]))])
            rows.append(dict(shift_log_nu=dnu, shift_log_amp=dR, expected_det=p.mean(), observed_det=obs,
                             median_lnu_true=np.median(x.lnu_exp + dnu), median_lnu_naive_expected=np.median(rep)))
    r = pd.DataFrame(rows)
    print(f'Galactic plane nu at the SMC positions: median log nu {x.lnu_exp.median():+.2f}; observed naive median log nu '
          f'(Whittle) {np.log10(x.nu_slf).median():+.2f}; observed detection above floor {obs:.2f} (N = {len(x)})')
    print(r.to_string(index=False, float_format=lambda v: f'{v:+.2f}'))
    r.to_csv(f'{HERE}/data_obs/floor_expected_detection.csv', index=False, float_format='%.3f')
    return r


if __name__ == '__main__' and os.environ.get('EXPECTED_DET'):
    expected_detection()
