#!/usr/bin/env python3
"""
nu_char metallicity test, unblocked by the Bowman & Dorn-Wallenstein 2022 GP <-> Lorentzian cross-calibration
(digitize_bdw22.py -> data_obs/bdw22_nuchar.csv; 30 Galactic O stars fitted with both methods).

Galactic Lorentzian nu_char (Bowman 2020 + Shen 2024, a0_ok) against the Bowman+2024 LMC and SMC GP values, after mapping GP
onto the Lorentzian scale in four ways: (a) none; (b) median offset; (c) inverse of the linear fit; (d) running median of
log(Lor/GP) as a function of log GP (nearest 10 stars). Residuals from the Galactic plane log nu = c0 + c1 log L + c2 log Teff,
restricted to the Galactic domain; bootstrap over the 30 calibration stars propagates the conversion uncertainty.
Also: the matched-Teff comparison of the paper (4.4 < log Teff < 4.65). Writes data_obs/nuchar_Z_test.csv.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
rng = np.random.default_rng(3)


def conversions(cal):
    lg, lb = np.log10(cal.nuchar_GP.values), np.log10(cal.nuchar_B20_tab.values)
    off = np.median(lb - lg)
    sl, ic = np.polyfit(lb, lg, 1)
    o = np.argsort(lg)

    def run(x):
        x = np.atleast_1d(x)
        out = np.empty_like(x, dtype=float)
        for i, v in enumerate(x):
            k = np.argsort(abs(lg - v))[:10]
            out[i] = v + np.median(lb[k] - lg[k])
        return out
    return {'none': lambda x: x, 'median offset': lambda x: x + off,
            'linear fit (inverse)': lambda x: (x - ic) / sl, 'running median (10 nearest)': run}


def main():
    E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
    mw = E[E['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic']) & np.isfinite(E.nuchar) & (E.nuchar > 0)]
    mc = E[E['sample'].isin(['Bowman2024_LMC', 'Bowman2024_SMC']) & np.isfinite(E.nuchar) & (E.nuchar > 0)].copy()
    cal = pd.read_csv(f'{HERE}/data_obs/bdw22_nuchar.csv').dropna(subset=['nuchar_GP', 'nuchar_B20_tab'])
    A = np.c_[np.ones(len(mw)), mw.lL, mw.lT]
    c = np.linalg.lstsq(A, np.log10(mw.nuchar), rcond=None)[0]
    rmw = np.log10(mw.nuchar) - A @ c
    lo_L, hi_L = np.percentile(mw.lL, [2, 98]); lo_T, hi_T = np.percentile(mw.lT, [2, 98])
    dom = mc[mc.lL.between(lo_L, hi_L) & mc.lT.between(lo_T, hi_T)]
    print(f'Galactic N = {len(mw)}; Magellanic GP stars {len(mc)}, in the Galactic domain {len(dom)} '
          f'(LMC {(dom["sample"] == "Bowman2024_LMC").sum()}, SMC {(dom["sample"] == "Bowman2024_SMC").sum()}); '
          f'calibration stars {len(cal)}, their log nu_B20 range {np.log10(cal.nuchar_B20_tab).min():.2f}..'
          f'{np.log10(cal.nuchar_B20_tab).max():.2f}, MC log nu_GP range {np.log10(dom.nuchar).min():.2f}..{np.log10(dom.nuchar).max():.2f}')
    rows = []
    for name in conversions(cal):
        for gal in ('LMC+SMC', 'LMC', 'SMC'):
            s = dom if gal == 'LMC+SMC' else dom[dom['sample'] == f'Bowman2024_{gal}']
            B = np.c_[np.ones(len(s)), s.lL, s.lT]
            def med(calib):
                f = conversions(calib)[name]
                return np.median(f(np.log10(s.nuchar.values)) - B @ c)
            m0 = med(cal)
            boots = [med(cal.iloc[rng.integers(0, len(cal), len(cal))]) for _ in range(300)]
            # also resample the Magellanic stars
            both = []
            for _ in range(300):
                cb = cal.iloc[rng.integers(0, len(cal), len(cal))]
                sb = s.iloc[rng.integers(0, len(s), len(s))]
                Bb = np.c_[np.ones(len(sb)), sb.lL, sb.lT]
                both.append(np.median(conversions(cb)[name](np.log10(sb.nuchar.values)) - Bb @ c))
            f = conversions(cal)[name]
            rr = f(np.log10(s.nuchar.values)) - B @ c
            rows.append(dict(conversion=name, sample=gal, n=len(s), median_residual_dex=m0,
                             lo_cal=np.percentile(boots, 16), hi_cal=np.percentile(boots, 84),
                             lo_all=np.percentile(both, 2.5), hi_all=np.percentile(both, 97.5),
                             mwu_p=mannwhitneyu(rr, rmw).pvalue))
    res = pd.DataFrame(rows)
    # paper's matched-Teff medians, with the median-offset conversion
    k_mw = mw.lT.between(4.4, 4.65); k_mc = mc.lT.between(4.4, 4.65)
    off = np.median(np.log10(cal.nuchar_B20_tab) - np.log10(cal.nuchar_GP))
    print(f'matched Teff 4.4-4.65: MW median {mw.nuchar[k_mw].median():.2f} (N={k_mw.sum()}), MC GP {mc.nuchar[k_mc].median():.2f} '
          f'(N={k_mc.sum()}), MC on Lorentzian scale (median offset {off:+.3f} dex) {10 ** off * mc.nuchar[k_mc].median():.2f}')
    res.to_csv(f'{HERE}/data_obs/nuchar_Z_test.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(res.to_string(index=False, float_format=lambda v: f'{v:+.3f}'))


if __name__ == '__main__':
    main()


def model_shift():
    """Model log nu_c(LMC grid) - log nu_c(MW grid) at the Magellanic star positions (v2, omega = 0, MS + post-MS)."""
    from scipy.interpolate import LinearNDInterpolator
    D2 = f'{HERE}/data_v2'
    E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
    mc = E[E['sample'].isin(['Bowman2024_LMC', 'Bowman2024_SMC']) & np.isfinite(E.nuchar)]
    mw = E[E['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic']) & np.isfinite(E.nuchar)]
    out = {}
    for z in ('MW', 'LMC'):
        df = pd.concat([pd.read_csv(f'{D2}/mesa_ms_{z}w0.0.csv'), pd.read_csv(f'{D2}/mesa_post_{z}w0.0.csv')])
        df = df[np.isfinite(df.logTeff) & np.isfinite(df.logLspec) & np.isfinite(df.FeCZ_nuc_d) & (df.FeCZ_nuc_d > 0)]
        f = LinearNDInterpolator(np.c_[df.logTeff / 0.03, df.logLspec / 0.08], np.log10(df.FeCZ_nuc_d))
        out[z] = (f(np.c_[mc.lT / 0.03, mc.lL / 0.08]), f(np.c_[mw.lT / 0.03, mw.lL / 0.08]))
    d = out['LMC'][0] - out['MW'][0]
    ok = np.isfinite(d)
    print(f'model log nu_c(LMC) - log nu_c(MW) at the {ok.sum()} Magellanic star positions with an FeCZ in both: '
          f'median {np.median(d[ok]):+.3f} dex, 16-84% [{np.percentile(d[ok], 16):+.3f}, {np.percentile(d[ok], 84):+.3f}]')
    d2 = out['LMC'][1] - out['MW'][1]
    ok2 = np.isfinite(d2)
    print(f'  same at the {ok2.sum()} Galactic star positions: median {np.median(d2[ok2]):+.3f} dex')


if __name__ == '__main__' and os.environ.get('MODEL_SHIFT'):
    model_shift()
