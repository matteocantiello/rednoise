#!/usr/bin/env python3
"""
Ingest Pedersen & Bildsten 2025 (MNRAS 539, 2742; Zenodo 15261328): SLF fits for 49 Cygnus OB stars and a refit, with the
same method, of 53 stars of the Bowman+2020 sample (B20).

Their alpha_0 is a power density at zero frequency (log ppm^2/uHz) and nu_char is in uHz; ours (Bowman 2020, Shen 2024) is
an amplitude-spectrum alpha_0 in umag and nu_char in 1/d. The B20 overlap gives the empirical mapping.

Writes data_obs/pb25_calibration.csv (overlap comparison and fitted relations) and data_obs/pb25_stars.csv (all 102 stars:
native PB25 values, the model-independent rms and nu_50, calibrated alpha0 in umag and nu_char in 1/d, log L_spec, log Teff).
"""
import os
import sys
import warnings
import numpy as np
import pandas as pd
from astropy.io import ascii
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/pedersen_bildsten2025/tables'
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
OUT = f'{HERE}/data_obs'
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/skill')
from kernel import star_key          # noqa: E402

UHZ_TO_D = 86400e-6                  # 1 uHz = 0.0864 1/d
ELL_SUN = 5777.0 ** 4 / 27400.0


def rd(name):
    """MRT table; the PB25 files declare numbers as A4/A5 strings, so coerce everything numeric-looking."""
    t = ascii.read(f'{RAW}/{name}', format='mrt').to_pandas()
    for c in t.columns:
        if c in ('Gaia', 'TIC', 'star', 'Group', 'FFI sectors', '2-min sectors', 'FFI', '2-min', 'star name', 'Variability'):
            continue
        v = pd.to_numeric(t[c], errors='coerce')
        if v.notna().sum() >= 0.5 * t[c].notna().sum():
            t[c] = v
    return t


def bisector(x, y):
    """OLS-bisector slope and intercept (Isobe et al. 1990): symmetric in x and y."""
    b1 = np.polyfit(x, y, 1)[0]
    b2 = 1.0 / np.polyfit(y, x, 1)[0]
    b = (b1 * b2 - 1 + np.sqrt((1 + b1 ** 2) * (1 + b2 ** 2))) / (b1 + b2)
    return b, np.mean(y) - b * np.mean(x)


def main():
    b1, b2 = rd('table_B1.mrt'), rd('table_B2.mrt')
    b1['TIC'] = b1.TIC.astype(str).str.replace('TIC ', '').str.strip()
    b2['TIC'] = b2.TIC.astype(str).str.strip()
    b20 = b1.merge(b2, on='TIC', how='inner')
    b20['key'] = b20.star.map(star_key)
    E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
    ours = E[E['sample'] == 'Bowman2020_Galactic'].copy()
    ours['key'] = ours.star.map(star_key)
    m = b20.merge(ours[['key', 'star', 'lT', 'lL', 'alpha0', 'nuchar', 'gamma']], on='key', how='inner')
    # failed / fixed fits have zero quoted uncertainty on alpha_0
    m['pb_ok'] = (m.e_log_alpha_0 > 0) | (m.E_log_alpha_0 > 0)
    print(f'B20 refit: {len(b20)} stars; matched to our Bowman 2020 entries: {len(m)} ({m.pb_ok.sum()} with nonzero PB25 errors)')

    g = m[m.pb_ok & np.isfinite(m.alpha0) & np.isfinite(m.nuchar)]
    x, y = g.log_alpha_0.values, np.log10(g.alpha0.values)
    s_amp, i_amp = bisector(x, y)
    res_amp = y - (i_amp + s_amp * x)
    xr = g.log_RMS.values
    s_rms, i_rms = bisector(xr, y)
    res_rms = y - (i_rms + s_rms * xr)
    lnu_pb = np.log10(g.nu_char.values * UHZ_TO_D)
    dnu = np.log10(g.nuchar.values) - lnu_pb
    lnu50 = np.log10(g.nu_50.values * UHZ_TO_D)
    s_n50, i_n50 = bisector(lnu50, np.log10(g.nuchar.values))
    cal = [
        dict(quantity='log alpha0[umag] vs log PSD alpha_0 [ppm^2/uHz]', n=len(g), slope=s_amp, intercept=i_amp,
             scatter_dex=np.std(res_amp), spearman=spearmanr(x, y).correlation),
        dict(quantity='log alpha0[umag] vs log RMS [ppm]', n=len(g), slope=s_rms, intercept=i_rms,
             scatter_dex=np.std(res_rms), spearman=spearmanr(xr, y).correlation),
        dict(quantity='log nu_char (Bowman, 1/d) - log nu_char (PB25, converted)', n=len(g), slope=np.nan,
             intercept=np.median(dnu), scatter_dex=1.4826 * np.median(np.abs(dnu - np.median(dnu))),
             spearman=spearmanr(lnu_pb, np.log10(g.nuchar)).correlation),
        dict(quantity='log nu_char (Bowman, 1/d) vs log nu_50 (PB25, converted)', n=len(g), slope=s_n50, intercept=i_n50,
             scatter_dex=np.std(np.log10(g.nuchar.values) - (i_n50 + s_n50 * lnu50)),
             spearman=spearmanr(lnu50, np.log10(g.nuchar)).correlation),
        dict(quantity='gamma (Bowman) - gamma (PB25)', n=len(g), slope=np.nan, intercept=np.median(g.gamma_y - g.gamma_x)
             if 'gamma_y' in g else np.median(g['gamma'] - g['gamma']), scatter_dex=np.nan, spearman=np.nan),
    ]
    cal = pd.DataFrame(cal)
    cal.to_csv(f'{OUT}/pb25_calibration.csv', index=False, float_format='%.4f')
    pd.set_option('display.width', 200)
    print(cal.to_string(index=False, float_format=lambda v: f'{v:.3f}'))

    # ---- all 102 PB25 stars, native + calibrated
    a1, a3 = rd('table_A1.mrt'), rd('table_A3.mrt')
    a1['Gaia'] = a1.Gaia.astype(str); a3['Gaia'] = a3.Gaia.astype(str)
    cyg = a1.merge(a3, on='Gaia', how='inner')
    for c in ('log_L_spec', 'log_alpha_0', 'e_log_alpha_0', 'E_log_alpha_0', 'nu_char', 'gamma', 'log_RMS', 'nu_50', 'w'):
        cyg[c] = pd.to_numeric(cyg[c], errors='coerce')
    cyg = pd.DataFrame(dict(star='Gaia EDR3 ' + cyg.Gaia, sample='PB25_CygOB', lT=cyg.log_Teff, lL=cyg.log_L_spec,
                            logL=cyg.log_L, M=cyg.M, log_psd0=cyg.log_alpha_0, e_psd=cyg.e_log_alpha_0,
                            E_psd=cyg.E_log_alpha_0, nu_char_uHz=cyg.nu_char, gamma=cyg.gamma, log_rms=cyg.log_RMS,
                            nu50_uHz=cyg.nu_50, w=cyg.w))
    # B20 refit: L_spec from our catalogue (their table gives classical L only)
    b20m = b20.merge(ours[['key', 'star', 'lT', 'lL']], on='key', how='left', suffixes=('', '_ours'))
    b20s = pd.DataFrame(dict(star=b20m.star, sample='PB25_B20refit', lT=b20m.lT.fillna(b20m.log_Teff), lL=b20m.lL,
                             logL=b20m.log_L, M=np.nan, log_psd0=b20m.log_alpha_0, e_psd=b20m.e_log_alpha_0,
                             E_psd=b20m.E_log_alpha_0, nu_char_uHz=b20m.nu_char, gamma=b20m.gamma, log_rms=b20m.log_RMS,
                             nu50_uHz=b20m.nu_50, w=b20m.w))
    allpb = pd.concat([cyg, b20s], ignore_index=True)
    allpb['pb_ok'] = (allpb.e_psd > 0) | (allpb.E_psd > 0)
    allpb['alpha0_cal_umag'] = 10 ** (i_amp + s_amp * allpb.log_psd0)
    allpb['nuchar_cal_d'] = allpb.nu_char_uHz * UHZ_TO_D * 10 ** cal.intercept[2]
    allpb.to_csv(f'{OUT}/pb25_stars.csv', index=False, float_format='%.5g')
    c = allpb[allpb['sample'] == 'PB25_CygOB']
    print(f"\nCygOB: {len(c)} stars, log L_spec {c.lL.min():.2f}-{c.lL.max():.2f} (median {c.lL.median():.2f}); "
          f"below log L_spec 3.16: {(c.lL < 3.16).sum()}; log Teff {c.lT.min():.2f}-{c.lT.max():.2f}; ok fits {c.pb_ok.sum()}")
    print(f"our primary sample below 3.16: {(E[E['sample'].isin(['Bowman2020_Galactic','Shen2024_Galactic'])].lL < 3.16).sum()}")


if __name__ == '__main__':
    main()
