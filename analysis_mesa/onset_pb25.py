#!/usr/bin/env python3
"""
Red-noise amplitude onset with the Pedersen & Bildsten 2025 stars (ingest_pb25.py).

(a) calibrated: primary sample (Bowman 2020 + Shen 2024, umag) + 49 Cygnus OB stars with alpha0 mapped to umag through the
    B20 overlap (log alpha0 = -0.296 + 0.505 log PSD; 0.32 dex scatter);
(b) internal: PB25's own homogeneous set (B20 refit with valid errors + Cygnus OB), with their native PSD alpha_0 and the
    model-independent rms, so no cross-method calibration is involved.
Hinge fit with a log Teff covariate, as in the paper (lmc_vmac_test.hinge_boot). Writes data_obs/onset_pb25.csv.
"""
import os
import numpy as np
import pandas as pd
from lmc_vmac_test import hinge_boot

HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
pb = pd.read_csv(f'{HERE}/data_obs/pb25_stars.csv')
E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
prim = E[E['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic']) & E.a0_ok & np.isfinite(E.alpha0)]
cyg = pb[(pb['sample'] == 'PB25_CygOB') & pb.pb_ok & np.isfinite(pb.lL)]
own = pb[pb.pb_ok & np.isfinite(pb.lL) & np.isfinite(pb.lT)]

sets = {
    'primary only (umag)': (prim.lL, prim.lT, np.log10(prim.alpha0)),
    'primary + CygOB calibrated (umag)': (pd.concat([prim.lL, cyg.lL]), pd.concat([prim.lT, cyg.lT]),
                                          pd.concat([np.log10(prim.alpha0), np.log10(cyg.alpha0_cal_umag)])),
    'PB25 internal: log PSD alpha_0': (own.lL, own.lT, own.log_psd0),
    'PB25 internal: log rms': (own.lL, own.lT, own.log_rms),
}
# robustness: drop the log L_spec < 2.0 star (2.5 Msun, 10.7 kK, alpha0 ~ 1.8e4 umag; not a massive-star SLF object) and,
# separately, the two B20 stars whose PB25 spectra look like coherent g modes (HD 37711, HD 27563: low nu_char, w < 1)
PULS = {'HD 37711', 'HD 27563'}
cyg2, own2 = cyg[cyg.lL >= 2.0], own[own.lL >= 2.0]
own3 = own2[~own2.star.isin(PULS)]
prim3 = prim[~prim.star.isin(PULS)]
sets.update({
    'primary + CygOB calibrated, log L >= 2': (pd.concat([prim.lL, cyg2.lL]), pd.concat([prim.lT, cyg2.lT]),
                                               pd.concat([np.log10(prim.alpha0), np.log10(cyg2.alpha0_cal_umag)])),
    'primary + CygOB calibrated, log L >= 2, no SPB-like': (pd.concat([prim3.lL, cyg2.lL]), pd.concat([prim3.lT, cyg2.lT]),
                                               pd.concat([np.log10(prim3.alpha0), np.log10(cyg2.alpha0_cal_umag)])),
    'PB25 internal rms, log L >= 2': (own2.lL, own2.lT, own2.log_rms),
    'PB25 internal rms, log L >= 2, no SPB-like': (own3.lL, own3.lT, own3.log_rms),
})
rows = []
for name, (x, T, y) in sets.items():
    x, T, y = map(np.asarray, (x, T, y))
    ok = np.isfinite(x) & np.isfinite(T) & np.isfinite(y)
    f = hinge_boot(x[ok], T[ok], y[ok])
    rows.append(dict(sample=name, **f))
res = pd.DataFrame(rows)
res.to_csv(f'{HERE}/data_obs/onset_pb25.csv', index=False, float_format='%.3f')
pd.set_option('display.width', 200)
print(res[['sample', 'n', 'n_below', 'h', 'h_lo', 'h_hi', 'dBIC', 'slope_below', 'slope_above']].to_string(
    index=False, float_format=lambda v: f'{v:.2f}'))
print('\npaper (alpha0, OB incl. LMC ppm stars, N = 222): hinge 3.16 [2.94, 3.44], dBIC +16.3, slopes -0.02 -> +1.66')
