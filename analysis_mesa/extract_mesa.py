#!/usr/bin/env python3
"""
Extract per-timestep FeCZ (and HeII/HeI/HI) properties from the rednoise MESA grid into
analysis tables (handoff §8.2).

Outputs (analysis_mesa/data/):
  mesa_ms.csv    main sequence of every track with a ZAMS: rows picked at 201 uniform
                 steps in fractional MS age tau (nearest model in age, no interpolation, so
                 zones switching on/off are not blended). tau is age-based and only defined
                 for tracks that reached the TAMS (Xc < 1e-3); for tracks still on the MS the
                 rows are picked uniformly in s = 1 - Xc/Xc,ZAMS and tau = NaN.
                 Column w_time = MS time represented by the row [yr].
  mesa_post.csv  post-MS part, <= 400 rows per track uniform in arc length on the sHRD
                 (for maps), with w_time = time represented by each row.

Conventions: log L_spec as in the paper (ell_sun = 5777^4/(274*100));
nu_c = 1/(2 pi t_c) [1/d] with t_c = alpha_MLT H_P,aver / v_aver (Cantiello+2021; this
is what run_star_extras writes as turnover_<zone>); alpha_MLT = 1.6 in the grid.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/models/grids')
import grid_io as g

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
os.makedirs(OUT, exist_ok=True)
XC_TAMS = 1e-3
N_MS = 201
N_POST = 400
RSUN, MSUN, DAY = 6.957e10, 1.989e33, 86400.0

ZONES = ['FeCZ', 'HeII', 'HeI', 'HI']


def zone_cols(d, z):
    """Per-zone quantities; NaN where the zone is absent (MESA writes 0 or -1e99)."""
    def pos(c):
        v = np.asarray(d[c], float) if c in d else np.full(len(d['model_number']), np.nan)
        return np.where(np.isfinite(v) & (v > 0), v, np.nan)
    out = {
        f'{z}_vmax_kms': pos(f'v_{z}_max') / 1e5,
        f'{z}_vaver_kms': pos(f'v_{z}_aver') / 1e5,
        f'{z}_FcF': pos(f'{z}_Fc_max'),
        f'{z}_tc_s': pos(f'turnover_{z}'),
        f'{z}_mach_max': pos(f'mach_{z}_max'),
        f'{z}_rho': pos(f'rho_{z}_aver'),
        f'{z}_hp_Rsun': pos(f'{z}_hp_aver') / RSUN,
        f'{z}_mass_Msun': pos(f'{z}_mass') / MSUN,
    }
    out[f'{z}_nuc_d'] = DAY / (2 * np.pi * out[f'{z}_tc_s'])
    R = 10**np.asarray(d['log_R'], float) * RSUN
    rt, rb = pos(f'{z}_r_top'), pos(f'{z}_r_bottom')
    out[f'{z}_depth_top_R'] = 1 - rt / R        # depth of the zone top below the surface, in R
    out[f'{z}_thick_R'] = (rt - rb) / R
    out[f'{z}_present'] = np.isfinite(out[f'{z}_vmax_kms']) & np.isfinite(out[f'{z}_FcF'])
    return out


def base_cols(d):
    lt, lg = np.asarray(d['log_Teff'], float), np.asarray(d['log_g'], float)
    out = dict(model_number=d['model_number'], age=d['star_age'], M=d['star_mass'],
               logTeff=lt, logg=lg, logL=d['log_L'], logR=d['log_R'],
               logLspec=g.spec_ell(lt, lg), Xc=d['center_h1'], Yc=d['center_he4'],
               vsurf_kms=d.get('surf_avg_v_rot', np.full(len(lt), np.nan)),
               Gamma_surf=d.get('surf_avg_Lrad_div_Ledd', np.full(len(lt), np.nan)),
               n_subsurf=d.get('subsurface_convective_regions', np.full(len(lt), np.nan)))
    for z in ZONES:
        out.update(zone_cols(d, z))
    n = len(lt)
    get = lambda c: np.asarray(d[c], float) if c in d else np.full(n, np.nan)
    pos = lambda c: np.where(get(c) > 0, get(c), np.nan)
    # extra inputs for the transfer models (transfer_models.py)
    out.update(rho_surf=pos('rho_surf'), cs_surf_kms=pos('photosphere_csound') / 1e5,
               FeCZ_cs_kms=pos('cs_FeCZ_aver') / 1e5, FeCZ_mach_aver=pos('mach_FeCZ_aver'),
               FeCZ_mach_ahp=pos('mach_FeCZ_aver_ahp'), FeCZ_rho_ahp=pos('rho_FeCZ_aver_ahp'),
               FeCZ_vsurf_C09_kms=pos('v_FeCZ_surf') / 1e5, FeCZ_beq_G=pos('b_FeCZ_aver'),
               core_vmax_kms=pos('v_max_core') / 1e5, core_vaver_kms=pos('v_aver_core') / 1e5,
               core_rho=pos('rho_aver_core'), core_hp_Rsun=pos('hp_aver_core') / RSUN,
               core_r_Rsun=pos('r_core') / RSUN, core_mach_max=pos('mach_max_core'),
               core_tc_s=pos('turnover_core'))
    # v averaged over the top alpha*H_P of the FeCZ, recovered from the Cantiello+09 surface
    # velocity v_surf = v_ahp * sqrt(M_ahp rho_ahp / rho_surf)
    out['FeCZ_vahp_kms'] = out['FeCZ_vsurf_C09_kms'] / np.sqrt(out['FeCZ_mach_ahp'] * out['FeCZ_rho_ahp'] / out['rho_surf'])
    # depth of the FeCZ top in pressure scale heights below the surface, from the radii at
    # which ln P - ln P_surf = N (N = 1..100), interpolated at r_top
    NS = np.array([1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 50, 100])
    rhp = np.vstack([get(f'r_hp_{k}') for k in NS]).T          # decreasing with N
    rtop = pos('FeCZ_r_top')
    nhp = np.full(n, np.nan)
    for i in np.where(np.isfinite(rtop) & np.all(np.isfinite(rhp), axis=1))[0]:
        r = rhp[i]
        if rtop[i] >= r[0]:
            nhp[i] = NS[0] * (10**out['logR'][i] * RSUN - rtop[i]) / max(10**out['logR'][i] * RSUN - r[0], 1.0)
        elif rtop[i] > r[-1]:
            nhp[i] = np.interp(-rtop[i], -r, NS)
    out['FeCZ_nHP_top'] = nhp
    return out


def rows(d, idx):
    return base_cols({k: np.asarray(v)[idx] for k, v in d.items()})


# optional: python3 extract_mesa.py MW/w0.0 [LMC/w0.2 ...] -> mesa_ms_<tag>.csv (subset only)
SEL = [tuple(a.split('/')) for a in sys.argv[1:]] or g.SUBGRIDS
TAG = '' if len(sys.argv) == 1 else '_' + '_'.join(a.replace('/', '') for a in sys.argv[1:])

ms_tabs, post_tabs = [], []
for z, w in SEL:
    tracks = g.load_subgrid(z, w, verbose=False)
    for m, d in sorted(tracks.items()):
        h, age = d['center_h1'], d['star_age']
        i0 = g.find_zams(h)
        if i0 == 0 or h[i0:].min() > h[i0] - 0.01:
            continue  # never reached the ZAMS
        it = g.find_tams(h, XC_TAMS)
        status = g.model_status(z, w, m)['status']
        if it is not None:
            t0, t1 = age[i0], age[it]
            tgrid = np.linspace(t0, t1, N_MS)
            idx = np.clip(np.searchsorted(age, tgrid), i0, it)
            tau = (age[idx] - t0) / (t1 - t0)
            wt = np.full(N_MS, (t1 - t0) / (N_MS - 1)); wt[[0, -1]] /= 2
        else:
            s = 1 - h[i0:] / h[i0]
            sgrid = np.linspace(0, s.max(), N_MS)
            idx = i0 + np.clip(np.searchsorted(s, sgrid), 0, len(s) - 1)
            tau = np.full(N_MS, np.nan)
            wt = np.gradient(age[idx])
        r = rows(d, idx)
        r.update(Z=z, w=float(w[1:]), Minit=m, status=status, tau=tau,
                 s_xc=1 - h[idx] / h[i0], w_time=wt, ms_complete=it is not None)
        ms_tabs.append(pd.DataFrame(r))

        if it is not None and it < len(h) - 5:
            lt = d['log_Teff'][it:]
            ls = g.spec_ell(lt, d['log_g'][it:])
            step = np.hypot(np.diff(lt) / 0.03, np.diff(ls) / 0.08)
            cum = np.r_[0, np.cumsum(np.nan_to_num(step))]
            sel = np.unique(np.minimum(np.searchsorted(cum, np.linspace(0, cum[-1], N_POST)), len(lt) - 1))
            idx = it + sel
            edges = np.r_[age[idx[0]], 0.5 * (age[idx[1:]] + age[idx[:-1]]), age[-1]]
            r = rows(d, idx)
            r.update(Z=z, w=float(w[1:]), Minit=m, status=status, t_since_tams=age[idx] - age[it],
                     w_time=np.diff(edges))
            post_tabs.append(pd.DataFrame(r))
    print(f'{z}/{w}: {len(tracks)} tracks', flush=True)

ms = pd.concat(ms_tabs, ignore_index=True)
post = pd.concat(post_tabs, ignore_index=True)
ms.to_csv(f'{OUT}/mesa_ms{TAG}.csv', index=False, float_format='%.6g')
post.to_csv(f'{OUT}/mesa_post{TAG}.csv', index=False, float_format='%.6g')
print(f'mesa_ms.csv: {len(ms)} rows, {ms.groupby(["Z", "w", "Minit"]).ngroups} tracks '
      f'({ms.groupby(["Z", "w", "Minit"]).ms_complete.first().sum()} with complete MS)')
print(f'mesa_post.csv: {len(post)} rows')
