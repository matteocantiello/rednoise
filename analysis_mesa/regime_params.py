"""Convective-regime parameters of the FeCZ from the 3D envelope literature, computed on a MESA grid.

For each profile (by default only the triggered ones, priority 10 in profiles.index) this finds the Fe opacity
peak (max kappa with 1e5 < T < 5e5 K) and evaluates:

  tau_Fe          optical depth at the Fe peak
  tau_crit        c P_rad / ((P_rad + P_gas) v_c) = c (1 - beta) / v_c   (Jiang+2015; Schultz+2022, 2023)
  tau_ratio       tau_Fe / tau_crit  (< 1: plumes radiatively lossy, convection inefficient)
  Gamma_Fe        kappa_Fe L(r) / (4 pi G M c)  (local Eddington factor at the Fe peak)
  Y_Fe            pseudo-Mach number L / (4 pi r^2 a T^4) / c_iso, c_iso = sqrt(P_gas/rho)  (Schultz+2020)
  mach_max        max v_c / c_s in the FeCZ (MLT)
  v_c             max MLT convective velocity in the FeCZ [km/s]   (profile conv_vel: use NON-ROTATING grids only)
  t_th_zone       thermal time of the FeCZ, int_zone c_p T dm / L  [s]  (c_p T = P delta/(rho grad_ad))
  t_th_above      thermal time of everything above the FeCZ bottom [s]
  t_th_peak       thermal time of the material above the Fe peak [s]
  t_th_top        thermal time of the material above the FeCZ top [s]
  nu_turn         v_c / (2 pi alpha H_P) at the Fe peak [1/d]
  nu_th_zone, nu_th_above, nu_th_peak, nu_th_top   1/(2 pi t_th) [1/d]

Usage (from analysis_mesa/):
    RN_GRID=../models/grids_v2 python3 regime_params.py MW/w0.0 [--all-profiles] [--threads 16]
Writes $RN_DATA (default data/)/regime_<Z><w>.csv
"""
import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models', 'grids'))
import grid_io as g                                                     # noqa: E402

C, G, A_RAD, MSUN, RSUN, LSUN = 2.99792458e10, 6.6743e-8, 7.5657e-15, 1.98841e33, 6.957e10, 3.828e33
ALPHA = 1.6
ELL_SUN = 5777.0 ** 4 / 27400.0
COLS = ['mass', 'temperature', 'logRho', 'pressure', 'opacity', 'radius', 'luminosity', 'mlt_mixing_type',
        'csound', 'tau', 'pgas_div_ptotal', 'log_conv_vel', 'grada', 'pressure_scale_height']


def read_profile(path):
    with open(path) as f:
        f.readline()
        hn = f.readline().split()
        hv = f.readline().split()
        f.readline(); f.readline()
        cols = f.readline().split()
    head = dict(zip(hn, hv))
    idx = [cols.index(c) for c in COLS]
    d = np.loadtxt(path, skiprows=6, usecols=idx)
    return head, {c: d[:, i] for i, c in enumerate(COLS)}


def regime(path):
    head, d = read_profile(path)
    T, kap, tau = d['temperature'], d['opacity'], d['tau']
    r = d['radius'] * RSUN
    Lr = d['luminosity'] * LSUN
    M = float(head['star_mass']) * MSUN
    Teff = float(head['Teff'])
    R = float(head['photosphere_r']) * RSUN
    out = dict(model_number=int(head['model_number']), Minit=float(head['initial_mass']),
               M=float(head['star_mass']), Xc=float(head['center_h1']), logTeff=np.log10(Teff),
               logL=np.log10(float(head['photosphere_L'])))
    out['logLspec'] = np.log10(Teff ** 4 / (G * M / R ** 2) / ELL_SUN)
    band = np.nonzero((T > 1e5) & (T < 5e5))[0]
    if len(band) == 0:
        return out
    k = band[np.argmax(kap[band])]
    beta = d['pgas_div_ptotal'][k]
    rho = 10 ** d['logRho'][k]
    Pg = beta * d['pressure'][k]
    out.update(tau_Fe=tau[k], beta_Fe=beta, kap_Fe=kap[k], T_Fe=T[k], r_Fe_R=r[k] / R,
               Gamma_Fe=kap[k] * Lr[k] / (4 * np.pi * G * M * C),
               Y_Fe=Lr[k] / (4 * np.pi * r[k] ** 2 * A_RAD * T[k] ** 4) / np.sqrt(Pg / rho))
    # FeCZ: convective cells in the Fe band, contiguous around the zone with the largest v_c
    conv = (d['mlt_mixing_type'] == 1) & (T > 1e5) & (T < 5e5)
    if conv.any():
        vc = np.where(conv, 10 ** d['log_conv_vel'], 0.0)
        j = int(np.argmax(vc))
        top = j
        while top > 0 and conv[top - 1]:
            top -= 1
        bot = j
        while bot < len(T) - 1 and conv[bot + 1]:
            bot += 1
        vmax = vc[j]
        dm = np.abs(np.diff(d['mass'] * MSUN, append=d['mass'][-1] * MSUN))   # cell masses (surface first)
        # c_p T = P delta / (rho grad_ad), delta = (4 - 3 beta)/beta for ideal gas + radiation (grid profiles
        # have no cp column; ionization terms in delta are neglected -- fine for an order-unity timescale)
        b_all = d['pgas_div_ptotal']
        cpT = d['pressure'] * (4 - 3 * b_all) / b_all / (10 ** d['logRho'] * d['grada'])
        out.update(v_c=vmax / 1e5, mach_max=np.max((vc / d['csound'])[top:bot + 1]),
                   tau_crit=C * (1 - beta) / vmax,
                   t_th_zone=np.sum(cpT[top:bot + 1] * dm[top:bot + 1]) / Lr[k],
                   t_th_above=np.sum(cpT[:bot + 1] * dm[:bot + 1]) / Lr[k],
                   t_th_peak=np.sum(cpT[:k] * dm[:k]) / Lr[k],       # material above the Fe peak
                   t_th_top=np.sum(cpT[:top] * dm[:top]) / Lr[k],    # material above the FeCZ top
                   nu_turn=vmax / (2 * np.pi * ALPHA * d['pressure_scale_height'][k] * RSUN) * 86400,
                   FeCZ_ncells=bot - top + 1, tau_top=tau[top], tau_bot=tau[bot])
        out['tau_ratio'] = out['tau_Fe'] / out['tau_crit']
        out['nu_th_zone'] = 86400 / (2 * np.pi * out['t_th_zone'])
        out['nu_th_above'] = 86400 / (2 * np.pi * out['t_th_above'])
        out['nu_th_peak'] = 86400 / (2 * np.pi * out['t_th_peak'])
        out['nu_th_top'] = 86400 / (2 * np.pi * max(out['t_th_top'], 1.0))
    else:
        out.update(tau_ratio=0.0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('subgrid')
    ap.add_argument('--all-profiles', action='store_true')
    ap.add_argument('--threads', type=int, default=16)
    a = ap.parse_args()
    z, w = a.subgrid.split('/')
    jobs = []
    for m in g.MASSES:
        logs = f'{g.model_dir(z, w, m)}/LOGS'
        try:
            idx = np.loadtxt(f'{logs}/profiles.index', skiprows=1, ndmin=2).astype(int)
        except (OSError, ValueError):
            continue
        for mn, prio, pn in idx:
            if a.all_profiles or prio >= 10:
                jobs.append(f'{logs}/profile{pn}.data')
    print(f'{a.subgrid}: {len(jobs)} profiles')
    with ThreadPoolExecutor(a.threads) as ex:
        rows = [r for r in ex.map(lambda p: _safe(p), jobs) if r]
    df = pd.DataFrame(rows).sort_values(['Minit', 'model_number'])
    out = os.environ.get('RN_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'))
    fn = f'{out}/regime_{z}{w}.csv'
    df.to_csv(fn, index=False, float_format='%.6g')
    print(f'wrote {fn}: {len(df)} rows, {df.Minit.nunique()} masses')


def _safe(p):
    try:
        return regime(p)
    except Exception as e:                                   # partially written profile, etc.
        print('skip', p, e)
        return None


if __name__ == '__main__':
    main()
