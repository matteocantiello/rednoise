#!/usr/bin/env python3
"""
Joint test of complete FeCZ -> surface scenarios: one physical chain must predict v_macro,
alpha0 and nu_char at once, each with its physical exponent (beta = 1): only a
normalisation (and a floor for v_macro, alpha0) is fitted per observable.

Also pushes each scenario's fitted predictions, evaluated at the observed stars, through the
same regressions used for the data (L_spec+Teff plane; tau+log M) to get the implied
sensitivity ratios d log alpha0 / d log v_macro, to compare with 5.3 / 6.4 / 2.2.

Writes data/scenarios_<Z>_w<w>.csv.
"""
import os, sys
import numpy as np
import pandas as pd
from test_transfer import OBSSETS, predict, V, A, N, fit_floor, fit_lin, bic, D, Z, W
import transfer_models as tm

SCEN = {
    'S0 MLT at the FeCZ':               ('v_c,max', 'F_c/F', '1/(2pi t_c)'),
    'S1 incoherent cells':              ('v_c,max', 'F_c/F x cells', '1/(2pi t_c)'),
    'S2 Cantiello+09 microturbulence':  ('C09 microturb', 'M_s^2 [C09 microturb]', 'v_s/(2pi H_s) [C09 microturb]'),
    'S3 M_x F_c waves, v saturated':    ('sat[wave M_x F_c, v_g=v]', 'M_x F_c/F', '1/(2pi t_c)'),
    'S3t M_t F_c waves, v saturated':   ('sat[wave M_t F_c, v_g=v]', 'M_t F_c/F', '1/(2pi t_c)'),
    'S3c as S3, cells-diluted flux':    ('sat[wave M_x F_c, v_g=v]', 'M_x F_c/F x cells', '1/(2pi t_c)'),
    'S4 core IGW (undamped)':           ('core IGW, v_g=c_s', 'core IGW flux/F*', 'core turnover'),
    # post hoc (2026-09-25): the best single predictor per observable on the full v2 MW grid; physically, the
    # 3D picture (Schultz+2023): photospheric v ~ FeCZ v, brightness ~ Mach-scaled convective flux
    'S5 v_c + M_t F_c/F + turnover':    ('v_c,max', 'M_t F_c/F', '1/(2pi t_c)'),
}
OBS_RATIO = {'logL': 5.31, 'tau': 6.37, 'logM': 2.19}
NOFECZ = ('core', 'dynamical', 'nu_max')


def evaluate(obs, cands, key, quad=False):
    o = OBSSETS[obs]
    Xs = predict(cands[key][0], not key.startswith(NOFECZ), o['x'], o['y'])
    sel = np.isfinite(Xs) & np.isfinite(o['v'])
    if obs == 'nuchar':
        sel &= Xs > 0
    y, X = o['v'][sel], Xs[sel]
    Bm = np.c_[np.ones(sel.sum()), o['y'][sel], o['x'][sel]]
    bres = y - Bm @ np.linalg.lstsq(Bm, y, rcond=None)[0]
    if obs == 'nuchar':
        p, res = fit_lin(X, y, 1.0); k = 1
        yhat = y - res
        par = dict(a=p[0])
    else:
        p, res = fit_floor(X, y, 1.0, quad=quad); k = 2
        yhat = y - res
        par = dict(floor=10**p[0], a=p[1])
    out = dict(n=int(sel.sum()), rms=np.std(res), dBIC=bic(res, k) - bic(bres, 3), **par)
    return out, sel, yhat


def slopes(o, sel, yhat):
    """Regress the scenario prediction as the data were: plane and tau+logM."""
    x, yL, tau, lM = o['x'][sel], o['y'][sel], o['tau'][sel], o['logM'][sel]
    Bm = np.c_[np.ones(len(yhat)), yL, x]
    cL = np.linalg.lstsq(Bm, yhat, rcond=None)[0]
    Em = np.c_[np.ones(len(yhat)), tau, lM]
    ok = np.isfinite(tau) & np.isfinite(lM)
    cE = np.linalg.lstsq(Em[ok], yhat[ok], rcond=None)[0]
    return dict(logL=cL[1], logTeff=cL[2], tau=cE[1], logM=cE[2])



def main():
    rows = []
    for name, (kv, ka, kn) in SCEN.items():
        r = dict(scenario=name, v_pred=kv, a_pred=ka, nu_pred=kn)
        tot = 0
        for obs, cands, key in (('vmac', V, kv), ('vmac_hot', V, kv), ('alpha0', A, ka), ('nuchar', N, kn)):
            try:
                ev, sel, yhat = evaluate(obs, cands, key)
            except Exception as e:
                print(name, obs, 'failed', e); continue
            for k2, v2 in ev.items():
                r[f'{obs}_{k2}'] = v2
            if obs in ('vmac', 'alpha0', 'nuchar'):
                tot += ev['dBIC']
            sl = slopes(OBSSETS[obs], sel, yhat)
            for k2, v2 in sl.items():
                r[f'{obs}_slope_{k2}'] = v2
        r['joint_dBIC'] = tot
        # bootstrap over stars (each observable resampled independently) for the implied ratios
        nb = int(os.environ.get('NBOOT_SCEN', '0'))
        if nb:
            rng = np.random.default_rng(3)
            br = {ax: [] for ax in ('logL', 'tau', 'logM')}
            base = {}
            for obs, cands, key in (('vmac', V, kv), ('alpha0', A, ka)):
                o = OBSSETS[obs]
                Xs = predict(cands[key][0], not key.startswith(NOFECZ), o['x'], o['y'])
                base[obs] = (o, np.where(np.isfinite(Xs) & np.isfinite(o['v']))[0], Xs)
            for _ in range(nb):
                sl = {}
                for obs, (o, idx, Xs) in base.items():
                    ii = rng.choice(idx, len(idx))
                    oo = {k: (v[ii] if isinstance(v, np.ndarray) else v) for k, v in o.items()}
                    p_, res_ = fit_floor(Xs[ii], oo['v'], 1.0)
                    sl[obs] = slopes(oo, np.ones(len(ii), bool), oo['v'] - res_)
                for ax in br:
                    br[ax].append(sl['alpha0'][ax] / sl['vmac'][ax])
            for ax in br:
                r[f'ratio_{ax}_lo'], r[f'ratio_{ax}_hi'] = np.percentile(br[ax], [16, 84])
        for ax in ('logL', 'tau', 'logM'):
            r[f'ratio_{ax}'] = r.get(f'alpha0_slope_{ax}', np.nan) / r.get(f'vmac_slope_{ax}', np.nan)
        rows.append(r)

    res = pd.DataFrame(rows)
    res.to_csv(f'{D}/scenarios_{Z}_w{W}.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    print('dBIC < 0: better than the empirical L_spec+Teff plane (3 params) with 1-2 params; beta = 1 everywhere\n')
    print(res[['scenario', 'vmac_dBIC', 'vmac_a', 'vmac_floor', 'vmac_hot_dBIC', 'alpha0_dBIC', 'alpha0_a', 'alpha0_floor',
               'nuchar_dBIC', 'nuchar_a', 'joint_dBIC']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    if 'ratio_logL_lo' in res:
        print('\nratio intervals (16-84%, bootstrap over stars):')
        print(res[['scenario'] + [f'ratio_{a}{e}' for a in ('logL', 'tau', 'logM') for e in ('', '_lo', '_hi')]].to_string(index=False, float_format=lambda v: f'{v:.2f}'))
    print('\nImplied slopes (predictions at the observed stars, regressed as the data) and ratios')
    print(f"observed: vmac dlogL +0.29 dlogT +0.10 dtau +0.27 dlogM +0.60 | alpha0 dlogL +1.55 dlogT -3.13 dtau +1.73 dlogM +1.32 "
          f"| nu dlogT +1.53 dlogL -0.23 dtau -0.33 | ratios {OBS_RATIO}")
    print(res[['scenario', 'vmac_slope_logL', 'vmac_slope_logTeff', 'vmac_slope_tau', 'vmac_slope_logM',
               'alpha0_slope_logL', 'alpha0_slope_logTeff', 'alpha0_slope_tau', 'alpha0_slope_logM',
               'nuchar_slope_logTeff', 'nuchar_slope_logL', 'nuchar_slope_tau',
               'ratio_logL', 'ratio_tau', 'ratio_logM']].to_string(index=False, float_format=lambda v: f'{v:+.2f}'))


if __name__ == '__main__':
    main()
