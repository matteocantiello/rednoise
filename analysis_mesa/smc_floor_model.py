#!/usr/bin/env python3
"""
Noise-floor model for the BLOeM SMC PSF light curves (Van Daele+2026) and injection-recovery tests.

Likelihood: Whittle (Anderson+1990) on the independent frequencies (step 1/T, 0.1-40 d-1) of the periodogram power
P = A^2, with model power
    S(nu) = [a0 / (1 + (nu/nc)^g)]^2 + [k * F_s(nu; Tmag)]^2 + Cw^2 ,
where F_s is an empirical red-noise floor: a semi-Lorentzian with sector-median shape (nu_f, g_f) and a median fractional
amplitude that scales with Tmag, measured from the stars Van Daele class as showing no significant variability. k is a
free rescaling with a lognormal prior (sigma_k = measured scatter), so the floor can be larger or smaller for each star.
Models: 'white' (Cw), 'floor' (floor + Cw), 'slf' (SLF + Cw; the naive fit, no floor), 'slf+floor' (all).
SLF detection above the floor: Delta BIC(floor -> slf+floor) > 10.

Steps:
 1. floor calibration from the non-variable stars (Whittle SLF + Cw fits of their prewhitened light curves);
 2. injection-recovery: synthetic SLF (random-phase sum of sinusoids with Rayleigh amplitudes following the semi-Lorentzian)
    injected into every non-variable light curve, over a grid of nu_in and a0_in / floor amplitude; fitted with 'slf' and
    'slf+floor'; the floor template leaves the injected star out;
 3. all SMC star-sectors fitted with every model; the Galactic SPOC light curves fitted with 'slf' (Whittle) for consistency.
Writes data_obs/floor_calib.csv, floor_injection.csv, floor_fits_smc.csv, floor_fits_mw.csv.
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault('VD_FMIN', '0.1'); os.environ.setdefault('VD_FMAX', '40')
import vandaele_fit as vf                 # noqa: E402
import smc_instrumental_check as sc       # noqa: E402
from astropy.timeseries import LombScargle  # noqa: E402

LC_SMC = '/mnt/ceph/users/mcantiello/rednoise/vandaele2026/zenodo_lightcurves/PSF'
LC_MW = '/mnt/ceph/users/mcantiello/rednoise/tess_mw/SPOC'
FMIN, FMAX = 0.1, 40.0
DBIC_DET = 10.0
NU_IN = [0.3, 0.7, 1.5, 3.0]
R_IN = [float(v) for v in os.environ.get('FLOOR_R_IN', '0.5,1,2,4,8,16').split(',')]
G_IN = 2.0


def power_spectrum(t, m):
    """Independent frequencies (step 1/T) and power in amplitude^2 units (A = sqrt(4 P_psd / N))."""
    T = t.max() - t.min()
    nu = np.arange(FMIN, FMAX, 1.0 / T)
    p = LombScargle(t, m, fit_mean=False, center_data=True, normalization='psd').power(nu, method='fast')
    return nu, 4 * p / len(t)


def lor(nu, a0, nc, g):
    return (a0 / (1 + (nu / nc) ** g)) ** 2


class Model:
    def __init__(self, nu, P, floor=None, sig_k=0.6):
        self.nu, self.P, self.floor, self.sig_k = nu, P, floor, sig_k   # floor: (Af, nuf, gf) in amplitude units

    def S(self, p, kind):
        cw2 = 10 ** (2 * p[-1])
        s = np.full_like(self.nu, cw2)
        i = 0
        if 'slf' in kind:
            s = s + lor(self.nu, 10 ** p[0], 10 ** p[1], p[2]); i = 3
        if 'floor' in kind:
            Af, nuf, gf = self.floor
            s = s + lor(self.nu, Af * 10 ** p[i], nuf, gf)
        return s

    def nll(self, p, kind):
        S = self.S(p, kind)
        v = np.sum(np.log(S) + self.P / S)
        if 'floor' in kind:
            lk = p[3] if 'slf' in kind else p[0]
            v += 0.5 * (lk / self.sig_k) ** 2
        return v

    def fit(self, kind):
        lcw0 = 0.5 * np.log10(np.median(self.P[self.nu > 0.8 * self.nu.max()]))
        la0 = 0.5 * np.log10(max(np.median(self.P[self.nu < 0.3]), 10 ** (2 * lcw0)))
        starts = []
        for lnc in (-0.6, -0.2, 0.3):
            for g in (1.5, 3.0):
                if kind == 'white':
                    starts = [[lcw0]]; break
                p = []
                if 'slf' in kind:
                    p += [la0, lnc, g]
                if 'floor' in kind:
                    p += [0.0]
                starts.append(p + [lcw0])
            if kind == 'white' or ('slf' not in kind):
                break
        if 'slf' not in kind and 'floor' in kind:
            starts = [[k0, lcw0] for k0 in (-0.5, 0.0, 0.5)]
        bounds = {'white': [(-6, 4)], 'floor': [(-3, 3), (-6, 4)], 'slf': [(-6, 5), (np.log10(0.02), np.log10(40)), (0.3, 10), (-6, 4)],
                  'slf+floor': [(-6, 5), (np.log10(0.02), np.log10(40)), (0.3, 10), (-3, 3), (-6, 4)]}[kind]
        best = None
        for p0 in starts:
            r = minimize(self.nll, p0, args=(kind,), method='L-BFGS-B', bounds=bounds)
            if best is None or r.fun < best.fun:
                best = r
        k = len(best.x)
        return best.x, best.fun, 2 * best.fun + k * np.log(len(self.nu))


def fit_all(t, m, floor=None, sig_k=0.6, kinds=('white', 'slf', 'floor', 'slf+floor')):
    nu, P = power_spectrum(t, m)
    M = Model(nu, P, floor, sig_k)
    out = {}
    for kd in kinds:
        if 'floor' in kd and floor is None:
            continue
        p, nll, bic = M.fit(kd)
        out[kd] = dict(p=p, bic=bic)
    return out


def summarise(o):
    r = {}
    for kd, v in o.items():
        r[f'bic_{kd}'] = v['bic']
        if 'slf' in kd:
            r[f'nu_{kd}'] = 10 ** v['p'][1]; r[f'a0_{kd}'] = 10 ** v['p'][0]; r[f'g_{kd}'] = v['p'][2]
        if 'floor' in kd:
            r[f'logk_{kd}'] = v['p'][3] if 'slf' in kd else v['p'][0]
        r[f'cw_{kd}'] = 10 ** v['p'][-1]
    if 'slf+floor' in o:
        r['dbic_above_floor'] = o['floor']['bic'] - o['slf+floor']['bic']
    r['dbic_slf_vs_white'] = o['white']['bic'] - o['slf']['bic']
    return r


# ---------------------------------------------------------------------------------------------------------------- 1. floor
def lc_path(gaia, sector):
    return f'{LC_SMC}/GAIA DR3 {gaia}_sector{sector}_PSFlc.txt'


def load_pw(path):
    t, m = vf.load(path)
    if t is None:
        return None, None, None
    fl = np.median(np.genfromtxt(path, invalid_raise=False)[:, 1])
    m, _ = vf.prewhiten(t, m, 'lin')
    return t, m, fl


def calib_one(args):
    gaia, sector, tmag = args
    t, m, fl = load_pw(lc_path(gaia, sector))
    if t is None:
        return None
    o = fit_all(t, m, kinds=('white', 'slf'))
    p = o['slf']['p']
    return dict(gaia=gaia, sector=sector, Tmag=tmag, flux=fl, a0=10 ** p[0], nu=10 ** p[1], g=p[2], cw=10 ** p[-1],
                dbic=o['white']['bic'] - o['slf']['bic'])


def floor_model(cal, exclude=None):
    """Per-sector shape (median nu, g) and log(a0/flux) = c0_s + c1 (Tmag - 13); returns function and sigma."""
    c = cal[cal.gaia != exclude] if exclude is not None else cal
    c = c[(c.dbic > 0) & (c.nu > 0.03)]
    y = np.log10(c.a0 / c.flux)
    X = np.c_[pd.get_dummies(c.sector).values.astype(float), c.Tmag - 13]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    sig = np.std(y - X @ beta)
    secs = list(pd.get_dummies(c.sector).columns)
    shape = c.groupby('sector').agg(nu=('nu', 'median'), g=('g', 'median'))

    def F(sector, tmag, flux):
        j = secs.index(sector)
        la = beta[j] + beta[-1] * (tmag - 13)
        return (10 ** la * flux, shape.loc[sector, 'nu'], shape.loc[sector, 'g'])
    return F, sig, beta, shape


# ---------------------------------------------------------------------------------------------------------- 2. injection
def synth(t, a0, nc, g, rng):
    T = t.max() - t.min()
    nu = np.arange(1.0 / T, FMAX, 1.0 / T)
    amp = a0 / (1 + (nu / nc) ** g)
    z = amp * (rng.normal(size=len(nu)) + 1j * rng.normal(size=len(nu))) / np.sqrt(2)
    out = np.zeros_like(t)
    for i in range(0, len(nu), 200):
        ph = 2 * np.pi * np.outer(t - t.min(), nu[i:i + 200])
        out += (np.abs(z[i:i + 200]) * np.cos(ph + np.angle(z[i:i + 200]))).sum(1)
    return out


def inject_one(args):
    gaia, sector, tmag, cal = args
    t, m, fl = load_pw(lc_path(gaia, sector))
    if t is None:
        return []
    F, sig, _, _ = floor_model(cal, exclude=gaia)
    Af, nuf, gf = F(sector, tmag, fl)
    rng = np.random.default_rng(gaia % 100000 + sector)
    rows = []
    for nu_in in NU_IN:
        for R in R_IN:
            a0_in = R * Af
            mm = m + synth(t, a0_in, nu_in, G_IN, rng)
            o = fit_all(t, mm, (Af, nuf, gf), sig)
            r = summarise(o)
            rows.append(dict(gaia=gaia, sector=sector, Tmag=tmag, nu_in=nu_in, R=R, a0_in=a0_in, Af=Af, nu_floor=nuf, **r))
    return rows


# ------------------------------------------------------------------------------------------------------------- 3. fits
def smc_one(args):
    gaia, sector, tmag, cls, cal = args
    t, m, fl = load_pw(lc_path(gaia, sector))
    if t is None:
        return None
    F, sig, _, _ = floor_model(cal, exclude=gaia if cls == 'none' else None)
    Af, nuf, gf = F(sector, tmag, fl)
    o = fit_all(t, m, (Af, nuf, gf), sig)
    return dict(gaia=gaia, sector=sector, Tmag=tmag, cls=cls, Af=Af, nu_floor=nuf, **summarise(o))


def mw_one(path):
    import re
    b = os.path.basename(path)
    g = re.search(r'(.+)_sector(\d+)_', b)
    t, m = vf.load(path)
    if t is None:
        return None
    m, _ = vf.prewhiten(t, m, 'lin')
    o = fit_all(t, m, kinds=('white', 'slf'))
    return dict(id=g.group(1), sector=int(g.group(2)), **summarise(o))


def main():
    cl = sc.classes()
    files = sorted(glob.glob(f'{LC_SMC}/*_PSFlc.txt'))
    import re
    idx = pd.DataFrame([dict(gaia=int(re.search(r'GAIA DR3 (\d+)_sector(\d+)', f).group(1)),
                             sector=int(re.search(r'GAIA DR3 (\d+)_sector(\d+)', f).group(2))) for f in files]).merge(cl, on='gaia')
    none = idx[idx.cls == 'none']
    step = os.environ.get('FLOOR_STEP', 'all')
    with Pool(int(os.environ.get('NPROC', 8))) as pool:
        if step in ('all', 'calib'):
            cal = pd.DataFrame([r for r in pool.map(calib_one, list(zip(none.gaia, none.sector, none.Tmag))) if r])
            cal.to_csv(f'{HERE}/data_obs/floor_calib.csv', index=False, float_format='%.5g')
            F, sig, beta, shape = floor_model(cal)
            print(f'floor calibration: {len(cal)} sectors; sigma(log a0/flux) = {sig:.2f}; Tmag slope {beta[-1]:+.2f}')
            print(shape.round(3).to_string())
        cal = pd.read_csv(f'{HERE}/data_obs/floor_calib.csv')
        if step in ('all', 'fits'):
            smc = pd.DataFrame([r for r in pool.map(smc_one, [(a, b, c, d, cal) for a, b, c, d in
                                                            zip(idx.gaia, idx.sector, idx.Tmag, idx.cls)]) if r])
            smc.to_csv(f'{HERE}/data_obs/floor_fits_smc.csv', index=False, float_format='%.5g')
            mw = pd.DataFrame([r for r in pool.map(mw_one, sorted(glob.glob(f'{LC_MW}/*_SPOC.txt'))) if r])
            mw.to_csv(f'{HERE}/data_obs/floor_fits_mw.csv', index=False, float_format='%.5g')
            print(smc.groupby('cls').agg(n=('gaia', 'size'), det_above_floor=('dbic_above_floor', lambda v: (v > DBIC_DET).mean()),
                                         det_naive=('dbic_slf_vs_white', lambda v: (v > DBIC_DET).mean())).round(2).to_string())
        if step in ('all', 'inject'):
            rows = pool.map(inject_one, [(a, b, c, cal) for a, b, c in zip(none.gaia, none.sector, none.Tmag)])
            inj = pd.DataFrame([r for rr in rows for r in rr])
            inj.to_csv(f'{HERE}/data_obs/floor_injection{os.environ.get("FLOOR_INJ_TAG", "")}.csv', index=False, float_format='%.5g')
            print(f'injections: {len(inj)}')


if __name__ == '__main__':
    main()
