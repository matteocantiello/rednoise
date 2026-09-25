#!/usr/bin/env python3
"""
Refit of the SLF variability of the BLOeM SMC stars (Van Daele+2026, arXiv 2605.15757) from their PSF TESS light curves
(Zenodo 20540863, on Ceph at /mnt/ceph/users/mcantiello/rednoise/vandaele2026).

Method, following their Sect. 6 as closely as documented:
- one sector at a time; flux -> mmag, 5-sigma clip, mean removed; no detrending;
- amplitude spectrum A(nu) from 1/(2T) to the Nyquist frequency, oversampling 10;
- least-squares fit of alpha(nu) = alpha0 / (1 + (nu/nu_char)^gamma) + C_w (their Eq. 1);
- prewhitening (their Period04 step, done by hand and without a stated stopping rule): here automatic, iteratively removing
  the highest peak while it exceeds SNR_PW times the current SLF model, up to 40 frequencies (sinusoid fit to the light curve);
- cuts: nu_char >= 0.15 d-1; BIC of the SLF model better than a flat white-noise line.
Two fit metrics are run: residuals in log A ('log') and in A ('lin'), since theirs is not stated.
Writes data_obs/vandaele_fits.csv (one row per star, sector, metric, prewhitening on/off).
"""
import os
import glob
import re
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, curve_fit
from astropy.timeseries import LombScargle
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
LC = os.environ.get('VD_LC', '/mnt/ceph/users/mcantiello/rednoise/vandaele2026/zenodo_lightcurves/PSF')
GLOB = os.environ.get('VD_GLOB', '*_PSFlc.txt')
COL = int(os.environ.get('VD_COL', 1))                # flux column (Galactic SPOC files: 1 = PDCSAP, 2 = SAP)
# optional degrade test: add white noise until alpha0/C_w matches 10**VD_MATCH_A0CW (e.g. the SMC median, 1.34)
MATCH = float(os.environ['VD_MATCH_A0CW']) if os.environ.get('VD_MATCH_A0CW') else None
SEED = 11
SNR_PW = 5.0
NPW_MAX = 40


def load(path):
    a = np.genfromtxt(path, invalid_raise=False)
    if a.ndim != 2 or a.shape[0] < 200:
        return None, None
    t, f = a[:, 0], a[:, COL]
    ok = np.isfinite(t) & np.isfinite(f) & (f > 0)
    t, f = t[ok], f[ok]
    m = f - np.median(f) if os.environ.get('VD_FLUX', '1') == '1' else -2.5 * np.log10(f / np.median(f)) * 1e3   # e-/s (as theirs) or mmag
    for _ in range(3):
        s = np.std(m - np.mean(m)); k = abs(m - np.mean(m)) < 5 * s
        t, m = t[k], m[k]
    return t, m - np.mean(m)


FMIN, FMAX = float(os.environ.get('VD_FMIN', 0.1)), float(os.environ.get('VD_FMAX', 40.0))


def amp_spectrum(t, m):
    # frequency range as in their periodogram figures (0.1-40 d-1); VD_FMIN=0 / VD_FMAX=0 give 1/(2T) and Nyquist
    T = t.max() - t.min()
    dt = np.median(np.diff(t))
    nu = np.arange(FMIN if FMIN > 0 else 0.5 / T, FMAX if FMAX > 0 else 0.5 / dt, 0.1 / T)
    p = LombScargle(t, m, fit_mean=False, center_data=True, normalization='psd').power(nu, method='fast')
    return nu, np.sqrt(4 * p / len(t))


def slf(nu, a0, nc, g, cw):
    return a0 / (1 + (nu / nc) ** g) + cw


def fit_slf(nu, A, metric):
    # log-spaced binning weight is not used: every periodogram point counts once, as in a plain least-squares fit
    p0 = [np.median(A[nu < 0.3]) if (nu < 0.3).any() else A[0], 1.0, 2.0, np.median(A[nu > 0.8 * nu.max()])]
    lo, hi = [1e-6, 0.02, 0.3, 1e-7], [1e5, 50, 30, 1e4]
    p0 = np.clip(p0, np.array(lo) * 1.01, np.array(hi) * 0.99)
    if metric == 'log':
        fun = lambda p: np.log10(slf(nu, *p)) - np.log10(A)
    else:
        fun = lambda p: slf(nu, *p) - A
    r = least_squares(fun, p0, bounds=(lo, hi), x_scale='jac', max_nfev=4000)
    res = r.fun
    n = len(nu)
    rss = np.sum(res ** 2)
    try:
        J = r.jac; cov = np.linalg.inv(J.T @ J) * rss / (n - 4); err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(4, np.nan)
    # white-noise-only comparison model
    if metric == 'log':
        c = np.mean(np.log10(A)); rss0 = np.sum((np.log10(A) - c) ** 2)
    else:
        rss0 = np.sum((A - np.mean(A)) ** 2)
    bic1 = n * np.log(rss / n) + 4 * np.log(n)
    bic0 = n * np.log(rss0 / n) + 1 * np.log(n)
    return r.x, err, bic0 - bic1


def prewhiten(t, m, metric):
    m = m.copy(); freqs = []
    for _ in range(NPW_MAX):
        nu, A = amp_spectrum(t, m)
        p, _, _ = fit_slf(nu, A, metric)
        snr = A / slf(nu, *p)
        i = np.argmax(snr)
        if snr[i] < SNR_PW:
            break
        f0 = nu[i]
        model = lambda tt, f, a, ph: a * np.sin(2 * np.pi * (f * tt + ph))
        try:
            q, _ = curve_fit(model, t, m, p0=[f0, A[i], 0.0])
            q2, _ = curve_fit(model, t, m, p0=[f0, A[i], 0.25])
            q = q if np.sum((m - model(t, *q)) ** 2) < np.sum((m - model(t, *q2)) ** 2) else q2
        except RuntimeError:
            break
        m = m - model(t, *q); freqs.append((q[0], abs(q[1])))
    return m, freqs


def one(path):
    base = os.path.basename(path)
    gid = re.search(r'(.+)_sector(\d+)_', base)
    sid, sector = gid.group(1), int(gid.group(2))
    g = re.match(r'GAIA DR3 (\d+)$', sid)
    gaia = int(g.group(1)) if g else -1
    t, m = load(path)
    if t is None:
        print('skipped (no usable data):', os.path.basename(path)); return []
    added = 0.0
    if MATCH is not None:
        nu, A = amp_spectrum(t, m)
        p, _, _ = fit_slf(nu, A, 'lin')
        cw_t = p[0] / 10 ** MATCH
        if cw_t > p[3]:
            added = np.sqrt(len(t) / np.pi) * np.sqrt(cw_t ** 2 - p[3] ** 2)
            m = m + np.random.default_rng(SEED + sector).normal(0, added, len(t))
    rows = []
    for metric in ('log', 'lin'):
        for pw in (False, True):
            mm, freqs = prewhiten(t, m, metric) if pw else (m, [])
            nu, A = amp_spectrum(t, mm)
            p, e, dbic = fit_slf(nu, A, metric)
            rows.append(dict(id=sid, gaia=gaia, sector=sector, noise_added=added, metric=metric, prewhitened=pw, n_pw=len(freqs),
                             alpha0=p[0], nuchar=p[1], gamma=p[2], Cw=p[3], e_alpha0=e[0], e_nuchar=e[1], e_gamma=e[2],
                             e_Cw=e[3], dBIC_vs_white=dbic, T_span=t.max() - t.min(), n_pts=len(t),
                             dt_min=np.median(np.diff(t)) * 1440))
    return rows


def main():
    files = sorted(glob.glob(f'{LC}/{GLOB}'))
    with Pool(8) as pool:
        out = pool.map(one, files)
    d = pd.DataFrame([r for rr in out for r in rr])
    d['pass_nu'] = d.nuchar >= 0.15
    d['pass_bic'] = d.dBIC_vs_white > 0
    d.to_csv(f'{HERE}/data_obs/vandaele_fits{os.environ.get("VD_TAG", "")}.csv', index=False, float_format='%.5g')
    g = d.groupby(['metric', 'prewhitened'])
    print(g.apply(lambda x: pd.Series(dict(n=len(x), pass_nu=(x.pass_nu).sum(), pass_both=(x.pass_nu & x.pass_bic).sum(),
                                           stars_both=x[x.pass_nu & x.pass_bic].id.nunique(),
                                           med_lognu=np.log10(x.nuchar[x.pass_nu]).median()))).to_string())


if __name__ == '__main__':
    main()
