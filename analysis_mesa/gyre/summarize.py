"""Summarize a non-adiabatic GYRE run: linewidths vs mode spacing, photometric response.

Usage: python3 summarize.py OUTDIR
Linewidth (FWHM of the power Lorentzian, cyclic) = 2|Im(sigma)|/(2 pi) = |Im(freq)|/pi ... with freq in d^-1:
Gamma_nu = 2 |Im(nu)| since sigma = 2 pi nu.  Overlap ratio = Gamma_nu / (spacing to the neighbouring mode).
"""
import sys
import h5py
import numpy as np

d = sys.argv[1]
with h5py.File(f'{d}/summary_nad.h5', 'r') as f:
    S = {k: f[k][()] for k in f.keys()}
with h5py.File(f'{d}/summary_ad.h5', 'r') as f:
    A = {k: f[k][()] for k in f.keys()}


def cplx(x):
    return x['re'] + 1j * x['im'] if x.dtype.names else x


nu = cplx(S['freq'])
nu_ad = cplx(A['freq']).real
l = S['l']
print(f'{d}: {len(nu)} non-adiabatic modes ({len(nu_ad)} adiabatic)')
print(' l   n_pg   nu_r[1/d]   Gamma[1/d]  spacing  Gamma/spacing  Q=nu/Gamma    f_T      |lagL/L|/|xi_r|   eta')
for ll in (1, 2, 3):
    s = np.nonzero(l == ll)[0]
    s = s[np.argsort(nu[s].real)]
    nr = nu[s].real
    gam = 2 * np.abs(nu[s].imag)
    sp = np.gradient(nr) if len(nr) > 2 else np.full_like(nr, np.nan)
    na = np.sum(A['l'] == ll)
    print(f'-- l={ll}: {len(s)} nad / {na} ad modes; median Gamma/spacing {np.median(gam / sp):.3g}')
    lagL = np.abs(cplx(S['lag_L_ref'][s]))
    xir = np.abs(cplx(S['xi_r_ref'][s]))
    for j, i in enumerate(s):
        if j % max(1, len(s) // 12) == 0 or j == len(s) - 1:
            print(f"{ll:2d} {S['n_pg'][i]:6d} {nr[j]:10.4f} {gam[j]:11.3e} {sp[j]:8.4f} {gam[j] / sp[j]:12.3e} "
                  f"{nr[j] / gam[j]:10.3g} {S['f_T'][i]:9.3g} {lagL[j] / max(xir[j], 1e-300):12.3g} {S['eta'][i]:9.3g}")
