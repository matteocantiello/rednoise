"""Where do GYRE non-adiabatic modes lose energy? Shares of the total work W from radial regions.

Usage: python3 where_damped.py DETAIL.h5 [...]
W < 0 means damping. Regions in x = r/R, plus the outer non-adiabatic layer (log T < 5.0) and the FeCZ
(log T in 5.1-5.35, roughly; see the envelope figure).
"""
import sys
import h5py
import numpy as np


def cplx(x):
    return x['re'] + 1j * x['im'] if x.dtype.names else x


for name in sys.argv[1:]:
    with h5py.File(name, 'r') as f:
        x, dW, T = f['x'][()], f['dW_dx'][()], f['T'][()]
        nu = f.attrs['freq']                       # scalars are attributes in detail files
        nu = nu['re'] + 1j * nu['im'] if getattr(nu, 'dtype', None) is not None and nu.dtype.names else complex(nu)
    w = np.r_[0, 0.5 * (dW[1:] + dW[:-1]) * np.diff(x)]
    Wt = w.sum()
    lT = np.log10(T)
    print(f'{name}: nu = {nu.real:.3f}/d, Q = {nu.real / (2 * abs(nu.imag)):.3g}, W_total = {Wt:.3g}, N = {len(x)}')
    regions = [('x < 0.3 (above the core)', x < 0.3), ('0.3 <= x < 0.9', (x >= 0.3) & (x < 0.9)),
               ('0.9 <= x < 0.97', (x >= 0.9) & (x < 0.97)), ('x >= 0.97 (FeCZ and above)', x >= 0.97),
               ('  of which log T 5.1-5.35 (FeCZ)', (lT > 5.1) & (lT < 5.35)),
               ('  of which log T < 5.0 (outer layer)', lT < 5.0)]
    for lab, s in regions:
        print(f'    {lab:38s} share of W: {w[s].sum() / Wt:8.3f}')
