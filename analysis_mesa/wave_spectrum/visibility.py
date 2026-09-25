"""Disk-integration weights for photometry.

Write the emergent-flux fluctuation over the surface as dF/F = sum_lm a_lm Y_lm, with every m of a degree
equally excited, and let sigma_l^2 be the surface variance carried by degree l. By the addition theorem the
variance of the disk-integrated dL/L is then independent of inclination:

    var(dL/L) = sum_l W_l sigma_l^2,   W_l = [ int_0^1 P_l(mu) h(mu) mu dmu / int_0^1 h(mu) mu dmu ]^2

h(mu) is the limb-darkening law; the default is Eddington, h = 1 + 3/2 mu. W_0 = 1; W_l falls steeply with l,
and odd l > 1 are strongly suppressed.
"""
import numpy as np
from numpy.polynomial import legendre


def weights(ells, limb=lambda mu: 1 + 1.5 * mu, n=4000):
    mu = (np.arange(n) + 0.5) / n                   # midpoint rule on [0, 1]
    hmu = limb(mu) * mu
    denom = hmu.sum()
    out = []
    for l in np.atleast_1d(ells).astype(int):
        c = np.zeros(l + 1)
        c[l] = 1
        out.append((legendre.legval(mu, c) * hmu).sum() / denom)
    return np.array(out) ** 2
