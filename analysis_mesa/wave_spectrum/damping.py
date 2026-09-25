"""Radiative damping of internal gravity waves between the FeCZ top and the photosphere.

Same discretisation as igw_tau in template_v2/src/run_star_extras.f90 (quasi-adiabatic WKB,
Zahn, Talon & Matias 1997, N_T = N), vectorised over (omega, ell):

  tau_rad = [l(l+1)]^(3/2) sum K N^3/omega^4 * min(sqrt(N^2/(N^2-omega^2)), 10) dr / r^3   (N > omega)
  tau_ev  = 2 sum k_h sqrt(1 - N^2/omega^2) dr,  k_h = sqrt(l(l+1))/r                     (N < omega)
  K = 16 sigma T^3 / (3 kappa rho^2 c_p)

Cells run from the photosphere (tau >= 2/3) down to the cell above the FeCZ top. The transmitted
wave energy flux is ~ exp(-tau_rad - tau_ev).

Caveat: quasi-adiabatic. It is least reliable just below the photosphere, where the thermal time
approaches the wave period.
"""
import numpy as np

from .mesa_io import SIGMA_SB

SECDAY = 86400.0


def path_arrays(prof, k_top):
    """Per-cell arrays along the path from the photosphere to the cell above cell k_top."""
    k = np.arange(max(prof.k_phot, 0), min(k_top, prof.nz - 1))
    dr = prof.r[k] - prof.r[k + 1]
    rmid = 0.5 * (prof.r[k] + prof.r[k + 1])
    N2 = 0.5 * (prof.N2[k] + prof.N2[k + 1])
    K = 16 * SIGMA_SB * prof.T[k] ** 3 / (3 * prof.kappa[k] * prof.rho[k] ** 2 * prof.cp[k])
    return dict(k=k, dr=dr, rmid=rmid, N2=N2, K=K)


def q_nonad(path, omega, ell):
    """Local non-adiabaticity q = K k_r^2 / omega (thermal diffusion rate on a wavelength / wave frequency).
    The quasi-adiabatic tau is valid only where q << 1."""
    omega, ell = np.broadcast_arrays(np.asarray(omega, float), np.asarray(ell, float))
    w2 = (omega ** 2)[..., None]
    kh2 = (ell * (ell + 1))[..., None] / path['rmid'] ** 2
    kr2 = kh2 * np.clip(path['N2'] / w2 - 1, 0, None)
    return path['K'] * kr2 / omega[..., None]


def tau(path, omega, ell, stop_at_q=None):
    """tau_rad, tau_ev for angular frequencies omega [rad/s] and degrees ell (broadcast together).

    stop_at_q: if set (e.g. 1.0), integrate only from the FeCZ top up to the first cell (going outward)
    where q_nonad >= stop_at_q: the "transition layer" above which the wave's temperature fluctuations
    are radiated within a wave period. There, wave energy becomes emergent-flux fluctuations (the signal)
    instead of being dissipated. Default None = the run_star_extras integral (whole path to tau = 2/3).
    """
    omega, ell = np.broadcast_arrays(np.asarray(omega, float), np.asarray(ell, float))
    w2 = (omega ** 2)[..., None]
    L2 = (ell * (ell + 1))[..., None]
    N2 = path['N2']
    prop = N2 > w2
    if stop_at_q is not None:
        # cells are ordered surface -> FeCZ top: keep only the cells below the deepest one with q >= stop_at_q
        hi = q_nonad(path, omega, ell) >= stop_at_q
        idx = np.arange(N2.size)
        deepest = np.where(hi, idx, -1).max(-1, keepdims=True)
        keep = idx > deepest
    else:
        keep = True
    prop = prop & keep
    fac = np.minimum(np.sqrt(np.where(prop, N2 / np.where(prop, N2 - w2, 1.0), 1.0)), 10.0)
    t_rad = np.where(prop, L2 ** 1.5 * path['K'] * np.clip(N2, 0, None) ** 1.5 / (w2 * w2) * fac
                     * path['dr'] / path['rmid'] ** 3, 0.0).sum(-1)
    t_ev = np.where(prop | ~np.broadcast_to(keep, prop.shape), 0.0, 2 * np.sqrt(L2) / path['rmid']
                    * np.sqrt(np.clip(1 - np.clip(N2, 0, None) / w2, 0, None)) * path['dr']).sum(-1)
    return t_rad, t_ev


def nu_damp(path, ell, lo=1e-10, hi=1e-2):
    """Cyclic frequency [1/d] at which tau_rad = 1 (bisection in log omega, as in run_star_extras)."""
    a, b = np.log(lo), np.log(hi)
    if tau(path, np.exp(b), ell)[0] > 1:
        return np.exp(b) / (2 * np.pi) * SECDAY
    if tau(path, np.exp(a), ell)[0] < 1:
        return np.exp(a) / (2 * np.pi) * SECDAY
    for _ in range(40):
        m = 0.5 * (a + b)
        if tau(path, np.exp(m), ell)[0] > 1:
            a = m
        else:
            b = m
    return np.exp(0.5 * (a + b)) / (2 * np.pi) * SECDAY


def tau_profile(path, omega, ell):
    """Cumulative tau_rad from the FeCZ top outwards (for locating where the damping accrues)."""
    w2, L2 = omega ** 2, ell * (ell + 1)
    N2 = path['N2']
    prop = N2 > w2
    fac = np.minimum(np.sqrt(np.where(prop, N2 / np.where(prop, N2 - w2, 1.0), 1.0)), 10.0)
    d = np.where(prop, L2 ** 1.5 * path['K'] * np.clip(N2, 0, None) ** 1.5 / w2 ** 2 * fac
                 * path['dr'] / path['rmid'] ** 3, 0.0)
    return np.cumsum(d[::-1])[::-1]   # tau accumulated between each cell and the FeCZ top


def tau_cumulative(path, omega, ell):
    """Quasi-adiabatic tau_rad accumulated from the FeCZ top out to each path cell, for arrays of
    omega [rad/s] and ell (broadcast). Returns shape (..., ncell); cell order as in path (surface first)."""
    omega, ell = np.broadcast_arrays(np.asarray(omega, float), np.asarray(ell, float))
    w2 = (omega ** 2)[..., None]
    L2 = (ell * (ell + 1))[..., None]
    N2 = path['N2']
    prop = N2 > w2
    fac = np.minimum(np.sqrt(np.where(prop, N2 / np.where(prop, N2 - w2, 1.0), 1.0)), 10.0)
    d = np.where(prop, L2 ** 1.5 * path['K'] * np.clip(N2, 0, None) ** 1.5 / (w2 * w2) * fac
                 * path['dr'] / path['rmid'] ** 3, 0.0)
    return np.flip(np.cumsum(np.flip(d, -1), -1), -1)


def nu_nonad(path, ell):
    """Cyclic frequency [1/d] below which waves of degree ell are non-adiabatic (q = K k_r^2/omega > 1)
    at each path cell. For N >> omega, q = 1 gives omega^3 = K k_h^2 N^2."""
    kh2 = ell * (ell + 1) / path['rmid'] ** 2
    w = (path['K'] * kh2 * np.clip(path['N2'], 0, None)) ** (1 / 3)
    return np.where(path['N2'] > 0, w / (2 * np.pi) * SECDAY, np.nan)
