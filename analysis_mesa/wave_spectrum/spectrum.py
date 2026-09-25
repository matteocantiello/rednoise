"""Surface photometric spectrum of FeCZ-excited waves, and the red-noise fit used on the observations.

    PSD(nu) ∝ sum_l W_l * R(omega) * [dL_w,l/dnu * exp(-tau_rad - tau_ev)] / L

- dL_w,l/dnu: source.Source (L_w = M_t F_c/F L, split in omega and ell)
- tau: damping.tau (FeCZ top -> photosphere; validated against run_star_extras)
- W_l: visibility.weights (disk integration)
- R(omega) = (omega/omega_c)^r: brightness response per unit wave flux. It is a placeholder for the
  non-adiabatic surface response (plan step 3); r = 0 means "emergent-flux fluctuations ∝ wave flux".
  The overall normalisation is not predicted yet, so only the SHAPE (nu_char, gamma) and the relative
  visible fraction are meaningful.

The observed amplitude spectrum is fitted as alpha(nu) = alpha0 / (1 + (nu/nu_char)^gamma) + C_w
(Bowman et al. 2019, 2020). The model has no white noise, so C_w = 0 and alpha(nu) ∝ sqrt(PSD).
"""
import numpy as np
from scipy.optimize import curve_fit

from . import damping, visibility
from .source import Source

SECDAY = 86400.0


def model_spectrum(prof, hist, nu=None, ells=np.arange(1, 41), r=0.0, stop_at_q=None, **src_kw):
    """Model PSD shape for one profile + matching history row.

    nu: cyclic frequencies [1/d] (default: TESS-like 0.01..25 d^-1, linear steps of 0.01)
    stop_at_q: None = damp all the way to the photosphere (quasi-adiabatic, run_star_extras); a number
               (e.g. 1) = damp only below the non-adiabatic transition layer (see damping.tau)
    Returns a dict with nu, psd (relative), per-ell psd, the visible fraction and diagnostics.
    """
    if nu is None:
        nu = np.arange(0.01, 25.0 + 1e-9, 0.01)
    omega = 2 * np.pi * nu / SECDAY
    k_top, k_bot = prof.fecz()
    path = damping.path_arrays(prof, k_top)
    src = Source.from_history(hist, **src_kw)

    dL = src.dL_dlnomega_dell(omega, ells)                      # (nu, ell) erg/s per ln omega
    t_rad, t_ev = damping.tau(path, omega[:, None], ells[None, :].astype(float), stop_at_q=stop_at_q)
    trans = np.exp(-(t_rad + t_ev))
    W = visibility.weights(ells)
    R = (omega / src.omega_c) ** r
    psd_l = W[None, :] * R[:, None] * dL * trans / nu[:, None] / (10 ** hist['log_L'] * 3.828e33)
    psd = psd_l.sum(1)

    dlnnu = np.gradient(np.log(nu))
    L_surf_vis = np.sum((W[None, :] * dL * trans).sum(1) * dlnnu)      # visible transmitted luminosity
    L_surf = np.sum((dL * trans).sum(1) * dlnnu)                        # transmitted luminosity, all ell
    return dict(nu=nu, psd=psd, psd_l=psd_l, ells=ells, W=W,
                L_w=src.L_w, f_trans=L_surf / src.L_w, f_vis=L_surf_vis / src.L_w,
                nu_c=src.omega_c / (2 * np.pi) * SECDAY, nu_peak=nu[np.argmax(psd)],
                nu_damp_l1=damping.nu_damp(path, 1.0), k_top=k_top, k_phot=prof.k_phot)


def red_noise(nu, alpha0, nu_char, gamma):
    return alpha0 / (1 + (nu / nu_char) ** gamma)


def fit_red_noise(nu, amp, nu_min=0.02, nu_max=25.0):
    """Fit alpha0, nu_char, gamma to an amplitude spectrum, in log amplitude (C_w = 0)."""
    s = (nu >= nu_min) & (nu <= nu_max) & (amp > 0)
    x, y = nu[s], np.log(amp[s])
    p0 = (amp[s].max(), nu[s][np.argmax(amp[s])] if np.argmax(amp[s]) > 0 else 1.0, 3.0)
    f = lambda x, a0, nc, g: np.log(red_noise(x, np.exp(a0), np.exp(nc), g))
    popt, pcov = curve_fit(f, x, y, p0=(np.log(p0[0]), np.log(max(p0[1], 0.05)), p0[2]), maxfev=20000)
    resid = y - f(x, *popt)
    return dict(alpha0=np.exp(popt[0]), nu_char=np.exp(popt[1]), gamma=popt[2], rms_dex=np.std(resid) / np.log(10))
