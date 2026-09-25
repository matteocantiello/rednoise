"""Excitation spectrum of the waves launched at the top of the FeCZ.

Total wave luminosity (S3t in REPORT_transfer.md; Goldreich & Kumar 1990, sharp interface):
    L_w = M_t * F_c,max/F * L,   M_t = Mach number over the top alpha H_P (history mach_FeCZ_aver_ahp)

Distribution, parametric (the parameters are the physics still to be pinned down; see WAVE_MODEL.md):
  frequency:  dL_w/dln(omega) ∝ (omega/omega_c)^(-a)  for omega > omega_c     (a = 13/2: GK, sharp interface)
                              ∝ (omega/omega_c)^(+b)  for omega < omega_c
  degree:     at frequency omega, waves are excited by eddies of size h(omega) = alpha H_P (omega/omega_c)^(-3/2)
              (Kolmogorov, omega > omega_c; h = alpha H_P below), up to ell_max(omega) = ell_eddy (omega/omega_c)^(3/2).
              Within 1 <= ell <= ell_max:  dL_w/d ell ∝ ell^p.
              p (the low-ell slope) is the key unknown: only ell <~ 5 survive damping and disk averaging,
              so the visible fraction scales roughly as ell_max^-(p+1).
"""
import numpy as np


class Source:
    def __init__(self, omega_c, ell_eddy, L_w, a=6.5, b=0.0, p=1.0):
        self.omega_c, self.ell_eddy, self.L_w = omega_c, ell_eddy, L_w
        self.a, self.b, self.p = a, b, p

    @classmethod
    def from_history(cls, h, **kw):
        """Source parameters from a template_v2 history row (L in erg/s from log_L)."""
        L = 10 ** h['log_L'] * 3.828e33
        L_w = h['mach_FeCZ_aver_ahp'] * h['FeCZ_Fc_max'] * L
        return cls(h['FeCZ_omega_c'], h['FeCZ_ell_eddy'], L_w, **kw)

    def phi(self, omega):
        """dL_w/dln(omega) / L_w. With b = 0 the flat part below omega_c is truncated at omega_c/100."""
        x = np.asarray(omega) / self.omega_c
        if self.b > 0:
            f = np.where(x > 1, x ** -self.a, x ** self.b)
            norm = 1.0 / self.a + 1.0 / self.b
        else:
            f = np.where(x > 1, x ** -self.a, np.where(x > 1e-2, 1.0, 0.0))
            norm = 1.0 / self.a + np.log(100.0)
        return f / norm

    def ell_max(self, omega):
        x = np.maximum(np.asarray(omega) / self.omega_c, 1.0)
        return np.maximum(self.ell_eddy * x ** 1.5, 1.0)

    def ell_weights(self, omega, ells):
        """Fraction of the wave luminosity at frequency omega carried by each integer ell (shape omega x ell)."""
        omega = np.atleast_1d(omega)
        ells = np.asarray(ells, float)
        lmax = self.ell_max(omega)[:, None]
        # normalise over all 1 <= ell <= ell_max, even though only the ells passed in are returned
        norm = np.array([np.sum(np.arange(1, int(m) + 1, dtype=float) ** self.p) for m in lmax[:, 0]])[:, None]
        return np.where(ells[None, :] <= lmax, ells[None, :] ** self.p, 0.0) / norm

    def dL_dlnomega_dell(self, omega, ells):
        """Wave luminosity per ln(omega) per ell [erg/s], shape (omega, ell)."""
        return self.L_w * self.phi(omega)[:, None] * self.ell_weights(omega, ells)
