"""Read MESA profiles/history and locate the FeCZ, mirroring template_v2/src/run_star_extras.f90.

Profiles supply the propagation region (T, rho, kappa, c_p, N^2, r, tau). The wave source (omega_c,
Mach number, convective flux) comes from the matching history row: template_v2 computes it from the
MLT velocity (mlt_vc), whereas the profile conv_vel includes rotational mixing.
"""
import os
import numpy as np

RSUN = 6.957e10          # cm (MESA r26 const_def)
LSUN = 3.828e33          # erg/s
SIGMA_SB = 5.670374419e-5

T_FECZ = (1e5, 5e5)      # run_star_extras: zone type set by T of its bottom cell
T_LIMIT = 1e6            # only regions above T = 1e6 K are examined
CONVECTIVE_MIXING = 1    # const_def: convective_mixing


def read_mesa_table(path):
    """Return (header dict, data dict of numpy arrays) for a MESA history/profile file."""
    with open(path) as f:
        lines = f.readlines()
    hnames, hvals = lines[1].split(), lines[2].split()
    header = {}
    for n, v in zip(hnames, hvals):
        try:
            header[n] = float(v)
        except ValueError:
            header[n] = v.strip('"')
    cols = lines[5].split()
    rows = [l.split() for l in lines[6:]]
    rows = [r for r in rows if len(r) == len(cols)]   # drop partially written lines
    arr = np.array(rows, dtype=float)
    return header, {c: arr[:, i] for i, c in enumerate(cols)}


class Profile:
    """One MESA profile in cgs units, surface first (cell 0 = surface), as in MESA."""

    def __init__(self, path):
        self.path = path
        self.header, d = read_mesa_table(path)
        self.model_number = int(self.header['model_number'])
        self.r = d['radius'] * RSUN
        self.T = d['temperature']
        self.rho = 10 ** d['logRho']
        self.P = d['pressure']
        self.kappa = d['opacity']
        self.cp = d['cp']
        self.N2 = d['brunt_N2']
        self.cs = d['csound']
        self.Hp = d['pressure_scale_height'] * RSUN
        self.tau = d['tau']
        self.L = d['luminosity'] * LSUN
        self.mlt_type = d['mlt_mixing_type'].astype(int)
        self.Lconv_frac = d['lum_conv_div_L']
        self.nz = len(self.r)
        self.R = self.r[0]
        self.k_phot = int(np.argmax(self.tau >= 2.0 / 3.0))   # first cell with tau >= 2/3
        self.Teff = (self.L[0] / (4 * np.pi * self.R**2 * SIGMA_SB)) ** 0.25

    def conv_regions(self):
        """Contiguous runs of MLT-convective cells above T_LIMIT: list of (top, bottom) cell indices."""
        above = np.nonzero(self.T < T_LIMIT)[0]
        n_limit = above.max() if len(above) else 0
        regions, top = [], None
        for k in range(n_limit + 1):
            if self.mlt_type[k] == CONVECTIVE_MIXING:
                if top is None:
                    top = k
                bottom = k
            elif top is not None:
                regions.append((top, bottom))
                top = None
        if top is not None:
            regions.append((top, bottom))
        return regions

    def fecz(self):
        """(top, bottom) cell indices of the FeCZ, or None. Several fragments: keep the max L_conv/L."""
        cands = [(t, b) for t, b in self.conv_regions() if T_FECZ[0] < self.T[b] < T_FECZ[1]]
        if not cands:
            return None
        return max(cands, key=lambda tb: self.Lconv_frac[tb[0]:tb[1] + 1].max())


class Track:
    """A MESA run directory: history plus its profiles, paired by model number."""

    def __init__(self, run_dir):
        self.run_dir = run_dir
        logs = os.path.join(run_dir, 'LOGS')
        _, self.hist = read_mesa_table(os.path.join(logs, 'history.data'))
        idx = np.loadtxt(os.path.join(logs, 'profiles.index'), skiprows=1, ndmin=2).astype(int)
        self.profile_paths = {m: os.path.join(logs, f'profile{p}.data') for m, _, p in idx}

    def history_row(self, model_number):
        """History values at a model number (last occurrence, i.e. after any retries/restarts)."""
        i = np.nonzero(self.hist['model_number'] == model_number)[0]
        if len(i) == 0:
            raise KeyError(f'model {model_number} not in history (buffered output not flushed yet?)')
        return {c: v[i[-1]] for c, v in self.hist.items()}

    def pairs(self, xc_max=0.7, xc_min=1e-3):
        """(Profile, history row) for every profile on the main sequence."""
        out = []
        for m, path in sorted(self.profile_paths.items()):
            if not os.path.exists(path):
                continue
            try:
                h = self.history_row(m)
            except KeyError:
                continue
            if xc_min < h['center_h1'] < xc_max:
                out.append((Profile(path), h))
        return out
