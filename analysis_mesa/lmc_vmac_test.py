#!/usr/bin/env python3
"""
Metallicity test with macroturbulence: does the v_macro onset hinge move from the Galaxy to the LMC as the v2 grids
predict?

Data
  MW : IACOB compilation (handoff macroturbulence_evol.csv, 832 stars) + the Holgado+2022 O stars not in it
       (J/A+A/665/A150, same IACOB-BROAD method).
  LMC: Serebriakova+2024 (J/A+A/692/A245; LMC values from Serebriakova+2023 and Gebruers+2022), galaxy assigned
       from coordinates. There are no stars in common with IACOB, so the v_macro scales are not tied; the hinge
       LOCATION is insensitive to a multiplicative offset, the absolute level is not.
Method (as in the paper, section 3.4): log v = c0 + c1 (x - h) + c2 max(0, x - h) + c3 log Teff, x = log L_spec,
  hinge h scanned on a 0.02-dex grid with >= 10 stars on each side; dBIC = BIC_single - BIC_broken; 16-84% from
  300 bootstrap resamples of the stars.
Model prediction: the same hinge fit applied to the model v_c,max (with the fitted floor, beta = 1) at each star's
  own position in its own galaxy's v2 grid (MS + post-MS models, non-rotating).

Writes $RN_DATA/lmc_vmac_test.csv and figures_v2/fig_lmc_vmac.{png,pdf}.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D2 = os.environ.get('RN_DATA', f'{HERE}/data_v2')
FIG = f'{HERE}/figures_v2'
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
RAW = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vizier'
sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/skill')
from kernel import vz_read_tsv, star_key          # noqa: E402  (skill helpers; see handoff Appendix A)
rng = np.random.default_rng(7)
NBOOT = 300


# ------------------------------------------------------------------ data
def load_mw():
    V = pd.read_csv(f'{OBS}/macroturbulence_evol.csv')
    V = V[np.isfinite(V.vmac) & (V.vmac > 0) & np.isfinite(V.logL_sp)]
    mw = pd.DataFrame(dict(name=V.key, lT=V.logTeff_sp, lL=V.logL_sp, vmac=V.vmac, vsini=V.vsini, src='IACOB'))
    h = vz_read_tsv(f'{RAW}/holgado22.tsv')
    h['k'] = h.Target.map(star_key)
    have = set(pd.read_csv(f'{OBS}/macroturbulence_catalog.csv').key.map(star_key))
    h = h[~h.k.isin(have)].copy()
    vm = pd.to_numeric(h.vmacro, errors='coerce')
    T = pd.to_numeric(h.Teff, errors='coerce') * 1e3                     # kK (ReadMe)
    L = pd.to_numeric(h.logLsp, errors='coerce')                         # log L_spec / L_sun (ReadMe G1)
    ok = np.isfinite(vm) & (vm > 0) & np.isfinite(T) & np.isfinite(L) & (h.l_logLsp.fillna('') == '') & \
         (h.l_Teff.fillna('') == '')
    add = pd.DataFrame(dict(name=h.Target[ok], lT=np.log10(T[ok]), lL=L[ok], vmac=vm[ok],
                            vsini=pd.to_numeric(h.vsini, errors='coerce')[ok], src='Holgado22'))
    return pd.concat([mw, add], ignore_index=True)


def load_lmc():
    s = vz_read_tsv(f'{RAW}/serebriakova24.tsv')
    ra, de = pd.to_numeric(s._RAJ2000, errors='coerce'), pd.to_numeric(s._DEJ2000, errors='coerce')
    lmc = (ra > 60) & (ra < 100) & (de < -63) & (de > -75)
    vm = pd.to_numeric(s.Vmacro, errors='coerce')
    ok = lmc & np.isfinite(vm) & (vm > 0)
    return pd.DataFrame(dict(name=s.Name[ok], lT=np.log10(s.Teff[ok].astype(float)), lL=s.Ls[ok].astype(float),
                             vmac=vm[ok], vsini=pd.to_numeric(s['Vsini-fourier'], errors='coerce')[ok],
                             src='Serebriakova24/' + s.Samp[ok]))


# ------------------------------------------------------------------ hinge fit
def bic(res, k):
    n = len(res)
    return n * np.log(np.sum(res ** 2) / n) + k * np.log(n)


def hinge_fit(x, T, y, grid=np.arange(2.2, 4.31, 0.02), nmin=10):
    A1 = np.c_[np.ones(len(y)), x, T]
    r1 = y - A1 @ np.linalg.lstsq(A1, y, rcond=None)[0]
    b1 = bic(r1, 3)
    best = (np.inf, np.nan, None)
    for h in grid:
        if (x < h).sum() < nmin or (x >= h).sum() < nmin:
            continue
        A = np.c_[np.ones(len(y)), x - h, np.maximum(0, x - h), T]
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        b = bic(y - A @ c, 5)
        if b < best[0]:
            best = (b, h, c)
    return dict(h=best[1], dBIC=b1 - best[0], slope_below=best[2][1] if best[2] is not None else np.nan,
                slope_above=(best[2][1] + best[2][2]) if best[2] is not None else np.nan)


def hinge_boot(x, T, y):
    f = hinge_fit(x, T, y)
    hs = []
    for _ in range(NBOOT):
        i = rng.integers(0, len(y), len(y))
        hs.append(hinge_fit(x[i], T[i], y[i])['h'])
    f['h_lo'], f['h_hi'] = np.nanpercentile(hs, [16, 84])
    f['n'] = len(y)
    f['n_below'] = int((x < f['h']).sum())
    return f


# ------------------------------------------------------------------ model prediction
class GridMap:
    """Model v_c,max [km/s] on (log Teff, log L_spec) from a v2 sub-grid extract (MS + post-MS)."""
    SX, SY = 0.03, 0.08

    def __init__(self, z):
        df = pd.concat([pd.read_csv(f'{D2}/mesa_ms_{z}w0.0.csv'), pd.read_csv(f'{D2}/mesa_post_{z}w0.0.csv')])
        df = df[np.isfinite(df.logTeff) & np.isfinite(df.logLspec)]
        self.pts = np.c_[df.logTeff / self.SX, df.logLspec / self.SY]
        v = df.FeCZ_vmax_kms.values
        ok = np.isfinite(v) & (v > 0)
        self.v = LinearNDInterpolator(self.pts[ok], np.log10(v[ok]))
        self.pres = LinearNDInterpolator(self.pts, df.FeCZ_present.values.astype(float))
        self.tree = cKDTree(self.pts)

    def __call__(self, lT, lL):
        P = np.c_[lT / self.SX, lL / self.SY]
        out = 10 ** self.v(P)
        pf = self.pres(P)
        out = np.where(np.isfinite(pf) & (pf < 0.5), 0.0, out)
        out[self.tree.query(P)[0] > 1.5] = np.nan
        return out


def floor_fit(X, y):
    """log y = log10(10^f + 10^a X) (beta = 1); returns the prediction."""
    from scipy.optimize import least_squares
    lx = np.log10(np.where(X > 0, X, np.nan))
    pos = X > 0
    m = lambda p: np.log10(10 ** p[0] + np.where(pos, 10 ** (p[1] + np.nan_to_num(lx)), 0.0))
    r = least_squares(lambda p: y - m(p), [np.percentile(y, 10), 0.0])
    return m(r.x), r.x


def main():
    mw, lmc = load_mw(), load_lmc()
    print(f'MW: {len(mw)} stars ({(mw.src == "Holgado22").sum()} new from Holgado+22); LMC: {len(lmc)} stars')
    Tlo, Thi = lmc.lT.min(), lmc.lT.max()
    rows = []
    samples = {'MW all': mw, f'MW, log Teff {Tlo:.2f}-{Thi:.2f} (LMC range)': mw[mw.lT.between(Tlo, Thi)],
               'LMC all': lmc}
    maps = {'MW': GridMap('MW'), 'LMC': GridMap('LMC')}
    for name, d in samples.items():
        x, T, y = d.lL.values, d.lT.values, np.log10(d.vmac.values)
        f = hinge_boot(x, T, y)
        rows.append(dict(sample=name, kind='observed', **f))
        z = 'LMC' if name.startswith('LMC') else 'MW'
        X = maps[z](T, x)
        sel = np.isfinite(X)
        yhat, p = floor_fit(X[sel], y[sel])
        fm = hinge_boot(x[sel], T[sel], yhat)
        rows.append(dict(sample=name, kind=f'model ({z} v2 grid, v_c,max + floor)', floor_kms=10 ** p[0], norm=10 ** p[1], **fm))
    res = pd.DataFrame(rows)
    res.to_csv(f'{D2}/lmc_vmac_test.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(res[['sample', 'kind', 'n', 'n_below', 'h', 'h_lo', 'h_hi', 'dBIC', 'slope_below', 'slope_above']].to_string(
        index=False, float_format=lambda v: f'{v:.2f}'))

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.2), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.88, wspace=0.08)
    edges = np.arange(1.0, 4.8, 0.2)
    for ax, (lab, d, col) in zip(axes, (('MW (IACOB + Holgado 2022)', mw, '#2a78d6'),
                                         ('LMC (Serebriakova 2023/2024)', lmc, '#eb6834'))):
        ax.scatter(d.lL, d.vmac, s=5, color=col, alpha=0.35, lw=0)
        b = []
        for a, c in zip(edges[:-1], edges[1:]):
            k = d.lL.between(a, c)
            if k.sum() >= 4:
                b.append((0.5 * (a + c), *np.percentile(d.vmac[k], [50, 16, 84])))
        b = np.array(b)
        ax.errorbar(b[:, 0], b[:, 1], yerr=[b[:, 1] - b[:, 2], b[:, 3] - b[:, 1]], fmt='o', color='k', ms=3.5, lw=0.8)
        key = 'MW all' if lab.startswith('MW') else 'LMC all'
        o = res[(res['sample'] == key) & (res.kind == 'observed')].iloc[0]
        m = res[(res['sample'] == key) & (res.kind != 'observed')].iloc[0]
        ax.axvspan(o.h_lo, o.h_hi, color=col, alpha=0.15, lw=0)
        ax.axvline(o.h, color=col, lw=1.2, label=f'observed hinge {o.h:.2f} [{o.h_lo:.2f}, {o.h_hi:.2f}]')
        ax.axvline(m.h, color='0.3', ls='--', lw=1.0, label=f'model-predicted hinge {m.h:.2f}')
        ax.set_yscale('log'); ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
        ax.set_title(f'{lab}, N = {len(d)}', loc='left')
        ax.legend(loc='lower right', frameon=False, fontsize=6.2)
    axes[0].set_ylabel(r'$v_{\rm macro}$ [km s$^{-1}$]')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}/fig_lmc_vmac.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
    print('saved', f'{FIG}/fig_lmc_vmac.png')


if __name__ == '__main__':
    main()
