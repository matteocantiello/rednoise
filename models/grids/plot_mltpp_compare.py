#!/usr/bin/env python3
"""
Differential test of MLT++ on FeCZ properties: MW/w0.0 (MLT++ off, production)
vs MW_mltpp/w0.0 (okay_to_reduce_gradT_excess = .true.), otherwise identical physics.

Tracks are compared on the main sequence at matched central hydrogen,
  s = 1 - X_c / X_c,ZAMS   (0 at ZAMS, 1 at TAMS),
which is defined for incomplete (running / failed) tracks too.

Outputs (figures/):
  mltpp_hrd.png          HRD + sHRD, both runs overlaid
  mltpp_fecz_tracks.png  FeCZ quantities vs s for selected masses
  mltpp_fecz_bias.png    MS-median and range of dlog Q = log Q(MLT++) - log Q(std) vs mass
  mltpp_fecz_sHRD.png    dlog Q mapped on the sHRD
and mltpp_fecz_bias.csv with the numbers of the bias figure.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
import grid_io as g

FIGDIR = f'{g.GRID}/figures'
os.makedirs(FIGDIR, exist_ok=True)
STD, MLT = ('MW', 'w0.0'), ('MW_mltpp', 'w0.0')

QUANT = [
    ('v_FeCZ_max',          r'$v_{\rm conv}^{\rm max}$ [cm/s]'),
    ('mach_FeCZ_max',       r'Mach$_{\rm max}$'),
    ('turnover_FeCZ',       r'$\tau_{\rm conv}$ [s]'),
    ('FeCZ_Fc_max',         r'$F_{c,\rm max}$ [erg/cm$^2$/s]'),
    ('b_FeCZ_max',          r'$B_{\rm eq}^{\rm max}$ [G]'),
    ('FeCZ_B_shutoff_conv', r'$B_{\rm shutoff}$ [G]'),
]
SHOW_MASSES = [10, 15, 20, 25, 30, 40, 50, 60]
S_GRID = np.linspace(0.0, 0.99, 100)
norm = LogNorm(8, 120)
cmap = plt.get_cmap('plasma_r')

std = g.load_subgrid(*STD)
mlt = g.load_subgrid(*MLT)


def ms_part(d):
    """Main-sequence slice and s = 1 - Xc/Xc,ZAMS; None if the track never reached ZAMS."""
    h = d['center_h1']
    i0 = g.find_zams(h)
    if i0 == 0 or h[i0:].min() > h[i0] - 0.01:
        return None  # still on the pre-MS (e.g. MLT++ M80-M120, dt collapsed there)
    it = g.find_tams(h)
    sl = slice(i0, it if it is not None else len(h))
    s = 1 - h[sl] / h[i0]
    return sl, s


def on_s(d, q):
    """Quantity q (log10) interpolated onto S_GRID over the MS; NaN outside the covered range."""
    r = ms_part(d)
    if r is None or q not in d:
        return np.full(len(S_GRID), np.nan)
    sl, s = r
    v = d[q][sl].astype(float)
    lv = np.where(np.isfinite(v) & (v > 0), np.log10(np.clip(v, 1e-300, None)), np.nan)
    ok = np.isfinite(lv)
    if ok.sum() < 5:
        return np.full(len(S_GRID), np.nan)
    order = np.argsort(s[ok])
    return np.interp(S_GRID, s[ok][order], lv[ok][order], left=np.nan, right=np.nan)


common = [m for m in g.MASSES if m in std and m in mlt
          and ms_part(std[m]) is not None and ms_part(mlt[m]) is not None]
print('masses with MS in both runs:', common)
skipped = [m for m in g.MASSES if m in mlt and ms_part(mlt[m]) is None]
print('MLT++ tracks without a MS yet:', skipped)

# ---------------------------------------------------------------- 1. HRD / sHRD
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
for ax, kind in zip(axes, ('HRD', 'sHRD')):
    for m in g.MASSES:
        for tr, ls, lw in ((std, '-', 0.9), (mlt, '--', 0.9)):
            if m not in tr or ms_part(tr[m]) is None:
                continue  # skip tracks still on the pre-MS
            d = tr[m]
            i0 = g.find_zams(d['center_h1'])
            y = d['log_L'] if kind == 'HRD' else g.spec_ell(d['log_Teff'], d['log_g'])
            ax.plot(d['log_Teff'][i0:], y[i0:], ls=ls, lw=lw, color=cmap(norm(max(m, 8))), alpha=0.85)
    ax.set_xlim(4.8, 3.45)
    ax.set_xlabel(r'$\log T_{\rm eff}$ [K]')
    ax.set_ylabel(r'$\log L/L_\odot$' if kind == 'HRD' else r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
    ax.grid(alpha=0.2, ls='--')
axes[0].legend(handles=[Line2D([], [], color='k', ls='-', label='MLT++ off (production)'),
                        Line2D([], [], color='k', ls='--', label='MLT++ on')], loc='lower left')
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cb = fig.colorbar(sm, ax=axes, pad=0.01, fraction=0.03)
cb.set_label(r'Initial mass [$M_\odot$]')
cb.set_ticks([8, 10, 15, 20, 30, 50, 80, 120]); cb.set_ticklabels(['≤8', '10', '15', '20', '30', '50', '80', '120'])
fig.suptitle(f'MW $Z=0.014$, non-rotating: MLT++ on vs off (MLT++ {len(mlt)} tracks, '
             f'{len(skipped)} still pre-ZAMS, not shown: ' + ', '.join(f'{m:g}' for m in skipped) + r' $M_\odot$)', fontsize=12, x=0.45)
fig.savefig(f'{FIGDIR}/mltpp_hrd.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)

# ---------------------------------------------------------------- 2. tracks vs s
masses = [m for m in SHOW_MASSES if m in common]
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), sharex=True)
for ax, (q, lab) in zip(axes.flat, QUANT):
    for m in masses:
        c = cmap(norm(m))
        ax.plot(S_GRID, on_s(std[m], q), color=c, lw=1.3, ls='-')
        ax.plot(S_GRID, on_s(mlt[m], q), color=c, lw=1.3, ls='--')
    ax.set_ylabel(r'$\log$ ' + lab)
    ax.grid(alpha=0.2, ls='--')
for ax in axes[1]:
    ax.set_xlabel(r'$1 - X_c/X_{c,\rm ZAMS}$  (ZAMS → TAMS)')
handles = [Line2D([], [], color=cmap(norm(m)), lw=2, label=f'{m:g} $M_\\odot$') for m in masses]
handles += [Line2D([], [], color='k', ls='-', label='MLT++ off'), Line2D([], [], color='k', ls='--', label='MLT++ on')]
fig.legend(handles=handles, loc='lower center', ncol=len(handles), fontsize=9, frameon=False,
           bbox_to_anchor=(0.5, -0.03))
fig.suptitle('FeCZ properties on the main sequence, MW non-rotating: MLT++ off (solid) vs on (dashed)',
             fontsize=12)
plt.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(f'{FIGDIR}/mltpp_fecz_tracks.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)

# ---------------------------------------------------------------- 3. bias vs mass
rows = []
fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for ax, (q, lab) in zip(axes.flat, QUANT):
    med, lo, hi, late = [], [], [], []
    for m in common:
        dl = on_s(mlt[m], q) - on_s(std[m], q)
        ok = np.isfinite(dl)
        if ok.sum() < 10:
            med.append(np.nan); lo.append(np.nan); hi.append(np.nan); late.append(np.nan)
            continue
        med.append(np.median(dl[ok])); lo.append(np.percentile(dl[ok], 5)); hi.append(np.percentile(dl[ok], 95))
        lt = ok & (S_GRID > 0.8)
        late.append(np.median(dl[lt]) if lt.sum() else np.nan)
        rows.append((q, m, med[-1], lo[-1], hi[-1], late[-1], ok.mean()))
    ax.axhline(0, color='0.5', lw=0.8)
    ax.fill_between(common, lo, hi, color='#2a78d6', alpha=0.2, lw=0, label='5–95% over MS')
    ax.plot(common, med, 'o-', color='#2a78d6', ms=4, lw=1.5, label='MS median')
    ax.plot(common, late, 's--', color='#eb6834', ms=3.5, lw=1.0, label=r'median for $s>0.8$ (late MS)')
    ax.set_xscale('log')
    ax.set_title(lab, fontsize=11)
    ax.set_ylabel(r'$\Delta\log$ (MLT++ − std)')
    ax.grid(alpha=0.2, ls='--')
for ax in axes[1]:
    ax.set_xlabel(r'Initial mass [$M_\odot$]')
    ax.set_xticks([5, 10, 20, 50]); ax.set_xticklabels(['5', '10', '20', '50'])
axes[0, 0].legend(fontsize=8, loc='best')
fig.suptitle('MLT++ bias on FeCZ properties along the main sequence (MW, non-rotating)', fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f'{FIGDIR}/mltpp_fecz_bias.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)
with open(f'{FIGDIR}/mltpp_fecz_bias.csv', 'w') as f:
    f.write('quantity,mass,dlog_median,dlog_p05,dlog_p95,dlog_median_late_MS,frac_MS_both_have_FeCZ\n')
    for r in rows:
        f.write(','.join(str(x) if isinstance(x, str) else f'{x:.4g}' for x in r) + '\n')

# ---------------------------------------------------------------- 4. dlog Q on the sHRD
# Pair models at the same (mass, s); place the difference at the MLT-off position.
fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
GT, GL = np.meshgrid(np.linspace(4.1, 4.75, 160), np.linspace(2.4, 4.6, 160))
for ax, (q, lab) in zip(axes.flat, QUANT):
    X, Y, D = [], [], []
    for m in common:
        sl, s = ms_part(std[m])
        d = std[m]
        o = np.argsort(s)
        x = np.interp(S_GRID, s[o], d['log_Teff'][sl][o])
        y = np.interp(S_GRID, s[o], g.spec_ell(d['log_Teff'], d['log_g'])[sl][o])
        dl = on_s(mlt[m], q) - on_s(std[m], q)
        X.append(x); Y.append(y); D.append(dl)
    X, Y, D = map(np.concatenate, (X, Y, D))
    ok = np.isfinite(D)
    zi = griddata((X[ok], Y[ok]), D[ok], (GT, GL), method='linear')
    lim = max(0.05, np.nanpercentile(np.abs(D[ok]), 98))
    cf = ax.contourf(GT, GL, zi, levels=np.linspace(-lim, lim, 21), cmap='RdBu_r',
                     norm=TwoSlopeNorm(0, -lim, lim), extend='both')
    for m in common:
        sl, _ = ms_part(std[m])
        d = std[m]
        ax.plot(d['log_Teff'][sl], g.spec_ell(d['log_Teff'], d['log_g'])[sl], color='0.3', lw=0.35)
        if m in (10, 15, 20, 30, 50, 60):
            ax.text(d['log_Teff'][sl][0], g.spec_ell(d['log_Teff'][sl][0], d['log_g'][sl][0]) - 0.05,
                    f'{m:g}', fontsize=6, color='0.3', ha='center', va='top')
    cb = fig.colorbar(cf, ax=ax, pad=0.02, fraction=0.05)
    cb.set_label(r'$\Delta\log$ (MLT++ − std)', fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_title(lab, fontsize=11)
    ax.set_xlim(4.75, 4.1); ax.set_ylim(2.4, 4.6)
for ax in axes[1]:
    ax.set_xlabel(r'$\log T_{\rm eff}$ [K]')
for ax in axes[:, 0]:
    ax.set_ylabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
fig.suptitle('MLT++ effect on FeCZ properties mapped on the sHRD (main sequence; red = MLT++ higher)',
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f'{FIGDIR}/mltpp_fecz_sHRD.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved mltpp_* figures to', FIGDIR)
