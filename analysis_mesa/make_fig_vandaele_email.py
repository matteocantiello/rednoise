#!/usr/bin/env python3
"""
Figures for the note to P. Van Daele (Vandaele_correspondence/): what we did with the BLOeM PSF light curves, and the issues.
Fig. 1 (reproduction): our per-sector refit vs values digitised from their Fig. 13 (nu_char, alpha0/C_w); our fit on identical
        data with two fit metrics; our Galactic SPOC fits vs Shen 2024 on the same sectors.
Fig. 2 (noise floor): nu_char and alpha0 of their 'no significant variability' stars vs their SLF stars (our fits);
        example periodograms; residuals from our Galactic nu_char plane.
Writes Vandaele_correspondence/fig{1,2}_*.png.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vandaele_fit as vf                 # noqa: E402
import mw_smc_nuchar as ms                # noqa: E402
import smc_instrumental_check as sc       # noqa: E402
OUT = '/mnt/home/mcantiello/work/rednoise/Vandaele_correspondence'
C1, C2, C3 = '#2a78d6', '#eb6834', '#555555'
plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.7})


def one2one(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], color='0.6', lw=0.8, ls='--', zorder=0)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')


def stats(x, y):
    d = y - x
    return f'N = {len(d)}, median Δ = {np.median(d):+.2f}, MAD = {np.median(abs(d - np.median(d))):.2f} dex'


def fig1():
    c = pd.read_csv(f'{HERE}/data_obs/vandaele_compare.csv')
    c = c[(c.metric == 'lin') & c.prewhitened]
    f = pd.read_csv(f'{HERE}/data_obs/vandaele_fits.csv')
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.4))
    fig.subplots_adjust(wspace=0.62)
    for a, q, lab, lim in ((ax[0], 'nu', r'$\log\nu_{\rm char}$ [d$^{-1}$]', (-1.2, 0.4)),
                           (ax[1], 'a', r'$\log\alpha_0/C_{\rm w}$', (0.0, 2.6))):
        s = c[c.qty == q]
        a.scatter(s.theirs, s.ours, s=10, c=C1, lw=0)
        one2one(a, *lim)
        a.set_xlabel(f'Van Daele+2026 (digitised Fig. 13)\n{lab}'); a.set_ylabel(f'our refit, {lab}')
        a.set_title(stats(s.theirs, s.ours), fontsize=6.5, loc='left')
    a = ax[2]
    p = f[f.prewhitened].pivot_table(index=['gaia', 'sector'], columns='metric', values='nuchar')
    p = p[(p.lin >= 0.15) & (p.log >= 0.15)]
    a.scatter(np.log10(p.lin), np.log10(p.log), s=8, c=C3, lw=0)
    one2one(a, -1.2, 1.0)
    a.set_xlabel(r'our fit, residuals in $A$: $\log\nu_{\rm char}$'); a.set_ylabel(r'our fit, residuals in $\log A$: $\log\nu_{\rm char}$')
    a.set_title('same SMC light curves, two fit metrics\n' + stats(np.log10(p.lin), np.log10(p.log)), fontsize=6.5, loc='left')
    a = ax[3]
    fmw = pd.read_csv(f'{HERE}/data_obs/vandaele_fits_mw.csv')
    m, _ = ms.names()
    g = ms.good(fmw, 'lin').copy(); g['Name'] = g.id.map(m)
    s2 = ms.vz_read_tsv(f'{ms.RAW}/shen24_t2.tsv')
    s2['Name'] = s2.Name.str.strip(); s2['sector'] = pd.to_numeric(s2.Sec, errors='coerce')
    s2['nu_pub'] = pd.to_numeric(s2['nu-char'], errors='coerce')
    x = g.merge(s2[['Name', 'sector', 'nu_pub']], on=['Name', 'sector'])
    a.scatter(np.log10(x.nu_pub), np.log10(x.nuchar), s=8, c=C2, lw=0)
    one2one(a, -1.2, 1.2)
    a.set_xlabel(r'Shen+2024 (same sector): $\log\nu_{\rm char}$'); a.set_ylabel(r'our fit, Galactic SPOC: $\log\nu_{\rm char}$')
    a.set_title('Galactic check of our pipeline\n' + stats(np.log10(x.nu_pub), np.log10(x.nuchar)), fontsize=6.5, loc='left')
    for i, a in enumerate(ax):
        a.text(-0.28, 1.02, 'abcd'[i], transform=a.transAxes, fontweight='bold', fontsize=10)
    fig.savefig(f'{OUT}/fig1_reproduction.png', dpi=200, bbox_inches='tight', facecolor='white')


def fig2():
    cl = sc.classes()
    f = pd.read_csv(f'{HERE}/data_obs/vandaele_fits.csv').merge(cl, on='gaia')
    g = f[(f.metric == 'lin') & f.prewhitened]
    ok = (g.nuchar >= 0.15) & (g.dBIC_vs_white > 0)
    fig, ax = plt.subplots(1, 4, figsize=(12, 3.3))
    fig.subplots_adjust(wspace=0.5)
    a = ax[0]
    bins = np.linspace(-1.0, 0.8, 19)
    for k, col, lab in (('SLF', C1, 'their "SLF" stars'), ('none', C2, 'their "no significant variability" stars')):
        s = g[ok & (g.cls == k)]
        a.hist(np.log10(s.nuchar), bins=bins, histtype='step', lw=1.6, color=col,
               label=f'{lab}\n({len(s)}/{(g.cls == k).sum()} sectors pass our cuts)')
    a.set_xlabel(r'our fit: $\log\nu_{\rm char}$ [d$^{-1}$]'); a.set_ylabel('sectors'); a.set_ylim(0, a.get_ylim()[1] * 1.45); a.legend(frameon=False, fontsize=6, loc='upper left')
    a = ax[1]
    for k, col in (('SLF', C1), ('none', C2)):
        s = g[ok & (g.cls == k)]
        a.scatter(s.Tmag, np.log10(s.alpha0), s=9, c=col, lw=0, alpha=0.8)
    a.set_xlabel('Tmag'); a.set_ylabel(r'our fit: $\log\alpha_0$ [e$^-$ s$^{-1}$]')
    a.set_title('amplitude of the fitted red noise', fontsize=7, loc='left')
    # example periodograms: a 'none' and an SLF star of similar Tmag, same sector
    a = ax[2]
    s = g[ok & (g.cls == 'none')]
    ex_n = s.iloc[int(np.argmin(abs(np.log10(s.nuchar.values) - np.log10(s.nuchar).median())))]
    s2 = g[ok & (g.cls == 'SLF') & (g.sector == ex_n.sector)]
    ex_s = s2.iloc[int(np.argmin(abs(s2.Tmag.values - ex_n.Tmag)))]
    for ex, col, lab, off in ((ex_n, C2, 'no variability', 1.0), (ex_s, C1, 'SLF', 10.0)):
        path = f'{vf.LC}/GAIA DR3 {int(ex.gaia)}_sector{int(ex.sector)}_PSFlc.txt'
        t, mm = vf.load(path)
        mm, _ = vf.prewhiten(t, mm, 'lin')
        nu, A = vf.amp_spectrum(t, mm)
        a.loglog(nu, A * off, color=col, lw=0.4, alpha=0.7)
        a.loglog(nu, vf.slf(nu, ex.alpha0, ex.nuchar, ex.gamma, ex.Cw) * off, color='k', lw=1.0)
        a.text(0.97, 0.14 if lab == 'SLF' else 0.06, f'BLOeM {ex.BLOeM} ({lab}), T = {ex.Tmag:.1f}, '
               rf'$\nu_{{\rm char}}$ = {ex.nuchar:.2f} d$^{{-1}}$', fontsize=5.8, color=col, transform=a.transAxes, ha='right', va='bottom')
    a.set_xlabel(r'frequency [d$^{-1}$]'); a.set_ylabel(r'amplitude [e$^-$ s$^{-1}$] (SLF star ×10)')
    a.set_title(f'sector {int(ex_n.sector)}; our fit in black', fontsize=7, loc='left'); a.set_ylim(bottom=a.get_ylim()[0] / 30)
    a = ax[3]
    fmw = pd.read_csv(f'{HERE}/data_obs/vandaele_fits_mw.csv')
    mw = ms.place_mw(ms.per_star(fmw, 'lin')); mw = mw[np.isfinite(mw.lT) & np.isfinite(mw.lL)]
    A_ = np.c_[np.ones(len(mw)), mw.lL, mw.lT]; cc = np.linalg.lstsq(A_, mw.lnu, rcond=None)[0]
    rmw = mw.lnu - A_ @ cc
    smc = ms.place_smc(ms.per_star(pd.read_csv(f'{HERE}/data_obs/vandaele_fits.csv'), 'lin'))
    smc['gaia'] = smc.id.str.replace('GAIA DR3 ', '').astype(int)
    smc = smc.merge(cl, on='gaia'); smc = smc[np.isfinite(smc.lT) & np.isfinite(smc.lL)]
    smc['res'] = smc.lnu - np.c_[np.ones(len(smc)), smc.lL, smc.lT] @ cc
    a.scatter(mw.lL, rmw, s=6, c='0.7', lw=0, label='Galactic (SPOC 2-min), our fit')
    for k, col, lab in (('SLF', C1, 'SMC "SLF"'), ('none', C2, 'SMC "no variability"')):
        s = smc[smc.cls == k]
        a.scatter(s.lL, s.res, s=12, c=col, lw=0, label=f'{lab}: median {s.res.median():+.2f}')
    a.axhline(0, color='k', lw=0.6)
    a.set_xlabel(r'$\log\mathcal{L}/\mathcal{L}_\odot$'); a.set_ylabel(r'$\log\nu_{\rm char}$ − Galactic plane($\mathcal{L}$, $T_{\rm eff}$)')
    a.legend(frameon=False, fontsize=5.8, loc='lower left')
    for i, a in enumerate(ax):
        a.text(-0.28, 1.02, 'abcd'[i], transform=a.transAxes, fontweight='bold', fontsize=10)
    fig.savefig(f'{OUT}/fig2_noise_floor.png', dpi=200, bbox_inches='tight', facecolor='white')



def fig3():
    """Noise-floor model: injection-recovery maps, and what Galactic-like SLF would look like in the SMC light curves."""
    s = pd.read_csv(f'{HERE}/data_obs/floor_injection_summary.csv')
    e = pd.read_csv(f'{HERE}/data_obs/floor_expected_detection.csv')
    x = pd.read_csv(f'{HERE}/data_obs/floor_expected_R.csv')
    fig, ax = plt.subplots(1, 4, figsize=(12, 3.3))
    fig.subplots_adjust(wspace=0.5)
    cols = {0.3: '#9b59b6', 0.7: C2, 1.5: C1, 3.0: '#2ca02c'}
    for nu_in, g in s.groupby('nu_in'):
        ax[0].plot(g.R, g.det_frac, 'o-', ms=3, lw=1, color=cols[nu_in], label=rf'$\nu_{{\rm in}}$ = {nu_in} d$^{{-1}}$')
        ax[1].plot(g.R, g.bias_naive, 'o-', ms=3, lw=1, color=cols[nu_in])
        ax[1].plot(g.R, g.bias_floor_det, 's--', ms=3, lw=0.8, color=cols[nu_in], alpha=0.7)
    ax[0].axhline(s[s.R == 0.1].det_frac.min(), color='0.6', ls=':', lw=0.8)
    ax[0].set_xscale('log'); ax[0].set_xlabel(r'injected $\alpha_0$ / floor amplitude ($R$)'); ax[0].set_ylabel('detected above floor (ΔBIC > 10)')
    ax[0].legend(frameon=False, fontsize=6); ax[0].set_title('injections into the 104 "no variability" sectors', fontsize=6.5, loc='left')
    ax[1].axhline(0, color='k', lw=0.6); ax[1].set_xscale('log'); ax[1].set_ylim(-1, 0.7)
    ax[1].set_xlabel(r'$R$'); ax[1].set_ylabel(r'recovered − injected $\log\nu_{\rm char}$')
    ax[1].set_title('solid: naive fit (SLF + white)\ndashed: SLF + floor + white (detected only)', fontsize=6.5, loc='left')
    a = ax[2]
    a.hist(x.logR_exp, bins=np.linspace(-1.5, 1.5, 16), color='0.6', histtype='stepfilled', alpha=0.6)
    a.axvline(0, color='k', lw=0.6)
    a.set_xlabel(r'expected $\log R$ for Galactic-like SLF'); a.set_ylabel('SMC star-sectors')
    a.set_title('Galactic α₀/flux (our SPOC fits) vs SMC floor', fontsize=6.5, loc='left')
    a = ax[3]
    obs_det, obs_nu = e.observed_det.iloc[0], np.log10(x.nu_slf).median()
    for dnu, mk in ((0.0, 'o'), (-0.2, 's'), (-0.4, '^'), (-0.6, 'D')):
        g = e[e.shift_log_nu == dnu]
        a.plot(g.expected_det, g.median_lnu_naive_expected, mk + '-', ms=4, lw=0.8, color=C1, alpha=0.4 + 0.15 * abs(dnu) / 0.2,
               label=rf'$\Delta\log\nu$ = {dnu:+.1f} (amplitude ×1 → ×0.1)')
    a.plot([obs_det], [obs_nu], '*', ms=12, color=C2, label='observed (SMC, 85 sectors)')
    a.set_xlabel('fraction detected above floor'); a.set_ylabel(r'median naive $\log\nu_{\rm char}$')
    a.legend(frameon=False, fontsize=5.5, loc='lower right')
    a.set_title('forward model: Galactic-like SLF + SMC floor', fontsize=6.5, loc='left')
    for i, a in enumerate(ax):
        a.text(-0.28, 1.02, 'abcd'[i], transform=a.transAxes, fontweight='bold', fontsize=10)
    fig.savefig(f'{OUT}/fig3_floor_model.png', dpi=200, bbox_inches='tight', facecolor='white')


if __name__ == '__main__':
    fig1(); fig2(); fig3()
    print('saved', OUT)
