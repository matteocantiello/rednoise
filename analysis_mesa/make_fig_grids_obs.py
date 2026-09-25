#!/usr/bin/env python3
"""
v2 grids (MW, LMC, SMC; non-rotating, MS + post-MS) against every dataset we have, per galaxy.

fig_shrd_grids_obs   3 x 3 sHRD maps (columns MW, LMC, SMC). Background: median model quantity in (log Teff, log L_spec)
                     cells; grey = models without an FeCZ; black contour = FeCZ present in half the models (dashed: the
                     other two metallicities). Points: observations on the SAME colour scale.
                     row 1  v_c,max (FeCZ)                vs v_macro (IACOB + Holgado 2022; LMC Serebriakova 2024)
                     row 2  v_c,max (FeCZ)                vs v_mic (Markova 2025, MW); LMC/SMC: stars with SLF light curves
                     row 3  k_nu * nu_c (FeCZ turnover)   vs nu_char (Galactic Lorentzian, Bowman 2024 GP LMC/SMC,
                                                          Bowman 2019b LMC); k_nu fitted on MW only and applied to all Z
fig_profiles_Z       hot stars (log Teff > 4.3) vs log L_spec: (a) velocities; (b) alpha0 vs k_a M_t F_c/F;
                     (c) nu_char vs k_nu nu_c; (d) fraction of MS models with an FeCZ per Z, with the log L_spec
                     distribution of each galaxy's stars
fig_fecz_window      the Magellanic samples against the MW/LMC/SMC FeCZ boundaries (Van Daele correspondence item 1)
Writes figures_v2/*.png|pdf and data_v2/grids_obs_summary.csv.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import binned_statistic_2d
from scipy.interpolate import LinearNDInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lmc_vmac_test import load_mw, load_lmc, vz_read_tsv, star_key, RAW     # noqa: E402
import vmic_test                                                          # noqa: E402

OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
D2, FIG = f'{HERE}/data_v2', f'{HERE}/figures_v2'
LSUN = np.log10(5777.0 ** 4 / 27400.0)
ZS = ('MW', 'LMC', 'SMC')
LS = {'MW': '-', 'LMC': '--', 'SMC': ':'}
TB = np.arange(3.90, 4.76, 0.0125)
LB = np.arange(2.0, 4.66, 0.03)
plt.rcParams.update({'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5, 'legend.fontsize': 6,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5})


# ------------------------------------------------------------------------------------------------------------ models
def models(z, w='0.0'):
    fs = [f'{D2}/mesa_ms_{z}w{w}.csv', f'{D2}/mesa_post_{z}w{w}.csv']
    if not all(os.path.exists(f) for f in fs):
        return None
    ms, po = pd.read_csv(fs[0]), pd.read_csv(fs[1])
    ms['phase'] = 'MS'; po['phase'] = 'post'
    d = pd.concat([ms, po], ignore_index=True)
    d = d[np.isfinite(d.logTeff) & np.isfinite(d.logLspec)]
    d['aprox'] = d.FeCZ_mach_ahp * d.FeCZ_FcF                       # S5 / S3t amplitude predictor M_t F_c/F
    return d


_GC = {}


def _grid(d):
    """Cell centres, and a mask of cells within ~one track spacing of a model (scaled distance < 1.2)."""
    key = id(d)
    if key not in _GC:
        from scipy.spatial import cKDTree
        X, Y = np.meshgrid(0.5 * (TB[1:] + TB[:-1]), 0.5 * (LB[1:] + LB[:-1]), indexing='ij')
        pts = np.c_[d.logTeff / 0.03, d.logLspec / 0.08]
        P = np.c_[X.ravel() / 0.03, Y.ravel() / 0.08]
        near = (cKDTree(pts).query(P)[0] < 1.2).reshape(X.shape)
        pres = LinearNDInterpolator(pts, d.FeCZ_present.astype(float).values)(P).reshape(X.shape)
        _GC[key] = (X, Y, P, near, np.where(near, pres, np.nan))
    return _GC[key]


def cellmap(d, col, log=True, fecz_only=True):
    """Model quantity interpolated between tracks (FeCZ models only), masked off the grid and where FeCZ is absent."""
    X, Y, P, near, pres = _grid(d)
    v = d[col].values.astype(float)
    ok = np.isfinite(v) & ((v > 0) if log else True)
    if fecz_only:
        ok &= d.FeCZ_present.values.astype(bool)
    vv = np.log10(v[ok]) if log else v[ok]
    f = LinearNDInterpolator(np.c_[d.logTeff.values[ok] / 0.03, d.logLspec.values[ok] / 0.08], vv)
    out = f(P).reshape(X.shape)
    out = np.where(near & (np.nan_to_num(pres) >= 0.5), out, np.nan)
    return 10 ** out if log else out


def presence(d):
    return _grid(d)[4]


def coverage(d):
    return _grid(d)[3].astype(float)


def interp(d, col, lT, lL, fecz_zero=True):
    v = d[col].values.astype(float)
    ok = np.isfinite(v) & (v > 0)
    f = LinearNDInterpolator(np.c_[d.logTeff[ok] / 0.03, d.logLspec[ok] / 0.08], np.log10(v[ok]))
    out = 10 ** f(np.c_[np.asarray(lT) / 0.03, np.asarray(lL) / 0.08])
    if fecz_zero:
        p = LinearNDInterpolator(np.c_[d.logTeff / 0.03, d.logLspec / 0.08], d.FeCZ_present.astype(float))(
            np.c_[np.asarray(lT) / 0.03, np.asarray(lL) / 0.08])
        out = np.where(np.isfinite(p) & (p < 0.5), 0.0, out)
    return out


def ms_lines(d):
    m = d[(d.phase == 'MS') & d.ms_complete.astype(bool)].sort_values('tau').groupby('Minit')
    return m.first().sort_index(), m.last().sort_index()


# ------------------------------------------------------------------------------------------------------ observations
def obs():
    o = {}
    o['vmac_MW'] = load_mw()
    o['vmac_LMC'] = load_lmc()
    v = vmic_test.load()
    o['vmic_MW'] = v[(v.lT > 4.0) & (v.vmic > 0) & ~v.pulsating]
    E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
    E = E[np.isfinite(E.lT) & np.isfinite(E.lL)]
    o['rn_MW'] = E[E['sample'].isin(['Bowman2020_Galactic', 'Shen2024_Galactic'])]
    o['rn_LMC'] = E[E['sample'].isin(['Bowman2024_LMC', 'Bowman2019b_LMC'])]
    o['rn_SMC'] = E[E['sample'] == 'Bowman2024_SMC']
    pb = pd.read_csv(f'{HERE}/data_obs/pb25_stars.csv')
    o['cyg'] = pb[(pb['sample'] == 'PB25_CygOB') & pb.pb_ok & np.isfinite(pb.lL)]
    # BLOeM (Bestenlehner+2025) and the Van Daele PSF light-curve stars
    b = vz_read_tsv(f'{RAW}/bestenlehner25_a1.tsv')
    for c in ('Teff', 'logg', 'logL', '_RAJ2000', '_DEJ2000'):
        b[c] = pd.to_numeric(b[c], errors='coerce')
    b['lT'] = np.log10(b.Teff * 1e3); b['lL'] = 4 * b.lT - b.logg - LSUN
    o['bloem'] = b[np.isfinite(b.lL)]
    try:
        import vandaele_compare as vc
        import smc_instrumental_check as sc
        g = vc.placed(); g = g[g.bidx >= 0].merge(sc.classes(), on='gaia')
        g['lL'] = 4 * g.lT - pd.to_numeric(g.logg, errors='coerce') - LSUN
        o['vd'] = g[np.isfinite(g.lL)]
    except Exception as e:
        print('Van Daele placement unavailable:', e)
    return o


# ------------------------------------------------------------------------------------------------------------ figure 1
def fig_maps(M, o, kn):
    fig, ax = plt.subplots(3, 3, figsize=(10.5, 10.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.9, bottom=0.05, top=0.96, wspace=0.05, hspace=0.12)
    vnorm = LogNorm(1, 100); nnorm = LogNorm(0.1, 10)
    X, Y = np.meshgrid(0.5 * (TB[1:] + TB[:-1]), 0.5 * (LB[1:] + LB[:-1]), indexing='ij')
    for j, z in enumerate(ZS):
        d = M.get(z)
        for i in range(3):
            a = ax[i, j]
            if d is None:
                a.text(0.5, 0.5, f'{z}: no v2 extract', transform=a.transAxes, ha='center'); continue
            cov = coverage(d) > 0
            pres = presence(d)
            if i < 2:
                val = cellmap(d, 'FeCZ_vmax_kms'); norm, cmap = vnorm, 'viridis'
            else:
                val = cellmap(d, 'FeCZ_nuc_d') * kn; norm, cmap = nnorm, 'magma'
            grey = np.where(cov & ~np.isfinite(val), 1.0, np.nan)
            a.pcolormesh(TB, LB, grey.T, cmap='Greys', vmin=0, vmax=3, shading='flat', rasterized=True)
            im = a.pcolormesh(TB, LB, val.T, cmap=cmap, norm=norm, shading='flat', alpha=0.85, rasterized=True)
            for z2 in ZS:
                d2 = M.get(z2)
                if d2 is None:
                    continue
                p2 = presence(d2)
                a.contour(X, Y, np.nan_to_num(p2, nan=0.0), levels=[0.5], colors='k' if z2 == z else '0.35',
                          linewidths=1.1 if z2 == z else 0.6, linestyles=LS[z2])
            za, ta = ms_lines(d)
            a.plot(za.logTeff, za.logLspec, color='w', lw=0.8); a.plot(ta.logTeff, ta.logLspec, color='w', lw=0.5, ls='--')
            # observations
            kw = dict(s=11, edgecolors='k', linewidths=0.25, zorder=5)
            if i == 0 and f'vmac_{z}' in o:
                s = o[f'vmac_{z}']; sc_ = a.scatter(s.lT, s.lL, c=s.vmac, cmap=cmap, norm=norm, **kw)
                a.set_title(f'{z}: model v$_{{c,\\max}}$ vs v$_{{\\rm macro}}$ (N = {len(s)})', loc='left')
            elif i == 0:
                a.set_title(f'{z}: model v$_{{c,\\max}}$ (no v$_{{\\rm macro}}$ sample)', loc='left')
            if i == 1:
                if z == 'MW':
                    s = o['vmic_MW']; a.scatter(s.lT, s.lL, c=s.vmic, cmap=cmap, norm=norm, marker='D', **kw)
                    a.set_title(f'MW: model v$_{{c,\\max}}$ vs v$_{{\\rm mic}}$ (Markova+25, N = {len(s)})', loc='left')
                else:
                    s = o[f'rn_{z}']
                    a.scatter(s.lT, s.lL, s=12, facecolors='none', edgecolors='k', linewidths=0.6, zorder=5,
                              label=f'Bowman 2019b/2024 ({len(s)})')
                    if z == 'SMC' and 'vd' in o:
                        v = o['vd']
                        a.scatter(v[v.cls == 'SLF'].lT, v[v.cls == 'SLF'].lL, marker='^', s=16, c='#eb6834', edgecolors='k',
                                  linewidths=0.3, zorder=6, label='Van Daele: SLF')
                        a.scatter(v[v.cls != 'SLF'].lT, v[v.cls != 'SLF'].lL, marker='v', s=16, facecolors='w',
                                  edgecolors='#eb6834', linewidths=0.6, zorder=6, label='Van Daele: other/none')
                    a.legend(loc='lower left', frameon=False)
                    a.set_title(f'{z}: model v$_{{c,\\max}}$; stars with SLF light curves', loc='left')
            if i == 2:
                s = o[f'rn_{z}']; s = s[np.isfinite(s.nuchar) & (s.nuchar > 0)]
                a.scatter(s.lT, s.lL, c=s.nuchar, cmap=cmap, norm=norm, marker='o' if z == 'MW' else 's', **kw)
                if z == 'MW':
                    c = o['cyg']; a.scatter(c.lT, c.lL, c=c.nuchar_cal_d, cmap=cmap, norm=norm, marker='P', **kw)
                a.set_title(f'{z}: {kn:.1f}×model ν$_c$ vs ν$_{{\\rm char}}$ (N = {len(s)})', loc='left')
            a.set_xlim(4.75, 3.95); a.set_ylim(2.0, 4.6)
            if j == 0:
                a.set_ylabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
            if i == 2:
                a.set_xlabel(r'$\log T_{\rm eff}$')
        if j == 2:
            for i, (norm, cmap, lab) in enumerate(((vnorm, 'viridis', 'velocity [km s$^{-1}$]'),
                                                   (vnorm, 'viridis', 'velocity [km s$^{-1}$]'),
                                                   (nnorm, 'magma', r'$\nu$ [d$^{-1}$]'))):
                cax = fig.add_axes([0.91, ax[i, 2].get_position().y0, 0.012, ax[i, 2].get_position().height])
                fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, label=lab)
    ax[0, 0].plot([], [], 'k-', lw=1.1, label='FeCZ in 50% of models (this Z)')
    for z2 in ZS:
        ax[0, 0].plot([], [], color='0.35', ls=LS[z2], lw=0.6, label=f'{z2} boundary')
    ax[0, 0].plot([], [], 'w-', label='ZAMS / TAMS')
    ax[0, 0].legend(loc='lower left', frameon=True, fontsize=5.5)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}/fig_shrd_grids_obs.{ext}', dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------ figure 2
def binned(x, y, edges, nmin=5):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (x >= a) & (x < b) & np.isfinite(y)
        if k.sum() >= nmin:
            out.append((0.5 * (a + b), *np.percentile(y[k], [50, 16, 84]), k.sum()))
    return np.array(out) if out else np.zeros((0, 5))


def model_profile(d, col, edges, hot=4.3, ms_only=True, fecz_zero=True):
    s = d[(d.logTeff > hot) & ((d.phase == 'MS') if ms_only else True)]
    v = s[col].values.astype(float)
    if fecz_zero:
        v = np.where(s.FeCZ_present.astype(bool), v, np.nan)
    return binned(s.logLspec.values, np.log10(np.where(v > 0, v, np.nan)), edges, nmin=3)


def fig_profiles(M, o, kn, ka, rows):
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.3))
    fig.subplots_adjust(wspace=0.32)
    e = np.arange(2.3, 4.5, 0.15)
    C = {'MW': '#2a78d6', 'LMC': '#eb6834', 'SMC': '#1baf7a'}
    hot = lambda s: s[s.lT > 4.3]
    # (a) velocities
    a = ax[0]
    for key, lab, mk, col in (('vmac_MW', 'v$_{\\rm macro}$ MW', 'o', C['MW']), ('vmac_LMC', 'v$_{\\rm macro}$ LMC', 'o', C['LMC']),
                              ('vmic_MW', 'v$_{\\rm mic}$ MW', 'D', '0.3')):
        s = hot(o[key]); y = s.vmac if 'vmac' in key else s.vmic
        b = binned(s.lL.values, np.log10(y.values), e)
        a.errorbar(b[:, 0], 10 ** b[:, 1], yerr=[10 ** b[:, 1] - 10 ** b[:, 2], 10 ** b[:, 3] - 10 ** b[:, 1]],
                   fmt=mk + '-', ms=3, lw=0.9, color=col, label=lab, capsize=0)
    for z, d in M.items():
        b = model_profile(d, 'FeCZ_vmax_kms', e)
        a.plot(b[:, 0], 10 ** b[:, 1], ls=LS[z], color=C[z], lw=1.6, label=f'model v$_{{c,\\max}}$ {z}')
    a.set_yscale('log'); a.set_ylim(0.5, 200); a.set_xlabel(r'$\log\mathcal{L}$'); a.set_ylabel('km s$^{-1}$')
    a.legend(frameon=False, fontsize=5.3, ncol=2, loc='lower right'); a.set_title('(a) velocities, log T$_{\\rm eff}$ > 4.3', loc='left')
    # (b) alpha0
    a = ax[1]
    s = hot(o['rn_MW']); s = s[s.a0_ok & (s.alpha0 > 0)]
    b = binned(s.lL.values, np.log10(s.alpha0.values), e)
    a.errorbar(b[:, 0], 10 ** b[:, 1], yerr=[10 ** b[:, 1] - 10 ** b[:, 2], 10 ** b[:, 3] - 10 ** b[:, 1]], fmt='o-', ms=3,
               color=C['MW'], lw=0.9, label='α$_0$ MW (Bowman 2020 + Shen 2024)')
    c = hot(o['cyg'])
    b = binned(c.lL.values, np.log10(c.alpha0_cal_umag.values), e, nmin=4)
    a.errorbar(b[:, 0], 10 ** b[:, 1], yerr=[10 ** b[:, 1] - 10 ** b[:, 2], 10 ** b[:, 3] - 10 ** b[:, 1]], fmt='P-', ms=3,
               color='0.4', lw=0.8, label='α$_0$ Cyg OB (PB25, calibrated)')
    for z, d in M.items():
        bm = model_profile(d, 'aprox', e)
        a.plot(bm[:, 0], ka * 10 ** bm[:, 1], ls=LS[z], color=C[z], lw=1.6, label=f'{ka:.2g}×$\\mathcal{{M}}_t F_c/F$ {z}')
    a.set_yscale('log'); a.set_ylim(3, 2e4); a.set_xlabel(r'$\log\mathcal{L}$'); a.set_ylabel('α$_0$ [µmag]')
    a.legend(frameon=False, fontsize=5.3); a.set_title('(b) red-noise amplitude', loc='left')
    # (c) nu_char
    a = ax[2]
    for z, mk in (('MW', 'o'), ('LMC', 's'), ('SMC', 's')):
        s = hot(o[f'rn_{z}']); s = s[s.nuchar > 0]
        b = binned(s.lL.values, np.log10(s.nuchar.values), e, nmin=4)
        if len(b):
            a.errorbar(b[:, 0] + {'MW': 0, 'LMC': 0.02, 'SMC': 0.04}[z], 10 ** b[:, 1],
                       yerr=[10 ** b[:, 1] - 10 ** b[:, 2], 10 ** b[:, 3] - 10 ** b[:, 1]], fmt=mk, ms=3.5, color=C[z],
                       lw=0.9, label=f'ν$_{{\\rm char}}$ {z}' + (' (GP/Lor.)' if z != 'MW' else ''))
    for z, d in M.items():
        bm = model_profile(d, 'FeCZ_nuc_d', e)
        a.plot(bm[:, 0], kn * 10 ** bm[:, 1], ls=LS[z], color=C[z], lw=1.6, label=f'{kn:.1f}×ν$_c$ {z}')
    a.set_yscale('log'); a.set_xlabel(r'$\log\mathcal{L}$'); a.set_ylabel(r'$\nu$ [d$^{-1}$]')
    a.legend(frameon=False, fontsize=5.3, ncol=2); a.set_title('(c) characteristic frequency', loc='left')
    # (d) FeCZ presence and where the data are
    a = ax[3]
    for z, d in M.items():
        s = d[(d.phase == 'MS') & (d.logTeff > 4.3)]
        b = binned(s.logLspec.values, s.FeCZ_present.astype(float).values, np.arange(2.3, 4.5, 0.05), nmin=3)
        a.plot(b[:, 0], b[:, 1], ls=LS[z], color=C[z], lw=1.6, label=f'{z} MS models with FeCZ')
    a2 = a.twinx()
    for key, z in (('rn_MW', 'MW'), ('rn_LMC', 'LMC'), ('rn_SMC', 'SMC')):
        s = hot(o[key])
        a2.hist(s.lL, bins=np.arange(2.3, 4.5, 0.1), histtype='step', color=C[z], lw=0.9, alpha=0.8)
    if 'vd' in o:
        v = hot(o['vd'])
        a2.hist(v.lL, bins=np.arange(2.3, 4.5, 0.1), histtype='step', color=C['SMC'], lw=0.9, ls='--')
    a2.set_ylabel('hot red-noise stars (histograms; dashed: Van Daele)', fontsize=6)
    a.set_ylim(-0.02, 1.05); a.set_xlabel(r'$\log\mathcal{L}$'); a.set_ylabel('fraction of MS models with FeCZ')
    a.legend(frameon=False, fontsize=5.5, loc='upper left'); a.set_title('(d) where the FeCZ exists vs where the data are', loc='left')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}/fig_profiles_Z.{ext}', dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------ figure 3
def fig_window(M, o, rows):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
    fig.subplots_adjust(wspace=0.22)
    X, Y = np.meshgrid(0.5 * (TB[1:] + TB[:-1]), 0.5 * (LB[1:] + LB[:-1]), indexing='ij')
    C = {'MW': '#2a78d6', 'LMC': '#eb6834', 'SMC': '#1baf7a'}
    for k, (a, z0) in enumerate(zip(ax, ('SMC', 'LMC'))):
        d = M.get(z0)
        if d is not None:
            pres = presence(d); cov = coverage(d) > 0
            nofe = np.where(cov & (np.nan_to_num(pres) < 0.5), 1.0, np.nan)
            a.pcolormesh(TB, LB, nofe.T, cmap='Greys', vmin=0, vmax=4, shading='flat', rasterized=True)
            za, ta = ms_lines(d)
            a.plot(za.logTeff, za.logLspec, color='k', lw=0.8); a.plot(ta.logTeff, ta.logLspec, color='k', lw=0.5, ls='--')
        for z2, d2 in M.items():
            a.contour(X, Y, np.nan_to_num(presence(d2), nan=0.0), levels=[0.5], colors=C[z2], linewidths=1.3, linestyles=LS[z2])
            a.plot([], [], color=C[z2], ls=LS[z2], label=f'{z2}: FeCZ in 50% of models')
        if z0 == 'SMC':
            b = o['bloem']; a.scatter(b.lT, b.lL, s=3, c='0.6', lw=0, label=f'BLOeM (Bestenlehner+25, {len(b)})', zorder=2)
            if 'vd' in o:
                v = o['vd']
                a.scatter(v[v.cls == 'SLF'].lT, v[v.cls == 'SLF'].lL, marker='^', s=22, c=C['SMC'], edgecolors='k', lw=0.3,
                          zorder=5, label=f'Van Daele PSF: SLF ({(v.cls == "SLF").sum()})')
                a.scatter(v[v.cls != 'SLF'].lT, v[v.cls != 'SLF'].lL, marker='v', s=22, facecolors='w', edgecolors=C['SMC'],
                          lw=0.7, zorder=5, label=f'Van Daele PSF: other/none ({(v.cls != "SLF").sum()})')
            s = o['rn_SMC']; a.scatter(s.lT, s.lL, marker='s', s=20, c='#9b59b6', edgecolors='k', lw=0.3, zorder=6,
                                       label=f'Bowman+2024 SMC ({len(s)})')
        else:
            s = o['vmac_LMC']; a.scatter(s.lT, s.lL, s=6, c='0.55', lw=0, label=f'Serebriakova+24 v$_{{\\rm macro}}$ ({len(s)})', zorder=2)
            s = o['rn_LMC']; a.scatter(s.lT, s.lL, marker='s', s=20, c='#9b59b6', edgecolors='k', lw=0.3, zorder=6,
                                       label=f'Bowman 2019b/2024 LMC ({len(s)})')
        a.set_xlim(4.75, 3.95); a.set_ylim(2.0, 4.6)
        a.set_xlabel(r'$\log T_{\rm eff}$'); a.set_ylabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
        a.set_title(f'{z0}: grey = {z0} models without an FeCZ (v2, ω = 0)', loc='left')
        a.legend(loc='lower left', frameon=True, fontsize=5.6)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}/fig_fecz_window.{ext}', dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def window_fractions(M, o):
    """Fraction of each Magellanic sample sitting where its own-Z models have no FeCZ (MS + post-MS, omega = 0)."""
    rows = []
    for z, keys in (('SMC', ('bloem', 'vd', 'rn_SMC')), ('LMC', ('vmac_LMC', 'rn_LMC'))):
        d = M.get(z)
        if d is None:
            continue
        P = LinearNDInterpolator(np.c_[d.logTeff / 0.03, d.logLspec / 0.08], d.FeCZ_present.astype(float))
        for k in keys:
            if k not in o:
                continue
            s = o[k]
            p = P(np.c_[s.lT / 0.03, s.lL / 0.08])
            ok = np.isfinite(p)
            r = dict(Z=z, sample=k, n=len(s), n_on_grid=int(ok.sum()), frac_no_fecz=float((p[ok] < 0.5).mean()))
            if k == 'vd':
                for cl in ('SLF', 'none'):
                    m = ok & (s.cls == cl).values
                    r[f'frac_no_fecz_{cl}'] = float((p[m] < 0.5).mean()) if m.any() else np.nan
                    r[f'n_{cl}'] = int(m.sum())
            rows.append(r)
    return rows


def main():
    M = {z: models(z) for z in ZS}
    M = {z: d for z, d in M.items() if d is not None}
    o = obs()
    mw = M['MW']
    # scale factors fitted on the MW (primary sample, where the MW models have an FeCZ); applied unchanged to LMC/SMC
    s = o['rn_MW']; s = s[(s.nuchar > 0)]
    nuc = interp(mw, 'FeCZ_nuc_d', s.lT.values, s.lL.values)
    ok = np.isfinite(nuc) & (nuc > 0)
    kn = 10 ** np.median(np.log10(s.nuchar.values[ok] / nuc[ok]))
    s2 = o['rn_MW']; s2 = s2[s2.a0_ok & (s2.alpha0 > 0)]
    ap = interp(mw, 'aprox', s2.lT.values, s2.lL.values)
    ok2 = np.isfinite(ap) & (ap > 0)
    ka = 10 ** np.median(np.log10(s2.alpha0.values[ok2] / ap[ok2]))
    rows = [dict(Z='MW', sample='scale', k_nu=kn, k_alpha=ka, n_nu=int(ok.sum()), n_alpha=int(ok2.sum()))]
    # per-galaxy offsets of the observations from the MW-scaled models, at the stars' own positions
    for z in M:
        for key, col, mcol, k in (('rn_' + z, 'nuchar', 'FeCZ_nuc_d', kn), ('vmac_' + z, 'vmac', 'FeCZ_vmax_kms', 1.0)):
            if key not in o:
                continue
            ss = o[key]; ss = ss[(ss[col] > 0) & (ss.lT > 4.3)]
            pr = interp(M[z], mcol, ss.lT.values, ss.lL.values) * k
            g = np.isfinite(pr) & (pr > 0)
            rows.append(dict(Z=z, sample=key, quantity=col, n=int(g.sum()),
                             median_log_obs_over_model=float(np.median(np.log10(ss[col].values[g] / pr[g]))) if g.any() else np.nan,
                             n_no_fecz=int((pr == 0).sum())))
    rows += window_fractions(M, o)
    fig_maps(M, o, kn)
    fig_profiles(M, o, kn, ka, rows)
    fig_window(M, o, rows)
    r = pd.DataFrame(rows)
    r.to_csv(f'{D2}/grids_obs_summary.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 220)
    print(r.to_string(index=False, float_format=lambda v: f'{v:.3f}'))


if __name__ == '__main__':
    main()
