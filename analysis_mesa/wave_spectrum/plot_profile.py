"""Diagnostic figures for one model profile: where waves are excited and damped, and what sets them.

Usage (from analysis_mesa/):
    python3 -m wave_spectrum.plot_profile [--run RUN_DIR] [--xc 0.5] [--tag NAME]
Writes figures/wave_envelope_<tag>.{png,pdf} and figures/propagation_<tag>.{png,pdf}.

Figure 1 (envelope, x = log T, surface on the right): opacity; nabla_rad vs nabla_ad; v_conv and c_s;
H_P; fluxes (convective, launched wave flux, and wave flux left after quasi-adiabatic damping);
envelope propagation diagram (N, S_1, nu_c, non-adiabatic frequency); cumulative tau; q = K k_r^2/omega.
Figure 2 (whole star, x = r/R): classic propagation diagram with the g-mode cavity.

Wave source: source.Source defaults (a = 6.5, p = 1), with L_w = M_t F_c/F L from the history row.
Damping: quasi-adiabatic (run_star_extras); the q panel shows where that approximation fails.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .mesa_io import Track
from . import damping as D
from .source import Source

SECDAY = 86400.0
# reference palette (dataviz skill): categorical slots 1-4 in fixed order; validated light mode
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
LS = ['-', '--', '-.', ':']                       # secondary encoding (aqua/yellow < 3:1 contrast)
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
ZONE = {'FeCZ': '#cde2fb', 'HeII': '#f3e3c3', 'HeI': '#f3e3c3', 'HI': '#f3e3c3', 'UNKNOWN': '#ebebe7'}
T_TYPES = [('HI', 3e3, 1.1e4), ('HeI', 1.1e4, 3.5e4), ('HeII', 3.5e4, 1e5), ('FeCZ', 1e5, 5e5)]

plt.rcParams.update({
    'font.size': 9, 'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'xtick.color': MUTED,
    'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.6,
    'legend.frameon': False, 'legend.fontsize': 8, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb',
})


def zone_type(T_bottom):
    for name, lo, hi in T_TYPES:
        if lo < T_bottom < hi:
            return name
    return 'UNKNOWN'


def label_end(ax, x, y, text, color, dx=0.0, where='left'):
    """Direct label at the end of a line (text in secondary ink; a short colored bar carries identity)."""
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.any():
        return
    i = np.nonzero(ok)[0][0 if where == 'left' else -1]
    ax.annotate(text, (x[i], y[i]), xytext=(4 if where == 'left' else -4, 0), textcoords='offset points',
                ha='left' if where == 'left' else 'right', va='center', fontsize=7.5, color=INK2,
                bbox=dict(boxstyle='square,pad=0.1', fc='#fcfcfb', ec='none', alpha=0.8))


def label_at(ax, x, y, x0, text, dy=4, va='bottom'):
    """Direct label placed on the curve at x = x0 (secondary ink, on a surface-colored backing)."""
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.any():
        return
    i = np.argmin(np.where(ok, np.abs(x - x0), np.inf))
    ax.annotate(text, (x[i], y[i]), xytext=(0, dy), textcoords='offset points', ha='center', va=va,
                fontsize=7.5, color=INK2, bbox=dict(boxstyle='square,pad=0.1', fc='#fcfcfb', ec='none', alpha=0.85))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default='/home/mcantiello/rednoise_tests/T1a_M20_w00')
    ap.add_argument('--xc', type=float, default=0.5)
    ap.add_argument('--tag', default=None)
    ap.add_argument('--outdir', default='figures')
    args = ap.parse_args()

    tr = Track(args.run)
    p, h = min(tr.pairs(), key=lambda x: abs(x[1]['center_h1'] - args.xc))
    M = h['star_mass']
    Mi = p.header.get('initial_mass', M)
    tag = args.tag or f"M{Mi:g}_Xc{h['center_h1']:.2f}"
    kt, kb = p.fecz()
    path = D.path_arrays(p, kt)
    src = Source.from_history(h)
    logT = np.log10(p.T)
    L = 10 ** h['log_L'] * 3.828e33
    nu_c = src.omega_c / (2 * np.pi) * SECDAY
    title = (f"{Mi:g} M$_\\odot$ (now {M:.2f}), X$_c$ = {h['center_h1']:.2f}, log T$_{{\\rm eff}}$ = {np.log10(p.Teff):.3f}, "
             f"model {p.model_number}  —  {os.path.basename(args.run.rstrip('/'))}")

    # x range: from a bit below the FeCZ bottom to just above the photosphere
    x_lo = np.log10(p.T[p.k_phot]) - 0.08
    x_hi = logT[kb] + 0.35
    env = (logT >= x_lo) & (logT <= x_hi)

    # ---------------- derived quantities on the full profile ----------------
    N2 = p.N2
    N_nu = np.where(N2 > 0, np.sqrt(np.clip(N2, 0, None)) / (2 * np.pi) * SECDAY, np.nan)
    S1_nu = np.sqrt(2.0) * p.cs / p.r / (2 * np.pi) * SECDAY
    S2_nu = np.sqrt(6.0) * p.cs / p.r / (2 * np.pi) * SECDAY
    hd, dd = {}, None
    from .mesa_io import read_mesa_table
    _, dd = read_mesa_table(p.path)
    grada, gradr = dd['grada'], dd['gradr']
    vconv = 10 ** dd['log_conv_vel'] / 1e5
    vconv = np.where((p.mlt_type == 1) & (vconv > 1e-3), vconv, np.nan)

    # wave flux vs depth: sum over omega, ell of the source, damped cumulatively from the FeCZ top
    nu = np.logspace(-2, 1.5, 90)
    om = 2 * np.pi * nu / SECDAY
    ells = np.arange(1, int(min(1500, 1.2 * src.ell_max(om.max()))) + 1)   # effectively all excited ell
    dL = src.dL_dlnomega_dell(om, ells)                                  # (nu, ell)
    tc = D.tau_cumulative(path, om[:, None], ells[None, :].astype(float))   # (nu, ell, cell)
    dln = np.gradient(np.log(nu))
    Lw_all = np.einsum('ij,ijk,i->k', dL, np.exp(-tc), dln) / L
    low = ells <= 3
    Lw_low = np.einsum('ij,ijk,i->k', dL[:, low], np.exp(-tc[:, low]), dln) / L
    Lw0 = src.L_w / L
    xpath = logT[path['k']]

    # zones (above T = 1e6 K)
    zones = [(t, b, zone_type(p.T[b])) for t, b in p.conv_regions()]

    def shade(ax, labels=False):
        for t, b, name in zones:
            ax.axvspan(logT[b], logT[t], color=ZONE.get(name, ZONE['UNKNOWN']), lw=0, zorder=0)
            if labels:
                ax.text(0.5 * (logT[b] + logT[t]), 1.02, name, transform=ax.get_xaxis_transform(),
                        ha='center', va='bottom', fontsize=8, color=INK2)
        ax.axvline(logT[kt], color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.axvline(np.log10(p.T[p.k_phot]), color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=1)

    # ---------------- Figure 1: envelope ----------------
    fig, axs = plt.subplots(8, 1, figsize=(7.2, 15.5), sharex=True,
                            gridspec_kw=dict(height_ratios=[1, 0.9, 1, 0.8, 1.3, 1.6, 1.2, 1.0], hspace=0.12))
    fig.subplots_adjust(top=0.965, bottom=0.04)
    fig.suptitle(title, fontsize=9.5, color=INK, y=0.99)
    xm = x_lo + 0.62 * (x_hi - x_lo)      # x position for direct labels (radiative layer above the FeCZ)
    xm = min(xm, logT[kt] - 0.06)
    for i, ax in enumerate(axs):
        shade(ax, labels=(i == 0))
        ax.text(0.005, 0.95, '(' + 'abcdefgh'[i] + ')', transform=ax.transAxes, va='top', fontsize=8.5, color=INK)

    ax = axs[0]
    ax.plot(logT[env], p.kappa[env], color=C[0])
    ax.set_ylabel('κ [cm² g⁻¹]')
    ax.text(logT[kt], 0.05, ' FeCZ top: waves launched', transform=ax.get_xaxis_transform(),
            fontsize=7.5, color=INK2, ha='left', va='bottom')
    ax.text(np.log10(p.T[p.k_phot]), 0.05, 'τ = 2/3 ', transform=ax.get_xaxis_transform(),
            fontsize=7.5, color=MUTED, ha='right', va='bottom')

    ax = axs[1]
    ax.plot(logT[env], gradr[env], color=C[0], ls=LS[0])
    ax.plot(logT[env], grada[env], color=C[1], ls=LS[1])
    ax.set_ylim(0, max(0.5, np.nanmax(np.where(env, gradr, 0)) * 1.15))
    ax.set_ylabel('∇')
    ax.legend(ax.lines[-2:], ['∇$_{\\rm rad}$', '∇$_{\\rm ad}$'], loc='upper right', ncol=2)

    ax = axs[2]
    ax.semilogy(logT[env], p.cs[env] / 1e5, color=C[1], ls=LS[1])
    ax.semilogy(logT[env], vconv[env], color=C[0], ls=LS[0])
    ax.set_ylabel('velocity [km s⁻¹]')
    ax.set_ylim(0.1, 300)
    ax.legend(ax.lines[-2:], ['c$_s$', 'v$_{\\rm conv}$ (MLT)'], loc='upper right', ncol=2)
    label_at(ax, logT[env], (p.cs / 1e5)[env], xm, 'c$_s$')
    ax.text(0.99, 0.06, f"FeCZ: $\\mathcal{{M}}_t$ = {h['mach_FeCZ_aver_ahp']:.3f}, v$_{{\\max}}$ = {h['v_FeCZ_max'] / 1e5:.1f} km/s",
            transform=ax.transAxes, ha='right', fontsize=7.5, color=INK2)

    ax = axs[3]
    ax.semilogy(logT[env], p.Hp[env] / p.R, color=C[0])
    ax.set_ylabel('H$_P$ / R')
    ax.text(0.99, 0.08, f"FeCZ: {kb - kt + 1} cells, top {np.log(p.P[kt] / p.P[p.k_phot]):.2f} ln P below τ = 2/3",
            transform=ax.transAxes, ha='right', fontsize=7.5, color=INK2)

    ax = axs[4]
    lc = np.where(p.Lconv_frac > 1e-12, p.Lconv_frac, np.nan)
    ax.semilogy(logT[env], lc[env], color=C[0], ls=LS[0])
    ax.semilogy(xpath, Lw_all, color=C[1], ls=LS[1])
    ax.semilogy(xpath, Lw_low, color=C[2], ls=LS[2])
    ax.axhline(Lw0, color=INK2, lw=0.8, ls=(0, (4, 3)))
    ax.set_ylim(1e-12, 1.5)
    ax.set_ylabel('flux / F$_*$')
    ax.legend(ax.lines[2:5], ['L$_{\\rm conv}$/L (FeCZ)', 'wave flux, all ℓ', 'wave flux, ℓ ≤ 3 (visible)'],
              loc='upper right', ncol=1)
    label_at(ax, xpath, Lw_all, logT[kt] - 0.05, 'all ℓ', dy=5)
    label_at(ax, xpath, Lw_low, logT[kt] - 0.05, 'ℓ ≤ 3', dy=-5, va='top')
    ax.text(logT[kt], Lw0 * 1.8, f' launched: L$_w$/L = $\\mathcal{{M}}_t$F$_c$/F = {Lw0:.1e}', fontsize=7.5, color=INK2)
    ax.text(0.99, 0.04, 'quasi-adiabatic damping (invalid where q > 1, panel h)', transform=ax.transAxes,
            ha='right', fontsize=7, color=MUTED)

    ax = axs[5]
    ax.semilogy(logT[env], N_nu[env], color=C[0], ls=LS[0])
    ax.semilogy(logT[env], S1_nu[env], color=C[1], ls=LS[1])
    ax.semilogy(xpath, D.nu_nonad(path, 1.0), color=C[2], ls=LS[2])
    ax.semilogy(xpath, D.nu_nonad(path, 5.0), color=C[3], ls=LS[3])
    ax.fill_between(logT[env], 1e-3, np.fmin(N_nu[env], S1_nu[env]), where=np.isfinite(N_nu[env]),
                    color=C[0], alpha=0.07, lw=0)
    ax.axhline(nu_c, color=INK2, lw=0.8, ls=(0, (4, 3)))
    x_s = x_lo + 0.25 * (x_hi - x_lo)    # near-surface side, clear of the FeCZ-top crossings
    ax.text(x_s, nu_c * 0.85, f'ν$_c$ = 1/(2πt$_c$) = {nu_c:.2f} d⁻¹ (FeCZ turnover)', fontsize=7.5,
            color=INK2, ha='center', va='top')
    ax.set_ylim(1e-2, 1e3)
    ax.set_ylabel('ν [d⁻¹]')
    ax.legend(ax.lines[2:6], ['N/2π', 'S$_1$/2π (Lamb, ℓ=1)', 'ν below which waves are non-adiabatic, ℓ=1',
               'same, ℓ=5'], loc='upper right', ncol=2, fontsize=7)
    label_at(ax, logT[env], N_nu[env], xm, 'N/2π')
    label_at(ax, xpath, D.nu_nonad(path, 1.0), x_s + 0.12, 'q = 1, ℓ=1', dy=-5, va='top')
    label_at(ax, xpath, D.nu_nonad(path, 5.0), x_s + 0.12, 'q = 1, ℓ=5', dy=5)
    ax.text(xm, 0.015, 'g-mode cavity (ℓ=1): ν < N and ν < S$_1$', ha='center', fontsize=7, color=MUTED)

    ax = axs[6]
    for j, nuj in enumerate((0.3, 1.0, 3.0)):
        t = D.tau_cumulative(path, 2 * np.pi * nuj / SECDAY, 1.0)
        ax.semilogy(xpath, np.where(t > 1e-6, t, np.nan), color=C[j], ls=LS[j], label=f'ν = {nuj:g} d⁻¹')
        label_at(ax, xpath, np.where(t > 1e-6, t, np.nan), xm, f'{nuj:g} d⁻¹')
    ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.88), ncol=1)
    ax.axhline(1.0, color=MUTED, lw=0.8)
    ax.set_ylim(1e-5, 1e4)
    ax.set_ylabel('τ$_{\\rm rad}$ (ℓ=1),\nfrom FeCZ top')

    ax = axs[7]
    for j, nuj in enumerate((0.3, 1.0, 3.0)):
        q = D.q_nonad(path, 2 * np.pi * nuj / SECDAY, 1.0)
        q = np.where(q > 0, q, np.nan)
        ax.semilogy(xpath, q, color=C[j], ls=LS[j], label=f'ν = {nuj:g} d⁻¹')
        label_at(ax, xpath, q, xm, f'{nuj:g} d⁻¹')
    ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.88), ncol=1)
    ax.axhline(1.0, color=MUTED, lw=0.8)
    ax.text(x_hi - 0.01, 0.5, 'above the line (q > 1):\nquasi-adiabatic damping fails', fontsize=7.5,
            color=INK2, va='top')
    ax.set_ylim(1e-5, 1e5)
    ax.set_ylabel('q = K k$_r^2$/ω (ℓ=1)')
    ax.set_xlabel('log T [K]   (surface →)')
    ax.set_xlim(x_hi, x_lo)

    os.makedirs(args.outdir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(args.outdir, f'wave_envelope_{tag}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Figure 2: propagation diagram, whole star ----------------
    x = p.r / p.R
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.set_title(title, fontsize=9, color=INK)
    for t, b, name in zones:
        ax.axvspan(x[b], x[t], color=ZONE.get(name, ZONE['UNKNOWN']), lw=0, zorder=0)
    core = np.nonzero(p.mlt_type == 1)[0]
    core = core[p.T[core] > 1e6]
    if len(core):
        ax.axvspan(0, x[core.min()], color='#ebebe7', lw=0, zorder=0)
        ax.text(0.5 * x[core.min()], 0.02, 'convective\ncore', ha='center', va='bottom', fontsize=7.5, color=INK2)
    ax.semilogy(x, N_nu, color=C[0], ls=LS[0])
    ax.semilogy(x, S1_nu, color=C[1], ls=LS[1])
    ax.semilogy(x, S2_nu, color=C[2], ls=LS[2])
    g1 = np.fmin(N_nu, S1_nu)
    ax.fill_between(x, 1e-4, g1, where=np.isfinite(N_nu), color=C[0], alpha=0.10, lw=0)
    ax.axhline(nu_c, color=INK2, lw=0.8, ls=(0, (4, 3)))
    ax.text(0.25, nu_c * 1.2, f'FeCZ turnover ν$_c$ = {nu_c:.2f} d⁻¹', fontsize=7.5, color=INK2)
    if h.get('turnover_core', 0) > 0:
        nu_core = 1 / (2 * np.pi * h['turnover_core']) * SECDAY
        ax.axhline(nu_core, color=MUTED, lw=0.8, ls=(0, (1, 2)))
        ax.text(0.25, nu_core * 1.2, f'core turnover 1/(2πt$_c$) = {nu_core:.2g} d⁻¹', fontsize=7.5, color=INK2)
    ax.text(0.6, 0.02, 'g-mode cavity (ℓ=1): ν < N and ν < S$_1$', ha='center', fontsize=7.5, color=INK2)
    ax.annotate('FeCZ', (x[kt], 0.93), xycoords=('data', 'axes fraction'), xytext=(-6, 0),
                textcoords='offset points', ha='right', fontsize=8, color=INK2)
    ax.legend(ax.lines[:3], ['N/2π (buoyancy)', 'S$_1$/2π (Lamb, ℓ=1)', 'S$_2$/2π (Lamb, ℓ=2)'],
              loc='upper center', ncol=3)
    label_at(ax, x, N_nu, 0.55, 'N/2π')
    label_at(ax, x, S1_nu, 0.72, 'S$_1$/2π', dy=-5, va='top')
    label_at(ax, x, S2_nu, 0.33, 'S$_2$/2π')
    ax.set_xlim(0, 1.0)
    ax.set_ylim(1e-3, 1e4)
    ax.set_xlabel('r / R')
    ax.set_ylabel('ν [d⁻¹]')
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(args.outdir, f'propagation_{tag}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {args.outdir}/wave_envelope_{tag}.png/.pdf and {args.outdir}/propagation_{tag}.png/.pdf')


if __name__ == '__main__':
    main()
