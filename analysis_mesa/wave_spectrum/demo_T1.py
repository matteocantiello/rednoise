"""Forward-model spectra on the local test track, with the red-noise fit, for a set of assumptions.

Usage (from analysis_mesa/):  python3 -m wave_spectrum.demo_T1 [--q Q ...] [--p P ...] [--run RUN_DIR] [--plot]
  --q: non-adiabatic transition threshold(s) (see damping.tau); 'none' = quasi-adiabatic to the photosphere
  --p: low-ell slope(s) of the source (see source.Source)
Also prints where the damping accrues and the local non-adiabaticity q (the diagnostic that showed the
quasi-adiabatic integral is invalid in the outer pressure scale height).
"""
import argparse
import numpy as np
from .mesa_io import Track
from . import damping as D, spectrum as S

ap = argparse.ArgumentParser()
ap.add_argument('--run', default='/home/mcantiello/rednoise_tests/T1a_M20_w00')
ap.add_argument('--q', nargs='+', default=['none', '1', '10', '100', '1000'])
ap.add_argument('--p', nargs='+', type=float, default=[1.0])
ap.add_argument('--a', type=float, default=6.5)
ap.add_argument('--plot', action='store_true', help='write figures/wave_spectrum_T1.png')
args = ap.parse_args()

tr = Track(args.run)
pairs = tr.pairs()

# where does the damping accrue? (one mid-MS profile, ell = 1)
p, h = min(pairs, key=lambda x: abs(x[1]['center_h1'] - 0.5))
path = D.path_arrays(p, p.fecz()[0])
lnP = np.log(p.P[path['k']] / p.P[p.k_phot])
print(f'model {p.model_number} (X_c = {h["center_h1"]:.2f}): FeCZ top {np.log(p.P[p.fecz()[0]] / p.P[p.k_phot]):.2f} '
      f'ln P above the photosphere')
for nu in (0.3, 1.0, 3.0):
    w = 2 * np.pi * nu / 86400
    cum = D.tau_profile(path, w, 1.0)
    dtau = np.r_[cum[:-1] - cum[1:], cum[-1]]
    q = D.q_nonad(path, w, 1.0)
    outer = lnP < 1
    print(f'  nu = {nu:3.1f}/d: tau = {cum[0]:9.3g}; share from the outer H_P {dtau[outer].sum() / cum[0]:5.1%} '
          f'(median q there {np.median(q[outer]):.3g}); share from cells with q > 1: {dtau[q > 1].sum() / cum[0]:5.1%}')

print('\n q      p   X_c  | nu_c  nu_damp1 | red-noise fit: nu_char  gamma  rms[dex] | f_trans  f_vis')
curves = {}
for qs in args.q:
    Q = None if qs == 'none' else float(qs)
    for pp in args.p:
        for p, h in pairs[::2]:
            m = S.model_spectrum(p, h, stop_at_q=Q, p=pp, a=args.a)
            amp = np.sqrt(m['psd'])
            try:
                f = S.fit_red_noise(m['nu'], amp)
                fit = f"{f['nu_char']:7.2f} {f['gamma']:6.2f} {f['rms_dex']:6.2f}"
            except Exception:
                fit = '   fit failed      '
            print(f"{qs:>5s} {pp:4.1f} {h['center_h1']:.2f} | {m['nu_c']:4.2f} {m['nu_damp_l1']:6.2f} | {fit:>36s} | "
                  f"{m['f_trans']:.1e} {m['f_vis']:.1e}")
            curves.setdefault((qs, pp), []).append((h['center_h1'], m['nu'], amp))

if args.plot:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, len(args.q), figsize=(3.2 * len(args.q), 3.2), sharey=True)
    for i, qs in enumerate(args.q):
        for xc, nu, amp in curves[(qs, args.p[0])]:
            ax[i].loglog(nu, amp / amp.max(), label=f'X_c={xc:.2f}')
        ax[i].set_title(f'q = {qs}')
        ax[i].set_xlabel('nu [1/d]')
        ax[i].set_ylim(1e-4, 2)
    ax[0].set_ylabel('amplitude (normalised)')
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig('figures/wave_spectrum_T1.png', dpi=130)
    print('wrote figures/wave_spectrum_T1.png')
