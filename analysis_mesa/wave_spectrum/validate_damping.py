"""Check the Python damping integral against the run_star_extras history columns (template_v2).

Usage (from analysis_mesa/):  python3 -m wave_spectrum.validate_damping [RUN_DIR]
Default RUN_DIR: the local 20 Msun test /home/mcantiello/rednoise_tests/T1a_M20_w00.
Result 2026-09-24: identical to all printed digits (tau at ell=1 and ell_eddy, nu_damp ell=1,5, FeCZ cells).
"""
import sys
from .mesa_io import Track
from . import damping as D

run = sys.argv[1] if len(sys.argv) > 1 else '/home/mcantiello/rednoise_tests/T1a_M20_w00'
tr = Track(run)
print('model  Xc   | FeCZ cells py mesa | tau(wc,l=1) py mesa | tau(wc,l_eddy) py mesa | nu_damp l=1 py mesa | l=5 py mesa')
for p, h in tr.pairs():
    kt, kb = p.fecz()
    path = D.path_arrays(p, kt)
    wc, le = h['FeCZ_omega_c'], h['FeCZ_ell_eddy']
    print(f"{p.model_number:5d} {h['center_h1']:.2f} | {kb - kt + 1:4d} {h['FeCZ_ncells']:4.0f} | "
          f"{D.tau(path, wc, 1.0)[0]:9.3g} {h['FeCZ_tau_rad_wc_l1']:9.3g} | "
          f"{D.tau(path, wc, le)[0]:9.3g} {h['FeCZ_tau_rad_wc']:9.3g} | "
          f"{D.nu_damp(path, 1.0):5.3f} {h['FeCZ_nu_damp_l1']:5.3f} | {D.nu_damp(path, 5.0):5.2f} {h['FeCZ_nu_damp_l5']:5.2f}")
