"""Run GYRE (non-adiabatic, l = 1-3) on a MESA pulse-data file.

Usage (any directory):
    python3 run_gyre.py MODEL.GYRE OUTDIR [--outer VACUUM|UNNO|JCD|ISOTHERMAL] [--branch E_NEG|F_NEG]
                        [--fmin 0.05] [--fmax 5] [--pts-per-mode 3] [--threads 8]

The number of frequency points per degree is set from the asymptotic g-mode period spacing
(Pi_0 / sqrt(l(l+1)), with Pi_0 computed from the model's N^2), so that each mode gets about
--pts-per-mode points of the INVERSE scan. Output: OUTDIR/{gyre.in, gyre.log, summary_*.h5, detail.*.h5}.
Keep OUTDIR on local disk (detail files are numerous).
"""
import argparse
import os
import subprocess
import numpy as np

GYRE_DIR = '/mnt/home/mcantiello/software/gyre-9.1.1'
SDK = '/mnt/home/mcantiello/mesasdk-26.6.1'
HERE = os.path.dirname(os.path.abspath(__file__))


def pi0_from_gyre_file(path):
    """Asymptotic Pi_0 = 2 pi^2 / int N dr/r [s] from a MESA GYRE-format file (cols: k r m L P T rho nabla N2 ...)."""
    d = np.loadtxt(path, skiprows=1)
    r, N2 = d[:, 1], d[:, 8]
    s = r > 0
    return 2 * np.pi ** 2 / np.trapezoid(np.sqrt(np.clip(N2[s], 0, None)) / r[s], r[s])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('outdir')
    ap.add_argument('--outer', default='VACUUM')
    ap.add_argument('--branch', default='E_NEG')
    ap.add_argument('--fmin', type=float, default=0.05)
    ap.add_argument('--fmax', type=float, default=5.0)
    ap.add_argument('--pts-per-mode', type=float, default=3.0)
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()

    model = os.path.abspath(a.model)
    pi0 = pi0_from_gyre_file(model)
    span = 86400.0 * (1 / a.fmin - 1 / a.fmax)          # period range [s]
    nf = {l: int(np.ceil(a.pts_per_mode * span * np.sqrt(l * (l + 1)) / pi0)) + 10 for l in (1, 2, 3)}
    print(f'Pi_0 = {pi0:.0f} s; expected modes per l: '
          + ', '.join(f'l={l}: {span * np.sqrt(l * (l + 1)) / pi0:.0f}' for l in (1, 2, 3))
          + f'; n_freq = {nf}')

    tpl = open(os.path.join(HERE, 'nad_lowl.in.template')).read()
    for k, v in {'@MODEL@': model, '@OUTER_BOUND@': a.outer, '@BRANCH@': a.branch, '@FMIN@': str(a.fmin),
                 '@FMAX@': str(a.fmax), '@NF1@': str(nf[1]), '@NF2@': str(nf[2]), '@NF3@': str(nf[3])}.items():
        tpl = tpl.replace(k, v)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, 'gyre.in'), 'w') as f:
        f.write(tpl)

    cmd = (f'export MESASDK_ROOT={SDK}; source $MESASDK_ROOT/bin/mesasdk_init.sh; '
           f'export GYRE_DIR={GYRE_DIR}; export OMP_NUM_THREADS={a.threads}; '
           f'cd {os.path.abspath(a.outdir)} && $GYRE_DIR/bin/gyre gyre.in > gyre.log 2>&1')
    r = subprocess.run(['bash', '-c', cmd])
    print('GYRE exit code', r.returncode, '->', os.path.join(a.outdir, 'gyre.log'))


if __name__ == '__main__':
    main()
