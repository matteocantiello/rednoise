"""Shared helpers for the rednoise MESA grid: model status and cached history loading.

Status is read from each model's rn.out (last "termination code" line), history from
LOGS/history.data. Parsed histories are cached as .npz on ceph, keyed by file size
and mtime, so re-plotting only re-reads models that changed.
"""
import os, re, glob, time
import numpy as np

# Grid to read: this directory (v1) by default; set RN_GRID to another grid tree (e.g. ../grids_v2).
# Each non-default grid gets its own cache directory, because cache files are keyed by sub-grid/mass only.
_HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.abspath(os.environ.get('RN_GRID', _HERE))
CACHE = '/mnt/ceph/users/mcantiello/rednoise/cache/history_npz'
if GRID != _HERE:
    CACHE += '_' + os.path.basename(GRID)

MASSES = [5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0,
          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
          21, 22, 23, 24, 25, 30, 40, 50, 60, 80, 100, 120]
SUBGRIDS = [(z, w) for z in ('MW', 'LMC', 'SMC') for w in ('w0.0', 'w0.2', 'w0.4', 'w0.6')]
SUBGRIDS += [('MW_mltpp', 'w0.0')]
SUBGRIDS = [(z, w) for z, w in SUBGRIDS if os.path.isdir(f'{GRID}/{z}/{w}')]   # only those set up
ZVAL = {'MW': 0.014, 'MW_mltpp': 0.014, 'LMC': 0.006, 'SMC': 0.002}
ZTITLE = {'MW': r'MW ($Z=0.014$)', 'LMC': r'LMC ($Z=0.006$)',
          'SMC': r'SMC ($Z=0.002$)', 'MW_mltpp': r'MW + MLT++ ($Z=0.014$)'}

# Spectroscopic luminosity, user's normalisation
ELL_SUN = 5777.0**4 / (274 * 100)


def mdir(mass):
    return f'M{int(mass)}' if mass >= 10 else f'M{mass:.1f}'


def model_dir(z, w, mass):
    return f'{GRID}/{z}/{w}/{mdir(mass)}'


def spec_ell(log_teff, log_g):
    return np.log10((10**np.asarray(log_teff))**4 / 10**np.asarray(log_g) / ELL_SUN)


def model_status(z, w, mass, running_hours=1.0):
    """Return dict(status, min_dt) from rn.out.

    status: 'done' (xa_central_lower_limit), 'running' (no termination code, rn.out
    updated within running_hours), 'crash' (Fortran ERROR STOP), 'stalled', 'nostart',
    or the raw termination code.
    min_dt: the min_timestep_limit MESA echoed on a min_timestep failure (else None).
    """
    out = f'{model_dir(z, w, mass)}/rn.out'
    if not os.path.exists(out):
        return {'status': 'nostart', 'min_dt': None}
    code, min_dt = None, None
    # Termination lines are at the end; read only the tail of large files
    with open(out, 'rb') as f:
        f.seek(0, 2)
        size = f.tell()
        f.seek(max(0, size - 200_000))
        tail = f.read().decode('ascii', 'replace')
    m = re.findall(r'termination code:\s*(\S+)', tail)
    if m:
        code = m[-1]
    d = re.findall(r'min_timestep_limit\s+([0-9.]+D[-+]\d+)', tail)
    if d:
        min_dt = float(d[-1].replace('D', 'E'))
    if code is None and re.search(r'ERROR STOP|Error termination', tail):
        return {'status': 'crash', 'min_dt': None}   # Fortran abort, no termination code
    if code is None:
        age_h = (time.time() - os.path.getmtime(out)) / 3600
        return {'status': 'running' if age_h < running_hours else 'stalled', 'min_dt': None}
    return {'status': 'done' if code == 'xa_central_lower_limit' else code, 'min_dt': min_dt}


BASIC_COLS = ['model_number', 'star_age', 'star_mass', 'log_abs_mdot', 'log_dt',
              'log_L', 'log_Teff', 'log_R', 'log_g', 'center_h1', 'center_he4',
              'surf_avg_v_rot', 'surf_avg_omega_div_omega_crit', 'surf_avg_v_div_v_crit',
              'surf_avg_Lrad_div_Ledd', 'max_L_rad_div_Ledd', 'surface_he4', 'surface_n14',
              'he_core_mass', 'subsurface_convective_regions', 'num_retries',
              'rho_surf', 'photosphere_csound', 'log_L_div_Ledd', 'F0_div_omega_c',
              'r_hp_1', 'r_hp_2', 'r_hp_3', 'r_hp_4', 'r_hp_5', 'r_hp_6', 'r_hp_7', 'r_hp_8',
              'r_hp_10', 'r_hp_15', 'r_hp_20', 'r_hp_30', 'r_hp_50', 'r_hp_100']
# every sub-surface convection-zone column (FeCZ, HeII, HeI, HI) plus the core-convection ones
CZ_PATTERN = re.compile(r'(FeCZ|HeII|HeI|HI)|_core$|core_')


def _history_columns(path):
    with open(path) as f:
        for _ in range(5):
            f.readline()
        allcols = f.readline().split()
        data_start = f.tell()
    return allcols, data_start


def _parse_rows(path, allcols, cols, start, end):
    """Parse complete rows in byte range [start, end) of history.data with awk."""
    import subprocess
    idx = [allcols.index(c) + 1 for c in cols]
    prog = f"NF=={len(allcols)} {{print " + ",".join(f"${i}" for i in idx) + "}"
    with open(path, 'rb') as f:
        f.seek(start)
        raw = f.read(end - start)
    txt = subprocess.run(['awk', prog], input=raw, capture_output=True, check=True).stdout
    arr = np.array(txt.decode().replace('D', 'E').split(), dtype=float)
    return arr.reshape(-1, len(cols))


def _dedupe(arr, cols):
    """After a restart from a photo MESA re-writes models already in the file:
    keep the last occurrence (walk backwards keeping strictly decreasing model_number)."""
    mn = arr[:, cols.index('model_number')]
    keep = np.zeros(len(mn), bool)
    cur = np.inf
    for i in range(len(mn) - 1, -1, -1):
        if mn[i] < cur:
            keep[i] = True
            cur = mn[i]
    return arr[keep]


def load_history(z, w, mass):
    """Return dict column -> array for one model (BASIC_COLS + convection-zone columns),
    or None if no usable history.

    Cache is incremental: the raw parsed rows and the byte offset already parsed are stored,
    so a growing history (running model) only has its new tail read from ceph. If the file
    shrank, its first row changed (model restarted from scratch) or the header changed,
    it is re-parsed in full.
    """
    path = f'{model_dir(z, w, mass)}/LOGS/history.data'
    if not os.path.exists(path):
        return None
    size = os.path.getsize(path)
    cfile = f'{CACHE}/{z}_{w}_{mdir(mass)}.npz'
    try:
        allcols, data_start = _history_columns(path)
        cols = [c for c in allcols if c in BASIC_COLS or CZ_PATTERN.search(c)]
        with open(path, 'rb') as f:  # first data row identifies the run (./rn rewrites the file)
            f.seek(data_start)
            first = f.readline()[:200].decode('ascii', 'replace')
        raw, offset = np.empty((0, len(cols))), data_start
        if os.path.exists(cfile):
            c = np.load(cfile, allow_pickle=False)
            if list(c['cols']) == cols and int(c['offset']) <= size and str(c['first']) == first:
                raw, offset = c['raw'], int(c['offset'])
        if offset < size:
            # only parse up to the last complete line
            with open(path, 'rb') as f:
                f.seek(max(offset, size - 1_000_000))
                last_nl = f.tell() + f.read().rfind(b'\n') + 1
            if last_nl > offset:
                new = _parse_rows(path, allcols, cols, offset, last_nl)
                raw = np.vstack([raw, new.astype(np.float64)])
                offset = last_nl
                os.makedirs(CACHE, exist_ok=True)
                tmp = cfile + '.tmp.npz'
                np.savez(tmp, raw=raw, offset=offset, cols=np.array(cols), first=first)
                os.replace(tmp, cfile)
    except Exception as e:
        print(f'  could not parse {path}: {e}')
        return None
    arr = _dedupe(raw, cols)
    if len(arr) < 5:
        return None
    return {c: arr[:, i] for i, c in enumerate(cols)}


def find_zams(h1):
    """ZAMS index: where center_h1 first drops 0.003 below its maximum."""
    mask = h1 < (h1.max() - 0.003)
    return int(np.argmax(mask)) if mask.any() else 0


def find_tams(h1, xlim=1e-4):
    """TAMS index: first point with center_h1 < xlim (None if not reached)."""
    idx = np.where(h1 < xlim)[0]
    return int(idx[0]) if len(idx) else None


def load_subgrid(z, w, masses=MASSES, verbose=True):
    """Load all histories of a sub-grid, in parallel threads (I/O bound)."""
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(12) as ex:
        res = list(ex.map(lambda m: load_history(z, w, m), masses))
    tracks = {m: d for m, d in zip(masses, res) if d is not None}
    if verbose:
        print(f'  {z}/{w}: {len(tracks)} tracks loaded')
    return tracks


def last_row(z, w, mass):
    """Last complete row of history.data as a dict (cheap: reads only the file tail)."""
    path = f'{model_dir(z, w, mass)}/LOGS/history.data'
    if not os.path.exists(path):
        return None
    allcols, data_start = _history_columns(path)
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        f.seek(max(data_start, size - 200_000))
        lines = f.read().decode('ascii', 'replace').split('\n')
    for l in reversed(lines):
        p = l.split()
        if len(p) == len(allcols):
            try:
                return dict(zip(allcols, (float(x.replace('D', 'E')) for x in p)))
            except ValueError:
                continue
    return None
