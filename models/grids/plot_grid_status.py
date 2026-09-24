#!/usr/bin/env python3
"""Status of every model in the rednoise MESA grid, read live from rn.out + history.data.

One row per sub-grid (Z, omega/omega_crit), one column per initial mass.
Also writes grid_status.txt (one line per model) next to the figure.
"""
import os
from concurrent.futures import ThreadPoolExecutor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import grid_io as g

XC_TAMS = 1e-3   # centre H below this = main sequence finished

STATUS = {  # key: (face color, glyph, legend label)
    'done':     ('#0ca30c', '✓', 'complete (He exhaustion)'),
    'running':  ('#2a78d6', '▶', 'running'),
    'crawl':    ('#6a3fb5', '…', 'running but stuck before ZAMS'),
    'fail_pms': ('#ec835a', 'T', 'failed after TAMS (MS usable)'),
    'fail_ms':  ('#d03b3b', '✗', 'failed on MS (value = $X_c$)'),
    'other':    ('#8a8a86', '?', 'stalled / not started'),
}


def classify(z, w, m):
    st = g.model_status(z, w, m)
    row = g.last_row(z, w, m)
    xc = row['center_h1'] if row else None
    s = st['status']
    if s == 'running' and row and row['star_age'] < 1e4 and row['model_number'] > 2e4:
        cls = 'crawl'   # tens of thousands of steps and still < 10^4 yr old: dt collapsed on the pre-MS
    elif s in ('done', 'running'):
        cls = s
    elif s in ('nostart', 'stalled'):
        cls = 'other'
    else:
        cls = 'fail_pms' if (xc is not None and xc < XC_TAMS) else 'fail_ms'
    return dict(z=z, w=w, m=m, cls=cls, code=s, min_dt=st['min_dt'], xc=xc,
                model=int(row['model_number']) if row else None,
                age=row['star_age'] if row else None)


jobs = [(z, w, m) for z, w in g.SUBGRIDS for m in g.MASSES]
with ThreadPoolExecutor(16) as ex:
    res = list(ex.map(lambda a: classify(*a), jobs))

# ---- text summary
with open(f'{g.GRID}/grid_status.txt', 'w') as f:
    f.write('# subgrid mass class termination_code min_timestep_limit model_number Xc star_age\n')
    for r in res:
        f.write(f"{r['z']}/{r['w']} {g.mdir(r['m'])} {r['cls']} {r['code']} {r['min_dt']} "
                f"{r['model']} {r['xc']} {r['age']}\n")

# ---- figure
rows = [f'{z}/{w}' for z, w in g.SUBGRIDS]
fig, ax = plt.subplots(figsize=(15, 6.8))
for r in res:
    i = rows.index(f"{r['z']}/{r['w']}")
    j = g.MASSES.index(r['m'])
    color, glyph, _ = STATUS[r['cls']]
    ax.add_patch(Rectangle((j + 0.04, i + 0.04), 0.92, 0.92, fc=color, ec='white', lw=0))
    txt = glyph
    if r['cls'] == 'fail_ms' and r['xc'] is not None:
        txt = f"{r['xc']:.3f}"[1:] if r['xc'] < 0.1 else f"{r['xc']:.2f}"[1:]
    ax.text(j + 0.5, i + 0.5, txt, ha='center', va='center', color='white',
            fontsize=6.5 if r['cls'] == 'fail_ms' else 9, fontweight='bold')
    # failures from before the 2026-09-23 fix (min_timestep_limit = 0.1 s): rerun candidates
    if r['min_dt'] is not None and r['min_dt'] > 1e-3:
        ax.plot(j + 0.86, i + 0.16, 'o', ms=3.2, mfc='black', mec='white', mew=0.5)

# row summaries on the right
for i, row in enumerate(rows):
    sub = [r for r in res if f"{r['z']}/{r['w']}" == row]
    n = {k: sum(r['cls'] == k for r in sub) for k in STATUS}
    ax.text(len(g.MASSES) + 0.3, i + 0.5,
            f"{n['done']:2d} done  {n['running'] + n['crawl']:2d} run  {n['fail_pms']:2d} post-TAMS  {n['fail_ms']:2d} MS",
            va='center', fontsize=8, family='monospace', color='0.2')

ax.set_xlim(0, len(g.MASSES))
ax.set_ylim(len(rows), 0)
ax.set_xticks([j + 0.5 for j in range(len(g.MASSES))])
ax.set_xticklabels([f'{m:g}' for m in g.MASSES], fontsize=8)
ax.set_yticks([i + 0.5 for i in range(len(rows))])
ax.set_yticklabels([r.replace('/w', r'  $\omega/\omega_c$=') for r in rows], fontsize=9)
for y in (4, 8, 12):
    ax.axhline(y, color='k', lw=1.2)
ax.set_xlabel(r'Initial mass ($M_\odot$)')
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

tot = {k: sum(r['cls'] == k for r in res) for k in STATUS}
n_old = sum(1 for r in res if r['cls'].startswith('fail') and r['min_dt'] and r['min_dt'] > 1e-3)
ax.set_title(f"MESA grid status — {len(res)} models: {tot['done']} complete, {tot['running'] + tot['crawl']} running, "
             f"{tot['fail_pms']} failed after TAMS, {tot['fail_ms']} failed on MS",
             fontsize=12, loc='left')

handles = [Patch(fc=STATUS[k][0], label=f'{STATUS[k][1]}  {STATUS[k][2]}') for k in STATUS if tot[k]]
handles.append(Line2D([], [], ls='', marker='o', ms=5, mfc='black', mec='white',
                      label=f'failed with old min_timestep_limit=0.1 s → rerun candidate ({n_old})'))
ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=3,
          fontsize=8.5, frameon=False)

plt.tight_layout()
outpath = f'{g.GRID}/grid_status.png'
fig.savefig(outpath, dpi=150, bbox_inches='tight', facecolor='white')
print(f'Saved to {outpath}')
print({k: v for k, v in tot.items()}, 'old-limit failures:', n_old)
