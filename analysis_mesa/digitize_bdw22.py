#!/usr/bin/env python3
"""
Bowman & Dorn-Wallenstein 2022 (A&A 668, A134; arXiv 2211.08347) give no per-star table. Their appendix summary figures
(figures/tess<TIC>_fitted_log_FT.png, one per star, from the arXiv source in handoff/data_survey_raw/bowman_dw2022) mark
nu_char from the Bowman+2020 amplitude-spectrum fit (lime dotted line) and from their celerite2 GP (SHO) fit (navy dotted line)
on a log-frequency axis. We read both line positions at the pixel level, calibrating the axis per image from its tick marks.

Validation: the lime value must reproduce the tabulated Bowman 2020 nu_char (J/A+A/640/A36 table A2).
Output: data_obs/bdw22_nuchar.csv; prints the GP <-> Lorentzian conversion.
"""
import os
import glob
import re
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/bowman_dw2022/figures'
B20 = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data/bowman2020_A36_tablea2.tsv'


def clusters(idx, gap=2):
    """Group sorted integer indices into runs; return (centre, size) per run."""
    if len(idx) == 0:
        return []
    out, start = [], [idx[0]]
    for i in idx[1:]:
        if i - start[-1] <= gap:
            start.append(i)
        else:
            out.append((np.mean(start), len(start))); start = [i]
    out.append((np.mean(start), len(start)))
    return out


def read_panel(path):
    im = np.asarray(Image.open(path).convert('RGB')).astype(int)
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    H = im.shape[0]
    dark = im.max(2) < 80
    # bottom panel frame: the full-width dark rows below the top panel are its top and bottom spines
    rows = dark.sum(1); lower = np.arange(int(0.38 * H), H)
    cand = lower[rows[lower] > 0.6 * im.shape[1]]
    top, bot = cand.min(), cand.max()
    # outward x ticks just below the bottom spine
    tk = np.where(dark[bot + 1:bot + 12].sum(0) >= 4)[0]
    cc = [c for c, n in clusters(tk)]
    left = min(c for c, n in clusters(np.where(dark[top:bot].sum(0) > 0.8 * (bot - top))[0]))
    cc = [c for c in cc if c > left + 5]
    gaps = np.diff(cc)
    # a decade tick is preceded by a small gap (9 -> 10) and followed by a large one (10 -> 20)
    dec = [cc[i] for i in range(1, len(cc) - 1) if gaps[i - 1] < 0.25 * gaps[i]]
    D = np.median(np.diff(dec))
    x1 = dec[0]
    if abs((x1 - left) - D) > 3:          # the first decade must sit one decade right of the 0.1 d-1 spine
        raise RuntimeError(f'axis calibration failed for {path}: left {left}, first decade {x1}, D {D}')
    band = slice(top + 2, bot - 2)
    lime = (g[band] > 200) & (r[band] < 120) & (b[band] < 120)
    navy = (b[band] > 90) & (r[band] < 40) & (g[band] < 40)
    out = {}
    for name, mask in (('B20', lime), ('GP', navy)):
        cs = mask.sum(0)
        runs = [(c, n) for c, n in clusters(np.where(cs > 0.15 * (bot - top))[0]) if c > left + 3]
        if len(runs) != 1:
            out[name] = np.nan; continue
        out[name] = 10 ** ((runs[0][0] - x1) / D)
    return out, D


def main():
    rows = []
    for p in sorted(glob.glob(f'{SRC}/tess0*_fitted_log_FT.png')):
        tic = int(re.search(r'tess0*(\d+)_fitted', p).group(1))
        v, D = read_panel(p)
        rows.append(dict(TIC=tic, nuchar_B20_fig=v['B20'], nuchar_GP=v['GP'], px_per_decade=D))
    d = pd.DataFrame(rows)
    b = pd.read_csv(B20, sep='\t', comment='#')
    b = b[pd.to_numeric(b.TIC, errors='coerce').notna()]
    b['TIC'] = b.TIC.astype(int); b['nuchar'] = b.nuchar.astype(float); b['Name'] = b.Name.str.strip()
    d = d.merge(b[['TIC', 'Name', 'nuchar']].rename(columns={'nuchar': 'nuchar_B20_tab'}), on='TIC', how='left')
    ok = np.isfinite(d.nuchar_B20_fig) & np.isfinite(d.nuchar_B20_tab)
    dv = np.log10(d.nuchar_B20_fig[ok] / d.nuchar_B20_tab[ok])
    print(f'{len(d)} panels; GP read in {np.isfinite(d.nuchar_GP).sum()}, B20 line in {np.isfinite(d.nuchar_B20_fig).sum()}, '
          f'matched to B20 table {d.nuchar_B20_tab.notna().sum()}')
    print(f'validation, log(B20 figure / B20 table): median {np.median(dv):+.3f}, max |.| {np.max(abs(dv)):.3f} dex (N = {ok.sum()})')
    g = np.isfinite(d.nuchar_GP) & d.nuchar_B20_tab.notna()
    lg, lb = np.log10(d.nuchar_GP[g]), np.log10(d.nuchar_B20_tab[g])
    sl, ic = np.polyfit(lb, lg, 1)
    print(f'GP vs Lorentzian (B20 table), N = {g.sum()}: median log(GP/B20) {np.median(lg - lb):+.3f}, '
          f'scatter {np.std(lg - lb):.3f} dex; fit log GP = {ic:+.3f} + {sl:.3f} log B20; Spearman {spearmanr(lg, lb).correlation:.2f}')
    d.to_csv(f'{HERE}/data_obs/bdw22_nuchar.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 200)
    print(d.sort_values('nuchar_B20_tab').to_string(index=False, float_format=lambda v: f'{v:.3f}'))


if __name__ == '__main__':
    main()
