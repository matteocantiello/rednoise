#!/usr/bin/env python3
"""
Digitise the per-sector SLF parameters of Van Daele+2026 from their Fig. (figs/confindance_intervals_and_pvalues_sectors.png,
arXiv source): log nu_char and log(alpha0/C_w) against log Teff (left) and log L (right), one colour per TESS sector.
Markers are found by colour (matplotlib tab10 at alpha ~0.7 over white), error bars removed by a morphological opening.
Points in the Teff and L panels are paired by equal ordinate within a sector, which recovers (log Teff, log L, value) per
star-sector. Axes are calibrated from the detected tick marks (see the pixel values below).
Writes data_obs/vandaele_digitised.csv.
"""
import os
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vandaele2026/src/figs/confindance_intervals_and_pvalues_sectors.png'
TAB10 = {27: (31, 119, 180), 28: (255, 127, 14), 67: (44, 160, 44), 68: (214, 39, 40)}
PANELS = {  # (row0, row1, col0, col1) inner frame
    'TL': (56, 1555, 290, 1721), 'TR': (56, 1555, 2013, 3445), 'BL': (1779, 3279, 290, 1721), 'BR': (1779, 3279, 2013, 3445)}
xT = lambda c: 4.6 - (c - 311.5) * 0.8 / (1668.5 - 311.5)
xL = lambda c: 4.0 + (c - 2111.5) * 2.0 / (3309.5 - 2111.5)
yNu = lambda r: -(r - 84.5) / (1542.5 - 84.5)
yA = lambda r: 2.0 - (r - 2078.5) * 1.5 / (2989.5 - 2078.5)


def markers(im, panel, sector, alpha=0.7):
    r0, r1, c0, c1 = PANELS[panel]
    sub = im[r0:r1, c0:c1]
    ref = alpha * np.array(TAB10[sector]) + (1 - alpha) * 255
    d_ref = np.sqrt(((sub - ref) ** 2).sum(-1))
    d_full = np.sqrt(((sub - np.array(TAB10[sector])) ** 2).sum(-1))       # overlaps darken towards the full colour
    mask = (d_ref < 45) | (d_full < 45)
    # other sectors' colours must be further away
    for s, c in TAB10.items():
        if s != sector:
            mask &= np.sqrt(((sub - (alpha * np.array(c) + (1 - alpha) * 255)) ** 2).sum(-1)) > d_ref
    mask = ndimage.binary_opening(mask, structure=np.ones((5, 5)))
    lab, n = ndimage.label(mask)
    out = []
    if n == 0:
        return out
    areas = ndimage.sum(mask, lab, range(1, n + 1))
    cms = ndimage.center_of_mass(mask, lab, range(1, n + 1))
    for a, (rr, cc) in zip(areas, cms):
        if a < 40:
            continue
        if panel == 'TL' and 73 <= rr + r0 <= 500 and 735 <= cc + c0 <= 1640:      # legend box
            continue
        out.append((rr + r0, cc + c0, a))
    return out


def main():
    im = np.asarray(Image.open(FIG).convert('RGB')).astype(float)
    rows = []
    for s in TAB10:
        for top, (pT, pL, fy) in (('nu', ('TL', 'TR', yNu)), ('a', ('BL', 'BR', yA))):
            mT, mL = markers(im, pT, s), markers(im, pL, s)
            medA = np.median([m[2] for m in mT + mL]) if mT + mL else np.nan
            for side, mm in (('T', mT), ('L', mL)):
                for rr, cc, a in mm:
                    rows.append(dict(sector=s, qty=top, side=side, y=fy(rr), x=(xT if side == 'T' else xL)(cc),
                                     area=a, n_merged=max(1, int(round(a / medA)))))
    d = pd.DataFrame(rows)
    print(d.groupby(['qty', 'sector', 'side']).agg(n=('y', 'size'), merged=('n_merged', lambda v: (v > 1).sum())).unstack())
    # pair Teff- and L-panel points of the same sector and quantity with equal ordinate
    pairs = []
    for (q, s), g in d.groupby(['qty', 'sector']):
        T, L = g[g.side == 'T'], g[g.side == 'L']
        used = set()
        for _, t in T.sort_values('y').iterrows():
            dy = abs(L.y - t.y)
            dy = dy[~dy.index.isin(used)]
            if len(dy) and dy.min() < 0.006:
                j = dy.idxmin(); used.add(j)
                pairs.append(dict(qty=q, sector=s, lT=t.x, lL=L.loc[j].x, val=t.y, dy=dy.min()))
    p = pd.DataFrame(pairs)
    print(p.groupby(['qty', 'sector']).size())
    d.to_csv(f'{HERE}/data_obs/vandaele_digitised_raw.csv', index=False, float_format='%.4f')
    p.to_csv(f'{HERE}/data_obs/vandaele_digitised.csv', index=False, float_format='%.4f')


if __name__ == '__main__':
    main()
