#!/usr/bin/env python3
"""
Two-panel geographic distribution map of L4.11 strains.
  (a) Pie charts by country showing sub-lineage proportions
  (b) Inferred dispersal routes (mugration) with arrows
"""
import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ── Configuration ──────────────────────────────────────────────────────────
STRAINS_CSV = '/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv'
OUTPUT_PNG = '/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/map_l411_geographic.png'

# L4.11.1 strain IDs
L4111_STRAINS = set("""
ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448 ERR4553470
ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568 ERR4553613
ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721 ERR4553747
ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834 ERR4553841
ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942 ERR4553947
ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739 SRR21661641
SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865 SRR29016881
SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138 SRR29017144
SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234 SRR29017242
SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583 SRR29017611
SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674 SRR29017683
SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528 SRR29440700
SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600 SRR35281602
SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186 SRR6797722 SRR6797801
""".split())

# Country centroids (lon, lat) — slightly adjusted to avoid overlaps
COUNTRY_COORDS = {
    'Peru': (-76.0, -10.0),
    'Bangladesh': (90.3, 23.7),
    'United States of America': (-100.0, 40.0),
    'Argentina': (-64.0, -34.6),
    'Italy': (12.5, 42.5),
    'Sweden': (16.0, 62.0),
    'South Africa': (25.0, -29.0),
    'Uganda': (32.3, 1.4),
    'Australia': (134.0, -25.3),
    'United Kingdom of Great Britain and Northern Ireland': (-4.0, 55.0),
    'Paraguay': (-58.4, -23.4),
    'Germany': (10.4, 51.2),
}

# Short display names
COUNTRY_SHORT = {
    'United States of America': 'USA',
    'United Kingdom of Great Britain and Northern Ireland': 'UK',
    'South Africa': 'S. Africa',
}

# Colors
C_L4111 = '#1976D2'   # Blue for L4.11.1
C_L4112 = '#E53935'   # Red for L4.11.2

# ── Load data ──────────────────────────────────────────────────────────────
with open(STRAINS_CSV, 'r') as f:
    reader = csv.DictReader(f)
    strains = list(reader)

country_counts = {}
for s in strains:
    country = s['Country'].strip().strip('"')
    if not country:
        continue
    sid = s['Id'].strip().strip('"')
    sublineage = 'L4.11.1' if sid in L4111_STRAINS else 'L4.11.2'
    if country not in country_counts:
        country_counts[country] = {'L4.11.1': 0, 'L4.11.2': 0}
    country_counts[country][sublineage] += 1

print("Country distribution:")
for c, v in sorted(country_counts.items(), key=lambda x: -(x[1]['L4.11.1'] + x[1]['L4.11.2'])):
    total = v['L4.11.1'] + v['L4.11.2']
    print(f"  {c:50s}  L4.11.1={v['L4.11.1']:3d}  L4.11.2={v['L4.11.2']:3d}  Total={total:3d}")

# ── Helper: draw pie chart on a cartopy axis ───────────────────────────────
def draw_pie(ax, lon, lat, sizes, colors, max_count, proj, label=None):
    """Draw a small pie chart at given lon/lat."""
    total = sum(sizes)
    if total == 0:
        return
    # Much smaller radius — scaled with sqrt
    radius = 1.8 * math.sqrt(total / max_count)
    radius = max(radius, 0.6)  # minimum visible size

    x, y = proj.transform_point(lon, lat, ccrs.PlateCarree())

    start = 0
    for size, color in zip(sizes, colors):
        if size == 0:
            continue
        theta1 = start * 360
        theta2 = (start + size / total) * 360
        wedge = mpatches.Wedge((x, y), radius * 1e6, theta1, theta2,
                               facecolor=color, edgecolor='white',
                               linewidth=0.4, transform=ax.transData,
                               zorder=10, alpha=0.85)
        ax.add_patch(wedge)
        start += size / total

    # Label with count
    if label is None:
        label = str(total)
    ax.text(x, y - radius * 1.4e6, label,
            transform=ax.transData, fontsize=4.5, ha='center', va='top',
            fontweight='bold', color='#333333', zorder=15,
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='none', alpha=0.7))

# ── Helper: draw arrow between two countries ───────────────────────────────
def draw_arrow(ax, lon1, lat1, lon2, lat2, proj, width=1.5, color='#333',
               style='->', label=None, bidirectional=False):
    """Draw a curved arrow between two lon/lat points."""
    x1, y1 = proj.transform_point(lon1, lat1, ccrs.PlateCarree())
    x2, y2 = proj.transform_point(lon2, lat2, ccrs.PlateCarree())

    arrow_style = '<->' if bidirectional else '->'

    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle=f'{arrow_style},head_width=5,head_length=4',
        connectionstyle='arc3,rad=0.15',
        linewidth=width, color=color, zorder=8,
        transform=ax.transData, alpha=0.7
    )
    ax.add_patch(arrow)

    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mx, my, label, transform=ax.transData,
                fontsize=5, ha='center', va='center', color=color,
                fontweight='bold', fontstyle='italic',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.8), zorder=12)

# ── Setup map features ─────────────────────────────────────────────────────
def setup_map(ax):
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='#f5f5f0', edgecolor='none')
    ax.add_feature(cfeature.OCEAN, facecolor='#e8f4f8')
    ax.add_feature(cfeature.BORDERS, linewidth=0.2, edgecolor='#bbbbbb')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#888888')

# ── Create two-panel figure ────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10), dpi=300)

proj = ccrs.Robinson()

# Panel A: Geographic distribution with pie charts
ax1 = fig.add_subplot(2, 1, 1, projection=proj)
setup_map(ax1)
ax1.set_title('(a) Geographic distribution by sub-lineage',
              fontsize=10, fontweight='bold', loc='left', pad=6)

max_count = max(v['L4.11.1'] + v['L4.11.2'] for v in country_counts.values())

for country, counts in country_counts.items():
    if country not in COUNTRY_COORDS:
        print(f"  ⚠ No coordinates for: {country}")
        continue
    lon, lat = COUNTRY_COORDS[country]
    total = counts['L4.11.1'] + counts['L4.11.2']
    short = COUNTRY_SHORT.get(country, country)
    lbl = f"{short} ({total})" if total >= 5 else str(total)
    draw_pie(ax1, lon, lat,
             [counts['L4.11.1'], counts['L4.11.2']],
             [C_L4111, C_L4112],
             max_count, proj, label=lbl)

# Legend
legend_patches = [
    mpatches.Patch(facecolor=C_L4111, edgecolor='white', label='L4.11.1 (n=91)'),
    mpatches.Patch(facecolor=C_L4112, edgecolor='white', label='L4.11.2 (n=737)'),
]
ax1.legend(handles=legend_patches, loc='lower left', fontsize=7,
           frameon=True, facecolor='white', edgecolor='#cccccc',
           framealpha=0.9, title='Sub-lineage', title_fontsize=8)

n_with_country = sum(v['L4.11.1'] + v['L4.11.2'] for v in country_counts.values())
ax1.text(0.99, 0.02,
         f'Strains with known country: {n_with_country}/828 ({n_with_country*100/828:.0f}%)',
         transform=ax1.transAxes, fontsize=6, ha='right', color='#888888')

# Panel B: Inferred dispersal routes
ax2 = fig.add_subplot(2, 1, 2, projection=proj)
setup_map(ax2)
ax2.set_title('(b) Inferred dispersal routes (mugration analysis)',
              fontsize=10, fontweight='bold', loc='left', pad=6)

# Draw dots for key countries
key_countries = {
    'Bangladesh': (C_L4111, 62),
    'Peru': (C_L4112, 410),
    'United States of America': ('#9C27B0', 58),
    'Argentina': (C_L4112, 32),
    'Uganda': (C_L4111, 4),
    'South Africa': (C_L4111, 5),
    'Italy': (C_L4112, 8),
    'Sweden': (C_L4112, 5),
}

for country, (color, count) in key_countries.items():
    if country not in COUNTRY_COORDS:
        continue
    lon, lat = COUNTRY_COORDS[country]
    x, y = proj.transform_point(lon, lat, ccrs.PlateCarree())
    size = 4 * math.sqrt(count)
    ax2.plot(x, y, 'o', markersize=size, color=color,
             markeredgecolor='white', markeredgewidth=0.5,
             transform=ax2.transData, zorder=10, alpha=0.8)
    short = COUNTRY_SHORT.get(country, country)
    ax2.text(x, y - size * 0.7e5, f'{short}\n(n={count})',
             transform=ax2.transData, fontsize=4.5, ha='center', va='top',
             fontweight='bold', color='#333',
             bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                       edgecolor='none', alpha=0.7), zorder=12)

# Draw dispersal arrows
bd = COUNTRY_COORDS['Bangladesh']
pe = COUNTRY_COORDS['Peru']
us = COUNTRY_COORDS['United States of America']
ar = COUNTRY_COORDS['Argentina']
ug = COUNTRY_COORDS['Uganda']

# Main ancestral dispersal: Bangladesh → Peru (via east, crossing Pacific/Indian Ocean)
draw_arrow(ax2, bd[0], bd[1], pe[0], pe[1], proj,
           width=2.5, color='#E65100', label='Ancestral\ndispersal')

# Peru ↔ USA (strongest: W=14.9)
draw_arrow(ax2, pe[0], pe[1], us[0], us[1], proj,
           width=2.0, color='#1565C0', label='W=14.9',
           bidirectional=True)

# Peru → Argentina (W=7.7)
draw_arrow(ax2, pe[0], pe[1], ar[0], ar[1], proj,
           width=1.5, color='#2E7D32', label='W=7.7')

# Bangladesh → Uganda (W=1.4)
draw_arrow(ax2, bd[0], bd[1], ug[0], ug[1], proj,
           width=1.0, color='#6A1B9A', label='W=1.4')

# Legend for panel B
from matplotlib.lines import Line2D
legend_items = [
    Line2D([0], [0], color='#E65100', linewidth=2.5, label='Ancestral dispersal (Bangladesh → Peru)'),
    Line2D([0], [0], color='#1565C0', linewidth=2.0, label='Peru ↔ USA (W=14.9)'),
    Line2D([0], [0], color='#2E7D32', linewidth=1.5, label='Peru → Argentina (W=7.7)'),
    Line2D([0], [0], color='#6A1B9A', linewidth=1.0, label='Bangladesh → Uganda (W=1.4)'),
]
ax2.legend(handles=legend_items, loc='lower left', fontsize=6,
           frameon=True, facecolor='white', edgecolor='#cccccc',
           framealpha=0.9, title='Migration corridors', title_fontsize=7)

ax2.text(0.99, 0.02,
         'W = symmetrised GTR transition rate from TreeTime mugration',
         transform=ax2.transAxes, fontsize=5.5, ha='right', color='#888888')

plt.tight_layout(h_pad=1.5)
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ Two-panel map saved to: {OUTPUT_PNG}")
