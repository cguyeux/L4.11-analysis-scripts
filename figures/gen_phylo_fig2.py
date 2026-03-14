#!/usr/bin/env python3
"""
Generate the rooted phylogenetic tree figure for L4.11 article.
- Prunes 3 basal strains (not truly basal phylogenetically)
- Roots the tree using midpoint rooting
- Colors tips: blue=L4.11.1, red=L4.11.2
- Circular layout
"""
import os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

# === Paths ===
TREE_FILE = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/T3_l411_only.nwk"
OUTPUT_PNG = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/phylo_tree_l411.png"
OUTPUT_PDF = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/phylo_tree_l411.pdf"

# === Strains to prune (not truly basal) ===
BASAL_TO_PRUNE = {'ERR3802166', 'ERR4195920', 'ERR4195919'}

# === L4.11.1 strains (91) ===
L4111 = set("""ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448
ERR4553470 ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568
ERR4553613 ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721
ERR4553747 ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834
ERR4553841 ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942
ERR4553947 ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739
SRR21661641 SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865
SRR29016881 SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138
SRR29017144 SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234
SRR29017242 SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583
SRR29017611 SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674
SRR29017683 SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528
SRR29440700 SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600
SRR35281602 SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186
SRR6797722 SRR6797801""".split())

# === Load tree ===
print("Loading tree...")
tree = Phylo.read(TREE_FILE, "newick")
tips_before = len(tree.get_terminals())
print(f"Tree has {tips_before} tips")

# === Prune basal strains ===
print("Pruning basal strains...")
for tip in tree.get_terminals():
    sra = tip.name.replace('L4.11_', '')
    if sra in BASAL_TO_PRUNE:
        tree.prune(tip)
        print(f"  Pruned: {tip.name}")

tips_after = len(tree.get_terminals())
print(f"Tree now has {tips_after} tips")

# === Midpoint rooting ===
print("Midpoint rooting...")
tree.root_at_midpoint()
tree.ladderize()

# === Assign colors ===
def get_sra(name):
    return name.replace('L4.11_', '')

colors = {}
for tip in tree.get_terminals():
    sra = get_sra(tip.name)
    if sra in L4111:
        colors[tip.name] = '#1f77b4'  # blue
    else:
        colors[tip.name] = '#d62728'  # red

n_blue = sum(1 for v in colors.values() if v == '#1f77b4')
n_red = sum(1 for v in colors.values() if v == '#d62728')
print(f"L4.11.1 (blue): {n_blue}, L4.11.2 (red): {n_red}")

# === Get tips in traversal order ===
tips_info = []
for tip in tree.get_terminals():
    depth = tree.distance(tip)
    tips_info.append((tip.name, depth))

n_tips = len(tips_info)
max_depth = max(d for _, d in tips_info)

# Assign angular positions
angles = {}
for i, (name, _) in enumerate(tips_info):
    angles[name] = 2 * np.pi * i / n_tips

# === Render circular phylogeny ===
print("Rendering circular phylogeny...")
fig, ax = plt.subplots(figsize=(14, 14), subplot_kw=dict(polar=True))

# Draw tips as colored dots
for name, depth in tips_info:
    angle = angles[name]
    r = depth / max_depth
    ax.plot(angle, r, 'o', color=colors[name], markersize=1.8, alpha=0.85)

# Draw branches recursively
def draw_clade(clade, parent_depth=0):
    terminals = clade.get_terminals()
    if clade.is_terminal():
        depth = tree.distance(clade)
        r = depth / max_depth
        angle = angles[clade.name]
        ax.plot([angle, angle], [parent_depth/max_depth, r],
                color=colors[clade.name], linewidth=0.3, alpha=0.6)
        return

    depth = tree.distance(clade)
    r = depth / max_depth
    r_parent = parent_depth / max_depth

    clade_tips = [t.name for t in terminals]
    clade_angles = [angles[t] for t in clade_tips if t in angles]
    if not clade_angles:
        return

    mean_angle = np.arctan2(
        np.mean([np.sin(a) for a in clade_angles]),
        np.mean([np.cos(a) for a in clade_angles])
    )

    # Radial line from parent to this node
    ax.plot([mean_angle, mean_angle], [r_parent, r],
            color='#333333', linewidth=0.3, alpha=0.5)

    # Arc connecting children
    child_angles = []
    for child in clade.clades:
        child_tips = [t.name for t in child.get_terminals()]
        child_clade_angles = [angles[t] for t in child_tips if t in angles]
        if child_clade_angles:
            child_mean = np.arctan2(
                np.mean([np.sin(a) for a in child_clade_angles]),
                np.mean([np.cos(a) for a in child_clade_angles])
            )
            child_angles.append(child_mean)

    if len(child_angles) >= 2:
        arc_angles = np.linspace(min(child_angles), max(child_angles), 50)
        ax.plot(arc_angles, [r]*50, color='#333333', linewidth=0.3, alpha=0.5)

    for child in clade.clades:
        child_tips = [t.name for t in child.get_terminals()]
        if child_tips:
            draw_clade(child, depth)

draw_clade(tree.root)

# Clean up
ax.set_ylim(0, 1.05)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['polar'].set_visible(False)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=10, label=f'L4.11.1 (n={n_blue})'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#d62728', markersize=10, label=f'L4.11.2 (n={n_red})'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=11,
          bbox_to_anchor=(1.1, -0.05))

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_PDF, bbox_inches='tight', facecolor='white')
print(f"Saved: {OUTPUT_PNG}")
print(f"Saved: {OUTPUT_PDF}")
print("Done!")
