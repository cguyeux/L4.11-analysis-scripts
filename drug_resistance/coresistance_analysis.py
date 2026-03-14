#!/usr/bin/env python3
"""
Co-resistance analysis for L4.11 strains.
- Resistance prevalence by sub-lineage
- Pairwise co-resistance matrix (Jaccard + counts)
- Fisher exact tests for resistance enrichment by sub-lineage
- Publication-quality heatmap (matplotlib/seaborn)
"""

import csv, re
import numpy as np
from collections import Counter
from scipy.stats import fisher_exact, chi2_contingency
import itertools

CSV = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv"
OUT_PNG = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/coresistance_heatmap.png"
OUT_TXT = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/stat_coresistance.txt"

DRUGS = ["Isoniazid", "Rifampicin", "Ethambutol", "Pyrazinamide"]
DRUG_SHORT = {"Isoniazid": "INH", "Rifampicin": "RIF", "Ethambutol": "EMB", "Pyrazinamide": "PZA"}

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

# ── Load data ────────────────────────────────────────────────────────────────
strains = []
with open(CSV) as f:
    for row in csv.DictReader(f):
        acc = row["Id"]
        sub = "L4.11.1" if acc in L4111 else "L4.11.2"
        dr = {d: (row.get(d, "N").strip().upper() == "Y") for d in DRUGS}
        strains.append({"id": acc, "sub": sub, **dr})

total = len(strains)
print(f"Total strains: {total}")

# ── Resistance prevalence ────────────────────────────────────────────────────
print(f"\n{'='*80}")
print(f"RESISTANCE PREVALENCE BY SUB-LINEAGE")
print(f"{'='*80}")
print(f"{'Drug':<15} {'Overall':>10} {'L4.11.1':>12} {'L4.11.2':>12} {'OR':>8} {'p(Fisher)':>12} {'Sig':>5}")
print("-" * 80)

results_dr = []
for sub_label in ["all", "L4.11.1", "L4.11.2"]:
    subset = [s for s in strains if sub_label == "all" or s["sub"] == sub_label]
    n = len(subset)
    for d in DRUGS:
        r = sum(1 for s in subset if s[d])
        if sub_label == "all":
            pass  # just count

n_l4111 = sum(1 for s in strains if s["sub"] == "L4.11.1")
n_l4112 = sum(1 for s in strains if s["sub"] == "L4.11.2")

for d in DRUGS:
    r_all = sum(1 for s in strains if s[d])
    r_l4111 = sum(1 for s in strains if s["sub"] == "L4.11.1" and s[d])
    r_l4112 = sum(1 for s in strains if s["sub"] == "L4.11.2" and s[d])
    
    # 2x2: resistant_yes/no × L4.11.1/L4.11.2
    a, b = r_l4111, r_l4112
    c, d2 = n_l4111 - r_l4111, n_l4112 - r_l4112
    table = np.array([[a, b], [c, d2]])
    or_val, p_fisher = fisher_exact(table)
    sig = "***" if p_fisher < 0.001 else "**" if p_fisher < 0.01 else "*" if p_fisher < 0.05 else "ns"
    
    pct_all = r_all / total * 100
    pct_l1 = r_l4111 / n_l4111 * 100
    pct_l2 = r_l4112 / n_l4112 * 100
    
    print(f"{DRUG_SHORT[d]:<15} {r_all:>4}/{total:<4} ({pct_all:4.1f}%)  {r_l4111:>3}/{n_l4111:<3} ({pct_l1:4.1f}%)  {r_l4112:>3}/{n_l4112:<3} ({pct_l2:4.1f}%)  {or_val:>8.2f} {p_fisher:>12.2e} {sig:>5}")
    results_dr.append({"drug": DRUG_SHORT[d], "r_all": r_all, "pct_all": pct_all,
                        "r_l1": r_l4111, "pct_l1": pct_l1, "r_l2": r_l4112, "pct_l2": pct_l2,
                        "OR": or_val, "p": p_fisher, "sig": sig})

# ── Resistance profiles ─────────────────────────────────────────────────────
print(f"\n{'='*80}")
print(f"RESISTANCE PROFILES")
print(f"{'='*80}")

for sub_label in ["L4.11.1", "L4.11.2", "All"]:
    subset = [s for s in strains if sub_label == "All" or s["sub"] == sub_label]
    n = len(subset)
    
    # Profile = tuple of resistant drugs
    profiles = Counter()
    for s in subset:
        r_drugs = tuple(DRUG_SHORT[d] for d in DRUGS if s[d])
        if r_drugs:
            profiles[r_drugs] += 1
        else:
            profiles[("Pan-susceptible",)] += 1
    
    print(f"\n  {sub_label} (n={n}):")
    for prof, cnt in profiles.most_common(10):
        label = "+".join(prof)
        mdr = "MDR" if "INH" in prof and "RIF" in prof else ""
        xdr = "XDR" if mdr and len(prof) >= 3 else ""
        tag = f"  [{xdr or mdr}]" if mdr else ""
        print(f"    {label:<30} {cnt:>5} ({cnt/n*100:5.1f}%){tag}")

# ── Co-resistance matrix (pairwise) ─────────────────────────────────────────
print(f"\n{'='*80}")
print(f"CO-RESISTANCE MATRICES")
print(f"{'='*80}")

drug_labels = [DRUG_SHORT[d] for d in DRUGS]
n_drugs = len(drug_labels)

for sub_label in ["All", "L4.11.1", "L4.11.2"]:
    subset = [s for s in strains if sub_label == "All" or s["sub"] == sub_label]
    n = len(subset)
    
    # Co-resistance count matrix
    co_matrix = np.zeros((n_drugs, n_drugs), dtype=int)
    for s in subset:
        for i, d1 in enumerate(DRUGS):
            for j, d2 in enumerate(DRUGS):
                if s[d1] and s[d2]:
                    co_matrix[i, j] += 1
    
    # Jaccard similarity matrix
    jaccard = np.zeros((n_drugs, n_drugs))
    for i in range(n_drugs):
        for j in range(n_drugs):
            union = sum(1 for s in subset if s[DRUGS[i]] or s[DRUGS[j]])
            if union > 0:
                jaccard[i, j] = co_matrix[i, j] / union
    
    print(f"\n  {sub_label} (n={n}) — Co-resistance counts:")
    print(f"    {'':>6}", end="")
    for dl in drug_labels:
        print(f" {dl:>6}", end="")
    print()
    for i, dl in enumerate(drug_labels):
        print(f"    {dl:>6}", end="")
        for j in range(n_drugs):
            print(f" {co_matrix[i,j]:>6}", end="")
        print()
    
    print(f"  Jaccard similarity:")
    print(f"    {'':>6}", end="")
    for dl in drug_labels:
        print(f" {dl:>6}", end="")
    print()
    for i, dl in enumerate(drug_labels):
        print(f"    {dl:>6}", end="")
        for j in range(n_drugs):
            print(f" {jaccard[i,j]:>6.3f}", end="")
        print()

# ── MDR/XDR/Pre-XDR summary ─────────────────────────────────────────────────
print(f"\n{'='*80}")
print(f"MDR/XDR SUMMARY")
print(f"{'='*80}")
for sub_label in ["All", "L4.11.1", "L4.11.2"]:
    subset = [s for s in strains if sub_label == "All" or s["sub"] == sub_label]
    n = len(subset)
    n_inh = sum(1 for s in subset if s["Isoniazid"])
    n_rif = sum(1 for s in subset if s["Rifampicin"])
    n_mdr = sum(1 for s in subset if s["Isoniazid"] and s["Rifampicin"])
    n_any = sum(1 for s in subset if any(s[d] for d in DRUGS))
    n_all4 = sum(1 for s in subset if all(s[d] for d in DRUGS))
    
    print(f"  {sub_label:>10} (n={n:>4}): Any-R={n_any:>4} ({n_any/n*100:.1f}%)  INH-R={n_inh:>4} ({n_inh/n*100:.1f}%)  RIF-R={n_rif:>4} ({n_rif/n*100:.1f}%)  MDR={n_mdr:>4} ({n_mdr/n*100:.1f}%)  Pan-R={n_all4:>4} ({n_all4/n*100:.1f}%)")

# ── Publication-quality heatmap ──────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(14, 5), dpi=200)
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.35)

panels = [
    ("All L4.11", "All"),
    ("L4.11.1 (n=%d)" % n_l4111, "L4.11.1"),
    ("L4.11.2 (n=%d)" % n_l4112, "L4.11.2"),
]

for idx, (title, sub_label) in enumerate(panels):
    ax = fig.add_subplot(gs[idx])
    subset = [s for s in strains if sub_label == "All" or s["sub"] == sub_label]
    n = len(subset)
    
    # Build co-resistance percentage matrix
    co_pct = np.zeros((n_drugs, n_drugs))
    for i, d1 in enumerate(DRUGS):
        for j, d2 in enumerate(DRUGS):
            co = sum(1 for s in subset if s[d1] and s[d2])
            co_pct[i, j] = co / n * 100
    
    # Heatmap
    cmap = "YlOrRd" if sub_label != "L4.11.1" else "YlGnBu"
    im = ax.imshow(co_pct, cmap=cmap, vmin=0, vmax=max(co_pct.max(), 1), aspect="equal")
    
    # Labels
    ax.set_xticks(range(n_drugs))
    ax.set_yticks(range(n_drugs))
    ax.set_xticklabels(drug_labels, fontsize=10, fontweight="bold")
    ax.set_yticklabels(drug_labels, fontsize=10, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    
    # Annotate cells
    for i in range(n_drugs):
        for j in range(n_drugs):
            val = co_pct[i, j]
            count = sum(1 for s in subset if s[DRUGS[i]] and s[DRUGS[j]])
            color = "white" if val > co_pct.max() * 0.6 else "black"
            ax.text(j, i, f"{val:.1f}%\n({count})", ha="center", va="center",
                    fontsize=8, color=color, fontweight="bold")
    
    plt.colorbar(im, ax=ax, shrink=0.8, label="Co-resistance (%)")

plt.suptitle("Drug co-resistance patterns in L4.11 sub-lineages", 
             fontsize=14, fontweight="bold", y=1.02)
plt.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
print(f"\nHeatmap: {OUT_PNG}")

# ── Save text results ────────────────────────────────────────────────────────
with open(OUT_TXT, "w") as f:
    f.write("Co-resistance analysis — L4.11 sub-lineages\n")
    f.write(f"N = {total} strains (L4.11.1: {n_l4111}, L4.11.2: {n_l4112})\n\n")
    f.write("Resistance prevalence:\n")
    for r in results_dr:
        f.write(f"  {r['drug']:<6} All:{r['pct_all']:5.1f}%  L4.11.1:{r['pct_l1']:5.1f}%  L4.11.2:{r['pct_l2']:5.1f}%  OR={r['OR']:.2f}  p={r['p']:.2e} {r['sig']}\n")

print(f"Results: {OUT_TXT}")
print("\nDone!")
