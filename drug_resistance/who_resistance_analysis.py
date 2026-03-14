#!/usr/bin/env python3
"""
Cross-reference WHO Catalogue v2.1 mutations with L4.11 report.json SNPs.
Uses gene+HGVS from report.json annotations to match WHO GARC1 mutations.
Produces a full resistance profile across all 15 antibiotics.
"""

import csv, json, os, re, sys
from collections import Counter, defaultdict
from pathlib import Path

WHO_CSV = "/home/christophe/Documents/codes/MTBC/TBannotator/data/WHO_catalogue_v2.csv"
BDD_DIR = "/home/christophe/Documents/codes/MTBC/TBannotator/BDD/L4.11"
OUT_TXT = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/who_resistance_profile.txt"
OUT_PNG = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/who_resistance_heatmap.png"

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

DRUG_NAMES = {
    "INH": "Isoniazid", "RIF": "Rifampicin", "EMB": "Ethambutol", "PZA": "Pyrazinamide",
    "STM": "Streptomycin", "ETH": "Ethionamide", "KAN": "Kanamycin", "AMI": "Amikacin",
    "CAP": "Capreomycin", "LEV": "Levofloxacin", "MXF": "Moxifloxacin",
    "BDQ": "Bedaquiline", "CFZ": "Clofazimine", "DLM": "Delamanid", "LZD": "Linezolid"
}

# ── 1. Parse WHO catalogue ──────────────────────────────────────────────────
print("Loading WHO catalogue v2.1...")
who_mutations = {}   # (gene, mutation_str) → {drug: prediction}
who_by_gene = defaultdict(list)

with open(WHO_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        drug = row["DRUG"]
        mut_raw = row["MUTATION"]    # e.g. "rpoB@S450L" or "katG@S315T"
        pred = row["PREDICTION"]     # R, F, U, S
        
        if "@" not in mut_raw:
            continue
        
        gene, mut = mut_raw.split("@", 1)
        
        # Normalize gene name aliases
        gene_norm = gene.strip()
        
        key = (gene_norm, mut.strip())
        if key not in who_mutations:
            who_mutations[key] = {}
        who_mutations[key][drug] = pred
        who_by_gene[gene_norm].append((mut.strip(), drug, pred))

print(f"  {len(who_mutations)} unique (gene, mutation) entries")
print(f"  {len(who_by_gene)} genes tracked")
print(f"  Genes: {sorted(who_by_gene.keys())}")

# ── 2. Map HGVS protein notation to WHO GARC mutation ────────────────────────
def hgvs_p_to_garc(hgvs_p):
    """Convert SnpEff HGVS protein notation to WHO GARC format.
    e.g. 'p.Ser315Thr' → 'S315T'
    """
    if not hgvs_p or hgvs_p == "?" or not hgvs_p.startswith("p."):
        return None
    
    aa3 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Ter': '*', 'Stop': '*'
    }
    
    s = hgvs_p[2:]  # remove "p."
    
    # Handle frameshift: p.Ile67fs → I67fs
    m = re.match(r'([A-Z][a-z]{2})(\d+)fs', s)
    if m:
        ref_aa = aa3.get(m.group(1), '?')
        return f"{ref_aa}{m.group(2)}fs"
    
    # Simple substitution: p.Ser315Thr → S315T
    m = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', s)
    if m:
        ref_aa = aa3.get(m.group(1), '?')
        alt_aa = aa3.get(m.group(3), '?')
        return f"{ref_aa}{m.group(2)}{alt_aa}"
    
    # Stop gain: p.Arg202* → R202*
    m = re.match(r'([A-Z][a-z]{2})(\d+)\*', s)
    if m:
        ref_aa = aa3.get(m.group(1), '?')
        return f"{ref_aa}{m.group(2)}*"
    
    return None

# Also map gene names between SnpEff and WHO names
GENE_ALIASES = {
    "Rv0678": "mmpR5", "mmpR5": "mmpR5",
    "rrs": "rrs", "rrl": "rrl",
}

def normalize_gene(gene_name):
    """Normalize gene name to match WHO catalogue."""
    return GENE_ALIASES.get(gene_name, gene_name)

# ── 3. Scan report.json files ────────────────────────────────────────────────
print(f"\nScanning BDD/L4.11/ report.json files...")

# Drug order for display
DRUG_ORDER = ["INH", "RIF", "EMB", "PZA", "STM", "ETH", "KAN", "AMI", 
              "CAP", "LEV", "MXF", "BDQ", "CFZ", "DLM", "LZD"]

strain_results = {}   # strain → {drug: set of (prediction, gene, mutation)}
total_found = 0
total_strains = 0

for strain_dir in sorted(os.listdir(BDD_DIR)):
    report_path = os.path.join(BDD_DIR, strain_dir, "NC_000962.3", "report.json")
    if not os.path.exists(report_path):
        continue
    
    total_strains += 1
    
    try:
        with open(report_path) as f:
            data = json.load(f)
    except:
        continue
    
    strain_drugs = defaultdict(set)
    
    for snp in data.get("snp", []):
        for ann in snp.get("annotations", []):
            gene = ann.get("gene_name", "")
            hgvs_p = ann.get("hgvs_p", "")
            
            gene_norm = normalize_gene(gene)
            garc = hgvs_p_to_garc(hgvs_p)
            
            if garc and gene_norm in who_by_gene:
                key = (gene_norm, garc)
                if key in who_mutations:
                    for drug, pred in who_mutations[key].items():
                        if pred in ("R", "F"):
                            strain_drugs[drug].add((pred, gene_norm, garc))
                            total_found += 1
    
    strain_results[strain_dir] = dict(strain_drugs)
    
    if total_strains % 100 == 0:
        print(f"  {total_strains} strains processed...", file=sys.stderr)

print(f"  {total_strains} strains scanned (with report.json)")
print(f"  {total_found} WHO R/F mutations found across all strains")

# ── 4. Summary statistics ────────────────────────────────────────────────────
print(f"\n{'='*90}")
print(f"WHO CATALOGUE v2.1 — RESISTANCE PREDICTIONS FOR L4.11 ({total_strains} strains)")
print(f"{'='*90}")

n_l1 = sum(1 for s in strain_results if s in L4111)
n_l2 = sum(1 for s in strain_results if s not in L4111)
print(f"  L4.11.1: {n_l1} strains  |  L4.11.2: {n_l2} strains")

print(f"\n{'Drug':<6} {'Full name':<20} {'All':>8} {'%':>6} {'L4.11.1':>8} {'%':>6} {'L4.11.2':>8} {'%':>6}")
print("-" * 80)

drug_summary = []
for drug in DRUG_ORDER:
    r_all = sum(1 for s, res in strain_results.items() if drug in res)
    r_l1 = sum(1 for s, res in strain_results.items() if s in L4111 and drug in res)
    r_l2 = sum(1 for s, res in strain_results.items() if s not in L4111 and drug in res)
    
    pct_all = r_all / total_strains * 100 if total_strains > 0 else 0
    pct_l1 = r_l1 / n_l1 * 100 if n_l1 > 0 else 0
    pct_l2 = r_l2 / n_l2 * 100 if n_l2 > 0 else 0
    
    name = DRUG_NAMES.get(drug, drug)
    print(f"{drug:<6} {name:<20} {r_all:>6}/{total_strains:<4} {pct_all:>5.1f}% {r_l1:>6}/{n_l1:<3} {pct_l1:>5.1f}% {r_l2:>6}/{n_l2:<3} {pct_l2:>5.1f}%")
    drug_summary.append({
        "drug": drug, "name": name,
        "r_all": r_all, "pct_all": pct_all,
        "r_l1": r_l1, "pct_l1": pct_l1,
        "r_l2": r_l2, "pct_l2": pct_l2
    })

# ── 5. Top mutations per drug ────────────────────────────────────────────────
print(f"\n{'='*90}")
print(f"TOP MUTATIONS PER DRUG (WHO R/F only)")
print(f"{'='*90}")

for drug in DRUG_ORDER:
    mut_counter = Counter()
    for s, res in strain_results.items():
        if drug in res:
            for pred, gene, mut in res[drug]:
                mut_counter[f"{gene}@{mut} [{pred}]"] += 1
    
    if not mut_counter:
        continue
    
    print(f"\n  {drug} ({DRUG_NAMES.get(drug, drug)}):")
    for mut, cnt in mut_counter.most_common(8):
        print(f"    {mut:<35} {cnt:>5}/{total_strains} ({cnt/total_strains*100:.1f}%)")

# ── 6. Heatmap ───────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(8, 6), dpi=200, 
                         gridspec_kw={"width_ratios": [1, 1], "wspace": 0.4})

for ax_idx, (title, sub_filter) in enumerate([
    (f"L4.11.1 (n={n_l1})", lambda s: s in L4111),
    (f"L4.11.2 (n={n_l2})", lambda s: s not in L4111),
]):
    ax = axes[ax_idx]
    subset = {s: r for s, r in strain_results.items() if sub_filter(s)}
    n_sub = len(subset)
    
    pcts = []
    labels = []
    for drug in DRUG_ORDER:
        r = sum(1 for s, res in subset.items() if drug in res)
        pct = r / n_sub * 100 if n_sub > 0 else 0
        pcts.append(pct)
        labels.append(drug)
    
    colors = []
    for pct in pcts:
        if pct > 50:
            colors.append("#dc2626")
        elif pct > 20:
            colors.append("#f59e0b")
        elif pct > 5:
            colors.append("#3b82f6")
        elif pct > 0:
            colors.append("#9ca3af")
        else:
            colors.append("#e5e7eb")
    
    bars = ax.barh(range(len(labels)), pcts, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9, fontweight="bold")
    ax.set_xlabel("Resistance prevalence (%)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(0, 105)
    ax.invert_yaxis()
    ax.grid(True, axis="x", alpha=0.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    for bar, pct, drug in zip(bars, pcts, labels):
        if pct > 0:
            r = sum(1 for s, res in subset.items() if drug in res)
            ax.text(pct + 1.5, bar.get_y() + bar.get_height()/2, 
                    f"{pct:.1f}% ({r})", va="center", fontsize=7.5, fontweight="bold")

plt.suptitle("WHO Catalogue v2 — Drug resistance predictions\n(based on SNP genotypes in report.json)",
             fontsize=12, fontweight="bold", y=1.03)
plt.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
print(f"\nHeatmap: {OUT_PNG}")

# Save text output
with open(OUT_TXT, "w") as f:
    f.write(f"WHO Catalogue v2.1 resistance profile — L4.11\n")
    f.write(f"Strains: {total_strains} (L4.11.1: {n_l1}, L4.11.2: {n_l2})\n\n")
    for ds in drug_summary:
        f.write(f"{ds['drug']:<6} All: {ds['pct_all']:5.1f}%  L4.11.1: {ds['pct_l1']:5.1f}%  L4.11.2: {ds['pct_l2']:5.1f}%\n")

# Copy script
import shutil
shutil.copy(__file__, "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/who_resistance_analysis.py")

print(f"Results: {OUT_TXT}")
print("\nDone!")
