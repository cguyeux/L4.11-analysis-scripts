#!/usr/bin/env python3
"""
Regenerate all figures from v2 data (analysis_v2.json).
1. Co-resistance heatmap (4 first-line drugs)
2. WHO resistance heatmap (all drugs)
3. pN/pS analysis by gene family

Uses: analysis_v2.json (from analyse_reports_v2.py)
"""
import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).parent
FIG_DIR = ROOT.parent / "article" / "figures"
ANALYSIS = ROOT / "analysis_v2.json"

# L4.11.1 strain IDs
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

FIRST_LINE = ["INH", "RIF", "EMB", "PZA"]
ALL_DRUGS = ["INH", "RIF", "EMB", "PZA", "STM", "LEV", "MXF", "ETH",
             "AMI", "KAN", "CAP", "BDQ", "CFZ", "DLM", "LZD"]


def load_data():
    with open(ANALYSIS) as f:
        data = json.load(f)
    for d in data:
        d['sub'] = 'L4.11.1' if d['strain_id'] in L4111 else 'L4.11.2'
    return data


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CO-RESISTANCE HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════
def make_coresistance_heatmap(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    drug_labels = FIRST_LINE
    n_drugs = len(drug_labels)
    n_l1 = sum(1 for d in data if d['sub'] == 'L4.11.1')
    n_l2 = sum(1 for d in data if d['sub'] == 'L4.11.2')

    fig = plt.figure(figsize=(14, 5), dpi=200)
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.35)

    panels = [
        (f"All L4.11 (n={len(data)})", "All"),
        (f"L4.11.1 (n={n_l1})", "L4.11.1"),
        (f"L4.11.2 (n={n_l2})", "L4.11.2"),
    ]

    for idx, (title, sub_label) in enumerate(panels):
        ax = fig.add_subplot(gs[idx])
        subset = [d for d in data if sub_label == "All" or d['sub'] == sub_label]
        n = len(subset)

        co_pct = np.zeros((n_drugs, n_drugs))
        co_cnt = np.zeros((n_drugs, n_drugs), dtype=int)
        for i, d1 in enumerate(drug_labels):
            for j, d2 in enumerate(drug_labels):
                co = sum(1 for s in subset
                         if d1 in s.get('who_resistance', {})
                         and d2 in s.get('who_resistance', {}))
                co_cnt[i, j] = co
                co_pct[i, j] = co / n * 100 if n > 0 else 0

        cmap = "YlGnBu" if sub_label == "L4.11.1" else "YlOrRd"
        im = ax.imshow(co_pct, cmap=cmap, vmin=0,
                       vmax=max(co_pct.max(), 1), aspect="equal")

        ax.set_xticks(range(n_drugs))
        ax.set_yticks(range(n_drugs))
        ax.set_xticklabels(drug_labels, fontsize=10, fontweight="bold")
        ax.set_yticklabels(drug_labels, fontsize=10, fontweight="bold")
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

        for i in range(n_drugs):
            for j in range(n_drugs):
                val = co_pct[i, j]
                count = co_cnt[i, j]
                color = "white" if val > co_pct.max() * 0.6 else "black"
                ax.text(j, i, f"{val:.1f}%\n({count})", ha="center",
                        va="center", fontsize=8, color=color, fontweight="bold")

        plt.colorbar(im, ax=ax, shrink=0.8, label="Co-resistance (%)")

    plt.suptitle("Drug co-resistance patterns in L4.11 sub-lineages (v2, WHO Catalogue)",
                 fontsize=14, fontweight="bold", y=1.02)
    out = FIG_DIR / "coresistance_heatmap.png"
    plt.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"✅ Co-resistance heatmap: {out}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. WHO RESISTANCE HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════
def make_who_heatmap(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    l4111_data = [d for d in data if d['sub'] == 'L4.11.1']
    l4112_data = [d for d in data if d['sub'] == 'L4.11.2']

    # Compute prevalence per drug per sub-lineage
    drug_labels = ALL_DRUGS
    prev_l1 = []
    prev_l2 = []
    for drug in drug_labels:
        r1 = sum(1 for d in l4111_data if drug in d.get('who_resistance', {}))
        r2 = sum(1 for d in l4112_data if drug in d.get('who_resistance', {}))
        prev_l1.append(r1 / len(l4111_data) * 100 if l4111_data else 0)
        prev_l2.append(r2 / len(l4112_data) * 100 if l4112_data else 0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=200)

    def color_for_pct(p):
        if p > 50:
            return '#d32f2f'  # red
        elif p > 20:
            return '#ff9800'  # orange
        elif p > 5:
            return '#1976d2'  # blue
        elif p > 0:
            return '#90a4ae'  # grey-blue
        else:
            return '#e0e0e0'  # light grey

    for ax, prev, label, n_sub in [(ax1, prev_l1, 'L4.11.1', len(l4111_data)),
                                    (ax2, prev_l2, 'L4.11.2', len(l4112_data))]:
        colors = [color_for_pct(p) for p in prev]
        bars = ax.barh(range(len(drug_labels)), prev, color=colors, edgecolor='white')
        ax.set_yticks(range(len(drug_labels)))
        ax.set_yticklabels(drug_labels, fontsize=10, fontweight='bold')
        ax.set_xlabel('Strains with resistance (%)', fontsize=11)
        ax.set_title(f'{label} (n={n_sub})', fontsize=13, fontweight='bold')
        ax.set_xlim(0, 100)
        ax.invert_yaxis()

        for i, (p, bar) in enumerate(zip(prev, bars)):
            if p > 0:
                ax.text(p + 1, i, f'{p:.1f}%', va='center', fontsize=9,
                        fontweight='bold')

        # Add grid
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d32f2f', label='>50%'),
        Patch(facecolor='#ff9800', label='20-50%'),
        Patch(facecolor='#1976d2', label='5-20%'),
        Patch(facecolor='#90a4ae', label='<5%'),
        Patch(facecolor='#e0e0e0', label='0%'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=5,
               fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.suptitle('WHO Catalogue v2 — Resistance predictions by sub-lineage (822 strains)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    out = FIG_DIR / "who_resistance_heatmap.png"
    plt.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"✅ WHO resistance heatmap: {out}")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. pN/pS ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════
def make_pnps_analysis():
    print("\n=== pN/pS Analysis ===")
    pnps_file = ROOT / "pnps_v2.json"
    if not pnps_file.exists():
        print(f"  ❌ {pnps_file} not found")
        return

    with open(pnps_file) as f:
        pnps_data = json.load(f)

    # Aggregate: count unique (gene, hgvs_p) syn and nonsyn across all strains
    gene_syn = Counter()
    gene_nonsyn = Counter()
    family_syn = Counter()
    family_nonsyn = Counter()

    VIRULENCE_PREFIXES = {
        "PE_PGRS": "PE/PPE", "PPE": "PE/PPE", "PE": "PE/PPE",
        "mce": "mce", "ecc": "ESX", "esp": "ESX", "esx": "ESX",
        "fad": "Lipid", "pks": "Lipid", "mmpL": "Lipid",
        "fadD": "Lipid", "fadE": "Lipid",
        "vap": "TA", "maz": "TA", "rel": "TA",
        "sig": "Sigma",
        "dos": "DosR", "dev": "DosR", "hsp": "DosR",
        "emb": "DR genes", "rpo": "DR genes", "kat": "DR genes",
        "gyr": "DR genes", "rps": "DR genes", "pnc": "DR genes",
        "eth": "DR genes", "inh": "DR genes",
    }

    def classify(gene):
        for prefix, fam in sorted(VIRULENCE_PREFIXES.items(), key=lambda x: -len(x[0])):
            if gene.startswith(prefix):
                return fam
        return "Other"

    # Collect unique variants across ALL strains
    variants_seen = set()
    for sid, muts in pnps_data.items():
        for m in muts:
            gene = m.get('gene', '')
            hgvs = m.get('hgvs_p', '') or m.get('hgvs_c', '')
            is_syn = m.get('is_syn', False)
            is_nonsyn = m.get('is_nonsyn', False)
            key = (gene, hgvs, is_syn, is_nonsyn)
            if key not in variants_seen:
                variants_seen.add(key)
                fam = classify(gene)
                if is_syn:
                    gene_syn[gene] += 1
                    family_syn[fam] += 1
                if is_nonsyn:
                    gene_nonsyn[gene] += 1
                    family_nonsyn[fam] += 1

    # Gene-level pN/pS
    all_nonsyn = sum(family_nonsyn.values())
    all_syn = sum(family_syn.values())
    overall = all_nonsyn / all_syn if all_syn > 0 else 0
    print(f"  Overall: pN={all_nonsyn}, pS={all_syn}, pN/pS={overall:.2f}")
    print(f"  Genes with ≥1 variant: {len(gene_syn) + len(gene_nonsyn)}")

    # Family-level
    print(f"\n  {'Family':<15} {'pN':>6} {'pS':>6} {'pN/pS':>8}")
    print(f"  {'-'*40}")
    families_sorted = sorted(set(list(family_syn.keys()) + list(family_nonsyn.keys())))
    results = []
    for fam in families_sorted:
        ns = family_nonsyn.get(fam, 0)
        s = family_syn.get(fam, 0)
        ratio = ns / s if s > 0 else float('inf')
        print(f"  {fam:<15} {ns:>6} {s:>6} {ratio:>8.2f}")
        results.append((fam, ns, s, ratio))

    # Save
    with open(ROOT / "pnps_v2_summary.txt", 'w') as f:
        f.write(f"pN/pS Analysis (v2, 822 strains)\n")
        f.write(f"Overall: pN={all_nonsyn}, pS={all_syn}, pN/pS={overall:.2f}\n\n")
        f.write(f"{'Family':<15} {'pN':>6} {'pS':>6} {'pN/pS':>8}\n")
        for fam, ns, s, ratio in results:
            f.write(f"{fam:<15} {ns:>6} {s:>6} {ratio:>8.2f}\n")
    print(f"\n  Saved: {ROOT / 'pnps_v2_summary.txt'}")


def main():
    print("Loading v2 analysis data...")
    data = load_data()
    print(f"  {len(data)} strains loaded")

    make_coresistance_heatmap(data)
    make_who_heatmap(data)
    make_pnps_analysis()

    print("\n✅ All figures regenerated!")


if __name__ == "__main__":
    main()
