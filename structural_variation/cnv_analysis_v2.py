#!/usr/bin/env python3
"""
CNV (Copy Number Variation) detection from TBannotator v2 gene coverage data.

For each strain, normalizes per-gene median coverage by the genome-wide median.
Genes with ratio ≥ 2.0 are candidate duplications; ratio ≤ 0.5 are candidate deletions.
Identifies lineage-specific CNV events shared by ≥50% of a sub-lineage.
"""
import json
import os
import sys
from collections import defaultdict, Counter
from pathlib import Path
from statistics import median, mean
import numpy as np

ROOT = Path(__file__).parent
REPORTS_DIR = ROOT / "reports_v2"
FIG_DIR = ROOT.parent / "article" / "figures"

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

DUP_THRESHOLD = 2.0   # ratio ≥ 2.0 → potential duplication
DEL_THRESHOLD = 0.3   # ratio ≤ 0.3 → potential deletion (stricter than missing_genes)
LINEAGE_FREQ = 0.50   # ≥ 50% of sub-lineage = lineage-associated CNV


def process_reports():
    """Extract per-gene coverage ratios for all strains."""
    files = sorted(REPORTS_DIR.glob("*.json"))
    print(f"Processing {len(files)} reports...")

    # gene -> {strain_id: ratio}
    all_ratios = {}
    gene_list = None
    strain_subs = {}

    for i, f in enumerate(files):
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(files)}...", file=sys.stderr)

        with open(f) as fh:
            report = json.load(fh)

        sid = report['strain_id']
        sub = 'L4.11.1' if sid in L4111 else 'L4.11.2'
        strain_subs[sid] = sub

        genes = report.get('genes', [])
        if not genes:
            continue

        # Get gene list from first report
        if gene_list is None:
            gene_list = [g['locus_tag'] for g in genes]

        # Compute genome-wide median coverage
        coverages = [g['median_coverage'] for g in genes if g['median_coverage'] > 0]
        if not coverages:
            continue
        genome_median = median(coverages)
        if genome_median < 5:  # skip very low-coverage strains
            continue

        # Compute per-gene ratio
        ratios = {}
        for g in genes:
            cov = g['median_coverage']
            ratio = cov / genome_median if genome_median > 0 else 0
            ratios[g['locus_tag']] = ratio

        all_ratios[sid] = ratios

    return all_ratios, strain_subs, gene_list


def find_cnv_events(all_ratios, strain_subs, gene_list):
    """Find genes with consistent CNV across sub-lineages."""
    n_l1 = sum(1 for s in strain_subs.values() if s == 'L4.11.1')
    n_l2 = sum(1 for s in strain_subs.values() if s == 'L4.11.2')

    print(f"\nStrains: {len(all_ratios)} total, L4.11.1={n_l1}, L4.11.2={n_l2}")
    print(f"Genes: {len(gene_list)}")
    print(f"Thresholds: DUP≥{DUP_THRESHOLD}×, DEL≤{DEL_THRESHOLD}×, lineage freq≥{LINEAGE_FREQ*100:.0f}%")

    results = []

    for gene in gene_list:
        dup_l1 = 0
        dup_l2 = 0
        del_l1 = 0
        del_l2 = 0
        ratios_l1 = []
        ratios_l2 = []

        for sid, ratios in all_ratios.items():
            r = ratios.get(gene, 0)
            sub = strain_subs[sid]
            if sub == 'L4.11.1':
                ratios_l1.append(r)
                if r >= DUP_THRESHOLD:
                    dup_l1 += 1
                if r <= DEL_THRESHOLD:
                    del_l1 += 1
            else:
                ratios_l2.append(r)
                if r >= DUP_THRESHOLD:
                    dup_l2 += 1
                if r <= DEL_THRESHOLD:
                    del_l2 += 1

        if not ratios_l1 or not ratios_l2:
            continue

        freq_dup_l1 = dup_l1 / len(ratios_l1)
        freq_dup_l2 = dup_l2 / len(ratios_l2)
        freq_del_l1 = del_l1 / len(ratios_l1)
        freq_del_l2 = del_l2 / len(ratios_l2)
        med_l1 = median(ratios_l1)
        med_l2 = median(ratios_l2)

        # Flag if any sub-lineage has ≥50% of strains with CNV
        is_interesting = (
            freq_dup_l1 >= LINEAGE_FREQ or freq_dup_l2 >= LINEAGE_FREQ or
            freq_del_l1 >= LINEAGE_FREQ or freq_del_l2 >= LINEAGE_FREQ
        )

        if is_interesting:
            cnv_type = []
            if freq_dup_l1 >= LINEAGE_FREQ:
                cnv_type.append(f"DUP_L1({freq_dup_l1*100:.0f}%)")
            if freq_dup_l2 >= LINEAGE_FREQ:
                cnv_type.append(f"DUP_L2({freq_dup_l2*100:.0f}%)")
            if freq_del_l1 >= LINEAGE_FREQ:
                cnv_type.append(f"DEL_L1({freq_del_l1*100:.0f}%)")
            if freq_del_l2 >= LINEAGE_FREQ:
                cnv_type.append(f"DEL_L2({freq_del_l2*100:.0f}%)")

            results.append({
                'gene': gene,
                'med_l1': med_l1,
                'med_l2': med_l2,
                'dup_l1': dup_l1, 'dup_l2': dup_l2,
                'del_l1': del_l1, 'del_l2': del_l2,
                'freq_dup_l1': freq_dup_l1, 'freq_dup_l2': freq_dup_l2,
                'freq_del_l1': freq_del_l1, 'freq_del_l2': freq_del_l2,
                'n_l1': len(ratios_l1), 'n_l2': len(ratios_l2),
                'cnv_type': ', '.join(cnv_type),
                'mean_l1': mean(ratios_l1),
                'mean_l2': mean(ratios_l2),
            })

    results.sort(key=lambda x: max(x['med_l1'], x['med_l2']), reverse=True)
    return results


def make_cnv_figure(results, all_ratios, strain_subs, gene_list):
    """Create a CNV visualization."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not results:
        print("No significant CNV events found.")
        return

    # --- Figure 1: Summary bar chart of CNV genes ---
    dup_genes = [r for r in results if 'DUP' in r['cnv_type']]
    del_genes = [r for r in results if 'DEL' in r['cnv_type']]

    print(f"\n{'='*80}")
    print(f"CNV RESULTS SUMMARY")
    print(f"{'='*80}")
    print(f"Genes with lineage-associated duplications: {len(dup_genes)}")
    print(f"Genes with lineage-associated deletions: {len(del_genes)}")

    # Print top duplications
    print(f"\n--- TOP DUPLICATIONS (ratio ≥ {DUP_THRESHOLD}) ---")
    print(f"{'Gene':<15} {'Type':<25} {'Med L1':>8} {'Med L2':>8} {'DupL1%':>8} {'DupL2%':>8}")
    print("-" * 80)
    for r in sorted(dup_genes, key=lambda x: max(x['med_l1'], x['med_l2']), reverse=True)[:30]:
        print(f"{r['gene']:<15} {r['cnv_type']:<25} {r['med_l1']:>8.2f} {r['med_l2']:>8.2f} "
              f"{r['freq_dup_l1']*100:>7.1f}% {r['freq_dup_l2']*100:>7.1f}%")

    # Print top deletions
    print(f"\n--- TOP DELETIONS (ratio ≤ {DEL_THRESHOLD}) ---")
    print(f"{'Gene':<15} {'Type':<25} {'Med L1':>8} {'Med L2':>8} {'DelL1%':>8} {'DelL2%':>8}")
    print("-" * 80)
    for r in sorted(del_genes, key=lambda x: min(x['med_l1'], x['med_l2']))[:30]:
        print(f"{r['gene']:<15} {r['cnv_type']:<25} {r['med_l1']:>8.2f} {r['med_l2']:>8.2f} "
              f"{r['freq_del_l1']*100:>7.1f}% {r['freq_del_l2']*100:>7.1f}%")

    # --- Figure: CNV ratio comparison between sub-lineages ---
    fig, axes = plt.subplots(1, 2, figsize=(16, max(6, len(dup_genes[:25]) * 0.35)), dpi=200)

    # Panel 1: Duplications
    ax = axes[0]
    top_dup = sorted(dup_genes, key=lambda x: max(x['med_l1'], x['med_l2']), reverse=True)[:25]
    if top_dup:
        genes_d = [r['gene'] for r in top_dup]
        y = range(len(genes_d))
        ax.barh([yi - 0.15 for yi in y], [r['med_l1'] for r in top_dup],
                height=0.3, color='#1976d2', label='L4.11.1', alpha=0.85)
        ax.barh([yi + 0.15 for yi in y], [r['med_l2'] for r in top_dup],
                height=0.3, color='#d32f2f', label='L4.11.2', alpha=0.85)
        ax.axvline(x=DUP_THRESHOLD, color='grey', linestyle='--', alpha=0.7, label=f'DUP threshold ({DUP_THRESHOLD}×)')
        ax.axvline(x=1.0, color='black', linestyle='-', alpha=0.3)
        ax.set_yticks(list(y))
        ax.set_yticklabels(genes_d, fontsize=8)
        ax.set_xlabel('Median coverage ratio (gene/genome)', fontsize=11)
        ax.set_title(f'Top {len(top_dup)} duplicated genes', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

    # Panel 2: Deletions
    ax = axes[1]
    top_del = sorted(del_genes, key=lambda x: min(x['med_l1'], x['med_l2']))[:25]
    if top_del:
        genes_d = [r['gene'] for r in top_del]
        y = range(len(genes_d))
        ax.barh([yi - 0.15 for yi in y], [r['med_l1'] for r in top_del],
                height=0.3, color='#1976d2', label='L4.11.1', alpha=0.85)
        ax.barh([yi + 0.15 for yi in y], [r['med_l2'] for r in top_del],
                height=0.3, color='#d32f2f', label='L4.11.2', alpha=0.85)
        ax.axvline(x=DEL_THRESHOLD, color='grey', linestyle='--', alpha=0.7, label=f'DEL threshold ({DEL_THRESHOLD}×)')
        ax.axvline(x=1.0, color='black', linestyle='-', alpha=0.3)
        ax.set_yticks(list(y))
        ax.set_yticklabels(genes_d, fontsize=8)
        ax.set_xlabel('Median coverage ratio (gene/genome)', fontsize=11)
        ax.set_title(f'Top {len(top_del)} deleted genes (by coverage)', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

    plt.suptitle('Copy Number Variation in L4.11 sub-lineages\n(gene coverage ratio: gene median / genome median)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIG_DIR / "cnv_analysis.png"
    plt.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\n✅ CNV figure: {out}")


def main():
    all_ratios, strain_subs, gene_list = process_reports()
    results = find_cnv_events(all_ratios, strain_subs, gene_list)
    make_cnv_figure(results, all_ratios, strain_subs, gene_list)

    # Save results
    out_json = ROOT / "cnv_v2.json"
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"✅ Results: {out_json}")

    out_txt = ROOT / "cnv_v2_summary.txt"
    with open(out_txt, 'w') as f:
        f.write(f"CNV Analysis (v2, {len(all_ratios)} strains)\n")
        f.write(f"DUP threshold: ≥{DUP_THRESHOLD}×, DEL threshold: ≤{DEL_THRESHOLD}×\n")
        f.write(f"Lineage frequency: ≥{LINEAGE_FREQ*100:.0f}%\n\n")
        dup = [r for r in results if 'DUP' in r['cnv_type']]
        dele = [r for r in results if 'DEL' in r['cnv_type']]
        f.write(f"Duplications: {len(dup)} genes\n")
        f.write(f"Deletions: {len(dele)} genes\n\n")
        for r in results:
            f.write(f"{r['gene']}\t{r['cnv_type']}\tmed_L1={r['med_l1']:.2f}\tmed_L2={r['med_l2']:.2f}\n")
    print(f"✅ Summary: {out_txt}")


if __name__ == "__main__":
    main()
