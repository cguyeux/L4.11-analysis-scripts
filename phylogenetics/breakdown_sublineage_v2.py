#!/usr/bin/env python3
"""
Sub-lineage breakdown of v2 analysis results.
Uses the L4.11.1 strain list from stat_country_sublineage.py and analysis_v2.json.

Produces: analysis_v2_sublineage.txt  (detailed comparison)
          analysis_v2_comparison.csv   (tabular data for article update)
"""
import json
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).parent

# L4.11.1 strain IDs (91 strains)
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


def fmt_pct(n, total):
    return f"{n}/{total} ({n/total*100:.1f}%)" if total > 0 else "0/0"


def main():
    # Load analysis results
    with open(ROOT / "analysis_v2.json") as f:
        data = json.load(f)

    # Split by sub-lineage
    l4111 = [d for d in data if d['strain_id'] in L4111]
    l4112 = [d for d in data if d['strain_id'] not in L4111]
    
    print(f"Total: {len(data)}")
    print(f"L4.11.1: {len(l4111)}")
    print(f"L4.11.2: {len(l4112)}")

    out = []
    out.append(f"V2 Sub-lineage Analysis ({len(data)} strains)")
    out.append(f"  L4.11.1: {len(l4111)}, L4.11.2: {len(l4112)}")
    out.append("=" * 70)

    # === Mapping stats ===
    out.append("\n### MAPPING STATS ###")
    for label, subset in [("L4.11.1", l4111), ("L4.11.2", l4112)]:
        depths = [d['mean_depth'] for d in subset if d['mean_depth'] > 0]
        covs = [d['coverage_pct'] for d in subset if d['coverage_pct'] > 0]
        out.append(f"  {label}: median depth = {statistics.median(depths):.1f}x " +
                   f"(range {min(depths):.1f}-{max(depths):.1f}), " +
                   f"median cov = {statistics.median(covs)*100:.1f}%")

    # === IS6110 ===
    out.append("\n### IS6110 ###")
    for label, subset in [("L4.11.1", l4111), ("L4.11.2", l4112)]:
        is6110s = [d['is6110_total'] for d in subset]
        out.append(f"  {label}: median = {statistics.median(is6110s):.0f} " +
                   f"(range {min(is6110s)}-{max(is6110s)})")

    # IS1081 for comparison
    for label, subset in [("L4.11.1", l4111), ("L4.11.2", l4112)]:
        is1081s = [d['is1081_total'] for d in subset]
        out.append(f"  {label} IS1081: median = {statistics.median(is1081s):.0f}")

    # === Missing genes ===
    out.append("\n### MISSING GENES ###")
    for label, subset in [("L4.11.1", l4111), ("L4.11.2", l4112)]:
        mg = [d['missing_genes_count'] for d in subset]
        out.append(f"  {label}: median = {statistics.median(mg):.0f} " +
                   f"(range {min(mg)}-{max(mg)})")

    # === Drug resistance (WHO) ===
    out.append("\n### DRUG RESISTANCE (WHO Catalogue v2) ###")
    drugs = ['INH', 'RIF', 'EMB', 'PZA', 'STM', 'LEV', 'MXF', 'ETH',
             'AMI', 'KAN', 'CAP', 'BDQ', 'CFZ', 'DLM', 'LZD']
    for drug in drugs:
        n_all = sum(1 for d in data if drug in d.get('who_resistance', {}))
        n_1 = sum(1 for d in l4111 if drug in d.get('who_resistance', {}))
        n_2 = sum(1 for d in l4112 if drug in d.get('who_resistance', {}))
        if n_all > 0:
            out.append(f"  {drug:>3}: ALL {fmt_pct(n_all, len(data)):>20}  " +
                       f"L4.11.1 {fmt_pct(n_1, len(l4111)):>20}  " +
                       f"L4.11.2 {fmt_pct(n_2, len(l4112)):>20}")

    # MDR (INH + RIF)
    out.append("\n  MDR (INH+RIF):")
    for label, subset in [("ALL", data), ("L4.11.1", l4111), ("L4.11.2", l4112)]:
        mdr = sum(1 for d in subset
                  if 'INH' in d.get('who_resistance', {})
                  and 'RIF' in d.get('who_resistance', {}))
        out.append(f"    {label}: {fmt_pct(mdr, len(subset))}")

    # === Key DR mutations ===
    out.append("\n### KEY DR MUTATIONS ###")
    key_mutations = [
        ('katG', 'S315T'), ('katG', 'S315N'),
        ('rpoB', 'S450L'), ('rpoB', 'D435Y'), ('rpoB', 'D435V'),
        ('embB', 'M306I'), ('embB', 'G406S'), ('embB', 'G406D'), ('embB', 'G406A'),
        ('rpsL', 'K43R'),
        ('gyrA', 'S91P'), ('gyrA', 'A90V'), ('gyrA', 'S95T'),
        ('pncA', 'D109N'),  # Actually panD
        ('rpoC', 'V483G'), ('rpoC', 'Q435P'), ('rpoC', 'L527V'),
        ('mmpR5', 'I67fs'),
        ('mmpL5', 'R202fs'),
        ('panD', 'D109N'),
        ('nat', 'G78D'), ('nat', 'D27N'),
    ]

    for gene, mut in key_mutations:
        n_all = 0
        n_1 = 0
        n_2 = 0
        for d in data:
            for dm in d.get('dr_mutations', []):
                if dm.get('gene') == gene and dm.get('garc') == mut:
                    n_all += 1
                    if d['strain_id'] in L4111:
                        n_1 += 1
                    else:
                        n_2 += 1
                    break
        if n_all > 0:
            out.append(f"  {gene} {mut}: ALL {fmt_pct(n_all, len(data)):>20}  " +
                       f"L4.11.1 {fmt_pct(n_1, len(l4111)):>20}  " +
                       f"L4.11.2 {fmt_pct(n_2, len(l4112)):>20}")

    # === Missing RD ===
    out.append("\n### KEY MISSING RD ###")
    key_rds = ['RD715', 'DS5', 'RD3', 'RD149']
    for rd in key_rds:
        n_1 = sum(1 for d in l4111 if rd in d.get('missing_rd', []))
        n_2 = sum(1 for d in l4112 if rd in d.get('missing_rd', []))
        out.append(f"  {rd}: L4.11.1 {fmt_pct(n_1, len(l4111)):>20}  " +
                   f"L4.11.2 {fmt_pct(n_2, len(l4112)):>20}")

    # === Total SNPs ===
    out.append("\n### TOTAL SNPs ###")
    for label, subset in [("L4.11.1", l4111), ("L4.11.2", l4112)]:
        snps = [d['total_snps'] for d in subset]
        out.append(f"  {label}: median = {statistics.median(snps):.0f} " +
                   f"(range {min(snps)}-{max(snps)})")

    # Print and save
    text = "\n".join(out)
    print(text)
    with open(ROOT / "analysis_v2_sublineage.txt", 'w') as f:
        f.write(text)
    print(f"\nSaved: {ROOT / 'analysis_v2_sublineage.txt'}")


if __name__ == "__main__":
    main()
