#!/usr/bin/env python3
"""Analyze IS and RD data across L4.11 sub-lineages — optimized version.
Uses jq for large report.json files, direct json for small is_report.json."""
import json, os, sys, subprocess
from collections import Counter
import numpy as np

BASE = '/home/christophe/Documents/codes/MTBC/TBannotator/BDD/L4.11'

L4111_IDS = set("""
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

strains = sorted(os.listdir(BASE))

# ══════════════════════════════════════════════════════════════════════════
# PART 1: IS ANALYSIS (small files, fast)
# ══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("IS6110 COPY NUMBER ANALYSIS")
print("=" * 70)

is6110 = {'L4.11.1': [], 'L4.11.2': []}
is_novel = {'L4.11.1': [], 'L4.11.2': []}
is_total = {'L4.11.1': [], 'L4.11.2': []}
n_is = {'L4.11.1': 0, 'L4.11.2': 0}

for strain in strains:
    isp = os.path.join(BASE, strain, 'NC_000962.3', 'is_report.json')
    if not os.path.exists(isp):
        continue
    sub = 'L4.11.1' if strain in L4111_IDS else 'L4.11.2'
    n_is[sub] += 1
    with open(isp) as f:
        data = json.load(f)
    ises = data.get('insertion_sequences', [])
    is_total[sub].append(len(ises))
    n6110 = sum(1 for e in ises if e['name'] == 'IS6110')
    nn = sum(1 for e in ises if e['name'] == 'IS6110' and not e['is_reference'])
    is6110[sub].append(n6110)
    is_novel[sub].append(nn)

for sub in ['L4.11.1', 'L4.11.2']:
    a = np.array(is6110[sub])
    nv = np.array(is_novel[sub])
    tot = np.array(is_total[sub])
    print(f"\n{sub} (n={n_is[sub]}):")
    print(f"  IS6110 total:    mean={a.mean():.1f} ± {a.std():.1f}, median={np.median(a):.0f}, range=[{a.min()}-{a.max()}]")
    print(f"  IS6110 novel:    mean={nv.mean():.1f} ± {nv.std():.1f}, median={np.median(nv):.0f}")
    print(f"  All IS elements: mean={tot.mean():.1f}")

# ══════════════════════════════════════════════════════════════════════════
# PART 2: RD + MISSING GENES (large files — use jq)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("MISSING RD & GENES ANALYSIS (via jq)")
print("=" * 70)

# Build list of report.json paths
report_files = []
for strain in strains:
    rp = os.path.join(BASE, strain, 'NC_000962.3', 'report.json')
    if os.path.exists(rp):
        report_files.append((strain, rp))

print(f"Processing {len(report_files)} report.json files with jq...")

missing_rd = {'L4.11.1': Counter(), 'L4.11.2': Counter()}
missing_genes_c = {'L4.11.1': Counter(), 'L4.11.2': Counter()}
n_rd = {'L4.11.1': 0, 'L4.11.2': 0}

for i, (strain, rp) in enumerate(report_files):
    sub = 'L4.11.1' if strain in L4111_IDS else 'L4.11.2'
    n_rd[sub] += 1
    
    # Extract missing_rd and missing_genes with jq
    result = subprocess.run(
        ['jq', '-c', '{missing_rd, missing_genes}', rp],
        capture_output=True, text=True, timeout=10
    )
    if result.returncode != 0:
        continue
    
    data = json.loads(result.stdout)
    for rd in data.get('missing_rd', []):
        missing_rd[sub][rd] += 1
    for g in data.get('missing_genes', []):
        missing_genes_c[sub][g] += 1
    
    if (i+1) % 50 == 0:
        print(f"  ...{i+1}/{len(report_files)}", file=sys.stderr)

print(f"\nL4.11.1: {n_rd['L4.11.1']} strains, L4.11.2: {n_rd['L4.11.2']} strains")

# Named RDs
all_rds = set(missing_rd['L4.11.1'].keys()) | set(missing_rd['L4.11.2'].keys())
named = sorted([r for r in all_rds if not r.startswith('CUS_GS_')])
custom = sorted([r for r in all_rds if r.startswith('CUS_GS_')])

print(f"\n--- Named RDs ({len(named)}) ---")
print(f"  {'RD':<30} {'L4.11.1':>20} {'L4.11.2':>20}")
print("  " + "-" * 72)
for rd in named:
    c1 = missing_rd['L4.11.1'].get(rd, 0)
    c2 = missing_rd['L4.11.2'].get(rd, 0)
    p1 = c1/n_rd['L4.11.1']*100 if n_rd['L4.11.1'] else 0
    p2 = c2/n_rd['L4.11.2']*100 if n_rd['L4.11.2'] else 0
    flag = ' ★' if abs(p1-p2) > 20 else ''
    print(f"  {rd:<28} {c1:3d}/{n_rd['L4.11.1']:3d} ({p1:5.1f}%)   {c2:3d}/{n_rd['L4.11.2']:3d} ({p2:5.1f}%){flag}")

# CUS_GS with significant differences
print(f"\n--- Custom regions (CUS_GS_*): {len(custom)} total ---")
diff_cus = []
for rd in custom:
    c1 = missing_rd['L4.11.1'].get(rd, 0)
    c2 = missing_rd['L4.11.2'].get(rd, 0)
    p1 = c1/n_rd['L4.11.1']*100 if n_rd['L4.11.1'] else 0
    p2 = c2/n_rd['L4.11.2']*100 if n_rd['L4.11.2'] else 0
    if abs(p1-p2) > 30:
        diff_cus.append((rd, c1, p1, c2, p2))

if diff_cus:
    print(f"  Regions with >30% difference:")
    for rd, c1, p1, c2, p2 in diff_cus:
        print(f"    {rd}: L4.11.1={p1:.0f}% vs L4.11.2={p2:.0f}%")
else:
    print("  No regions with >30% difference between sub-lineages")

# Missing genes
print(f"\n--- Missing Genes ---")
all_genes = set(missing_genes_c['L4.11.1'].keys()) | set(missing_genes_c['L4.11.2'].keys())
print(f"  Total distinct: {len(all_genes)}")
print(f"  {'Gene':<15} {'L4.11.1':>20} {'L4.11.2':>20}")
print("  " + "-" * 57)
for gene in sorted(all_genes):
    c1 = missing_genes_c['L4.11.1'].get(gene, 0)
    c2 = missing_genes_c['L4.11.2'].get(gene, 0)
    p1 = c1/n_rd['L4.11.1']*100 if n_rd['L4.11.1'] else 0
    p2 = c2/n_rd['L4.11.2']*100 if n_rd['L4.11.2'] else 0
    flag = ' ★' if abs(p1-p2) > 20 else ''
    print(f"  {gene:<13} {c1:3d}/{n_rd['L4.11.1']:3d} ({p1:5.1f}%)   {c2:3d}/{n_rd['L4.11.2']:3d} ({p2:5.1f}%){flag}")
