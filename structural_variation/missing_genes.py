#!/usr/bin/env python3
"""Analyse missing_genes and gene coverage from report.json: potential deletions."""

import json, os, re
from collections import Counter, defaultdict

BDD = "BDD/L4.11"
GFF3 = "resources/NC_000962.3.gff3"

L4111 = set("ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448 ERR4553470 ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568 ERR4553613 ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721 ERR4553747 ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834 ERR4553841 ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942 ERR4553947 ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739 SRR21661641 SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865 SRR29016881 SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138 SRR29017144 SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234 SRR29017242 SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583 SRR29017611 SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674 SRR29017683 SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528 SRR29440700 SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600 SRR35281602 SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186 SRR6797722 SRR6797801".split())

# Load GFF3 for gene annotation
gene_products = {}
with open(GFF3) as f:
    for line in f:
        if line.startswith("#") or "\t" not in line:
            continue
        parts = line.strip().split("\t")
        if parts[2] != "CDS":
            continue
        attrs = parts[8]
        lt = re.search(r'locus_tag=([^;]+)', attrs)
        prod = re.search(r'product=([^;]+)', attrs)
        gene = re.search(r'gene=([^;]+)', attrs)
        if lt:
            tag = lt.group(1)
            p = prod.group(1).replace("%2C", ",") if prod else ""
            g = gene.group(1) if gene else ""
            gene_products[tag] = (g, p)

# ── 1. MISSING GENES ──
l1_missing = Counter()
l2_missing = Counter()
l1_missing_per_strain = []
l2_missing_per_strain = []
total_l1 = total_l2 = 0

# ── 2. GENES WITH LOW COVERAGE (percent_missing > 50%) ──
l1_low_cov = Counter()  # locus -> count of strains with >50% missing
l2_low_cov = Counter()

for strain_id in sorted(os.listdir(BDD)):
    rpath = os.path.join(BDD, strain_id, "NC_000962.3", "report.json")
    if not os.path.exists(rpath):
        continue
    sl = "L4.11.1" if strain_id in L4111 else "L4.11.2"
    if sl == "L4.11.1":
        total_l1 += 1
    else:
        total_l2 += 1
    
    try:
        with open(rpath) as f:
            data = json.load(f)
    except:
        continue
    
    # Missing genes
    mg = data.get("missing_genes", [])
    if sl == "L4.11.1":
        l1_missing_per_strain.append(len(mg))
        for g in mg:
            l1_missing[g] += 1
    else:
        l2_missing_per_strain.append(len(mg))
        for g in mg:
            l2_missing[g] += 1
    
    # Genes with high percent_missing (potential partial deletions)
    for gene_entry in data.get("genes", []):
        pct = gene_entry.get("percent_missing", 0)
        lt = gene_entry.get("locus_tag", "")
        if pct > 50:
            if sl == "L4.11.1":
                l1_low_cov[lt] += 1
            else:
                l2_low_cov[lt] += 1

print(f"Processed: L4.11.1={total_l1}, L4.11.2={total_l2}")

# ── RESULTS ──
print("\n" + "=" * 100)
print("1. MISSING GENES COUNT PER STRAIN")
print("=" * 100)
if l1_missing_per_strain:
    print(f"  L4.11.1: mean={sum(l1_missing_per_strain)/len(l1_missing_per_strain):.1f}, "
          f"range={min(l1_missing_per_strain)}-{max(l1_missing_per_strain)}")
    for n, c in sorted(Counter(l1_missing_per_strain).items()):
        print(f"    {n:3d} missing genes: {c:3d} strains")
if l2_missing_per_strain:
    print(f"  L4.11.2: mean={sum(l2_missing_per_strain)/len(l2_missing_per_strain):.1f}, "
          f"range={min(l2_missing_per_strain)}-{max(l2_missing_per_strain)}")
    for n, c in sorted(Counter(l2_missing_per_strain).items()):
        if c >= 3:
            print(f"    {n:3d} missing genes: {c:3d} strains")

print("\n" + "=" * 100)
print("2. MISSING GENES BY FREQUENCY (in >=3 strains)")
print("=" * 100)
all_missing = set(l1_missing.keys()) | set(l2_missing.keys())
# Sort by total frequency
freq_list = []
for g in all_missing:
    n1, n2 = l1_missing.get(g, 0), l2_missing.get(g, 0)
    if n1 + n2 >= 3:
        gname, product = gene_products.get(g, ("", ""))
        freq_list.append((g, gname, product, n1, n2))

freq_list.sort(key=lambda x: -(x[3]+x[4]))
print(f"{'Locus':12s} {'Gene':10s} {'L1':>5s} {'L2':>5s} {'Total':>6s}  {'Product'}")
print("-" * 100)
for lt, gn, prod, n1, n2 in freq_list:
    flag = ""
    if n1 > 0 and n2 == 0:
        flag = " ◄L1"
    elif n2 > 0 and n1 == 0:
        flag = " ◄L2"
    elif n1 > 0.8 * total_l1 and n2 > 0.8 * total_l2:
        flag = " ★SHARED"
    print(f"{lt:12s} {gn:10s} {n1:5d} {n2:5d} {n1+n2:6d}  {prod[:55]}{flag}")

# ── 3. LINEAGE-SPECIFIC MISSING GENES ──
print("\n" + "=" * 100)
print("3. LINEAGE-SPECIFIC MISSING GENES")
print("=" * 100)
print("\n  == L4.11.1 ONLY (missing in >=3 L1 strains, 0 L2) ==")
for lt, gn, prod, n1, n2 in freq_list:
    if n1 >= 3 and n2 == 0:
        pct = n1 * 100 // total_l1
        print(f"    {lt:12s} {gn:10s} {n1:3d}/{total_l1} ({pct}%)  {prod[:55]}")

print("\n  == L4.11.2 ONLY (missing in >=3 L2 strains, 0 L1) ==")
for lt, gn, prod, n1, n2 in freq_list:
    if n2 >= 3 and n1 == 0:
        pct = n2 * 100 // total_l2
        print(f"    {lt:12s} {gn:10s} {n2:3d}/{total_l2} ({pct}%)  {prod[:55]}")

# ── 4. LOW COVERAGE GENES (>50% missing, freq >= 10 strains) ──
print("\n" + "=" * 100)
print("4. GENES WITH >50% COVERAGE LOSS (freq >= 10)")
print("=" * 100)
all_low = set(l1_low_cov.keys()) | set(l2_low_cov.keys())
low_list = []
for g in all_low:
    n1, n2 = l1_low_cov.get(g, 0), l2_low_cov.get(g, 0)
    if n1 + n2 >= 10:
        gname, product = gene_products.get(g, ("", ""))
        low_list.append((g, gname, product, n1, n2))

low_list.sort(key=lambda x: -(x[3]+x[4]))
print(f"{'Locus':12s} {'Gene':10s} {'L1':>5s} {'L2':>5s}  {'Product'}")
print("-" * 100)
for lt, gn, prod, n1, n2 in low_list:
    flag = ""
    if n1 > total_l1 * 0.5 and n2 < total_l2 * 0.1:
        flag = " ◄L1-specific"
    elif n2 > total_l2 * 0.5 and n1 < total_l1 * 0.1:
        flag = " ◄L2-specific"
    elif n1 > total_l1 * 0.5 and n2 > total_l2 * 0.5:
        flag = " ★SHARED"
    print(f"{lt:12s} {gn:10s} {n1:5d} {n2:5d}  {prod[:55]}{flag}")

# ── 5. LIPID/CELL ENVELOPE GENES CHECK ──
print("\n" + "=" * 100)
print("5. LIPID/ENVELOPE GENES IN MISSING OR LOW-COVERAGE")
print("=" * 100)
lipid_keywords = ["lipid", "lipase", "acyl", "fatty", "myco", "wax", "pks", "fas",
                   "mce", "mmp", "pim", "mbt", "pps", "fadD", "fadE", "lpp"]
for lt, gn, prod, n1, n2 in freq_list + low_list:
    if any(kw.lower() in (prod + gn).lower() for kw in lipid_keywords):
        if n1 + n2 >= 3:
            print(f"  {lt:12s} {gn:10s} L1={n1:3d} L2={n2:3d}  {prod[:60]}")
