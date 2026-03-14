#!/usr/bin/env python3
"""Map IS6110 insertion positions to genes using GFF3 annotation."""

import os, re

# GFF3 file for H37Rv
GFF3 = "resources/NC_000962.3.gff3"

# IS6110 positions of interest (from previous analysis)
positions_of_interest = {
    # L4.11.2 exclusive core
    1960285: ("L4.11.2 exclusive", "366/368 (99%), 0/40 L1"),
    2038789: ("L4.11.2 exclusive", "360/368 (97%), 0/40 L1"),
    3549142: ("L4.11.2 exclusive", "350/368 (95%), 0/40 L1"),
    # L4.11.2 quasi-exclusive core 
    888836:  ("L4.11.2 quasi-excl", "355/368 (96%), 6/40 L1"),
    1996100: ("L4.11.2 quasi-excl", "352/368 (96%), 8/40 L1"),
    4212926: ("L4.11.2 quasi-excl", "366/368 (99%), 5/40 L1"),
    # Shared core
    1987702: ("SHARED core", "40/40 L1, 320/368 L2"),
    3120522: ("SHARED core", "40/40 L1, 360/368 L2"),
    # L4.11.1 specific
    1543432: ("L4.11.1 specific", "22/40 (55%), 0/368 L2"),
    # Additional L4.11.1 positions from the fingerprints
    889020:  ("L4.11.1 minor", "3/40 L1"),
    850086:  ("L4.11.1 minor", "2/40 L1"),
    2635576: ("L4.11.1 minor", "2/40 L1"),
    # L4.11.2 extra positions (from fingerprint patterns)
    1891:    ("L4.11.2 sub-clade", "31/368 L2"),
}

# Parse GFF3 to get gene annotations
genes = []  # (start, end, strand, locus_tag, gene_name, product)
with open(GFF3) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] not in ("gene", "CDS"):
            continue
        
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        
        locus_tag = ""
        gene_name = ""
        product = ""
        
        m = re.search(r'locus_tag=([^;]+)', attrs)
        if m:
            locus_tag = m.group(1)
        m = re.search(r'Name=([^;]+)', attrs)
        if m:
            gene_name = m.group(1)
        elif re.search(r'gene=([^;]+)', attrs):
            gene_name = re.search(r'gene=([^;]+)', attrs).group(1)
        m = re.search(r'product=([^;]+)', attrs)
        if m:
            product = m.group(1).replace("%2C", ",").replace("%3B", ";")
        
        if parts[2] == "gene":
            genes.append((start, end, strand, locus_tag, gene_name, product))

# For CDS, update product if missing
# Actually, let's build a product lookup from CDS
cds_products = {}
with open(GFF3) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        attrs = parts[8]
        lt = ""
        prod = ""
        m = re.search(r'locus_tag=([^;]+)', attrs)
        if m:
            lt = m.group(1)
        m = re.search(r'product=([^;]+)', attrs)
        if m:
            prod = m.group(1).replace("%2C", ",").replace("%3B", ";")
        if lt and prod:
            cds_products[lt] = prod

# Sort genes by start position
genes.sort(key=lambda x: x[0])

print(f"Loaded {len(genes)} genes from GFF3")
print(f"Loaded {len(cds_products)} CDS product annotations")

# Map each IS6110 position to gene(s)
WINDOW = 1355  # IS6110 is ~1355 bp long

print("\n" + "=" * 120)
print("IS6110 INSERTION POSITIONS → GENE MAPPING")
print("=" * 120)
print(f"{'Position':>12s}  {'Category':20s}  {'In Gene?':10s}  {'Locus':12s}  {'Gene':10s}  {'Product'}")
print("-" * 120)

for pos in sorted(positions_of_interest.keys()):
    cat, detail = positions_of_interest[pos]
    
    # Find genes overlapping with [pos, pos+1355]
    is_start = pos
    is_end = pos + WINDOW
    
    hit_genes = []
    upstream_gene = None
    downstream_gene = None
    
    for g_start, g_end, strand, locus, gname, prod in genes:
        # Direct overlap
        if g_start <= is_end and g_end >= is_start:
            product = cds_products.get(locus, prod)
            hit_genes.append((locus, gname, product, strand, g_start, g_end))
        # Track nearest upstream/downstream
        if g_end < is_start:
            upstream_gene = (locus, gname, cds_products.get(locus, prod), strand, g_start, g_end)
        if g_start > is_end and downstream_gene is None:
            downstream_gene = (locus, gname, cds_products.get(locus, prod), strand, g_start, g_end)
    
    if hit_genes:
        for locus, gname, product, strand, gs, ge in hit_genes:
            loc_type = "INSIDE" if gs <= is_start and ge >= is_end else "OVERLAP"
            print(f"{pos:>12d}  {cat:20s}  {loc_type:10s}  {locus:12s}  {gname:10s}  {product[:60]}")
    else:
        # Intergenic
        up_info = f"{upstream_gene[0]}({upstream_gene[1]})" if upstream_gene else "?"
        dn_info = f"{downstream_gene[0]}({downstream_gene[1]})" if downstream_gene else "?"
        up_dist = pos - upstream_gene[4] if upstream_gene else 0
        dn_dist = downstream_gene[4] - pos if downstream_gene else 0
        print(f"{pos:>12d}  {cat:20s}  INTERGENIC  between {up_info} (+{up_dist}bp) and {dn_info} (-{dn_dist}bp)")
    
    # Print detail on second line
    print(f"{'':>12s}  → {detail}")
    print()
