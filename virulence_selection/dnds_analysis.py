#!/usr/bin/env python3
"""Per-gene dN/dS ratio analysis across L4.11 strains.

Counts unique non-synonymous (dN: missense + stop_gained + frameshift) 
and synonymous (dS) variants per gene across the lineage.
Reports pN/pS ratio (proportion-based) for genes with sufficient data.
"""

import json, os, re
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = '"NC_000962.3:1700209:G:A"'

# Gene family classifier
FAMILIES = {
    "PE/PPE": re.compile(r'^(PE\d|PPE\d|PE_PGRS)'),
    "mce": re.compile(r'^(mce\d|yrbE)'),
    "ESX/T7SS": re.compile(r'^(ecc[A-E]|esp[A-K]|esx[A-Z]|mycP|EsxN)'),
    "lipid": re.compile(r'^(fad[A-Z]\d|pks\d|mas$|mmp[SL]|kas[AB])'),
    "drug_target": re.compile(r'^(rpoB|rpoC|katG|inhA|embB|gyrA|gyrB|rpsL|pncA|ethA)'),
    "sigma": re.compile(r'^sig[A-M]'),
    "TA_system": re.compile(r'^(vap[BC]|maz[EF]|rel[BE]|hig[AB]|par[DE])'),
}

def get_family(gene):
    for fam, pat in FAMILIES.items():
        if pat.match(gene):
            return fam
    return "other"

# Collect unique variants per gene (use SPDI as unique key)
gene_ns = defaultdict(set)   # gene -> set of unique non-synonymous SPDIs
gene_syn = defaultdict(set)  # gene -> set of unique synonymous SPDIs
gene_ns_l1 = defaultdict(set)
gene_syn_l1 = defaultdict(set)
gene_ns_l2 = defaultdict(set)
gene_syn_l2 = defaultdict(set)

total = 0
for sid in sorted(os.listdir(BDD)):
    rpath = os.path.join(BDD, sid, "NC_000962.3", "report.json")
    if not os.path.exists(rpath):
        continue
    total += 1
    
    try:
        with open(rpath) as f:
            raw = f.read()
        lin = "L4.11.1" if L4111_SPDI in raw else "L4.11.2"
        data = json.loads(raw)
        
        for snp in data.get("snp", []):
            spdi = snp.get("spdi", "")
            anns = snp.get("annotations", [])
            if not anns:
                continue
            ann = anns[0]
            gene = ann.get("gene_name", "")
            effects = ann.get("annotation", [])
            
            if not gene or "-" in gene:  # skip intergenic
                continue
            
            for eff in effects:
                if eff in ("missense_variant", "stop_gained", "frameshift_variant"):
                    gene_ns[gene].add(spdi)
                    if lin == "L4.11.1":
                        gene_ns_l1[gene].add(spdi)
                    else:
                        gene_ns_l2[gene].add(spdi)
                    break
                elif eff == "synonymous_variant":
                    gene_syn[gene].add(spdi)
                    if lin == "L4.11.1":
                        gene_syn_l1[gene].add(spdi)
                    else:
                        gene_syn_l2[gene].add(spdi)
                    break
    except:
        pass

print(f"Processed: {total} strains\n")

all_genes = set(gene_ns.keys()) | set(gene_syn.keys())
print(f"Genes with variants: {len(all_genes)}")
print(f"Total unique non-synonymous sites: {sum(len(v) for v in gene_ns.values())}")
print(f"Total unique synonymous sites: {sum(len(v) for v in gene_syn.values())}")

# Compute pN/pS per gene
results = []
for gene in all_genes:
    ns = len(gene_ns.get(gene, set()))
    syn = len(gene_syn.get(gene, set()))
    ns_l1 = len(gene_ns_l1.get(gene, set()))
    syn_l1 = len(gene_syn_l1.get(gene, set()))
    ns_l2 = len(gene_ns_l2.get(gene, set()))
    syn_l2 = len(gene_syn_l2.get(gene, set()))
    
    total_var = ns + syn
    if total_var < 2:
        continue
    
    # pN/pS ratio (avoid div by 0)
    ratio = ns / syn if syn > 0 else float('inf') if ns > 0 else 0
    ratio_l1 = ns_l1 / syn_l1 if syn_l1 > 0 else (float('inf') if ns_l1 > 0 else 0)
    ratio_l2 = ns_l2 / syn_l2 if syn_l2 > 0 else (float('inf') if ns_l2 > 0 else 0)
    
    fam = get_family(gene)
    
    results.append({
        "gene": gene,
        "family": fam,
        "ns": ns, "syn": syn, "total": total_var,
        "ratio": ratio,
        "ns_l1": ns_l1, "syn_l1": syn_l1,
        "ns_l2": ns_l2, "syn_l2": syn_l2,
        "ratio_l1": ratio_l1,
        "ratio_l2": ratio_l2,
    })

# ── 1. GLOBAL pN/pS BY GENE FAMILY ──
print(f"\n{'=' * 100}")
print("1. AGGREGATE pN/pS BY GENE FAMILY")
print(f"{'=' * 100}")

fam_stats = defaultdict(lambda: {"ns": 0, "syn": 0, "n_genes": 0})
for r in results:
    f = r["family"]
    fam_stats[f]["ns"] += r["ns"]
    fam_stats[f]["syn"] += r["syn"]
    fam_stats[f]["n_genes"] += 1

print(f"\n  {'Family':15s} {'Genes':>6s} {'dN':>6s} {'dS':>6s} {'pN/pS':>8s} {'Interpretation'}")
print("  " + "-" * 75)
for fam in ["PE/PPE", "mce", "ESX/T7SS", "lipid", "drug_target", "sigma", "TA_system", "other"]:
    s = fam_stats.get(fam, {"ns": 0, "syn": 0, "n_genes": 0})
    if s["n_genes"] == 0:
        continue
    ratio = s["ns"] / s["syn"] if s["syn"] > 0 else float('inf')
    interp = ""
    if ratio > 2.0:
        interp = "◄ POSITIVE selection"
    elif ratio > 1.0:
        interp = "◆ relaxed/diversifying"
    else:
        interp = "purifying"
    print(f"  {fam:15s} {s['n_genes']:6d} {s['ns']:6d} {s['syn']:6d} {ratio:8.2f} {interp}")

# ── 2. TOP GENES WITH HIGH pN/pS ──
print(f"\n{'=' * 100}")
print("2. TOP GENES WITH HIGHEST pN/pS (min 3 variants)")
print(f"{'=' * 100}")

high_pnps = sorted([r for r in results if r["total"] >= 3 and r["ratio"] > 1.0],
                    key=lambda x: (-x["ratio"], -x["total"]))
print(f"\n  {len(high_pnps)} genes with pN/pS > 1.0")
print(f"\n  {'Gene':20s} {'Family':>12s} {'dN':>5s} {'dS':>5s} {'pN/pS':>8s} {'Total':>6s}")
print("  " + "-" * 65)
for r in high_pnps[:40]:
    rat = f"{r['ratio']:.2f}" if r['ratio'] != float('inf') else "inf"
    print(f"  {r['gene']:20s} {r['family']:>12s} {r['ns']:5d} {r['syn']:5d} {rat:>8s} {r['total']:6d}")

# ── 3. GENES UNDER PURIFYING SELECTION (low pN/pS) ──
print(f"\n{'=' * 100}")
print("3. GENES UNDER STRONG PURIFYING SELECTION (pN/pS < 0.5, min 5 var)")
print(f"{'=' * 100}")
purifying = sorted([r for r in results if r["total"] >= 5 and r["ratio"] < 0.5],
                   key=lambda x: (x["ratio"], -x["total"]))
print(f"\n  {len(purifying)} genes")
print(f"\n  {'Gene':20s} {'Family':>12s} {'dN':>5s} {'dS':>5s} {'pN/pS':>8s}")
print("  " + "-" * 55)
for r in purifying[:20]:
    print(f"  {r['gene']:20s} {r['family']:>12s} {r['ns']:5d} {r['syn']:5d} {r['ratio']:8.2f}")

# ── 4. VIRULENCE GENES DETAIL ──
print(f"\n{'=' * 100}")
print("4. VIRULENCE GENE FAMILIES — PER-GENE pN/pS")
print(f"{'=' * 100}")

for fam in ["PE/PPE", "mce", "ESX/T7SS", "lipid", "sigma", "TA_system"]:
    fam_genes = sorted([r for r in results if r["family"] == fam], 
                       key=lambda x: (-x["ratio"], -x["total"]))
    if not fam_genes:
        continue
    print(f"\n  --- {fam} ---")
    print(f"  {'Gene':20s} {'dN':>5s} {'dS':>5s} {'pN/pS':>8s}  {'dN_L1':>5s} {'dS_L1':>5s}  {'dN_L2':>5s} {'dS_L2':>5s}")
    print("  " + "-" * 75)
    for r in fam_genes:
        rat = f"{r['ratio']:.2f}" if r['ratio'] != float('inf') else "inf"
        print(f"  {r['gene']:20s} {r['ns']:5d} {r['syn']:5d} {rat:>8s}"
              f"  {r['ns_l1']:5d} {r['syn_l1']:5d}  {r['ns_l2']:5d} {r['syn_l2']:5d}")

# ── 5. GENOME-WIDE SUMMARY ──
print(f"\n{'=' * 100}")
print("5. GENOME-WIDE SUMMARY")
print(f"{'=' * 100}")
total_ns = sum(r["ns"] for r in results)
total_syn = sum(r["syn"] for r in results)
total_ns_l1 = sum(r["ns_l1"] for r in results)
total_syn_l1 = sum(r["syn_l1"] for r in results)
total_ns_l2 = sum(r["ns_l2"] for r in results)
total_syn_l2 = sum(r["syn_l2"] for r in results)

print(f"\n  Genome-wide pN/pS: {total_ns}/{total_syn} = {total_ns/total_syn:.3f}")
print(f"  L4.11.1:           {total_ns_l1}/{total_syn_l1} = {total_ns_l1/total_syn_l1:.3f}" if total_syn_l1 else "")
print(f"  L4.11.2:           {total_ns_l2}/{total_syn_l2} = {total_ns_l2/total_syn_l2:.3f}" if total_syn_l2 else "")
