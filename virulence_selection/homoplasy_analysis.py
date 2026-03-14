#!/usr/bin/env python3
"""Homoplasy analysis: find mutations present in both sub-lineages
at intermediate frequencies, suggesting independent emergence.

A putative homoplasy is an SPDI that:
- Is NOT ancestral (not fixed >90% in both sub-lineages)
- IS present in both sub-lineages at >5% frequency
- This pattern suggests convergent evolution under selective pressure
"""

import os, sys, json
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = "NC_000962.3:1700209:G:A"

# ── 1. LOAD SPDI DATA AND CLASSIFY ──
print("Loading SPDI data...", file=sys.stderr)
spdi_l1 = defaultdict(int)  # spdi -> count in L4.11.1
spdi_l2 = defaultdict(int)  # spdi -> count in L4.11.2
n1 = 0
n2 = 0

for sid in sorted(os.listdir(BDD)):
    spath = os.path.join(BDD, sid, "NC_000962.3", "spdi.txt")
    if not os.path.exists(spath):
        continue
    
    with open(spath) as f:
        spdis = set(line.strip() for line in f if line.strip())
    
    is_l1 = L4111_SPDI in spdis
    if is_l1:
        n1 += 1
        for s in spdis:
            spdi_l1[s] += 1
    else:
        n2 += 1
        for s in spdis:
            spdi_l2[s] += 1

print(f"L4.11.1: {n1} strains, L4.11.2: {n2} strains", file=sys.stderr)

# ── 2. IDENTIFY HOMOPLASIES ──
all_spdis = set(spdi_l1.keys()) | set(spdi_l2.keys())
print(f"Total unique SPDIs: {len(all_spdis)}", file=sys.stderr)

# Homoplasy criteria:
# - Present in >5% of BOTH sub-lineages
# - NOT fixed (>90%) in BOTH sub-lineages (these are ancestral)
MIN_FREQ = 0.05
MAX_ANCESTRAL = 0.90

homoplasies = []
for spdi in all_spdis:
    c1 = spdi_l1.get(spdi, 0)
    c2 = spdi_l2.get(spdi, 0)
    f1 = c1 / n1
    f2 = c2 / n2
    
    # Present in both at >5%
    if f1 >= MIN_FREQ and f2 >= MIN_FREQ:
        # Not ancestral (both >90%)
        if not (f1 > MAX_ANCESTRAL and f2 > MAX_ANCESTRAL):
            homoplasies.append({
                "spdi": spdi,
                "c1": c1, "c2": c2,
                "f1": f1, "f2": f2,
            })

print(f"Homoplasies found: {len(homoplasies)}", file=sys.stderr)

# ── 3. ANNOTATE WITH GENE INFO ──
# Load gene annotations from one report.json
print("Loading gene annotations...", file=sys.stderr)
spdi_to_gene = {}
spdi_to_effect = {}
spdi_to_hgvs = {}

# Sample a few report.json files for annotation
sample_sids = []
for sid in sorted(os.listdir(BDD)):
    rpath = os.path.join(BDD, sid, "NC_000962.3", "report.json")
    if os.path.exists(rpath):
        sample_sids.append(sid)
    if len(sample_sids) >= 20:
        break

for sid in sample_sids:
    rpath = os.path.join(BDD, sid, "NC_000962.3", "report.json")
    try:
        with open(rpath) as f:
            data = json.load(f)
        for snp in data.get("snp", []):
            spdi = snp.get("spdi", "")
            if spdi and spdi not in spdi_to_gene:
                anns = snp.get("annotations", [])
                if anns:
                    ann = anns[0]
                    spdi_to_gene[spdi] = ann.get("gene_name", "intergenic")
                    effects = ann.get("annotation", [])
                    spdi_to_effect[spdi] = effects[0] if effects else ""
                    spdi_to_hgvs[spdi] = ann.get("hgvs_p", "") or ""
    except:
        pass

print(f"Annotated {len(spdi_to_gene)} SPDIs", file=sys.stderr)

# ── 4. RESULTS ──
print(f"\n{'=' * 120}")
print(f"HOMOPLASY ANALYSIS: Mutations present in both L4.11.1 (n={n1}) and L4.11.2 (n={n2})")
print(f"Criteria: >5% in both sub-lineages, NOT fixed (>90%) in both")
print(f"{'=' * 120}")
print(f"\nTotal homoplasies: {len(homoplasies)}")

# Categorize
nonsyn = [h for h in homoplasies if spdi_to_effect.get(h["spdi"], "") in 
          ("missense_variant", "stop_gained", "frameshift_variant")]
syn = [h for h in homoplasies if spdi_to_effect.get(h["spdi"], "") == "synonymous_variant"]
interg = [h for h in homoplasies if "intergenic" in spdi_to_gene.get(h["spdi"], "intergenic")]

print(f"  Non-synonymous: {len(nonsyn)}")
print(f"  Synonymous: {len(syn)}")
print(f"  Intergenic/unannotated: {len(interg)}")

# ── 5. TOP NON-SYNONYMOUS HOMOPLASIES ──
print(f"\n{'=' * 120}")
print("TOP NON-SYNONYMOUS HOMOPLASIES (sorted by combined frequency)")
print(f"{'=' * 120}")

nonsyn_sorted = sorted(nonsyn, key=lambda x: -(x["f1"] + x["f2"]))
print(f"\n  {'Gene':20s} {'HGVS_p':20s} {'Effect':22s} {'L4.11.1':>15s} {'L4.11.2':>15s} {'SPDI'}")
print("  " + "-" * 130)

for h in nonsyn_sorted[:50]:
    gene = spdi_to_gene.get(h["spdi"], "?")
    hgvs = spdi_to_hgvs.get(h["spdi"], "")
    effect = spdi_to_effect.get(h["spdi"], "?")
    flag = ""
    # Flag if highly asymmetric (one lineage much higher)
    if h["f1"] > 0.5 and h["f2"] < 0.3:
        flag = " ◄ L4.11.1-enriched"
    elif h["f2"] > 0.5 and h["f1"] < 0.3:
        flag = " ◄ L4.11.2-enriched"
    elif h["f1"] > 0.3 and h["f2"] > 0.3:
        flag = " ★ TRUE CONVERGENCE"
    
    print(f"  {gene:20s} {hgvs:20s} {effect:22s} "
          f"{h['c1']:3d}/{n1} ({h['f1']*100:5.1f}%) "
          f"{h['c2']:3d}/{n2} ({h['f2']*100:5.1f}%){flag}")

# ── 6. GENE-LEVEL SUMMARY ──
print(f"\n{'=' * 120}")
print("GENES WITH MOST HOMOPLASTIC MUTATIONS")
print(f"{'=' * 120}")

gene_counts = defaultdict(int)
for h in nonsyn:
    gene = spdi_to_gene.get(h["spdi"], "?")
    gene_counts[gene] += 1

print(f"\n  {'Gene':20s} {'Homoplasies':>12s}")
print("  " + "-" * 35)
for gene, cnt in sorted(gene_counts.items(), key=lambda x: -x[1])[:30]:
    print(f"  {gene:20s} {cnt:12d}")

# ── 7. VIRULENCE GENE HOMOPLASIES ──
print(f"\n{'=' * 120}")
print("VIRULENCE GENE HOMOPLASIES")
print(f"{'=' * 120}")

import re
vir_patterns = {
    "PE/PPE": re.compile(r'^(PE|PPE)'),
    "mce": re.compile(r'^(mce|yrbE)'),
    "ESX": re.compile(r'^(ecc|esp|esx|mycP)'),
    "lipid": re.compile(r'^(fad|pks|mmp[SL]|mas$)'),
    "drug": re.compile(r'^(rpoB|rpoC|katG|inhA|embB|gyrA|pncA|ethA)'),
}

for fam, pat in vir_patterns.items():
    fam_homo = [h for h in nonsyn if pat.match(spdi_to_gene.get(h["spdi"], ""))]
    if fam_homo:
        print(f"\n  --- {fam} ({len(fam_homo)} homoplasies) ---")
        for h in sorted(fam_homo, key=lambda x: -(x["f1"] + x["f2"])):
            gene = spdi_to_gene.get(h["spdi"], "?")
            hgvs = spdi_to_hgvs.get(h["spdi"], "")
            flag = "★" if h["f1"] > 0.3 and h["f2"] > 0.3 else ""
            print(f"    {gene:20s} {hgvs:20s} "
                  f"L1={h['c1']:2d}/{n1} ({h['f1']*100:5.1f}%) "
                  f"L2={h['c2']:2d}/{n2} ({h['f2']*100:5.1f}%) {flag}")
