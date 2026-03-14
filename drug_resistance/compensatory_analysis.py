#!/usr/bin/env python3
"""Compensatory mutation analysis: rpoB/rpoC/rpoA co-occurrence.

Known compensatory mutations in rpoC (Comas et al. 2012, Song et al. 2014):
- rpoC: V483G, V483A, N698S, I491T, L516P, A172V, etc.
Known rpoA compensatory:
- rpoA: various positions

This script identifies ALL non-synonymous mutations in rpoB, rpoC, rpoA
and analyzes co-occurrence patterns.
"""

import json, os, sys
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = '"NC_000962.3:1700209:G:A"'

# Known compensatory mutations in rpoC (literature)
KNOWN_COMPENSATORY_RPOC = {
    "p.Val483Gly", "p.Val483Ala", "p.Asn698Ser", "p.Ile491Thr",
    "p.Leu516Pro", "p.Ala172Val", "p.Trp484Gly", "p.Asp485Tyr",
    "p.Asp485Asn", "p.Asn698His", "p.Gly594Glu", "p.Val1121Ala",
    "p.Ile491Val", "p.Pro1040Leu",
}

# Known major rpoB RRDR mutations
KNOWN_RPOB_RRDR = {
    "p.Ser450Leu", "p.His445Tyr", "p.His445Asp", "p.His445Asn",
    "p.Asp435Val", "p.Asp435Tyr", "p.Leu430Pro", "p.Ser450Trp",
    "p.His445Arg", "p.His445Leu", "p.Gln432Lys",
}

results = []
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
        
        rpoB_muts = []
        rpoC_muts = []
        rpoA_muts = []
        
        for snp in data.get("snp", []):
            anns = snp.get("annotations", [])
            if not anns:
                continue
            ann = anns[0]
            gene = ann.get("gene_name", "")
            effect = ann.get("annotation", [""])[0] if ann.get("annotation") else ""
            hgvs_p = ann.get("hgvs_p", "") or ""
            
            if gene == "rpoB" and ("missense" in effect or "stop" in effect or "frameshift" in effect):
                rpoB_muts.append(hgvs_p)
            elif gene == "rpoC" and ("missense" in effect or "stop" in effect or "frameshift" in effect):
                rpoC_muts.append(hgvs_p)
            elif gene == "rpoA" and ("missense" in effect or "stop" in effect or "frameshift" in effect):
                rpoA_muts.append(hgvs_p)
        
        results.append({
            "sid": sid,
            "lineage": lin,
            "rpoB": rpoB_muts,
            "rpoC": rpoC_muts,
            "rpoA": rpoA_muts,
            "has_RRDR": any(m in KNOWN_RPOB_RRDR for m in rpoB_muts),
            "has_comp_rpoC": any(m in KNOWN_COMPENSATORY_RPOC for m in rpoC_muts),
        })
    except Exception as e:
        pass

n = len(results)
n1 = sum(1 for r in results if r["lineage"] == "L4.11.1")
n2 = n - n1
print(f"Processed: {total} strains (L1={n1}, L2={n2})\n")

# ── 1. rpoB mutations ──
print("=" * 100)
print("1. rpoB MUTATIONS (rifampicin resistance)")
print("=" * 100)
rpoB_counts = defaultdict(lambda: {"L4.11.1": 0, "L4.11.2": 0})
for r in results:
    for m in r["rpoB"]:
        rpoB_counts[m][r["lineage"]] += 1

print(f"\n  {'Mutation':30s} {'L4.11.1':>15s} {'L4.11.2':>15s} {'Known RRDR'}")
print("  " + "-" * 75)
for m in sorted(rpoB_counts, key=lambda x: -(rpoB_counts[x]["L4.11.1"] + rpoB_counts[x]["L4.11.2"])):
    c = rpoB_counts[m]
    rrdr = "★ RRDR" if m in KNOWN_RPOB_RRDR else ""
    print(f"  {m:30s} {c['L4.11.1']:3d}/{n1} ({c['L4.11.1']*100/n1:5.1f}%) "
          f"{c['L4.11.2']:3d}/{n2} ({c['L4.11.2']*100/n2:5.1f}%) {rrdr}")

# ── 2. rpoC mutations ──
print(f"\n{'=' * 100}")
print("2. rpoC MUTATIONS (potential compensatory)")
print("=" * 100)
rpoC_counts = defaultdict(lambda: {"L4.11.1": 0, "L4.11.2": 0})
for r in results:
    for m in r["rpoC"]:
        rpoC_counts[m][r["lineage"]] += 1

print(f"\n  {'Mutation':30s} {'L4.11.1':>15s} {'L4.11.2':>15s} {'Known comp.'}")
print("  " + "-" * 75)
for m in sorted(rpoC_counts, key=lambda x: -(rpoC_counts[x]["L4.11.1"] + rpoC_counts[x]["L4.11.2"])):
    c = rpoC_counts[m]
    comp = "★ KNOWN" if m in KNOWN_COMPENSATORY_RPOC else ""
    print(f"  {m:30s} {c['L4.11.1']:3d}/{n1} ({c['L4.11.1']*100/n1:5.1f}%) "
          f"{c['L4.11.2']:3d}/{n2} ({c['L4.11.2']*100/n2:5.1f}%) {comp}")

# ── 3. rpoA mutations ──
print(f"\n{'=' * 100}")
print("3. rpoA MUTATIONS (potential compensatory)")
print("=" * 100)
rpoA_counts = defaultdict(lambda: {"L4.11.1": 0, "L4.11.2": 0})
for r in results:
    for m in r["rpoA"]:
        rpoA_counts[m][r["lineage"]] += 1

if rpoA_counts:
    print(f"\n  {'Mutation':30s} {'L4.11.1':>15s} {'L4.11.2':>15s}")
    print("  " + "-" * 65)
    for m in sorted(rpoA_counts, key=lambda x: -(rpoA_counts[x]["L4.11.1"] + rpoA_counts[x]["L4.11.2"])):
        c = rpoA_counts[m]
        print(f"  {m:30s} {c['L4.11.1']:3d}/{n1} ({c['L4.11.1']*100/n1:5.1f}%) "
              f"{c['L4.11.2']:3d}/{n2} ({c['L4.11.2']*100/n2:5.1f}%)")
else:
    print("\n  No rpoA mutations found")

# ── 4. CO-OCCURRENCE ANALYSIS ──
print(f"\n{'=' * 100}")
print("4. rpoB + rpoC CO-OCCURRENCE")
print("=" * 100)

has_rpoB = sum(1 for r in results if r["rpoB"])
has_rpoC = sum(1 for r in results if r["rpoC"])
has_both = sum(1 for r in results if r["rpoB"] and r["rpoC"])
has_rrdr = sum(1 for r in results if r["has_RRDR"])
has_rrdr_comp = sum(1 for r in results if r["has_RRDR"] and r["has_comp_rpoC"])
has_rrdr_any_rpoc = sum(1 for r in results if r["has_RRDR"] and r["rpoC"])

print(f"\n  Strains with any rpoB mutation:      {has_rpoB}/{n} ({has_rpoB*100/n:.1f}%)")
print(f"  Strains with any rpoC mutation:      {has_rpoC}/{n} ({has_rpoC*100/n:.1f}%)")
print(f"  Strains with both rpoB + rpoC:       {has_both}/{n} ({has_both*100/n:.1f}%)")
print(f"  Strains with RRDR rpoB mutation:     {has_rrdr}/{n} ({has_rrdr*100/n:.1f}%)")
print(f"  RRDR + known compensatory rpoC:      {has_rrdr_comp}/{has_rrdr} ({has_rrdr_comp*100/has_rrdr:.1f}%)" if has_rrdr else "")
print(f"  RRDR + any rpoC mutation:            {has_rrdr_any_rpoc}/{has_rrdr} ({has_rrdr_any_rpoc*100/has_rrdr:.1f}%)" if has_rrdr else "")

# ── 5. BY LINEAGE ──
print(f"\n{'=' * 100}")
print("5. CO-OCCURRENCE BY SUB-LINEAGE")
print("=" * 100)
for lin in ["L4.11.1", "L4.11.2"]:
    ss = [r for r in results if r["lineage"] == lin]
    nn = len(ss)
    rB = sum(1 for r in ss if r["rpoB"])
    rC = sum(1 for r in ss if r["rpoC"])
    both = sum(1 for r in ss if r["rpoB"] and r["rpoC"])
    rrdr = sum(1 for r in ss if r["has_RRDR"])
    comp = sum(1 for r in ss if r["has_RRDR"] and r["has_comp_rpoC"])
    print(f"\n  {lin} (n={nn}):")
    print(f"    rpoB mutant:    {rB}/{nn} ({rB*100/nn:.1f}%)")
    print(f"    RRDR mutant:    {rrdr}/{nn} ({rrdr*100/nn:.1f}%)")
    print(f"    rpoC mutant:    {rC}/{nn} ({rC*100/nn:.1f}%)")
    print(f"    RRDR + any rpoC:{sum(1 for r in ss if r['has_RRDR'] and r['rpoC'])}/{rrdr}" if rrdr else "")
    print(f"    RRDR + known comp:{comp}/{rrdr}" if rrdr else "")

# ── 6. DETAILED CO-OCCURRENCE TABLE ──
print(f"\n{'=' * 100}")
print("6. STRAIN-LEVEL rpoB/rpoC CO-OCCURRENCE (RRDR strains only)")
print("=" * 100)
rrdr_strains = [r for r in results if r["has_RRDR"]]
print(f"\n  {'Strain':20s} {'Lineage':>8s} {'rpoB':>25s} {'rpoC':>30s}")
print("  " + "-" * 90)
for r in sorted(rrdr_strains, key=lambda x: (x["lineage"], -len(x["rpoC"]))):
    rpoB_str = ",".join(r["rpoB"][:2])
    rpoC_str = ",".join(r["rpoC"][:3]) if r["rpoC"] else "—"
    comp_flag = " ★" if r["has_comp_rpoC"] else ""
    print(f"  {r['sid']:20s} {r['lineage']:>8s} {rpoB_str:>25s} {rpoC_str:>30s}{comp_flag}")
