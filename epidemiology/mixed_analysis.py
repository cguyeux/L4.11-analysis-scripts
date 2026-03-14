#!/usr/bin/env python3
"""Mixed infection analysis: detect intra-host heterogeneity via allelic ratios."""

import json, os, sys
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = '"NC_000962.3:1700209:G:A"'

# Thresholds
HET_THRESHOLD = 0.1   # ref/(alt+ref) > 10% = heterozygous site
MIXED_MIN_SITES = 10  # Strains with >= 10 het sites flagged

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
        
        is_l1 = L4111_SPDI in raw
        data = json.loads(raw)
        
        snps = data.get("snp", [])
        n_snp = len(snps)
        
        # Count heterozygous-like sites
        het_sites = []
        for s in snps:
            ac = s.get("alt_count", 0) or 0
            rc = s.get("ref_count", 0) or 0
            dp = ac + rc
            if dp >= 10 and rc > 0:
                ratio = rc / dp
                if ratio >= HET_THRESHOLD:
                    het_sites.append({
                        "spdi": s.get("spdi", ""),
                        "alt": ac,
                        "ref": rc,
                        "depth": dp,
                        "ratio": ratio,
                        "ann": s.get("annotations", [{}])[0].get("gene_name", "") if s.get("annotations") else ""
                    })
        
        # Allelic ratio stats
        ratios = [h["ratio"] for h in het_sites]
        near_50 = sum(1 for r in ratios if 0.3 <= r <= 0.7)
        
        results.append({
            "sid": sid,
            "lineage": "L4.11.1" if is_l1 else "L4.11.2",
            "n_snp": n_snp,
            "n_het": len(het_sites),
            "n_near50": near_50,
            "het_pct": len(het_sites) * 100 / n_snp if n_snp else 0,
            "het_sites": het_sites,
        })
        
    except Exception as e:
        print(f"Error {sid}: {e}", file=sys.stderr)

print(f"Processed: {total} strains\n")

# ── 1. OVERVIEW ──
print("=" * 100)
print("1. HETEROZYGOSITY OVERVIEW")
print("=" * 100)
het_counts = sorted([r["n_het"] for r in results])
n = len(results)
print(f"  Het sites per strain: median={het_counts[n//2]}, "
      f"mean={sum(het_counts)/n:.1f}, range={het_counts[0]}-{het_counts[-1]}")
print(f"  Strains with 0 het sites: {sum(1 for h in het_counts if h == 0)}")
print(f"  Strains with >=10 het sites: {sum(1 for h in het_counts if h >= 10)}")
print(f"  Strains with >=20 het sites: {sum(1 for h in het_counts if h >= 20)}")
print(f"  Strains with >=50 het sites: {sum(1 for h in het_counts if h >= 50)}")

# Distribution
print("\n  Distribution:")
bins = [(0,1,"0"), (1,5,"1-4"), (5,10,"5-9"), (10,20,"10-19"), (20,50,"20-49"), (50,999,">=50")]
for lo, hi, label in bins:
    cnt = sum(1 for h in het_counts if lo <= h < hi)
    bar = "█" * (cnt * 40 // n) if n else ""
    print(f"    {label:8s}: {cnt:4d} ({cnt*100/n:5.1f}%) {bar}")

# ── 2. BY SUB-LINEAGE ──
print("\n" + "=" * 100)
print("2. BY SUB-LINEAGE")
print("=" * 100)
for lin in ["L4.11.1", "L4.11.2"]:
    ss = [r for r in results if r["lineage"] == lin]
    hc = sorted([r["n_het"] for r in ss])
    nn = len(ss)
    print(f"\n  {lin} (n={nn}):")
    print(f"    Het sites: median={hc[nn//2]}, mean={sum(hc)/nn:.1f}, max={hc[-1]}")
    print(f"    Mixed candidates (>=10): {sum(1 for h in hc if h >= 10)}")

# ── 3. MIXED INFECTION CANDIDATES ──
print("\n" + "=" * 100)
print("3. MIXED INFECTION CANDIDATES (>=10 heterozygous sites)")
print("=" * 100)
mixed = sorted([r for r in results if r["n_het"] >= MIXED_MIN_SITES], 
               key=lambda x: -x["n_het"])
print(f"\n  {len(mixed)} candidates")
print(f"\n  {'Strain':20s} {'SNPs':>5s} {'Het':>4s} {'~50%':>4s} {'%Het':>6s} {'Lineage':>8s}")
print("  " + "-" * 55)
for r in mixed:
    print(f"  {r['sid']:20s} {r['n_snp']:5d} {r['n_het']:4d} {r['n_near50']:4d} "
          f"{r['het_pct']:6.1f} {r['lineage']:>8s}")

# ── 4. TOP AFFECTED GENES ──
if mixed:
    print("\n" + "=" * 100)
    print("4. MOST FREQUENTLY HETEROZYGOUS GENES (across mixed candidates)")
    print("=" * 100)
    gene_counts = defaultdict(int)
    for r in mixed:
        for h in r["het_sites"]:
            g = h["ann"] or "intergenic"
            gene_counts[g] += 1
    
    for gene, cnt in sorted(gene_counts.items(), key=lambda x: -x[1])[:25]:
        print(f"    {gene:20s}: {cnt}")

# ── 5. NEAR 50:50 ANALYSIS ──
print("\n" + "=" * 100)
print("5. STRAINS WITH MANY NEAR-50:50 SITES (strong mixed signal)")
print("=" * 100)
strong = sorted([r for r in results if r["n_near50"] >= 5], 
                key=lambda x: -x["n_near50"])
if strong:
    print(f"\n  {len(strong)} strains with >=5 sites where ref ratio is 30-70%")
    for r in strong[:15]:
        print(f"  {r['sid']:20s} het={r['n_het']:3d} near50={r['n_near50']:3d} {r['lineage']}")
else:
    print("\n  No strains with >=5 near-50:50 sites")
