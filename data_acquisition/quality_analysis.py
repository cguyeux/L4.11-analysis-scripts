#!/usr/bin/env python3
"""Quality analysis v2: mapping_stats + covered_bases_percent + BioProject."""

import json, os
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = '"NC_000962.3:1700209:G:A"'

# BioProject from SRA metadata JSON if available
bp_map = {}
sra_meta = "data/sra_metadata"
if os.path.exists(sra_meta):
    for fn in os.listdir(sra_meta):
        if fn.endswith('.json'):
            try:
                with open(os.path.join(sra_meta, fn)) as f:
                    m = json.load(f)
                sid = m.get("run_accession", m.get("acc", fn.replace('.json','')))
                bp = m.get("bioproject", m.get("study_accession", ""))
                if sid and bp:
                    bp_map[sid] = bp
            except:
                pass

# Try to get BioProject from bioproject_analysis_L4.11 artifact
bp_artifact = "/home/christophe/.gemini/antigravity/brain/820bdbdc-4ce6-47a4-9944-b0abde9f4ddf/bioproject_analysis_L4.11.md"
if os.path.exists(bp_artifact):
    with open(bp_artifact) as f:
        content = f.read()

# Known BioProject prefixes 
known_bp = {
    "ERR1034": "PRJEB6273",   # CRyPTIC project
    "ERR2184": "PRJEB26924",  # Damien Foundation  
    "ERR2513": "PRJEB27979",
    "ERR2514": "PRJEB27979",
    "ERR3287": "PRJEB31443",  # CRyPTIC
    "ERR4553": "PRJEB38938",  # Bangladesh (KiteFin)
    "ERR4811": "PRJEB40777",  # CRyPTIC
    "ERR4812": "PRJEB40777",  # CRyPTIC
    "ERR4813": "PRJEB40777",  # CRyPTIC
    "ERR1328": "PRJEB33659",
    "SRR1049": "PRJNA234517",
    "SRR1062": "PRJNA234517",
    "SRR1140": "PRJNA234517",
    "SRR1080": "PRJNA595834", # TANDEM/Peru
    "SRR2166": "PRJNA749058", # Lima lipid
    "SRR2901": "PRJNA633333",
    "SRR2904": "PRJNA633333",
    "SRR3044": "PRJNA633333",
    "SRR3528": "PRJNA803760",
    "SRR4423": "PRJNA344785",
    "SRR6650": "PRJNA376914",
    "SRR6797": "PRJNA414593",
    "SRR2905": "PRJNA633333",
}

def get_bp(sid):
    """Get BioProject from known prefixes."""
    if sid in bp_map:
        return bp_map[sid]
    for prefix, bp in known_bp.items():
        if sid.startswith(prefix):
            return bp
    return "other"

stats = []
for sid in sorted(os.listdir(BDD)):
    rpath = os.path.join(BDD, sid, "NC_000962.3", "report.json")
    if not os.path.exists(rpath):
        continue
    try:
        with open(rpath) as f:
            raw = f.read()
        is_l1 = L4111_SPDI in raw
        data = json.loads(raw)
        ms = data.get("mapping_stats", {})
        
        stats.append({
            "sid": sid,
            "lineage": "L4.11.1" if is_l1 else "L4.11.2",
            "bp": get_bp(sid),
            "mean_depth": ms.get("mean_depth", 0),
            "covered_pct": ms.get("covered_bases_percent", 0),
            "mean_mapq": ms.get("mean_mapq", 0),
            "mean_baseq": ms.get("mean_baseq", 0),
            "reads": ms.get("reads_count", 0),
            "covered_bases": ms.get("covered_bases_count", 0),
            "n_missing": len(data.get("missing_genes", [])),
        })
    except:
        pass

print(f"Processed: {len(stats)} strains\n")

# ── 1. OVERALL ──
print("=" * 100)
print("1. OVERALL QUALITY METRICS")
print("=" * 100)
depths = sorted([s["mean_depth"] for s in stats])
covp = sorted([s["covered_pct"] for s in stats])
n = len(stats)
print(f"  Mean depth:   median={depths[n//2]:.1f}x, mean={sum(depths)/n:.1f}x, "
      f"range={depths[0]:.1f}-{depths[-1]:.1f}x")
print(f"  Covered %:    median={covp[n//2]:.2f}%, mean={sum(covp)/n:.2f}%, "
      f"range={covp[0]:.2f}-{covp[-1]:.2f}%")

# Depth buckets
print("\n  Depth distribution:")
bins = [(0,30,"<30x (low)"), (30,50,"30-50x"), (50,100,"50-100x"), 
        (100,200,"100-200x"), (200,500,"200-500x"), (500,9999,">500x")]
for lo, hi, label in bins:
    cnt = sum(1 for d in depths if lo <= d < hi)
    bar = "█" * (cnt * 40 // n) if n else ""
    print(f"    {label:15s}: {cnt:4d} ({cnt*100/n:5.1f}%) {bar}")

# ── 2. BY LINEAGE ──
print("\n" + "=" * 100)
print("2. BY SUB-LINEAGE")
print("=" * 100)
for lin in ["L4.11.1", "L4.11.2"]:
    ss = [s for s in stats if s["lineage"] == lin]
    d = sorted([s["mean_depth"] for s in ss])
    c = sorted([s["covered_pct"] for s in ss])
    m = sorted([s["n_missing"] for s in ss])
    nn = len(ss)
    print(f"\n  {lin} (n={nn}):")
    print(f"    Depth:    median={d[nn//2]:.1f}x, range={d[0]:.1f}-{d[-1]:.1f}x")
    print(f"    Covered:  median={c[nn//2]:.2f}%, range={c[0]:.2f}-{c[-1]:.2f}%")
    print(f"    Missing:  median={m[nn//2]}, range={m[0]}-{m[-1]}")

# ── 3. BY BIOPROJECT ──
print("\n" + "=" * 100)
print("3. BY BIOPROJECT")
print("=" * 100)
bp_g = defaultdict(list)
for s in stats:
    bp_g[s["bp"]].append(s)

print(f"\n{'BioProject':18s} {'n':>4s} {'MedDepth':>9s} {'MinD':>6s} {'MaxD':>6s} "
      f"{'MedCov%':>8s} {'MedMiss':>8s} {'L1':>3s} {'L2':>4s}")
print("-" * 80)
for bp, ss in sorted(bp_g.items(), key=lambda x: -len(x[1])):
    if len(ss) < 3:
        continue
    d = sorted([s["mean_depth"] for s in ss])
    c = sorted([s["covered_pct"] for s in ss])
    m = sorted([s["n_missing"] for s in ss])
    nn = len(ss)
    l1 = sum(1 for s in ss if s["lineage"] == "L4.11.1")
    print(f"{bp:18s} {nn:4d} {d[nn//2]:9.1f} {d[0]:6.1f} {d[-1]:6.1f} "
          f"{c[nn//2]:8.2f} {m[nn//2]:8d} {l1:3d} {nn-l1:4d}")

# ── 4. LOW-QUALITY ──
print("\n" + "=" * 100)
print("4. LOW-QUALITY STRAINS (depth < 30x or covered < 95%)")
print("=" * 100)
low = [s for s in stats if s["mean_depth"] < 30 or s["covered_pct"] < 95]
low.sort(key=lambda x: x["mean_depth"])
print(f"\n  {len(low)} strains below threshold")
print(f"\n  {'Strain':20s} {'Depth':>7s} {'Cov%':>7s} {'Miss':>5s} {'Lineage':>8s} {'BioProject'}")
print("  " + "-" * 75)
for s in low:
    print(f"  {s['sid']:20s} {s['mean_depth']:7.1f} {s['covered_pct']:7.2f} {s['n_missing']:5d} "
          f"{s['lineage']:>8s} {s['bp']}")
