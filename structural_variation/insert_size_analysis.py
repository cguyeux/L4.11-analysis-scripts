#!/usr/bin/env python3
"""
Insert size analysis by BioProject from v2 reports.
Extracts quality.insert_size.peak for each strain,
groups by BioProject and sub-lineage.
"""
import json, os, csv, sys
from collections import defaultdict
from pathlib import Path
from statistics import median, mean, stdev

REPORTS_DIR = Path("/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/reports_v2")
STRAINS_CSV = Path("/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv")
FIG_DIR = Path("/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures")

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

# Load BioProject mapping from strains.csv
strain_bioproject = {}
with open(STRAINS_CSV) as f:
    for row in csv.DictReader(f):
        strain_bioproject[row["Id"]] = row.get("BioProject", "Unknown")

# Extract insert size + QC from reports
data = []
files = sorted(REPORTS_DIR.glob("*.json"))
print(f"Processing {len(files)} reports...")

for i, f in enumerate(files):
    if (i+1) % 200 == 0:
        print(f"  {i+1}/{len(files)}...", file=sys.stderr)
    
    with open(f) as fh:
        report = json.load(fh)
    
    sid = report["strain_id"]
    sub = "L4.11.1" if sid in L4111 else "L4.11.2"
    bp = strain_bioproject.get(sid, "Unknown")
    
    quality = report.get("quality", {})
    insert_size = quality.get("insert_size", {})
    peak = insert_size.get("peak", 0)
    
    before = quality.get("before_filtering", {})
    after = quality.get("after_filtering", {})
    
    dup = quality.get("duplication", {})
    dup_rate = dup.get("rate", 0) if isinstance(dup, dict) else 0
    
    read_len = after.get("read1_mean_length", 0)
    total_bases = after.get("total_bases", 0)
    gc = after.get("gc_content", 0)
    q30 = after.get("q30_rate", 0)
    total_reads = after.get("total_reads", 0)
    
    mapping = report.get("mapping_stats", {})
    depth = mapping.get("mean_depth", 0)
    coverage = mapping.get("covered_bases_percent", 0)
    
    data.append({
        "sid": sid, "sub": sub, "bp": bp,
        "insert_peak": peak, "read_len": read_len,
        "total_bases_M": total_bases / 1e6,
        "gc": gc, "q30": q30, "dup_rate": dup_rate,
        "depth": depth, "coverage": coverage,
        "total_reads_M": total_reads / 1e6,
    })

print(f"Loaded {len(data)} strains")

# ── Group by BioProject ──
bp_groups = defaultdict(list)
for d in data:
    bp_groups[d["bp"]].append(d)

# Sort by size
bp_sorted = sorted(bp_groups.items(), key=lambda x: -len(x[1]))

print(f"\n{'='*120}")
print(f"INSERT SIZE & QC BY BIOPROJECT")
print(f"{'='*120}")
print(f"{'BioProject':<20} {'N':>4} {'Sub':>8} {'Insert(med)':>12} {'ReadLen':>8} {'Depth(med)':>11} {'GC(med)':>8} {'Q30(med)':>9} {'Dup(med)':>9} {'Bases(med)':>11}")
print("-" * 120)

for bp, strains in bp_sorted:
    n = len(strains)
    subs = set(d["sub"] for d in strains)
    sub_str = "/".join(sorted(subs))
    
    inserts = [d["insert_peak"] for d in strains if d["insert_peak"] > 0]
    readlens = [d["read_len"] for d in strains if d["read_len"] > 0]
    depths = [d["depth"] for d in strains if d["depth"] > 0]
    gcs = [d["gc"] for d in strains if d["gc"] > 0]
    q30s = [d["q30"] for d in strains if d["q30"] > 0]
    dups = [d["dup_rate"] for d in strains if d["dup_rate"] > 0]
    bases = [d["total_bases_M"] for d in strains if d["total_bases_M"] > 0]
    
    ins_med = f"{median(inserts):.0f}" if inserts else "N/A"
    rl_med = f"{median(readlens):.0f}" if readlens else "N/A"
    dep_med = f"{median(depths):.1f}x" if depths else "N/A"
    gc_med = f"{median(gcs)*100:.1f}%" if gcs else "N/A"
    q30_med = f"{median(q30s)*100:.1f}%" if q30s else "N/A"
    dup_med = f"{median(dups)*100:.1f}%" if dups else "N/A"
    bas_med = f"{median(bases):.0f}M" if bases else "N/A"
    
    print(f"{bp:<20} {n:>4} {sub_str:>8} {ins_med:>12} {rl_med:>8} {dep_med:>11} {gc_med:>8} {q30_med:>9} {dup_med:>9} {bas_med:>11}")

# ── Figure ──
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 1. Insert size distribution by sub-lineage
fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=200)

# Panel A: Insert size histogram by sub-lineage
ax = axes[0, 0]
l1_inserts = [d["insert_peak"] for d in data if d["sub"] == "L4.11.1" and d["insert_peak"] > 0]
l2_inserts = [d["insert_peak"] for d in data if d["sub"] == "L4.11.2" and d["insert_peak"] > 0]
ax.hist(l1_inserts, bins=40, alpha=0.7, color="#1976d2", label=f"L4.11.1 (n={len(l1_inserts)})", density=True)
ax.hist(l2_inserts, bins=40, alpha=0.7, color="#d32f2f", label=f"L4.11.2 (n={len(l2_inserts)})", density=True)
ax.set_xlabel("Insert size peak (bp)", fontsize=11)
ax.set_ylabel("Density", fontsize=11)
ax.set_title("A. Insert size by sub-lineage", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(axis='y', alpha=0.3)

# Panel B: Insert size by BioProject (top 15)
ax = axes[0, 1]
top_bps = [bp for bp, strains in bp_sorted[:15] if bp != "Unknown"]
bp_inserts = []
bp_labels = []
bp_colors = []
for bp in top_bps:
    ins = [d["insert_peak"] for d in bp_groups[bp] if d["insert_peak"] > 0]
    if ins:
        bp_inserts.append(ins)
        subs = set(d["sub"] for d in bp_groups[bp])
        color = "#1976d2" if subs == {"L4.11.1"} else "#d32f2f" if subs == {"L4.11.2"} else "#9c27b0"
        bp_labels.append(f"{bp}\n(n={len(ins)})")
        bp_colors.append(color)

parts = ax.violinplot(bp_inserts, showmedians=True, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(bp_colors[i])
    pc.set_alpha(0.7)
parts['cmedians'].set_color('black')
ax.set_xticks(range(1, len(bp_labels) + 1))
ax.set_xticklabels(bp_labels, fontsize=7, rotation=45, ha='right')
ax.set_ylabel("Insert size peak (bp)", fontsize=11)
ax.set_title("B. Insert size by BioProject", fontsize=13, fontweight="bold")
ax.grid(axis='y', alpha=0.3)

# Panel C: Read length vs depth
ax = axes[1, 0]
for sub, color, marker in [("L4.11.1", "#1976d2", "o"), ("L4.11.2", "#d32f2f", "s")]:
    subset = [d for d in data if d["sub"] == sub]
    ax.scatter([d["read_len"] for d in subset], [d["depth"] for d in subset],
               c=color, alpha=0.3, s=10, marker=marker, label=sub)
ax.set_xlabel("Read length (bp)", fontsize=11)
ax.set_ylabel("Mean depth (×)", fontsize=11)
ax.set_title("C. Read length vs sequencing depth", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# Panel D: GC content vs Q30
ax = axes[1, 1]
for sub, color, marker in [("L4.11.1", "#1976d2", "o"), ("L4.11.2", "#d32f2f", "s")]:
    subset = [d for d in data if d["sub"] == sub]
    ax.scatter([d["gc"]*100 for d in subset], [d["q30"]*100 for d in subset],
               c=color, alpha=0.3, s=10, marker=marker, label=sub)
ax.set_xlabel("GC content (%)", fontsize=11)
ax.set_ylabel("Q30 rate (%)", fontsize=11)
ax.set_title("D. GC content vs Q30 rate", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.axvline(x=65.6, color='grey', linestyle='--', alpha=0.5, label='H37Rv GC=65.6%')

plt.suptitle("Sequencing Quality Control — 822 L4.11 strains (TBannotator v2)",
             fontsize=14, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.95])
out = FIG_DIR / "insert_size_qc.png"
plt.savefig(out, bbox_inches="tight", facecolor="white")
plt.close()
print(f"\n✅ Figure: {out}")
