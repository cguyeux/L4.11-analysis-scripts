#!/usr/bin/env python3
"""Transmission clustering: pairwise SNP distances and cluster identification.

Approach: 
- Load SPDI sets for all 828 strains
- Compute pairwise symmetric differences (= SNP distance)
- Identify clusters at <12 SNP threshold using single-linkage
"""

import os, sys
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = "NC_000962.3:1700209:G:A"
THRESHOLD = 12  # SNP distance threshold for recent transmission

# ── 1. LOAD ALL SPDI SETS ──
print("Loading SPDI sets...", file=sys.stderr)
strains = {}
lineage = {}

for sid in sorted(os.listdir(BDD)):
    spath = os.path.join(BDD, sid, "NC_000962.3", "spdi.txt")
    if not os.path.exists(spath):
        continue
    with open(spath) as f:
        spdis = set(line.strip() for line in f if line.strip())
    strains[sid] = spdis
    lineage[sid] = "L4.11.1" if L4111_SPDI in spdis else "L4.11.2"

n = len(strains)
print(f"Loaded {n} strains", file=sys.stderr)
sids = sorted(strains.keys())

# ── 2. PAIRWISE DISTANCES ──
print("Computing pairwise distances...", file=sys.stderr)
close_pairs = []  # (sid1, sid2, dist)
dist_histogram = defaultdict(int)

total_pairs = n * (n - 1) // 2
done = 0

for i in range(n):
    s1 = strains[sids[i]]
    for j in range(i + 1, n):
        s2 = strains[sids[j]]
        # Symmetric difference = sites different between the two
        dist = len(s1.symmetric_difference(s2))
        
        # Bin the distance
        if dist < 5:
            dist_histogram["0-4"] += 1
        elif dist < 12:
            dist_histogram["5-11"] += 1
        elif dist < 25:
            dist_histogram["12-24"] += 1
        elif dist < 50:
            dist_histogram["25-49"] += 1
        elif dist < 100:
            dist_histogram["50-99"] += 1
        else:
            dist_histogram["100+"] += 1
        
        if dist < THRESHOLD:
            close_pairs.append((sids[i], sids[j], dist))
        
        done += 1
        if done % 50000 == 0:
            print(f"  {done}/{total_pairs} ({done*100/total_pairs:.1f}%)", file=sys.stderr)

# ── 3. DISTANCE DISTRIBUTION ──
print(f"\n{'=' * 100}")
print("1. PAIRWISE SNP DISTANCE DISTRIBUTION")
print(f"{'=' * 100}")
print(f"\n  Total pairs: {total_pairs:,}")
for label in ["0-4", "5-11", "12-24", "25-49", "50-99", "100+"]:
    cnt = dist_histogram.get(label, 0)
    pct = cnt * 100 / total_pairs
    bar = "█" * int(pct * 2)
    print(f"  {label:>8s}: {cnt:8,} ({pct:6.2f}%) {bar}")

# ── 4. CLUSTERS AT <12 SNPs ──
print(f"\n{'=' * 100}")
print(f"2. CLOSE PAIRS (<{THRESHOLD} SNPs)")
print(f"{'=' * 100}")
print(f"\n  {len(close_pairs)} pairs below threshold")

# Single-linkage clustering
clusters = {}  # sid -> cluster_id
cluster_id = 0
for s1, s2, d in close_pairs:
    c1 = clusters.get(s1)
    c2 = clusters.get(s2)
    if c1 is None and c2 is None:
        clusters[s1] = cluster_id
        clusters[s2] = cluster_id
        cluster_id += 1
    elif c1 is not None and c2 is None:
        clusters[s2] = c1
    elif c2 is not None and c1 is None:
        clusters[s1] = c2
    elif c1 != c2:
        # Merge clusters
        old = c2
        for s in clusters:
            if clusters[s] == old:
                clusters[s] = c1

# Group by cluster
cluster_groups = defaultdict(list)
for sid, cid in clusters.items():
    cluster_groups[cid].append(sid)

# Remove singletons (shouldn't be any with our approach)
cluster_groups = {k: v for k, v in cluster_groups.items() if len(v) >= 2}

print(f"\n{'=' * 100}")
print(f"3. TRANSMISSION CLUSTERS (single-linkage, <{THRESHOLD} SNPs)")
print(f"{'=' * 100}")
print(f"\n  {len(cluster_groups)} clusters encompassing {len(clusters)} strains")

# Strains involved in clusters
clustered_strains = set(clusters.keys())
l1_clustered = sum(1 for s in clustered_strains if lineage[s] == "L4.11.1")
l2_clustered = sum(1 for s in clustered_strains if lineage[s] == "L4.11.2")
l1_total = sum(1 for s in sids if lineage[s] == "L4.11.1")
l2_total = sum(1 for s in sids if lineage[s] == "L4.11.2")

print(f"\n  Clustering rate:")
print(f"    L4.11.1: {l1_clustered}/{l1_total} ({l1_clustered*100/l1_total:.1f}%)")
print(f"    L4.11.2: {l2_clustered}/{l2_total} ({l2_clustered*100/l2_total:.1f}%)")
print(f"    Total:   {len(clustered_strains)}/{n} ({len(clustered_strains)*100/n:.1f}%)")

# Detail each cluster
print(f"\n  {'Cluster':>8s} {'Size':>5s} {'L4.11.1':>8s} {'L4.11.2':>8s} {'Max dist':>8s}")
print("  " + "-" * 45)

for cid in sorted(cluster_groups, key=lambda x: -len(cluster_groups[x])):
    members = sorted(cluster_groups[cid])
    l1 = sum(1 for s in members if lineage[s] == "L4.11.1")
    l2 = len(members) - l1
    # Max distance within cluster
    max_d = 0
    for p1, p2, d in close_pairs:
        if p1 in members and p2 in members:
            max_d = max(max_d, d)
    print(f"  C{cid:>7d} {len(members):5d} {l1:8d} {l2:8d} {max_d:8d}")

# ── 5. LARGEST CLUSTERS DETAIL ──
print(f"\n{'=' * 100}")
print("4. LARGEST CLUSTERS (top 10)")
print(f"{'=' * 100}")

sorted_clusters = sorted(cluster_groups.items(), key=lambda x: -len(x[1]))
for cid, members in sorted_clusters[:10]:
    members = sorted(members)
    l1 = sum(1 for s in members if lineage[s] == "L4.11.1")
    l2 = len(members) - l1
    print(f"\n  Cluster C{cid} (n={len(members)}, L1={l1}, L2={l2}):")
    for s in members[:20]:
        print(f"    {s} ({lineage[s]})")
    if len(members) > 20:
        print(f"    ... and {len(members)-20} more")

# ── 6. CROSS-LINEAGE PAIRS ──
print(f"\n{'=' * 100}")
print("5. CROSS-LINEAGE CLOSE PAIRS")
print(f"{'=' * 100}")
cross = [(s1, s2, d) for s1, s2, d in close_pairs 
         if lineage[s1] != lineage[s2]]
intra_l1 = [(s1, s2, d) for s1, s2, d in close_pairs 
            if lineage[s1] == "L4.11.1" and lineage[s2] == "L4.11.1"]
intra_l2 = [(s1, s2, d) for s1, s2, d in close_pairs 
            if lineage[s1] == "L4.11.2" and lineage[s2] == "L4.11.2"]

print(f"\n  Intra-L4.11.1 pairs: {len(intra_l1)}")
print(f"  Intra-L4.11.2 pairs: {len(intra_l2)}")
print(f"  Cross-lineage pairs: {len(cross)}")
if cross:
    print("\n  Cross-lineage pair details:")
    for s1, s2, d in sorted(cross, key=lambda x: x[2]):
        print(f"    {s1} ({lineage[s1]}) <-> {s2} ({lineage[s2]}) = {d} SNPs")
