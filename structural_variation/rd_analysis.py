#!/usr/bin/env python3
"""RD-typing validation: extract large_rd and missing_rd profiles for L4.11."""

import json, os
from collections import Counter, defaultdict

BDD = "BDD/L4.11"

L4111 = set("ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448 ERR4553470 ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568 ERR4553613 ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721 ERR4553747 ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834 ERR4553841 ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942 ERR4553947 ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739 SRR21661641 SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865 SRR29016881 SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138 SRR29017144 SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234 SRR29017242 SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583 SRR29017611 SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674 SRR29017683 SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528 SRR29440700 SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600 SRR35281602 SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186 SRR6797722 SRR6797801".split())

# Per-strain data
total_l1 = total_l2 = 0
l1_large_rd = defaultdict(lambda: Counter())  # rd_name -> {present/absent/partial: count}
l2_large_rd = defaultdict(lambda: Counter())
l1_missing_rd = Counter()
l2_missing_rd = Counter()

# Fingerprint profiles
l1_profiles = Counter()
l2_profiles = Counter()

# First pass: explore schema from one sample
sample_found = False

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
    
    # Print schema from first strain
    if not sample_found:
        sample_found = True
        print("=== SCHEMA EXPLORATION (first strain) ===")
        for rd in data.get("large_rd", [])[:3]:
            print(f"  large_rd entry: {rd}")
        for rd in data.get("missing_rd", [])[:3]:
            print(f"  missing_rd entry: {rd}")
        print()
    
    # Extract large_rd
    large_rd_dict = {}
    for rd in data.get("large_rd", []):
        name = rd if isinstance(rd, str) else rd.get("name", rd.get("rd_name", str(rd)))
        status = "present"
        if isinstance(rd, dict):
            status = rd.get("status", rd.get("type", "present"))
            pct = rd.get("percent_coverage", rd.get("pct", 100))
            if pct < 50:
                status = "deleted"
        large_rd_dict[name] = status
        if sl == "L4.11.1":
            l1_large_rd[name][status] += 1
        else:
            l2_large_rd[name][status] += 1
    
    # Extract missing_rd
    missing_rds = set()
    for rd in data.get("missing_rd", []):
        name = rd if isinstance(rd, str) else rd.get("name", rd.get("rd_name", str(rd)))
        missing_rds.add(name)
        if sl == "L4.11.1":
            l1_missing_rd[name] += 1
        else:
            l2_missing_rd[name] += 1
    
    # Create profile signature
    profile = tuple(sorted(missing_rds))
    if sl == "L4.11.1":
        l1_profiles[profile] += 1
    else:
        l2_profiles[profile] += 1

print(f"Processed: L4.11.1={total_l1}, L4.11.2={total_l2}")

# ── 1. LARGE RD SUMMARY ──
print("\n" + "=" * 100)
print("1. LARGE RD STATUS")
print("=" * 100)
all_rd = sorted(set(l1_large_rd.keys()) | set(l2_large_rd.keys()))
for rd in all_rd:
    l1_statuses = dict(l1_large_rd[rd])
    l2_statuses = dict(l2_large_rd[rd])
    l1_str = ", ".join(f"{s}={n}" for s, n in sorted(l1_statuses.items()))
    l2_str = ", ".join(f"{s}={n}" for s, n in sorted(l2_statuses.items()))
    print(f"  {rd:25s}  L1: {l1_str:30s}  L2: {l2_str}")

# ── 2. MISSING RD ──
print("\n" + "=" * 100)
print("2. MISSING RD (Regions of Difference with insufficient coverage)")
print("=" * 100)
all_missing = sorted(set(l1_missing_rd.keys()) | set(l2_missing_rd.keys()),
                      key=lambda x: -(l1_missing_rd.get(x, 0) + l2_missing_rd.get(x, 0)))
print(f"{'RD':25s} {'L1':>8s} {'L2':>8s}  Specificity")
print("-" * 80)
for rd in all_missing:
    n1, n2 = l1_missing_rd.get(rd, 0), l2_missing_rd.get(rd, 0)
    flag = ""
    if n1 > total_l1 * 0.8 and n2 > total_l2 * 0.8:
        flag = "★ SHARED (lineage-fixed)"
    elif n1 > total_l1 * 0.8 and n2 < total_l2 * 0.1:
        flag = "◄ L4.11.1 specific"
    elif n2 > total_l2 * 0.8 and n1 < total_l1 * 0.1:
        flag = "◄ L4.11.2 specific"
    elif n1 > total_l1 * 0.5 and n2 < total_l2 * 0.1:
        flag = "◄ L4.11.1 major"
    elif n2 > total_l2 * 0.5 and n1 < total_l1 * 0.1:
        flag = "◄ L4.11.2 major"
    print(f"{rd:25s} {n1:4d}/{total_l1:3d}  {n2:4d}/{total_l2:3d}  {flag}")

# ── 3. RD PROFILES ──
print("\n" + "=" * 100)
print("3. MISSING RD PROFILES (unique combinations)")
print("=" * 100)
print(f"\nL4.11.1: {len(l1_profiles)} unique RD profiles from {total_l1} strains")
for prof, n in l1_profiles.most_common(5):
    print(f"  n={n:3d}: {prof}")

print(f"\nL4.11.2: {len(l2_profiles)} unique RD profiles from {total_l2} strains")
for prof, n in l2_profiles.most_common(5):
    print(f"  n={n:3d}: {prof}")

# ── 4. CONSENSUS PROFILE ──
print("\n" + "=" * 100)
print("4. CONSENSUS RD PROFILE")
print("=" * 100)
# RDs missing in >90% of ALL strains
total = total_l1 + total_l2
print("\n  RDs absent in >90% of ALL L4.11 strains:")
for rd in all_missing:
    n = l1_missing_rd.get(rd, 0) + l2_missing_rd.get(rd, 0)
    if n > total * 0.9:
        print(f"    {rd:25s}: {n}/{total} ({n*100//total}%)")

print("\n  RDs absent in >90% of L4.11.1 only:")
for rd in all_missing:
    n1 = l1_missing_rd.get(rd, 0)
    n2 = l2_missing_rd.get(rd, 0)
    if n1 > total_l1 * 0.9 and n2 < total_l2 * 0.1:
        print(f"    {rd:25s}: L1={n1}/{total_l1}, L2={n2}/{total_l2}")

print("\n  RDs absent in >90% of L4.11.2 only:")
for rd in all_missing:
    n1 = l1_missing_rd.get(rd, 0)
    n2 = l2_missing_rd.get(rd, 0)
    if n2 > total_l2 * 0.9 and n1 < total_l1 * 0.1:
        print(f"    {rd:25s}: L1={n1}/{total_l1}, L2={n2}/{total_l2}")
