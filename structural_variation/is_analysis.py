#!/usr/bin/env python3
"""IS6110 in silico fingerprinting: compare L4.11.1 vs L4.11.2 IS profiles."""

import json, os
from collections import Counter, defaultdict

BDD = "BDD/L4.11"

L4111 = set("ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448 ERR4553470 ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568 ERR4553613 ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721 ERR4553747 ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834 ERR4553841 ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942 ERR4553947 ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739 SRR21661641 SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865 SRR29016881 SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138 SRR29017144 SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234 SRR29017242 SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583 SRR29017611 SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674 SRR29017683 SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528 SRR29440700 SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600 SRR35281602 SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186 SRR6797722 SRR6797801".split())

# Collect IS data per strain
strain_is = {}  # strain -> {is_name: [positions]}
all_is_names = Counter()
total_l1 = total_l2 = 0

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
    
    is_data = defaultdict(list)
    for entry in data.get("insertion_sequences", []):
        name = entry.get("name", "")
        pos = entry.get("position", 0)
        is_ref = entry.get("is_reference", False)
        is_data[name].append((pos, is_ref))
        all_is_names[name] += 1
    
    strain_is[strain_id] = dict(is_data)

print(f"Processed: L4.11.1={total_l1}, L4.11.2={total_l2}, Total={total_l1+total_l2}")

# ── 1. IS INVENTORY ──
print("\n" + "=" * 80)
print("1. INVENTAIRE DES SÉQUENCES D'INSERTION")
print("=" * 80)
for name, count in all_is_names.most_common():
    # Count per sublineage
    l1_strains = sum(1 for s, d in strain_is.items() if name in d and s in L4111)
    l2_strains = sum(1 for s, d in strain_is.items() if name in d and s not in L4111)
    l1_avg = sum(len(d.get(name, [])) for s, d in strain_is.items() if s in L4111) / max(total_l1, 1)
    l2_avg = sum(len(d.get(name, [])) for s, d in strain_is.items() if s not in L4111) / max(total_l2, 1)
    print(f"  {name:15s}: L1 present={l1_strains:3d}/{total_l1} (avg={l1_avg:.1f} copies)  "
          f"L2 present={l2_strains:3d}/{total_l2} (avg={l2_avg:.1f} copies)")

# ── 2. IS6110 COPY NUMBER DISTRIBUTION ──
print("\n" + "=" * 80)
print("2. IS6110 COPY NUMBER DISTRIBUTION")
print("=" * 80)
for sl_name, sl_ids in [("L4.11.1", L4111), ("L4.11.2", None)]:
    copies = []
    for s, d in strain_is.items():
        if sl_ids is not None and s not in sl_ids:
            continue
        if sl_ids is None and s in L4111:
            continue
        n_is6110 = len(d.get("IS6110", []))
        copies.append(n_is6110)
    
    if copies:
        cn_dist = Counter(copies)
        print(f"\n{sl_name} (n={len(copies)}):")
        print(f"  Range: {min(copies)}-{max(copies)}, Mean={sum(copies)/len(copies):.1f}, Median={sorted(copies)[len(copies)//2]}")
        for cn in sorted(cn_dist):
            bar = "█" * cn_dist[cn]
            print(f"  {cn:3d} copies: {cn_dist[cn]:4d} {bar}")

# ── 3. IS6110 POSITION ANALYSIS ──
print("\n" + "=" * 80)
print("3. IS6110 POSITIONS BY SUBLINEAGE")
print("=" * 80)

# Collect all IS6110 positions
l1_positions = Counter()
l2_positions = Counter()
for s, d in strain_is.items():
    if "IS6110" not in d:
        continue
    for pos, is_ref in d["IS6110"]:
        if s in L4111:
            l1_positions[pos] += 1
        else:
            l2_positions[pos] += 1

# Core positions (present in >80% of strains)
print("\nCore IS6110 positions (>80% of strains in sublineage):")
print(f"\n  L4.11.1 core (>{int(total_l1*0.8)} strains):")
for pos, n in sorted(l1_positions.items(), key=lambda x: x[0]):
    if n > total_l1 * 0.8:
        also_l2 = l2_positions.get(pos, 0)
        shared = "SHARED" if also_l2 > total_l2 * 0.8 else f"L2={also_l2}"
        print(f"    pos={pos:>10d}  n={n:3d}/{total_l1}  ({shared})")

print(f"\n  L4.11.2 core (>{int(total_l2*0.8)} strains):")
for pos, n in sorted(l2_positions.items(), key=lambda x: x[0]):
    if n > total_l2 * 0.8:
        also_l1 = l1_positions.get(pos, 0)
        shared = "SHARED" if also_l1 > total_l1 * 0.8 else f"L1={also_l1}"
        print(f"    pos={pos:>10d}  n={n:3d}/{total_l2}  ({shared})")

# L4.11.1-specific positions
print("\nIS6110 positions specific to L4.11.1 (>50% L1, <5% L2):")
for pos, n in sorted(l1_positions.items(), key=lambda x: -x[1]):
    if n > total_l1 * 0.5 and l2_positions.get(pos, 0) < total_l2 * 0.05:
        print(f"    pos={pos:>10d}  L1={n:3d}/{total_l1} ({n*100//total_l1}%)  L2={l2_positions.get(pos,0):3d}/{total_l2}")

# L4.11.2-specific positions
print("\nIS6110 positions specific to L4.11.2 (>50% L2, <5% L1):")
for pos, n in sorted(l2_positions.items(), key=lambda x: -x[1]):
    if n > total_l2 * 0.5 and l1_positions.get(pos, 0) < total_l1 * 0.05:
        print(f"    pos={pos:>10d}  L2={n:3d}/{total_l2} ({n*100//total_l2}%)  L1={l1_positions.get(pos,0):3d}/{total_l1}")

# ── 4. IS6110 FINGERPRINT PATTERN ──
print("\n" + "=" * 80)
print("4. IS6110 FINGERPRINT PATTERNS (unique position combinations)")
print("=" * 80)
# Create a fingerprint per strain (sorted tuple of positions)
l1_fingerprints = Counter()
l2_fingerprints = Counter()
for s, d in strain_is.items():
    if "IS6110" not in d:
        fp = tuple()
    else:
        fp = tuple(sorted(pos for pos, _ in d["IS6110"]))
    if s in L4111:
        l1_fingerprints[fp] += 1
    else:
        l2_fingerprints[fp] += 1

print(f"\nL4.11.1: {len(l1_fingerprints)} unique fingerprints from {total_l1} strains")
for fp, n in l1_fingerprints.most_common(5):
    print(f"  n={n:3d}: {len(fp)} copies at {fp[:5]}{'...' if len(fp)>5 else ''}")

print(f"\nL4.11.2: {len(l2_fingerprints)} unique fingerprints from {total_l2} strains")
for fp, n in l2_fingerprints.most_common(5):
    print(f"  n={n:3d}: {len(fp)} copies at {fp[:5]}{'...' if len(fp)>5 else ''}")

# ── 5. OTHER IS ELEMENTS ──
print("\n" + "=" * 80)
print("5. IS1081 COPY NUMBER COMPARISON")
print("=" * 80)
for sl_name, sl_ids in [("L4.11.1", L4111), ("L4.11.2", None)]:
    copies = []
    for s, d in strain_is.items():
        if sl_ids is not None and s not in sl_ids:
            continue
        if sl_ids is None and s in L4111:
            continue
        n = len(d.get("IS1081", []))
        copies.append(n)
    if copies:
        print(f"  {sl_name}: Mean={sum(copies)/len(copies):.1f}, Range={min(copies)}-{max(copies)}")
