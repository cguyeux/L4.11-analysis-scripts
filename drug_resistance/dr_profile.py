#!/usr/bin/env python3
"""Complete DR mutational profile: all drug resistance genes across L4.11."""

import json, os, csv
from collections import Counter, defaultdict

BDD = "BDD/L4.11"
STRAINS_CSV = "articles/L4.11/résultats/strains.csv"

L4111 = set("ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448 ERR4553470 ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568 ERR4553613 ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721 ERR4553747 ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834 ERR4553841 ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942 ERR4553947 ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739 SRR21661641 SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865 SRR29016881 SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138 SRR29017144 SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234 SRR29017242 SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583 SRR29017611 SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674 SRR29017683 SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528 SRR29440700 SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600 SRR35281602 SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186 SRR6797722 SRR6797801".split())

# Drug resistance genes organized by drug/category
DR_GENES = {
    # First-line
    "INH":  ["katG", "inhA", "fabG1", "ahpC", "kasA", "ndh", "nat"],
    "RIF":  ["rpoB", "rpoA", "rpoC"],  # includes compensatory
    "EMB":  ["embB", "embC", "embA", "embR", "ubiA"],
    "PZA":  ["pncA", "panD", "rpsA", "clpC1"],
    # Second-line
    "FQ":   ["gyrA", "gyrB"],
    "AG":   ["rrs", "rpsL", "eis", "gidB"],
    "BDQ":  ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"],
    "LZD":  ["rrl", "rplC"],
    "ETH":  ["ethA", "ethR", "mshA"],
    "DCS":  ["alr", "ddl", "cycA"],
}

# Flatten
ALL_DR_GENES = set()
GENE_TO_DRUG = {}
for drug, genes in DR_GENES.items():
    for g in genes:
        ALL_DR_GENES.add(g)
        GENE_TO_DRUG[g] = drug

# Load phenotypes
strain_info = {}
with open(STRAINS_CSV) as f:
    for r in csv.DictReader(f):
        strain_info[r["Id"]] = {
            "INH": r["Isoniazid"], "RIF": r["Rifampicin"],
            "EMB": r["Ethambutol"], "PZA": r["Pyrazinamide"],
            "sl": "L4.11.1" if r["Id"] in L4111 else "L4.11.2"
        }

# Parse all reports
mutations = defaultdict(lambda: {"L4.11.1": [], "L4.11.2": []})
strain_dr_muts = defaultdict(list)  # strain -> list of (drug, gene, mut)
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
    
    for snp in data.get("snp", []):
        for ann in snp.get("annotations", []):
            gname = ann.get("gene_name", "")
            if gname not in ALL_DR_GENES:
                continue
            impact = ann.get("impact", "")
            if impact == "LOW":  # skip synonymous for DR analysis
                continue
            hgvs_p = ann.get("hgvs_p") or ""
            hgvs_c = ann.get("hgvs_c") or ""
            drug = GENE_TO_DRUG.get(gname, "?")
            key = f"{drug}|{gname}|{hgvs_p or hgvs_c}|{impact}"
            mutations[key][sl].append(strain_id)
            strain_dr_muts[strain_id].append((drug, gname, hgvs_p or hgvs_c))

print(f"Processed: L4.11.1={total_l1}, L4.11.2={total_l2}")
print(f"Unique non-synonymous DR mutations: {len(mutations)}")

# ── 1. ALL DR MUTATIONS BY FREQUENCY ──
print("\n" + "=" * 100)
print("1. ALL DR GENE MUTATIONS (non-synonymous, freq >= 3)")
print("=" * 100)
print(f"{'Drug':5s} {'Gene':8s} {'Mutation':30s} {'Impact':10s} {'L1':>5s} {'L2':>5s} {'Total':>6s}")
print("-" * 100)

for key, by_sl in sorted(mutations.items(), key=lambda x: -(len(x[1]["L4.11.1"])+len(x[1]["L4.11.2"]))):
    n1, n2 = len(by_sl["L4.11.1"]), len(by_sl["L4.11.2"])
    total = n1 + n2
    if total < 3:
        continue
    parts = key.split("|")
    drug, gene, mut, impact = parts[0], parts[1], parts[2], parts[3]
    # Flag if lineage-specific
    flag = ""
    if n1 > 0 and n2 == 0:
        flag = " ◄L1"
    elif n2 > 0 and n1 == 0:
        flag = " ◄L2"
    elif n1 > 0.8 * total_l1 and n2 > 0.8 * total_l2:
        flag = " ★FIXED"
    print(f"{drug:5s} {gene:8s} {mut:30s} {impact:10s} {n1:5d} {n2:5d} {total:6d}{flag}")

# ── 2. COMPENSATORY MUTATIONS (rpoA, rpoC) ──
print("\n" + "=" * 100)
print("2. COMPENSATORY MUTATIONS (rpoA, rpoC) — ALL frequencies")
print("=" * 100)
for key, by_sl in sorted(mutations.items()):
    parts = key.split("|")
    if parts[1] in ("rpoA", "rpoC"):
        n1, n2 = len(by_sl["L4.11.1"]), len(by_sl["L4.11.2"])
        print(f"  {parts[1]:6s} {parts[2]:30s} {parts[3]:10s} L1={n1:3d} L2={n2:3d}")

# ── 3. rpoB MUTATIONS ──
print("\n" + "=" * 100)
print("3. rpoB MUTATIONS (all)")
print("=" * 100)
for key, by_sl in sorted(mutations.items(), key=lambda x: -(len(x[1]["L4.11.1"])+len(x[1]["L4.11.2"]))):
    parts = key.split("|")
    if parts[1] == "rpoB":
        n1, n2 = len(by_sl["L4.11.1"]), len(by_sl["L4.11.2"])
        print(f"  {parts[2]:30s} {parts[3]:10s} L1={n1:3d} L2={n2:3d} Total={n1+n2}")

# ── 4. DRUG-LEVEL SUMMARY ──
print("\n" + "=" * 100)
print("4. DRUG-LEVEL SUMMARY (strains with ≥1 non-syn mutation in corresponding gene)")
print("=" * 100)
for drug_cat, gene_list in DR_GENES.items():
    l1_hit = set()
    l2_hit = set()
    for key, by_sl in mutations.items():
        parts = key.split("|")
        if parts[0] == drug_cat:
            for s in by_sl["L4.11.1"]:
                l1_hit.add(s)
            for s in by_sl["L4.11.2"]:
                l2_hit.add(s)
    l1p = len(l1_hit)/total_l1*100 if total_l1 else 0
    l2p = len(l2_hit)/total_l2*100 if total_l2 else 0
    print(f"  {drug_cat:5s} ({', '.join(gene_list):40s}): L1={len(l1_hit):3d}/{total_l1} ({l1p:5.1f}%)  L2={len(l2_hit):3d}/{total_l2} ({l2p:5.1f}%)")

# ── 5. LINEAGE-SPECIFIC PATTERNS ──
print("\n" + "=" * 100)
print("5. LINEAGE-SPECIFIC DR MUTATIONS (present in one sublineage only)")
print("=" * 100)
print("\n  == L4.11.1 ONLY ==")
for key, by_sl in sorted(mutations.items(), key=lambda x: -len(x[1]["L4.11.1"])):
    if len(by_sl["L4.11.1"]) >= 3 and len(by_sl["L4.11.2"]) == 0:
        parts = key.split("|")
        n = len(by_sl["L4.11.1"])
        print(f"    {parts[0]:5s} {parts[1]:8s} {parts[2]:30s} n={n:3d}/{total_l1} ({n/total_l1*100:.0f}%)")

print("\n  == L4.11.2 ONLY ==")
for key, by_sl in sorted(mutations.items(), key=lambda x: -len(x[1]["L4.11.2"])):
    if len(by_sl["L4.11.2"]) >= 3 and len(by_sl["L4.11.1"]) == 0:
        parts = key.split("|")
        n = len(by_sl["L4.11.2"])
        print(f"    {parts[0]:5s} {parts[1]:8s} {parts[2]:30s} n={n:3d}/{total_l2} ({n/total_l2*100:.0f}%)")

# ── 6. katG vs inhA pathway ──
print("\n" + "=" * 100)
print("6. INH RESISTANCE PATHWAY: katG vs inhA promoter")
print("=" * 100)
for sl_name, sl_total in [("L4.11.1", total_l1), ("L4.11.2", total_l2)]:
    katG_strains = set()
    inhA_strains = set()
    for key, by_sl in mutations.items():
        parts = key.split("|")
        if parts[1] == "katG":
            for s in by_sl[sl_name]:
                katG_strains.add(s)
        elif parts[1] in ("inhA", "fabG1"):
            for s in by_sl[sl_name]:
                inhA_strains.add(s)
    both = katG_strains & inhA_strains
    only_k = katG_strains - inhA_strains
    only_i = inhA_strains - katG_strains
    print(f"  {sl_name}: katG-only={len(only_k)}, inhA-only={len(only_i)}, both={len(both)}, neither={sl_total-len(katG_strains|inhA_strains)}")
