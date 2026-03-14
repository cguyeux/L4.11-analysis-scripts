#!/usr/bin/env python3
"""Analysis #8: Convergent non-synonymous mutations in virulence genes.

Scans report.json SNP data for:
- PE/PPE family mutations
- mce operon mutations  
- ESX/Type VII secretion system mutations
- Other virulence-associated genes
Compares frequency between L4.11.1 and L4.11.2.
"""

import json, os, re
from collections import defaultdict

BDD = "BDD/L4.11"
L4111_SPDI = '"NC_000962.3:1700209:G:A"'

# Virulence gene families to track
VIRULENCE_PATTERNS = {
    "PE/PPE": re.compile(r'\b(PE\d|PPE\d|PE_PGRS)', re.I),
    "mce": re.compile(r'\bmce\d|yrbE', re.I),
    "ESX/T7SS": re.compile(r'\becc[A-E]|esp[A-K]|esx[A-Z]|mycP|EsxN|esxN|PE35|PPE68', re.I),
    "DosR regulon": re.compile(r'\bdos[RST]|hsp[X]|acr|narK|fdxA|pfkB|Rv207[89]', re.I),
    "lipid metabolism": re.compile(r'\bfad[A-Z]\d|pks\d|mas$|mmp[SL]|kas[AB]|acs[A-Z]|fbi[ABC]', re.I),
    "cell wall": re.compile(r'\bpim[A-F]|mmp[SL]|emb[A-C]|glf|rml[A-D]|wbb[L]', re.I),
    "toxin-antitoxin": re.compile(r'\bvap[BC]|maz[EF]|rel[BE]|hig[AB]|par[DE]', re.I),
    "sigma factors": re.compile(r'\bsig[A-M]|rsk[A]|rsd[A]', re.I),
    "two-component": re.compile(r'\bsens[XR]|mpr[AB]|pho[PRS]|dev[RS]|nar[LS]|trcp[RS]|kdp[DE]', re.I),
}

# Collect mutations by gene family
mutations_by_family = defaultdict(lambda: defaultdict(lambda: {"L4.11.1": 0, "L4.11.2": 0, "gene": "", "effect": ""}))
strain_counts = {"L4.11.1": 0, "L4.11.2": 0}
all_missense = defaultdict(lambda: {"L4.11.1": 0, "L4.11.2": 0, "gene": "", "effect": ""})

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
        strain_counts[lin] += 1
        
        data = json.loads(raw)
        snps = data.get("snp", [])
        
        for snp in snps:
            anns = snp.get("annotations", [])
            if not anns:
                continue
            ann = anns[0]
            effect = ann.get("annotation", [""])[0] if ann.get("annotation") else ""
            gene = ann.get("gene_name", "")
            hgvs_p = ann.get("hgvs_p", "") or ""
            spdi = snp.get("spdi", "")
            
            # Only non-synonymous
            if "missense" not in effect and "stop_gained" not in effect and "frameshift" not in effect:
                continue
            
            # Check against virulence patterns
            for family, pattern in VIRULENCE_PATTERNS.items():
                if pattern.search(gene):
                    key = f"{gene}:{hgvs_p}" if hgvs_p else f"{gene}:{spdi}"
                    mutations_by_family[family][key][lin] += 1
                    mutations_by_family[family][key]["gene"] = gene
                    mutations_by_family[family][key]["effect"] = hgvs_p or effect
                    break
            
            # Track ALL recurrent missense (>5% in either lineage)
            key = f"{gene}:{hgvs_p}" if hgvs_p else f"{gene}:{spdi}"
            all_missense[key][lin] += 1
            all_missense[key]["gene"] = gene
            all_missense[key]["effect"] = hgvs_p or effect
            
    except Exception as e:
        pass

print(f"Processed: {total} strains (L1={strain_counts['L4.11.1']}, L2={strain_counts['L4.11.2']})\n")

# ── VIRULENCE FAMILIES ──
for family in VIRULENCE_PATTERNS:
    muts = mutations_by_family[family]
    if not muts:
        continue
    
    print("=" * 100)
    print(f"  {family} MUTATIONS")
    print("=" * 100)
    
    # Sort by total count
    sorted_muts = sorted(muts.items(), key=lambda x: -(x[1]["L4.11.1"] + x[1]["L4.11.2"]))
    
    n1 = strain_counts["L4.11.1"]
    n2 = strain_counts["L4.11.2"]
    
    print(f"\n  {'Mutation':40s} {'L4.11.1':>15s} {'L4.11.2':>15s} {'Specificity'}")
    print("  " + "-" * 90)
    
    for key, counts in sorted_muts[:30]:
        c1 = counts["L4.11.1"]
        c2 = counts["L4.11.2"]
        pct1 = c1 * 100 / n1 if n1 else 0
        pct2 = c2 * 100 / n2 if n2 else 0
        
        if c1 + c2 < 3:  # skip very rare
            continue
        
        spec = ""
        if pct1 > 80 and pct2 < 5:
            spec = "◄ L4.11.1 specific"
        elif pct2 > 80 and pct1 < 5:
            spec = "◄ L4.11.2 specific"
        elif pct1 > 50 and pct2 > 50:
            spec = "★ SHARED"
        elif pct1 > 20 or pct2 > 20:
            spec = "◆ recurrent"
            
        print(f"  {key:40s} {c1:3d}/{n1} ({pct1:5.1f}%) {c2:3d}/{n2} ({pct2:5.1f}%) {spec}")
    print()

# ── HIGHLY RECURRENT NON-SYNONYMOUS (>30% in either) ──
print("\n" + "=" * 100)
print("  HIGHLY RECURRENT NON-SYNONYMOUS MUTATIONS (>30% in either lineage)")
print("=" * 100)

n1 = strain_counts["L4.11.1"]
n2 = strain_counts["L4.11.2"]

recurrent = []
for key, counts in all_missense.items():
    c1 = counts["L4.11.1"]
    c2 = counts["L4.11.2"]
    pct1 = c1 * 100 / n1 if n1 else 0
    pct2 = c2 * 100 / n2 if n2 else 0
    if pct1 > 30 or pct2 > 30:
        recurrent.append((key, counts, pct1, pct2))

recurrent.sort(key=lambda x: -(x[2] + x[3]))
print(f"\n  {len(recurrent)} mutations found")
print(f"\n  {'Mutation':40s} {'L4.11.1':>15s} {'L4.11.2':>15s}")
print("  " + "-" * 75)
for key, counts, pct1, pct2 in recurrent:
    c1 = counts["L4.11.1"]
    c2 = counts["L4.11.2"]
    spec = ""
    if pct1 > 80 and pct2 < 5:
        spec = "◄ L1"
    elif pct2 > 80 and pct1 < 5:
        spec = "◄ L2"
    elif pct1 > 50 and pct2 > 50:
        spec = "★ SHARED"
    print(f"  {key:40s} {c1:3d}/{n1} ({pct1:5.1f}%) {c2:3d}/{n2} ({pct2:5.1f}%) {spec}")
