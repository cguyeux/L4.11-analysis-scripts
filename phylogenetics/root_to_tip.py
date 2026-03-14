#!/usr/bin/env python3
"""Root-to-tip regression using SRA RunInfo for dates.
Uses smaller Entrez batches and esummary instead of efetch.
"""

import os, sys, time
from Bio import Phylo, Entrez
from io import StringIO
import numpy as np
from scipy import stats
import xml.etree.ElementTree as ET

Entrez.email = "christophe.guyeux@univ-fcomte.fr"

# ── 1. LOAD TREE ──
print("Loading tree...", file=sys.stderr)
with open("articles/L4.11/résultats/T3.raxml.bestTree.cleaned") as f:
    tree_str = f.read().strip()
tree = Phylo.read(StringIO(tree_str), "newick")
tree.root_at_midpoint()

leaf_distances = {}
for clade in tree.get_terminals():
    leaf_distances[clade.name] = tree.distance(clade)
print(f"Loaded {len(leaf_distances)} leaves", file=sys.stderr)

# ── 2. GET DATES VIA ESUMMARY (lighter than efetch) ──
sra_accs = [a for a in leaf_distances if a.startswith(('SRR', 'ERR', 'DRR'))]
dates = {}
batch_size = 50  # smaller batches

for i in range(0, len(sra_accs), batch_size):
    batch = sra_accs[i:i+batch_size]
    query = " OR ".join(f"{acc}[ACCN]" for acc in batch)
    
    try:
        handle = Entrez.esearch(db="sra", term=query, retmax=batch_size+10)
        result = Entrez.read(handle)
        handle.close()
        
        ids = result.get("IdList", [])
        if not ids:
            continue
        
        # Use esummary - lighter
        handle = Entrez.esummary(db="sra", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        
        for summary in summaries:
            # Parse the ExpXml field for run accession and date
            exp_xml = summary.get("ExpXml", "")
            runs_xml = summary.get("Runs", "")
            
            # Extract run accession
            for acc in batch:
                if acc in runs_xml:
                    # Extract collection date from ExpXml  
                    # Look for date patterns in biosample attributes
                    import re
                    # Try collection_date from summary
                    date_match = re.search(r'collection.date.*?(\d{4})', exp_xml)
                    if not date_match:
                        date_match = re.search(r'(\d{4})', exp_xml)
                    if date_match:
                        year = int(date_match.group(1))
                        if 1990 <= year <= 2026:
                            dates[acc] = float(year)
        
        print(f"  Batch {i//batch_size + 1}/{len(sra_accs)//batch_size + 1}: "
              f"{len(dates)} dates total", file=sys.stderr)
        time.sleep(0.4)  # Rate limit
        
    except Exception as e:
        print(f"  Error batch {i//batch_size + 1}: {e}", file=sys.stderr)
        time.sleep(1)

print(f"\nTotal dates extracted: {len(dates)}", file=sys.stderr)

# Save cache
dates_file = "articles/L4.11/résultats/dates.tsv"
with open(dates_file, 'w') as f:
    for acc, yr in sorted(dates.items()):
        f.write(f"{acc}\t{yr}\n")

# ── 3. REGRESSION ──
paired = [(dates[a], leaf_distances[a], a) for a in dates if a in leaf_distances]
print(f"{len(paired)} paired points", file=sys.stderr)

if len(paired) < 10:
    print(f"Only {len(paired)} paired points - insufficient for regression")
    sys.exit(1)

years = np.array([p[0] for p in paired])
dists = np.array([p[1] for p in paired])

slope, intercept, r_value, p_value, std_err = stats.linregress(years, dists)
r_squared = r_value ** 2
tmrca = -intercept / slope if slope != 0 else float('nan')

print(f"\n{'=' * 80}")
print("ROOT-TO-TIP REGRESSION")
print(f"{'=' * 80}")
print(f"\n  n = {len(paired)}")
print(f"  Date range: {min(years):.0f} – {max(years):.0f}")
print(f"  R² = {r_squared:.4f}")
print(f"  p = {p_value:.2e}")
print(f"  Clock rate: {slope:.2e} subst/site/year")
print(f"  TMRCA: {tmrca:.1f}")

# Sub-lineage
L4111_SPDI = "NC_000962.3:1700209:G:A"
l1 = set()
for sid in os.listdir("BDD/L4.11"):
    sp = os.path.join("BDD/L4.11", sid, "NC_000962.3", "spdi.txt")
    if os.path.exists(sp):
        with open(sp) as f:
            if L4111_SPDI in f.read():
                l1.add(sid)

l1_p = [(y, d) for y, d, a in paired if a in l1]
l2_p = [(y, d) for y, d, a in paired if a not in l1]

for name, sub in [("L4.11.1", l1_p), ("L4.11.2", l2_p)]:
    if len(sub) >= 3:
        sy, sd = zip(*sub)
        sl, si, rv, pv, se = stats.linregress(sy, sd)
        print(f"\n  {name}: n={len(sub)}, R²={rv**2:.4f}, p={pv:.2e}, "
              f"rate={sl:.2e}, TMRCA={-si/sl:.1f}" if sl else "")
