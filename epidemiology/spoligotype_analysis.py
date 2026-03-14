#!/usr/bin/env python3
"""
Spoligotype pattern analysis for L4.11 strains.
1. Extract and classify Spol43 profiles
2. Match to SIT numbers via SIT.csv
3. Compare spoligotype clusters with SNP sub-lineage
4. Generate publication-quality visualization
"""

import csv, re, pickle
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

CSV = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv"
SIT_CSV = "/tmp/SIT.csv"
SIT2LINEAGE = "/home/christophe/Documents/codes/MTBC/TBannotator/data2/sit_to_lineage.pkl"
SIT2SRA = "/home/christophe/Documents/codes/MTBC/TBannotator/data2/sit_to_sra.pkl"
OUT_PNG = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/spoligotype_analysis.png"
OUT_TXT = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/spoligotype_analysis.txt"

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

PROTO = set("SRR3675589 ERR13289532 SRR35281598 SRR35281599 SRR35281596 SRR35281602 SRR35281600".split())

def spol_to_binary(spol_str):
    """Convert ■□ string to binary string (1/0)."""
    return spol_str.replace("■", "1").replace("□", "0")

def binary_to_octal(binary_43):
    """Convert 43-bit spoligotype binary to octal."""
    # Pad to 45 bits (15 octal digits) by adding 2 zeros at the end
    padded = binary_43 + "00"
    octals = []
    for i in range(0, 45, 3):
        octals.append(str(int(padded[i:i+3], 2)))
    return "".join(octals)

def assign_subgroup(acc):
    if acc in PROTO:
        return "Proto-L4.11.1"
    elif acc in L4111:
        return "L4.11.1"
    else:
        return "L4.11.2"

# ── 1. Load strains and extract spoligotypes ─────────────────────────────────
print("Loading strains...")
strains = []
with open(CSV) as f:
    for row in csv.DictReader(f):
        acc = row["Id"]
        spol43 = row.get("Spol43", "").strip()
        spol98 = row.get("Spol98", "").strip()
        country = row.get("Country", "").strip() or "Unknown"
        sub = assign_subgroup(acc)
        if spol43:
            binary = spol_to_binary(spol43)
            octal = binary_to_octal(binary)
            strains.append({
                "id": acc, "sub": sub, "country": country,
                "spol43": spol43, "binary43": binary, "octal": octal,
            })

print(f"  {len(strains)} strains with Spol43 data")

# ── 2. Build SIT lookup from SIT.csv ─────────────────────────────────────────
print("Loading SIT database...")
# SIT.csv: IsoNumber,Nb Strains,Spoligotype Binary,Spoligotype Octal,...,SIT,...,Clade,...
sit_db = {}  # binary → (SIT, clade)
sit_octal = {}  # octal → (SIT, clade)

# In SIT.csv, 'o' = spacer present (■/1), 'n' = spacer absent (□/0)
# Wait - let me verify the convention
# From earlier output: SIT 1 has "oooooooooooooooooooooooooooooooooonnnnnnnnn" = Beijing
# Beijing has spacers 1-34 absent, 35-43 present → so 'o' = absent, 'n' = present?
# No, Beijing has spacers 1-34 ABSENT → 'o' = absent (0), 'n' = present (1)?
# Actually in standard spoligotyping: Beijing shows NO hybridization for spacers 1-34
# Let me check: the octal for SIT 1 (Beijing) is 000000000003771
# 000000000003771 in binary (3 bits each): 000 000 000 000 000 000 000 000 000 011 111 111 001
# That means spacers 1-27 = 0, 28-34 = 1, 35-39 = 1, 40-43 = 01
# Wait, octal 000000000003771:
# Last digit '1' → 001 (but it's only 1 bit for last = spacer 43 present)
# Standard: 3771 for Beijing = spacers 35-43: 011 111 111 001
# Let's just match on octal

with open(SIT_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        sit = row.get("SIT", "").strip()
        clade = row.get("Clade", "").strip()
        octal_val = row.get("Spoligotype Octal", "").strip()
        binary_val = row.get("Spoligotype Binary", "").strip()
        
        if sit and octal_val:
            # Convert binary: o=0 (absent), n=1 (present)
            if binary_val:
                bin_conv = binary_val.replace("n", "1").replace("o", "0")
                key = (sit, clade)
                if octal_val not in sit_octal:
                    sit_octal[octal_val] = (sit, clade)

# Also build from binary matching
print(f"  {len(sit_octal)} unique octal patterns in SIT database")

# ── 3. Match strains to SIT ──────────────────────────────────────────────────
print("\nMatching L4.11 strains to SIT numbers...")

# Our Spol43: ■=1 (present), □=0 (absent)
# SIT binary: n=present, o=absent
# Our octal uses same convention (1=present) 
# But SIT octal uses: o=absent(0), n=present(1)? Let me reverse check.
# 
# SIT 358 from earlier: binary = nnnoonnnnnnnnnnnnnnnnnnnnnnnnnnnoooonnnnnnn
# Converting n→1, o→0: 111001111111111111111111111111110000111111111
# Wait that's the SIT convention for spacer notation.
# But our strains use ■=present(1), □=absent(0)
# SRR35281599 matched SIT358: ■■■□□■■■■■■■■■■■■■■■■■■■■■■■■■■■□□□□■■■■■■■
# Binary: 11100111111111111111111111111111000011111 (wait, only 43 chars)
# SIT 358 binary: nnnoonnnnnnnnnnnnnnnnnnnnnnnnnnnoooonnnnnnn (43 chars)
# n→1,o→0:         11100111111111111111111111111111000011111111

# Wait, SIT convention is INVERTED: 'n' = NO hybridization = spacer ABSENT
# and 'o' = hybridization OBSERVED = spacer PRESENT
# So: n=0, o=1
# SIT 358: nnnoonnnnnnnnnnnnnnnnnnnnnnnnnnnoooonnnnnnn
# n=0, o=1: 0001000000000000000000000000000011110000000
# Octal: 717777777760771? That doesn't match...
#
# Actually let me check SIT 1 (Beijing): 
# Binary: oooooooooooooooooooooooooooooooooonnnnnnnnn (43 chars)
# If o=1,n=0: 1111111111111111111111111111111110000000000 → not Beijing
# If o=0,n=1: 0000000000000000000000000000000001111111111 → Beijing pattern!
# Octal 000000000003771 matches o=0,n=1 convention
# So SIT: o=0 (absent), n=1 (present) — SAME as our ■=1,□=0
# Wait no: o=absent=0, n=present=1
# But Beijing has NO spacers 1-34, so they should be 0 → o=0 correct!
# And spacers 35-43 present → n=1... wait that's wrong semantically
# n means "negative" hybridization = NO signal = spacer ABSENT = 0
# o means "positive/observed" hybridization = spacer PRESENT = 1
# Let me verify: Beijing pattern has spacers 35-43 PRESENT only
# SIT1: oooooooooooooooooooooooooooooooooo nnnnnnnnn  
#        1-34: no signal (absent)            35-43: signal (present)?
# No! Standard Beijing: only spacers 35-43 are PRESENT
# If o at positions 1-34 means absent, then o=absent? But 'o' usually means 'observed'
#
# The octal 000000000003771 for Beijing:
# Binary (3 bits each, 15 groups): 
# 000 000 000 000 000 000 000 000 000 003 771
# Wait, octal digits 0-7. 3=011, 7=111, 1=001
# So: 000 000 000 000 000 000 000 000 000 011 111 111 001
# =    0   0   0   0   0   0   0   0   0  0 1 1 1 1 1 1 1 1 0 0 1
# That's... 45 bits for 43 spacers. Last 2 bits padding.
# Spacers: 000000000000000000000000000011111111001
# No that's 39 bits... Let me just match by octal string.

# Let me just compute our octals correctly and match.
# Our binary: 1=present, 0=absent (same as SIT n convention... wait)
# Actually our ■=present=1 and SIT's n=... 
# I need to check empirically. 
# SRR35281599 Spol43: ■■■□□■■■■■■■■■■■■■■■■■■■■■■■■■■■□□□□■■■■■■■
# binary: 11100111111111111111111111111111000011111 (wait 43 chars?)
# Let me count: ■■■□□■■■■■■■■■■■■■■■■■■■■■■■■■■■□□□□■■■■■■■
# 1110011111111111111111111111111100001111111 = 43 chars? Let me count
# ■■■ □□ ■■■■■■■■■■■■■■■■■■■■■■■■■■ □□□□ ■■■■■■■
# 3 + 2 + 26 + 4 + 7 = 42... let me just count directly

# OK let me just use the octal from the script and match
# Our convention: ■=1 (present), □=0 (absent)
# SIT convention in SIT.csv: from the octal 717777777760771 for SIT 358
# If we convert that octal to binary:
# 7=111, 1=001, 7=111, etc.
# 111 001 111 111 111 111 111 111 111 110 000 111 111 001
# = 11100111111111111111111111111111000011111100 1
# That's 43 bits (14*3 + 1 = 43)
# Matches: 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 1
# Hmm that's 44... 
# Standard octal spoligotype: 15 octal digits, each = 3 bits → 45 bits
# Last 2 bits are padding (always 0)
# So real spoligotype = first 43 bits

# Our SRR35281599: ■■■□□■■■■■■■■■■■■■■■■■■■■■■■■■■■□□□□■■■■■■■
# Let me count character by character:
# ■(1) ■(2) ■(3) □(4) □(5) ■(6) ■(7) ■(8) ■(9) ■(10)
# ■(11) ■(12) ■(13) ■(14) ■(15) ■(16) ■(17) ■(18) ■(19) ■(20)
# ■(21) ■(22) ■(23) ■(24) ■(25) ■(26) ■(27) ■(28) ■(29) ■(30)
# ■(31) ■(32) □(33) □(34) □(35) □(36) ■(37) ■(38) ■(39) ■(40)
# ■(41) ■(42) ■(43) = 43 characters ✓
# Binary: 11100111111111111111111111111111100001111111
# Hmm 43 chars: 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1
# Octal (pad to 45): 11100111111111111111111111111111100001111111 + 00
#                   111 001 111 111 111 111 111 111 111 111 000 011 111 110 0
# = 7 1 7 7 7 7 7 7 7 7 0 3 7 6 0
# Hmm...  wait pad with 00 at end:
# 111 001 111 111 111 111 111 111 111 111 000 011 111 1100
# Groups of 3: 111 | 001 | 111 | 111 | 111 | 111 | 111 | 111 | 111 | 111 | 000 | 011 | 111 | 110 | 0
# Last group needs padding: 0 → 000
# So: 7 1 7 7 7 7 7 7 7 7 0 3 7 6 0
# = 717777777703760
# But SIT 358 octal is 717777777760771... these don't match.

# I think the issue is the o/n convention is OPPOSITE to ■/□
# Let me try: n=0 (negative/absent), o=1 (observed/present)
# SIT 358 binary: nnnoonnnnnnnnnnnnnnnnnnnnnnnnnnnoooonnnnnnn
# n=0,o=1:         00010000000000000000000000000000111100000000  (43 chars? count: nnn=3, oo=5, nnn...=31+, oooo=35-38, nnn...=39-43)
# Actually: n n n o o n n n n n n n n n n n n n n n n n n n n n n n n n n n o o o o n n n n n n n
# = 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0
# Octal: pad → 00011000000000000000000000000000011110000000 + 00
# Groups: 000 110 000 000 000 000 000 000 000 000 000 000 111 000 000
# = 0 6 0 0 0 0 0 0 0 0 0 0 7 0 0 = 060000000000700
# That also doesn't match 717777777760771

# I think the convention must be: n=1, o=0 (opposite of what semantics suggest)
# SIT 358: nnnoonnnnnnnnnnnnnnnnnnnnnnnnnnnoooonnnnnnn
# n=1,o=0:  1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1
# OUR SRR35281599: 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1
# MATCH! ✓ So n=1 (filled=present), o=0 (open/absent)!
# And our ■=1, □=0. So ■≡n and □≡o. Convention confirmed.
# Octal: 11100111111111111111111111111111100001111111 + 00  (pad to 45)
# 111 001 111 111 111 111 111 111 111 111 000 011 111 110 0→000
# Hmm let me recount carefully:
# pos: 1234567890123456789012345678901234567890123
# bin: 11100111111111111111111111111111100001111111
# That's 41 chars... Let me use actual code.

# OK, I'll do this properly in the script. The SIT.csv uses Spoligotype Octal directly.
# I need to compute our strains' octal the SAME way and match.

matched = 0
unmatched = 0

for s in strains:
    # Compute octal from our binary43
    b = s["binary43"]   # 43 chars of 0/1
    # Standard: pad to 45 with 2 trailing 0s, then group by 3
    padded = b + "00"
    octal_digits = []
    for i in range(0, 45, 3):
        octal_digits.append(str(int(padded[i:i+3], 2)))
    our_octal = "".join(octal_digits)
    s["octal"] = our_octal
    
    if our_octal in sit_octal:
        sit_num, clade = sit_octal[our_octal]
        s["SIT"] = sit_num
        s["clade"] = clade
        matched += 1
    else:
        s["SIT"] = "Orphan"
        s["clade"] = "Unknown"
        unmatched += 1

print(f"  Matched to known SIT: {matched}")
print(f"  Orphan patterns: {unmatched}")

# ── 4. Spoligotype profile summary ──────────────────────────────────────────
print(f"\n{'='*90}")
print(f"SPOLIGOTYPE PROFILES — L4.11 ({len(strains)} strains)")
print(f"{'='*90}")

# Distribution by SIT and sub-group
sit_counts = defaultdict(lambda: {"total": 0, "L4.11.1": 0, "L4.11.2": 0, "Proto-L4.11.1": 0, "clade": ""})
for s in strains:
    sit = s["SIT"]
    sit_counts[sit]["total"] += 1
    sit_counts[sit][s["sub"]] += 1
    sit_counts[sit]["clade"] = s["clade"]

print(f"\n{'SIT':<12} {'Clade':<15} {'Total':>6} {'L4.11.1':>8} {'Proto':>6} {'L4.11.2':>8} {'Sample pattern':>50}")
print("-" * 110)

# Sort by total count (descending)
for sit, cnts in sorted(sit_counts.items(), key=lambda x: -x[1]["total"]):
    if cnts["total"] >= 2:
        # Find one example pattern
        example = next(s["spol43"] for s in strains if s["SIT"] == sit)
        print(f"{sit:<12} {cnts['clade']:<15} {cnts['total']:>6} {cnts['L4.11.1']:>8} {cnts['Proto-L4.11.1']:>6} {cnts['L4.11.2']:>8}  {example[:43]}")

# Orphan count
orphan_pattern_counts = Counter(s["spol43"] for s in strains if s["SIT"] == "Orphan")
print(f"\nOrphan patterns (unique): {len(orphan_pattern_counts)}")
for pat, cnt in orphan_pattern_counts.most_common(10):
    examples = [s for s in strains if s["spol43"] == pat]
    subs = Counter(s["sub"] for s in examples)
    print(f"  {pat[:43]}  n={cnt}  {dict(subs)}")

# ── 5. Concordance spoligotype × SNP sub-lineage ────────────────────────────
print(f"\n{'='*90}")
print(f"CONCORDANCE: SIT × SNP Sub-lineage")
print(f"{'='*90}")

# Major SITs by sub-lineage
for sub in ["L4.11.1", "Proto-L4.11.1", "L4.11.2"]:
    sub_strains = [s for s in strains if s["sub"] == sub]
    sit_dist = Counter(s["SIT"] for s in sub_strains)
    n = len(sub_strains)
    print(f"\n  {sub} (n={n}):")
    for sit, cnt in sit_dist.most_common(10):
        clade = next((s["clade"] for s in sub_strains if s["SIT"] == sit), "?")
        print(f"    SIT {sit:<10} ({clade:<12}): {cnt:>5} ({cnt/n*100:5.1f}%)")

# ── 6. Spacer deletion analysis ─────────────────────────────────────────────
print(f"\n{'='*90}")
print(f"SPACER PRESENCE/ABSENCE — per sub-lineage")
print(f"{'='*90}")

for sub in ["L4.11.1", "Proto-L4.11.1", "L4.11.2", "All"]:
    sub_strains = [s for s in strains if sub == "All" or s["sub"] == sub]
    n = len(sub_strains)
    if n == 0:
        continue
    
    # Spacer presence frequency
    presence = np.zeros(43)
    for s in sub_strains:
        b = s["binary43"]
        for i in range(min(43, len(b))):
            if b[i] == "1":
                presence[i] += 1
    
    pct = presence / n * 100
    print(f"\n  {sub} (n={n}):")
    # Show spacers with <100% and >0% presence
    variable = [(i+1, pct[i]) for i in range(43) if 0 < pct[i] < 100]
    absent = [(i+1, pct[i]) for i in range(43) if pct[i] == 0]
    present = [(i+1, pct[i]) for i in range(43) if pct[i] == 100]
    print(f"    Always present ({len(present)} spacers): {[sp for sp, _ in present]}")
    print(f"    Always absent ({len(absent)} spacers): {[sp for sp, _ in absent]}")
    print(f"    Variable ({len(variable)} spacers):")
    for sp, p in variable:
        print(f"      Spacer {sp:>2}: {p:5.1f}%")

# ── 7. Publication figure ────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(14, 7), dpi=200)
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.2], hspace=0.35, wspace=0.3)

# Panel A: SIT distribution pie chart for each sub-lineage
for idx, (sub, color_main) in enumerate([("L4.11.1", "#2563eb"), ("L4.11.2", "#dc2626")]):
    ax = fig.add_subplot(gs[0, idx])
    sub_strains = [s for s in strains if s["sub"] == sub or (sub == "L4.11.1" and s["sub"] == "Proto-L4.11.1")]
    sit_dist = Counter(s["SIT"] for s in sub_strains)
    n = len(sub_strains)
    
    # Top SITs + Other
    top_sits = sit_dist.most_common(6)
    top_total = sum(c for _, c in top_sits)
    other = n - top_total
    
    labels = [f"SIT{s}" if s != "Orphan" else "Orphan" for s, _ in top_sits]
    sizes = [c for _, c in top_sits]
    if other > 0:
        labels.append("Other")
        sizes.append(other)
    
    colors_pie = plt.cm.Set3(np.linspace(0, 1, len(labels)))
    wedges, texts, autotexts = ax.pie(sizes, labels=None, autopct=lambda p: f'{p:.0f}%' if p > 3 else '',
                                       colors=colors_pie, startangle=90,
                                       textprops=dict(fontsize=7))
    ax.legend(wedges, [f"{l} ({s})" for l, s in zip(labels, sizes)],
              fontsize=7, loc="center left", bbox_to_anchor=(-0.35, 0.5))
    title = f"{sub}" if sub == "L4.11.2" else "L4.11.1 + Proto"
    ax.set_title(f"{title} (n={n})", fontsize=11, fontweight="bold")

# Panel B: Spacer presence heatmap
ax = fig.add_subplot(gs[1, :])

# Build matrix: rows = sub-lineage, cols = 43 spacers
sub_labels = ["L4.11.1", "Proto-L4.11.1", "L4.11.2"]
matrix = np.zeros((3, 43))
for row_idx, sub in enumerate(sub_labels):
    sub_strains = [s for s in strains if s["sub"] == sub]
    n = len(sub_strains)
    if n == 0:
        continue
    for s in sub_strains:
        b = s["binary43"]
        for i in range(min(43, len(b))):
            if b[i] == "1":
                matrix[row_idx, i] += 1
    matrix[row_idx] /= n
    matrix[row_idx] *= 100

im = ax.imshow(matrix, cmap="RdYlGn", vmin=0, vmax=100, aspect="auto")
ax.set_yticks(range(3))
ax.set_yticklabels(sub_labels, fontsize=10, fontweight="bold")
ax.set_xticks(range(43))
ax.set_xticklabels(range(1, 44), fontsize=6, rotation=90)
ax.set_xlabel("Spacer position", fontsize=10)
ax.set_title("Spacer presence frequency (%) by sub-lineage", fontsize=11, fontweight="bold")
plt.colorbar(im, ax=ax, shrink=0.6, label="Presence (%)")

# Annotate cells with < 100% but > 0%
for i in range(3):
    for j in range(43):
        val = matrix[i, j]
        if 0 < val < 100:
            ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=5, fontweight="bold",
                    color="black" if val > 50 else "white")

plt.suptitle("Spoligotype analysis of L4.11 sub-lineages", fontsize=13, fontweight="bold")
plt.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
print(f"\nFigure: {OUT_PNG}")

# Save text results
with open(OUT_TXT, "w") as f:
    f.write(f"Spoligotype analysis — L4.11\n")
    f.write(f"Strains with Spol43: {len(strains)}\n")
    f.write(f"Matched to SIT: {matched}, Orphan: {unmatched}\n\n")
    for sit, cnts in sorted(sit_counts.items(), key=lambda x: -x[1]["total"]):
        if cnts["total"] >= 2:
            f.write(f"SIT {sit:<10} ({cnts['clade']:<12}): total={cnts['total']}, L4.11.1={cnts['L4.11.1']}, Proto={cnts['Proto-L4.11.1']}, L4.11.2={cnts['L4.11.2']}\n")

import shutil
shutil.copy(__file__, "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/spoligotype_analysis.py")

print(f"Results: {OUT_TXT}")
print("\nDone!")
