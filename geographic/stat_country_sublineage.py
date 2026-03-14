#!/usr/bin/env python3
"""
Formal statistical tests for country × sub-lineage associations in L4.11.
- Overall Chi² test on the contingency table
- Per-country Fisher exact tests + Odds Ratios with 95% CI
- Results exported as LaTeX table fragment
"""

import csv, re
import numpy as np
from collections import Counter
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact
import math

CSV = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv"
OUT = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/stat_country_sublineage.txt"

# ── L4.11.1 list (91 strains) ────────────────────────────────────────────────
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

# ── Load data ────────────────────────────────────────────────────────────────
strains = []
with open(CSV) as f:
    for row in csv.DictReader(f):
        acc = row["Id"]
        country = row.get("Country", "").strip()
        if not country:
            country = "Unknown"
        sub = "L4.11.1" if acc in L4111 else "L4.11.2"
        strains.append({"id": acc, "country": country, "sub": sub})

total = len(strains)
n_l4111 = sum(1 for s in strains if s["sub"] == "L4.11.1")
n_l4112 = total - n_l4111
print(f"Total strains: {total}")
print(f"  L4.11.1: {n_l4111}")
print(f"  L4.11.2: {n_l4112}")

# ── Build contingency table ──────────────────────────────────────────────────
country_counts = Counter(s["country"] for s in strains)
# Only keep countries with ≥3 strains for meaningful stats
countries_min = [c for c, n in country_counts.most_common() if n >= 3]

print(f"\nCountries with ≥3 strains: {len(countries_min)}")

# Full contingency table (for Chi²)
# Rows = countries, Cols = [L4.11.1, L4.11.2]
table_countries = []
table_data = []
for country in countries_min:
    n1 = sum(1 for s in strains if s["country"] == country and s["sub"] == "L4.11.1")
    n2 = sum(1 for s in strains if s["country"] == country and s["sub"] == "L4.11.2")
    table_countries.append(country)
    table_data.append([n1, n2])

contingency = np.array(table_data)
print(f"\nContingency table ({len(table_countries)} countries × 2 sub-lineages):")
print(f"{'Country':<45} {'L4.11.1':>8} {'L4.11.2':>8} {'Total':>8}")
print("-" * 75)
for i, c in enumerate(table_countries):
    t = contingency[i].sum()
    print(f"{c:<45} {contingency[i,0]:>8} {contingency[i,1]:>8} {t:>8}")
print("-" * 75)
print(f"{'TOTAL':<45} {contingency[:,0].sum():>8} {contingency[:,1].sum():>8} {contingency.sum():>8}")

# ── Overall Chi² test ────────────────────────────────────────────────────────
chi2, p_chi2, dof, expected = chi2_contingency(contingency)
print(f"\n══════════════════════════════════════════════")
print(f"OVERALL CHI-SQUARED TEST")
print(f"  χ² = {chi2:.2f}, df = {dof}, p = {p_chi2:.2e}")
print(f"  Cramér's V = {np.sqrt(chi2 / (contingency.sum() * (min(contingency.shape) - 1))):.4f}")
print(f"══════════════════════════════════════════════")

# ── Per-country Fisher exact tests + OR ──────────────────────────────────────
print(f"\n{'Country':<40} {'n(L4.11.1)':>10} {'n(L4.11.2)':>10} {'OR':>8} {'95% CI':>18} {'p(Fisher)':>12} {'Sig':>5}")
print("=" * 110)

results = []
for country in countries_min:
    # 2×2 table: Country_yes/no × L4.11.1/L4.11.2
    a = sum(1 for s in strains if s["country"] == country and s["sub"] == "L4.11.1")
    b = sum(1 for s in strains if s["country"] == country and s["sub"] == "L4.11.2")
    c = n_l4111 - a
    d = n_l4112 - b
    
    table_2x2 = np.array([[a, b], [c, d]])
    
    # Fisher exact test
    odds_ratio, p_fisher = fisher_exact(table_2x2)
    
    # Compute 95% CI for OR using log transform
    if a > 0 and b > 0 and c > 0 and d > 0:
        log_or = math.log(odds_ratio)
        se_log_or = math.sqrt(1/a + 1/b + 1/c + 1/d)
        ci_low = math.exp(log_or - 1.96 * se_log_or)
        ci_high = math.exp(log_or + 1.96 * se_log_or)
    elif a == 0 or d == 0:
        ci_low, ci_high = 0.0, float('inf')
    else:
        ci_low, ci_high = float('inf'), float('inf')
    
    # Significance stars
    if p_fisher < 0.001:
        sig = "***"
    elif p_fisher < 0.01:
        sig = "**"
    elif p_fisher < 0.05:
        sig = "*"
    else:
        sig = "ns"
    
    # Bonferroni correction
    p_bonf = min(p_fisher * len(countries_min), 1.0)
    if p_bonf < 0.001:
        sig_bonf = "***"
    elif p_bonf < 0.01:
        sig_bonf = "**"
    elif p_bonf < 0.05:
        sig_bonf = "*"
    else:
        sig_bonf = "ns"
    
    ci_str = f"[{ci_low:.2f}–{ci_high:.2f}]" if ci_high < 1e6 else "[0.00–∞]"
    print(f"{country:<40} {a:>10} {b:>10} {odds_ratio:>8.2f} {ci_str:>18} {p_fisher:>12.2e} {sig:>5}")
    
    results.append({
        "country": country, "n1": a, "n2": b, "total": a+b,
        "OR": odds_ratio, "ci_low": ci_low, "ci_high": ci_high,
        "p_fisher": p_fisher, "p_bonf": p_bonf, "sig": sig, "sig_bonf": sig_bonf
    })

# ── Summary of significant associations ─────────────────────────────────────
print(f"\n══════════════════════════════════════════════")
print(f"SIGNIFICANT ASSOCIATIONS (Bonferroni-corrected)")
print(f"══════════════════════════════════════════════")
for r in sorted(results, key=lambda x: x["p_fisher"]):
    if r["sig_bonf"] != "ns":
        direction = "enriched in L4.11.1" if r["OR"] > 1 else "enriched in L4.11.2"
        ci_str = f"[{r['ci_low']:.2f}–{r['ci_high']:.2f}]" if r['ci_high'] < 1e6 else "[0.00–∞]"
        print(f"  {r['country']:<35} OR={r['OR']:>7.2f} {ci_str:>18}  p={r['p_fisher']:.2e} (Bonf: {r['p_bonf']:.2e}) {r['sig_bonf']}  → {direction}")

# ── Save results ─────────────────────────────────────────────────────────────
with open(OUT, "w") as f:
    f.write(f"Statistical tests: country × sub-lineage association in L4.11\n")
    f.write(f"N = {total} strains (L4.11.1: {n_l4111}, L4.11.2: {n_l4112})\n\n")
    f.write(f"Overall Chi²: χ² = {chi2:.2f}, df = {dof}, p = {p_chi2:.2e}\n")
    f.write(f"Cramér's V = {np.sqrt(chi2 / (contingency.sum() * (min(contingency.shape) - 1))):.4f}\n\n")
    f.write(f"Per-country Fisher exact tests + OR:\n")
    for r in sorted(results, key=lambda x: x["p_fisher"]):
        ci_str = f"[{r['ci_low']:.2f}-{r['ci_high']:.2f}]" if r['ci_high'] < 1e6 else "[0.00-Inf]"
        f.write(f"  {r['country']:<35} n={r['total']:>4} ({r['n1']:>3}/{r['n2']:>3})  OR={r['OR']:>7.2f} {ci_str:>18}  p={r['p_fisher']:.2e} (Bonf: {r['p_bonf']:.2e}) {r['sig_bonf']}\n")

print(f"\nResults saved to: {OUT}")
print("\nDone!")
