#!/usr/bin/env python3
"""
Publication-quality Bayesian Skyline Plot from TreeTime output.
Shows effective population size (Ne) over time for L4.11 ONLY.
Uses the pruned tree (825 L4.11 strains, no outgroup).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SKYLINE = "/tmp/treetime_skyline_l411_v2/skyline.tsv"
OUT_PNG = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/skyline_plot.png"

# Parse skyline data
dates, ne, lo, hi = [], [], [], []
with open(SKYLINE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) >= 4:
            dates.append(float(parts[0]))
            ne.append(float(parts[1]))
            lo.append(float(parts[2]))
            hi.append(float(parts[3]))

dates = np.array(dates)
ne = np.array(ne)
lo = np.array(lo)
hi = np.array(hi)

print(f"Skyline data: {len(dates)} time points")
print(f"Date range: {dates[0]:.0f} to {dates[-1]:.0f}")
print(f"Ne range: {ne.min():.1e} to {ne.max():.1e}")

# ── Publication figure ───────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=200,
                                gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.3})

# Panel A: Full history (log scale)
ax1.fill_between(dates, lo, hi, alpha=0.2, color="#3b82f6", label="95% CI")
ax1.plot(dates, ne, color="#1d4ed8", linewidth=2.5, label="$N_e$ (effective pop. size)")
ax1.set_yscale("log")
ax1.set_xlabel("Year (CE/BCE)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Effective population size ($N_e$)", fontsize=11, fontweight="bold")
ax1.set_title("(a) BSP — L4.11 full history (825 strains)", fontsize=11, fontweight="bold")
ax1.grid(True, alpha=0.15, which="both")
ax1.legend(fontsize=9, loc="upper left")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Key epoch annotations
ax1.axvline(x=0, color="gray", linestyle="--", alpha=0.3, linewidth=0.8)
ax1.axvline(x=-10000, color="#059669", linestyle=":", alpha=0.4, linewidth=1)
ax1.text(-10000, ne.max() * 1.5, "Neolithic\ntransition", fontsize=7, 
         color="#059669", ha="center", fontweight="bold")

# Panel B: Recent expansion (last 25000 years = most interesting part)
recent_mask = dates > -25000
if recent_mask.sum() > 2:
    ax2.fill_between(dates[recent_mask], lo[recent_mask], hi[recent_mask], 
                     alpha=0.2, color="#3b82f6")
    ax2.plot(dates[recent_mask], ne[recent_mask], color="#1d4ed8", linewidth=2.5,
             marker="o", markersize=4)
    ax2.set_yscale("log")
    ax2.set_xlabel("Year (CE/BCE)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("$N_e$", fontsize=11, fontweight="bold")
    ax2.set_title("(b) Recent expansion (last 25,000 years)", fontsize=11, fontweight="bold")
    ax2.grid(True, alpha=0.15, which="both")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    # Annotate expansion phases
    ax2.axvline(x=-10000, color="#059669", linestyle=":", alpha=0.4)
    ax2.text(-10000, lo[recent_mask].min() * 2, "Neolithic", fontsize=7, 
             color="#059669", ha="center", fontweight="bold", rotation=90)
    
    # Peak annotation
    recent_ne = ne[recent_mask]
    recent_dates = dates[recent_mask]
    max_idx = np.argmax(recent_ne)
    peak_date = recent_dates[max_idx]
    peak_ne = recent_ne[max_idx]
    ax2.annotate(f"Peak $N_e$\n≈ {peak_ne:.1e}\n({peak_date:.0f} CE)", 
                 xy=(peak_date, peak_ne),
                 xytext=(peak_date - 6000, peak_ne * 0.2),
                 fontsize=8, fontweight="bold", color="#dc2626",
                 arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.5))
    
    # Mark 2015
    ax2.axvline(x=2015, color="#dc2626", linestyle=":", alpha=0.5)
    ax2.text(2015, lo[recent_mask][-1], "2015", fontsize=8, 
             color="#dc2626", ha="left", fontweight="bold")
    
    # Shade exponential growth phase
    growth_mask = (recent_dates >= -22000) & (recent_dates <= -1400)
    if growth_mask.sum() > 1:
        ax2.axvspan(-22000, -1400, alpha=0.05, color="#f59e0b")
        ax2.text(-12000, recent_ne.max() * 0.5, "~300× expansion", fontsize=8,
                 color="#b45309", ha="center", fontstyle="italic")

plt.suptitle("Bayesian Skyline Plot — L4.11 demographic history\n(TreeTime coalescent skyline, 825 L4.11 strains, 143 dated tips)",
             fontsize=12, fontweight="bold", y=1.04)
plt.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
print(f"\nFigure saved: {OUT_PNG}")

# Copy results
import shutil
for f in ["skyline.tsv", "timetree.nexus", "dates.tsv"]:
    src = f"/tmp/treetime_skyline_l411_v2/{f}"
    dst = f"/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/{f}"
    shutil.copy(src, dst)
    print(f"Copied: {f}")

shutil.copy(__file__, "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/skyline_plot.py")
print("\nDone!")
