# Analysis Scripts for the L4.11 Article

Scripts used in:

> **Genomic characterisation of *Mycobacterium tuberculosis* lineage L4.11
> reveals two geographically distinct sub-lineages with contrasting
> drug-resistance profiles.**
>
> Guyeux C., Leclercq C., Revon G., Sola C.

---

## Requirements

| Package | Tested version | Purpose |
|---------|---------------|---------|
| Python | ≥ 3.10 | All scripts |
| BioPython | ≥ 1.83 | Tree parsing, Entrez queries |
| matplotlib | ≥ 3.8 | Figures |
| numpy / scipy | ≥ 1.26 / 1.12 | Statistics, distances |
| scikit-learn | ≥ 1.4 | PCoA, clustering |
| cartopy | ≥ 0.22 | Geographic maps |
| plotly | ≥ 5.18 | Interactive plots |
| ete3 | ≥ 3.1.3 | Tree manipulation |
| statsmodels | ≥ 0.14 | Statistical tests |

Install all dependencies:

```bash
pip install biopython matplotlib numpy scipy scikit-learn cartopy plotly ete3 statsmodels
```

## Input Data

All scripts expect the following data files in `../résultats/`:

- `reports_v2/` — Individual TBannotator v2 JSON reports (one per strain)
- `analysis_v2.json` — Aggregated analysis data
- `T3_l411_only.nwk` — RAxML-NG maximum-likelihood tree (Newick)
- `treetime_output/` — TreeTime molecular dating output

The raw sequencing data are available from the NCBI Sequence Read Archive
under the BioProject accessions listed in Table S1.

---

## Directory Structure

### `data_acquisition/` — Data Retrieval and QC

| Script | Description |
|--------|-------------|
| `download_reports_v2.py` | Batch download of TBannotator v2 reports via API |
| `analyse_reports_v2.py` | Unified parsing of 822 reports → `analysis_v2.json` |
| `quality_analysis.py` | Mapping statistics, coverage, per-BioProject QC |

### `phylogenetics/` — Lineage Definition and Sub-lineage Structure

| Script | Description |
|--------|-------------|
| `breakdown_sublineage_v2.py` | Sub-lineage assignment (L4.11.1 vs L4.11.2) from SPDI markers |
| `l4111_marker.py` | Utility: identify L4.11.1 strains from report.json |
| `heatmap_distances.py` | Pairwise patristic distance heatmap (Figure 3) |
| `root_to_tip.py` | Root-to-tip regression for temporal signal |
| `pcoa_l411.py` | PCoA ordination by sub-lineage |
| `pcoa_l411_intra.py` | PCoA ordination: L4.11.2 vs core L4.11.1 vs Proto-L4.11.1 |

### `geographic/` — Geographic Distribution

| Script | Description |
|--------|-------------|
| `map_l411_v2.py` | Two-panel world map of L4.11 strain distribution (Figure 4) |
| `stat_country_sublineage.py` | Fisher exact tests for country × sub-lineage associations |

### `drug_resistance/` — Drug Resistance Profiling

| Script | Description |
|--------|-------------|
| `dr_profile.py` | Complete DR mutational profile across all drug-resistance genes |
| `coresistance_analysis.py` | Co-resistance pattern analysis (INH+RIF, MDR, XDR) |
| `who_resistance_analysis.py` | Cross-reference with WHO Catalogue v2.1 mutations |
| `compensatory_analysis.py` | Compensatory mutation analysis: rpoB/rpoC/rpoA co-occurrence |

### `structural_variation/` — IS6110, RD, and CNV Analysis

| Script | Description |
|--------|-------------|
| `is_analysis.py` | IS6110 in silico fingerprinting: L4.11.1 vs L4.11.2 profiles |
| `is_gene_map.py` | Map IS6110 insertion positions to genes (GFF3) |
| `is_rd_analysis.py` | Combined IS and RD analysis across sub-lineages |
| `rd_analysis.py` | RD-typing validation (large_rd and missing_rd profiles) |
| `missing_genes.py` | Gene deletions from coverage data |
| `cnv_analysis_v2.py` | Copy-number variation detection from gene coverage |
| `insert_size_analysis.py` | Insert size distribution by BioProject |

### `virulence_selection/` — Virulence and Selection Pressure

| Script | Description |
|--------|-------------|
| `convergent_analysis.py` | Convergent non-synonymous mutations in virulence genes |
| `dnds_analysis.py` | Per-gene dN/dS ratio analysis |
| `homoplasy_analysis.py` | Homoplasy detection: independently emerged mutations in both sub-lineages |

### `epidemiology/` — Transmission and Epidemiological Analysis

| Script | Description |
|--------|-------------|
| `transmission_clustering.py` | Pairwise SNP distance and transmission cluster identification |
| `mixed_analysis.py` | Mixed infection detection via allelic ratios |
| `spoligotype_analysis.py` | Spoligotype pattern analysis |

### `figures/` — Figure Generation

| Script | Description |
|--------|-------------|
| `gen_phylo_fig2.py` | Circular phylogenetic tree (Figure 2) |
| `regenerate_figures_v2.py` | Regenerate all figures from `analysis_v2.json` |
| `skyline_plot.py` | Bayesian Skyline Plot from TreeTime output |

---

## Reproducibility

Each script can be run independently. The typical workflow is:

```bash
# 1. Download reports (requires TBannotator API access)
cd data_acquisition
python download_reports_v2.py

# 2. Aggregate into unified dataset
python analyse_reports_v2.py

# 3. Run any analysis script
cd ../drug_resistance
python dr_profile.py
```

Most scripts read from `../résultats/reports_v2/` or `../résultats/analysis_v2.json`
and produce output figures and tables to `../article/figures/` or `stdout`.

## Citation

If you use these scripts, please cite:

```
Guyeux C., Leclercq C., Revon G., Sola C. (2026).
Genomic characterisation of Mycobacterium tuberculosis lineage L4.11
reveals two geographically distinct sub-lineages with contrasting
drug-resistance profiles.
```

## License

These scripts are released under the MIT License.
