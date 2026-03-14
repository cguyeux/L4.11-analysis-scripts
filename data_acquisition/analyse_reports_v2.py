#!/usr/bin/env python3
"""
Unified analysis of 822 L4.11 v2 reports.
Extracts: IS6110, missing genes, missing RDs, mapping stats, SNPs per drug gene.

Output: articles/L4.11/résultats/analysis_v2.json (consolidated per-strain data)
        articles/L4.11/résultats/analysis_v2_summary.txt (human-readable summary)

Run: python3 analyse_reports_v2.py
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).parent
REPORTS_DIR = ROOT / "reports_v2"
STRAINS_CSV = ROOT / "strains.csv"
WHO_CSV = ROOT.parent.parent.parent / "data" / "WHO_catalogue_v2.csv"
OUT_JSON = ROOT / "analysis_v2.json"
OUT_SUMMARY = ROOT / "analysis_v2_summary.txt"

# Drug resistance genes of interest
DR_GENES = {
    "rpoB", "katG", "inhA", "embB", "embC", "embA", "ubiA",
    "pncA", "rpsL", "rrs", "gyrA", "gyrB",
    "mmpR5", "Rv0678", "mmpL5",
    "rpoC", "rpoA",
    "ethA", "ethR",
    "fabG1",  # inhA promoter region
    "panD",
    "eis",
    "rplC",  # linezolid
    "atpE",  # bedaquiline
    "ddn",   # delamanid
    "fbiA", "fbiB", "fbiC", "fgd1",  # delamanid
    "tlyA",  # capreomycin
    "nat",   # arylamine N-acetyltransferase
}

# Virulence gene families (prefix-based)
VIRULENCE_PREFIXES = {
    "PE_PGRS": "PE/PPE", "PPE": "PE/PPE", "PE": "PE/PPE",
    "mce": "mce", "ecc": "ESX", "esp": "ESX", "esx": "ESX",
    "fad": "lipid", "pks": "lipid", "mmpL": "lipid",
    "fadD": "lipid", "fadE": "lipid",
    "vap": "TA", "maz": "TA", "rel": "TA",
    "sig": "sigma",
    "dos": "DosR", "dev": "DosR", "hsp": "DosR",
}


def classify_gene(gene_name):
    """Classify a gene into a virulence family."""
    if not gene_name:
        return None
    for prefix, family in sorted(VIRULENCE_PREFIXES.items(), key=lambda x: -len(x[0])):
        if gene_name.startswith(prefix):
            return family
    return None


def load_who_catalogue(path):
    """Load WHO v2 catalogue into dict: (gene, mutation) -> {drug: prediction}"""
    cat = defaultdict(dict)
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mut = row.get('MUTATION', '')
            drug = row.get('DRUG', '')
            pred = row.get('PREDICTION', '')
            if '@' in mut and drug:
                gene, change = mut.split('@', 1)
                cat[(gene, change)][drug] = pred
    return dict(cat)


def hgvs_to_garc(hgvs_p):
    """Convert SnpEff HGVS protein notation to GARC-like notation.
    E.g. p.Ser450Leu -> S450L
    """
    if not hgvs_p or not hgvs_p.startswith('p.'):
        return None

    # 3-letter to 1-letter
    aa3to1 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Ter': '*', 'fs': 'fs',
    }

    s = hgvs_p[2:]  # Remove 'p.'

    # Handle frameshifts
    if 'fs' in s:
        # e.g. p.Ile67fs -> I67fs
        import re
        m = re.match(r'([A-Z][a-z]{2})(\d+)fs', s)
        if m:
            aa = aa3to1.get(m.group(1), '?')
            return f'{aa}{m.group(2)}fs'
        return None

    # Handle stop_gained: p.Arg305*
    if s.endswith('*'):
        import re
        m = re.match(r'([A-Z][a-z]{2})(\d+)\*', s)
        if m:
            aa = aa3to1.get(m.group(1), '?')
            return f'{aa}{m.group(2)}*'
        return None

    # Standard missense: p.Ser450Leu -> S450L
    import re
    m = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', s)
    if m:
        aa1 = aa3to1.get(m.group(1), '?')
        aa2 = aa3to1.get(m.group(3), '?')
        return f'{aa1}{m.group(2)}{aa2}'

    return None


def analyse_report(filepath, who_cat):
    """Extract key data from a single v2 report."""
    with open(filepath) as f:
        d = json.load(f)

    strain_id = d.get('strain_id', filepath.stem)
    result = {'strain_id': strain_id}

    # Mapping stats
    ms = d.get('mapping_stats', {})
    result['mean_depth'] = ms.get('mean_depth', 0)
    result['coverage_pct'] = ms.get('covered_bases_percent', 0)

    # IS6110 counting
    is_all = d.get('insertion_sequences', [])
    is6110 = [x for x in is_all if x.get('name') == 'IS6110']
    is6110_novel = [x for x in is6110 if not x.get('is_reference', True)]
    result['is6110_total'] = len(is6110)
    result['is6110_novel'] = len(is6110_novel)
    result['is6110_positions'] = [x.get('position') for x in is6110]
    result['is_total'] = len(is_all)

    # IS1081 counting (for comparison)
    is1081 = [x for x in is_all if x.get('name') == 'IS1081']
    result['is1081_total'] = len(is1081)

    # Missing genes
    mg = d.get('missing_genes', [])
    result['missing_genes_count'] = len(mg)
    result['missing_genes'] = mg

    # Missing RD
    mrd = d.get('missing_rd', [])
    result['missing_rd_count'] = len(mrd)
    result['missing_rd'] = mrd

    # SNP analysis for drug resistance and virulence
    snps = d.get('snp', [])
    result['total_snps'] = len(snps)

    dr_mutations = []
    all_gene_mutations = []  # For pN/pS
    who_resistance = defaultdict(list)  # drug -> list of mutations

    for snp in snps:
        spdi = snp.get('spdi', '')
        for ann in snp.get('annotations', []):
            gene = ann.get('gene_name', '')
            hgvs_p = ann.get('hgvs_p')
            hgvs_c = ann.get('hgvs_c', '')
            impact = ann.get('impact', '')
            annotation_types = ann.get('annotation', [])
            locus = ann.get('gene_locus_tag', '')
            prot_pos = ann.get('protein_position')

            if not gene or gene == 'null':
                continue

            # Track all coding mutations for pN/pS
            is_syn = 'synonymous_variant' in annotation_types
            is_nonsyn = any(a in annotation_types for a in [
                'missense_variant', 'stop_gained', 'frameshift_variant',
                'start_lost', 'stop_lost', 'inframe_deletion', 'inframe_insertion'
            ])
            if is_syn or is_nonsyn:
                all_gene_mutations.append({
                    'gene': gene, 'locus': locus, 'spdi': spdi,
                    'hgvs_p': hgvs_p, 'hgvs_c': hgvs_c,
                    'is_syn': is_syn, 'is_nonsyn': is_nonsyn,
                    'family': classify_gene(gene),
                })

            # Drug resistance check
            clean_gene = gene.replace('gene-', '')
            # Map Rv0678 -> mmpR5
            gene_alias = {'Rv0678': 'mmpR5', 'Rv0206c': 'mmpL5'}
            check_gene = gene_alias.get(clean_gene, clean_gene)

            if check_gene in DR_GENES:
                garc = hgvs_to_garc(hgvs_p) if hgvs_p else None
                mut_info = {
                    'gene': check_gene, 'hgvs_p': hgvs_p, 'hgvs_c': hgvs_c,
                    'garc': garc, 'spdi': spdi, 'impact': impact,
                }

                # Check WHO catalogue
                if garc:
                    key = (check_gene, garc)
                    if key in who_cat:
                        for drug, pred in who_cat[key].items():
                            if pred in ('R', 'F'):  # Resistant or Fail
                                who_resistance[drug].append({
                                    'gene': check_gene,
                                    'mutation': garc,
                                    'prediction': pred,
                                })
                        mut_info['who_grading'] = who_cat[key]

                dr_mutations.append(mut_info)

    result['dr_mutations'] = dr_mutations
    result['who_resistance'] = dict(who_resistance)
    result['pnps_mutations'] = all_gene_mutations

    # Quality
    q = d.get('quality', {})
    af = q.get('after_filtering', {})
    result['total_reads'] = af.get('total_reads', 0)
    result['q30_rate'] = af.get('q30_rate', 0)

    return result


def main():
    print("Loading WHO catalogue v2...")
    who_cat = load_who_catalogue(WHO_CSV)
    print(f"  {len(who_cat)} unique gene@mutation entries")

    # Load strain metadata
    with open(STRAINS_CSV) as f:
        meta = {row['Id']: row for row in csv.DictReader(f)}

    # List available reports
    reports = sorted(REPORTS_DIR.glob("*.json"))
    print(f"Reports to analyse: {len(reports)}")

    results = []
    errors = []

    for i, rp in enumerate(reports, 1):
        try:
            r = analyse_report(rp, who_cat)
            # Add metadata
            sid = r['strain_id']
            if sid in meta:
                r['country'] = meta[sid].get('Country', '')
                r['bioproject'] = meta[sid].get('Bioproject', '')
                r['sublineage'] = ''  # Will be set based on tree
            results.append(r)
        except Exception as e:
            errors.append({'file': str(rp), 'error': str(e)})

        if i % 100 == 0 or i == len(reports):
            print(f"  [{i}/{len(reports)}] processed", flush=True)

    print(f"\nProcessed: {len(results)} | Errors: {len(errors)}")
    if errors:
        for e in errors[:5]:
            print(f"  ERROR: {e['file']}: {e['error']}")

    # Save consolidated results (stripped of heavy fields for JSON size)
    # Save lightweight version (no missing_genes lists, no pnps_mutations)
    light = []
    for r in results:
        lr = {k: v for k, v in r.items()
              if k not in ('missing_genes', 'pnps_mutations', 'is6110_positions')}
        light.append(lr)

    with open(OUT_JSON, 'w') as f:
        json.dump(light, f, indent=1)
    print(f"\nSaved lightweight JSON: {OUT_JSON} ({os.path.getsize(OUT_JSON)/1024:.0f} KB)")

    # Save full data for pN/pS analysis separately
    pnps_data = {}
    for r in results:
        pnps_data[r['strain_id']] = r.get('pnps_mutations', [])
    with open(ROOT / "pnps_v2.json", 'w') as f:
        json.dump(pnps_data, f)
    print(f"Saved pN/pS data: {ROOT / 'pnps_v2.json'}")

    # Save IS6110 positions for IS analysis
    is_data = {}
    for r in results:
        is_data[r['strain_id']] = r.get('is6110_positions', [])
        # Restore from original result
        for orig in results:
            if orig['strain_id'] == r['strain_id']:
                is_data[r['strain_id']] = [
                    x.get('position') for x in
                    json.load(open(REPORTS_DIR / f"{r['strain_id']}.json")).get('insertion_sequences', [])
                    if x.get('name') == 'IS6110'
                ] if False else r.get('is6110_positions', [])  # use cached
                break

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    depths = [r['mean_depth'] for r in results if r['mean_depth'] > 0]
    print(f"\nMapping: median depth = {sorted(depths)[len(depths)//2]:.1f}x")
    print(f"  range: {min(depths):.1f}x - {max(depths):.1f}x")

    is6110s = [r['is6110_total'] for r in results]
    print(f"\nIS6110: median = {sorted(is6110s)[len(is6110s)//2]}")
    print(f"  range: {min(is6110s)} - {max(is6110s)}")

    miss_genes = [r['missing_genes_count'] for r in results]
    print(f"\nMissing genes: median = {sorted(miss_genes)[len(miss_genes)//2]}")

    # Drug resistance summary
    dr_counts = Counter()
    for r in results:
        for drug in r.get('who_resistance', {}):
            dr_counts[drug] += 1
    print(f"\nWHO resistance (strains with ≥1 R/F mutation):")
    for drug, n in dr_counts.most_common(20):
        print(f"  {drug}: {n} ({n/len(results)*100:.1f}%)")

    # Save summary to text
    with open(OUT_SUMMARY, 'w') as f:
        f.write(f"V2 Analysis Summary ({len(results)} strains)\n")
        f.write(f"{'='*60}\n\n")
        f.write(f"Mapping: median depth = {sorted(depths)[len(depths)//2]:.1f}x\n")
        f.write(f"IS6110: median = {sorted(is6110s)[len(is6110s)//2]}\n")
        f.write(f"Missing genes: median = {sorted(miss_genes)[len(miss_genes)//2]}\n\n")
        f.write("WHO resistance:\n")
        for drug, n in dr_counts.most_common(20):
            f.write(f"  {drug}: {n} ({n/len(results)*100:.1f}%)\n")

    print(f"\nSummary saved: {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
