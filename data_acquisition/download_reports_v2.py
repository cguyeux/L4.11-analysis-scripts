#!/usr/bin/env python3
"""
Download all v2 reports for L4.11 strains via the TBannotator v2 API.
Saves to articles/L4.11/résultats/reports_v2/<SRA_ID>.json

Usage:
    python3 download_reports_v2.py [--workers 4] [--resume]
"""
import csv
import json
import os
import sys
import time
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Config
BASE_URL = "https://darthos.freeboxos.fr/mcp/download/report"
STRAINS_CSV = Path(__file__).parent / "strains.csv"
OUTPUT_DIR = Path(__file__).parent / "reports_v2"
MAX_RETRIES = 3
TIMEOUT = 120  # seconds per request


def download_report(sra_id: str, output_dir: Path, retry: int = 0) -> dict:
    """Download a single report.json for a given SRA ID."""
    outfile = output_dir / f"{sra_id}.json"
    
    # Skip if already downloaded and valid
    if outfile.exists() and outfile.stat().st_size > 100:
        try:
            with open(outfile) as f:
                data = json.load(f)
            if data.get("strain_id") == sra_id:
                return {"id": sra_id, "status": "skipped", "size": outfile.stat().st_size}
        except (json.JSONDecodeError, KeyError):
            pass  # Re-download corrupt file
    
    url = f"{BASE_URL}/{sra_id}"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            data = resp.read()
        
        # Validate JSON
        parsed = json.loads(data)
        if parsed.get("strain_id") != sra_id:
            # Some reports may have different strain_id format
            pass
        
        with open(outfile, "wb") as f:
            f.write(data)
        
        return {"id": sra_id, "status": "ok", "size": len(data), "snps": len(parsed.get("snp", []))}
    
    except Exception as e:
        if retry < MAX_RETRIES:
            time.sleep(2 ** retry)  # Exponential backoff
            return download_report(sra_id, output_dir, retry + 1)
        return {"id": sra_id, "status": "error", "error": str(e)}


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Download v2 reports for L4.11")
    parser.add_argument("--workers", type=int, default=4, help="Parallel downloads")
    parser.add_argument("--resume", action="store_true", help="Skip already downloaded")
    args = parser.parse_args()
    
    # Read strain IDs
    with open(STRAINS_CSV) as f:
        reader = csv.DictReader(f)
        sra_ids = [row["Id"] for row in reader]
    
    print(f"Strains to download: {len(sra_ids)}")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Check already downloaded
    existing = {f.stem for f in OUTPUT_DIR.glob("*.json")}
    if args.resume:
        remaining = [s for s in sra_ids if s not in existing]
        print(f"Already downloaded: {len(existing)}, remaining: {len(remaining)}")
    else:
        remaining = sra_ids
    
    if not remaining:
        print("All reports already downloaded!")
        return
    
    # Download with progress
    ok = 0
    errors = []
    t0 = time.time()
    
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download_report, sid, OUTPUT_DIR): sid for sid in remaining}
        
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            if result["status"] == "error":
                errors.append(result)
                sym = "✗"
            elif result["status"] == "skipped":
                ok += 1
                sym = "→"
            else:
                ok += 1
                sym = "✓"
            
            elapsed = time.time() - t0
            rate = i / elapsed if elapsed > 0 else 0
            eta = (len(remaining) - i) / rate if rate > 0 else 0
            
            print(f"  [{i}/{len(remaining)}] {sym} {result['id']} "
                  f"({result.get('size', 0) / 1024:.0f} KB) "
                  f"[{rate:.1f}/s, ETA {eta:.0f}s]",
                  flush=True)
    
    # Summary
    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Downloaded: {ok} | Errors: {len(errors)} | Time: {elapsed:.0f}s")
    
    if errors:
        print("\nFailed downloads:")
        for e in errors:
            print(f"  {e['id']}: {e.get('error', '?')}")
        
        with open(OUTPUT_DIR / "errors.txt", "w") as f:
            for e in errors:
                f.write(f"{e['id']}\t{e.get('error', '?')}\n")
    
    # Final validation
    total_valid = 0
    for f in OUTPUT_DIR.glob("*.json"):
        try:
            with open(f) as fh:
                json.load(fh)
            total_valid += 1
        except:
            pass
    
    print(f"\nValid reports: {total_valid}/{len(sra_ids)}")


if __name__ == "__main__":
    main()
