#!/usr/bin/env python3
"""Utility module: identify L4.11.1 strains from report.json using SPDI marker.

L4.11.1 is defined by the core-exclusive SPDI NC_000962.3:1700209:G:A.
This replaces the hardcoded strain ID list with dynamic detection.
"""

import os

# L4.11.1 defining marker
L4111_SPDI = "NC_000962.3:1700209:G:A"

def get_l4111_strains(bdd_dir="BDD/L4.11"):
    """Return set of strain IDs that carry the L4.11.1 marker SPDI."""
    l4111 = set()
    for strain_id in os.listdir(bdd_dir):
        rpath = os.path.join(bdd_dir, strain_id, "NC_000962.3", "report.json")
        if not os.path.exists(rpath):
            continue
        # Fast text search instead of JSON parsing (much faster)
        with open(rpath) as f:
            content = f.read()
        if f'"spdi":"{L4111_SPDI}"' in content or f'"spdi": "{L4111_SPDI}"' in content:
            l4111.add(strain_id)
    return l4111

if __name__ == "__main__":
    strains = get_l4111_strains()
    print(f"L4.11.1 strains (n={len(strains)}): {sorted(strains)[:5]}...")
