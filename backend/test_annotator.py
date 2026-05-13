"""Test GRN annotation with bacteriorhodopsin."""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from main import GRNAnnotator

# Bacteriorhodopsin (1C3W) - K216 should be at 7.50
BR_SEQ = (
    "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLL"
    "GYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVW"
    "WAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDV"
    "SAKV GFGLILLRSRAIFGEAEAPEPSAGDGAAATSD"
).replace(" ", "")

if __name__ == "__main__":
    print("Initializing annotator...")
    ann = GRNAnnotator()
    print(f"Refs: {len(ann.ref_seqs)}, Strict GRNs: {len(ann.grns_str_strict)}")

    print(f"\nAnnotating BR ({len(BR_SEQ)} aa)...")
    result = ann.annotate(BR_SEQ, "bacteriorhodopsin")
    print(f"Annotated {len(result)} / {len(BR_SEQ)} positions\n")

    # Check Schiff base
    sb = [r for r in result if r["grn"] == "7.50"]
    if sb:
        print(f"7.50 (Schiff base): {sb[0]['residue']}{sb[0]['position']}")
        assert sb[0]["residue"] == "K", f"Expected K, got {sb[0]['residue']}"
        print("PASS: Schiff base lysine correct!\n")

    # Conserved midpoints
    print("TM midpoints (x.50):")
    for tm in range(1, 8):
        pos = f"{tm}.50"
        m = [r for r in result if r["grn"] == pos]
        print(f"  {pos}: {m[0]['residue']}{m[0]['position']}" if m else f"  {pos}: missing")

    # First/last 10
    print("\nFirst 10:")
    for r in result[:10]:
        print(f"  {r['residue']}{r['position']:>4d} -> {r['grn']}")
    print(f"\nLast 10:")
    for r in result[-10:]:
        print(f"  {r['residue']}{r['position']:>4d} -> {r['grn']}")
