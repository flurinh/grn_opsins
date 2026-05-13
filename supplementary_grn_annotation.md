# Supplementary Methods: Runtime GRN annotation via the `protos` webservice

The annotation procedure is distinct from the reference table itself. The reference table — our core contribution — is the manual curation of 129 microbial rhodopsins spanning the family's structural and functional diversity, with canonical anchors (x.50) placed in each transmembrane helix by retinal proximity and sequence conservation. The annotation procedure is the application that consumes that table: a propagation tool that transfers the curated numbering onto a new sequence. Its quality is bounded by the reference: queries close to a reference entry yield reliable results; distant queries still return an annotation, but with low coverage and alignment scores that should prompt cautious interpretation. At runtime the algorithm has no knowledge of retinal, structure, or conservation; it transfers GRNs by sequence alignment to the closest reference, then fills unassigned positions by fixed geometric rules. Because the reference covers pumps, cation and anion channels, sensory rhodopsins, heliorhodopsins, and viral and synthetic variants, most queries find a useful neighbour — broadening the reference table directly improves annotation behaviour.

Each query is compared against the 129 reference sequences with MMseqs2 to identify its closest reference, then globally aligned to that reference with BLOSUM62. Reference selection is per-query, so a channelrhodopsin is typically annotated against a channelrhodopsin, a proton pump against a pump, and so on: no single entry acts as a global pivot, and annotation error is bounded by local rather than global sequence distance.

Annotation is transferred in two stages. First, a seed is built column by column: every query residue aligned to a reference residue carrying a canonical GRN inherits that GRN. This stage alone annotates most of the query when the reference is close, leaving three kinds of residue unassigned — termini, canonical GRNs lost to alignment gaps, and residues in insertions or inter-helical loops. The second stage fills these by fixed rules following Isberg et al. (2015). Terminal residues receive N- or C-tail labels. Unassigned canonical GRNs are placed by pivoting from the nearest assigned GRN in the same helix; a position that cannot be filled is left unassigned, and the resulting skipped number signals a 3_10-helix-like compression. An unassigned stretch whose two flanks both lie inside the same helix is an insertion and receives a three-decimal suffix on the nearest flank, marking a π-helix-like extension.

Inter-helical loop residues receive labels of the form HxHy.<position>, where Hx is the nearer transmembrane helix, Hy is the other flanking helix, and <position> is a three-digit counter growing outward from the Hx boundary. The first digit identifies the side of the loop the residue sits on, not the loop as a whole: a residue in the TM6–TM7 loop is labelled 67.<position> on the TM6 side and 76.<position> on the TM7 side. The boundary between the two sides is the midpoint of the unassigned stretch, so every residue receives exactly one label with no overlap or gap. For example, the TM6–TM7 loop of HsBR (PDB 1c3w) is labelled G192=67.001, S193=67.002, E194=67.003, G195=67.004, A196=67.005 on the TM6 side and G197=76.004, I198=76.003, V199=76.002, P200=76.001 on the TM7 side. E194 — discussed in the HsBR literature as both a proton-release-group member and a possible distorted-TM6 residue — is carried as 67.003: whether it belongs to TM6 or the loop is a secondary-structure question, not a numbering question, and the label is fixed regardless.

The service returns a table of query residues keyed by canonical GRN, with per-query metadata reporting the selected reference, the alignment score, and the coverage fraction. Coverage fraction is the user-facing confidence indicator; low values flag queries whose fold diverges substantially from every reference. Extending the service to a new clade is a matter of adding a curated entry to the reference table, not of modifying the algorithm. Expressed as pseudocode:

```
INPUT:  query_sequences     # dict: query_id -> amino-acid sequence
        reference_table     # 129 curated references x canonical GRN columns

for query_id, query_seq in query_sequences:
    reference_id = mmseqs2_top_hit(query_seq, reference_sequences)
    alignment    = biopython_align(query_seq, reference_sequences[reference_id],
                                   matrix=BLOSUM62)
    ref_grn      = grn_map_of_reference_row(reference_table, reference_id)

    # stage 1: seed transfer
    seed = {}
    for query_pos, ref_pos in aligned_pairs(alignment):
        grn = ref_grn.get(ref_pos)
        if grn is not None:
            seed[query_pos] = grn

    # stage 2: expansion
    assigned = seed
    assigned.update(n_tail_labels(query_seq, seed))
    assigned.update(c_tail_labels(query_seq, seed))

    for grn in canonical_grns_not_in(assigned):
        pivot     = nearest_assigned_grn_in_same_helix(grn, assigned)
        candidate = position_of(pivot) + offset(grn, pivot)
        if candidate is unassigned and within query:
            assigned[candidate] = grn
        # else leave unassigned (skipped number = 3_10-like compression)

    for interval in contiguous_unassigned_intervals(assigned):
        if both_flanks_in_same_helix(interval):
            for rank, pos in enumerate(interval, start=1):
                assigned[pos] = flank.grn + 0.001 * rank    # e.g. 3.52 -> 3.521
        else:
            for pos in n_half(interval):
                assigned[pos] = f"{Hx}{Hy}.{distance_to_Hx:03d}"
            for pos in c_half(interval):
                assigned[pos] = f"{Hy}{Hx}.{distance_to_Hy:03d}"

    emit(query_id, assigned)
```
