"""
Phase 0b: rules-compliance sequence audit (structures-free).
Rule 4: designs must be ">=10 aa different from known binders (scFv-2, scFv-4, VHB937 CDRs)".
Reference sequences from the addendum Appendix (Szykowska 2021 scFv-2/-4; Novartis MOR44698,
WO2020079580A1 — a lead clone of the VHB937 patent family, NOT confirmed to be VHB937 itself).

For each of the 141 designs, vs each reference (full scFv-2, full scFv-4, each CDR):
  - Smith-Waterman local alignment (BLOSUM62): max local %identity, aligned length, score
  - Longest common substring (exact match) — the interpretable rules number
  - Global Levenshtein to the full scFvs
Flag: LCS >= 5 aa to any CDR, OR >= 70% local identity over >= 10 aa to any reference.
"""
import pandas as pd
from pathlib import Path
from Bio.Align import PairwiseAligner, substitution_matrices
import Levenshtein

ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"

SCFV2 = "SMEVQLLESGGGLVQPGGSLRLSCAASGFTFYSSYMGWVRQAPGKGLEWVSYISSSGSSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARVGGYYSWGNGIDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYGVYYPFTFGQGTKLEIK"
SCFV4 = "SMEVQLLESGGGLVQPGGSLRLSCAASGFTFSYYYMGWVRQAPGKGLEWVSGISPSSGYTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARYYYGYYYSHMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRSGLHTFGQGTKLEIK"
CDRS = {
    "scFv2_H1": "YSSYMG", "scFv2_H2": "YISSSGSST", "scFv2_H3": "VGGYYSWGNGIDY",
    "scFv2_L1": "RASQSISSYLN", "scFv2_L2": "AASSLQS", "scFv2_L3": "QQYGVYYPFT",
    "scFv4_H1": "SYYYMG", "scFv4_H2": "GISPSSGYT", "scFv4_H3": "YYYGYYYSHMDY",
    "scFv4_L1": "RASQSISSYLN", "scFv4_L2": "AASSLQS", "scFv4_L3": "QQSRSGLHT",
    "MOR44698_H1": "GYTFTGYHMS", "MOR44698_H2": "VINPVSGNTVYAQKFQG", "MOR44698_H3": "IPSYTYAFDY",
    "MOR44698_L1": "RASQDISNYLA", "MOR44698_L2": "RASSLQS", "MOR44698_L3": "FQYRHMPSQT",
}
FULL = {"scFv2_full": SCFV2, "scFv4_full": SCFV4}

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "local"
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

def lcs(a, b):
    """longest common substring length (exact)."""
    m, n = len(a), len(b)
    dp = [0]*(n+1); best = 0
    for i in range(1, m+1):
        prev = 0
        for j in range(1, n+1):
            tmp = dp[j]
            if a[i-1] == b[j-1]:
                dp[j] = prev + 1; best = max(best, dp[j])
            else:
                dp[j] = 0
            prev = tmp
    return best

def local_identity(design, ref):
    """best local alignment; return (%id over aligned block, aligned length, score)."""
    aln = aligner.align(design, ref)[0]
    # count identities and aligned columns
    idc = 0; length = 0
    for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
        block = e1 - s1
        length += block
        idc += sum(1 for k in range(block) if design[s1+k] == ref[s2+k])
    pct = 100*idc/length if length else 0
    return pct, length, aln.score

df = pd.read_parquet(ROOT / "data/designs.parquet")
rows = []
for _, r in df.iterrows():
    seq = r.sequence; rec = {"design_id": r.design_id, "name": r["name"], "cohort": r.cohort,
                             "team": r.team, "length": len(seq)}
    # full scFv: local identity + Levenshtein
    for name, ref in FULL.items():
        pct, length, score = local_identity(seq, ref)
        rec[f"{name}_locid"] = round(pct, 1); rec[f"{name}_loclen"] = length
        rec[f"{name}_lev"] = Levenshtein.distance(seq, ref)
    # CDRs: LCS + local identity
    max_cdr_lcs = 0; max_cdr_lcs_which = ""
    for name, ref in CDRS.items():
        l = lcs(seq, ref)
        rec[f"{name}_lcs"] = l
        if l > max_cdr_lcs: max_cdr_lcs = l; max_cdr_lcs_which = name
    rec["max_cdr_lcs"] = max_cdr_lcs; rec["max_cdr_lcs_which"] = max_cdr_lcs_which
    # flags
    flag_lcs = max_cdr_lcs >= 5
    flag_locid = any(rec[f"{n}_locid"] >= 70 and rec[f"{n}_loclen"] >= 10 for n in FULL)
    rec["FLAG_lcs_ge5_cdr"] = flag_lcs
    rec["FLAG_locid_ge70_over10"] = flag_locid
    rec["FLAG_any"] = flag_lcs or flag_locid
    rows.append(rec)

out = pd.DataFrame(rows)
out.to_csv(RES / "rules_compliance.csv", index=False)
print(f"n designs: {len(out)}")
print(f"\nFLAG_lcs_ge5_cdr:   {int(out.FLAG_lcs_ge5_cdr.sum())}")
print(f"FLAG_locid_ge70_over10: {int(out.FLAG_locid_ge70_over10.sum())}")
print(f"FLAG_any:           {int(out.FLAG_any.sum())}")
print(f"\nmax CDR-LCS distribution:")
print(out.max_cdr_lcs.value_counts().sort_index())
print(f"\nmax local identity to full scFv-2/-4 (top 8 designs):")
out["max_full_locid"] = out[["scFv2_full_locid","scFv4_full_locid"]].max(axis=1)
print(out.sort_values("max_full_locid",ascending=False).head(8)[
    ["design_id","name","cohort","length","scFv2_full_locid","scFv4_full_locid","max_cdr_lcs","max_cdr_lcs_which"]].to_string(index=False))
if out.FLAG_any.any():
    print("\n=== FLAGGED designs ===")
    print(out[out.FLAG_any][["design_id","name","cohort","max_full_locid","max_cdr_lcs","max_cdr_lcs_which"]].to_string(index=False))
print("\nwrote results/rules_compliance.csv")
