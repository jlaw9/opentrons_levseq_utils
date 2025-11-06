""" Script to map DNA variants from LevSeq to amino acid changes
Author: Jeff Law
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

pd.options.mode.chained_assignment = None  # default='warn'


# map the DNA variants to the amino acid changes
def map_dna_to_aa(dna_seq, seq_idx=None):
    seq_idx_str = f"'{seq_idx}'" if seq_idx is not None else ""
    if dna_seq is None or pd.isnull(dna_seq):
        return np.nan
    aa_seq = []
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if codon in codontab:
            aa = codontab[codon]
            if aa == '*':
                break
            aa_seq.append(aa)
        else:
            print(f"Unknown codon: {codon}")
            aa_seq.append('X')
    if aa == "*" and i < len(dna_seq) - 3:
        print(f"seq idx {seq_idx_str}: stop codon encountered at position {i} in sequence {seq_idx_str} ({i//3} AA in sequence).")
    aa_seq.append('*')
        #print(f"Warning: Sequence is longer than expected. Remaining sequence: {dna_seq[i:]}")
    return ''.join(aa_seq)


def apply_dna_mutations(mutations, offset=0):
    # Create a list to hold the mutated sequence
    mutated_seq = list(refseq)

    # Apply the mutation
    for mut in mutations.split('_'):
        ref_bp, pos, new_bp = mut[0], mut[1:-1], mut[-1]
        if 'DEL' in mut or 'INS' in mut:
            pos = mut[1:-3]
        pos = int(pos) - 1  # Convert to 0-based index
        # if a barcode was included at the start, then remove that
        pos -= offset
        if pos >= len(refseq):
            print(f"WARNING: {pos = } > {len(refseq) = } for '{mutations}'. Skipping")
            return np.nan

        if refseq[pos] != ref_bp:
            print(f"WARNING: '{refseq[pos]}' != '{ref_bp}'. "
                  f"Reference base  does not match at position {pos + 1}. Mutations: {mutations}")
            return ''.join(mutated_seq)
        # assert refseq[pos] == ref_bp, f"'{refseq[pos]}' != '{ref_bp}'. Reference base  does not match at position {pos + 1}. Mutations: {mutations}"

        try:
            # need to handle deletions and insertions
            # e.g., T551DEL, C651INS
            if "DEL" in mut:
                # Deletion: remove the base
                mutated_seq.pop(pos)
            elif "INS" in mut:
                # Insertion: copy the current base (it always accompanies a substitution mutation)
                # e.g., C651INS, C651T
                mutated_seq.insert(pos, ref_bp)
            else:
                # Substitution: replace the base
                mutated_seq[pos] = new_bp
        except IndexError:
            print(f"Error when handling mutations: '{mutations}'")
            return np.nan

    return ''.join(mutated_seq)


def get_aa_mutations(ref_seq, mut_seq):
    """ List the mutations between two amino acid sequences. """
    # Convert the sequences to lists of amino acids
    ref_aa = list(ref_seq)
    mut_aa = list(mut_seq)

    # TODO better handle insertions and deletions
    if len(ref_aa) != len(mut_aa):
        return np.nan

    # Find the positions of the mutations
    mutations = []
    for i in range(len(mut_aa)):
        if ref_aa[i] != mut_aa[i]:
            mutations += [f"{ref_aa[i]}{i+1}{mut_aa[i]}"]
            #mutations.append((i+1, ref_aa[i], mut_aa[i]))

    return '_'.join(mutations)


# codon mapping
codontab = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',    # serine
    'TTC': 'F', 'TTT': 'F',    # Fenilalanine
    'TTA': 'L', 'TTG': 'L',    # Leucine
    'TAC': 'Y', 'TAT': 'Y',    # Tirosine
    'TAA': '*', 'TAG': '*',    # Stop
    'TGC': 'C', 'TGT': 'C',    # Cisteine
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofane
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucine
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',    # Proline
    'CAC': 'H', 'CAT': 'H',    # Histidine
    'CAA': 'Q', 'CAG': 'Q',    # Glutamine
    'CGA': 'R', 'CGC': 'R',    # arginine
    'CGG': 'R', 'CGT': 'R',    # arginine
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',    # Treonine
    'AAC': 'N', 'AAT': 'N',    # asparagine
    'AAA': 'K', 'AAG': 'K',    # lisine
    'AGC': 'S', 'AGT': 'S',    # serine
    'AGA': 'R', 'AGG': 'R',    # arginine
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',    # valine
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',    # alanine
    'GAC': 'D', 'GAT': 'D',    # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',    # Glutamic Acid
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'     # glicine
}


if __name__ == "__main__":
    # setup the argument parser
    parser = argparse.ArgumentParser(description="Map LevSeq DNA variants to amino acid changes")
    parser.add_argument("--var_file", type=str, help="Path to the LevSeq variants file (variants.csv)")
    parser.add_argument("--run_name", type=str, default=None, help="Optional run name to rename output files")
    parser.add_argument("--offset", type=int, default=0, help="Optional offset to apply to mutation positions (e.g., barcode length)")
    args = parser.parse_args()  

    offset = args.offset
    # run_dir = Path(args.run_dir)

    var_file = Path(args.var_file)
    if not var_file.exists():
        raise FileNotFoundError(f"Variants file not found: {var_file}")
    print(f"Processing '{var_file}'")
    df = pd.read_csv(var_file)
    print(f"\t{len(df)} rows")
    print(df.head(2))

    print("")
    print(f"# unique reference sequences: {len(df['refseq'].unique())}")
    refseq = df['refseq'].unique()[0]
    aa_refseq = map_dna_to_aa(refseq)
    # print(aa_refseq)

    df_mut = df[df['P value'] < 0.05]
    print(f"# rows with a p-value < 0.05: {len(df_mut)}")
    # print(df_mut.Variant.value_counts())

    print("")
    print("Applying mutations")
    df_mut["AA_seq"] = df_mut.apply(lambda row: 
                                    map_dna_to_aa(
                                        apply_dna_mutations(row['Variant'], offset=offset), 
                                        seq_idx=row.name), 
                                    axis=1)
    df_mut = df_mut.dropna(subset="AA_seq")
    df_mut = df_mut[~df_mut.AA_seq.str.contains("X")]
    print(f"# rows that successfully applied sequence mutations and mapped to AA: {len(df_mut)}")

    print("")
    print("Extracting AA mutations")
    df_mut["AA_mutations"] = df_mut.AA_seq.apply(lambda x: get_aa_mutations(aa_refseq, x))
    print(df_mut.AA_mutations.value_counts())
    print("")

    # if a run name is given, use that
    rename_file = var_file
    if args.run_name is not None:
        rename_file = str(var_file).replace("variants.csv", f"{args.run_name}_levseq_variants.csv")
        print(f"Copying the variants file to include the run name: {rename_file}")
        df.to_csv(rename_file, index=False)

    out_file = str(rename_file).replace(".csv", "_aa_mut.csv")
    print(f"Writing AA mutations file: {out_file}")
    df_mut.to_csv(out_file, index=False)
