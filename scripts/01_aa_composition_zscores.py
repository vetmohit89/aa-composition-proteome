#!/usr/bin/env python3
"""
Calculate amino acid percentage and z-score for each protein in a FASTA file.

Outputs:
  <prefix>_aa_percentages.tsv
  <prefix>_aa_zscores.tsv

Filtering:
  - skips empty sequences
  - skips sequences containing non-standard amino acids
    (anything outside ACDEFGHIKLMNPQRSTVWY)
"""

import argparse
import math
import re
import sys
from Bio import SeqIO

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
VALID_AA_SET = set(AMINO_ACIDS)


def is_usable_protein_sequence(seq):
    """
    Return True if the sequence is a usable protein sequence:
    - not empty
    - contains only the 20 standard amino acids
    """
    if not seq:
        return False
    return set(seq).issubset(VALID_AA_SET)


def extract_gene_symbol(description):
    """
    Extract gene symbol from UniProt header if present.
    Example: GN=TP53
    """
    match = re.search(r"\bGN=(\S+)", description)
    return match.group(1) if match else "NA"


def calculate_percentages(seq):
    """
    Calculate amino acid percentages for a protein sequence.
    """
    length = len(seq)
    return [(seq.count(aa) / length) * 100 for aa in AMINO_ACIDS]


def main():
    parser = argparse.ArgumentParser(
        description="Calculate amino acid percentages and z-scores from a FASTA file."
    )
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument(
        "--output-prefix",
        default="proteome",
        help="Prefix for output files (default: proteome)"
    )
    args = parser.parse_args()

    proteins = []

    print("Reading FASTA and calculating amino acid percentages...")

    processed = 0
    skipped_empty = 0
    skipped_nonstandard = 0

    for record in SeqIO.parse(args.fasta_file, "fasta"):
        seq = str(record.seq).upper().replace("*", "")
        if not seq:
            skipped_empty += 1
            continue

        if not is_usable_protein_sequence(seq):
            skipped_nonstandard += 1
            continue

        length = len(seq)
        gene_symbol = extract_gene_symbol(record.description)
        percentages = calculate_percentages(seq)

        proteins.append({
            "Protein_ID": record.id,
            "Gene_Symbol": gene_symbol,
            "Length": length,
            "Percentages": percentages
        })

        processed += 1
        if processed % 1000 == 0:
            print(f"  Processed {processed} proteins...")

    if not proteins:
        print("No usable protein sequences found.", file=sys.stderr)
        sys.exit(1)

    print(f"Total usable proteins processed: {len(proteins)}")
    print(f"Skipped empty sequences: {skipped_empty}")
    print(f"Skipped non-standard sequences: {skipped_nonstandard}")

    # Calculate z-score parameters for each amino acid
    means = []
    std_devs = []

    for aa_index in range(len(AMINO_ACIDS)):
        values = [p["Percentages"][aa_index] for p in proteins]
        mean = sum(values) / len(values)
        variance = sum((x - mean) ** 2 for x in values) / len(values)
        std_dev = math.sqrt(variance)

        means.append(mean)
        std_devs.append(std_dev)

    # Write amino acid percentage file
    percent_file = f"{args.output_prefix}_aa_percentages.tsv"
    print(f"Writing percentages to {percent_file}...")
    with open(percent_file, "w") as f:
        header = ["Protein_ID", "Gene_Symbol", "Length"] + AMINO_ACIDS
        f.write("\t".join(header) + "\n")

        for protein in proteins:
            row = [
                protein["Protein_ID"],
                protein["Gene_Symbol"],
                str(protein["Length"])
            ] + [f"{x:.3f}" for x in protein["Percentages"]]
            f.write("\t".join(row) + "\n")

    # Write amino acid z-score file
    zscore_file = f"{args.output_prefix}_aa_zscores.tsv"
    print(f"Writing z-scores to {zscore_file}...")
    with open(zscore_file, "w") as f:
        header = ["Protein_ID", "Gene_Symbol", "Length"] + AMINO_ACIDS
        f.write("\t".join(header) + "\n")

        for protein in proteins:
            zscores = []
            for i, value in enumerate(protein["Percentages"]):
                if std_devs[i] == 0:
                    z = 0.0
                else:
                    z = (value - means[i]) / std_devs[i]
                zscores.append(f"{z:.3f}")

            row = [
                protein["Protein_ID"],
                protein["Gene_Symbol"],
                str(protein["Length"])
            ] + zscores
            f.write("\t".join(row) + "\n")

    print("Done.")
    print(f"Created:\n  {percent_file}\n  {zscore_file}")


if __name__ == "__main__":
    main()
