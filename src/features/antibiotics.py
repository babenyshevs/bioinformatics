from typing import Dict, List


def dna_to_rna(dna: str) -> str:
    """
    Transcribe a DNA sequence into RNA.

    Args:
        dna (str): DNA sequence consisting of 'A', 'T', 'C', and 'G'.

    Returns:
        str: RNA sequence obtained by replacing 'T' with 'U'.

    Example:
        >>> dna_to_rna("ATCGATCGATCG")
        'AUCGAUCGAUCG'
    """
    return dna.replace("T", "U")


def translate_rna_to_protein(rna: str) -> str:
    """
    Translate an RNA sequence into a protein sequence (amino acids).

    Args:
        rna (str): RNA sequence.

    Returns:
        str: Protein sequence (amino acids).

    Example:
        >>> translate_rna_to_protein("AUCGAUCGAUCG")
        'IDRIDR'
    """
    codon_table = {
        "UUU": "F",
        "UUC": "F",
        "UUA": "L",
        "UUG": "L",
        "UCU": "S",
        "UCC": "S",
        "UCA": "S",
        "UCG": "S",
        "UAU": "Y",
        "UAC": "Y",
        "UAA": "*",
        "UAG": "*",
        "UGU": "C",
        "UGC": "C",
        "UGA": "*",
        "UGG": "W",
        "CUU": "L",
        "CUC": "L",
        "CUA": "L",
        "CUG": "L",
        "CCU": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAU": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGU": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AUU": "I",
        "AUC": "I",
        "AUA": "I",
        "AUG": "M",
        "ACU": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAU": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGU": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GUU": "V",
        "GUC": "V",
        "GUA": "V",
        "GUG": "V",
        "GCU": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAU": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGU": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    protein = []
    for i in range(0, len(rna), 3):
        codon = rna[i : i + 3]
        amino_acid = codon_table.get(codon, "")
        if amino_acid == "*":
            break
        protein.append(amino_acid)

    return "".join(protein)


def protein_to_dna(protein: str) -> str:
    """
    Reverse-translate a protein sequence into a DNA sequence.

    Args:
        protein (str): Protein sequence consisting of single-letter amino acid codes.

    Returns:
        str: DNA sequence obtained by translating each amino acid to all possible corresponding codons.

    Raises:
        ValueError: If an invalid amino acid code is encountered that does not correspond to any codon.

    Example:
        >>> protein_to_dna("IDRIDR")
        'AUCGAUCGAUCG'
    """
    codon_table = {
        "F": ["UUU", "UUC"],
        "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "Y": ["UAU", "UAC"],
        "*": ["UAA", "UAG", "UGA"],
        "C": ["UGU", "UGC"],
        "W": ["UGG"],
        "P": ["CCU", "CCC", "CCA", "CCG"],
        "H": ["CAU", "CAC"],
        "Q": ["CAA", "CAG"],
        "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "I": ["AUU", "AUC", "AUA"],
        "M": ["AUG"],
        "T": ["ACU", "ACC", "ACA", "ACG"],
        "N": ["AAU", "AAC"],
        "K": ["AAA", "AAG"],
        "V": ["GUU", "GUC", "GUA", "GUG"],
        "A": ["GCU", "GCC", "GCA", "GCG"],
        "D": ["GAU", "GAC"],
        "E": ["GAA", "GAG"],
        "G": ["GGU", "GGC", "GGA", "GGG"],
    }

    dna_sequences = []

    for amino_acid in protein:
        possible_codons = codon_table.get(amino_acid, [])
        if not possible_codons:
            raise ValueError(f"Invalid amino acid code '{amino_acid}'")

        dna_sequences.append(possible_codons)

    # Generate all combinations of codons
    import itertools

    all_combinations = list(itertools.product(*dna_sequences))

    # Join each combination into a single DNA sequence
    dna_sequences = ["".join(combination) for combination in all_combinations]

    return " ".join(dna_sequences)
