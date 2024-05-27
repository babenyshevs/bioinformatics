import math
from typing import Dict, List, Set, Tuple

from src.features.replication_origin import hamming_distance, neighbors


def motif_enumeration(dna: List[str], k: int, d: int) -> Set[str]:
    """
    Find all (k, d)-motifs in a collection of strings.

    Args:
    - dna (List[str]): A collection of strings.
    - k (int): Length of the motifs.
    - d (int): Maximum number of mismatches allowed.

    Returns:
    - Set[str]: Set of all (k, d)-motifs found in DNA.

    Example:
    >>> dna = ["ATCGTAGC", "AAGCTAGC", "GTAGCTAC"]
    >>> k = 4
    >>> d = 1
    >>> motif_enumeration(Dna, k, d)
    {'GTAC', 'ATAG', 'ATGC', 'GTAG', 'CTAC', 'ATAC', 'CTAG', 'ATCG', 'GCTA', 'ATCA', 'ATGA'}
    """
    patterns = set()
    for i in range(len(dna[0]) - k + 1):
        pattern = dna[0][i : i + k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            if all(
                any(
                    hamming_distance(neighbor, dna_string[j : j + k]) <= d
                    for j in range(len(dna_string) - k + 1)
                )
                for dna_string in dna
            ):
                patterns.add(neighbor)
    return patterns


def count_motif_percent(motifs: List[str]) -> Dict[int, Dict[str, float]]:
    """
    Counts the percentage of each nucleotide ('A', 'C', 'G', 'T') at each position in the given list of motifs.

    Args:
        motifs (List[str]): A list of DNA motifs, where each motif is a string consisting of nucleotides.

    Returns:
        Dict[int, Dict[str, float]]: A dictionary where the keys represent the position in the motifs, and the values are dictionaries containing the percentage of each nucleotide at that position.

    Example:
        >>> motifs = ['ATG', 'CCG', 'AAG']
        >>> count_motif_percent(motifs)
        {
            0: {'A': 0.6666666666666666, 'C': 0.3333333333333333, 'G': 0.0, 'T': 0.0},
            1: {'A': 0.0, 'C': 0.6666666666666666, 'G': 0.3333333333333333, 'T': 0.0},
            2: {'A': 0.6666666666666666, 'C': 0.0, 'G': 0.0, 'T': 0.3333333333333333}
        }
    """
    count = {}
    columns = []
    for i in range(len(motifs[0])):
        columns.append([motif[i] for motif in motifs])
    for i in range(len(columns)):
        count[i] = {
            "A": columns[i].count("A") / len(columns[i]),
            "C": columns[i].count("C") / len(columns[i]),
            "G": columns[i].count("G") / len(columns[i]),
            "T": columns[i].count("T") / len(columns[i]),
        }

    return count


def motif_entropy(motifs: List[str]) -> float:
    """
    Calculate the motif entropy of a list of motifs.

    Args:
    motifs (List[str]): A list of DNA motifs represented as strings.

    Returns:
    float: The calculated motif entropy.

    Example:
    >>> motifs = ['ATG', 'GTC', 'CCA']
    >>> motif_entropy(motifs)
    1.584962500721156
    """
    entropy = 0
    percents = count_motif_percent(motifs)
    for i in range(len(percents)):
        for nucleotide in percents[i]:
            if percents[i][nucleotide] != 0:
                entropy += percents[i][nucleotide] * math.log2(percents[i][nucleotide])
    return -entropy


def median_string(k: int, dna: List[str]) -> str:
    """
    Finds the median string of length k that minimizes the total Hamming distance to a set of DNA sequences.

    Args:
    - k (int): The length of the desired median string.
    - dna (List[str]): A list of DNA sequences as strings.

    Returns:
    - str: The median string of length k that minimizes the total Hamming distance to the input DNA sequences.

    Example:
    >>> dna_sequences = [
    ...     "TCGGGGTTTTT",
    ...     "CCGGTGACTTAC",
    ...     "ACGGGGATTTTC",
    ...     "TTGGGGACTTTT",
    ...     "AAGGGGACTTCC",
    ...     "TTGGGGACTTCC",
    ...     "TCGGGGATTCAT",
    ...     "TCGGGGATTCCT",
    ...     "TAGGGGAACTAC",
    ...     "TCGGGTATAACC"
    ... ]
    >>> median_string(3, dna_sequences)
    'GGG'
    """
    distance = float("inf")
    pattern = ""

    for i in range(len(dna[0]) - k + 1):
        current_pattern = dna[0][i : i + k]
        current_distance = sum(
            min(hamming_distance(current_pattern, seq[i : i + k]) for seq in dna)
            for i in range(len(dna[0]) - k + 1)
        )

        if current_distance < distance:
            distance = current_distance
            pattern = current_pattern

    return pattern


def most_probable_kmer(text: str, k: int, profile: List[List[float]]) -> str:
    """
    Finds the Profile-most probable k-mer in a given text.

    Args:
        text: A string representing the input DNA sequence.
        k: An integer specifying the length of the k-mer.
        profile: A 4 × k matrix representing the profile of the DNA sequence.

    Returns:
        A Profile-most probable k-mer in the text.

    Example:
        >>> text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
        >>> k = 5
        >>> profile = [
        ...     [0.2, 0.2, 0.3, 0.2, 0.3],
        ...     [0.4, 0.3, 0.1, 0.5, 0.1],
        ...     [0.3, 0.3, 0.5, 0.2, 0.4],
        ...     [0.1, 0.2, 0.1, 0.1, 0.2]
        ... ]
        >>> profile_most_probable_kmer(text, k, profile)
        'CCGAG'
    """
    max_probability = -1
    most_probable_kmer = ""

    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        probability = 1
        for j in range(k):
            if kmer[j] == "A":
                probability *= profile[0][j]
            elif kmer[j] == "C":
                probability *= profile[1][j]
            elif kmer[j] == "G":
                probability *= profile[2][j]
            elif kmer[j] == "T":
                probability *= profile[3][j]
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer


def profile_most_probable_kmer(text: str, k: int, profile: dict) -> str:
    """
    Find the most probable k-mer in a DNA sequence based on a given profile.

    Args:
    - text (str): A DNA sequence.
    - k (int): The length of the k-mer.
    - profile (dict): A dictionary representing the profile of nucleotides.

    Returns:
    - str: The most probable k-mer in the given sequence.

    Example:
    >>> text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    >>> k = 5
    >>> profile = {'A': [0.2, 0.2, 0.3, 0.2, 0.3],
                  'C': [0.4, 0.3, 0.1, 0.5, 0.1],
                  'G': [0.3, 0.3, 0.5, 0.2, 0.4],
                  'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
    >>> profile_most_probable_kmer(text, k, profile)
    'CCGAG'
    """
    max_probability = -1.0
    most_probable_kmer = ""

    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        probability = 1.0

        for j, nucleotide in enumerate(kmer):
            probability *= profile[nucleotide][j]

        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer


def compute_pr(sequence: str, profile: Dict[str, List[float]]) -> float:
    """
    Compute the probability of a given sequence according to the profile matrix.

    Args:
    sequence (str): The sequence for which probability needs to be computed.
    profile (Dict[str, List[float]]): The profile matrix containing probabilities for each nucleotide.

    Returns:
    float: The probability of the sequence according to the profile matrix.

    # Example usage
        Profile = {
            "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
            "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
            "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
            "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0],
        }

        sequence = "AAGTTC"
        result = compute_pr(sequence, Profile)
    """
    probability = 1.0
    for i, nucleotide in enumerate(sequence):
        probability *= profile[nucleotide][i]
    return probability


def form_profile(motifs: List[str]) -> dict:
    """
    Form a profile matrix based on a list of DNA motifs.

    Args:
    - motifs (List[str]): A list of DNA motifs.

    Returns:
    - dict: A dictionary representing the profile matrix.

    Example:
    >>> motifs = ["ACGT", "ACGT", "ACGT"]
    >>> form_profile(motifs)
    {'A': [1.0, 1.0, 1.0, 1.0], 'C': [1.0, 1.0, 1.0, 1.0], 'G': [1.0, 1.0, 1.0, 1.0], 'T': [1.0, 1.0, 1.0, 1.0]}
    """
    k = len(motifs[0])
    profile = {nucleotide: [0] * k for nucleotide in "ACGT"}
    t = len(motifs)

    for i in range(k):
        column = [motif[i] for motif in motifs]
        for nucleotide in "ACGT":
            count = column.count(nucleotide)
            profile[nucleotide][i] = count / t

    return profile


def score_motifs(motifs: List[str]) -> int:
    """
    Calculate the score of a list of DNA motifs.

    Args:
    - motifs (List[str]): A list of DNA motifs.

    Returns:
    - int: The score of the motifs.

    Example:
    >>> motifs = ["ACGT", "ACGT", "ACGT"]
    >>> score_motifs(motifs)
    0
    """
    k = len(motifs[0])
    t = len(motifs)
    score = 0

    for i in range(k):
        column = [motif[i] for motif in motifs]
        max_count = max(column.count(nucleotide) for nucleotide in "ACGT")
        score += t - max_count

    return score


def greedy_motif_search(dna: List[str], k: int, t: int) -> List[str]:
    """
    Perform greedy motif search to find the best motifs from the given DNA sequences.

    Args:
    - dna (List[str]): A list of DNA sequences.
    - k (int): The length of the motif.
    - t (int): The number of DNA sequences.

    Returns:
    - List[str]: The best motifs found by greedy motif search algorithm.

    Example:
    ```python
    # Example usage:
    dna_sequences = [
        "GGCGTTCAGGCA",
        "AAGAATCAGTCA",
        "CAAGGAGTTCGC",
        "CACGTCAATCAC",
        "CAATAATATTCG"
    ]
    k = 3
    t = 5
    best_motifs = greedy_motif_search(dna_sequences, k, t)
    print(best_motifs)
    # Output:
    # ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
    ```
    """
    best_motifs = [seq[:k] for seq in dna]

    for i in range(len(dna[0]) - k + 1):
        motif1 = dna[0][i : i + k]
        motifs = [motif1]

        for j in range(1, t):
            profile = form_profile(motifs)
            motif_i = profile_most_probable_kmer(dna[j], k, profile)
            motifs.append(motif_i)

        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs

    return best_motifs


def profile_matrix_with_pseudocounts(motifs: List[str]) -> Dict[str, List[float]]:
    """
    Construct a profile matrix with pseudocounts based on a list of motifs.

    Args:
        motifs (List[str]): A list of DNA motifs.

    Returns:
        Dict[str, List[float]]: The profile matrix with pseudocounts, represented as a dictionary
        where each key is a nucleotide ('A', 'C', 'G', 'T') and each value is a list of probabilities
        corresponding to each position in the motifs.

    Example:
        >>> motifs = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
        >>> profile_matrix_with_pseudocounts(motifs)
        {'A': [0.2, 0.2, 0.0, 0.6, 0.2, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2],
         'C': [0.4, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.2, 0.4, 0.2, 0.2, 0.4],
         'G': [0.2, 0.8, 1.0, 0.2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0],
         'T': [0.2, 0.0, 0.0, 0.2, 0.6, 0.0, 1.0, 0.6, 0.4, 0.4, 0.6, 0.4]}
    """

    k = len(motifs[0])
    profile = {"A": [1] * k, "C": [1] * k, "G": [1] * k, "T": [1] * k}
    total_motifs = len(motifs) + 4  # pseudocounts for each nucleotide

    for i in range(k):
        for motif in motifs:
            profile[motif[i]][i] += 1

    for nucleotide in profile:
        for i in range(k):
            profile[nucleotide][i] /= total_motifs

    return profile


def greedy_motif_search_with_pseudocounts(dna: List[str], k: int, t: int) -> List[str]:
    """
    Implement a greedy motif search algorithm with pseudocounts.

    Args:
        dna (List[str]): A list of DNA sequences.
        k (int): The length of the motifs to search for.
        t (int): The number of DNA sequences in the input.

    Returns:
        List[str]: The best motifs found by the greedy motif search algorithm.

    Example:
        >>> dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
        >>> k = 3
        >>> t = 5
        >>> greedy_motif_search_with_pseudocounts(dna, k, t)
        ['TTC', 'ATC', 'TTC', 'TTC', 'TTC']
    """

    best_motifs = [string[:k] for string in dna]

    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i : i + k]]

        for j in range(1, t):
            profile = profile_matrix_with_pseudocounts(motifs)
            motifs.append(profile_most_probable_kmer(dna[j], k, profile))

        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs

    return best_motifs


def distance_between_pattern_and_strings(Pattern: str, Dna: List[str]) -> int:
    """
    Computes the total Hamming distance between a pattern and a collection of strings.

    Args:
        Pattern (str): The pattern to compare against.
        Dna (List[str]): A collection of strings.

    Returns:
        int: The sum of distances between the pattern and each string in Dna.

    Example:
        >>> distance_between_pattern_and_strings("AAA", ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG"])
        5
    """
    k = len(Pattern)
    distance = 0
    for Text in Dna:
        HammingDistance = float("inf")
        for i in range(len(Text) - k + 1):
            Pattern_prime = Text[i : i + k]
            if HammingDistance > hamming_distance(Pattern, Pattern_prime):
                HammingDistance = hamming_distance(Pattern, Pattern_prime)
        distance += HammingDistance
    return distance
