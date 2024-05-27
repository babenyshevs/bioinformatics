from collections import defaultdict
from typing import List


def hamming_distance(str1: str, str2: str) -> int:
    """
    Compute the Hamming distance between two strings.

    The Hamming distance between two strings of equal length is the number
    of positions at which the corresponding symbols are different.

    Args:
        str1 (str): The first input string.
        str2 (str): The second input string.

    Returns:
        int: The Hamming distance between the two input strings.

    Raises:
        ValueError: If the input strings are of different lengths.

    Example:
        >>> string1 = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
        >>> string2 = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
        >>> hamming_distance(string1, string2)
        20
    """
    # Check if the lengths of the strings are equal
    if len(str1) != len(str2):
        raise ValueError("Strings must be of equal length")

    # Initialize the distance counter
    distance = 0

    # Iterate through each character in the strings and compare
    for char1, char2 in zip(str1, str2):
        # If characters are different, increment the distance counter
        if char1 != char2:
            distance += 1

    return distance


def skew(genome: str, i: int) -> int:
    """
    Calculate the skew for a given genome up to position i.

    Skew is defined as the difference between the total number of occurrences
    of G and the total number of occurrences of C in the first i nucleotides of Genome.

    Args:
        genome (str): The genome sequence.
        i (int): The position up to which to calculate the skew.

    Returns:
        int: The skew value up to position i in the genome.

    # Example usage:
        genome = "CATTCCAGTACTTCATGATGGCGTGAAGA"
        position = 15
        print("Skew at position", position, ":", skew(genome, position))
    """
    # Initialize the skew counter
    skew_value = 0

    # Iterate through the genome up to position i
    for j in range(i):
        # Increment the skew counter based on the nucleotide at position j
        if genome[j] == "G":
            skew_value += 1
        elif genome[j] == "C":
            skew_value -= 1

    return skew_value


def skew_values(genome: str) -> list[int]:
    """
    Compute all values of Skewi for a given genome.

    Args:
        genome (str): The genome sequence.

    Returns:
        list: A space-separated string of all Skewi values for i ranging from 0 to len(genome).
    """
    # Initialize the skew counter and the list to store skew values
    skew_value = 0
    skew_list = [0]  # Start with skew value 0 at position 0

    # Iterate through the genome
    for nucleotide in genome:
        # Update the skew value based on the current nucleotide
        if nucleotide == "G":
            skew_value += 1
        elif nucleotide == "C":
            skew_value -= 1
        skew_list.append(skew_value)

    return skew_list


def min_skew_indices(genome: str) -> list:
    """
    Return the index or list of indices of minimum values in the skew list.

    Args:
        genome (str): The genome sequence.

    Returns:
        list: A list of indices of minimum values in the skew list.
    """
    # Calculate skew values for the given genome
    skew_list = skew_values(genome)

    # Find the minimum skew value
    min_skew = min(skew_list)

    # Find indices of minimum skew values
    min_indices = [i for i, skew in enumerate(skew_list) if skew == min_skew]

    return min_indices


def approximate_pattern_matching(pattern: str, text: str, d: int) -> list[int]:
    """
    Find all starting positions where Pattern appears as a substring of Text with at most d mismatches.

    Args:
    - pattern (str): The pattern string to search for.
    - text (str): The text string to search within.
    - d (int): Maximum number of mismatches allowed.

    Returns:
    - list: list of integers specifying all starting positions where Pattern appears in Text with at most d mismatches.

    Example:
        pattern = "CTTGATCAT"
        text = "CTTGATCATCTTGATCATCTTGATCAT"
        d = 1
        positions = approximate_pattern_matching(pattern, text, d)

    # Test the function
        pattern = "CTTGATCAT"
        text = "CTTGATCATCTTGATCATCTTGATCAT"
        d = 1
        positions = approximate_pattern_matching(pattern, text, d)

    """
    positions = []

    pattern_len = len(pattern)
    text_len = len(text)

    # Iterate through the text string
    for i in range(text_len - pattern_len + 1):
        # Extract the substring of length pattern_len starting at position i
        substring = text[i : i + pattern_len]

        # Calculate the Hamming distance between the pattern and the extracted substring
        distance = hamming_distance(pattern, substring)

        # Check if the Hamming distance is less than or equal to d
        if distance <= d:
            positions.append(i)

    return positions


def reverse_complement(pattern: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        pattern (str): The DNA sequence.

    Returns:
        str: The reverse complement of the input DNA sequence.

    Example:
        >>> reverse_complement("ATGC")
        'GCAT'
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[base] for base in reversed(pattern))


def count_mismatches(pattern: str, text: str) -> int:
    """
    Count the number of mismatches between a pattern and a text.

    Args:
        pattern (str): The pattern to compare.
        text (str): The text to compare against.

    Returns:
        int: The count of mismatches between the pattern and the text.

    Example:
        >>> count_mismatches("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")
        7
    """
    count = 0
    for i in range(len(pattern)):
        if pattern[i] != text[i]:
            count += 1
    return count


def generate_neighbors(pattern: str, d: int) -> List[str]:
    """
    Generate all approximate neighbors of a pattern.

    Args:
        pattern (str): The pattern to generate neighbors for.
        d (int): Maximum number of mismatches allowed.

    Returns:
        List[str]: A list of all approximate neighbors of the pattern.

    Example:
        >>> generate_neighbors("ACGT", 1)
        ['ACGT', 'CCGT', 'GCGT', 'TCGT', 'ATGT', 'AGGT', 'ACCT', 'ACAT', 'ACGA', 'ACGC']
    """
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]
    neighborhood = []
    suffix_neighbors = generate_neighbors(pattern[1:], d)
    for neighbor in suffix_neighbors:
        if count_mismatches(pattern[1:], neighbor) < d:
            for nucleotide in ["A", "C", "G", "T"]:
                neighborhood.append(nucleotide + neighbor)
        else:
            neighborhood.append(pattern[0] + neighbor)
    return neighborhood


def freq_rc_kmers_in_neighbour(text: str, k: int, d: int) -> List[str]:
    """
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.

    Args:
        text (str): The DNA string to search in.
        k (int): Length of the k-mers.
        d (int): Maximum number of mismatches allowed.

    Returns:
        List[str]: A list of the most frequent k-mers with mismatches and their reverse complements.

    Example:
        >>> frequent_words_with_mismatches_and_reverse_complements("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
        ['ACAT', 'ATGT', 'ATGC', 'CATG', 'CATC', 'CTAC', 'CTAG', 'GATG', 'GATC', 'GCAT', 'GCAC', 'TACG', 'TCAT', 'TCAG', 'TGCA']
    """
    patterns = defaultdict(int)
    for i in range(len(text) - k + 1):
        pattern = text[i : i + k]
        neighborhood = generate_neighbors(pattern, d)
        for neighbor in neighborhood:
            patterns[neighbor] += 1
            rc_neighbor = reverse_complement(neighbor)
            patterns[rc_neighbor] += 1
    max_count = max(patterns.values())
    frequent_patterns = [pattern for pattern, count in patterns.items() if count == max_count]
    return frequent_patterns


def freq_kmers_in_neighbour(text: str, k: int, d: int) -> List[str]:
    """
    Find the most frequent k-mers with mismatches in a text.

    Args:
        text (str): The text to search in.
        k (int): Length of the k-mers.
        d (int): Maximum number of mismatches allowed.

    Returns:
        List[str]: A list of the most frequent k-mers with mismatches.

    Example:
        >>> frequent_words_with_mismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
        ['ATGT', 'GATG', 'ATGC']
    """
    patterns = defaultdict(int)
    for i in range(len(text) - k + 1):
        pattern = text[i : i + k]
        neighborhood = generate_neighbors(pattern, d)
        for neighbor in neighborhood:
            patterns[neighbor] += 1
    max_count = max(patterns.values())
    frequent_patterns = [pattern for pattern, count in patterns.items() if count == max_count]
    return frequent_patterns


def suffix(pattern: str) -> str:
    """
    Returns the suffix of a given pattern.

    Args:
        pattern (str): The input pattern.

    Returns:
        str: The suffix of the input pattern.

    Example:
        >>> suffix("ACGT")
        'CGT'
    """
    return pattern[1:]


def first_symbol(pattern: str) -> str:
    """
    Returns the first symbol of a given pattern.

    Args:
        pattern (str): The input pattern.

    Returns:
        str: The first symbol of the input pattern.

    Example:
        >>> first_symbol("ACGT")
        'A'
    """
    return pattern[0]


def neighbors(pattern: str, d: int) -> set:
    """
    Generates a set of closely related patterns based on a given pattern and a maximum distance.

    Args:
        pattern (str): The input pattern.
        d (int): The maximum distance allowed.

    Returns:
        set: A set of closely related patterns.

    Example:
        >>> neighbors("ACG", 1)
        {'ACT', 'AGG', 'ACA', 'ACG', 'ACC', 'TCG', 'GCG', 'CCG', 'ATG', 'AAG', 'AGA', 'AC T', 'AC C', 'AC A', 'AC G'}
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}

    neighborhood = set()
    suffix_neighbors = neighbors(suffix(pattern), d)
    for text in suffix_neighbors:
        if hamming_distance(suffix(pattern), text) < d:
            for nucleotide in ["A", "C", "G", "T"]:
                neighborhood.add(nucleotide + text)
        else:
            neighborhood.add(first_symbol(pattern) + text)

    return neighborhood
