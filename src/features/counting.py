def pattern_count(Text: str, pattern: str) -> int:
    """
    Counts the number of times a pattern appears in a given text.

    Parameters:
    Text (str): The text in which to search for the pattern.
    pattern (str): The pattern to search for.

    Returns:
    int: The number of times the pattern appears in the text.

    Example:
    >>> text = "GCGCG"
    >>> pattern = "GCG"
    >>> pattern_count(text, pattern)
    2
    """
    count = 0
    text_length = len(Text)
    pattern_length = len(pattern)

    for i in range(text_length - pattern_length + 1):
        if Text[i : i + pattern_length] == pattern:
            count += 1

    return count


def frequency_table(text: str, k: int) -> dict[str, int]:
    """
    Generates a frequency table for all k-mers in the given text.

    Args:
    - text (str): The input text.
    - k (int): The length of k-mers to consider.

    Returns:
    - dict[str, int]: A dictionary mapping k-mers to their frequencies.
    """
    freq_map = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i : i + k]
        if pattern not in freq_map:
            freq_map[pattern] = 1
        else:
            freq_map[pattern] += 1
    return freq_map


def max_map(freq_map: dict[str, int]) -> int:
    """
    Finds the maximum value in the frequency map.

    Args:
    - freq_map (dict): A dictionary mapping strings to integers.

    Returns:
    - int: The maximum value in the frequency map.
    """
    if not freq_map:
        return 0
    return max(freq_map.values())


def frequent_words(text: str, k: int) -> list[str]:
    """
    Finds all frequent k-mers in the given text.

    Args:
    - text (str): The input text.
    - k (int): The length of k-mers to consider.

    Returns:
    - list[str]: An array of frequent k-mers.
    """
    frequent_patterns = []
    freq_map = frequency_table(text, k)
    max_freq = max_map(freq_map)
    for pattern, freq in freq_map.items():
        if freq == max_freq:
            frequent_patterns.append(pattern)
    return frequent_patterns


def reverse_complement(dna_string: str) -> str:
    """
    Find the reverse complement of a DNA string.

    Args:
    - dna_string (str): A DNA string.

    Returns:
    - str: The reverse complement of the input DNA string.

    Example:
        dna_string = "ATCGATCG"
        reverse_complement_string = reverse_complement(dna_string)
    """

    # Dictionary to store the complement of each nucleotide
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    # Initialize an empty string to store the reverse complement
    reverse_complement_string = ""

    # Iterate through the characters of the input DNA string in reverse order
    for nucleotide in reversed(dna_string):
        # Find the complement of each nucleotide and append it to the reverse complement string
        reverse_complement_string += complement[nucleotide]

    return reverse_complement_string


def pattern_matching(pattern: str, genome: str) -> str:
    """
    Find all starting positions where Pattern appears as a substring of Genome.

    Args:
    - pattern (str): The pattern string to search for.
    - genome (str): The genome string to search within.

    Returns:
    - str: A space-separated string of integers specifying all starting positions where Pattern appears in Genome.

    Example:
        pattern = "CTTGATCAT"
        genome = "CTTGATCATCTTGATCATCTTGATCAT"
        positions = pattern_matching(pattern, genome)
    """
    positions = []
    pattern_len = len(pattern)
    genome_len = len(genome)

    # Iterate through the genome string
    for i in range(genome_len - pattern_len + 1):
        # Check if the substring of length pattern_len starting at position i matches the pattern
        if genome[i : i + pattern_len] == pattern:
            positions.append(str(i))

    return " ".join(positions)


def find_clumps(text: str, k: int, L: int, t: int, allow_duplicates: bool = True) -> list:
    """
    Find clumps of k-mers in a genome where a k-mer appears at least t times
    within a window of length L.

    Args:
    - text (str): The genome text.
    - k (int): The length of k-mers.
    - L (int): The length of the window.
    - t (int): The minimum number of times a k-mer should appear in the window.
    - allow_duplicates (bool): Whether to allow counting a k-mer more than once.

    Returns:
    - list: A list of unique k-mers that form clumps in the genome.
    """

    def frequency_table(s: str, k: int) -> dict:
        freq_map = {}
        for i in range(len(s) - k + 1):
            kmer = s[i : i + k]
            freq_map[kmer] = freq_map.get(kmer, 0) + 1
        return freq_map

    patterns = set() if not allow_duplicates else []

    freq_map = frequency_table(text[:L], k)
    for pattern, count in freq_map.items():
        if count >= t:
            if allow_duplicates:
                patterns.add(pattern)
            else:
                patterns.add(pattern)

    for i in range(1, len(text) - L + 1):
        prev_kmer = text[i - 1 : i - 1 + k]
        new_kmer = text[i + L - k : i + L]
        freq_map[prev_kmer] -= 1
        freq_map[new_kmer] = freq_map.get(new_kmer, 0) + 1
        if freq_map[new_kmer] >= t:
            if allow_duplicates:
                patterns.add(new_kmer)
            elif new_kmer not in patterns:
                patterns.add(new_kmer)

    return list(patterns)
