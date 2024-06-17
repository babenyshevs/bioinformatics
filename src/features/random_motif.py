import random
from typing import Dict, List


def random_motif(dna: List[str], k: int, t: int) -> List[str]:
    """
    Generate a list of random k-mers, one from each DNA string in dna.

    Args:
        dna (List[str]): List of DNA strings.
        k (int): Length of k-mer.
        t (int): Number of DNA strings.

    Returns:
        List[str]: List of random k-mers.

    Example:
        >>> dna = ["AGCTAGC", "TTCGAGT", "GCTACGT"]
        >>> random_motif(dna, 3, 3)
        ['AGC', 'CGA', 'CTA']
    """
    random_motifs = []
    for string in dna:
        random_num = random.randint(0, len(string) - k)
        random_motifs.append(string[random_num : random_num + k])
    return random_motifs


def profile(motifs: List[str]) -> Dict[str, List[float]]:
    """
    Create a profile matrix with pseudocounts from a list of motifs.

    Args:
        motifs (List[str]): List of k-mers.

    Returns:
        Dict[str, List[float]]: Profile matrix with pseudocounts.

    Example:
        >>> motifs = ["AGC", "TGC", "ATC"]
        >>> profile(motifs)
        {'A': [0.5, 0.25, 0.25], 'C': [0.25, 0.5, 1.0], 'G': [0.25, 0.25, 0.25], 'T': [0.25, 0.25, 0.25]}
    """
    profile_pseudo = {}
    t = len(motifs)
    k = len(motifs[0])
    for symbol in "ACGT":
        profile_pseudo[symbol] = []
        for i in range(k):
            column = [motif[i] for motif in motifs]
            profile_pseudo[symbol].append((column.count(symbol) + 1) / (t + 4))
    return profile_pseudo


def profile_most_probable_kmer(text: str, k: int, profile: Dict[str, List[float]]) -> str:
    """
    Find the k-mer in a string that is most probable given a profile matrix.

    Args:
        text (str): DNA string.
        k (int): Length of k-mer.
        profile (Dict[str, List[float]]): Profile matrix.

    Returns:
        str: Most probable k-mer.

    Example:
        >>> text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
        >>> profile = {'A': [0.2, 0.4, 0.3, 0.1], 'C': [0.4, 0.3, 0.1, 0.2], 'G': [0.3, 0.1, 0.5, 0.4], 'T': [0.1, 0.2, 0.1, 0.3]}
        >>> profile_most_probable_kmer(text, 4, profile)
        'CCTA'
    """
    max_prob = -1
    best_kmer = text[0:k]
    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile[nucleotide][j]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer


def score(motifs: List[str]) -> int:
    """
    Calculate the score of a set of motifs.

    Args:
        motifs (List[str]): List of k-mers.

    Returns:
        int: Score of the motifs.

    Example:
        >>> motifs = ["AGC", "TGC", "ATC"]
        >>> score(motifs)
        3
    """
    total_score = 0
    k = len(motifs[0])
    t = len(motifs)
    for i in range(k):
        column = [motif[i] for motif in motifs]
        max_count = max(column.count(x) for x in set(column))
        total_score += t - max_count
    return total_score


def motifs(profile: Dict[str, List[float]], dna: List[str]) -> List[str]:
    """
    Generate the most probable motifs for each DNA string given a profile matrix.

    Args:
        profile (Dict[str, List[float]]): Profile matrix.
        dna (List[str]): List of DNA strings.

    Returns:
        List[str]: List of most probable k-mers.

    Example:
        >>> profile = {'A': [0.2, 0.4, 0.3, 0.1], 'C': [0.4, 0.3, 0.1, 0.2], 'G': [0.3, 0.1, 0.5, 0.4], 'T': [0.1, 0.2, 0.1, 0.3]}
        >>> dna = ["ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", "ACGCGGCTCTGAAAATAGCCGCCATACGGACATATCCCTTGTGTATCCAG"]
        >>> motifs(profile, dna)
        ['CCTA', 'GCGC']
    """
    t = len(dna)
    k = len(profile["A"])
    profile_most_probable_kmers = []
    for i in range(t):
        most_probable_kmer = profile_most_probable_kmer(dna[i], k, profile)
        profile_most_probable_kmers.append(most_probable_kmer)
    return profile_most_probable_kmers


def randomized_motif_search(dna: List[str], k: int, t: int) -> List[str]:
    """
    Perform the Randomized Motif Search algorithm to find the best motifs.

    Args:
        dna (List[str]): List of DNA strings.
        k (int): Length of k-mer.
        t (int): Number of DNA strings.

    Returns:
        List[str]: Best motifs found.

    Example:
        >>> dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
        >>> randomized_motif_search(dna, 8, 5)
        ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TCTCGGGG', 'TCTCGGGG']
    """
    motifs_random = random_motif(dna, k, t)
    best_motifs = motifs_random
    while True:
        for i in range(1000):
            current_profile = profile(motifs_random)
            motifs_random = motifs(current_profile, dna)
            if score(motifs_random) < score(best_motifs):
                best_motifs = motifs_random
            else:
                return best_motifs


def repeated_randomized_motif_search(dna: List[str], k: int, t: int) -> List[str]:
    """
    Repeat the Randomized Motif Search algorithm multiple times to find the best motifs.

    Args:
        dna (List[str]): List of DNA strings.
        k (int): Length of k-mer.
        t (int): Number of DNA strings.

    Returns:
        List[str]: Best motifs found over multiple runs.

    Example:
        >>> dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
        >>> repeated_randomized_motif_search(dna, 8, 5)
        ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TCTCGGGG', 'TCTCGGGG']
    """
    best_motifs_overall = []
    best_score_overall = float("inf")
    for i in range(1000):
        current_motifs = randomized_motif_search(dna, k, t)
        current_score = score(current_motifs)
        if current_score < best_score_overall:
            best_motifs_overall = current_motifs
            best_score_overall = current_score
    return best_motifs_overall


def gibbs_sampler(dna: List[str], k: int, t: int, n: int) -> List[str]:
    """
    Perform the Gibbs Sampling algorithm to find the best motifs.

    Args:
        dna (List[str]): List of DNA strings.
        k (int): Length of k-mer.
        t (int): Number of DNA strings.
        n (int): Number of iterations.

    Returns:
        List[str]: Best motifs found.

    Example:
        >>> dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
        >>> gibbs_sampler(dna, 8, 5, 100)
        ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TCTCGGGG', 'TCTCGGGG']
    """

    def weighted_random_choice(kmer_probs: Dict[str, float]) -> str:
        """
        Select a k-mer based on the probabilities from the profile.

        Args:
            kmer_probs (Dict[str, float]): Probabilities of k-mers.

        Returns:
            str: Chosen k-mer.
        """
        k_mers = list(kmer_probs.keys())
        probabilities = list(kmer_probs.values())
        return random.choices(k_mers, probabilities, k=1)[0]

    motifs = random_motif(dna, k, t)
    best_motifs = motifs[:]
    for _ in range(n):
        i = random.randint(0, t - 1)
        excluded_motifs = motifs[:i] + motifs[i + 1 :]
        current_profile = profile(excluded_motifs)

        # Calculate probabilities of k-mers in dna[i]
        kmer_probs = {}
        for j in range(len(dna[i]) - k + 1):
            kmer = dna[i][j : j + k]
            prob = 1
            for l, nucleotide in enumerate(kmer):
                prob *= current_profile[nucleotide][l]
            kmer_probs[kmer] = prob

        # Normalize probabilities
        total_prob = sum(kmer_probs.values())
        for kmer in kmer_probs:
            kmer_probs[kmer] /= total_prob

        motifs[i] = weighted_random_choice(kmer_probs)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]

    return best_motifs


def repeated_gibbs_sampler(dna: List[str], k: int, t: int, n: int) -> List[str]:
    """
    Repeat the Gibbs Sampler algorithm multiple times to find the best motifs.

    Args:
        dna (List[str]): List of DNA strings.
        k (int): Length of k-mer.
        t (int): Number of DNA strings.
        n (int): Number of iterations per run.

    Returns:
        List[str]: Best motifs found over multiple runs.

    Example:
        >>> dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
        >>> repeated_gibbs_sampler(dna, 8, 5, 100)
        ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TCTCGGGG', 'TCTCGGGG']
    """
    best_motifs_overall = None
    best_score_overall = float("inf")
    for _ in range(20):  # 20 random starts
        current_motifs = gibbs_sampler(dna, k, t, n)
        current_score = score(current_motifs)
        if current_score < best_score_overall:
            best_motifs_overall = current_motifs
            best_score_overall = current_score
    return best_motifs_overall
