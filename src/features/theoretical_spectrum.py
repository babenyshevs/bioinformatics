def amino_acid_mass_dict() -> dict:
    """
    Returns a dictionary mapping amino acid symbols to their integer masses.

    Example:
    >>> amino_acid_mass_dict()
    {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113,
     'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147,
     'R': 156, 'Y': 163, 'W': 186}
    """
    mass_table = """
    G 57
    A 71
    S 87
    P 97
    V 99
    T 101
    C 103
    I 113
    L 113
    N 114
    D 115
    K 128
    Q 128
    E 129
    M 131
    H 137
    F 147
    R 156
    Y 163
    W 186"""
    mass = mass_table.split()
    return {mass[i]: int(mass[i + 1]) for i in range(0, len(mass), 2)}


def linear_spectrum(peptide: str) -> list:
    """
    Generate the linear spectrum of a peptide.

    Args:
    - peptide: A string representing an amino acid sequence.

    Returns:
    - A sorted list of integers representing the linear spectrum.

    Example:
    >>> linear_spectrum("LEQN")
    [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
    """
    mass_dict = amino_acid_mass_dict()
    n = len(peptide)
    prefix_mass = [0]
    for i in range(n):
        prefix_mass.append(prefix_mass[i] + mass_dict[peptide[i]])
    l_spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            l_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(l_spectrum)


def cyclic_spectrum(peptide: str) -> list:
    """
    Generate the cyclic spectrum of a peptide.

    Args:
    - peptide: A string representing an amino acid sequence.

    Returns:
    - A sorted list of integers representing the cyclic spectrum.

    Example:
    >>> cyclic_spectrum("LEQN")
    [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
    """
    mass_dict = amino_acid_mass_dict()
    n = len(peptide)
    prefix_mass = [0]
    for i in range(n):
        prefix_mass.append(prefix_mass[i] + mass_dict[peptide[i]])
    peptide_mass = prefix_mass[n]
    c_spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            c_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < n:
                c_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(c_spectrum)


if __name__ == "__main__":
    peptide = input().strip()
    c_spectrum = linear_spectrum(peptide)
    print(" ".join(map(str, c_spectrum)))
