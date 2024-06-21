from collections import Counter
from typing import List

# def calculate_matching_score(
#     theoretical_spectrum: List[float], experimental_spectrum: List[float]
# ) -> float:
#     """
#     Calculates a matching score between a theoretical spectrum and an experimental spectrum.

#     Parameters:
#     - theoretical_spectrum (List[float]): Theoretical spectrum as a list of peak masses.
#     - experimental_spectrum (List[float]): Experimental spectrum as a list of peak masses.

#     Returns:
#     - float: Matching score between 0 and 1.

#     Example:
#     >>> theoretical_spectrum = [100.0, 200.0, 300.0]
#     >>> experimental_spectrum = [100.0, 300.0]
#     >>> score = calculate_matching_score(theoretical_spectrum, experimental_spectrum)
#     >>> print(score)
#     0.6666666666666666
#     """
#     # Define a matching criterion, e.g., mass tolerance
#     mass_tolerance = 0.1  # hypothetical mass tolerance

#     # Initialize variables for score calculation
#     matched_peaks = 0
#     total_intensity_theoretical = sum(theoretical_spectrum)
#     total_intensity_experimental = sum(experimental_spectrum)

#     # Loop through experimental spectrum peaks
#     for exp_peak in experimental_spectrum:
#         # Check for a matching peak in theoretical spectrum
#         for theo_peak in theoretical_spectrum:
#             if abs(theo_peak - exp_peak) <= mass_tolerance:
#                 matched_peaks += 1
#                 break  # Matched, move to next experimental peak

#     # Calculate score
#     score = matched_peaks / len(experimental_spectrum)  # Simple matching ratio

#     return score


# def calculate_matching_count(
#     theoretical_spectrum: List[float], experimental_spectrum: List[float]
# ) -> int:
#     """
#     Calculates the count of matching peaks between a theoretical spectrum and an experimental spectrum.

#     Parameters:
#     - theoretical_spectrum (List[float]): Theoretical spectrum as a list of peak masses.
#     - experimental_spectrum (List[float]): Experimental spectrum as a list of peak masses.

#     Returns:
#     - int: Count of matching peaks between the theoretical and experimental spectra.

#     Example:
#     >>> theoretical_spectrum = [100.0, 200.0, 300.0]
#     >>> experimental_spectrum = [100.0, 300.0]
#     >>> count = calculate_matching_count(theoretical_spectrum, experimental_spectrum)
#     >>> print(count)
#     2
#     """
#     # Convert lists to sets for efficient intersection
#     theoretical_set = set(theoretical_spectrum)
#     experimental_set = set(experimental_spectrum)

#     # Find the intersection to get matching peaks
#     matching_peaks = theoretical_set.intersection(experimental_set)

#     # Count the number of matching peaks
#     matching_count = len(matching_peaks)

#     return matching_count


###################################################################################################


def find_max_multiplicity_element(spectrum):
    """
    Finds the element with the largest multiplicity in the spectral convolution of a given spectrum.

    Parameters:
    - spectrum (list of float): The input spectrum as a list of mass values.

    Returns:
    - float: The element with the largest multiplicity in the spectral convolution.

    Example:
    >>> spectrum = [100.0, 200.0, 300.0]
    >>> result = find_max_multiplicity_element(spectrum)
    >>> print(result)
    200.0
    """
    # Calculate spectral convolution
    convoluted_spectrum = []
    n = len(spectrum)

    for i in range(n):
        for j in range(i + 1, n):
            convoluted_spectrum.append(abs(spectrum[i] - spectrum[j]))

    # Count frequencies of each element in the spectral convolution
    freq_counter = Counter(convoluted_spectrum)

    # Find the element with the maximum frequency (multiplicity)
    max_element = max(freq_counter, key=freq_counter.get)

    return max_element


####################################################################################################


import sys
from collections import Counter

# filename = sys.argv[1]

# with open(filename) as file:
#     data = []
#     for line in file:
#         data.append(line)
# peptide = data[0][:-1]
# exp_spec = map(int, data[1].split())


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


def linear_spectrum(peptide, aa_mass):
    aa_mass = amino_acid_mass_dict()
    prefix_mass = [0 for x in range(len(peptide) + 1)]
    for i in range(1, len(peptide) + 1):
        prefix_mass[i] = prefix_mass[i - 1] + int(aa_mass[str(peptide[i - 1])])
    linear = []
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linear.append(prefix_mass[j] - prefix_mass[i])
    linear.append(0)
    return sorted(linear)


def cyclic_spectrum(peptide, aa_mass):
    aa_mass = amino_acid_mass_dict()
    prefix_mass = [0 for x in range(len(peptide) + 1)]
    for i in range(1, len(peptide) + 1):
        prefix_mass[i] = prefix_mass[i - 1] + int(aa_mass[str(peptide[i - 1])])
    cyclic = []
    peptide_mass = prefix_mass[-1]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            cyclic.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cyclic.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    cyclic.append(0)
    return sorted(cyclic)


def scoring(peptide, exp_spec, type="linear"):
    aa_mass = amino_acid_mass_dict()
    if type == "cyclic":
        theo_spec = cyclic_spectrum(peptide, aa_mass)
    else:
        theo_spec = linear_spectrum(peptide, aa_mass)
    theo_count = Counter(theo_spec)
    exp_count = Counter(exp_spec)
    count = 0
    for mass in exp_count.keys():
        if mass in theo_count.keys():
            count += min(exp_count[mass], theo_count[mass])
    return count
