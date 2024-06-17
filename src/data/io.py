from typing import List


def read_text_from_file(file_path: str) -> dict:
    """
    Reads text from a specified file and returns a dictionary.

    Parameters:
    file_path (str): The path to the text file.

    Returns:
    dict: A dictionary where each line of the text file is a key.

    Example:
    >>> text_dict = read_text_from_file("example.txt")
    >>> print(text_dict)
    {'param1': 'This is line 1', 'param2': 'This is line 2', ...}
    """
    text_dict = {}
    with open(file_path, "r") as file:
        lines = file.readlines()
        for idx, line in enumerate(lines, start=1):
            key = f"param_{idx}"
            text_dict[key] = line.strip()
    return text_dict


def print_list(lst):
    """
    Print the elements of a list as a string separated by space.

    Args:
    lst (list): The list of values to print.

    # Example usage:
        my_list = [1, 2, 3, 4, 5]
        print_list_as_string(my_list)
    """
    print(" ".join(map(str, lst)))


def string_to_list(input_string: str) -> List[str]:
    """
    Convert a given string into a list of strings, using a space as the separator.

    Args:
        input_string (str): The input string containing substrings separated by spaces.

    Returns:
        List[str]: A list of substrings separated by spaces from the input string.

    Example:
        >>> input_string = "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
        >>> string_to_list(input_string)
        ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    """
    return input_string.split()
