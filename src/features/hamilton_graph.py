from collections import defaultdict, deque
from typing import Dict, List, Tuple


class Graph:
    def __init__(self, adjacency_list: Dict[str, List[str]]):
        """
        Initializes the graph with the given adjacency list.

        :param adjacency_list: A dictionary where keys are vertex names and values are lists of adjacent vertices.

        Example usage:
        adjacency_list = {
            "A": ["B", "C"],
            "B": ["A", "C", "D"],
            "C": ["A", "B", "D"],
            "D": ["B", "C", "E"],
            "E": ["D"],
        }

        graph = Graph(adjacency_list)
        print("Adjacency List Representation of the Graph:")
        print(graph)
        print("\nHamiltonian Cycle:")
        hamiltonian_cycle = ["A", "B", "D", "E", "D", "C", "A"]
        graph.display_hamiltonian_cycle(hamiltonian_cycle)

        """
        self.adjacency_list = adjacency_list

    def display_hamiltonian_cycle(self, cycle: List[str]) -> None:
        """
        Displays the Hamiltonian cycle.

        :param cycle: A list of vertices representing the Hamiltonian cycle.
        """
        if self.is_valid_hamiltonian_cycle(cycle):
            print(" -> ".join(cycle))
        else:
            print("The provided cycle is not a valid Hamiltonian cycle.")

    def is_valid_hamiltonian_cycle(self, cycle: List[str]) -> bool:
        """
        Validates if the given cycle is a valid Hamiltonian cycle.

        :param cycle: A list of vertices representing the Hamiltonian cycle.
        :return: True if the cycle is a valid Hamiltonian cycle, False otherwise.
        """
        if len(cycle) != len(self.adjacency_list) + 1:
            return False

        visited = set()
        for i in range(len(cycle) - 1):
            if cycle[i] in visited or cycle[i + 1] not in self.adjacency_list[cycle[i]]:
                return False
            visited.add(cycle[i])

        return cycle[0] == cycle[-1] and len(visited) == len(self.adjacency_list)

    def __str__(self) -> str:
        """
        Provides a string representation of the graph.

        :return: A string representing the adjacency list of the graph.
        """
        return "\n".join(
            [
                f"{vertex}: {', '.join(neighbors)}"
                for vertex, neighbors in self.adjacency_list.items()
            ]
        )


####################################################################################################


def is_n_universal(binary_string: str, n: int = 3) -> bool:
    """
    Check if the given binary string is n-universal linear binary string.

    A binary string is n-universal if every possible n-bit binary string
    appears as a substring.

    Args:
        binary_string (str): The binary string to check.
        n (int): The length of substrings to check for universality. (defaults to 3)

    Returns:
        bool: True if the string is n-universal, False otherwise.

    Example:
        >>> is_n_universal("000111010011", 3)
        True
        >>> is_n_universal("000110", 3)
        False
    """
    # List of all possible n-bit binary strings
    required_substrings = {format(i, f"0{n}b") for i in range(2**n)}

    # Generate all n-bit substrings from the given binary string
    string_length = len(binary_string)
    substrings = {binary_string[i : i + n] for i in range(string_length - n + 1)}

    # Check if all required substrings are present
    return required_substrings.issubset(substrings)


####################################################################################################


def construct_de_bruijn_graph(
    kmers: List[str],
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Constructs the de Bruijn graph from a collection of k-mers.

    Args:
        kmers (List[str]): A list of k-mers.

    Returns:
        Tuple[Dict[str, List[str]], Dict[str, List[str]]]: A tuple containing two dictionaries:
            - adjacency_list: adjacency list of the de Bruijn graph.
            - incoming_edges: dictionary mapping nodes to their incoming edges.
    """
    adjacency_list = defaultdict(list)
    incoming_edges = defaultdict(list)

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_list[prefix].append(suffix)
        incoming_edges[suffix].append(prefix)

    return adjacency_list, incoming_edges


def count_indegree(node: str, incoming_edges: Dict[str, List[str]]) -> int:
    """
    Counts the indegree of a given node in the de Bruijn graph.

    Args:
        node (str): The node whose indegree is to be counted.
        incoming_edges (Dict[str, List[str]]): Dictionary mapping nodes to their incoming edges.

    Returns:
        int: The indegree of the node.

    # Example usage
    kmers = ["ATG", "TGA", "GAT", "ATC", "TCG", "CGA", "GAA", "AAT", "ATG", "TGC", "GCA"]

    # Step 1: Construct the de Bruijn graph
    adjacency_list, incoming_edges = construct_de_bruijn_graph(kmers)

    # Step 2: Count the indegree of "GTA"
    node = "GTA"
    indegree_gta = count_indegree(node, incoming_edges)

    print(f"The indegree of {node} is {indegree_gta}")
    """
    return len(incoming_edges[node])


####################################################################################################


def find_eulerian_path(adjacency_list: Dict[str, List[str]]) -> List[str]:
    """
    Finds an Eulerian path in a directed graph represented by an adjacency list.

    Args:
        adjacency_list (Dict[str, List[str]]): The adjacency list of the graph.

    Returns:
        List[str]: The Eulerian path as a list of nodes.
    """

    def balance(graph):
        """
        Returns the balance of each node: outdegree - indegree.
        """
        outdegree = defaultdict(int)
        indegree = defaultdict(int)

        for node in graph:
            outdegree[node] += len(graph[node])
            for neighbor in graph[node]:
                indegree[neighbor] += 1

        balance = {
            node: outdegree[node] - indegree[node] for node in set(outdegree) | set(indegree)
        }
        return balance

    balance_map = balance(adjacency_list)
    start_node = [node for node in balance_map if balance_map[node] == 1][0]

    # Hierholzer's algorithm to find the Eulerian path
    path = []
    stack = [start_node]
    current_path = []

    while stack:
        current = stack[-1]
        if adjacency_list[current]:
            next_node = adjacency_list[current].pop()
            stack.append(next_node)
        else:
            current_path.append(stack.pop())

    return current_path[::-1]


####################################################################################################


def reconstruct_string_from_kmers(kmers: List[str]) -> str:
    """
    Reconstructs a linear string from the composition of 4-mers.

    Args:
        kmers (List[str]): A list of 4-mers.

    Returns:
        str: The reconstructed string.

    # Example usage:
    kmers = ["ATGC", "TGCA", "GCAT", "CATG"]
    reconstructed_string = reconstruct_string_from_kmers(kmers)
    print(reconstructed_string)  # Output should be "ATGCATG"

    """
    adjacency_list, _ = construct_de_bruijn_graph(kmers)
    eulerian_path = find_eulerian_path(adjacency_list)

    # Reconstruct the string from the Eulerian path
    reconstructed_string = eulerian_path[0]
    for node in eulerian_path[1:]:
        reconstructed_string += node[-1]

    return reconstructed_string


####################################################################################################


def min_edges_to_balance_graph(adj_list: dict) -> int:
    """
    Calculate the minimum number of edges to add to make each node in the graph balanced.

    :param adj_list: Adjacency list of the graph as a dictionary where keys are nodes and
                     values are lists of adjacent nodes.
    :return: Minimum number of edges needed to balance the graph.

    Example:
    >>> adj_list = {
    ...     0: [1, 2],
    ...     1: [2],
    ...     2: [0],
    ...     3: [2, 1]
    ... }
    >>> min_edges_to_balance_graph(adj_list)
    2
    """
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    # Calculate in-degree and out-degree
    for node, neighbors in adj_list.items():
        out_degree[node] += len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1

    # Ensure all nodes are considered (even if they have no outgoing edges)
    all_nodes = set(adj_list.keys()).union(
        set(n for neighbors in adj_list.values() for n in neighbors)
    )

    imbalance = 0
    for node in all_nodes:
        out_minus_in = out_degree[node] - in_degree[node]
        if out_minus_in > 0:
            imbalance += out_minus_in

    return imbalance


####################################################################################################


# def parse_kmers(kmers: List[str]) -> List[Tuple[str, str]]:
#     """
#     Parse the (3,1)-mers into pairs of k-mers.

#     Args:
#         kmers (List[str]): List of (3,1)-mers in the format (ACC|ATA).

#     Returns:
#         List[Tuple[str, str]]: List of k-mer pairs.
#     """
#     return [tuple(kmer.split("|")) for kmer in kmers]


# def build_adjacency_list(pairs: List[Tuple[str, str]]) -> Tuple[dict, dict, dict]:
#     """
#     Build the adjacency list for the given k-mer pairs.

#     Args:
#         pairs (List[Tuple[str, str]]): List of k-mer pairs.

#     Returns:
#         Tuple[dict, dict, dict]: Adjacency list, in-degree dictionary, and out-degree dictionary.
#     """
#     adjacency_list = defaultdict(list)
#     in_degree = defaultdict(int)
#     out_degree = defaultdict(int)

#     for prefix, suffix in pairs:
#         adjacency_list[prefix].append(suffix)
#         out_degree[prefix] += 1
#         in_degree[suffix] += 1

#     return adjacency_list, in_degree, out_degree


# def find_start_node(adjacency_list: dict, in_degree: dict, out_degree: dict) -> str:
#     """
#     Find the start node for the Eulerian path.

#     Args:
#         adjacency_list (dict): The adjacency list.
#         in_degree (dict): In-degree count of nodes.
#         out_degree (dict): Out-degree count of nodes.

#     Returns:
#         str: The start node.
#     """
#     start = None
#     for node in adjacency_list:
#         if out_degree[node] > in_degree[node]:
#             return node
#         if out_degree[node] == in_degree[node]:
#             start = node
#     return start


# def eulerian_path(adjacency_list: dict, start_node: str) -> List[str]:
#     """
#     Find the Eulerian path in the given graph.

#     Args:
#         adjacency_list (dict): The adjacency list.
#         start_node (str): The starting node for the Eulerian path.

#     Returns:
#         List[str]: The Eulerian path.
#     """
#     path = []
#     stack = [start_node]
#     current_path = []

#     while stack:
#         current_node = stack[-1]
#         if adjacency_list[current_node]:
#             next_node = adjacency_list[current_node].pop()
#             stack.append(next_node)
#         else:
#             current_path.append(stack.pop())

#     return current_path[::-1]


# def reconstruct_string(path: List[str]) -> str:
#     """
#     Reconstruct the original string from the Eulerian path.

#     Args:
#         path (List[str]): The Eulerian path.

#     Returns:
#         str: The reconstructed string.
#     """
#     result = path[0]
#     for node in path[1:]:
#         result += node[-1]
#     return result


# def find_linear_string(kmers: List[str]) -> str:
#     """
#     Find the single linear string from the list of (3,1)-mers.

#     Args:
#         kmers (List[str]): List of (3,1)-mers in the format (ACC|ATA).

#     Returns:
#         str: The single linear string.

#     Example:
#         >>> find_linear_string(["ACC|CGA", "CGA|ATA", "ATA|TAA"])
#         'ACCGATAA'
#     """
#     pairs = parse_kmers(kmers)
#     adjacency_list, in_degree, out_degree = build_adjacency_list(pairs)
#     start_node = find_start_node(adjacency_list, in_degree, out_degree)
#     path = eulerian_path(adjacency_list, start_node)
#     return reconstruct_string(path)


# def reconstruct_genome_from_reads(reads: List[str]) -> str:
#     """
#     Reconstructs a single linear genome string from given (3,1)-mer composition reads.

#     Args:
#     - reads (List[str]): List of paired reads in the format "(3-mer|3-mer)".

#     Returns:
#     - str: Reconstructed genome string.

#     Example:
#     >>> reads = [
#     ...     "(ACC|ATA)",
#     ...     "(ACT|ATT)",
#     ...     "(ATA|TGA)",
#     ...     "(ATT|TGA)",
#     ...     "(CAC|GAT)",
#     ...     "(CCG|TAC)",
#     ...     "(CGA|ACT)",
#     ...     "(CTG|AGC)",
#     ...     "(CTG|TTC)",
#     ...     "(GAA|CTT)",
#     ...     "(GAT|CTG)",
#     ...     "(GAT|CTG)",
#     ...     "(TAC|GAT)",
#     ...     "(TCT|AAG)",
#     ...     "(TGA|GCT)",
#     ...     "(TGA|TCT)",
#     ...     "(TTC|GAA)"
#     ... ]
#     >>> reconstruct_genome_from_reads(reads)
#     'CACGATTGAGCTGAAGCTGAGCACCTGAGCTTCTGACTTGAACTGAGCTGATCTGATCTGAACTTGA'
#     """
#     # Initialize variables
#     read_pairs = [read.strip("()").split("|") for read in reads]
#     genome = read_pairs[0][0]  # Start with the first read of the first pair

#     # Iterate over the read pairs to reconstruct the genome
#     for i in range(len(read_pairs)):
#         genome += read_pairs[i][1][-1]  # Append the last character of the second read in each pair

#     return genome
