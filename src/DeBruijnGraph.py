from typing import Union
from BioFileProcessor import BioFileProcessor
import log as log
import copy  
import numpy as np
import py4cytoscape as p4c
import networkx
import pandas as pd

class DeBruijnGraph:
  """
  De Bruijn Graph.
  
  Can take a list of sequences, create k-mers, add them to the graph, resolve bubbles statistically,
  and finally build contigs based on the De Bruijn graph using the Eulerian path algorithm.
  """
  
  def __init__(self, k: int, filepath: Union[str, None] = None, percentile: int = 0):
    """
    Construct an instance of DeBruijnGraph, initializing the adjacency list, in-degrees and out-degrees dictionaries, and k value.

    :param k: K-mer length.
    :type k: int
    :param filepath: (optional) Path to the .fasta/.fastq file to initialize the graph.
    :type filepath: str
    :param percentile: (optional, default 0) Minimum percentile of the edges to keep.
    :type percentile: int
    """
    self.adj_list: dict[str, dict[str, int]] = {}
    self.in_degrees: dict[str, int] = {}
    self.out_degrees: dict[str, int] = {}
    self.k: int = k
    log.initialize_file_log()

    if filepath is not None: 
      self.build_graph_from_file(filepath)
      self.delete_weak_edges(percentile)
      log.logger().info(f'Graph for file [{filepath}] initialized successfully')
    else:
      log.logger().info(f'Empty graph initialized successfully')
          
  def build_graph_from_file(self, filepath: str) -> None:
    """
    Builds the graph based on the reads contained in the given .fasta/.fastq file.

    :param filepath: The path to the .fasta/.fastq file.
    :type filepath: str
    """
    log.logger().info(f'Building the DeBruijnGraph for file [{filepath}]...')
    processor = BioFileProcessor(filepath)
    for read in processor.iter_filtered_reads(self.k):
      self.add_read(read)
    log.logger().info(f'DeBruijnGraph for file [{filepath}] built')

  def add_read(self, read: str) -> None:
    """
    Adds a read to the graph, building its k-mers.

    :param read: Read to add to the graph.
    :type read: str
    """
    self.add_kmers_to_graph(self.create_kmers(read))

  def create_kmers(self, read: str) -> list[str]:
    """
    Builds k-mers for a given read.

    :param read: Read to divide into k-mers.
    :type read: str
    :return: List of k-mers.
    :rtype: list[str]
    """
    kmers: list[str] = []
    for i in range(len(read) - self.k + 1):
      kmers.append(read[i : i + self.k])
    return kmers

  def add_kmers_to_graph(self, kmers: list[str]) -> None:
    """
    Add a list of k-mers to the graph.

    :param kmers: The k-mers to add to the graph.
    :type kmers: list[str]
    """
    for i in range(len(kmers) - 1):
      prefix = kmers[i]
      suffix = kmers[i + 1]
      # Update in-degrees for suffix node
      if suffix not in self.in_degrees:
        self.in_degrees[suffix] = 0
      self.in_degrees[suffix] += 1
      # Update out-degrees for prefix node
      if prefix not in self.out_degrees:
        self.out_degrees[prefix] = 0
      self.out_degrees[prefix] += 1
      # Update adjacency list for prefix node
      if prefix in self.adj_list:
        if suffix in self.adj_list[prefix]:
          self.adj_list[prefix][suffix] += 1
        else:
          self.adj_list[prefix][suffix] = 1
      else:
        self.adj_list[prefix] = {suffix: 1}

  def delete_weak_edges(self, percentile: int) -> None:
    """
    Remove edges with coverage less than a given percentile of the coverage of all edges.

    :param percentile: Minimum percentile of the edges to keep.
    :type percentile: int
    """
    log.logger().info('Deleting weak edges...')
    # Collect all coverages in a list
    all_coverages = [coverage for edges in self.adj_list.values() for coverage in edges.values()]
    # Compute the coverage percentile
    coverage_threshold = np.percentile(all_coverages, percentile)
    for node, edges in self.adj_list.items():
      for edge, coverage in edges.items():
        if coverage < coverage_threshold:
          edges.pop(edge)
          # Update in-degrees for edge
          self.in_degrees[edge] -= 1
          if self.in_degrees[edge] == 0:
            self.in_degrees.pop(edge)
          # Update out-degrees for node
          self.out_degrees[node] -= 1
          if self.out_degrees[node] == 0:
            self.out_degrees.pop(node)
      if not edges:
        self.adj_list.pop(node)
    log.logger().info('Weak edges deleted successfully')

  def get_start_nodes(self) -> list[str]:
    """
    Finds the starting nodes in the graph.

    :return: A list of starting nodes.
    :rtype: list[str]
    """
    log.logger().info('Finding starting nodes...')
    start_nodes: list[str] = []
    for node in self.adj_list:
      if self.in_degrees.get(node, 0) < self.out_degrees.get(node, 0):
        start_nodes.append(node)
    log.logger().info(f'Starting nodes found: {len(start_nodes)}')
    return start_nodes

  def get_eulerian_paths(self) -> tuple[dict[str, list[str]], dict[str, str]]:
    """
    Finds the Eulerian path for each starting node.

    :return: Two dictionaries containing the starting nodes as keys and their corresponding longest paths and sequences as values.
    :rtype: tuple[dict[str, list[str]], dict[str, str]]
    """
    start_nodes = self.get_start_nodes()
    print(len(start_nodes))
    if start_nodes:
      eulerian_paths: dict[str, list[str]] = {start_node: [] for start_node in start_nodes}
      eulerian_sequences: dict[str, str] = {start_node: "" for start_node in start_nodes}
      # Avoid modifying the original graph
      unused_edges = copy.deepcopy(self.adj_list)
      log.logger().info(f'Starting the Eulerian trails research...')
      for start_node in start_nodes:
        # Initialize cycle
        cycle = [start_node]
        node = start_node
        while True:
          # Ensure node exists in unused_edges
          while unused_edges.get(node, []):
            # Use popitem to remove an arbitrary (key, value) pair
            next_node, _ = unused_edges[node].popitem()
            cycle.append(next_node)
            node = next_node
          # If there are unused edges, choose a new start node
          if any(unused_edges.values()):
            unused_nodes = [n for n in cycle if unused_edges.get(n)]
            if not unused_nodes:
              break
            node = unused_nodes[0]
          else:
            # No unused edges are left, break the loop
            break
        # Assemble sequence from path
        sequence = cycle[0] + ''.join(kmer[-1] for kmer in cycle[1:])
        # Update path for start node
        eulerian_paths[start_node] = cycle
        # Update sequence for start node
        eulerian_sequences[start_node] = sequence
        log.logger().info(f'Found path of {len(sequence)} bases. Continuing the research...')
      log_sequences = '\n'.join(f'[{len(sequence)}]: sequence' for sequence in eulerian_sequences.values())
      log.logger().info('Eulerian trails research finished.')
      log.logger().info(f'Sequences found: {log_sequences}')
      return eulerian_paths, eulerian_sequences
    return {}, {}

  def visualize_with_cytoscape(self, show_paths: bool = True) -> None:
    """
    Visualize the graph with py4cytoscape. Eulerian paths are highlighted if show_paths is true.

    :param show_paths: (optional, default: True) Enable path highlighting.
    :type show_paths: bool
    """
    edges = {
      'sources': [],
      'targets': [],
      'weights': []
    }
    for source, targets in self.adj_list.items():
      for target, weight in targets.items():
        if target in self.adj_list.keys():
          edges['sources'].append(source)
          edges['targets'].append(target)
          edges['weights'].append(weight)

    edges = pd.DataFrame({'source': edges['sources'], 'target': edges['targets'], 'interaction': ['connected with'], 'weight': edges['weights']})
    nodes = pd.DataFrame({'id': list(self.adj_list.keys())})
    # Check if all sources are in targets and if all targets are in sources
    # nodes = list(self.adj_list.keys())
    # graph = networkx.Graph()
    # graph.add_nodes_from(nodes)
    # graph.add_edges_from(edges)
    # Create a new network in Cytoscape from our data
    # p4c.create_network_from_networkx(graph, 'De Bruijn Graph')
    log.logger().info('Creating Cytoscape network from dataframes...')
    p4c.create_network_from_data_frames(nodes, edges, title="De Bruijn Graph", collection="De Bruijn Graph Collection")
    p4c.layout_network('force-directed')
    p4c.create_view()
    log.logger().info('Cytoscape network initialized.')
    if show_paths:
      log.logger().info(f'Setting custom path appearance...')
      eulerian_paths, eulerian_sequences = self.get_eulerian_paths()
      # Apply the Eulerian path style to the appropriate edges
      for path in eulerian_paths.values():
        for _ in range(len(path) - 1):
          edge_name = ''
          p4c.set_edge_line_style_default('SOLID', edge_name)
          p4c.set_edge_line_width_default(2, edge_name)
          p4c.set_edge_opacity_default(200, edge_name)
          p4c.set_edge_color_default('#FF0000', edge_name)
      log.logger().info(f'Path appearance set')
    log.logger().info('Cytoscape network should be available on the client')
