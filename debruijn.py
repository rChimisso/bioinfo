from igraph import Graph
import igraph as IGraph
import matplotlib.pyplot as plt
import random as rand


class Node:
  def __init__(self, label):
    self.label = label
    self.indegree = 0
    self.outdegree = 0

class DeBruijnGraph:
  def __init__(self):
    self.vertices = {}
    self.edges = {}

  def add_vertex(self, v):
    self.vertices[v] = Node(v)
    self.edges[v] = []
      
  def add_indegree(self, v, n = 1):
    self.vertices[v].indegree += n
  
  def add_outdegree(self, v, n = 1):
    self.vertices[v].outdegree += n
  
  def sub_indegree(self, v, n = 1):
    self.vertices[v].indegree -= n
  
  def sub_outdegree(self, v, n = 1):
    self.vertices[v].outdegree -= n
      
  def add_edge(self, v1, v2):
    if v1 not in self.vertices.keys():
      self.add_vertex(v1)
    self.add_outdegree(v1)
    
    if v2 not in self.vertices.keys():
      self.add_vertex(v2)
    self.add_indegree(v2)
    
    self.edges[v1] += [v2]
  
  def find_null_indegree_vertices(self):
    null_indegree_vertices = []

    for v in self.vertices.keys():
      if self.vertices[v].indegree == 0:
        null_indegree_vertices += [v]

    return null_indegree_vertices
  
  def find_null_outdegree_vertices(self):
    null_outdegree_vertices = []

    for v in self.vertices.keys():
      if self.vertices[v].outdegree == 0:
        null_outdegree_vertices += [v]
            
    return null_outdegree_vertices
  
  def build_contigs(self):
    starts = self.find_null_indegree_vertices()
    contigs = []

    for start in starts:
      contig = start
      current_vertex = start

      while len(self.edges[current_vertex]) > 0:
        next_vertex = self.edges[current_vertex][0]
        del self.edges[current_vertex][0]
        contig += next_vertex[-1]
        current_vertex = next_vertex

      print(contig)

      contigs += [contig]

    return contigs

def build_de_Bruijn_graph(reads, k):
  g = DeBruijnGraph()

  for read in reads:
    i = 0
    for i in range(len(read) - k):
      v1 = read[i:i+k]
      v2 = read[i+1:i+k+1]
      g.add_edge(v1, v2)

  print("Graph built")

  return g

def plot_de_Bruijn_graph(g: DeBruijnGraph):
  graph = Graph.ListDict(g.edges, directed=True, vertex_name_attr='label')
  layout = graph.layout("kk")
  #fig, ax = plt.subplots()
  #plt.bbox(3000,3000)
  IGraph.plot(
    graph,
    layout = layout,
    #vertex_size=1.3,
    vertex_color=['yellow'],
    #label_dist = 200,
    #edge_width=[2],
    #edge_color=['black'],
    #margin = [200],
    #target = ax,
    #edge_arrow_size = 0.02,
    #edge_arrow_width = 2
    bbox = (1000, 1000)
  )
