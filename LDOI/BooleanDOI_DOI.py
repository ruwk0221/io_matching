'''
This file contains functions to calculate the Logic Domain of Influence of a set of node states.
The algorithm is described in a paper titled "Target Control in Logical Models Using the Domain of Influence of Nodes"
to be publisehd in the Frontiers of Physiology.
The code is written by Gang Yang, Department of Physics, Penn State University.
'''
import networkx as nx
from collections import deque
from collections import defaultdict

def Child(G, node):
  '''
  Return all child nodes for a given node
  '''
  return [x[1] for x in G.edges() if x[0]==node]

def Parent(G,node):
  '''
  Return all parent nodes for a given node
  '''
  return [x[0] for x in G.edges() if x[1]==node]


def Negation(node):
  '''
  Return the complementary node of a given node in the expanded network generated by Jorge's program.
  That is take the negation of the given node.
  This does not apply to composite node.
  '''
  if '.' not in node:
    return str(-int(node))
  else:
    return str(-float(node))

def Negation_in_expanded(node):
  '''
  Return the complementary node of a given node in the expanded network.
  That is take the negation of the given node.
  This does not apply to composite node.
  '''
  if '~' in node:
    return node[1:]
  else:
    return '~'+node



def truncated_node_of_influence_BFS(G, source, composite_node_delimiter='_', Negation_func=Negation_in_expanded):
  '''
  Algorithm to find the Logic Domain of Influecen(LDOI) of a set of node states on expanded network
  The algorithm is a modified Breadth first search on the expanded network.
  The difference from a regular BFS is described as follows:
	  1.Include composite node only if all parents node are included
	  2.Include source node if revisited
	  3.Truncate descendants to avoid conflict against given intervention during BFS

  Parameters
  ----------
  G         : the given expanded network as an DiGraph
  source    : the list of node states that we are interested in to find LDOI
              source node can either be a single node or a list of nodes,
              DO NOT INCLUDE A NODE AND ITS NEGATION FOR SOURCE
  composite_node_delimiter='_' : as name suggested, used to parse a composite node
                                 e.g. a composite node 'A_B' means A AND B for the default value
  Negation_func= : the function to calculate the complementary node of a given node on the expanded network

  Returns
  -------
  LDOI in the format of a set,
  complementary list for the visited part of the search
  a list containing any found conflict name to the source during the search process
  a list containing all the composite nodes that is not included in the LDOI

  Notice LDOI does not include composite node in the formal definition,
  however, the algorithm returns it for completeness.

  '''
  queue=deque()  #frontier of the search
  visited=set()  #visited part of the search
  complist=[]    #complementary list for the visited part of the search
  waiting=[]     #waiting list for composite nodes
  conflict=[]    #any visited and truncated node that is the negation of any source node
  cpnode_count=defaultdict(int)    #a dict to record the visited time for a composite node
  #Modified BFS
  queue.extend(source)
  while len(queue)!=0:
    node=queue.popleft()
    if node not in visited and Negation_func(node) not in visited:
      visited.add(node)
      complist.append(node)
      children=Child(G,node)
      for child in children:
        if composite_node_delimiter not in child:
          queue.append(child)
        elif cpnode_count[child] == G.in_degree(child) - 1:
          queue.append(child)
        else:
          cpnode_count[child] += 1
    elif Negation_func(node) in visited:
      conflict.append(node)
    elif node in source:  #if the node is already the source, L won't include it, we demand CL to include it to see whether the source node appears again
      complist.append(node)
  assert visited==set(complist)
  # Modified return output wating->cpnode_count
  return set(complist[len(source):]),complist,conflict,cpnode_count


