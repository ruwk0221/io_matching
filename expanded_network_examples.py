import networkx as nx
from LDOI import BooleanDOI_processing as BDOIp
import jReversion as jR

NetworkName = 'models/RC_mut'
# BooleanRuleFileName = 'Models/TLGLNetwork_jijoo.booleannet'
BooleanRuleFileName = NetworkName + '.txt'
Prefix, Suffix = 'n', 'n'

TempI = jR.initialize(BooleanRuleFileName, Prefix, Suffix)
Mapping = TempI['mapping']
InverseMapping = TempI['inverse_mapping']
GRead = TempI['Gread']
ReadNodes = TempI['read_nodes']
GExpanded = BDOIp.Get_expanded_network(GRead, prefix=Prefix, suffix=Suffix)

Gcopy = nx.DiGraph()
Gcopy.add_edges_from(GRead.edges())
Gcopy = nx.relabel_nodes(Gcopy, ReadNodes)
GExpanded = nx.relabel_nodes(GExpanded, InverseMapping)

# path = NetworkName + '_topology.gml'
# with open(path, 'wb') as f:
#     nx.write_gml(G=Gcopy, path=f)

path = NetworkName + '_expanded.gml'
with open(path, 'wb') as f:
    nx.write_gml(G=GExpanded, path=f)

