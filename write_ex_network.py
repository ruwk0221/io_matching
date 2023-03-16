import jReversion as jR
import networkx as nx
from LDOI import BooleanDOI_processing as BDOIp

# ModelList = [jR.toy_model(14)]
# ModelList = [jR.erbb(), jR.bauer(), jR.grieco(),
#              jR.prostate(), jR.colitis_jijoo(), jR.blt(), jR.saez(), jR.fibroblast()]
# ModelList = [jR.tlgl(), jR.emt()]
# ModelList = [jR.toy_model(x) for x in range(1, 9)]
ModelList =[{
    'BooleanRule_filename': 'models/Grieco_robust_br1.txt',
    'network_name': 'grieco_br1',
    # 'num_inputs':,
    # 'num_input_conditions':,
    # 'input_conditions':,
    # 'output_nodes':

}

]
for Model in ModelList:
    BooleanRuleFileName = Model['BooleanRule_filename']
    NetworkName = Model['network_name']
    #
    # NumInputs = Model['num_inputs']
    # NumInputConditions = Model['num_input_conditions']
    #
    # InputConditions = Model['input_conditions']
    #
    # OutputNodes = Model['output_nodes']

    # Set parameters
    # Note the node name for Gread is the index (integer), one can encode the nodename by adding prefix and suffix
    # If the node name from the input file is not this simple, one need to create a dictionary to record the mapping
    Prefix, Suffix = 'n', 'n'

    TempI = jR.initialize(BooleanRuleFileName, Prefix, Suffix)
    Mapping = TempI['mapping']
    InverseMapping = TempI['inverse_mapping']
    GRead = TempI['Gread']
    ReadNodes = TempI['read_nodes']

    try:
        GExpanded = nx.read_gml('networks/' + NetworkName + '_expanded_network.gml')
    except IOError:
        GExpanded = BDOIp.Get_expanded_network(GRead, prefix=Prefix, suffix=Suffix)
        nx.write_gml(GExpanded, 'networks/' + NetworkName + '_expanded_network.gml')

    Gtest = nx.DiGraph()
    Gtest.add_edges_from(GExpanded.edges())
    Gtest = nx.relabel_nodes(Gtest, InverseMapping)

    path = 'networks/' + NetworkName + '_ex_topology.gml'
    with open(path, 'wb') as f:
        nx.write_gml(G=Gtest, path=f)

