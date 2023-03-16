import jReversion as jR
from LDOI import BooleanDOI_processing as BDOIp
import networkx as nx
import par_helper as ph
import sys

models = ['bortezomib',
    # 'igvh',
    'apoptosis',
    # 'aurora',
    'bt474_long',
    'bt474_short',
    # 'cd4t',
    'colitis',
    'death',
    # 'egfr',
    # 'erbb',
    # 'fa_brca',
    # 'fa_check',
    'hcc1954_long',
    'hcc1954_short',
    'hgf',
    'mammalian',
    # 'mammalian_2006',
    'mapk',
    'oxidative',
    # 'pro_inflammatory',
    'fibroblasts',
    'skbr3_long',
    'skbr3_short',
    'tlgl_2008',
    'tlgl_2011',
    # 'tlgl_2011_reduced',
    # 'prostate',
    'migration']
num_worker = 5

for model in models:
    print(model)
    Prefix, Suffix = 'n', 'n'
    TEMP = jR.cellcollective(model, Prefix, Suffix)

    BooleanRuleFileName = TEMP['BooleanRule_filename']
    NetworkName = TEMP['network_name']

    NumInputs = TEMP['num_inputs']
    NumInputConditions = TEMP['num_input_conditions']

    InputConditions = TEMP['input_conditions']

    OutputNodes = TEMP['output_nodes']
    InputNodes = TEMP['input_nodes']

    # print(model + '- Nodes: ' + str(len(TEMP['Gread'])) + ' Input:' + str(len(InputNodes)) + ' Output:' + str(len(OutputNodes)))

    Mapping = TEMP['mapping']
    InverseMapping = TEMP['inverse_mapping']
    GRead = TEMP['Gread']
    ReadNodes = TEMP['read_nodes']
    # if TEMP['reduction_required']:
    #     GRead = BDOIp.Get_reduced_network(GRead)

    # form expanded network
    # GExpanded = BDOIp.Get_expanded_network(GRead, prefix=Prefix, suffix=Suffix)
    # nx.write_gml(GExpanded, NetworkName + '_expanded_network.gml')
    try:
        GExpanded = nx.read_gml('networks/' + NetworkName + '_expanded_network.gml')
    except IOError:
        # GExpanded = BDOIp.Get_expanded_network(GRead, prefix=Prefix, suffix=Suffix)
        GExpanded = ph.par_get_expanded_network(GRead, prefix=Prefix, suffix=Suffix, worker=num_worker)
        nx.write_gml(GExpanded, 'networks/' + NetworkName + '_expanded_network.gml')

    TempGIOW = jR.get_input_output_relation(GExpanded, Mapping, InverseMapping, InputConditions, OutputNodes,
                                            constant_nodes=[])
    LDOIs = TempGIOW['LDOIs']
    GeneLDOIs = TempGIOW['gene_LDOIs']
    Conflicts = TempGIOW['conflicts']
    GeneConflicts = TempGIOW['gene_conflicts']
    IORelation = TempGIOW['io_relation']
    GRemained = TempGIOW['G_remained']

    R0 = jR.identifying_r0_mutations_ldoi(GExpanded, OutputNodes, Mapping, InverseMapping, InputConditions, IORelation)
    # DM = jR.identifying_disconnecting_mutations(GExpanded, InputNodes, OutputNodes, Mapping, InverseMapping)
    R1 = jR.identifying_r1_mutations_ldoi(GExpanded, OutputNodes, Mapping, InverseMapping, InputConditions, IORelation)
    IIDC = jR.identifying_input_independent_canalizing_mutations(GExpanded, OutputNodes, Mapping, InverseMapping)
    UN = jR.identifying_input_unreachable_nodes(GExpanded, OutputNodes, Mapping, InverseMapping, InputConditions)
    RN = jR.identifying_rn_mutations(GExpanded, Mapping, InverseMapping, InputNodes, OutputNodes, InputConditions,
                                     IORelation,
                                     R0['ineffective_mutations'], R1['r1_mutations'])

    with open('data/' + NetworkName + '_IO.txt', 'w') as f:
        f.write(str(IORelation[0].tolist()))

    CutOffRange = [3, 4, 5]
    GRead.remove_edges_from(list(nx.selfloop_edges(GRead)))
    FFLTest = jR.ffl_test_inter(GRead, ReadNodes, CutOffRange)
    FBLTest = jR.fbl_test(GRead, ReadNodes)

    with open('data/' + NetworkName + '_table_for_original_network.tsv', 'w') as f:
        TMP = '\t'.join(['FFL' + str(CUT) for CUT in CutOffRange])
        f.write('node\tClass\tCanalizing\tUnreachableEffective\t' + TMP + '\tFBL' + '\n')
        NodeList = set(ReadNodes.values())
        NodeList.difference_update(InputNodes)
        NodeList.difference_update(OutputNodes)
        for NODE in NodeList:
            negNODE = '~' + NODE
            canalizing = False
            unreachable = False
            FFL = '\t'.join([str(x) for x in FFLTest[NODE].values()])
            if R0['ineffective'][NODE] and R0['ineffective'][negNODE]:
                nodeClass = 'C0'
            elif NODE in RN['rn_mutations']:
                if RN['rn_mutations'][NODE] == 'R1':
                    nodeClass = 'C1'
                    canalizing = IIDC['iid_canalizing'][NODE]
                    unreachable = UN['input_unreachable'][NODE]
                else:
                    nodeClass = 'C2'
                    canalizing = IIDC['iid_canalizing'][NODE]
                    unreachable = UN['input_unreachable'][NODE]
            elif negNODE in RN['rn_mutations']:
                if RN['rn_mutations'][negNODE] == 'R1':
                    nodeClass = 'C1'
                    canalizing = IIDC['iid_canalizing'][NODE]
                    unreachable = UN['input_unreachable'][NODE]
            else:
                nodeClass = 'C3'
                canalizing = IIDC['iid_canalizing'][NODE]
                unreachable = UN['input_unreachable'][NODE]
            f.write(NODE + '\t' + nodeClass + '\t' + str(canalizing) + '\t'
                    + str(unreachable) + '\t' + FFL + '\t' + str(FBLTest[NODE]) + '\n')