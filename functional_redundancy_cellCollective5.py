import sys
import os
# sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname('__file__'))))
import jReversion as jR
import par_helper as ph
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, balanced_accuracy_score
import networkx as nx
from statannot import add_stat_annotation
import csv

num_worker = 5
networkModel =['bortezomib',
                # 'igvh',
                'apoptosis',
                # 'aurora',
                'bt474_long',
                # 'bt474_short',
                # 'cd4t',
                'colitis',
                'death',
                # 'egfr',
                # 'erbb',
                # 'fa_brca',
                # 'fa_check',
                'hcc1954_long',
                # 'hcc1954_short',
                'hgf',
                'mammalian',
                # 'mammalian_2006',
                'mapk',
                'oxidative',
                # 'pro_inflammatory',
                # 'fibroblasts',
                'skbr3_long',
                # 'skbr3_short',
                'tlgl_2008',
                'tlgl_2011',
                # 'tlgl_2011_reduced',
                # 'prostate',
                'migration']
# networkModel = ['mammalian']
io_pathway_counts = dict()
for Model in networkModel:
    Prefix, Suffix = 'n', 'n'
    TEMP = jR.cellcollective(Model, Prefix, Suffix, directory='')

    BooleanRuleFileName = TEMP['BooleanRule_filename']
    NetworkName = TEMP['network_name']
    print(NetworkName)
    NumInputs = TEMP['num_inputs']
    NumInputConditions = TEMP['num_input_conditions']

    InputConditions = TEMP['input_conditions']

    OutputNodes = TEMP['output_nodes']
    InputNodes = TEMP['input_nodes']

    Mapping = TEMP['mapping']
    InverseMapping = TEMP['inverse_mapping']
    GRead = TEMP['Gread']
    ReadNodes = TEMP['read_nodes']

    InverseReadNodes = dict()
    for key,val in ReadNodes.items():
        InverseReadNodes[val] = key
    #
    # try:
    #     GExpanded = nx.read_gml('networks/' + NetworkName + '_expanded_network.gml')
    # except IOError:
    #     # GExpanded = BDOIp.Get_expanded_network(GRead, prefix=Prefix, suffix=Suffix)
    #     GExpanded = ph.par_get_expanded_network(GRead, prefix=Prefix, suffix=Suffix, worker=num_worker)
    #     nx.write_gml(GExpanded, 'networks/' + NetworkName + '_expanded_network.gml')

    # TempGIOW = jR.get_input_output_relation(GExpanded, Mapping, InverseMapping, InputConditions, OutputNodes,
    #                                         constant_nodes=[])
    # LDOIs = TempGIOW['LDOIs']
    # GeneLDOIs = TempGIOW['gene_LDOIs']
    # # Conflicts = TempGIOW['conflicts']
    # # GeneConflicts = TempGIOW['gene_conflicts']
    # IORelation = TempGIOW['io_relation']
    # # LDOIs_union = set()
    # # for idx in range(NumInputConditions):
    # #     LDOIs_union = LDOIs_union.union(LDOIs[0, idx])
    # # LDOIs_subnetwork = nx.subgraph(GExpanded, LDOIs_union)
    io_pathway_count = 0
    io_count = 0
    LDOI_union = set()
    # for LDOI in LDOIs[0,:]:
    #     LDOI_union = LDOI_union.union(LDOI)
    # LDOI_union = LDOI_union.union([Mapping[x] for x in InputConditions.reshape((2*NumInputs,))])
    # GSub = nx.subgraph(GExpanded, LDOI_union)
    # for i in range(NumInputConditions):
    #     # G_copied = G_expanded.copy()
    #     BinaryBit = np.binary_repr(i, width=NumInputs)
    #     InputCondition = [InputConditions[x, int(BinaryBit[x])] for x in range(NumInputs)]
    #     # source_nodes = input_condition + constant_nodes
    #     Sources = set([x for x in InputCondition])
    #     # for node in source:
    #     #     G_copied.remove_edges_from(list(G_copied.in_edges(node)))
    #     #     G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(node))))
    #
    #     Targets = set([x for x in IORelation[0, i]])
    Sources = InputNodes
    Targets = OutputNodes

    for s in Sources:
        for t in Targets:
            if nx.has_path(G=GRead, source=InverseReadNodes[s.replace('~', '')], target=InverseReadNodes[t.replace('~', '')]):
                io_count += 1
                len_shortest_path = nx.shortest_path_length(G=GRead, source=InverseReadNodes[s.replace('~', '')], target=InverseReadNodes[t.replace('~', '')])
                for p in nx.all_simple_paths(G=GRead, source=InverseReadNodes[s.replace('~', '')], target=InverseReadNodes[t.replace('~', '')]):
                    io_pathway_count += 1

    if io_count == 0:
        io_pathway_counts[Model] = str(-1)
    else:
        io_pathway_counts[Model] = str(io_pathway_count/float(io_count))

    # GRemained = TempGIOW['G_remained']
    #     CutOffRange = [3, 4, 5, 6, 7]
    # GRead.remove_edges_from(list(nx.selfloop_edges(GRead)))
    # num_cycle = 0
    # for cycle in nx.simple_cycles(G=GRead):
    #     num_cycle += 1
    # feedback_counts[network_name] = num_cycle

with open('io_pathway_counts_cellCollective7.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(io_pathway_counts.keys())
    writer.writerow(io_pathway_counts.values())

