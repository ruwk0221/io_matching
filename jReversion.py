import networkx as nx
import numpy as np
from LDOI import BooleanDOI_DOI as BDOId, BooleanDOI_TargetControl as BDOItc, BooleanDOI_processing as BDOIp
import itertools
import random

def cellcollective(model_name, prefix, suffix, directory=''
                   , boolean_rule_filename=''
                   ):
    if boolean_rule_filename == '':
        boolean_rule_filename = directory + 'models/' + model_name + '.txt'
    network_name = model_name
    num_inputs = 0
    input_nodes = []
    lines = []

    # lines.append('# Boolean rules\n')
    # lines.append('# Inputs\n')
    with open(directory + 'cellCollective/' + model_name + '/expr/external_components.ALL.txt', 'r') as ecf:
        for line in ecf:
            input_node = ' '.join(line.split())
            lines.append(input_node + '*= ' + input_node + '\n')
            num_inputs += 1
            input_nodes.append(input_node)

    # lines.append('# Internal nodes\n')
    with open(directory + 'cellCollective/' + model_name + '/expr/expressions.ALL.txt', 'r') as exf:
        for line in exf:
            line = ' '.join(line.split())
            logic = line.replace(' =', '*=')
            lines.append(logic)

    num_input_conditions = 2 ** num_inputs
    input_conditions = np.ndarray((num_inputs, 2), dtype=object)
    for idx, input_node in enumerate(input_nodes):
        input_conditions[idx, 0] = '~' + input_node
        input_conditions[idx, 1] = input_node



    Gread, read_nodes = BDOIp.form_network(lines, sorted_nodename=False)

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    read_nodes_dict = {}
    for i, node in enumerate(read_nodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node
        read_nodes_dict[i] = node
    nonself_Gread = Gread.copy()
    nonself_Gread.remove_edges_from(list(nx.selfloop_edges(Gread)))
    output_node_idx = [idx for idx in Gread.nodes if nonself_Gread.out_degree(idx) == 0]
    output_nodes = []
    output_nodes.extend([read_nodes[idx] for idx in output_node_idx])
    output_nodes.extend(['~'+read_nodes[idx] for idx in output_node_idx])
    input_nodes = input_conditions.reshape(num_inputs * 2, ).tolist()

    output = {"BooleanRule_filename": boolean_rule_filename,
              "network_name": network_name,
              "num_inputs": num_inputs,
              "num_input_conditions": num_input_conditions,
              "input_conditions": input_conditions,
              "output_nodes": output_nodes,
              "input_nodes": input_nodes,
              "reduction_required": False,
              'Gread': Gread,
              'mapping': mapping,
              'inverse_mapping': inverse_mapping,
              'read_nodes': read_nodes_dict}

    return output



def get_input_output_relation(G_expanded, mapping, inverse_mapping, input_conditions,
                              output_nodes, constant_nodes=[]):
    num_inputs, _ = input_conditions.shape
    num_input_conditions = 2 ** num_inputs

    results = np.ndarray((1, num_input_conditions), dtype=object)
    conflicts = np.ndarray((1, num_input_conditions), dtype=object)
    LDOIs = np.ndarray((1, num_input_conditions), dtype=object)
    gene_LDOIs = np.ndarray((1, num_input_conditions), dtype=object)
    gene_conflicts = np.ndarray((1, num_input_conditions), dtype=object)
    G_remained = np.ndarray((1, num_input_conditions), dtype=object)
    G_LDOI = np.ndarray((1, num_input_conditions), dtype=object)
    cpnode_counts = np.ndarray((1, num_input_conditions), dtype=object)

    for i in range(num_input_conditions):
        G_copied = G_expanded.copy()
        temp = np.binary_repr(i, width=num_inputs)
        input_condition = [input_conditions[x, int(temp[x])] for x in range(num_inputs)]
        source_nodes = input_condition + constant_nodes
        source = set([mapping[x] for x in source_nodes])
        for node in source:
            G_copied.remove_edges_from(list(G_copied.in_edges(node)))
            G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(node))))

        output = set([mapping[x] for x in output_nodes])
        LDOI_raw, _, conflict, cpnode = BDOId.truncated_node_of_influence_BFS(G_copied, source)

        LDOI = set([x for x in LDOI_raw if x.find('_') < 0])
        LDOIs[0, i] = LDOI_raw.copy()
        output.intersection_update(LDOI)
        results[0, i] = sorted([inverse_mapping[x] for x in output])
        gene_LDOIs[0, i] = sorted([inverse_mapping[x] for x in LDOI])
        conflict = set([x for x in conflict if x.find('_') < 0])
        conflicts[0, i] = conflict
        gene_conflicts[0, i] = [inverse_mapping[x] for x in conflict]
        cpnode_counts[0, i] = cpnode
        G_copied.remove_nodes_from(LDOI_raw)
        G_copied.remove_nodes_from(conflict)
        G_remained[0, i] = G_copied.copy()
        LDOI_raw.update(cpnode.keys())
        G_LDOI[0, i] = nx.subgraph(G_expanded, LDOI_raw)

    output = {'LDOIs': LDOIs,
              'gene_LDOIs': gene_LDOIs,
              'conflicts': conflicts,
              'gene_conflicts': gene_conflicts,
              'io_relation': results,
              'G_remained': G_remained,
              'G_LDOI': G_LDOI,
              'cpnode_counts': cpnode_counts}
    return output



def identifying_r0_mutations_ldoi(G_expanded, output_nodes, mapping, inverse_mapping, input_conditions, wt_io):
    num_inputs, _ = input_conditions.shape
    num_input_conditions = 2 ** num_inputs
    node_list = np.array(list(inverse_mapping))
    num_nodes = len(node_list)

    results = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    conflicts = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    LDOIs = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    gene_LDOIs = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    gene_conflicts = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    ineffective = dict()
    mut_io = dict()
    for ni, node in enumerate(node_list):
        G_copied = G_expanded.copy()
        G_copied.remove_edges_from(list(G_copied.in_edges(node)))
        G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(node))))
        for i in range(num_input_conditions):
            temp = np.binary_repr(i, width=num_inputs)
            input_condition = [input_conditions[x, int(temp[x])] for x in range(num_inputs)]
            mutation_nodes = [inverse_mapping[node]]
            neg_mutation_node = BDOItc.Negation_in_expanded(mutation_nodes[0])

            source_nodes = input_condition + mutation_nodes
            source = set([mapping[x] for x in source_nodes])
            output = set([mapping[x] for x in output_nodes])
            LDOI_raw, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(G_copied, source)

            LDOI = set([x for x in LDOI_raw if x.find('_') < 0])
            LDOIs[ni, i] = LDOI.copy()
            output.intersection_update(LDOI)
            results[ni, i] = sorted([inverse_mapping[x] for x in output])
            gene_LDOIs[ni, i] = sorted([inverse_mapping[x] for x in LDOIs[ni, i]])
            conflict = set([x for x in conflict if x.find('_') < 0])
            conflicts[ni, i] = conflict
            gene_conflicts[ni, i] = [inverse_mapping[x] for x in conflict]
        ineffective[inverse_mapping[node]] = (str(wt_io[0].tolist()) == str(results[ni].tolist()))
        mut_io[inverse_mapping[node]] = str(results[ni].tolist())
    # irreversible_nodes = [inverse_mapping[x] for x in node_list[irreversibility.values()]]
    ineffective_mutations = [k for k, v in ineffective.items() if v]
    effective_mutations = [k for k, v in ineffective.items() if not v]

    output = {'ineffective': ineffective,
              'ineffective_mutations': ineffective_mutations,
              'effective_mutations': effective_mutations,
              'mut_io': mut_io}
    return output




def initialize(file, prefix, suffix):
    with open(file, 'rU') as f:
        lines = f.readlines()
    Gread, read_nodes = BDOIp.form_network(lines, sorted_nodename=False)

    # calculate the mapping from string nodename to index
    mapping = {}  # nodename to number index
    inverse_mapping = {}  # number index to nodename
    read_nodes_dict = {}
    inverse_read_nodes_dict = {}
    for i, node in enumerate(read_nodes):
        index = prefix + str(i) + suffix
        mapping[node] = index
        inverse_mapping[index] = node
        mapping['~' + node] = '~' + index
        inverse_mapping['~' + index] = '~' + node
        read_nodes_dict[i] = node
        inverse_read_nodes_dict[node] = i
    output = {'Gread': Gread,
              'mapping': mapping,
              'inverse_mapping': inverse_mapping,
              'read_nodes': read_nodes_dict,
              'inverse_read_nodes': inverse_read_nodes_dict}
    return output


def all_reachability(G, source_list, output_list):
    reachability= np.ndarray((1, len(source_list), len(output_list)), dtype=bool)
    for ns, s in enumerate(source_list):
        for no, o in enumerate(output_list):
            reachability[0, ns, no] = nx.has_path(G, s, o)
    return reachability


def identifying_disconnecting_mutations(G_expanded, input_nodes, output_nodes, mapping, inverse_mapping):
    num_inputs = len(input_nodes)
    num_outputs = len(output_nodes)
    # num_input_conditions = 2 ** num_inputs
    node_list = np.array(list(inverse_mapping))
    num_nodes = len(node_list)

    results = np.ndarray((num_nodes, num_inputs, num_outputs), dtype=object)
    # conflicts = np.ndarray((num_nodes, num_inputs), dtype=object)
    # LDOIs = np.ndarray((num_nodes, num_inputs), dtype=object)
    # gene_LDOIs = np.ndarray((num_nodes, num_inputs), dtype=object)
    # gene_conflicts = np.ndarray((num_nodes, num_inputs), dtype=object)
    disconnecting = dict()

    source = set([mapping[x] for x in input_nodes])
    output = set([mapping[x] for x in output_nodes])
    wt_reachability = all_reachability(G_expanded, source, output)

    for ni, node in enumerate(node_list):
        G_copied = G_expanded.copy()

        # G_copied.remove_edges_from(G_copied.in_edges(node))
        G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(node))))
        results[ni] = all_reachability(G_copied, source, output)
        disconnecting[inverse_mapping[node]] = (str(wt_reachability[0].tolist()) != str(results[ni].tolist()))
    # disconnecting_nodes = [inverse_mapping[x] for x in node_list[disconnecting.values()]]
    disconnecting_mutations = [k for k, v in disconnecting.items() if v]

    output = {'disconnecting': disconnecting,
              'disconnecting_mutations': disconnecting_mutations}
    return output


def identifying_r1_mutations_ldoi(G_expanded, output_nodes, mapping, inverse_mapping, input_conditions, wt_io):
    num_inputs, _ = input_conditions.shape
    num_input_conditions = 2 ** num_inputs
    node_list = np.array(list(inverse_mapping))
    num_nodes = len(node_list)

    results_mut = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    results_neu = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    conflicts = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    LDOIs = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    gene_LDOIs = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    gene_conflicts = np.ndarray((num_nodes, num_input_conditions), dtype=object)
    r1 = dict()
    for ni, node in enumerate(node_list):
        G_copied = G_expanded.copy()
        G_copied.remove_edges_from(list(G_copied.in_edges(node)))
        G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(node))))

        # G_copied.remove_edges_from(G_copied.in_edges(BDOId.Negation_in_expanded(node)))
        for i in range(num_input_conditions):
            temp = np.binary_repr(i, width=num_inputs)
            input_condition = [input_conditions[x, int(temp[x])] for x in range(num_inputs)]
            # mutation_nodes = [inverse_mapping[node]]
            # neg_mutation_node = BDOItc.Negation_in_expanded(mutation_nodes[0])
            # temp = copy.deepcopy(constant_nodes)
            # if neg_mutation_node in constant_nodes:  # avoiding conflict
            #     temp.remove(neg_mutation_node)

            source_nodes = input_condition
            source_mut = set([mapping[x] for x in source_nodes])
            source_mut.add(node)
            output = set([mapping[x] for x in output_nodes])
            LDOI_mut, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(G_copied, source_mut)
            LDOIs[ni, i] = LDOI_mut.copy()
            source_neu = set([mapping[x] for x in source_nodes])
            source_neu.add(BDOId.Negation_in_expanded(node))
            LDOI_neu, _, _, _ = BDOId.truncated_node_of_influence_BFS(G_copied, source_neu)
            LDOI_mut = set([x for x in LDOI_mut if x.find('_') < 0])
            output_mut = output.intersection(LDOI_mut)
            results_mut[ni, i] = sorted([inverse_mapping[x] for x in output_mut])
            LDOI_neu = set([x for x in LDOI_neu if x.find('_') < 0])
            output_neu = output.intersection(LDOI_neu)
            results_neu[ni, i] = sorted([inverse_mapping[x] for x in output_neu])
            gene_LDOIs[ni, i] = sorted([inverse_mapping[x] for x in LDOI_mut])
            conflict = set([x for x in conflict if x.find('_') < 0])
            conflicts[ni, i] = conflict
            gene_conflicts[ni, i] = [inverse_mapping[x] for x in conflict]
        # disconnecting[inverse_mapping[node]] = (str(wt_io[0].tolist()) != str(results[ni].tolist()))
        r1[inverse_mapping[node]] = ((str(wt_io[0].tolist()) != str(results_mut[ni].tolist()))
                                              and (str(wt_io[0].tolist()) == str(results_neu[ni].tolist())))


    # irreversible_nodes = [inverse_mapping[x] for x in node_list[irreversibility.values()]]
    r1_mutations = [k for k, v in r1.items() if v]

    output = {'r1': r1,
              'r1_mutations': r1_mutations
              }
    return output



def identifying_input_independent_canalizing_mutations(G_expanded, output_nodes, mapping, inverse_mapping):
    # num_inputs, _ = input_conditions.shape
    # num_input_conditions = 2 ** num_inputs
    node_list = np.array(list(inverse_mapping))
    num_nodes = len(node_list)

    results = np.ndarray((num_nodes, ), dtype=object)
    conflicts = np.ndarray((num_nodes,), dtype=object)
    LDOIs = np.ndarray((num_nodes, ), dtype=object)
    gene_LDOIs = np.ndarray((num_nodes, ), dtype=object)
    gene_conflicts = np.ndarray((num_nodes, ), dtype=object)
    iid_canalizing = dict()
    for ni, node in enumerate(node_list):
        G_copied = G_expanded.copy()
        G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(node))))

        source = [node]
        # source = set([mapping[x] for x in source_nodes])
        output = set([mapping[x] for x in output_nodes])
        LDOI_raw, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(G_copied, source)

        LDOI = set([x for x in LDOI_raw if x.find('_') < 0])
        LDOIs[ni] = LDOI.copy()
        output.intersection_update(LDOI)
        results[ni] = sorted([inverse_mapping[x] for x in output])
        gene_LDOIs[ni] = sorted([inverse_mapping[x] for x in LDOIs[ni]])
        conflict = set([x for x in conflict if x.find('_') < 0])
        conflicts[ni] = conflict
        gene_conflicts[ni] = [inverse_mapping[x] for x in conflict]
        # disconnecting[inverse_mapping[node]] = (str(wt_io[0].tolist()) != str(results[ni].tolist()))
        iid_canalizing[inverse_mapping[node]] = len(output) > 0
    # irreversible_nodes = [inverse_mapping[x] for x in node_list[irreversibility.values()]]
    iid_canalizing_mutations = [k for k, v in iid_canalizing.items() if v]

    output = {'iid_canalizing': iid_canalizing,
              'iid_canalizing_mutations': iid_canalizing_mutations}
    return output


def identifying_input_unreachable_nodes(G_expanded, output_nodes, mapping, inverse_mapping, input_conditions):
    num_inputs, _ = input_conditions.shape
    num_input_conditions = 2 ** num_inputs
    node_list = np.array(list(inverse_mapping))
    num_nodes = len(node_list)

    results = np.ndarray((num_input_conditions, ), dtype=object)
    # conflicts = np.ndarray((num_input_conditions, ), dtype=object)
    LDOIs = np.ndarray((num_input_conditions, ), dtype=object)
    gene_LDOIs = np.ndarray((num_input_conditions, ), dtype=object)
    # gene_conflicts = np.ndarray((num_input_conditions, ), dtype=object)
    input_unreachable = dict()

    G_copied = G_expanded.copy()
    input_reachable = set()
    for i in range(num_input_conditions):
        temp = np.binary_repr(i, width=num_inputs)
        input_condition = [input_conditions[x, int(temp[x])] for x in range(num_inputs)]
        # mutation_nodes = [inverse_mapping[node]]
        # neg_mutation_node = BDOItc.Negation_in_expanded(mutation_nodes[0])
        # temp = copy.deepcopy(constant_nodes)
        # if neg_mutation_node in constant_nodes:  # avoiding conflict
        #     temp.remove(neg_mutation_node)

        source_nodes = input_condition
        source = set([mapping[x] for x in source_nodes])
        output = set([mapping[x] for x in output_nodes])
        LDOI_raw, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(G_copied, source)

        LDOI = set([x for x in LDOI_raw if x.find('_') < 0])
        LDOIs[i] = LDOI.copy()
        # output.intersection_update(LDOI)
        # results[i] = [inverse_mapping[x] for x in output]
        gene_LDOIs[i] = sorted([inverse_mapping[x] for x in LDOIs[i]])
        input_reachable = input_reachable.union(set(LDOIs[i]))
        # conflict = set([x for x in conflict if x.find('_') < 0])
        # conflicts[ni, i] = conflict
        # gene_conflicts[ni, i] = [inverse_mapping[x] for x in conflict]
    # disconnecting[inverse_mapping[node]] = (str(wt_io[0].tolist()) != str(results[ni].tolist()))

    for node in node_list:
        input_unreachable[inverse_mapping[node]] = node not in input_reachable
    # irreversible_nodes = [inverse_mapping[x] for x in node_list[irreversibility.values()]]
    input_unreachable_nodes = [k for k, v in input_unreachable.items() if v]

    output = {'input_unreachable': input_unreachable,
              'input_unreachable_nodes': input_unreachable_nodes}
    return output


def ineffective_test_ldoi(G_expanded, output_nodes, mapping, inverse_mapping, input_conditions, wt_io, test_set):
    num_inputs, _ = input_conditions.shape
    num_input_conditions = 2 ** num_inputs
    # node_list = np.array(list(inverse_mapping))
    # num_nodes = len(node_list)
    results = np.ndarray((1, num_input_conditions), dtype=object)
    # conflicts = np.ndarray((1, num_input_conditions), dtype=object)
    # LDOIs = np.ndarray((1, num_input_conditions), dtype=object)
    # gene_LDOIs = np.ndarray((1, num_input_conditions), dtype=object)
    # gene_conflicts = np.ndarray((1, num_input_conditions), dtype=object)
    G_copied = G_expanded.copy()
    mutation_nodes =[]
    for node in test_set:
        G_copied.remove_edges_from(list(G_copied.in_edges(mapping[node])))
        G_copied.remove_edges_from(list(G_copied.in_edges(BDOId.Negation_in_expanded(mapping[node]))))
        mutation_nodes.append(node)

    for i in range(num_input_conditions):
        temp = np.binary_repr(i, width=num_inputs)
        input_condition = [input_conditions[x, int(temp[x])] for x in range(num_inputs)]

        # neg_mutation_node = BDOItc.Negation_in_expanded(mutation_nodes[0])

        source_nodes = input_condition + mutation_nodes
        source = set([mapping[x] for x in source_nodes])
        output = set([mapping[x] for x in output_nodes])
        LDOI_raw, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(G_copied, source)

        LDOI = set([x for x in LDOI_raw if x.find('_') < 0])
        # LDOIs_temp = LDOI.copy()
        output.intersection_update(LDOI)
        results[0, i] = sorted([inverse_mapping[x] for x in output])
        # gene_LDOIs_temp = [inverse_mapping[x] for x in LDOIs_temp]
        # conflict = set([x for x in conflict if x.find('_') < 0])
        # conflicts_temp = conflict
        # gene_conflicts_temp = [inverse_mapping[x] for x in conflict]
    ineffective = (str(wt_io[0].tolist()) == str(results[0].tolist()))
    return ineffective


def identifying_rn_mutations(g_expanded, mapping, inverse_mapping, input_nodes, output_nodes, input_conditions, wt_io,
                             ineffective_mutations, r1_mutations):
    candidates = set(mapping)
    candidates.difference_update(set(input_nodes))
    candidates.difference_update(set(output_nodes))
    candidates.difference_update(set(ineffective_mutations))
    candidates.difference_update(set(r1_mutations))

    rn_mutations = dict()
    rn_solutions = dict()

    #initialize
    # for im in ineffective_mutations:
    #     rn_mutations[im] = "IM"
    #     rn_solutions[im] = "IM"

    for r1 in r1_mutations:
        rn_mutations[r1] = "R1"
        rn_solutions[r1] = [BDOId.Negation_in_expanded(r1)]

    rn = r1_mutations
    n = 1
    while len(rn) > 0:
        n += 1
        candidates.difference_update(set(rn))
        rn = []

        for node in candidates:
            # ex_node = mapping[node]
            ldoi_raw, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(g_expanded, [mapping[node]])
            ldoi = set([x for x in ldoi_raw if x.find('_') < 0])
            ldoi_gene = [inverse_mapping[x] for x in ldoi]
            ldoi_rn = [x for x in ldoi_gene if x in rn_mutations]

            if len(ldoi_rn) > 0:
                # neg_ldoi_rn = [BDOId.Negation_in_expanded(x) for x in ldoi_rn]
                # g_copied = g_expanded.copy()
                # g_copied.remove_edges_from(g_copied.in_edges(mapping[node]))
                # g_copied.remove_edges_from(g_copied.in_edges(BDOId.Negation_in_expanded(mapping[node])))
                # for target in neg_ldoi_r1:
                #     g_copied.remove_edges_from(g_copied.in_edges(target))
                #     g_copied.remove_edges_from(g_copied.in_edges(BDOId.Negation_in_expanded(target)))
                ineffectivity = False
                k = 0
                while (k < len(ldoi_rn)) and (not ineffectivity):
                    for subset in itertools.combinations(ldoi_rn, k+1):
                        control = []
                        for tt in subset:
                            control += rn_solutions[tt]
                        source = control + [node]
                        ineffectivity = ineffective_test_ldoi(g_expanded, output_nodes, mapping, inverse_mapping,
                                                              input_conditions, wt_io, source)
                        if ineffectivity:
                            rn_mutations[node] = "R" + str(n)
                            rn_solutions[node] = control
                            # candidates.remove(node)
                            rn.append(node)
                            # rn_mutations[BDOId.Negation_in_expanded(node)] = "R" + str(n)
                            # tmp = rn_solutions[node] + [node]
                            # rn_solutions[BDOId.Negation_in_expanded(node)] = tmp
                            # rn.append(BDOId.Negation_in_expanded(node))
                            break
                    k += 1
        candidates.difference_update(set(rn))
        for node in candidates:
            if (BDOId.Negation_in_expanded(node) in rn_mutations) and not (node in rn_mutations):
                neg_node = BDOId.Negation_in_expanded(node)
                rn_mutations[node] = rn_mutations[neg_node]
                tmp = rn_solutions[neg_node] + [neg_node]
                rn_solutions[node] = tmp
                rn.append(node)

    for node in set(mapping):
        if not (node in rn_mutations):
            rn_solutions[node] = 'NOT RN'

    return {'rn_mutations': rn_mutations,
            'rn_solutions': rn_solutions}


def identifying_r2_mutations(g_expanded, mapping, inverse_mapping, input_nodes, output_nodes, input_conditions, wt_io,
                             ineffective_mutations, r1_mutations):
    candidates = set(mapping)
    candidates.difference_update(set(input_nodes))
    candidates.difference_update(set(output_nodes))
    candidates.difference_update(set(ineffective_mutations))
    candidates.difference_update(set(r1_mutations))

    r2_mutations = set()
    r2_solutions = dict()

    for node in candidates:
        # ex_node = mapping[node]
        ldoi_raw, _, conflict, _ = BDOId.truncated_node_of_influence_BFS(g_expanded, [mapping[node]])
        ldoi = set([x for x in ldoi_raw if x.find('_') < 0])
        ldoi_gene = [inverse_mapping[x] for x in ldoi]
        ldoi_r1 = [x for x in ldoi_gene if x in r1_mutations]

        if len(ldoi_r1) > 0:
            neg_ldoi_r1 = [BDOId.Negation_in_expanded(x) for x in ldoi_r1]
            # g_copied = g_expanded.copy()
            # g_copied.remove_edges_from(g_copied.in_edges(mapping[node]))
            # g_copied.remove_edges_from(g_copied.in_edges(BDOId.Negation_in_expanded(mapping[node])))
            # for target in neg_ldoi_r1:
            #     g_copied.remove_edges_from(g_copied.in_edges(target))
            #     g_copied.remove_edges_from(g_copied.in_edges(BDOId.Negation_in_expanded(target)))
            ineffectivity = False
            k = 0
            while (k < len(neg_ldoi_r1)) and (not ineffectivity):
                for subset in itertools.combinations(neg_ldoi_r1, k+1):
                    source = list(subset) + [node]
                    ineffectivity = ineffective_test_ldoi(g_expanded, output_nodes, mapping, inverse_mapping, input_conditions, wt_io, source)
                    if ineffectivity:
                        r2_mutations.add(node)
                        r2_solutions[node] = list(subset)
                        # r2_mutations.add(BDOId.Negation_in_expanded(node))
                        # tmp = r2_solutions[node] + [node]
                        # r2_solutions[BDOId.Negation_in_expanded(node)] = tmp

                        break
                k += 1

    # for node in r2_mutations:
    #     r2_mutations.add(BDOId.Negation_in_expanded(node))
    #     tmp = r2_solutions[node] + [node]
    #     r2_solutions[BDOId.Negation_in_expanded(node)] = tmp

    for node in set(mapping):
        if not (node in r2_mutations):
            neg_node = BDOId.Negation_in_expanded(node)
            if neg_node in r2_mutations:
                r2_mutations.add(node)
                tmp = r2_solutions[neg_node] + [neg_node]
                r2_solutions[node] = tmp
            else:
                r2_solutions[node] = 'NOT R2'
        # if not (node in r2_mutations):
        #     r2_solutions[node] = 'NOT R2'

    return {'r2_mutations': r2_mutations,
            'r2_solutions': r2_solutions}

def test_r2(g_expanded, mapping, inverse_mapping, output_nodes, input_conditions, wt_io,
            r1_mutation):
    """

    :param g_expanded: nx.DiGraph
    :param mapping: dict()
    :param inverse_mapping: dict()
    :param input_nodes:
    :param output_nodes:
    :param input_conditions:
    :param wt_io:
    :param r1_mutation:
    """
    children = g_expanded.successors(mapping[r1_mutation])
    children = [child for child in children if child.find('_') < 0]
    test_set = [BDOId.Negation_in_expanded(child) for child in children]
    test_set.append(mapping[r1_mutation])
    test_set = set([inverse_mapping[m] for m in test_set])

    return ineffective_test_ldoi(g_expanded, output_nodes, mapping, inverse_mapping, input_conditions, wt_io, test_set)


def test_r2_ver2(g_expanded, mapping, inverse_mapping, output_nodes, input_conditions, wt_io,
                 r1_mutation, r1_mutations=[]):

    """

    :param r1_mutations:
    :param g_expanded: nx.DiGraph
    :param mapping: dict()
    :param inverse_mapping: dict()
    :param input_nodes:
    :param output_nodes:
    :param input_conditions:
    :param wt_io:
    :param r1_mutation:
    """
    children = g_expanded.successors(mapping[r1_mutation])
    children = set([child for child in children if child.find('_') < 0]).intersection([mapping[r1] for r1 in r1_mutations])
    test_set = [BDOId.Negation_in_expanded(child) for child in children]
    test_set.append(mapping[r1_mutation])
    test_set = set([inverse_mapping[m] for m in test_set])

    return ineffective_test_ldoi(g_expanded, output_nodes, mapping, inverse_mapping, input_conditions, wt_io, test_set)


def random_combination_set(elements, combination_length, set_length):
    combinations = set()
    while len(combinations) < set_length:
        combinations.add(tuple(sorted(random.sample(elements, combination_length))))

    return combinations


def synchronous_update(g_read, bool_state, mutation={}):
    """

    :type g_read: object
    :type bool_state: String
    :type mutation: dict->(int, String)
    :return next_bool_state: String
    """
    next_state = ""
    for node_idx in sorted(g_read.nodes):
        temp = ""
        for regulator_idx in g_read.nodes[node_idx]['update_nodes']:
            temp += bool_state[regulator_idx]
        next_state += str(g_read.nodes[node_idx]['update_rules'][temp])

    for key, val in mutation.items():
        next_state = next_state[:key] + val + next_state[key + 1:]

    return next_state


def asynchronous_update(g_read, bool_state, mutation={}):
    next_synch_state = synchronous_update(g_read, bool_state, mutation)
    next_asynch_states = set()
    for idx in range(len(g_read)):
        next_state = bool_state[:idx] + next_synch_state[idx] + bool_state[idx+1:]
        next_asynch_states.add(next_state)

    return next_asynch_states

#
# def gen_rand_bin_state(bit_num, state_num):
#     """
#
#     :param bit_num: int
#     :param state_num: int
#     :return: set of String
#     """
#     bin_states = set()
#     while len(bin_states) <= state_num:
#         state = np.array(np.random.random_sample((bit_num, )) < 0.5, dtype=int)
#         state = [str(s) for s in state]
#         str_state = "".join(state)
#         bin_states.add(str_state)
#     return bin_states
#
#
# def expands_bin_states(bin_states, input_idx, input_values):
#     """
#
#     :param bin_states: set->String
#     :param input_idx: list->int
#     :param input_values: dict->(int, String)
#     :return:
#     """
#     input_idx.sort()
#     ex_bin_states = set()
#     for state in bin_states:
#         for idx in input_idx:
#             state = state[:idx] + input_values[idx] + state[idx:]
#         ex_bin_states.add(state)
#
#     return ex_bin_states



def randomize_regulator(g_read, mapping, input_nodes, output_nodes, p=0.2):
    g_copy = g_read.copy()
    input_nodes_idx = set([int(mapping[x].replace('n', '').replace('~','')) for x in input_nodes])
    output_nodes_idx = set([int(mapping[x].replace('n', '').replace('~', '')) for x in output_nodes])
    internal_nodes_idx = set(range(len(g_copy))).difference(input_nodes_idx.union(output_nodes_idx))

    for ni in internal_nodes_idx.union(output_nodes_idx):
        if np.random.rand() < p:
            num_regulator = len(g_copy.nodes[ni]['update_nodes'])
            temp = internal_nodes_idx.copy()
            temp = temp.union(input_nodes_idx)
            temp.discard(ni)
            new_regulator = random.sample(temp, num_regulator)
            g_copy.nodes[ni]['update_nodes'] = new_regulator
            g_copy.remove_edges_from(list(g_copy.in_edges(ni)))
            g_copy.add_edges_from(zip(new_regulator, np.tile([ni], num_regulator)))

    return g_copy


def randomize_target(g_read, mapping, input_nodes, output_nodes, p=0.2):
    g_copy = g_read.copy()
    input_nodes_idx = set([int(mapping[x].replace('n', '').replace('~','')) for x in input_nodes])
    output_nodes_idx = set([int(mapping[x].replace('n', '').replace('~', '')) for x in output_nodes])
    internal_nodes_idx = set(range(len(g_copy))).difference(input_nodes_idx.union(output_nodes_idx))

    for ni in internal_nodes_idx.union(output_nodes_idx):
        if np.random.rand() < p:
            num_regulator = len(g_copy.nodes[ni]['update_nodes'])
            temp = internal_nodes_idx.copy()
            temp = temp.union(input_nodes_idx)
            temp.discard(ni)
            new_regulator = random.sample(temp, num_regulator)
            # g_copy.nodes[ni]['update_nodes'] = new_regulator
            g_copy.remove_edges_from(list(g_copy.out_edges(ni)))
            g_copy.add_edges_from(zip(np.tile([ni], num_regulator), new_regulator))

    return g_copy


def identify_ffl(g_read, read_nodes, cutoff=None):
    num_node = len(g_read)
    mx_num_sp = np.zeros((num_node, num_node), dtype=int)
    mx_min_sp = np.zeros((num_node, num_node), dtype=int)
    mx_max_sp = np.zeros((num_node, num_node), dtype=int)
    for source in read_nodes.keys():
        for target in read_nodes.keys():
            num_sp = 0
            min_sp = num_node
            max_sp = 0
            if source is not target:
                for path in nx.all_simple_paths(g_read, source, target, cutoff):
                    num_sp += 1
                    length = len(path)
                    min_sp = min([length, min_sp])
                    max_sp = max([length, max_sp])
            mx_num_sp[source, target] = num_sp
            mx_min_sp[source, target] = min_sp
            mx_max_sp[source, target] = max_sp

    mx_ffl = mx_num_sp >= 2

    output = {'mx_num_sp': mx_num_sp,
              'mx_min_sp': mx_min_sp,
              'mx_max_sp': mx_max_sp,
              'mx_ffl': mx_ffl}

    return output


def identify_fbl(g_read, read_nodes):
    num_node = len(g_read)
    mx_num_sp = np.zeros((num_node, ), dtype=int)
    mx_min_sp = np.zeros((num_node, ), dtype=int)
    mx_max_sp = np.zeros((num_node, ), dtype=int)
    for source in read_nodes.keys():

        num_sp = 0
        min_sp = num_node
        max_sp = 0

        for path in nx.all_simple_paths(g_read, source, source):
            num_sp += 1
            length = len(path)
            min_sp = min([length, min_sp])
            max_sp = max([length, max_sp])
        mx_num_sp[source] = num_sp
        mx_min_sp[source] = min_sp
        mx_max_sp[source] = max_sp

    mx_fbl = mx_num_sp >= 1

    output = {'mx_num_sp': mx_num_sp,
              'mx_min_sp': mx_min_sp,
              'mx_max_sp': mx_max_sp,
              'mx_fbl': mx_fbl}

    return output


def node_deletion_effects_on_ffl(g_read, read_nodes, cutoff=None):
    # num_node = len(g_read)
    node_deletion_effects = dict()
    for mut in read_nodes.keys():
        g_mut = g_read.copy()
        g_mut.remove_node(mut)
        g_mut.add_node(mut)
        node_deletion_effects[mut] = identify_ffl(g_mut, read_nodes, cutoff)

    return node_deletion_effects


def node_deletion_effects_on_fbl(g_read, read_nodes):
    # num_node = len(g_read)
    node_deletion_effects = dict()
    for mut in read_nodes.keys():
        g_mut = g_read.copy()
        g_mut.remove_node(mut)
        g_mut.add_node(mut)
        node_deletion_effects[mut] = identify_fbl(g_mut, read_nodes)

    return node_deletion_effects


def ffl_test_inter(g_read, read_nodes, cut_off_range=[None]):
    result = dict()
    for node in read_nodes.keys():
        result[read_nodes[node]] = dict()

    for cut_off in cut_off_range:
        ffl = identify_ffl(g_read, read_nodes, cut_off)
        eff = node_deletion_effects_on_ffl(g_read, read_nodes, cut_off)

        for node in read_nodes.keys():

            test1 = ffl['mx_num_sp'] != eff[node]['mx_num_sp']
            # test1 = abs(test1) == 1
            # test2 = np.logical_and(FFL['mx_ffl'], EFF[node]['mx_num_sp'] > 0)
            test2 = np.logical_and(test1, ffl['mx_ffl'])
            # test = np.logical_and(test3, EFF[node]['mx_num_sp'] > 0)

            test = eff[node]['mx_num_sp']
            test = np.delete(test, node, axis=0)
            test = np.delete(test, node, axis=1)

            test2 = np.delete(test2, node, axis=0)
            test2 = np.delete(test2, node, axis=1)

            if test2.any():
                result[read_nodes[node]][cut_off] = (test[test2] > 0).all()
            else:
                result[read_nodes[node]][cut_off] = False

            # test1 = ffl['mx_num_sp'] != eff[node]['mx_num_sp']
            # # test1 = abs(test1) == 1
            # # test2 = np.logical_and(FFL['mx_ffl'], EFF[node]['mx_num_sp'] > 0)
            # test3 = np.logical_and(test1, ffl['mx_ffl'])
            # test = np.logical_and(test3, eff[node]['mx_num_sp'] > 0)
            #
            # # test = np.delete(test, node, axis=0)
            # # test = np.delete(test, node, axis=1)
            # result[ReadNodes[node]][cut_off] = np.sum(test) == np.sum(test3)

    return result


def ffl_test_disc(g_read, read_nodes, cut_off_range=[None]):
    result = dict()
    for node in read_nodes.keys():
        result[read_nodes[node]] = dict()

    for cut_off in cut_off_range:
        ffl = identify_ffl(g_read, read_nodes, cut_off)
        eff = node_deletion_effects_on_ffl(g_read, read_nodes, cut_off)

        for node in read_nodes.keys():

            test1 = ffl['mx_num_sp'] != eff[node]['mx_num_sp'] # node deletion 에 의해 sp num 이 바뀐 pair
            # test1 = abs(test1) == 1
            # test2 = np.logical_and(FFL['mx_ffl'], EFF[node]['mx_num_sp'] > 0)
            test2 = np.logical_and(test1, ffl['mx_ffl']) # node deletion 에 의해 바뀐 pair 중 원래 ffl 인 pair
            # test = np.logical_and(test3, EFF[node]['mx_num_sp'] > 0)

            test = eff[node]['mx_num_sp'] # node deletion 후 sp num
            test = np.delete(test, node, axis=0)
            test = np.delete(test, node, axis=1)

            test2 = np.delete(test2, node, axis=0)
            test2 = np.delete(test2, node, axis=1)

            if test2.any(): # ffl 에 변화가 있다
                result[read_nodes[node]][cut_off] = (test[test2] == 0).any() # ffl 중 하나라도 disconnecting 된다
            else:
                result[read_nodes[node]][cut_off] = False

            # test1 = ffl['mx_num_sp'] != eff[node]['mx_num_sp']
            # # test1 = abs(test1) == 1
            # # test2 = np.logical_and(FFL['mx_ffl'], EFF[node]['mx_num_sp'] > 0)
            # test3 = np.logical_and(test1, ffl['mx_ffl'])
            # test = np.logical_and(test3, eff[node]['mx_num_sp'] > 0)
            #
            # # test = np.delete(test, node, axis=0)
            # # test = np.delete(test, node, axis=1)
            # result[ReadNodes[node]][cut_off] = np.sum(test) == np.sum(test3)

    return result


def ffl_test(g_read, read_nodes, cut_off_range=[None]):
    result = dict()
    for node in read_nodes.keys():
        result[read_nodes[node]] = dict()

    for cut_off in cut_off_range:
        ffl = identify_ffl(g_read, read_nodes, cut_off)
        eff = node_deletion_effects_on_ffl(g_read, read_nodes, cut_off)

        for node in read_nodes.keys():

            test1 = ffl['mx_num_sp'] != eff[node]['mx_num_sp'] # node deletion 에 의해 sp num 이 바뀐 pair
            test2 = np.logical_and(test1, ffl['mx_ffl']) # node deletion 에 의해 바뀐 pair 중 원래 ffl 인 pair

            test = eff[node]['mx_num_sp'] # node deletion 후 sp num
            test = np.delete(test, node, axis=0)
            test = np.delete(test, node, axis=1)

            test2 = np.delete(test2, node, axis=0)
            test2 = np.delete(test2, node, axis=1)

            if test2.any(): # ffl 에 변화가 있다
                if (test[test2]==0).any(): # ffl 중 하나라도 disconnecting 된다
                    result[read_nodes[node]][cut_off] = 'Critical'
                elif (test[test2] > 0).all(): # ffl에 변화 있지만 connectivity 유지
                    result[read_nodes[node]][cut_off] = 'Redundant'
                else:
                    result[read_nodes[node]][cut_off] = 'Unexpected'
            else:
                result[read_nodes[node]][cut_off] = 'Unexpected'

    return result


def fbl_test(g_read, read_nodes):
    result = dict()
    fbl = identify_fbl(g_read, read_nodes)
    for node in read_nodes.keys():
        eff = node_deletion_effects_on_fbl(g_read, read_nodes)
        test = np.logical_and(fbl['mx_fbl'], eff[node]['mx_fbl'])
        test = np.delete(test, node)
        test2 = np.delete(fbl['mx_fbl'], node)
        result[read_nodes[node]] = np.sum(test2) - np.sum(test)

    return result






