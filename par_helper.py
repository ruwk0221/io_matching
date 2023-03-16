# import sys
# import os
# sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname('__file__'))))
import jReversion as jR
from LDOI import BooleanDOI_processing as BDOIp
# from LDOI import qm
# import pandas as pd
# import numpy as np
import networkx as nx
import concurrent.futures
# from collections import deque


def par_get_func(Gread, prefix, suffix, equal_sign, node):
    ON_list=[x for x in Gread.nodes[node]['update_rules'].keys() if Gread.nodes[node]['update_rules'][x]==1]
    rule=BDOIp.Getfunc(Gread.nodes[node]['update_nodes'],node,ON_list,prefix=prefix,suffix=suffix,equal_sign=equal_sign)
    OFF_list=[x for x in Gread.nodes[node]['update_rules'].keys() if Gread.nodes[node]['update_rules'][x]==0]
    negation_rule='~'+BDOIp.Getfunc(Gread.nodes[node]['update_nodes'],node,OFF_list,prefix=prefix,suffix=suffix,equal_sign=equal_sign)

    return { 'rule': rule,
             'negation_rule': negation_rule}


def par_get_func2(Gread, prefix, suffix, equal_sign, node_list):
    rules = []
    negation_rules = []
    for node in node_list:
        ON_list=[x for x in Gread.nodes[node]['update_rules'].keys() if Gread.nodes[node]['update_rules'][x]==1]
        rule=BDOIp.Getfunc(Gread.nodes[node]['update_nodes'],node,ON_list,prefix=prefix,suffix=suffix,equal_sign=equal_sign)
        OFF_list=[x for x in Gread.nodes[node]['update_rules'].keys() if Gread.nodes[node]['update_rules'][x]==0]
        negation_rule='~'+BDOIp.Getfunc(Gread.nodes[node]['update_nodes'],node,OFF_list,prefix=prefix,suffix=suffix,equal_sign=equal_sign)
        rules.append(rule)
        negation_rules.append(negation_rule)

    return { 'rules': rules,
             'negation_rules': negation_rules}


def par_get_expanded_network(Gread, prefix='n', suffix='n', equal_sign='*=', worker=10):
    '''
    paralleized by jijoo @ 210121
    Return the expanded network for a given Boolean network model.
    The Boolean network model is a DiGraph object in the output format of form_network().
    The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
    The Boolean rules will first be converted to a disjuctive normal form before generating the expanded network.

    Parameters
    ----------
    Gread     : the given Boolean network model
    prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
    suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
              e.g. node name '1' will become 'n1n' in the returned result
    equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format

    Returns
    -------
    The expanded network for the given Boolean network model.


    '''
    G_expand = nx.DiGraph()
    rules = []
    # first write rules for negation nodes
    negation_rules = []
    expanded_nodes = set()

    with concurrent.futures.ProcessPoolExecutor(max_workers=worker) as executor:
        futures = {executor.submit(par_get_func, Gread.copy(), prefix, suffix, equal_sign, node): node for node in Gread.nodes()}

    for future in concurrent.futures.as_completed(futures):
        node = futures[future]
        temp = future.result()
        temp_rule = temp['rule']
        temp_negation_rule = temp['negation_rule']
        rules.append(temp_rule)
        negation_rules.append(temp_negation_rule)
        expanded_nodes.add(temp_rule.split('*=')[0])
        expanded_nodes.add(temp_negation_rule.split('*=')[0])

    # then for each line in the rules, construct Boolean network
    composite_nodes = []
    rules.extend(negation_rules)
    for line in rules:
        child, update_rule = line.split(equal_sign)
        update_rule = update_rule.strip()
        if update_rule[0] == '(' and update_rule[-1] == ')':
            update_rule = update_rule[1:-1]
        # single parent situation
        if child[0] == '~':
            normal_child = child[1:]
        else:
            normal_child = child[:]
        normal_child = normal_child[len(prefix):len(normal_child) - len(suffix)]
        # deal with source node situation
        if not Gread.nodes[int(normal_child)]['update_nodes']:
            G_expand.add_node(child)  # maybe this do not need to be done
        else:
            if 'or' in update_rule:
                parents = update_rule.split(' or ')
            else:
                parents = [update_rule]
            parents.sort()
            for parent in parents:
                parent = parent.replace('not ', '~').replace('(', '').replace(')', '')
                if 'and' in parent:
                    composite_node = parent.replace(' and ', '_')
                    composite_nodes.append(composite_node)
                    G_expand.add_edge(composite_node, child)
                    for component in composite_node.split('_'):
                        G_expand.add_edge(component, composite_node)
                else:
                    G_expand.add_edge(parent, child)
    return G_expand.copy()


def par_get_expanded_subnetwork2(Gread, rules, prefix, suffix, equal_sign):
    G_expand = nx.DiGraph()
    # composite_nodes = []

    for line in rules:
        child, update_rule = line.split(equal_sign)
        update_rule = update_rule.strip()
        if update_rule[0] == '(' and update_rule[-1] == ')':
            update_rule = update_rule[1:-1]
        # single parent situation
        if child[0] == '~':
            normal_child = child[1:]
        else:
            normal_child = child[:]
        normal_child = normal_child[len(prefix):len(normal_child) - len(suffix)]
        # deal with source node situation
        if not Gread.nodes[int(normal_child)]['update_nodes']:
            G_expand.add_node(child)  # maybe this do not need to be done
        else:
            if 'or' in update_rule:
                parents = update_rule.split(' or ')
            else:
                parents = [update_rule]
            parents.sort()
            for parent in parents:
                parent = parent.replace('not ', '~').replace('(', '').replace(')', '')
                if 'and' in parent:
                    composite_node = parent.replace(' and ', '_')
                    # composite_nodes.append(composite_node)
                    G_expand.add_edge(composite_node, child)
                    for component in composite_node.split('_'):
                        G_expand.add_edge(component, composite_node)
                else:
                    G_expand.add_edge(parent, child)

    return {'G_expand': G_expand,
            # 'composite_nodes': composite_nodes,
            }


def par_get_expanded_network2(Gread, prefix='n', suffix='n', equal_sign='*=', worker=10):
    '''
    paralleized by jijoo @ 210121
    Return the expanded network for a given Boolean network model.
    The Boolean network model is a DiGraph object in the output format of form_network().
    The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
    The Boolean rules will first be converted to a disjuctive normal form before generating the expanded network.

    Parameters
    ----------
    Gread     : the given Boolean network model
    prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
    suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
              e.g. node name '1' will become 'n1n' in the returned result
    equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format

    Returns
    -------
    The expanded network for the given Boolean network model.


    '''
    G_expand = nx.DiGraph()
    rules = []
    # first write rules for negation nodes
    negation_rules = []
    expanded_nodes = set()
    chunk_size = int(len(Gread) / worker) + 1

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i+n]

    # node_chunks = list(chunks(list(Gread.nodes()), chunk_size))
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker) as executor:
        futures = {executor.submit(par_get_func2, Gread.copy(), prefix, suffix, equal_sign, chunk): chunk for chunk in chunks(list(Gread.nodes()), chunk_size)}

    for future in concurrent.futures.as_completed(futures):
        node = futures[future]
        temp = future.result()
        temp_rules = temp['rules']
        temp_negation_rules = temp['negation_rules']
        for temp_rule in temp_rules:
            rules.append(temp_rule)
            expanded_nodes.add(temp_rule.split('*=')[0])

        for temp_negation_rule in temp_negation_rules:
            negation_rules.append(temp_negation_rule)
            expanded_nodes.add(temp_negation_rule.split('*=')[0])

    # then for each line in the rules, construct Boolean network
    composite_nodes = []
    rules.extend(negation_rules)
    chunk_size = int(len(rules) / worker) + 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker) as executor:
        futures = {executor.submit(par_get_expanded_subnetwork2, Gread.copy(), rule_chunk, prefix, suffix, equal_sign): rule_chunk for rule_chunk in chunks(rules, chunk_size)}

    for future in concurrent.futures.as_completed(futures):
        temp = future.result()
        # composite_nodes.extend(temp['composite_nodes'])
        G_expand = nx.compose(G_expand, temp['G_expand'])

    return G_expand.copy()


def par_get_expanded_network3(Gread, prefix='n', suffix='n', equal_sign='*=', worker=10):
    '''
    paralleized by jijoo @ 210121
    Return the expanded network for a given Boolean network model.
    The Boolean network model is a DiGraph object in the output format of form_network().
    The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
    The Boolean rules will first be converted to a disjuctive normal form before generating the expanded network.

    Parameters
    ----------
    Gread     : the given Boolean network model
    prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
    suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
              e.g. node name '1' will become 'n1n' in the returned result
    equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format

    Returns
    -------
    The expanded network for the given Boolean network model.


    '''
    G_expand = nx.DiGraph()
    rules = []
    # first write rules for negation nodes
    # negation_rules = []
    # expanded_nodes = set()
    # chunk_size = int(len(Gread) / worker) + 1

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i+n]

    negation_rules = []
    expanded_nodes = set()
    for node in Gread.nodes():
        ON_list = [x for x in Gread.nodes[node]['update_rules'].keys() if Gread.nodes[node]['update_rules'][x] == 1]
        rule = BDOIp.Getfunc(Gread.nodes[node]['update_nodes'], node, ON_list, prefix=prefix, suffix=suffix,
                       equal_sign=equal_sign)
        OFF_list = [x for x in Gread.nodes[node]['update_rules'].keys() if Gread.nodes[node]['update_rules'][x] == 0]
        negation_rule = '~' + BDOIp.Getfunc(Gread.nodes[node]['update_nodes'], node, OFF_list, prefix=prefix, suffix=suffix,
                                      equal_sign=equal_sign)
        rules.append(rule)
        negation_rules.append(negation_rule)
        expanded_nodes.add(rule.split('*=')[0])
        expanded_nodes.add(negation_rule.split('*=')[0])

    # then for each line in the rules, construct Boolean network
    composite_nodes = []
    rules.extend(negation_rules)
    chunk_size = int(len(rules) / worker) + 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker) as executor:
        futures = {executor.submit(par_get_expanded_subnetwork2, Gread.copy(), rule_chunk, prefix, suffix, equal_sign): rule_chunk for rule_chunk in chunks(rules, chunk_size)}

    for future in concurrent.futures.as_completed(futures):
        temp = future.result()
        # composite_nodes.extend(temp['composite_nodes'])
        G_expand = nx.compose(G_expand, temp['G_expand'])

    return G_expand.copy()


