import jReversion as jR
import networkx as nx
import numpy as np


def int2state(int_st, n):
    """ Convert integer to state
    """
    fstr = "{0:0%db}" % (n)
    str_st = fstr.format(int_st)
    return str_st


def gen_bin_state(bit_num):
    """

    :param bit_num: int
    :return: set of String
    """
    bin_states = set()
    for state_int in range(2 ** bit_num):
        state = int2state(state_int, bit_num)
        # state = np.array(np.random.random_sample((bit_num, )) < 0.5, dtype=int)
        # state = [str(s) for s in state]
        # str_state = "".join(state)
        bin_states.add(state)
    return bin_states


def get_transition_map(g_read, mutation={}):
    state_transition_map = nx.DiGraph()
    bin_states = gen_bin_state(len(g_read))
    # bin_states = jR.gen_rand_bin_state(len(state_transition_map), 2**len(state_transition_map))
    while len(bin_states) > 0:
        current_state = bin_states.pop()
        next_synch_state = jR.synchronous_update(g_read, current_state, mutation)
        if not state_transition_map.has_node(next_synch_state):
            bin_states.add(next_synch_state)

        state_transition_map.add_edge(current_state, next_synch_state)
    return state_transition_map


def async_get_transition_map(g_read, mutation={}):
    state_transition_map = nx.DiGraph()
    bin_states = gen_bin_state(len(g_read))
    # bin_states = jR.gen_rand_bin_state(len(state_transition_map), 2**len(state_transition_map))
    while len(bin_states) > 0:
        current_state = bin_states.pop()
        for next_asynch_state in jR.asynchronous_update(g_read, current_state, mutation):

            if not state_transition_map.has_node(next_asynch_state):
                bin_states.add(next_asynch_state)

            state_transition_map.add_edge(current_state, next_asynch_state)

    return state_transition_map



# Model = jR.toy_model(111)
# BooleanRuleFileName = Model['BooleanRule_filename']
# NetworkName = Model['network_name']
#
# NumInputs = Model['num_inputs']
# NumInputConditions = Model['num_input_conditions']
#
# InputConditions = Model['input_conditions']
#
# OutputNodes = Model['output_nodes']
# InputNodes = Model['input_nodes']
#
# # Set parameters
# # Note the node name for Gread is the index (integer), one can encode the nodename by adding prefix and suffix
# # If the node name from the input file is not this simple, one need to create a dictionary to record the mapping
# Prefix, Suffix = 'n', 'n'
NetworkName = 'RC_ori'
# BooleanRuleFileName = 'Models/TLGLNetwork_jijoo.booleannet'
BooleanRuleFileName = 'models/'+ NetworkName + '.txt'
Prefix, Suffix = 'n', 'n'

TempI = jR.initialize(BooleanRuleFileName, Prefix, Suffix)


# TempI = jR.initialize(BooleanRuleFileName, Prefix, Suffix)
Mapping = TempI['mapping']
InverseMapping = TempI['inverse_mapping']
GRead = TempI['Gread']
ReadNodes = TempI['read_nodes']

nx.write_gml(G=async_get_transition_map(GRead, mutation={}), path='data/' + NetworkName + '_transition_map.gml')
nx.write_gml(G=async_get_transition_map(GRead, mutation={0:'1'}), path='data/' + NetworkName + '_A1_transition_map.gml')
nx.write_gml(G=async_get_transition_map(GRead, mutation={0:'1', 1:'0'}), path='data/' + NetworkName + '_A1B0_transition_map.gml')
nx.write_gml(G=async_get_transition_map(GRead, mutation={1:'0'}), path='data/' + NetworkName + '_B0_transition_map.gml')
nx.write_gml(G=async_get_transition_map(GRead, mutation={1:'1'}), path='data/' + NetworkName + '_B1_transition_map.gml')
nx.write_gml(G=async_get_transition_map(GRead, mutation={0:'1', 1:'1'}), path='data/' + NetworkName + '_A1B1_transition_map.gml')

