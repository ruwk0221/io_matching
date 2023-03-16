import networkx as nx
import numpy as np
import jReversion as jR

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


def gen_rand_bin_state(bit_num, state_num):
    """

    :param bit_num: int
    :param state_num: int
    :return: set of String
    """
    bin_states = set()
    while len(bin_states) < state_num:
        state = np.array(np.random.random_sample((bit_num, )) < 0.5, dtype=int)
        state = [str(s) for s in state]
        str_state = "".join(state)
        bin_states.add(str_state)
    return bin_states


def expands_bin_states(bin_states, input_idx, input_values):
    """

    :param bin_states: set->String
    :param input_idx: list->int
    :param input_values: dict->(int, String)
    :return:
    """
    input_idx.sort()
    ex_bin_states = set()
    for state in bin_states:
        for idx in input_idx:
            state = state[:idx] + input_values[idx] + state[idx:]
        ex_bin_states.add(state)

    return ex_bin_states

#
# def identifying_asynchronous_attractors(g_read, bin_states, mutation):
#     """
#
#     :param g_read: networkx.Digraph
#     :param bin_states: set
#     :param mutation: dict
#     :return:
#     """
#     states_in_attractors = set()
#     state_transition_map = nx.DiGraph()
#     while len(bin_states) > 0:
#
#         trace = set()
#         trace.add(bin_states.pop())
#         while len(trace) > 0:
#             current_state = trace.pop()
#
#             next_synch_state = synchronous_update(g_read, current_state, mutation)
#             next_asynch_states = set()
#             for idx in range(len(g_read)):
#                 next_state = current_state[:idx] + next_synch_state[idx] + current_state[idx+1:]
#                 next_asynch_states.add(next_state)
#
#             if current_state in next_asynch_states:
#                 state_transition_map.add_edge(current_state, current_state)
#                 next_asynch_states.remove(current_state)
#
#             temp = True
#             while temp and len(next_asynch_states) > 0:
#                 next_asynch_state = next_asynch_states.pop()
#                 if not state_transition_map.has_node(next_asynch_state):
#                     trace.add(next_asynch_state)
#                     temp = False
#
#                 state_transition_map.add_edge(current_state, next_asynch_state)
#             #
#             #
#             # for next_asynch_state in next_asynch_states:
#             #     if not state_transition_map.has_node(next_asynch_state):
#             #         trace.add(next_asynch_state)
#             #
#             #     state_transition_map.add_edge(current_state, next_asynch_state)
#
#     attractors = nx.attracting_components(state_transition_map)
#     for att in attractors:
#         for state in att:
#             states_in_attractors.add(state)
#
#         # state_transition_map.clear()
#
#     output = {"states_in_attractors": states_in_attractors}
#     return output
#
#
def identifying_synchronous_attractors_basins(g_read, bin_states, mutation={}):
    """

    :param g_read: networkx.Digraph
    :param bin_states: set
    :return:
    """
    basin = dict()
    attractors = set()
    type = dict()
    states_in_attractors = dict()

    while len(bin_states) > 0:
        state_trajectory = nx.DiGraph()
        current_state = bin_states.pop()
        next_synch_state = synchronous_update(g_read, current_state, mutation)
        # state_trajectory.add_edge(current_state, next_synch_state)
        while not state_trajectory.has_node(next_synch_state):
            # bin_states.add(next_synch_state)
            state_trajectory.add_edge(current_state, next_synch_state)
            current_state = next_synch_state
            next_synch_state = synchronous_update(g_read, current_state, mutation)
        state_trajectory.add_edge(current_state, next_synch_state)

        attractor = list(nx.attracting_components(state_trajectory))[0]
        assert len(attractor) > 0

        if len(attractor) == 1:
            state = attractor.pop()
            type[state] = 'point'
        else: # cyclic attractor
            states = list(attractor)
            states.sort()
            state = states[0]
            type[state] = 'cyclic'

        if state in attractors:
            basin[state] += 1
        else:
            attractors.add(state)
            basin[state] = 1
            if len(attractor) > 1:
                states_in_attractors[state] = list(attractor)
            else: # len(attractor) of a point attractor is 0 because of the result of pop()
                states_in_attractors[state] = [state]



    # state_transition_map.clear()
    output = {"attractors": attractors,
              "states_in_attractor": states_in_attractors,
              "basin_of_attractor": basin,
              "type_of_attractor": type}
    return output


def get_primary_point_attractor(attractors, basin_of_attractor, type_of_attractor):
    basin_of_point = dict()
    point_attractors = set()
    for att in attractors:
        if type_of_attractor[att] == 'point':
            basin_of_point[att] = basin_of_attractor[att]
            point_attractors.add(att)

    # maximum_basin = max(basin_of_point.values())
    primary_point_attractor = get_primary_attractor(point_attractors, basin_of_point)
    # for att in attractors:
    #     if basin_of_point[att] == maximum_basin:
    #         primary_point_attractor.add(att)

    return primary_point_attractor


def get_primary_attractor(attractors, basin_of_attractor):

    maximum_basin = max(basin_of_attractor.values())
    primary_attractor = set()
    for att in attractors:
        if basin_of_attractor[att] == maximum_basin:
            primary_attractor.add(att)

    return primary_attractor


# def test_robustness():
#     #single mutation에 대해 (또는 mutation 여러개에 대해) primary attractor 가 바뀌는지,
#     # 특히 cyclic attractor 인 경우, 길이의 변화나 point attractor 로 바뀌는 지 등에 대한 다양한 변화 고려 필요
#     # robustness 를 하나의 값으로 제시하는 것은 어려울 수도 있고, 변화 distribution 을 보여줘야 할 수도.
#     # 1. primary attractor(s) 가 여전히 primary attractor(s) 인가?
#     # 2. primary attractor(s) 가 여전히 attractor(s) 인가?
#     # 3. attractor(s) 가 여전히 attractor(s) 인가?
#     return

def robustness_primary_attractor(g_read, state_num, point_only=False):
    bit_num = len(g_read)
    initial_states = gen_rand_bin_state(bit_num, state_num)
    wt_result = identifying_synchronous_attractors_basins(g_read, initial_states.copy())
    if point_only:
        wt_primary_attractor = get_primary_point_attractor(wt_result['attractors'], wt_result['basin_of_attractor'],
                                                           wt_result['type_of_attractor'])
    else:
        wt_primary_attractor = get_primary_attractor(wt_result['attractors'], wt_result['basin_of_attractor'])

    robustness = dict()

    for node in g_read.nodes():
        # initial_states = gen_rand_bin_state(bit_num, state_num)
        oe_result = identifying_synchronous_attractors_basins(g_read, initial_states.copy(), {node: "1"})
        ko_result = identifying_synchronous_attractors_basins(g_read, initial_states.copy(), {node: "0"})
        if point_only:
            oe_primary_attractor = get_primary_point_attractor(oe_result['attractors'],
                                                               oe_result['basin_of_attractor'],
                                                               oe_result['type_of_attractor'])
            ko_primary_attractor = get_primary_point_attractor(ko_result['attractors'],
                                                               ko_result['basin_of_attractor'],
                                                               ko_result['type_of_attractor'])
        else:
            oe_primary_attractor = get_primary_attractor(oe_result['attractors'], oe_result['basin_of_attractor'])
            ko_primary_attractor = get_primary_attractor(ko_result['attractors'], ko_result['basin_of_attractor'])
        if point_only:
            wt_temp = set([x[:node] + x[node+1:] for x in wt_primary_attractor])
            oe_temp = set([x[:node] + x[node+1:] for x in oe_primary_attractor])
            ko_temp = set([x[:node] + x[node+1:] for x in ko_primary_attractor])
        else:
            wt_temp = []
            for att in wt_primary_attractor:
                wt_temp += wt_result['states_in_attractor'][att]
            wt_temp = set([x[:node] + x[node+1:] for x in wt_temp])
            oe_temp = []
            for att in oe_primary_attractor:
                oe_temp += oe_result['states_in_attractor'][att]
            oe_temp = set([x[:node] + x[node + 1:] for x in oe_temp])
            ko_temp = []
            for att in ko_primary_attractor:
                ko_temp += ko_result['states_in_attractor'][att]
            ko_temp = set([x[:node] + x[node + 1:] for x in ko_temp])

        oe_inter = wt_temp.intersection(oe_temp)
        ko_inter = wt_temp.intersection(ko_temp)

        oe_robustness = len(oe_inter) / float(len(wt_temp))
        ko_robustness = len(ko_inter) / float(len(wt_temp))

        robustness[node] = {'oe': oe_robustness, 'ko':ko_robustness, 'noa_wt':len(wt_primary_attractor)}

    return robustness


def par_robustness_primary_attractor(model, directory, state_num, point_only=False):
    prefix, suffix = 'n', 'n'
    temp = jR.cellcollective(model, prefix, suffix, directory)
    g_read = temp['Gread']
    return {'robustness': robustness_primary_attractor(g_read, state_num, point_only),
            'model': model}


def gen_perturbed_initial_states(initial_state, num_of_perturbation=1):
    perturbed_states = set()
    if num_of_perturbation == 1:
        for idx in range(len(initial_state)):
            if initial_state[idx] == '0':
                perturbed_states.add(initial_state[:idx] + '1' + initial_state[idx+1:])
            elif initial_state[idx] == '1':
                perturbed_states.add(initial_state[:idx] + '0' + initial_state[idx + 1:])

    return perturbed_states

def robustness_primary_attractor_perturbation(g_read, state_num, num_of_perturbation=1, point_only=False):
    bit_num = len(g_read)
    initial_states = gen_rand_bin_state(bit_num, state_num)
    wt_result = identifying_synchronous_attractors_basins(g_read, initial_states.copy())
    if point_only:
        wt_primary_attractor = get_primary_point_attractor(wt_result['attractors'], wt_result['basin_of_attractor'],
                                                           wt_result['type_of_attractor'])
    else:
        wt_primary_attractor = get_primary_attractor(wt_result['attractors'], wt_result['basin_of_attractor'])

    perturbed_states = set()
    # for state in initial_states:
    #     perturbed_states = perturbed_states.union(gen_perturbed_initial_states(state, num_of_perturbation))

    for attractor in wt_primary_attractor:
        if wt_result['type_of_attractor'][attractor] == 'point':
            perturbed_states = perturbed_states.union(gen_perturbed_initial_states(attractor, num_of_perturbation))
        else:
            for state in wt_result['states_in_attractor'][attractor]:
                perturbed_states = perturbed_states.union(gen_perturbed_initial_states(state, num_of_perturbation))

    robustness = dict()

    perturbed_result = identifying_synchronous_attractors_basins(g_read, perturbed_states.copy())
    if point_only:
        perturbed_primary_attractor = get_primary_point_attractor(perturbed_result['attractors'],
                                                           perturbed_result['basin_of_attractor'],
                                                           perturbed_result['type_of_attractor'])

    else:
        perturbed_primary_attractor = get_primary_attractor(perturbed_result['attractors'],
                                                            perturbed_result['basin_of_attractor'])


    inter = wt_primary_attractor.intersection(perturbed_primary_attractor)

    robustness = len(inter) / float(len(wt_primary_attractor))

    return robustness

def par_robustness_primary_attractor_perturbation(model, directory, state_num, num_of_perturbation=1, point_only=False):
    prefix, suffix = 'n', 'n'
    temp = jR.cellcollective(model, prefix, suffix, directory)
    g_read = temp['Gread']
    return {'robustness': robustness_primary_attractor_perturbation(g_read, state_num, num_of_perturbation, point_only),
            'model': model}


def robustness_initial_perturbation(g_read, state_num, num_of_perturbation=1, point_only=False):
    bit_num = len(g_read)
    initial_states = gen_rand_bin_state(bit_num, state_num)
    robustness = 0
    for state in initial_states:
        wt_result = identifying_synchronous_attractors_basins(g_read, set([state]))

        perturbed_states = gen_perturbed_initial_states(state, num_of_perturbation)

        for p_state in perturbed_states:
            perturbed_result = identifying_synchronous_attractors_basins(g_read, set([p_state]))
            if list(wt_result['attractors'])[0] == list(perturbed_result['attractors'])[0]:
                robustness += 1


    if num_of_perturbation == 1:
        robustness = float(robustness) / (state_num * bit_num)

    return robustness


def par_robustness_initial_perturbation(model, directory, state_num, num_of_perturbation=1, point_only=False):
    prefix, suffix = 'n', 'n'
    temp = jR.cellcollective(model, prefix, suffix, directory)
    g_read = temp['Gread']
    return {'robustness': robustness_initial_perturbation(g_read, state_num, num_of_perturbation, point_only),
            'model': model}

