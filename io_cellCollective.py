import jReversion as jR
from LDOI import BooleanDOI_processing as BDOIp
import networkx as nx

models = ['bortezomib',
    # 'igvh',
    'apoptosis',
    # 'aurora',
    'bt474_long',
    'bt474_short',
    'cd4t',
    'colitis',
    'death',
    'egfr',
    'erbb',
    'fa_brca',
    'fa_check',
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

for model in models:
    Prefix, Suffix = 'n', 'n'
    TEMP = jR.cellcollective(model, Prefix, Suffix)

    BooleanRuleFileName = TEMP['BooleanRule_filename']
    NetworkName = TEMP['network_name']

    NumInputs = TEMP['num_inputs']
    NumInputConditions = TEMP['num_input_conditions']

    InputConditions = TEMP['input_conditions']

    OutputNodes = TEMP['output_nodes']
    InputNodes = TEMP['input_nodes']

    print(model + '- Nodes: ' + str(len(TEMP['Gread'])) + ' Input:' + str(len(InputNodes)) + ' Output:' + str(len(OutputNodes)))