# model_name_list = ['BT474', 'cardiac', 'CD4_T', 'death_receptor', 'lymphopoesis', 'macrophage',
#                    'SASP', 'yeast_apoptosis']
model_name_list = ['BT474']

for model_name in model_name_list:
    with open('models/' + model_name + '.txt', 'w') as f:
        f.write('# Boolean rules\n')
        f.write('# Inputs\n')
        with open('cellcollective/' + model_name + '/expr/external_components.ALL.txt', 'r') as ecf:
            for line in ecf:
                input_node = line.rstrip('\n')
                f.write(input_node + '*= ' + input_node + '\n')

        f.write('# Internal nodes\n')
        with open('cellcollective/' + model_name + '/expr/expressions.ALL.txt', 'r') as exf:
            for line in exf:
                logic = line.replace(' =', '*=')
                f.write(logic)




