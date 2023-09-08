

configurations = {
    'prefix': '',
    'pseudo_dir': '../',
    'qe_atoms_file': 'QE_ATOMS',
    'qe_cell_file': 'QE_CELL',
    'kgrid_coarse': [1,1,1],
    'qshift_wfnq': [0,0,0],
    'kgrid_fine': [1,1,1],
    'kpoints_file_scf': 'KPOINTS_FILE_SCF',
    'kpoints_file_wfn': 'KPOINTS_FILE_WFN',
    'kpoints_file_wfnq': 'KPOINTS_FILE_WFNq',
    'kpoints_file_wfn_fi': 'KPOINTS_FILE_WFN_fi',
    'kpoints_file_wfnq_fi': 'KPOINTS_FILE_WFNq_fi',
    'nval': 0,
    'ncond_wfn': 0,
    'ncond_wfnq': 0,
    'ncond_wfnco': 0,
    'ncond_wfnfi': 0,
    'nmin_gw': 0,
    'nmax_gw': 0,
    'nval_kernel': 0,
    'ncond_kernel': 0,
    'nval_bse': 0,
    'ncond_bse': 0,
    'dft_ekin_cutoff': 40,
    'qshift_bse': [0,0,0],
    'epsilon_cutoff': 10, 
    'truncation_scheme': ''
}
        
import configparser
import ast

config_from_file = configparser.ConfigParser()
config_from_file.read("workflow_input")

print('Reading input file workflow_input')

for option, value in config_from_file['VARS'].items():
    print(f'Option: {option}, Value: {value}')
    # check if value is the same kind of variable as in default configurations
    value = ast.literal_eval(value)
    print('!!!!') 
    # print(value)
    # print(type(value))
    # print(type(configurations[option]))
    print('before mod', configurations[option])

    if type(value) == type(configurations[option]):

        if type(value) == list: # check if lists have the same size!
            if len(value) == len(configurations[option]):
                configurations[option] = value
            else:
                print(f'Expected {len(configurations[option])} values but got {len(value)} for variable {option}')
        else: # if not a list just load the data
            configurations[option] = value
                
    else:
        # checking if the user wrote one value that was suposed to be float but it is int
        if type(value) == int and type(configurations[option]) == float: 
            configurations[option] = float(value)
        else:
            print(f'Variable {option} in input file is {type(value)} but should be {type(configurations[option])}')

    print('after mod', configurations[option])