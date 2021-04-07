input_path = "/home/jovyan/work/Tobias Rether/Input/"
output_path = "/home/jovyan/work/Tobias Rether/Output/"

constants = {'f': 96485.33212, 'e': 1.602176634e-19, 'r': 8.314462618, 'na': 6.022140e23}

active_fuel_components = ['CH4', 'CO', 'H2', 'NH3']

steam_reforming = {'reactants': [['H2', 'CO', 'CH4', 'CO2', 'H2O'], [0, 0, 1, 0, 1]],
                   'products': [['H2', 'CO', 'CH4', 'CO2', 'H2O'], [3, 1, 0, 0, 0]]}

reactions_reactants = {'H2': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [1, 0, 0, 0.5, 0, 0, 0, 0]],
                       'CO': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 1, 0, 0.5, 0, 0, 0, 0]],
                       'CH4': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 1, 2, 0, 0, 0, 0]],
                       'NH3': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 3, 0, 0, 4, 0]],
                       'WGS': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 3, 0, 0, 4, 0]]}

reactions_products = {'H2': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 0, 0, 1, 0, 0]],
                      'CO': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 0, 1, 0, 0, 0]],
                      'CH4': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 0, 1, 2, 0, 0]],
                      'NH3': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 0, 0, 6, 0, 2]],
                      'WGS': [['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2'], [0, 0, 0, 3, 0, 0, 4, 0]]}