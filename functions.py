__author__ = "Tobias Rether"
__copyright__ = "Copyright 2021, Siemens Energy"
__department__ = "SE GP G SV TI EN RW PRS COE"
__email__ = "tobias.rether@siemens-energy.com"
__status__ = "Demo"

from input_parameters import*
# import plotly.offline as py '
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# Packages to use for processing e.g. numpy for numerical processing,
# cantera for electrochemical numerical calculations etc.
# import os
# import csv

import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import scipy as sy
import scipy.optimize as so
import xlsxwriter as Excel
import cantera as ct
# from datetime import date, time, datetime, timedelta
# from dateutil import tz
import sympy as sym
from sympy.solvers import solve
from sympy import Symbol
from scipy.optimize import fsolve
# import math
# import time


def composition_string(gas_composition_input_format):
    ''' Conversion of species-molar fraction string to cantera input format.

    Type: Function

    Input:
    gas_composition_input_format - gas object composition in format
    [[species_1, ... , species_n], [molar_fraction_1, ... , molar_fraction_n]]

    Output:
    gas_composition_output_format - gas object composition in format
    species_1:molar_fraction_1, ... , species_n:molar_fraction_n '''

    gas_composition_output_format = "self"

    if len(gas_composition_input_format[0]) == 1:

        gas_composition_output_format = gas_composition_input_format[0][0] + ":" + str(gas_composition_input_format[1][0])

    elif len(gas_composition_input_format[0]) > 1:

        for counter in range(len(gas_composition_input_format[0])):

            if counter == 0:
                gas_composition_output_format = gas_composition_input_format[0][counter] + ":" \
                                                + str(gas_composition_input_format[1][counter]) + ", "

            elif 0 < counter < len(gas_composition_input_format[0])-1:
                gas_composition_output_format = gas_composition_output_format + gas_composition_input_format[0][counter] \
                                                + ":" + str(gas_composition_input_format[1][counter]) + ", "

            elif counter == len(gas_composition_input_format[0])-1:
                gas_composition_output_format = gas_composition_output_format + gas_composition_input_format[0][counter] \
                                                + ":" + str(gas_composition_input_format[1][counter])

    return gas_composition_output_format

def gas_object(gas_composition: object, t: object, p: object) -> object:
    ''' Creation of Cantera Gas Object Phase based on given gas composition, Temperature, Pressure
        based on gri30 and Nasa Polynoms 9

        Type: Function

        Sub functions: composition_string

        Input:
        gas_composition - gas object composition in format
        [[species_1, ... , species_n], [molar_fraction_1, ... , molar_fraction_n]]
        t - Temperature in K
        p - pressure in Pa

        Output:
        gas_phase - Cantera gas phase object calculated w/ NASA7 Polynoms from gri30 input file
                    based on given composition, temperature and pressure from Input
        '''
    gas_phase = ct.Solution('gri30.xml')
    components = composition_string(gas_composition)
    gas_phase.TPX = t, p, components

    return gas_phase


def object_to_composition(gas_phase):
    # Decompile species and molar fractions of gas phase to input format in two lists,
    # one for species, second for molar fractions
    species_list = []
    molar_fraction_list = []
    for counter in range(len(gas_phase.species_names)):

        species = str(gas_phase.species_name(counter))
        molar_fraction = gas_phase.X[counter]

        if molar_fraction > 0:
            species_list.append(species)
            molar_fraction_list.append(molar_fraction)

    composition = []
    composition.append(species_list)
    composition.append(molar_fraction_list)

    return composition



def sensible_enthalpy(gas_phase, ref_type = None):
    t_actual = gas_phase.TP[0]
    p_actual = gas_phase.TP[1]

    t_ref_0 = 273.15
    t_ref_15 = 288.15
    t_ref_25 = 298.15
    p_ref = 101325

    h_meas_mass = gas_phase.enthalpy_mass
    h_meas_mole = gas_phase.enthalpy_mole

    h_ref_mass = None
    h_ref_mole = None

    if ref_type is None or ref_type == 1:

        gas_phase.TP = t_ref_0, p_ref
        h_ref_mass = gas_phase.enthalpy_mass
        h_ref_mole = gas_phase.enthalpy_mole

    elif ref_type == 2:

        gas_phase.TP = t_ref_15, p_ref
        h_ref_mass = gas_phase.enthalpy_mass
        h_ref_mole = gas_phase.enthalpy_mole

    elif ref_type == 3:

        gas_phase.TP = t_ref_25, p_ref
        h_ref_mass = gas_phase.enthalpy_mass
        h_ref_mole = gas_phase.enthalpy_mole

    h_sens_mass = h_meas_mass - h_ref_mass
    h_sens_mole = h_meas_mole - h_ref_mole

    gas_phase.TP = t_actual, p_actual

    return h_sens_mass, h_sens_mole, h_meas_mass, h_meas_mole


def fuel_to_air_ratio(oxidator, fuel, phi, m_fuel_start = None):

    m_oxidator = 1
    m_fuel = None

    if m_fuel_start is None:
        m_fuel = 1
    elif m_fuel_start is not None:
        m_fuel = m_fuel_start

    fuel_composition = composition_string(fuel)
    oxidator_phase = gas_object(oxidator, 288.15, 101325)
    fuel_phase = gas_object(fuel, 288.15, 101325)

    oxidator_bulk = ct.Quantity(oxidator_phase, mass = m_oxidator)
    fuel_bulk = ct.Quantity(fuel_phase, mass = m_fuel)

    mix = oxidator_bulk + fuel_bulk
    phi_start = mix.get_equivalence_ratio()

    phi_target_value = phi
    delta_phi = abs(phi_target_value - phi_start)

    counter = 0

    while delta_phi > 0.01:

        counter += 1
        # print(counter)
        factor = 1 + 1 / (counter)

        if (phi_target_value - phi_start) > 0:
            m_fuel = factor * m_fuel

        elif (phi_target_value - phi_start) < 0:
            m_fuel = m_fuel / factor

        oxidator_phase = gas_object(oxidator, 288.15, 101325)
        fuel_phase = gas_object(fuel, 288.15, 101325)

        oxidator_bulk = ct.Quantity(oxidator_phase, mass=m_oxidator)
        fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel)

        mix = oxidator_bulk + fuel_bulk

        phi_start = mix.get_equivalence_ratio()
        delta_phi = abs(phi_target_value - phi_start)

    fuel_to_air = m_fuel / m_oxidator

    return fuel_to_air


def heating_value(fuel):
    #Returns the LHV and HHV for the specified fuel
    gas = ct.Solution('gri30.cti')
    gas.TP = 298, ct.one_atm
    gas.set_equivalence_ratio(1.0, fuel, 'O2:1.0')
    h1 = gas.enthalpy_mass
    Y_fuel = gas[fuel].Y[0]

    # complete combustion products
    Y_products = {'CO2': gas.elemental_mole_fraction('C'),
                  'H2O': 0.5 * gas.elemental_mole_fraction('H'),
                  'N2': 0.5 * gas.elemental_mole_fraction('N')}

    gas.TPX = None, None, Y_products
    Y_H2O = gas['H2O'].Y[0]
    h2 = gas.enthalpy_mass
    LHV = -(h2 - h1) / Y_fuel

    water = ct.Water()
    # Set liquid water state, with vapor fraction x = 0
    water.TQ = 298, 0
    h_liquid = water.h
    # Set gaseous water state, with vapor fraction x = 1
    water.TQ = 298, 1
    h_gas = water.h

    HHV = -(h2 - h1 + (h_liquid - h_gas) * Y_H2O) / Y_fuel

    return LHV, HHV


def heating_value_mix(fuel_mix):
    t = 298.15
    p = 101325

    # start= time.time()

    fuel = gas_object(fuel_mix, t, p)
    LHV_average = 0
    HHV_average = 0

    for counter in range(len(fuel_mix[0])):
        LHV, HHV = heating_value(fuel_mix[0][counter])
        index = fuel.species_index(fuel_mix[0][counter])
        mass_fraction = fuel.Y[index]
        LHV_average = LHV_average + LHV * mass_fraction
        HHV_average = HHV_average + HHV * mass_fraction

    LHV_average = LHV_average / 1000
    HHV_average = HHV_average / 1000

    # end = time.time()
    # time_span= end-start
    # print(time_span)
    return LHV_average, HHV_average

def fuel_to_air_ratio2(oxidizer_phase, fuel_phase, phi):

    fuel_string = composition_string(object_to_composition(fuel_phase))

    gas = ct.Solution('gri30.yaml')
    gas.set_equivalence_ratio(phi, fuel=fuel_string, oxidizer='O2')

    index = gas.species_index('O2')
    oxygen_molar_fraction = gas.X[index]
    oxygen_mass_fraction = gas.Y[index]

    index = oxidizer_phase.species_index('O2')
    oxidizer_oxygen_molar_fraction = oxidizer_phase.X[index]
    oxidizer_oxygen_mass_fraction = oxidizer_phase.Y[index]

    fuel_to_air_mass_ratio = (1-oxygen_mass_fraction)/(oxygen_mass_fraction/oxidizer_oxygen_mass_fraction)

    return fuel_to_air_mass_ratio

def create_state_dataframe(gas_bulk, position):

    t_ref = 298.15
    p_ref = 101325

    species = ['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2']

    header = ['position', 'T actual [K]', 'p actual [Pa]', 'roh actual [kg/m3]', 'm [kg/s]', 'V actual [m3/s]',
              'n [kmol/s]', 'T ref [K]', 'p ref [Pa]', 'roh ref [kg/m3]', 'V ref [Nm3/s]', 'h in [kJ/kg]',
              's in [kJ/kgK]', 'cp mass [kJ/kgK]', 'cv mass[kJ/kgK]', 'kappa [-]',
              'H2 mol percentage [mol%]', 'CO mol percentage [mol%]', 'CH4 mol percentage [mol%]',
              'O2 mol percentage [mol%]', 'CO2 mol percentage [mol%]', 'H2O mol percentage [mol%]',
              'NH3 mol percentage [mol%]', 'N2 mol percentage [mol%]', 'H2 weight percentage [weight%]',
              'CO weight percentage [weight%]',
              'CH4 weight percentage [weight%]', 'O2 weight percentage [weight%]',
              'CO2 weight percentage [weight%]', 'H2O weight percentage [weight%]',
              'NH3 weight percentage [weight%]', 'N2 weight percentage [weight%]']

    t_act = float(gas_bulk.phase.T)
    p_act = float(gas_bulk.phase.P)
    roh_act = float(gas_bulk.phase.density_mass)
    massflow = float(gas_bulk.mass)
    volume_flow_act = float(gas_bulk.volume)
    molar_flow = float(gas_bulk.moles)

    enthalpy_mole_act = float(gas_bulk.phase.enthalpy_mole)
    enthalpy_mass_act = float(gas_bulk.phase.enthalpy_mass)
    entropy_mole_act = float(gas_bulk.phase.entropy_mole)
    entropy_mass_act = float(gas_bulk.phase.entropy_mass)

    cp_mass = float(gas_bulk.phase.cp_mass)
    cv_mass = float(gas_bulk.phase.cv_mass)
    kappa = float(cp_mass/cv_mass)

    cp_moles = float(gas_bulk.phase.cp_mole)
    cv_moles = float(gas_bulk.phase.cv_mole)

    gas_bulk_ref = ct.Quantity(gas_bulk.phase, mass=gas_bulk.mass)
    gas_bulk_ref.TP = t_ref, p_ref
    roh_ref = float(gas_bulk_ref.phase.density_mass)
    volume_flow_ref = float(gas_bulk_ref.volume)
    enthalpy_mole_ref = float(gas_bulk_ref.phase.enthalpy_mole)
    enthalpy_mass_ref = float(gas_bulk_ref.phase.enthalpy_mass)
    entropy_mole_ref = float(gas_bulk_ref.phase.entropy_mole)
    entropy_mass_ref = float(gas_bulk_ref.phase.entropy_mass)

    enthalpy = float(enthalpy_mass_act-enthalpy_mass_ref)
    entropy = float(entropy_mass_act-entropy_mass_ref)

    value = [position, t_act, p_act, roh_act, massflow, volume_flow_act, molar_flow, t_ref, p_ref, roh_ref,
             volume_flow_ref, enthalpy, entropy, cp_mass, cv_mass, kappa]

    for counter in range(len(species)):

        index = gas_bulk.phase.species_index(species[counter])

        value.append(float(gas_bulk.phase[index].X[0]))

    for counter in range(len(species)):
        index = gas_bulk.phase.species_index(species[counter])

        value.append(float(gas_bulk.phase[index].Y[0]))

    df = pd.DataFrame(data=value, index=header)

    return df

def lsv_to_compressor_map():
    file_path = input('Enter Filepath:') # e.g. C:/Users/Guest/Desktop
    file_name = input('Enter File-Name: ') # File name including file type e.g. Test.xlsx
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    if '/' in file_path:
        file_path = file_path + '/' + file_name
    else:
        file_path = file_path + '\\' + file_name
        file_path = file_path.replace('\\','/')
    if '"' in file_path:
        file_path = file_path.replace('"','')
    if not os.path.isfile(file_path):
        print('File does not exist.')
        exit()
    else:
        data = pd.read_excel(file_path)
        df = pd.DataFrame(data)
        print(df)
    val1 = input('Enter LSV: ')
    val2 = input('Enter TEB: ')
    try:
        a=df[(df['LSV'] == int(val1)) & (df['TEB'] == int(val2))].index[0]
    except IndexError:
        print('LSV and TEB incorrect.')
        exit()
    data1 = [{'MV1': df.at[int(a),'MV1'],'MENEXT%': df.at[int(a),'MENEXT%'],'PIV': df.at[int(a),'PIV'],'ETAVI': df.at[int(a),'ETAVI']}]
    index = 'LSV'+': '+str(val1)+'; '+'TEB'+': '+str(val2)
    df1 = pd.DataFrame(data1, index=[index])
    print(df1)

def igv_to_lsv():

    return 1

################# SOFC functions ###########################

def fuel_species_fractions(fuel_phase, fuel_active_species_vector):
    ''' Determination of molar fraction of single active fuel species related to sum of all molar fractions of
        active fuel species in fuel phase object defined in list of fuel_active_species_vector

        Type: Function

        Sub functions: -

        Input:
        fuel_phase - Fuel Cantera Phase Object
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]

        Output:
        molar_fraction_active_species - List of molar fractions of active fuel species in order of given list
                                        'fuel_active_species_vector'
    '''
    molar_fraction_all = 0
    molar_fraction_active_species = []

    for counter in range(len(fuel_active_species_vector)):
        active_species = fuel_active_species_vector[counter]
        index = fuel_phase.species_index(active_species)
        molar_fraction = fuel_phase.X[index]
        molar_fraction_all = molar_fraction_all + molar_fraction

    for counter in range(len(fuel_active_species_vector)):
        active_species = fuel_active_species_vector[counter]
        index = fuel_phase.species_index(active_species)
        molar_fraction = fuel_phase.X[index]
        molar_fraction_active_species.append(molar_fraction / molar_fraction_all)

    return molar_fraction_active_species


def fuel_species_charge(fuel_active_species_vector):
    ''' Determination of electron charge for oxidation of each list element of fuel active species list

        Type: Function

        Sub functions: -

        Input:
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]

        Output:
        fuel_active_species_charge_vector - List of charge of each list element/species in 'fuel_active_species_vector'
                                            in corresponding order

    '''

    fuel_active_species_charge_vector = []
    for counter in range(len(fuel_active_species_vector)):
        active_species = fuel_active_species_vector[counter]
        reactants = reactions_reactants[active_species]

        index = reactants[0].index(active_species)

        fuel_species_factor = reactants[1][index]

        index = reactants[0].index("O2")

        Oxygen_factor = reactants[1][index]

        fuel_active_species_charge_vector.append(Oxygen_factor / fuel_species_factor * 4)

    return fuel_active_species_charge_vector


def i_fraction(fuel_active_species_vector, fuel_active_species_molar_fraction_vector,
               fuel_active_species_charge_vector):
    ''' Determination of current fraction based in active fuel species charge and molar fraction refered to all
        fuel active species

            Type: Function

            Sub functions: -

            Input:
            fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]
            fuel_active_species_molar_fraction_vector - List of molar fractions of active fuel species in order
                                                        of given list 'fuel_active_species_vector'
            fuel_active_species_charge_vector - List of charge of each list element/species in
                                                'fuel_active_species_vector' in corresponding order

            Output:
            i_fraction_vector - List of current fractions of each list element/species
                                in 'fuel_active_species_vector' in corresponding order
        '''

    i_sum_fractions = 0
    i_fraction_vector = []

    for counter in range(len(fuel_active_species_vector)):
        i_sum_fractions = i_sum_fractions + float(fuel_active_species_molar_fraction_vector[counter]) * \
                          float(fuel_active_species_charge_vector[counter])

    for counter in range(len(fuel_active_species_vector)):
        i_fraction_vector.append(float(fuel_active_species_molar_fraction_vector[counter]) *
                      float(fuel_active_species_charge_vector[counter])/i_sum_fractions)

    return i_fraction_vector


def fuel_active_molar_fraction(fuel_phase, fuel_active_species_vector):
    ''' Determination molar fraction of sum of fuel active species in fuel phase object

            Type: Function

            Sub functions: -

            Input:
            fuel_phase - Fuel Cantera Phase Object
            fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]

            Output:
            fuel_active_species_molar_fraction - Molar fraction of summarized species in 'fuel phase' given by
                                                 'fuel_active_species_vector'
    '''

    fuel_active_species_molar_fraction = 0
    for counter in range(len(fuel_active_species_vector)):
        fuel_active_species_molar_fraction = fuel_active_species_molar_fraction + \
                                             fuel_phase[fuel_active_species_vector[counter]].X
    return fuel_active_species_molar_fraction[0]


def fuel_active_attribute_vector(fuel_phase, fuel_active_species):
    ''' Determination molar fraction of sum of fuel active species in fuel phase object

        Type: Function

        Sub functions: fuel_species_fractions, fuel_species_charge, i_fraction

        Input:
        fuel_phase - Fuel Cantera Phase Object
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]

        Output:
        fuel_active_attributes - List of attributes consist of list of active fuel species,
                                 list of active fuel species molar fractions,
                                 list of active fuel species charge,
                                 list of active fuel species current fractions
    '''
    fuel_active_attributes = []

    mol_fractions = fuel_species_fractions(fuel_phase, fuel_active_species)
    charge = fuel_species_charge(fuel_active_species)
    i_frac = i_fraction(fuel_active_species, mol_fractions, charge)

    fuel_active_attributes.append(fuel_active_species)
    fuel_active_attributes.append(mol_fractions)
    fuel_active_attributes.append(charge)
    fuel_active_attributes.append(i_frac)

    return fuel_active_attributes


def fuelcell_species_conversion_rates(i_fraction_absolute, reaction_type):
    ''' Determination of molar conversion of single active fuel species for given fraction of fuel cell absolute
        current

        Type: Function

        Sub functions: -

        Input:
        i_fraction_absolute - Partistion/fraction of fuel cell current for given species/reaction type
        reaction_type - reaction type of fuel active single species by complete oxidation

        Output:
        [n_species_all, n_species_anode, n_species_cathode]:

            - n_species_all - list of all moles consumed (negative) or produces (positive) by corresponding
                              absolute current fraction of given species in [kmol]
            - n_species_anode - list of all moles consumed (negative) or produces (positive) by corresponding
                                absolute current fraction of given species on anode side in [kmol]
            - n_species_cathode - list of all moles consumed (negative) or produces (positive) by corresponding
                                  absolute current fraction of given species on cathode side in [kmol]
    '''

    faraday = constants['f']

    reactants = reactions_reactants[reaction_type]

    products = reactions_products[reaction_type]

    index = reactants[0].index(reaction_type)

    fuel_species_factor = reactants[1][index]

    index = reactants[0].index("O2")

    Oxygen_factor = reactants[1][index]

    z = Oxygen_factor / fuel_species_factor * 4

    n_species_fuel_reactant = i_fraction_absolute / (z * faraday)

    n_species_oxygen_reactant = n_species_fuel_reactant * Oxygen_factor / fuel_species_factor
    n_species_anode_moles = []
    n_species_all_moles = []
    n_species_anode = []
    n_species_all = []

    for counter in range(len(products[1])):
        n_species_moles_anode = 0
        n_species_moles_all = 0
        if products[0][counter] == "O2":
            n_species_moles_anode = 0
            n_species_moles_all = n_species_fuel_reactant * (products[1][counter] - reactants[1][counter]) / 1000

        elif products[0][counter] != "O2":
            n_species_moles_anode = n_species_fuel_reactant * (products[1][counter] - reactants[1][counter]) / 1000
            n_species_moles_all = n_species_fuel_reactant * (products[1][counter] - reactants[1][counter]) / 1000

        n_species_anode_moles.append(n_species_moles_anode)
        n_species_all_moles.append(n_species_moles_all)

    n_species_anode.append(products[0])
    n_species_anode.append(n_species_anode_moles)
    n_species_all.append(products[0])
    n_species_all.append(n_species_all_moles)

    n_species_cathode = []
    n_species_cathode_name = ["O2"]
    n_species_cathode.append(n_species_cathode_name)
    n_species_cathode.append(n_species_oxygen_reactant)

    return [n_species_all, n_species_anode, n_species_cathode]


def molar_conversion(fuel_active_species, fuel_phase, i_absolute):
    ''' Determination of molar conversion of single active fuel species for given fraction of fuel cell absolute
        current

        Type: Function

        Sub functions: fuel_active_attribute_vector, i_fraction

        Input:
        fuel_phase - Fuel Cantera Phase Object
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]
        i_absolute - Absolute current in [A]


        Output:
        [anode_consumption, cathode_consumption]:

            - anode_consumption - list of all moles consumed (negative) or produces (positive) by
                                  absolute current of all active fuel species on anode side [kmol]
            - cathode_consumption - list of all moles consumed (negative) or produces (positive) of oxygen
                                    on cathode side [kmol]
    '''

    fuel_active_species_vector = fuel_active_attribute_vector(fuel_phase, fuel_active_species)[0]
    fuel_active_species_molar_fraction_vector = fuel_active_attribute_vector(fuel_phase, fuel_active_species)[1]
    fuel_active_species_charge_vector = fuel_active_attribute_vector(fuel_phase, fuel_active_species)[2]

    i_frac = i_fraction(fuel_active_species_vector, fuel_active_species_molar_fraction_vector,
                       fuel_active_species_charge_vector)

    components = ['H2', 'CO', 'CH4', 'O2', 'CO2', 'H2O', 'NH3', 'N2']

    component_names = []
    component_moles_anode = []
    component_moles_cathode = []
    anode_consumption = [components]
    cathode_consumption = [["O2"]]

    for n in range(len(components)):

        moles_components_anode = 0
        component_names.append(components[n])

        for counter in range(len(fuel_active_species_vector)):
            reaction_type = fuel_active_species_vector[counter]
            i_fraction_absolute = i_frac[counter]*i_absolute
            conv_vector = fuelcell_species_conversion_rates(i_fraction_absolute, reaction_type)

            moles_components_anode = moles_components_anode + conv_vector[1][1][n]

        component_moles_anode.append(moles_components_anode)
    anode_consumption.append(component_moles_anode)

    moles_components_cathode = 0

    for counter in range(len(fuel_active_species_vector)):
        reaction_type = fuel_active_species_vector[counter]
        i_fraction_absolute = i_frac[counter]*i_absolute
        conv_vector = fuelcell_species_conversion_rates(i_fraction_absolute, reaction_type)

        moles_components_cathode = moles_components_cathode + conv_vector[0][1][3]

    component_moles_cathode.append(moles_components_cathode)
    cathode_consumption.append(component_moles_cathode)
    return [anode_consumption, cathode_consumption]

def anode_inlet_bulk(fuel_phase, fuel_active_species_vector, i_absolute, fu):
    ''' Determination of quantity of anode bulk by given fuel_phase object, active fuel species, fuel cell current,
        and fuel utilization factor and create fuel bulk object as result

        Type: Function

        Sub functions: fuel_active_molar_fraction, molar_conversion

        Input:
        fuel_phase - Fuel Cantera Phase Object
        fuel_active_species_vector_vector - list of fuel active species in format [species_1, ... , species_n]
        i_absolute - Absolute current in [A]
        fu - Fuel utilzation Factor [-]


        Output:
        fuel_inlet_bulk - Cantera fuel bulk object generated by calculated quantity and fuel phase object, entering
                          anode inlet
    '''

    molar_fraction = fuel_active_molar_fraction(fuel_phase, fuel_active_species_vector)
    molar_conv = molar_conversion(fuel_active_species_vector, fuel_phase, i_absolute)
    molar_consumption_fuel = molar_conv[0]

    mole = 0

    for counter in range(len(fuel_active_species_vector)):
        active_species = fuel_active_species_vector[counter]
        index = molar_consumption_fuel[0].index(active_species)
        mole = mole + molar_consumption_fuel[1][index]

    factor = -mole/(fu*molar_fraction)
    fuel_inlet_bulk = ct.Quantity(fuel_phase, moles=factor)

    return fuel_inlet_bulk

def anode_outlet_bulk(anode_inlet_bulk, fuel_active_species_vector, i_absolute):
    ''' Determination of quantity of anode outlet bulk by given anode outlet bulk object, active fuel species,
        fuel cell current

        Type: Function

        Sub functions: molar_conversion

        Input:
        anode_inlet_bulk - Fuel Cantera Phase Object
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]
        i_absolute - Absolute current in [A]


        Output:
        fuel_outlet_bulk - Cantera fuel bulk object generated by calculated quantity and fuel phase object, exit
                          anode outlet
        '''

    molar_change_gas = molar_conversion(fuel_active_species_vector, anode_inlet_bulk, i_absolute)
    species_name = []
    species_mole = []
    gas = []

    for counter in range(len(anode_inlet_bulk.phase.species_names)):

        species_name.append(anode_inlet_bulk.phase.species_names[counter])
        species_mole.append(float(anode_inlet_bulk.phase.X[counter]*anode_inlet_bulk.moles))

    change_species = molar_change_gas[0][0]
    change_moles = molar_change_gas[0][1]

    for counter in range(len(species_name)):
        for j in range(len(change_species)):

            if species_name[counter] == change_species[j]:

                species_mole[counter] = float(species_mole[counter] + change_moles[j])

    gas.append(species_name)
    gas.append(species_mole)
    moles_outlet = 0

    for counter in range(len(gas[1])):
        moles_outlet = moles_outlet + gas[1][counter]

    t = anode_inlet_bulk.phase.T
    p = anode_inlet_bulk.phase.P

    outlet_gas = gas_object(gas, t, p)

    fuel_outlet_bulk = ct.Quantity(outlet_gas, moles=moles_outlet)

    return fuel_outlet_bulk


def cathode_inlet_bulk(anode_inlet_phase, cathode_inlet_phase, fuel_active_species_vector, i_absolute, airnumber):
    ''' Determination of quantity of cathode inlet bulk by given anode inlet phase object,
        cathode inlet phase object, fuel_active_species_vector, absolute fuel cell current and airnumber

        Type: Function

        Sub functions: molar_conversion

        Input:
        anode_inlet_phase - Fuel Cantera Phase Object
        cathode_inlet_phase - Air/oxidator Cantera Phase Object
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]
        i_absolute - Absolute current in [A]
        airnumber - Excess of air compared to stoechiometric oxidation process


        Output:
        cathode_inlet_bulk - Cantera air/oxidator bulk object generated by calculated quantity and cathode phase object,
                             entering cathode inlet
    '''

    molar_change_gas = molar_conversion(fuel_active_species_vector, anode_inlet_phase, i_absolute)
    moles_oxygen_min = -molar_change_gas[1][1][0]

    moles_air = moles_oxygen_min/0.21*airnumber

    cathode_inlet_bulk = ct.Quantity(cathode_inlet_phase, moles=moles_air)
    return cathode_inlet_bulk


def cathode_outlet_bulk(anode_inlet_phase, cathode_inlet_phase, fuel_active_species_vector, i_absolute, airnumber):
    ''' Determination of quantity of cathode inlet bulk by given anode inlet phase object,
        cathode inlet phase object, fuel_active_species_vector, absolute fuel cell current and airnumber

        Type: Function

        Sub functions: molar_conversion

        Input:
        anode_inlet_phase - Fuel Cantera Phase Object
        cathode_inlet_phase - Air/oxidator Cantera Phase Object
        fuel_active_species_vector - list of fuel active species in format [species_1, ... , species_n]
        i_absolute - Absolute current in [A]
        airnumber - Excess of air compared as ratio to stoechiometric oxidation process


        Output:
        cathode_outlet_bulk - Cantera air/oxidator bulk object generated by calculated quantity and cathode phase object,
                             exit cathode outlet
        '''
    air_inlet = cathode_inlet_bulk(anode_inlet_phase, cathode_inlet_phase, fuel_active_species_vector, i_absolute, airnumber)
    molar_change_gas = molar_conversion(fuel_active_species_vector, anode_inlet_phase, i_absolute)

    species_name = []
    species_mole = []
    gas = []

    for counter in range(len(air_inlet.phase.species_names)):
        species_name.append(air_inlet.phase.species_names[counter])
        species_mole.append(float(air_inlet.phase.X[counter] * air_inlet.moles))

    change_species = molar_change_gas[1][0]
    change_moles = molar_change_gas[1][1]

    for counter in range(len(species_name)):
        for j in range(len(change_species)):

            if species_name[counter] == change_species[j]:
                species_mole[counter] = float(species_mole[counter] + change_moles[j])

    gas.append(species_name)
    gas.append(species_mole)
    moles_outlet = 0

    for counter in range(len(gas[1])):
        moles_outlet = moles_outlet + gas[1][counter]

    t = cathode_inlet_phase.T
    p = cathode_inlet_phase.P

    outlet_gas = gas_object(gas, t, p)

    cathode_outlet_bulk = ct.Quantity(outlet_gas, moles=moles_outlet)

    return cathode_outlet_bulk


def oxygen_to_current(oxidator_phase, oxidator_quantity, airnumber=1, quant="mass"):
    # Function to determine absolute current [A] to convert oxygen partition given by composition of oxidator, airnumber and quantity of oxidator and electron charges of oxygen molecules

    # oxidator_phase - gas object of oxidator containing composition, partition, temperature and pressure
    # oxidator_quantity - Quantity of oxidator in kmol or kg
    # quant - quantity unit of oxidator_quantity

    # faraday - Faraday constont of 96485.33 [As/mol]
    # z_oxygen - electron charge per molecule/atom of Oxygen
    # oxidator_bulk -  Bulk of gas object containing oxygen delivered to cathode side of sofc
    # n_oxydator - moles of oxydator [mole]
    # n_oxygen - moles of oxygen in oxidatorbulk [mole]
    # i_absolute - current of converted oxygen partition in oxidator represanted by airnumber [A]

    faraday = constants['f']
    z_oxygen = 4

    if quant == "mass":

        oxidator_bulk = ct.Quantity(oxidator_phase, mass=oxidator_quantity)

    elif quant == "moles":

        oxidator_bulk = ct.Quantity(oxidator_phase, moles=oxidator_quantity)

    n_oxidator = oxidator_bulk.moles
    n_oxygen = oxidator_phase['O2'].X * n_oxidator

    i_absolute = 1000 * n_oxygen[0] * z_oxygen * faraday / airnumber

    return i_absolute, n_oxygen[0], n_oxidator


def recirculation(fuel_phase, s_to_c_ratio, fu_system_min, fu_stack, t_reformer, p_reformer):
    # Function to determine necessary recirculation rate of SOFC outlet/flue gas to reach minimal system fuel utilization by given fuel stack fuel utilization
    # or reach steam to carbon ratio for complete steam reforming of methane fraction incluting complte water-gas-shift reaction and excess of steam to enabling
    # numerical calculation of Cell Voltage by Nernst Equation - recommended s_to_c_ratio = 2.5

    # fuel phase -  fuel phase object contents composition in molar fractions of species, temperature in K, pressure in Pa etc., restricted to species CH4, H2, CO2, N2
    # s_to_c_ratio - steam to carbon ratio for complete steam reforming of methane fraction incluting complte water-gas-shift reaction and excess of steam to enabling
    # numerical calculation of Cell Voltage by Nernst Equation - recommended s_to_c_ratio = 2.5
    # fu_system_min - given minimum system fuel utlization factor
    # fu_stack - given stack fuel utlization factor

    # Minimum recirculation rate of SOFC Outlet/Flue Gas to reach SOFC system fuel utilization "fu_system_min"

    rr_min = (fu_stack - fu_system_min) / (fu_system_min * (fu_stack - 1))

    fuel_bulk_moles = 1

    # Recirculation system structured in 5 positions:
    # 1. Fuel Side upstream/before reformer
    # 2. Recirculation branch upstream/before reformer, splitted by recirculation ratio (rr) from SOFC outlet flue gas
    # 3. Reformate side downstream/after reformer and upstream/before SOFC Inlet with mix of converted methane and non converted species by
    # steam reforming of methane and water ags shift reaction from recirculation bramch
    # 4. Flue Gas downstream/after SOFC w/ electrochemical conversion of active fuel species H2 to H2O applying fuel utilization "fu_stack"
    # 5. Recirculation of SOFC flue gas

    # Defintion of variables for algebraic equation system
    h2_mole_2, co2_mole_2, n2_mole_2, h2o_mole_2 = sym.symbols("h2_mole_2 co2_mole_2 n2_mole_2 h2o_mole_2")
    h2_mole_3, co2_mole_3, n2_mole_3, h2o_mole_3 = sym.symbols("h2_mole_3 co2_mole_3 n2_mole_3 h2o_mole_3")
    h2_mole_4, co2_mole_4, n2_mole_4, h2o_mole_4, rr = sym.symbols("h2_mole_3 co2_mole_3 n2_mole_3 h2o_mole_3 rr")

    # Molar fraction of fuel species at Position 1
    ch4_mole_1 = float(fuel_phase['CH4'].X)
    h2_mole_1 = float(fuel_phase['H2'].X)
    co2_mole_1 = float(fuel_phase['CO2'].X)
    n2_mole_1 = float(fuel_phase['N2'].X)

    if ch4_mole_1 > 0:

        # Known Molar fraction of recirculation bulk given by steam to carbon ratio
        ch4_mole_2 = 0
        h2o_mole_2 = s_to_c_ratio * ch4_mole_1

        # Determination of Reformate composition expressed by varibles and conversion rates for complete steam reforming and water gas shift reaction
        ch4_mole_3 = 0
        h2o_mole_3 = h2o_mole_2 - 2 * ch4_mole_1
        h2_mole_3 = 4 * ch4_mole_1 + h2_mole_2 + h2_mole_1
        co2_mole_3 = ch4_mole_1 + co2_mole_2 + co2_mole_1
        n2_mole_3 = n2_mole_1 + n2_mole_2

        # Determination of SOFC Flue Gas/outlet composition expressed by varibles and conversion rates based on active SOFC fuel species (in this case H2) and stack fuel utilization
        ch4_mole_4 = 0
        h2o_mole_4 = fu_stack * h2_mole_3 + h2o_mole_3
        h2_mole_4 = (1 - fu_stack) * h2_mole_3
        co2_mole_4 = co2_mole_3
        n2_mole_4 = n2_mole_3

        # Recirculation of SOFC Flue Gas/outlet with recirculation rate (rr)
        ch4_mole_5 = 0
        h2o_mole_5 = rr * h2o_mole_4
        h2_mole_5 = rr * h2_mole_4
        co2_mole_5 = rr * co2_mole_4
        n2_mole_5 = rr * n2_mole_4

        # Equation system to determine variables for molar fraction in recurcalation branch and recirculation rate to keep given steam to carbon ratio (s_to_c_ratio)

        eq = [sym.Eq(h2_mole_5, h2_mole_2), sym.Eq(h2o_mole_5, h2o_mole_2), sym.Eq(co2_mole_5, co2_mole_2),
              sym.Eq(n2_mole_5, n2_mole_2)]

        # Solver for equation system to determine molar fractions H2, H2O, N2, CO2 in recircalation branch and recirculation rate to keep given steam to carbon ratio (s_to_c_ratio)
        # and fulfill boundary conditions given by single equations for closed recirculation loop

        recirculation_co2 = sym.solve(eq)[1][co2_mole_2]
        recirculation_h2 = sym.solve(eq)[1][h2_mole_2]
        recirculation_h2o = h2o_mole_2
        recirculation_n2 = sym.solve(eq)[1][n2_mole_2]
        rr = sym.solve(eq)[1][rr]

        recirculation_moles = recirculation_h2 + recirculation_h2o + recirculation_co2 + recirculation_n2

        # Determination of molar fractions and create composition vector in order to create phase objects for Recirculation branch, sofc inlet/reformer outlet branch, and SOFC outlet branch

        if rr > rr_min:

            sofc_inlet_co2 = ch4_mole_1 + recirculation_co2 + co2_mole_1
            sofc_inlet_h2 = 4 * ch4_mole_1 + recirculation_h2 + h2_mole_1
            sofc_inlet_h2o = recirculation_h2o - 2 * ch4_mole_1
            sofc_inlet_n2 = n2_mole_1 + recirculation_n2

            sofc_inlet_moles = sofc_inlet_h2 + sofc_inlet_h2o + sofc_inlet_co2 + sofc_inlet_n2

            sofc_outlet_co2 = sofc_inlet_co2
            sofc_outlet_h2 = (1 - fu_stack) * sofc_inlet_h2
            sofc_outlet_h2o = fu_stack * sofc_inlet_h2 + sofc_inlet_h2o
            sofc_outlet_n2 = sofc_inlet_n2

            recirculation_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                         [recirculation_h2, recirculation_h2o, recirculation_co2, recirculation_n2]]
            recirculation_phase = gas_object(recirculation_composition, t_reformer, p_reformer)
            recirculation_bulk = ct.Quantity(recirculation_phase, moles=recirculation_moles)

            sofc_inlet_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                      [sofc_inlet_h2, sofc_inlet_h2o, sofc_inlet_co2, sofc_inlet_n2]]
            sofc_inlet_phase = gas_object(sofc_inlet_composition, t_reformer, p_reformer)

            sofc_outlet_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                       [sofc_outlet_h2, sofc_outlet_h2o, sofc_outlet_co2, sofc_outlet_n2]]
            sofc_outlet_phase = gas_object(sofc_outlet_composition, t_reformer, p_reformer)
            sofc_inlet_bulk = ct.Quantity(sofc_inlet_phase, moles=sofc_inlet_moles)

            fu_system = sym.symbols("fu_system")

            eq_fu_system = sym.Eq((fu_stack - fu_system) / (fu_system * (fu_stack - 1)), rr)

            fu_system = sym.solve(eq_fu_system)[0]

        elif rr < rr_min:

            # Defintion of variables for algebraic equation system
            rr = rr_min
            h2_mole_2, co2_mole_2, n2_mole_2, h2o_mole_2, s_to_c_ratio = sym.symbols(
                "h2_mole_2 co2_mole_2 n2_mole_2 h2o_mole_2 s_to_c_ratio")
            h2_mole_3, co2_mole_3, n2_mole_3, h2o_mole_3 = sym.symbols("h2_mole_3 co2_mole_3 n2_mole_3 h2o_mole_3")
            h2_mole_4, co2_mole_4, n2_mole_4, h2o_mole_4 = sym.symbols("h2_mole_3 co2_mole_3 n2_mole_3 h2o_mole_3")

            # Molar fraction of fuel species at Position 1
            ch4_mole_1 = float(fuel_phase['CH4'].X)
            h2_mole_1 = float(fuel_phase['H2'].X)
            co2_mole_1 = float(fuel_phase['CO2'].X)
            n2_mole_1 = float(fuel_phase['N2'].X)

            # Known Molar fraction of recirculation bulk given by steam to carbon ratio
            ch4_mole_2 = 0
            # h2o_mole_2 = s_to_c_ratio*ch4_mole_1

            # Determination of Reformate composition expressed by varibles and conversion rates for complete steam reforming and water gas shift reaction
            ch4_mole_3 = 0
            h2o_mole_3 = h2o_mole_2 - 2 * ch4_mole_1
            h2_mole_3 = 4 * ch4_mole_1 + h2_mole_2 + h2_mole_1
            co2_mole_3 = ch4_mole_1 + co2_mole_2 + co2_mole_1
            n2_mole_3 = n2_mole_1 + n2_mole_2

            # Determination of SOFC Flue Gas/outlet composition expressed by varibles and conversion rates based on active SOFC fuel species (in this case H2) and stack fuel utilization
            ch4_mole_4 = 0
            h2o_mole_4 = fu_stack * h2_mole_3 + h2o_mole_3
            h2_mole_4 = (1 - fu_stack) * h2_mole_3
            co2_mole_4 = co2_mole_3
            n2_mole_4 = n2_mole_3

            # Recirculation of SOFC Flue Gas/outlet with recirculation rate (rr)
            ch4_mole_5 = 0
            h2o_mole_5 = rr * h2o_mole_4
            h2_mole_5 = rr * h2_mole_4
            co2_mole_5 = rr * co2_mole_4
            n2_mole_5 = rr * n2_mole_4

            # Equation system to determine variables for molar fraction in recurcalation branch and recirculation rate to keep given steam to carbon ratio (s_to_c_ratio)

            eq = [sym.Eq(h2_mole_5, h2_mole_2), sym.Eq(h2o_mole_5, h2o_mole_2), sym.Eq(co2_mole_5, co2_mole_2),
                  sym.Eq(n2_mole_5, n2_mole_2)]

            # Solver for equation system to determine molar fractions H2, H2O, N2, CO2 in recircalation branch and recirculation rate to keep given steam to carbon ratio (s_to_c_ratio)
            # and fulfill boundary conditions given by single equations for closed recirculation loop

            recirculation_co2 = sym.solve(eq)[co2_mole_2]
            recirculation_h2 = sym.solve(eq)[h2_mole_2]
            recirculation_h2o = sym.solve(eq)[h2o_mole_2]
            recirculation_n2 = sym.solve(eq)[n2_mole_2]

            recirculation_moles = recirculation_h2 + recirculation_h2o + recirculation_co2 + recirculation_n2

            s_to_c_ratio = recirculation_h2o / ch4_mole_1

            sofc_inlet_co2 = ch4_mole_1 + recirculation_co2 + co2_mole_1
            sofc_inlet_h2 = 4 * ch4_mole_1 + recirculation_h2 + h2_mole_1
            sofc_inlet_h2o = recirculation_h2o - 2 * ch4_mole_1
            sofc_inlet_n2 = n2_mole_1 + recirculation_n2

            sofc_inlet_moles = sofc_inlet_h2 + sofc_inlet_h2o + sofc_inlet_co2 + sofc_inlet_n2

            sofc_outlet_co2 = sofc_inlet_co2
            sofc_outlet_h2 = (1 - fu_stack) * sofc_inlet_h2
            sofc_outlet_h2o = fu_stack * sofc_inlet_h2 + sofc_inlet_h2o
            sofc_outlet_n2 = sofc_inlet_n2

            recirculation_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                         [recirculation_h2, recirculation_h2o, recirculation_co2, recirculation_n2]]
            recirculation_phase = gas_object(recirculation_composition, t_reformer, p_reformer)
            recirculation_bulk = ct.Quantity(recirculation_phase, moles=recirculation_moles)

            sofc_inlet_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                      [sofc_inlet_h2, sofc_inlet_h2o, sofc_inlet_co2, sofc_inlet_n2]]
            sofc_inlet_phase = gas_object(sofc_inlet_composition, t_reformer, p_reformer)
            sofc_inlet_bulk = ct.Quantity(sofc_inlet_phase, moles=sofc_inlet_moles)

            sofc_outlet_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                       [sofc_outlet_h2, sofc_outlet_h2o, sofc_outlet_co2, sofc_outlet_n2]]
            sofc_outlet_phase = gas_object(sofc_outlet_composition, t_reformer, p_reformer)

            fu_system = fu_system_min
            rr = rr_min

    elif ch4_mole_1 == 0.0:

        h2_mole_2, co2_mole_2, n2_mole_2, h2o_mole_2 = sym.symbols("h2_mole_2 co2_mole_2 n2_mole_2 h2o_mole_2")
        h2_mole_3, co2_mole_3, n2_mole_3, h2o_mole_3 = sym.symbols("h2_mole_3 co2_mole_3 n2_mole_3 h2o_mole_3")
        h2_mole_4, co2_mole_4, n2_mole_4, h2o_mole_4, rr = sym.symbols("h2_mole_3 co2_mole_3 n2_mole_3 h2o_mole_3 rr")

        # Molar fraction of fuel species at Position 1
        ch4_mole_1 = float(fuel_phase['CH4'].X)
        h2_mole_1 = float(fuel_phase['H2'].X)
        co2_mole_1 = float(fuel_phase['CO2'].X)
        n2_mole_1 = float(fuel_phase['N2'].X)
        rr = rr_min

        # h2o_mole_2 = s_to_c_ratio*ch4_mole_1

        # Determination of Reformate composition expressed by varibles and conversion rates for complete steam reforming and water gas shift reaction
        ch4_mole_3 = 0
        h2o_mole_3 = h2o_mole_2
        h2_mole_3 = h2_mole_2 + h2_mole_1
        co2_mole_3 = co2_mole_2 + co2_mole_1
        n2_mole_3 = n2_mole_1 + n2_mole_2

        # Determination of SOFC Flue Gas/outlet composition expressed by varibles and conversion rates based on active SOFC fuel species (in this case H2) and stack fuel utilization

        h2o_mole_4 = fu_stack * h2_mole_3 + h2o_mole_3
        h2_mole_4 = (1 - fu_stack) * h2_mole_3
        co2_mole_4 = co2_mole_3
        n2_mole_4 = n2_mole_3

        # Recirculation of SOFC Flue Gas/outlet with recirculation rate (rr)

        h2o_mole_5 = rr * h2o_mole_4
        h2_mole_5 = rr * h2_mole_4
        co2_mole_5 = rr * co2_mole_4
        n2_mole_5 = rr * n2_mole_4

        # Equation system to determine variables for molar fraction in recurcalation branch and recirculation rate to keep given steam to carbon ratio (s_to_c_ratio)

        eq = [sym.Eq(h2_mole_5, h2_mole_2), sym.Eq(h2o_mole_5, h2o_mole_2), sym.Eq(co2_mole_5, co2_mole_2),
              sym.Eq(n2_mole_5, n2_mole_2)]

        # Solver for equation system to determine molar fractions H2, H2O, N2, CO2 in recircalation branch and recirculation rate to keep given steam to carbon ratio (s_to_c_ratio)
        # and fulfill boundary conditions given by single equations for closed recirculation loop

        # print(eq[2])
        recirculation_co2 = sym.solve(eq)[co2_mole_2]
        # print(recirculation_co2)
        recirculation_h2 = sym.solve(eq)[h2_mole_2]
        recirculation_h2o = sym.solve(eq)[h2o_mole_2]
        recirculation_n2 = sym.solve(eq)[n2_mole_2]

        recirculation_moles = recirculation_h2 + recirculation_h2o + recirculation_co2 + recirculation_n2

        sofc_inlet_co2 = recirculation_co2 + co2_mole_1
        sofc_inlet_h2 = recirculation_h2 + h2_mole_1
        sofc_inlet_h2o = recirculation_h2o
        sofc_inlet_n2 = n2_mole_1 + recirculation_n2

        sofc_inlet_moles = sofc_inlet_h2 + sofc_inlet_h2o + sofc_inlet_co2 + sofc_inlet_n2

        sofc_outlet_co2 = sofc_inlet_co2
        sofc_outlet_h2 = (1 - fu_stack) * sofc_inlet_h2
        sofc_outlet_h2o = fu_stack * sofc_inlet_h2 + sofc_inlet_h2o
        sofc_outlet_n2 = sofc_inlet_n2

        recirculation_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                     [recirculation_h2, recirculation_h2o, recirculation_co2, recirculation_n2]]
        recirculation_phase = gas_object(recirculation_composition, t_reformer, p_reformer)
        recirculation_bulk = ct.Quantity(recirculation_phase, moles=recirculation_moles)

        sofc_inlet_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                  [sofc_inlet_h2, sofc_inlet_h2o, sofc_inlet_co2, sofc_inlet_n2]]
        sofc_inlet_phase = gas_object(sofc_inlet_composition, t_reformer, p_reformer)
        sofc_inlet_bulk = ct.Quantity(sofc_inlet_phase, moles=sofc_inlet_moles)

        sofc_outlet_composition = [['H2', 'H2O', 'CO2', 'N2'],
                                   [sofc_outlet_h2, sofc_outlet_h2o, sofc_outlet_co2, sofc_outlet_n2]]
        sofc_outlet_phase = gas_object(sofc_outlet_composition, t_reformer, p_reformer)

        fu_system = fu_system_min

    phase_vector = [recirculation_phase, sofc_inlet_phase, sofc_outlet_phase]
    parameter_vector = [s_to_c_ratio, rr, fu_system, fu_stack]

    return phase_vector, parameter_vector


def Nernst_Voltage_species(cathode_phase, anode_phase, t, p):
    # print(anode_phase['H2'].X, anode_phase['H2O'].X)
    h2_composition = [['H2'], [1]]
    co_composition = [['CO'], [1]]
    ch4_composition = [['CH4'], [1]]

    h2o_composition = [['H2O'], [1]]
    co2_composition = [['CO2'], [1]]
    o2_composition = [['O2'], [1]]

    t0 = 298.15
    p0 = 101325

    h2_phase_standard = gas_object(h2_composition, t, p0)
    co_phase_standard = gas_object(co_composition, t, p0)
    ch4_phase_standard = gas_object(ch4_composition, t, p0)

    h2o_phase_standard = gas_object(h2o_composition, t, p0)
    co2_phase_standard = gas_object(co2_composition, t, p0)
    o2_phase_standard = gas_object(o2_composition, t, p0)

    g_p0_h2 = h2_phase_standard.gibbs_mole / 1000000  # Free Enthalpy Hydrogen in kJ/mole
    g_p0_co = co_phase_standard.gibbs_mole / 1000000  # Free Enthalpy Carbonmonooxide in kJ/mole
    g_p0_ch4 = ch4_phase_standard.gibbs_mole / 1000000  # Free Enthalpy Methane in kJ/mole

    g_p0_h2o = h2o_phase_standard.gibbs_mole / 1000000  # Free Enthalpy Steam in kJ/mole
    g_p0_co2 = co2_phase_standard.gibbs_mole / 1000000  # Free Enthalpy Carbondioxide in kJ/mole
    g_p0_o2 = o2_phase_standard.gibbs_mole / 1000000  # Free Enthalpy Oxygen in kJ/mole

    coeff_h2_hydrogen = -reactions_reactants['H2'][1][0]
    coeff_o2_hydrogen = -reactions_reactants['H2'][1][3]
    coeff_h2o_hydrogen = reactions_products['H2'][1][5]

    coeff_co_carbonmonooxide = -reactions_reactants['CO'][1][1]
    coeff_o2_carbonmonooxide = -reactions_reactants['CO'][1][3]
    coeff_co2_carbonmonooxide = reactions_products['CO'][1][4]

    coeff_ch4_methane = -reactions_reactants['CH4'][1][2]
    coeff_o2_methane = -reactions_reactants['CH4'][1][3]
    coeff_co2_methane = reactions_products['CH4'][1][4]
    coeff_h20_methane = reactions_products['CH4'][1][5]

    g_p0_h2_oxidation = reactions_products['H2'][1][5] * g_p0_h2o - reactions_reactants['H2'][1][0] * g_p0_h2 - \
                        reactions_reactants['H2'][1][3] * g_p0_o2
    g_p0_co_oxidation = reactions_products['CO'][1][4] * g_p0_co2 - reactions_reactants['CO'][1][1] * g_p0_co - \
                        reactions_reactants['CO'][1][3] * g_p0_o2
    g_p0_ch4_oxidation = reactions_products['CH4'][1][4] * g_p0_co2 + reactions_products['CH4'][1][5] * g_p0_h2o - \
                         reactions_reactants['CH4'][1][2] * g_p0_ch4 - reactions_reactants['CH4'][1][3] * g_p0_o2

    f = constants["f"]
    r = constants["r"]

    p_h2 = anode_phase['H2'].X * p / p0
    p_h2o = anode_phase['H2O'].X * p / p0
    p_o2 = cathode_phase['O2'].X * p / p0

    p_co = anode_phase['CO'].X * p / p0
    p_ch4 = anode_phase['CH4'].X * p / p0
    p_co2 = anode_phase['CO2'].X * p / p0

    z_h2 = 4 * reactions_reactants['H2'][1][3]
    z_co = 4 * reactions_reactants['CO'][1][3]
    z_ch4 = 4 * reactions_reactants['CH4'][1][3]

    if p_h2 > 0 and p_o2 > 0 and p_h2o > 0:
        K_hydrogen = float(p_h2 ** coeff_h2_hydrogen * p_o2 ** coeff_o2_hydrogen * p_h2o ** coeff_h2o_hydrogen)
        u0_h2 = -g_p0_h2_oxidation / (z_h2 * f) * 1000 - r * t / (z_h2 * f) * (np.log(K_hydrogen))

    elif p_h2 == 0 or p_o2 == 0 or p_h2o == 0:
        K_hydrogen = None
        u0_h2 = None
    #############################################################################################################################

    if p_co > 0 and p_o2 > 0 and p_co2 > 0:
        K_carbonmonooxide = float(
            p_co ** coeff_co_carbonmonooxide * p_o2 ** coeff_o2_carbonmonooxide * p_co2 ** coeff_co2_carbonmonooxide)
        u0_co = -g_p0_co_oxidation / (z_co * f) * 1000 - r * t / (z_h2 * f) * np.log(K_carbonmonooxide)

    elif p_co == 0 or p_o2 == 0 or p_co2 == 0:
        K_carbonmonooxide = None
        u0_co = None

    ##############################################################################################################################
    if p_ch4 > 0 and p_o2 > 0 and p_co2 > 0 and p_h2o > 0:
        K_methane = float(
            p_ch4 ** coeff_ch4_methane * p_o2 ** coeff_o2_methane * p_h2o ** coeff_h20_methane * p_co2 ** coeff_co2_methane)
        u0_ch4 = -g_p0_ch4_oxidation / (z_ch4 * f) * 1000 - r * t / (z_ch4 * f) * np.log(K_methane)

    elif p_ch4 == 0 or p_o2 == 0 or p_co2 == 0 or p_h2o == 0:
        K_methane = None
        u0_ch4 = None

    return [u0_h2, u0_co, u0_ch4], [g_p0_h2_oxidation, g_p0_co_oxidation, g_p0_ch4_oxidation], [K_hydrogen,
                                                                                                K_carbonmonooxide,
                                                                                                K_methane]


def cell_voltage_local(cathode_phase, anode_phase, t, p, i):
    voltage_vector = Nernst_Voltage_species(cathode_phase, anode_phase, t, p)[0]

    nernst_voltage_hydrogen = voltage_vector[0]
    nernst_voltage_carbonmonooxide = voltage_vector[1]
    nernst_voltage_methane = voltage_vector[2]

    # ASR = 8,157E-9 Ohm*m² * exp(69 kJ / Rm / T)
    r = constants['r']
    asr = 10000 * 8.157 * 10 ** (-9) * np.exp(69005 / (r * t))

    if nernst_voltage_hydrogen != None and nernst_voltage_carbonmonooxide != None and nernst_voltage_methane != None:

        i_h2, i_co, i_ch4, u = sym.symbols("i_h2 i_co i_ch4 u")

        equations = [
            sym.Eq(nernst_voltage_hydrogen - i_h2 * asr, u),
            sym.Eq(nernst_voltage_carbonmonooxide - i_co * asr, u),
            sym.Eq(nernst_voltage_methane - i_ch4 * asr, u),
            sym.Eq(i_h2 + i_co + i_ch4, i)
        ]
        result = sym.solve(equations)

    elif nernst_voltage_hydrogen != None and nernst_voltage_carbonmonooxide != None and nernst_voltage_methane == None:

        i_h2, i_co, u = sym.symbols("i_h2 i_co u")

        equations = [
            sym.Eq(nernst_voltage_hydrogen - i_h2 * asr, u),
            sym.Eq(nernst_voltage_carbonmonooxide - i_co * asr, u),
            sym.Eq(i_h2 + i_co, i)
        ]
        result = sym.solve(equations)

    elif nernst_voltage_hydrogen != None and nernst_voltage_carbonmonooxide == None and nernst_voltage_methane != None:

        i_h2, i_ch4, u = sym.symbols("i_h2 i_ch4 u")

        equations = [
            sym.Eq(nernst_voltage_hydrogen - i_h2 * asr, u),
            sym.Eq(nernst_voltage_methane - i_ch4 * asr, u),
            sym.Eq(i_h2 + i_ch4, i)
        ]
        result = sym.solve(equations)

    elif nernst_voltage_hydrogen == 0 and nernst_voltage_carbonmonooxide != None and nernst_voltage_methane != None:

        i_co, i_ch4, u = sym.symbols("i_co i_ch4 u")

        equations = [
            sym.Eq(nernst_voltage_carbonmonooxide - i_co * asr, u),
            sym.Eq(nernst_voltage_methane - i_ch4 * asr, u),
            sym.Eq(i_co + i_ch4, i)
        ]
        result = sym.solve(equations)


    elif nernst_voltage_hydrogen != None and nernst_voltage_carbonmonooxide == None and nernst_voltage_methane == None:

        i_h2, u = sym.symbols("i_h2 u")

        equations = [
            sym.Eq(nernst_voltage_hydrogen - i_h2 * asr, u),
            sym.Eq(i_h2, i)
        ]
        result = sym.solve(equations)

    elif nernst_voltage_hydrogen == None and nernst_voltage_carbonmonooxide != None and nernst_voltage_methane == None:

        i_co, u = sym.symbols("i_co u")

        equations = [
            sym.Eq(nernst_voltage_carbonmonooxide - i_co * asr, u),
            sym.Eq(i_co, i)
        ]
        result = sym.solve(equations)

    elif nernst_voltage_hydrogen == None and nernst_voltage_carbonmonooxide == None and nernst_voltage_methane != None:

        i_ch4, u = sym.symbols("i_ch4 u")

        equations = [
            sym.Eq(nernst_voltage_methane - i_ch4 * asr, u),
            sym.Eq(i_ch4, i)
        ]
        result = sym.solve(equations)

    u = result[u]

    return u