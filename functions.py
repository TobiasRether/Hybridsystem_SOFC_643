__author__ = "Tobias Rether"
__copyright__ = "Copyright 2021, Siemens Energy"
__department__ = "SE GP G SV TI EN RW PRS COE"
__email__ = "tobias.rether@siemens-energy.com"
__status__ = "Demo"


# import plotly.offline as py '
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# Packages to use for processing e.g. numpy for numerical processing,
# cantera for electrochemical numerical calculations etc.
# import os
# import csv

import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy as sy
# import scipy.optimize as so
# import xlsxwriter as Excel
import cantera as ct
# from datetime import date, time, datetime, timedelta
# from dateutil import tz
# import sympy as sym
# from sympy.solvers import solve
# from sympy import Symbol
# from scipy.optimize import fsolve
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

        Subfunctions: composition_string

        Input:
        gas_composition - gas object composition in format
        [[species_1, ... , species_n], [molar_fraction_1, ... , molar_fraction_n]]
        t - Temperature in K
        p - pressure in Pa

        Output:
        gas_phase - Cantera gas phase object calculated w/ NASA7 Polynoms from gri30 input file
                    based on given composition, temperature and pressure from Input
                    :rtype: object
        '''
    gas_phase = ct.Solution('gri30.xml')
    components = composition_string(gas_composition)
    gas_phase.TPX = t, p, components

    return gas_phase


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


def append_dataframe(dataframe1, dataframe2):
    df = pd.concat([dataframe1, dataframe2], axis=1)

    return df