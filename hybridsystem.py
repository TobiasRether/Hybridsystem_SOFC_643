# from functions import *
from components import *

# Packages to use for processing e.g. numpy for numerical processing,
# cantera for electrochemical numerical calculations etc.
# import os
# import csv
# import yaml
import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy as sy
# import scipy.optimize as so
# import xlsxwriter as Excel
# from scipy import interpolate
# import cantera as ct
# from datetime import date, time, datetime, timedelta
# from dateutil import tz
# import sympy as sym
# from sympy.solvers import solve
# from sympy import Symbol
import time


# from scipy.optimize import fsolve
# import math
# import time


def gasturbine3(__compressor_gas_in_composition, __fuel_phase, __pv1, __tv1, __mv1, __piv, __etav, __menext, __m_fuel,
                __tt1, __dp_cc, __m_cooler, __p_booster, __t_hex_out, __eta_cc, __pt2, __etat, __controller):
    """ Generic numeric Model of  Gas Turbine w/ external cooling system

        Type: Function

        Sub functions: compressor, combustion_chamber3, turbine, create_state_dataframe, cantera function

        Input:
        __compressor_gas_in_composition - Gas object composition at compressor inlet in format
        [[species_1, ... , species_n], [molar_fraction_1, ... , molar_fraction_n]]
        __fuel_phase - Fuel Combustion Chamber Inlet Cantera Phase Object
        __pv1 - Compressor inlet pressure in [Pa]
        __tv1 - Compressor inlet temperature in [K]
        __mv1 - Compressor inlet mass flow in [kg/s]
        __piv - Compressor pressure ratio [-]
        __etav - Compressor isentropic efficiency decimal unit e.g. 0.85
        __menext - Cooling air equivalent in [%]
        __m_fuel - Fuel Mass Flow [kg/s]
        __tt1 - Turbine Inlet Temperature [K]
        __dp_cc - Differential Pressure Combustion Chamber [Pa]
        __m_cooler - Mass Flow in Cooling Air Line [kg/s]
        __p_booster - Power of Cooling Air Booster [kW]
        __t_hex_out - Temperature of Cooling Air after Heat Exchanger [K]
        __eta_cc - Efficiency Combustion Chamber []
        __pt2 - Pressure Turbine Inlet [Pa]
        __etat - Efficiency Turbine Section [-]
        __controller - controller value:
            1 - Fuel Mass Flow Given Value, TT1 calculated Value
            2 - TT1 Given Value, Fuel Mass Flow calculated Value


    Output:

        bulk_vector - vector of following Cantera Objects:

            1 - __compressor_in_bulk - Cantera Object representing Compressor Inlet bulk
            2 - __compressor_in_bulk - Cantera Object representing Compressor Inlet bulk
            3 - __combustor_in_bulk - Cantera Object representing Combustor Inlet bulk
            4 - __combustor_in_bulk - Cantera Object representing Combustor Outlet bulk
            5 - fuel_bulk - Cantera Object representing Fuel Combustor Inlet bulk
            6 - cooler_in_bulk - Cantera Object representing Cooler Inlet bulk
            7 - cooler_out_bulk - Cantera Object representing Cooler Outlet bulk
            8 - __turbine_inlet_bulk - Cantera Object representing Turbine Inlet bulk
            9 - turbine_outlet_bulk - Cantera Object representing Turbine Outlet bulk

        attribute_vector:
            gasturbine_attributes:

                gasturbine_gross_power - Gross mechanical coupling power at gear [kW]
                compressor_power - Compressor mechanical power [kW]
                turbine_power - Turbine mechanical power
                gasturbine_fuel_power - Fuel Power of gas turbine [kW]
                gasturbine_gross_efficiency - Gross Power efficiency at coupling gear [-]
                compressor_pressure_ratio - Compressor pressure ratio [-]
                turbine_inlet_temperature - Turbine Inlet Temperature [K]

            compressor_section_attributes:

                __mv1 - Compressor inlet mass flow in [kg/s]
                mv1_eq - Equivalent compressor inlet mass flow in [kg/s]
                __pv1 - Compressor inlet pressure in [Pa]
                pv2 - Compressor outlet pressure in [Pa]
                __piv - Compressor pressure ratio [-]
                __tv1 - Compressor inlet temperature in [K]
                tv2 - Compressor outlet temperature in [K]
                tv2_is - Compressor outlet temperature for isentropic compression in [K]
                compressor_power - Mechanical power for given compressor mass flow in [kg/s]
                compressor_power_equivalent - - Mechanical power for equivalent compressor mass flow in [kg/s]

            combustor_attributes:

                __m_fuel - Fuel mass flow [kg/s]
                lhv - Lower heat value [kJ/kgK]
                fuel_heat_power - Fuel Heat Power [kJ]
                pt1 - Combustor inlet pressure [Pa]
                __tt1 - Combustor outlet temperature [K]

            turbine_section_attributes:

                mt1 - Turbine Inlet Mass flow [kg/s]
                pt1 - Turbine Inlet Pressure [Pa]
                tt1 - Turbine Inlet Temperature [K]
                __pt2 - Turbine Inlet Pressure [Pa]
                tt2_is - Isentropic Turbine Outlet Temperature [K]
                tt2 - Turbine Outlet Temperature [K]
                turbine_power - turbine power [kW]

        gas_properties - dataframe w/ relevant thermodynamic data for compressor inlet, compressor outlet,
                        combustor inlet, combustor outlet, fuel, cooler inlet, cooler outlet, turbine outlet,
                        turbine outlet
        """

    compressor_section = compressor(__compressor_gas_in_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)

    __compressor_in_bulk = compressor_section[0][0]
    __compressor_out_bulk = compressor_section[0][1]
    __combustor_in_bulk = compressor_section[0][1]

    df1 = compressor_section[2]

    combustor = combustion_chamber3(__compressor_in_bulk, __combustor_in_bulk, __fuel_phase, __menext, __m_cooler,
                                    __p_booster, __t_hex_out, __eta_cc, __dp_cc, __tt1, __m_fuel, __controller)

    __cooler_in_bulk = combustor[0][2]
    __cooler_out_bulk = combustor[0][3]
    __fuel_bulk = combustor[0][4]
    __combustor_out_bulk = combustor[0][5]

    __turbine_inlet_bulk = combustor[0][5]

    df2 = combustor[2]

    turbine_section = turbine(__turbine_inlet_bulk, __pt2, __etat)

    __turbine_outlet_bulk = turbine_section[0][1]

    df3 = turbine_section[2]

    bulk_vector = [__compressor_in_bulk, __compressor_out_bulk, __combustor_in_bulk, __combustor_out_bulk, __fuel_bulk,
                   __cooler_in_bulk, __cooler_out_bulk, __turbine_inlet_bulk, __turbine_outlet_bulk]

    compressor_section_attributes = compressor_section[1]
    combustor_attributes = combustor[1]
    turbine_section_attributes = turbine_section[1]

    compressor_power = compressor_section_attributes[8]
    turbine_power = turbine_section_attributes[6]
    gasturbine_gross_power = turbine_power + compressor_power
    gasturbine_fuel_power = combustor_attributes[2]
    gasturbine_gross_efficiency = gasturbine_gross_power / gasturbine_fuel_power
    compressor_pressure_ratio = compressor_section_attributes[4]
    turbine_inlet_temperature = combustor_attributes[4]

    gasturbine_attributes = [gasturbine_gross_power, compressor_power, turbine_power, gasturbine_fuel_power,
                             gasturbine_gross_efficiency, compressor_pressure_ratio, turbine_inlet_temperature]

    attribute_vector = [gasturbine_attributes, compressor_section_attributes, combustor_attributes,
                        turbine_section_attributes]

    gas_properties = pd.concat([df1, df2, df3], axis=1)

    return bulk_vector, attribute_vector, gas_properties


# help(sofc3)


def sofc_gt643_hybridsystem(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext, split_factor_bypass,
                            fuel_phase_sofc, cathode_massflow_max, t_reformer, t_sofc, airnumber, fu_stack,
                            fu_system_min, s_to_c_ratio, dp_sofc):

    compressor_inlet_bulk = compressor(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)[0][0]
    compressor_outlet_bulk = compressor(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)[0][1]

    compressor_outlet_mass = compressor_outlet_bulk.mass
    sofc_oxidator_mass = (1-split_factor_bypass)*compressor_outlet_mass
    bypass_mass = split_factor_bypass*compressor_outlet_mass

    sofc_oxidator_inlet_bulk = ct.Quantity(compressor_outlet_bulk.phase, mass=sofc_oxidator_mass)
    bypass_bulk = ct.Quantity(compressor_outlet_bulk.phase, mass=bypass_mass)

    sofc_outlet_bulk = sofc3(fuel_phase_sofc, sofc_oxidator_inlet_bulk.phase, cathode_massflow_max,
                             sofc_oxidator_inlet_bulk.mass, t_reformer, t_sofc, airnumber, fu_stack, fu_system_min,
                             s_to_c_ratio, dp_sofc)[0][11]
    return sofc_outlet_bulk


air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
fuel_composition = [['CH4'], [1]]

pv1 = 101325
tv1 = 288.15
mv1 = 1

m_fuel = 0.02
t_fuel = 300
p_fuel = 2501325

piv = 7
etav = 0.88
menext = 2.4
split_factor_bypass = 0.5
fuel_phase_sofc = gas_object(fuel_composition, t_fuel, p_fuel)
cathode_massflow_max = 1
t_reformer = 650 + 273.15
t_sofc = 850 + 273.15
airnumber = 2
fu_stack = 0.6
fu_system_min = 0.8
s_to_c_ratio = 2.05
dp_sofc = 0

result = sofc_gt643_hybridsystem(air_composition, pv1, tv1, mv1, piv, etav, menext, split_factor_bypass,
                                 fuel_phase_sofc, cathode_massflow_max, t_reformer, t_sofc, airnumber, fu_stack,
                                 fu_system_min, s_to_c_ratio, dp_sofc)

result.phase()