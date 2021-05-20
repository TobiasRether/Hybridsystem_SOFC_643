from functions import *
from components import *
import csv

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
                            fu_system_min, s_to_c_ratio, dp_sofc, __m_cooler, __p_booster,__t_hex_out,fuel_phase_gt,
                            __eta_cc, __dp_cc,__tt1, __m_fuel, __controller, __pt2, __etat):

    compressor_inlet_bulk = compressor(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)[0][0]
    compressor_outlet_bulk = compressor(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)[0][1]

    compressor_outlet_mass = compressor_outlet_bulk.mass

    cooler_inlet_mass = (1-split_factor_bypass)*compressor_outlet_mass

    if cooler_inlet_mass <= 20:
        cooler_inlet_mass = cooler_inlet_mass
    elif cooler_inlet_mass >20:
        m_add_comb = cooler_inlet_mass-20
        cooler_inlet_mass = 20

    cooler_inlet_bulk = ct.Quantity(compressor_outlet_bulk.phase, mass=cooler_inlet_mass)
    cooler_outlet_bulk = cooler(cooler_inlet_bulk)[1]

    if cooler_outlet_bulk.mass > 15.4:
        m_add_comb_2 = cooler_outlet_bulk.mass - 15.4
        cooler_outlet_bulk.mass = 15.4
    elif cooler_outlet_bulk.mass <= 15.4:
        cooler_outlet_bulk.mass = cooler_outlet_bulk.mass

    p2_b = cooler_outlet_bulk.phase.P * 1.22
    eta_booster = 0.85

    booster_outlet_bulk = booster(cooler_outlet_bulk,p2_b,eta_booster)[0]

    sofc_oxidator_mass = booster_outlet_bulk.mass
    bypass_mass = split_factor_bypass*compressor_outlet_mass + m_add_comb + m_add_comb_2

    sofc_oxidator_inlet_bulk = ct.Quantity(booster_outlet_bulk.phase, mass=sofc_oxidator_mass)
    bypass_bulk = ct.Quantity(compressor_outlet_bulk.phase, mass=bypass_mass)

    sofc_outlet_bulk = sofc3(fuel_phase_sofc, sofc_oxidator_inlet_bulk.phase, cathode_massflow_max,
                             sofc_oxidator_inlet_bulk.mass, t_reformer, t_sofc, airnumber, fu_stack, fu_system_min,
                             s_to_c_ratio, dp_sofc)[0][11]
    combustor_inlet_oxidator_bulk = sofc_outlet_bulk + bypass_bulk
    combustor_outlet_flue_bulk = combustion_chamber3(compressor_inlet_bulk,combustor_inlet_oxidator_bulk, fuel_phase_gt,
                                                     __menext, __m_cooler, __p_booster,__t_hex_out, __eta_cc, __dp_cc,
                                                     __tt1, __m_fuel, __controller)[0][5]
    turbine_outlet_bulk = turbine(combustor_outlet_flue_bulk, __pt2, __etat)[0][1]



    gas_prop_comp = compressor(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)[2]
    gas_prop_cool = cooler(cooler_inlet_bulk)[3]
    gas_prop_boo = booster(cooler_outlet_bulk, p2_b, eta_booster)[2]
    gas_prop_comb = combustion_chamber3(compressor_inlet_bulk,combustor_inlet_oxidator_bulk, fuel_phase_gt,
                                                     __menext, __m_cooler, __p_booster,__t_hex_out, __eta_cc, __dp_cc,
                                                     __tt1, __m_fuel, __controller)[2]
    gas_prop_turb = turbine(combustor_outlet_flue_bulk, __pt2, __etat)[2]
    gas_prop_sofc = sofc3(fuel_phase_sofc, sofc_oxidator_inlet_bulk.phase, cathode_massflow_max,
                             sofc_oxidator_inlet_bulk.mass, t_reformer, t_sofc, airnumber, fu_stack, fu_system_min,
                             s_to_c_ratio, dp_sofc)[5]
    df_list = [gas_prop_comp, gas_prop_cool, gas_prop_boo, gas_prop_comb, gas_prop_turb, gas_prop_sofc]
    multiple_dfs(df_list, 'Bulk','Results.xlsx',1)

    attr_comp = pd.DataFrame (compressor(__air_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)[1],
                              index = ['MV1','MV1 EQ', 'PV1', 'PV2','PIV','TV1','TV2', 'TV2 IS','Compressor Power',
                                       'Compressor Power Equivalent'])
    attr_cool = pd.DataFrame(cooler(cooler_inlet_bulk)[2], index = ['M1', 'T2', 'P2'])
    attr_comb = pd.DataFrame(combustion_chamber3(compressor_inlet_bulk,combustor_inlet_oxidator_bulk, fuel_phase_gt,
                                                     __menext, __m_cooler, __p_booster,__t_hex_out, __eta_cc, __dp_cc,
                                                     __tt1, __m_fuel, __controller)[1], index=['M Fuel', 'LHV',
                                                                                               'Fuel Heat Power',
                                                                                               'PT1','TT1'])
    attr_turb = pd.DataFrame(turbine(combustor_outlet_flue_bulk, __pt2, __etat)[1], index = ['MT1', 'PT1', 'TT1', 'PT2',
                             'TT2 IS', 'TT2', 'Turbine Power'])
    attr_boo = pd.DataFrame(booster(cooler_outlet_bulk, p2_b, eta_booster)[1], index = ['M1','T2','T2 IS',
                                                                                        'Booster Power'])
    df_list_1 = [attr_comp,attr_cool, attr_boo, attr_comb, attr_turb]
    multiple_dfs_1(df_list_1, 'Attribute Vector', 'Results.xlsx', 1)

    perf_sofc = pd.DataFrame(sofc3(fuel_phase_sofc, sofc_oxidator_inlet_bulk.phase, cathode_massflow_max,
                             sofc_oxidator_inlet_bulk.mass, t_reformer, t_sofc, airnumber, fu_stack, fu_system_min,
                             s_to_c_ratio, dp_sofc)[2], index=['p_absolute', 'i_absolute', 'fuel_add_power',
                                                               'reformate_power','eta_add_fuel', 'eta_reformate',
                                                               'lhv_fuel', 'lhv_reformate'])
    heat_sofc = pd.DataFrame(sofc3(fuel_phase_sofc, sofc_oxidator_inlet_bulk.phase, cathode_massflow_max,
                             sofc_oxidator_inlet_bulk.mass, t_reformer, t_sofc, airnumber, fu_stack, fu_system_min,
                             s_to_c_ratio, dp_sofc)[1], index=['q_cathode', 'q_anode_fuel', 'q_reformer',
                                                               'q_steam_reforming[0]', 'q_reformer_anode',
                                                               'q_sofc_ohmic_losses','q_losses', 'q_deficite'])
    elec_sofc = pd.DataFrame(sofc3(fuel_phase_sofc, sofc_oxidator_inlet_bulk.phase, cathode_massflow_max,
                             sofc_oxidator_inlet_bulk.mass, t_reformer, t_sofc, airnumber, fu_stack, fu_system_min,
                             s_to_c_ratio, dp_sofc)[3], index=['u', 'u0', 'u_inlet', 'u_outlet', 'u0_inlet',
                                                               'u0_outlet', 'i_max', 'nominal_factor', 'i_nominal',
                                                               'i_absolute_max', 'i_absolute', 'i', 'a'])
    df_list_2 = [perf_sofc, heat_sofc, elec_sofc]
    multiple_dfs_1(df_list_2, 'SOFC', 'Results.xlsx', 1)

    eta_inverter = 0.94
    eta_gear = 0.99
    eta_geno_gt = 0.985
    list(attr_comp.columns.values)
    gross_turbo_comp_power = attr_comp.at['Compressor Power',attr_comp.columns[0]]+attr_turb.at['Turbine Power',
                                                                                                attr_turb.columns[0]]
    netto_turbo_comp_power = gross_turbo_comp_power * eta_gear * eta_geno_gt
    booster_power = attr_boo.at['Booster Power', attr_boo.columns[0]]
    gross_sofc_elec_power = perf_sofc.at['p_absolute',perf_sofc.columns[0]]
    netto_sofc_elec_power = gross_sofc_elec_power * eta_inverter
    sofc_fuel_power = perf_sofc.at['fuel_add_power',perf_sofc.columns[0]]
    comb_fuel_power = attr_comb.at['Fuel Heat Power', attr_comb.columns[0]]
    gross_eta_system = (gross_turbo_comp_power+booster_power+gross_sofc_elec_power)/(sofc_fuel_power+comb_fuel_power)
    netto_eta_system = (netto_turbo_comp_power+booster_power+netto_sofc_elec_power)/(sofc_fuel_power+comb_fuel_power)
    efficiencies = pd.DataFrame([eta_inverter, eta_gear, eta_geno_gt, gross_eta_system,netto_eta_system],
                                index=['Eta Inverter', 'Eta Gear','Eta Geno GT','Gross Eta System',
                                       'Netto Eta System'])
    power = pd.DataFrame ([gross_turbo_comp_power,netto_turbo_comp_power,booster_power, gross_sofc_elec_power,
                           netto_sofc_elec_power,sofc_fuel_power, comb_fuel_power], index=['Gross Turbo Components '
                                                                                           'Power','Netto Turbo'
                                                                                                   'Components Power',
                                                                                           'Booster Power',
                                                                                           'Gross SOFC Elec Power',
                                                                                           'Netto SOFC Power',
                                                                                           'SOFC Fuel Power',
                                                                                           'Combustion Chamber Fuel'
                                                                                           'Power'])
    df_list_3 = [efficiencies, power]
    multiple_dfs_1(df_list_3, 'Power', 'Results.xlsx',1)

    return compressor_inlet_bulk, turbine_outlet_bulk


#help(combustion_chamber3)
air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
fuel_composition = [['CH4'], [1]]

pv1 = 101325
tv1 = 288.15
mv1 = 189

m_fuel = 0.02
t_fuel = 300
p_fuel = 2501325

piv = 16
etav = 0.88
menext = 2.4
split_factor_bypass = 0.5
fuel_phase_sofc = gas_object(fuel_composition, t_fuel, p_fuel)
cathode_massflow_max = 180
t_reformer = 650 + 273.15
t_sofc = 850 + 273.15
airnumber = 2
fu_stack = 0.6
fu_system_min = 0.8
s_to_c_ratio = 2.05
dp_sofc = 200000
__m_cooler = 15
__p_booster = 600
__t_hex_out = 200 + 273.15
fuel_phase_gt = gas_object(fuel_composition, t_fuel, p_fuel)
__eta_cc = 0.998
__dp_cc = 200
__tt1 = 850 + 273.15
__m_fuel = 0.02
__controller = 2
__pt2 = pv1+ 2000
__etat = 0.88
fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

result = sofc_gt643_hybridsystem(air_composition, pv1, tv1, mv1, piv, etav, menext, split_factor_bypass,
                                 fuel_phase_sofc, cathode_massflow_max, t_reformer, t_sofc, airnumber, fu_stack,
                                 fu_system_min, s_to_c_ratio, dp_sofc,__m_cooler, __p_booster,__t_hex_out,fuel_phase_gt,
                                 __eta_cc, __dp_cc,__tt1, __m_fuel, __controller, __pt2, __etat)

#result[0].phase()
#result[1].phase()

#help(gasturbine3)

#result = gasturbine3(air_composition,fuel_phase,pv1,tv1,mv1,piv,etav,menext,__m_fuel,__tt1,__dp_cc,__m_cooler,
#                     __p_booster,__t_hex_out, __eta_cc, __pt2, __etat, __controller)

#print(result[1][0][4])

#help(sofc3)