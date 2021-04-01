from functions import*

# Packages to use for processing e.g. numpy for numerical processing,
# cantera for electrochemical numerical calculations etc.
# import os
# import csv
# import yaml
import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
# import scipy as sy
# import scipy.optimize as so
# import xlsxwriter as Excel
from scipy import interpolate
import cantera as ct
# from datetime import date, time, datetime, timedelta
# from dateutil import tz
import sympy as sym
from sympy.solvers import solve
from sympy import Symbol
import time
# from scipy.optimize import fsolve
# import math
# import time

def compressor(__gas_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext=0):
    """ Generic numeric Model of compressor parametrized by input values
        __gas_composition, __pv1, __tv1, __mv1, __piv, __menext

    Type: Function

    Sub functions: gas_object, create_state_dataframe, cantera function

    Input:
    __gas_composition_input_format - Gas object composition in format
    [[species_1, ... , species_n], [molar_fraction_1, ... , molar_fraction_n]]
    __pv1 - Compressor inlet pressure in [Pa]
    __tv1 - Compressor inlet temperature in [K]
    __mv1 - Compressor inlet mass flow in [kg/s]
    __piv - Compressor pressure ratio [-]
    __etav - Compressor isentropic efficiency decimal unit e.g. 0.85
    __menext - Cooling air equivalent in [%]

    Output:

    bulk_vector:

        compressor_inlet_bulk - Cantera gas quantity object based on phase object and quantity scalar (e.g. kg) for
                                medium at compressor inlet
        compressor_outlet_bulk - Cantera gas quantity object based on phase object and quantity scalar (e.g. kg) for
                                 medium at compressor inlet

    attribute_vector:
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

    gas_properties - dataframe w/ relevant thermodynamic data for compressor inlet, compressor outlet
    """

    m_eq = Symbol('m_eq')

    md = __menext / 100

    eq1 = sym.Eq((__mv1 / m_eq - 1), md)

    mv1_eq = solve(eq1, m_eq)[0]

    gas_phase = gas_object(__gas_composition, __tv1, __pv1)

    compressor_inlet_bulk = ct.Quantity(gas_phase, mass=__mv1)

    compressor_outlet_bulk = ct.Quantity(gas_phase, mass=__mv1)

    pv2 = __piv * __pv1
    sv1 = compressor_outlet_bulk.phase.s
    hv1 = compressor_outlet_bulk.phase.h

    compressor_outlet_bulk.SP = sv1, pv2
    tv2_is = compressor_outlet_bulk.phase.T
    hv2_is = compressor_outlet_bulk.phase.h
    hv2 = (hv2_is - hv1) / __etav + hv1

    compressor_outlet_bulk.HP = hv2, pv2
    tv2 = compressor_outlet_bulk.phase.T
    compressor_power = - __mv1 * (hv2 - hv1) / 1000
    compressor_power_equivalent = - mv1_eq * (hv2 - hv1) / 1000

    bulk_vector = [compressor_inlet_bulk, compressor_outlet_bulk]

    attribute_vector = [__mv1, mv1_eq, __pv1, pv2, __piv, __tv1, tv2, tv2_is, compressor_power, compressor_power_equivalent]

    df1 = create_state_dataframe(compressor_inlet_bulk, "Compressor Inlet")
    df2 = create_state_dataframe(compressor_outlet_bulk, "Compressor Outlet")
    gas_properties = pd.concat([df1, df2], axis=1)

    return bulk_vector, attribute_vector, gas_properties


def combustion_chamber3(__compressor_in_bulk, __combustor_in_bulk, __fuel_phase, __menext, __m_cooler, __p_booster,
                        __t_hex_out, __eta_cc, __dp_cc, __tt1, __m_fuel, __controller):

    """ Generic numeric Model of combustion chamber for Gas Turbine w/ external ccoling systen
     parametrized by input values __compressor_in_bulk, __combustor_in_bulk, __fuel_phase, __menext, __m_cooler, __p_booster,
                        __t_hex_out, __eta_cc, __dp_cc, __tt1, __m_fuel, __controller

    Type: Function

    Sub functions: gas_object, sensible_enthalpy, fuel_to_air_ratio2, object_to_composition, heating_value_mix
    create_state_dataframe, cantera function

    Input:
    __compressor_in_bulk - Compressor Inlet Cantera Quantity Object
    __combustor_in_bulk - Compressor Outlet/Oxidizer Combustion Chamber Inlet Cantera Quantity Object
    __fuel_phase - Fuel Combustion Chamber Inlet Cantera Phase Object
    __menext - Cooling air equivalent in [%]
    __m_cooler - Mass Flow in Cooling Air Line [kg/s]
    __p_booster - Power of Cooling Air Booster [kW]
    __t_hex_out - Temperature of Cooling Air after Heat Exchanger [K]
    __eta_cc - Efficiency Combustion Chamber []
    __dp_cc - Differential Pressure Combustion Chamber [Pa]
    __tt1 - Turbine Inlet Temperature [K]
    __m_fuel - Fuel Mass Flow [kg/s]
    __controller - controller value:
        1 - Fuel Mass Flow Given Value, TT1 calculated Value
        2 - TT1 Given Value, Fuel Mass Flow calculated Value


    Output:

    bulk_vector - vector of following Cantera Objects:

        1 - __compressor_in_bulk - Cantera Object representing Compressor Inlet bulk
        2 - __combustor_in_bulk - Cantera Object representing Compressor Outlet bulk
        3 - cooler_in_bulk - Cantera Object representing Cooler Inlet bulk
        4 - cooler_out_bulk - Cantera Object representing Cooler Outlet bulk
        5 - fuel_bulk - Cantera Object representing Fuel Combustor Inlet bulk
        6 - flue_gas_bulk - Cantera Object representing Flue Combustor Outlet bulk

    attribute_vector:
        __m_fuel - Fuel mass flow [kg/s]
        lhv - Lower heat value [kJ/kgK]
        fuel_heat_power - Fuel Heat Power [kJ]
        pt1 - Combustor inlet pressure [Pa]
        __tt1 - Combustor outlet temperature [K]

    gas_properties - dataframe w/ relevant thermodynamic data for compressor inlet, compressor outlet
    """

    hv1 = sensible_enthalpy(__compressor_in_bulk.phase, ref_type=1)[0]/1000
    mv1 = __compressor_in_bulk.mass

    hv2 = sensible_enthalpy(__combustor_in_bulk.phase, ref_type=1)[0]/1000
    mv2 = __combustor_in_bulk.mass
    pv2 = __combustor_in_bulk.phase.P

    cooler_in_bulk = ct.Quantity(__combustor_in_bulk.phase, mass=__m_cooler)
    cooler_out_bulk = ct.Quantity(__combustor_in_bulk.phase, mass=__m_cooler)
    cooler_out_bulk.TP = __t_hex_out, pv2

    h_hex_out = sensible_enthalpy(cooler_out_bulk.phase, ref_type=1)[0]/1000

    phi = 1

    fuel_to_air_ratio_stoech = fuel_to_air_ratio2(__combustor_in_bulk.phase, __fuel_phase, phi)
    __m_fuel_stoechiometric = float(mv2)*float(fuel_to_air_ratio_stoech)

    fuel_bulk = ct.Quantity(__fuel_phase, mass=__m_fuel_stoechiometric)

    enthalpy = []
    temperature = []
    n = 100

    for i in range(n+1):

        fuel_bulk.mass = __m_fuel_stoechiometric*i/n
        flue_gas_bulk = __combustor_in_bulk + fuel_bulk
        flue_gas_bulk.equilibrate('HP')
        enthalpy.append(sensible_enthalpy(flue_gas_bulk.phase, ref_type=1)[0]/1000)
        temperature.append(flue_gas_bulk.phase.T)

    interpolation_t_h_flue_gas = interpolate.interp1d(temperature, enthalpy)
    interpolation_h_t_flue_gas = interpolate.interp1d(enthalpy, temperature)

    m_eq = Symbol('m_eq')

    md = __menext / 100

    eq1 = sym.Eq((mv1 / m_eq - 1), md)

    mv1_eq = solve(eq1, m_eq)[0]

    m_ext = 0
    m_w = 0
    h_w = 0

    h_sens_fuel = sensible_enthalpy(__fuel_phase, ref_type=1)[0]/1000

    fuel_composition = object_to_composition(__fuel_phase)
    lhv = heating_value_mix(fuel_composition)[0]

    flue_gas_bulk = 'self'

    if __controller == 1:

        fuel_bulk.mass = __m_fuel
        ht1 = Symbol('ht1')

        eq3 = sym.Eq(((mv1 * hv1) + mv1_eq * (hv2 - hv1) - m_ext * hv2 + m_w * h_w - __m_cooler * (hv2 - h_hex_out) +
                      __m_fuel * __eta_cc * (lhv + h_sens_fuel) + __p_booster) / (mv1 + __m_fuel + m_w), ht1)

        ht1_calc = float(solve(eq3, ht1)[0])

        __tt1 = interpolation_h_t_flue_gas(ht1_calc)

        flue_gas_bulk = __combustor_in_bulk + fuel_bulk
        flue_gas_bulk.equilibrate('HP')

        pt1 = pv2 - __dp_cc
        flue_gas_bulk.TP = __tt1, pt1

    elif __controller == 2:

        __m_fuel = Symbol('__m_fuel')

        ht1 = interpolation_t_h_flue_gas(__tt1)

        eq2 = sym.Eq(((mv1*hv1) + mv1_eq*(hv2-hv1) - m_ext*hv2 + m_w * h_w - __m_cooler * (hv2 - h_hex_out) +
                      __m_fuel*__eta_cc*(lhv + h_sens_fuel) + __p_booster)/(mv1 + __m_fuel + m_w), ht1)

        __m_fuel_calc = solve(eq2, __m_fuel)[0]
        fuel_bulk.mass = __m_fuel_calc
        flue_gas_bulk = __combustor_in_bulk + fuel_bulk
        flue_gas_bulk.equilibrate('HP')

        pt1 = pv2 - __dp_cc
        flue_gas_bulk.TP = __tt1, pt1

    __m_fuel = fuel_bulk.mass

    fuel_heat_power = lhv*__m_fuel

    bulk_vector = [__compressor_in_bulk, __combustor_in_bulk, cooler_in_bulk, cooler_out_bulk, fuel_bulk, flue_gas_bulk]

    attribute_vector = [__m_fuel, lhv, fuel_heat_power, pt1, __tt1]

    df1 = create_state_dataframe(__combustor_in_bulk, "Oxidator Combustion Chamber Inlet")
    df2 = create_state_dataframe(cooler_in_bulk, "Cooler Inlet")
    df3 = create_state_dataframe(cooler_out_bulk, "Cooler Outlet")
    df4 = create_state_dataframe(fuel_bulk, "Fuel Combustion Chamber Inlet")
    df5 = create_state_dataframe(flue_gas_bulk, "Flue Gas Combustion Chamber Outlet")
    gas_properties = pd.concat([df1, df2, df3, df4, df5], axis=1)

    return bulk_vector, attribute_vector, gas_properties


def combustion_2(pv1, tv1, mv1, __dp_cc, eta_cc, compressor_outlet_phase, fuel_phase, menext=0, m_fuel=None, tt1=None,
                 phi=None, controller=None):
    t1 = time.time()

    # compressor_outlet_phase - compressor outlet phase object contents compostion in molar fractions of species, temperature in K, pressure in Pa etc.
    # fuel_phase - fuel phase object contents composition in molar fractions of species, temperature in K, pressure in Pa etc.
    # mv1 - Mass Flow of oxidator enters compressor inlet in kg/s
    # __dp_cc - Combustion Chamber pressure loss in bar
    # eta_cc - Combustion chamber efficiency by normalization of radiation/heat losses
    # m_fuel - Fuel Mass Flow (only needed if Controller is set to 1 - Fuel mass controller)
    # tt1 - Combustion Chamber Exhaust Tempertaure (only needed if Controller is set to 2 - Exhaust temperature controller)
    # phi - Combustion Chamber airnumeber/equivalent ratio control (only needed if Controller is set to 3 - Phi controller)
    #
    # controller - control principle to determine combustion chamber operation point:
    # None - No controller select
    # 1 - Mass Flow Controller -  determination of operation point by given m_fuel, t3 and phi will be ignored
    # 2 - Temperature Exhaust Controller - determination of operation point by given combustion outlet temperature t3 as target value
    # 3 - Phi Controller - determination of operation point by given equivalence (air to fuel mass ratio) phi

    if compressor_outlet_phase == None:
        print("Bitte Oxidator als Cantera Phase-Objekt (oxidator_phase) angeben")
        print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")
    if fuel_phase == None:
        print("Bitte Brennstoff als Cantera Phase-Objekt (fuel_phase) angeben")
        print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")
    if mv1 == None:
        print("Bitte Oxidatormassenstrom als Fließkommazahl oder Integer (m_oxi) in kg/s angeben")
        print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")
    if m_fuel == None:
        pass
    if __dp_cc == None:
        print("Bitte Brennkammerdruckverlust als Fließkommazahl oder Integer (__dp_cc) in mbar angeben")
        print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")
    if controller == None:
        print("Bitte Regelung (controller) als Integer angeben:\n"
              "1 - Mass Flow Controller -  determination of operation point by given m_fuel, t3 and phi will be ignored\n"
              "2 - Temperature Exhaust Controller - determination of operation point by given combustion outlet temperature t3 as target value\n"
              "3 - Phi Controller - determination of operation point by given equivalence (air to fuel mass ratio) phi")
        print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")

    if controller == 2 and tt1 == None:
        print("Bitte Brennkammeraustrittstemperatur t3 bei Controllerstellung 2 in kg/s angeben")
        print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")

    fuel_composition = object_to_composition(fuel_phase)
    compressor_outlet_composition = object_to_composition(compressor_outlet_phase)

    compressor_inlet_phase = gas_object(compressor_outlet_composition, tv1, pv1)

    tv2 = compressor_outlet_phase.T
    pv2 = compressor_outlet_phase.P
    pt1 = pv2 - __dp_cc

    hv1 = sensible_enthalpy(compressor_inlet_phase, ref_type=1)[0] / 1000
    hv2 = sensible_enthalpy(compressor_outlet_phase, ref_type=1)[0] / 1000
    h_sens_fuel = sensible_enthalpy(fuel_phase, ref_type=1)[0] / 1000

    m_eq = Symbol('m_eq')
    md = menext / 100

    eq1 = sym.Eq((mv1 / m_eq - 1), md)
    mv1_eq = solve(eq1, m_eq)[0]

    h_water = 0
    __m_cooler = 0
    __p_booster = 0

    m_water = 0
    m_ex_discharge = 0
    m_leakage = 0

    m_fuel_target_value = m_fuel
    phi_target_value = phi
    t_fuel = fuel_phase.T
    p_fuel = fuel_phase.P

    h_cooler_out = 0

    lhv = heating_value_mix(fuel_composition)[0]
    # lhv=50000

    # identify limits of combustion chamber operation tempoperature limits by Inlet Temperature as minimum and stoeciometric operation point as maximum

    # Minimum Temperature at phi=0

    tt1_min = tv2

    # Maximum Temperature at phi=1 (stoechiometric combustion)
    # calculating w/ fuel mass flow at phi=1
    phi = 1
    m_fuel_stoech = fuel_to_air_ratio(compressor_outlet_composition, fuel_composition, phi) * mv1_eq

    t2 = time.time()

    compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)

    fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

    compressor_outlet_bulk = ct.Quantity(compressor_outlet_phase, mass=mv1_eq)

    fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel_stoech)

    stoechiometric_mix = compressor_outlet_bulk + fuel_bulk
    stoechiometric_mix.equilibrate('HP')
    tt1_max = stoechiometric_mix.phase.T

    #################################################################################################
    # Calculation of exhaust flue gas properties based on fuel mass flow target value
    if controller == 1:
        # print(mv1)
        if m_fuel == None:
            print("Bitte Brennstoffmassenstrom m_fuel bei Controllerstellung 1 in kg/s angeben")
            print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")

        compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        compressor_outlet_bulk = ct.Quantity(compressor_outlet_phase, mass=mv1)
        fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel_target_value)

        flue_gas_bulk = compressor_outlet_bulk + fuel_bulk
        phi = flue_gas_bulk.get_equivalence_ratio()
        flue_gas_bulk.equilibrate('HP')

        tt1 = flue_gas_bulk.TP[0]
        flue_gas_bulk.TP = tt1, pt1

    t3 = time.time()
    #################################################################################################
    # Calculation of exhaust flue gas properties based on reverse air numer/equivalent ratio target value
    if controller == 3:

        if phi == None:
            print("Bitte Äquivalenzverhältnis phi bei Controllerstellung 3 angeben")
            print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")

        compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        m_fuel_target_value = phi_target_value * m_fuel_stoech
        # print(phi_target_value, m_fuel_stoech)
        compressor_outlet_bulk = ct.Quantity(oxidator_phase, mass=mv1_eq)
        fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel_target_value)

        flue_gas_bulk = oxidator_bulk + fuel_bulk
        phi = flue_gas_bulk.get_equivalence_ratio()
        flue_gas_bulk.equilibrate('HP')

        tt1 = flue_gas_bulk.TP[0]
        flue_gas_bulk.TP = tt1, pt1

    ##################################################################################################
    t4 = time.time()

    # Calculation of exhaust flue gas properties based on temperature target value
    if controller == 2:

        if tt1 == None:
            print("Bitte Brennkammeraustrittstemperatur t3 bei Controllerstellung 2 in kg/s angeben")
            print("combustion(oxidator_phase, fuel_phase, m_oxi, __dp_cc, m_fuel, t3, phi, controller)")

        if tt1 < tt1_min:
            print("Temperaturzielwert Brennkammeraustritt ist kleiner als Brennkammerientrittstemperatur")

        if tt1 > tt1_max:
            print(
                "Temperaturzielwert Brennkammeraustritt ist größer maximal erreichbare Temperatur bei stöchiometrischer Verbrennung")

        if tt1_min < tt1 < tt1_max:

            ht1 = Symbol('ht1')
            m_fuel = Symbol('m_fuel')

            mv2 = mv1 - m_ex_discharge

            m_fuel_ini = 0.35 * m_fuel_stoech
            compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)

            compressor_outlet_bulk_eq = ct.Quantity(compressor_outlet_phase, mass=mv1_eq)
            fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel_ini)
            # print(tt1, pt1)
            flue_gas_bulk_eq = compressor_outlet_bulk_eq + fuel_bulk
            flue_gas_bulk_eq.equilibrate('HP')
            flue_gas_bulk_eq.TP = tt1, pt1
            ht1 = sensible_enthalpy(flue_gas_bulk_eq.phase, ref_type=1)[0] / 1000

            delta_h = 1

            # print('mv1: ' + str(mv1))
            # print('mv1_eq: ' + str(mv1_eq))
            # print('pv1: ' + str(pv1))
            # print('tv1: ' + str(tv1))
            # print('hv1: ' + str(hv1))

            # print('tv2: ' + str(tv2))
            # print('pv2: ' + str(pv2))
            # print('hv2: ' + str(hv2))

            # print('lhv: ' + str(lhv))
            # print('h_sens_fuel: ' + str(h_sens_fuel))

            # print('eta_cc: ' + str(eta_cc))

            while delta_h > 0.1:
                compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)
                fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

                energy_balance_cc = sym.Eq(
                    (mv1 + m_fuel + m_water - m_ex_discharge - m_leakage) * ht1 + hv2 * __m_cooler + (1 - eta_cc) * (
                                lhv + h_sens_fuel) * m_fuel, (mv1_eq - m_ex_discharge) * hv2 + (
                                lhv + h_sens_fuel) * m_fuel + h_cooler_out * __m_cooler + h_water * m_water + (
                                mv1 - mv1_eq) * hv1 + __p_booster)

                m_fuel_ini = solve(energy_balance_cc, m_fuel)[0]
                # print(m_fuel_ini)
                compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)
                compressor_outlet_bulk_eq = ct.Quantity(compressor_outlet_phase, mass=mv1_eq)
                fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel_ini)

                flue_gas_bulk = compressor_outlet_bulk_eq + fuel_bulk
                flue_gas_bulk.equilibrate('HP')
                flue_gas_bulk.TP = tt1, pt1
                ht1_j = sensible_enthalpy(flue_gas_bulk.phase, ref_type=1)[0] / 1000

                delta_h = abs(ht1_j - ht1)
                ht1 = ht1_j

            # print('ht1: ' + str(ht1))
            m_fuel_target_value = m_fuel_ini
            fuel_bulk = ct.Quantity(fuel_phase, mass=m_fuel_target_value)
            compressor_outlet_phase = gas_object(compressor_outlet_composition, tv2, pv2)
            compressor_outlet_bulk = ct.Quantity(compressor_outlet_phase, mass=mv2)
            flue_gas_bulk = compressor_outlet_bulk + fuel_bulk
            flue_gas_bulk.equilibrate('HP')
            flue_gas_bulk.TP = tt1, pt1

    tt1 = flue_gas_bulk.phase.T
    pt1 = flue_gas_bulk.phase.P
    mt1 = flue_gas_bulk.mass

    t4 = time.time()

    df1 = create_state_dataframe(fuel_bulk, "Fuel")
    df2 = create_state_dataframe(compressor_outlet_bulk, "Combustion Chamber Inlet")
    df3 = create_state_dataframe(flue_gas_bulk, "Combustion Chamber Outlet")
    gas_properties = pd.concat([df1, df2, df3], axis=1)

    ht1 = sensible_enthalpy(flue_gas_bulk.phase, ref_type=1)[0] / 1000

    heat_power = fuel_bulk.mass * lhv
    t5 = time.time()

    return flue_gas_bulk, tt1, pt1, mt1, fuel_bulk.mass, lhv, heat_power, gas_properties


def turbine(__turbine_inlet_bulk, __pt2, __etat):
    """ Generic numeric Model of turbine for Gas Turbine parametrized by input values
        __turbine_inlet_bulk, __pt2, __etat

     Type: Function

     Sub functions: create_state_dataframe, cantera function

     Input:
     __turbine_inlet_bulk - Turbine Inlet Cantera Quantity Object
     __pt2 - Pressure Turbine Inlet [Pa]
     __etat - Efficiency Turbine Section []

     Output:

     bulk_vector - vector of following Cantera Objects:

         1 - __turbine_inlet_bulk - Cantera Object representing Turbine Inlet bulk
         2 - turbine_outlet_bulk - Cantera Object representing Turbine Outlet bulk

     attribute_vector:

         mt1 - Turbine Inlet Mass flow [kg/s]
         pt1 - Turbine Inlet Pressure [Pa]
         tt1 - Turbine Inlet Temperature [K]
         __pt2 - Turbine Inlet Pressure [Pa]
         tt2_is - Isentropic Turbine Outlet Temperature [K]
         tt2 - Turbine Outlet Temperature [K]
         turbine_power - turbine power [kW]

     gas_properties - dataframe w/ relevant thermodynamic data for compressor inlet, compressor outlet
     """

    tt1 = __turbine_inlet_bulk.phase.T
    pt1 = __turbine_inlet_bulk.phase.P
    mt1 = __turbine_inlet_bulk.mass
    st1 = __turbine_inlet_bulk.phase.s
    ht1 = __turbine_inlet_bulk.phase.h

    turbine_outlet_bulk = ct.Quantity(__turbine_inlet_bulk.phase, mass=mt1)
    turbine_outlet_bulk.SP = st1, __pt2
    tt2_is = turbine_outlet_bulk.phase.T
    ht2_is = turbine_outlet_bulk.phase.h

    ht2 = (ht2_is - ht1) * __etat + ht1

    turbine_outlet_bulk.HP = ht2, __pt2
    tt2 = turbine_outlet_bulk.phase.T

    turbine_power = - mt1 * (ht2 - ht1) / 1000

    df1 = create_state_dataframe(__turbine_inlet_bulk, "Turbine Inlet")
    df2 = create_state_dataframe(turbine_outlet_bulk, "Turbine Outlet")

    bulk_vector = [__turbine_inlet_bulk, turbine_outlet_bulk]

    attribute_vector = [mt1, pt1, tt1, __pt2, tt2_is, tt2, turbine_power]

    gas_properties = pd.concat([df1, df2], axis=1)

    return bulk_vector, attribute_vector, gas_properties