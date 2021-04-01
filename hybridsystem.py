from functions import*
from components import*

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


def gasturbine3(__gas_composition, __fuel_phase, __pv1, __tv1, __mv1, __piv, __etav, __menext, __m_fuel, __tt1, __dp_cc,
                __m_cooler, __p_booster, __t_hex_out, __eta_cc, __pt2, __etat, __controller):

    compressor_section = compressor(__gas_composition, __pv1, __tv1, __mv1, __piv, __etav, __menext)

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

    bulk_vector = [__compressor_in_bulk, __compressor_out_bulk, __combustor_in_bulk, __compressor_out_bulk, __fuel_bulk,
                   __cooler_in_bulk, __cooler_out_bulk, __turbine_inlet_bulk, __turbine_outlet_bulk]

    compressor_section_attributes = compressor_section[1]
    combustor_attributes = combustor[1]
    turbine_section_attributes = turbine_section[1]

    compressor_power = compressor_section_attributes[8]
    turbine_power = turbine_section_attributes[6]
    gasturbine_gross_power = turbine_power + compressor_power
    gasturbine_fuel_power = combustor_attributes[2]
    compressor_pressure_ratio = compressor_section_attributes[4]
    turbine_inlet_temperature = combustor_attributes[4]

    gasturbine_attributes = [gasturbine_gross_power, compressor_power, turbine_power, gasturbine_fuel_power]

    attribute_vector = [gasturbine_attributes, compressor_section_attributes, combustor_attributes,
                        turbine_section_attributes, compressor_pressure_ratio, turbine_inlet_temperature]

    gas_properties = pd.concat([df1, df2, df3], axis=1)

    return bulk_vector, attribute_vector, gas_properties

# compressor input


__gas_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
__pv1 = 101325
__tv1 = 288.15
__mv1 = 191
__piv = 16.5
__etav = 0.88
__menext = 2.4

# combustor input

fuel_composition = [['CH4'], [1]]
__m_fuel = 3.0
t_fuel = 300
p_fuel = 2501325

__fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

__tt1 = 1130 + 273.15
__dp_cc = 20000

__m_cooler = 16.74
__p_booster = 386
__t_hex_out = 198 + 273.15
__eta_cc = 0.998

__controller = 2

# turbine input

__pt2 = 103325
__etat = 0.88


result = gasturbine3(__gas_composition, __fuel_phase, __pv1, __tv1, __mv1, __piv, __etav, __menext, __m_fuel, __tt1,
                     __dp_cc, __m_cooler, __p_booster, __t_hex_out, __eta_cc, __pt2, __etat, __controller)

print(result[1])
