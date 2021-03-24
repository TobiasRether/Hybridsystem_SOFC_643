from functions import*

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
import cantera as ct
# from datetime import date, time, datetime, timedelta
# from dateutil import tz
import sympy as sym
from sympy.solvers import solve
from sympy import Symbol
# from scipy.optimize import fsolve
# import math
# import time


def compressor_2(gas, p1, t1, m1, pi, is_eff, menext=0):
    # gas - gas composition at inlet as tuple for species and fraction of species
    # p1 - pressure at compressor inlet
    # t1 - temperature at compressor inlet
    # m1 -  mass flow of gas at compressor inlet
    # is_eff - isentropic efficiency
    # menext - cooling air equivalent in [%]

    # create_state_dataframe(gas_bulk, position) to file gas properties @ assigned positions

    # Calculation of equivalent cooling air mass flow

    m_eq = Symbol('m_eq')

    md = menext / 100

    eq1 = sym.Eq((m1 / m_eq - 1), md)

    m1_eq = solve(eq1, m_eq)[0]

    ########################################################################################
    gas = gas_object(gas, t1, p1)

    compressor_inlet_bulk = ct.Quantity(gas, mass=m1)
    # t1 = compressor_inlet_bulk.phase.T
    # p1 = compressor_inlet_bulk.phase.P
    # m1 = compressor_inlet_bulk.mass
    df1 = create_state_dataframe(compressor_inlet_bulk, "Compressor Inlet")

    p2 = pi * p1

    s1 = compressor_inlet_bulk.phase.s
    h1 = compressor_inlet_bulk.phase.h

    compressor_inlet_bulk.SP = s1, p2
    t2_is = compressor_inlet_bulk.phase.T
    h2_is = compressor_inlet_bulk.phase.h

    h2 = (h2_is - h1) / is_eff + h1

    compressor_inlet_bulk.HP = h2, p2
    t2 = compressor_inlet_bulk.phase.T

    compressor_outlet_bulk = ct.Quantity(compressor_inlet_bulk.phase, mass=compressor_inlet_bulk.mass)

    df2 = create_state_dataframe(compressor_outlet_bulk, "Compressor Outlet")

    gas_properties = pd.concat([df1, df2], axis=1)

    compressor_power = - m1 * (h2 - h1) / 1000
    compressor_power_equivalent = - m1_eq * (h2 - h1) / 1000
    # print(h2-h1)
    # print(compressor_power, compressor_power_equivalent)
    # print(((164010.125002+compressor_power_equivalent)/164010.125002)*100)

    return compressor_outlet_bulk, m1, p2, t2_is, t2, compressor_power, gas_properties


gas_composition = [['N2'], [1]]
tv1 = 288.15
pv1 = 101325
mv1 = 1
pi = 7
is_eff = 0.85
menext = 0

test = compressor_2(gas_composition, pv1, tv1, mv1, pi, is_eff, menext)
