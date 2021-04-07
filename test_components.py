from unittest import TestCase
from components import*
import time


class Test(TestCase):
    def test_compressor(self):

        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        pv1 = 101325
        tv1 = 288.15
        mv1 = 1
        piv = 7
        etav = 0.88
        menext = 2.4

        time_start = time.time()
        result = compressor(air_composition, pv1, tv1, mv1, piv, etav, menext)
        time_end = time.time()
        print(time_end - time_start)
        print(result[0][0].phase.report())
        print(result[0][1].phase.report())

    def test_combustion_chamber3(self):
        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        fuel_composition = [['CH4'], [1]]

        t_compressor_in = 15 + 273.15
        p_compressor_in = 101325

        t_compressor_out = 404 + 273.15
        p_compressor_out = 1604000

        m_fuel = 3.47
        t_fuel = 300
        p_fuel = 2501325

        compressor_in_phase = gas_object(air_composition, t_compressor_in, p_compressor_in)
        compressor_out_phase = gas_object(air_composition, t_compressor_out, p_compressor_out)
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        m_compressor_out = 187.4
        m_compressor_in = 187.4
        tt1 = 1130 + 273.15
        dp = 20000

        menext = 2.36

        m_cooler = 16.74
        p_booster = 386
        t_hex_out = 198 + 273.15
        eta_cc = 0.998

        compressor_in_bulk = ct.Quantity(compressor_in_phase, mass=m_compressor_in)
        compressor_out_bulk = ct.Quantity(compressor_out_phase, mass=m_compressor_out)

        controller = 1

        time_start = time.time()
        result = combustion_chamber3(compressor_in_bulk, compressor_out_bulk, fuel_phase, menext, m_cooler, p_booster,
                                     t_hex_out, eta_cc, dp, tt1, m_fuel, controller)
        time_end = time.time()

        delta_time = time_end - time_start

        print(delta_time)

    def test_combustion_2(self):

        pv1 = 101325
        tv1 = 288.15
        mv1 = 180
        pv2 = 1601325
        tv2 = 688.15
        dp_cc = 20000
        eta_cc = 99.8

        t_fuel = 288.15
        p_fuel = 101325
        m_fuel = 3.6

        menext = 2.4

        tt1 = 1130 + 273.15
        phi = 0.35
        controller = 1

        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        fuel_composition = [['CH4'], [1]]

        compressor_outlet_phase = gas_object(air_composition, tv2, pv2)
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        time_start = time.time()
        result = combustion_2(pv1, tv1, mv1, dp_cc, eta_cc, compressor_outlet_phase, fuel_phase, menext, m_fuel, tt1,
                     phi, controller)
        time_end = time.time()

        print(time_end-time_start)

    def test_turbine(self):
        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        fuel_composition = [['CH4'], [1]]

        t_compressor_in = 15 + 273.15
        p_compressor_in = 101325

        t_compressor_out = 404 + 273.15
        p_compressor_out = 1604000

        m_fuel = 3.47
        t_fuel = 300
        p_fuel = 2501325

        compressor_in_phase = gas_object(air_composition, t_compressor_in, p_compressor_in)
        compressor_out_phase = gas_object(air_composition, t_compressor_out, p_compressor_out)
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        m_compressor_out = 187.4
        m_compressor_in = 187.4
        tt1 = 1130 + 273.15
        dp = 20000

        menext = 2.36

        m_cooler = 16.74
        p_booster = 386
        t_hex_out = 198 + 273.15
        eta_cc = 0.998

        compressor_in_bulk = ct.Quantity(compressor_in_phase, mass=m_compressor_in)
        compressor_out_bulk = ct.Quantity(compressor_out_phase, mass=m_compressor_out)

        controller = 2

        combustor = combustion_chamber3(compressor_in_bulk, compressor_out_bulk, fuel_phase, menext, m_cooler, p_booster,
                                     t_hex_out, eta_cc, dp, tt1, m_fuel, controller)

        __etat = 0.88
        __pt2 = 101325
        __turbine_inlet_bulk = combustor[0][5]

        result = turbine(__turbine_inlet_bulk, __pt2, __etat)

        print(result[1])

    def test_sofc3(self):
        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        fuel_composition = [['CH4'], [1]]

        t_compressor_out = 404 + 273.15
        p_compressor_out = 1604000

        cathode_massflow_max = 20
        cathode_massflow = 10

        t_fuel = 300
        p_fuel = 2501325

        compressor_out_phase = gas_object(air_composition, t_compressor_out, p_compressor_out)
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        t_reformer = 650 + 273.15
        t_sofc = 850 + 273.15

        airnumber = 3
        fu_stack = 0.6

        fu_system_min = 0.8
        s_to_c_ratio =2.05

        dp_sofc = 100000

        result = sofc3(fuel_phase, compressor_out_phase, cathode_massflow_max, cathode_massflow, t_reformer,
                       t_sofc, airnumber, fu_stack, fu_system_min, s_to_c_ratio, dp_sofc)

        print(result)