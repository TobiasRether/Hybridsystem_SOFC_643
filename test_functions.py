from unittest import TestCase
from functions import*


class Test(TestCase):
    def test_composition_string(self):

        gas_composition_input_format1 = [['N2', 'O2', 'AR'], [1, 1, 1]]
        gas_composition_output_format1 = 'N2:1, O2:1, AR:1'
        self.assertEqual(composition_string(gas_composition_input_format1), gas_composition_output_format1)

    def test_fuel_species_fractions(self):

        fuel_composition1 = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [2, 1, 1, 1, 1]]
        fuel_composition2 = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [1, 1, 1, 1, 1]]
        fuel_composition3 = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase1 = gas_object(fuel_composition1, t_fuel, p_fuel)
        fuel_phase2 = gas_object(fuel_composition2, t_fuel, p_fuel)
        fuel_phase3 = gas_object(fuel_composition3, t_fuel, p_fuel)
        molar_fraction1 = [0.5, 0.25, 0.25]
        molar_fraction2 = [0.333333333333333333, 0.333333333333333333333, 0.3333333333333333333333]
        molar_fraction3 = [0, 0.5, 0.5]

        fuel_active_species = ['CH4', 'H2', 'CO']
        self.assertEqual(fuel_species_fractions(fuel_phase1, fuel_active_species), molar_fraction1)
        self.assertEqual(fuel_species_fractions(fuel_phase2, fuel_active_species), molar_fraction2)
        self.assertEqual(fuel_species_fractions(fuel_phase3, fuel_active_species), molar_fraction3)

    def test_fuel_species_charge(self):

        fuel_active_species = ['CH4', 'H2', 'CO']
        charges = [['CH4', 'H2', 'CO'], [8, 2, 2]]
        self.assertEqual(fuel_species_charge(fuel_active_species), charges[1])

    def test_i_fraction(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [1, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        fuel_active_species_molar_fraction_vector = fuel_species_fractions(fuel_phase, fuel_active_species_vector)
        fuel_active_species_charge_vector = fuel_species_charge(fuel_active_species_vector)

        result = i_fraction(fuel_active_species_vector, fuel_active_species_molar_fraction_vector,
                            fuel_active_species_charge_vector)

        print(result)

    def test_fuel_active_molar_fraction(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition1 = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        fuel_composition2 = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [1, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase1 = gas_object(fuel_composition1, t_fuel, p_fuel)
        fuel_phase2 = gas_object(fuel_composition2, t_fuel, p_fuel)

        self.assertEqual(fuel_active_molar_fraction(fuel_phase1, fuel_active_species_vector), 0.5)
        self.assertEqual(fuel_active_molar_fraction(fuel_phase2, fuel_active_species_vector), 0.6)

    def test_fuel_active_attribute_vector(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        result = fuel_active_attribute_vector(fuel_phase, fuel_active_species_vector)
        print(result)

    def test_fuelcell_species_conversion_rates(self):

        i_fraction_absolute1 = 50000
        i_fraction_absolute2 = 500000
        reaction_type = 'H2'

        result1 = fuelcell_species_conversion_rates(i_fraction_absolute1, reaction_type)
        result2 = fuelcell_species_conversion_rates(i_fraction_absolute2, reaction_type)
        print(result1)
        print(result2)

    def test_molar_conversion(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        i_absolute = 50000

        result = molar_conversion(fuel_active_species_vector, fuel_phase, i_absolute)
        print(result)

    def test_anode_inlet_bulk(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        i_absolute = 500000

        fu = 0.8

        result = anode_inlet_bulk(fuel_phase, fuel_active_species_vector, i_absolute, fu)
        print(result.phase())
        print(result.moles)

    def test_anode_outlet_bulk(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        i_absolute = 500000

        fu = 0.8

        anode_in = anode_inlet_bulk(fuel_phase, fuel_active_species_vector, i_absolute, fu)

        result = anode_outlet_bulk(anode_in, fuel_active_species_vector, i_absolute)
        print(result.phase())
        print(result.moles)

    def test_cathode_inlet_bulk(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        t_air = 288.15
        p_air = 101325
        air_phase = gas_object(air_composition, t_air, p_air)

        i_absolute = 500000

        airnumber = 2

        result = cathode_inlet_bulk(fuel_phase, air_phase, fuel_active_species_vector, i_absolute, airnumber)

        print(result.phase())
        print(result.moles)

    def test_cathode_outlet_bulk(self):
        fuel_active_species_vector = ['CH4', 'H2', 'CO']

        fuel_composition = [['CH4', 'H2', 'CO', 'CO2', 'N2'], [0, 1, 1, 1, 1]]
        t_fuel = 288.15
        p_fuel = 101325
        fuel_phase = gas_object(fuel_composition, t_fuel, p_fuel)

        air_composition = [['N2', 'O2', 'AR'], [0.78, 0.21, 0.01]]
        t_air = 288.15
        p_air = 101325
        air_phase = gas_object(air_composition, t_air, p_air)

        i_absolute = 500000

        airnumber = 2

        result = cathode_outlet_bulk(fuel_phase, air_phase, fuel_active_species_vector, i_absolute, airnumber)

        print(result.phase())
        print(result.moles)