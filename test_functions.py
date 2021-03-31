from unittest import TestCase
from functions import*


class Test(TestCase):
    def test_composition_string(self):

        gas_composition_input_format1 = [['N2', 'O2', 'AR'], [1, 1, 1]]
        gas_composition_output_format1 = 'N2:1, O2:1, AR:1'
        self.assertEqual(composition_string(gas_composition_input_format1), gas_composition_output_format1)

