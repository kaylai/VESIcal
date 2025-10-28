import unittest
import VESIcal as v
import pandas as pd
import pathlib

# Allow unittest to find the file
TEST_FILE = pathlib.Path(__file__).parent.joinpath("ImportTest.xlsx")

class TestImportExcel(unittest.TestCase):
    def assertDataframeEqual(self, a, b, msg):
        """
        Creates a new type of unittest to assert that pd DataFrames
        are equal, inheriting from pandas testing routine
        """
        try:
            pd._testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        # Add assertDataframeEqual to unittest test cases
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

        self.data_dict = {'SiO2':  [47.95, 47.95],
                     'TiO2':  [1.67, 1.67],
                     'Al2O3': [17.32, 17.32],
                     'FeO':   [10.24, 10.24],
                     'Fe2O3': [0.1, 0.1],
                     'MgO':   [5.76, 5.76],
                     'CaO':   [10.93, 10.93],
                     'Na2O':  [3.45, 3.45],
                     'K2O':   [1.99, 1.99],
                     'P2O5':  [0.51, 0.51],
                     'MnO':   [0.1, 0.1],
                     'H2O':   [2.0, 2.0],
                     'CO2':   [0.1, 0.1],
                     'Notes': ['Normal sample', 'Duplicate sample']}
        
        self.df = pd.DataFrame(self.data_dict, 
                               index=['test_samp',
                                      'test_samp-duplicate-1'])

        self.myfile = v.BatchFile(TEST_FILE)

    def test_ImportExcel(self):
        self.assertEqual(self.df, self.myfile.get_data(), 
                         'DataFrames are different')