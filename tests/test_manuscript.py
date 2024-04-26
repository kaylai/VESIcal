import unittest
import VESIcal as v
import numpy as np
import pathlib
import pandas as pd
from pandas.testing import assert_frame_equal


# class TestSplitWeight(unittest.TestCase):
#     """
#     A class to define a method for asserting two pandas dataframes are equivalent.
#     From: https://stackoverflow.com/questions/38839402/how-to-use-assert-frame-equal-in-unittest
#     """
#     def assertDataframeEqual(self, a, b, msg):
#         try:
#             pd_testing.assert_frame_equal(a, b)
#         except AssertionError as e:
#             raise self.failureException(msg) from e

#     def setUp(self):
#         self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

#     def test_allZero(self):
#         self.assertEqual(pd.DataFrame([0,0,0,0]), pd.DataFrame([0,0,0,0]))

def print_msg_box(msg, indent=1, width=None, title=None):
    """Print message-box with optional title."""
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    print("\n")
    print(box)

# Allow unittest to find the file
TEST_FILE = pathlib.Path(__file__).parent.joinpath("ManuscriptBatchFile.xlsx")
PICKLE_DISSOLVED = pathlib.Path(__file__).parent.joinpath("manuscript_dissolved.p")
PICKLE_EQFLUID = pathlib.Path(__file__).parent.joinpath("manuscript_eqfluid.p")
PICKLE_EQFLUID_WTEMPS = pathlib.Path(__file__).parent.joinpath("manuscript_eqfluid_wTemps.p")
PICKLE_EQFLUID_WT = pathlib.Path(__file__).parent.joinpath("manuscript_eqfluid_wt.p")
PICKLE_EQFLUID_MOL = pathlib.Path(__file__).parent.joinpath("manuscript_eqfluid_mol.p")
PICKLE_SATPS = pathlib.Path(__file__).parent.joinpath("manuscript_satPs.p")
PICKLE_SATPS_WTEMPS = pathlib.Path(__file__).parent.joinpath("manuscript_satPs_wTemps.p")
PICKLE_ISOBARS = pathlib.Path(__file__).parent.joinpath("manuscript_isobars.p")
PICKLE_ISOPLETHS = pathlib.Path(__file__).parent.joinpath("manuscript_isopleths.p")
PICKLE_CLOSED_DF = pathlib.Path(__file__).parent.joinpath("manuscript_closed.p")
PICKLE_OPEN_DF = pathlib.Path(__file__).parent.joinpath("manuscript_open_df.p")
PICKLE_HALF_DF = pathlib.Path(__file__).parent.joinpath("manuscript_half_df.p")
PICKLE_EXSOLVED_DF = pathlib.Path(__file__).parent.joinpath("manuscript_exsolved_df.p")
PICKLE_START2000_DF = pathlib.Path(__file__).parent.joinpath("manuscript_start2000_df.p")    

class TestManuscriptCalculations(unittest.TestCase):
    """
    Some of the "correct" values hard coded here are not exactly the same as those in the
    manuscript. However, they are very very close and likely are different because of changes
    to normalization routines made after publication. I've noted original values in commented
    text below where they differ from currently accepted values.

    Current values last updated March 2024 by Kayla Iacovino. 
    Some degassing paths updated March 2024 by Simon Matthews (see comments below).
    """
    def setUp(self):
        # import any pickled objects
        self.pickle_dissolved = pd.read_pickle(PICKLE_DISSOLVED)
        self.pickle_eqfluid = pd.read_pickle(PICKLE_EQFLUID)
        self.pickle_eqfluid_wTemps = pd.read_pickle(PICKLE_EQFLUID_WTEMPS)
        self.pickle_eqfluid_wt = pd.read_pickle(PICKLE_EQFLUID_WT)
        self.pickle_eqfluid_mol = pd.read_pickle(PICKLE_EQFLUID_MOL)
        self.pickle_satPs = pd.read_pickle(PICKLE_SATPS)
        self.pickle_satPs_wTemps = pd.read_pickle(PICKLE_SATPS_WTEMPS)
        self.pickle_isobars = pd.read_pickle(PICKLE_ISOBARS)
        self.pickle_isopleths = pd.read_pickle(PICKLE_ISOPLETHS)
        self.pickle_closed_df = pd.read_pickle(PICKLE_CLOSED_DF)
        self.pickle_open_df = pd.read_pickle(PICKLE_OPEN_DF)
        self.pickle_half_df = pd.read_pickle(PICKLE_HALF_DF)
        self.pickle_exsolved_df = pd.read_pickle(PICKLE_EXSOLVED_DF)
        self.pickle_start2000_df = pd.read_pickle(PICKLE_START2000_DF)

        # example data file
        self.myfile = v.BatchFile(TEST_FILE)

        # mysample used in ms
        self.mysample = v.Sample({'SiO2':  77.3,
                                  'TiO2':   0.08, 
                                  'Al2O3': 12.6,
                                  'Fe2O3':  0.207,
                                  'Cr2O3':  0.0,
                                  'FeO':    0.473,
                                  'MnO':    0.0,
                                  'MgO':    0.03,
                                  'NiO':    0.0,
                                  'CoO':    0.0,
                                  'CaO':    0.43,
                                  'Na2O':   3.98,
                                  'K2O':    4.88,
                                  'P2O5':   0.0,
                                  'H2O':    6.5,
                                  'CO2':    0.05})
        
        # values are same as in ms unless noted!
        self.mysample_satP = {'SaturationP_bars': 2960.0,  # changed from 2960 in ms due to rounding
                              'FluidMass_grams': 0.0018160337487088,
                              'FluidProportion_wt': 0.0018160337487087978,
                              'XH2O_fl': 0.838064480487942,
                              'XCO2_fl': 0.161935519512058}

        # mybasalt sample used in ms
        self.mybasalt = v.Sample({'SiO2': 47, 
                    'TiO2': 1.01,  
                    'Al2O3': 17.46,
                    'Fe2O3': 0.89,
                    'FeO': 7.18,
                    'MgO': 7.63,
                    'CaO': 12.44,
                    'Na2O': 2.65,
                    'K2O': 0.03,
                    'P2O5': 0.08, 
                    'CO2': 0.1})
        
        # normalize in three ways
        self.mybasalt_std = self.mybasalt.get_composition(normalization="standard",
                                                          asSampleClass=True)
        self.mybasalt_add = self.mybasalt.get_composition(normalization="additionalvolatiles",
                                                          asSampleClass=True)
        self.mybasalt_fix = self.mybasalt.get_composition(normalization="fixedvolatiles",
                                                          asSampleClass=True)
        
        self.manuscript_vals_norm_test = [1848.03, 1906.55, 1848.27, 1848.26]
        """
        Actual values from manuscript:
        Raw Saturation Pressure = 1848.031831425599 bars
        standard Saturation Pressure = 1906.5453789627868 bars
        additionalvolatiles Saturation Pressure = 1848.2673972122493 bars
        fixedvolatiles Saturation Pressure = 1848.2611364359402 bars
        """
        
        self.sample_10star = v.Sample({'SiO2':      47.9600,
                                        'TiO2':     0.7800,
                                        'Al2O3':    18.7700,
                                        'Fe2O3':    0.0000,
                                        'Cr2O3':    0.0000,
                                        'FeO':      10.9200,
                                        'MnO':      0.1500,
                                        'MgO':      6.8600,
                                        'NiO':      0.0000,
                                        'CoO':      0.0000,
                                        'CaO':      12.2300,
                                        'Na2O':     1.9500,
                                        'K2O':      0.2100,
                                        'P2O5':     0.1700,
                                        'H2O':      4.5000,
                                        'CO2':      0.0479 })
        
        self.tenstar_verbose_diss_vol = {'H2O_liq': 2.69337030386574,
                                        'CO2_liq': 0.0638470697446765,
                                        'XH2O_fl': 0.500053772677827,
                                        'XCO2_fl': 0.499946227322173,
                                        'FluidProportion_wt': 0.19679860317100492}
        """
        Actual values from manuscript:
        {'H2O_liq': 2.69352739399806,
        'CO2_liq': 0.0638439414375309,
        'XH2O_fl': 0.500092686493868,
        'XCO2_fl': 0.499907313506132,
        'FluidProportion_wt': 0.18407321260435108}
        """

        self.tenstar_eqfluid = {'CO2': 0.0052346121751066,
                                'H2O': 0.994765387824893}
        """
        Actual values from manuscript:
        {'CO2': 0.00528661429366132,
         'H2O': 0.994713385706339}
        """
    
    def test_degassing_paths(self):
        print_msg_box("TestManuscript \ndegassing_paths")
        """Calculate open, closed, and closed + 2 wt% initial vapor"""
        closed_df = v.calculate_degassing_path(sample=self.sample_10star, temperature=1200).result
        assert_frame_equal(closed_df, self.pickle_closed_df, atol=1e-4)

        open_df = v.calculate_degassing_path(sample=self.sample_10star, temperature=1200,
                                             fractionate_vapor=1.0).result
        assert_frame_equal(open_df, self.pickle_open_df, atol=1e-4)

        half_df = v.calculate_degassing_path(sample=self.sample_10star, temperature=1200,
                                             fractionate_vapor=0.5).result
        assert_frame_equal(half_df, self.pickle_half_df, atol=5e-4)

        exsolved_df = v.calculate_degassing_path(sample=self.sample_10star, temperature=1200,
                                                 init_vapor=2.0).result
        assert_frame_equal(exsolved_df, self.pickle_exsolved_df, atol=1e-4)

        """Calculate closed-system degassing starting from a pressure of 2000 bars"""
        start2000_df = v.calculate_degassing_path(sample=self.sample_10star, temperature=1200,
                                                  pressure=2000.0).result
        assert_frame_equal(start2000_df, self.pickle_start2000_df, atol=1e-4)
    
    def test_dissolved_batch(self):
        print_msg_box("TestManuscript \ndissolved_batch")
        result = self.myfile.calculate_dissolved_volatiles(temperature=900.0, pressure=2000.0,
                                                           X_fluid=1, print_status=True)
        assert_frame_equal(result, self.pickle_dissolved, atol=1e-4)
    
    def test_dissolved_volatiles_sample_10star(self):
        print_msg_box("TestManuscript \ndissolved_volatiles_sample_10star")
        result = v.calculate_dissolved_volatiles(sample=self.sample_10star, temperature=900.0,
                                                 pressure=2000.0, X_fluid=0.5, verbose=True).result
        params = list(result.keys())
        for param in params:
            self.assertAlmostEqual(result[param], self.tenstar_verbose_diss_vol[param], places=4)
    
    def test_eqfluid_batch(self):
        print_msg_box("TestManuscript \neqfluid_batch")
        result = self.myfile.calculate_equilibrium_fluid_comp(temperature=900.0, pressure=1000.0)
        assert_frame_equal(result, self.pickle_eqfluid, atol=1e-4)
    
    def test_eqfluid_batch_wTemps(self):
        print_msg_box("TestManuscript \neqfluid_batch_wTemps")
        result = self.myfile.calculate_equilibrium_fluid_comp(temperature='Temp', pressure='Press')
        assert_frame_equal(result, self.pickle_eqfluid_wTemps, atol=1e-4)
    
    def test_eqfluid_molfrac_to_wt(self):
        print_msg_box("TestManuscript \neqfluid_molfrac_to_wt")
        result = v.fluid_molfrac_to_wt(self.pickle_eqfluid)
        assert_frame_equal(result, self.pickle_eqfluid_wt, atol=1e-4)
    
    def test_eqfluid_wt_to_molfrac(self):
        print_msg_box("TestManuscript \neqfluid_wt_to_molfrac")
        result = v.fluid_wt_to_molfrac(self.pickle_eqfluid_wt)
        assert_frame_equal(result, self.pickle_eqfluid_mol, atol=1e-4)
    
    def test_equilibrium_fluid_sample_10star(self):
        print_msg_box("TestManuscript \nequilibrium_fluid_sample_10star")
        result = v.calculate_equilibrium_fluid_comp(sample=self.sample_10star, temperature=900.0,
                                                    pressure=100.0).result
        self.assertAlmostEqual(result['H2O'], self.tenstar_eqfluid['H2O'], places=4)
        self.assertAlmostEqual(result['CO2'], self.tenstar_eqfluid['CO2'], places=4)
    
    def test_isobars_and_isopleths(self):
        print_msg_box("TestManuscript \nisobars_and_isopleths")
        result_isobars, result_isopleths = v.calculate_isobars_and_isopleths(
                                            sample=self.sample_10star, 
                                            temperature=1200.0,
                                            pressure_list=[1000.0, 2000.0, 3000.0],
                                            isopleth_list=[0.25,0.5,0.75]).result
        assert_frame_equal(result_isobars, self.pickle_isobars, atol=1e-4)
        assert_frame_equal(result_isopleths, self.pickle_isopleths, atol=1e-4)
    
    def test_manuscript_plots(self):
        """
        This simply tests that plots created in the manuscript don't throw an error.
        """
        print_msg_box("TestManuscript \ntest_manuscript_plots")
        v.calib_plot()
        v.plot(isobars=self.pickle_isobars, isopleths=self.pickle_isopleths)
        v.plot(degassing_paths=[self.pickle_open_df, self.pickle_half_df, self.pickle_closed_df,
                                self.pickle_exsolved_df],
               degassing_path_labels=["Open", "Half", "Closed", "Exsolved"])
        v.plot(degassing_paths=[self.pickle_exsolved_df, self.pickle_start2000_df],
               degassing_path_labels=["Exsolved", "2000 bars"])
        v.plot(isobars=self.pickle_isobars, 
               isopleths=self.pickle_isopleths, 
               degassing_paths=[self.pickle_open_df, self.pickle_closed_df], 
               degassing_path_labels=["Open System", "Closed System"])

    def test_satP_different_normalizations(self):
        print_msg_box("TestManuscript \nsatP_different_normalization")
        for basalt, msval in zip([self.mybasalt, self.mybasalt_std, self.mybasalt_add,
                                  self.mybasalt_fix], self.manuscript_vals_norm_test):
            result = v.calculate_saturation_pressure(sample=basalt, temperature=1200,
                                                     model="IaconoMarziano").result
            self.assertAlmostEqual(result, msval, places=2)
    
    def test_saturation_pressure_batch(self):
        print_msg_box("TestManuscript \nsaturation_pressure_batch")
        result = self.myfile.calculate_saturation_pressure(temperature=925.0)
        assert_frame_equal(result, self.pickle_satPs, atol=1e-4)
    
    def test_saturation_pressure_batch_wTemps(self):
        print_msg_box("TestManuscript \nsaturation_pressure_batch_wTemps")
        result = self.myfile.calculate_saturation_pressure(temperature="Temp")
        assert_frame_equal(result, self.pickle_satPs_wTemps, atol=1e-4)
    
    def test_saturation_pressure_mysample(self):
        print_msg_box("TestManuscript \nsaturation_pressure_mysample")
        result = v.calculate_saturation_pressure(sample=self.mysample, temperature=925.0,
                                                 verbose=True).result
        params = list(result.keys())
        for param in params:
            self.assertAlmostEqual(result[param], self.mysample_satP[param], places=4)
