import unittest
import pandas as pd
import numpy as np
import tempfile
import os

from VESIcal import sample_class, batchfile, calculate_classes, utils

class MockCalculate(calculate_classes.Calculate):
    def __init__(self, result):
        self.result = result

class TestSaveResults(unittest.TestCase):

    def setUp(self):
        # scalar
        self.scalar = 42
        self.np_scalar = np.float64(3.14159)
        self.calc_scalar = MockCalculate(result=42)
        
        # dict
        self.dict_obj = {"SiO2": 50.2, "MgO": 7.5}
        self.calc_dict = MockCalculate(result={"A": 1, "B": 2})
        
        # tuple
        self.tup_obj = (44, 45, 46)
        
        # pandas objects
        self.df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        self.series = pd.Series(self.dict_obj)
        self.calc_df = MockCalculate(result=pd.DataFrame({"col": [10, 20]}))
        self.calc_series = MockCalculate(result=pd.Series({"Z": 99}))
        
        # VESIcal objects
        # create temporary place to store written files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = self.temp_dir.name
        
        # Create a temporary CSV file to import to BatchFile object
        df = pd.DataFrame({
            'Lable': ['Test1'],
            'SiO2': [72.4],
            'Al2O3': [13.1],
            'FeO': [2.4]
        })
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_filename = f.name        
        
        # create objects of supported types for testing
        # VESIcal objects
        self.sample = sample_class.Sample({"SiO2": 70, "Al2O3": 15})
        self.batch = batchfile.BatchFile(temp_filename, label="Sample")
        self.calc_batch = MockCalculate(result=self.batch)
        self.calc_sample = MockCalculate(result=self.sample)

        os.remove(temp_filename)
                
        # mixed list
        self.mixed_data = [
            self.scalar,
            self.np_scalar,
            self.dict_obj,
            self.tup_obj,
            self.df,
            self.series,
            self.batch,
            self.sample,
            self.calc_scalar,
            self.calc_dict,
            self.calc_df,
            self.calc_series,
            self.calc_batch,
            self.calc_sample,
        ]

    def tearDown(self):
        self.temp_dir.cleanup()
        
    def test_save_csv_single_sheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "output.csv")
            utils.save_results(path, self.mixed_data, filetype="csv", mode="single_sheet")
            self.assertTrue(os.path.isfile(path))
            df = pd.read_csv(path)
            self.assertGreaterEqual(len(df), len(self.mixed_data))

    def test_save_csv_multi_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = os.path.join(tmp, "out.csv")
            utils.save_results(base, self.mixed_data, filetype="csv", mode="multi_file")
            csv_files = [f for f in os.listdir(tmp) if f.endswith(".csv")]
            self.assertEqual(len(csv_files), len(self.mixed_data))

    def test_save_excel_single_sheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.xlsx")
            utils.save_results(path, self.mixed_data, filetype="excel", mode="single_sheet")
            self.assertTrue(os.path.isfile(path))
            df = pd.read_excel(path)
            self.assertGreaterEqual(len(df), len(self.mixed_data))

    def test_save_excel_multi_sheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.xlsx")
            utils.save_results(path, self.mixed_data, filetype="excel", mode="multi_sheet")
            self.assertTrue(os.path.isfile(path))
            xl = pd.ExcelFile(path)
            self.assertEqual(len(xl.sheet_names), len(self.mixed_data))

    def test_save_excel_multi_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = os.path.join(tmp, "out.xlsx")
            utils.save_results(base, self.mixed_data, filetype="excel", mode="multi_file")
            xlsx_files = [f for f in os.listdir(tmp) if f.endswith(".xlsx")]
            self.assertEqual(len(xlsx_files), len(self.mixed_data))

    def test_save_with_descriptions(self):
        descriptions = [f"Item {i}" for i in range(len(self.mixed_data))]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "desc.csv")
            utils.save_results(path, self.mixed_data, filetype="csv", mode="single_sheet", descriptions=descriptions)
            df = pd.read_csv(path)
            self.assertIn("Description", df.columns)
            self.assertGreaterEqual(df["Description"].nunique(), len(self.mixed_data) // 2)  # flexible match

    def test_description_mismatch_raises(self):
        bad_descriptions = ["only one"]
        with self.assertRaises(ValueError):
            utils.save_results("test.csv", self.mixed_data, filetype="csv", descriptions=bad_descriptions)

    def test_invalid_mode_for_csv_raises(self):
        with self.assertRaises(ValueError):
            utils.save_results("bad.csv", self.mixed_data, filetype="csv", mode="multi_sheet")

    def test_invalid_filetype_raises(self):
        with self.assertRaises(ValueError):
            utils.save_results("bad.txt", self.mixed_data, filetype="txt")

    def test_invalid_mode_excel_raises(self):
        with self.assertRaises(ValueError):
            utils.save_results("bad.xlsx", self.mixed_data, filetype="excel", mode="nonsense_mode")
