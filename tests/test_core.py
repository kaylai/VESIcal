import unittest
from VESIcal import core
from VESIcal import sample_class
import pandas as pd
    
class TestOxideCationConversions(unittest.TestCase):
    def setUp(self):
        self.oxide_data_dict_keys = [
            'mass',
            'cation',
            'cation_num',
            'oxygen_num',
            'cation_charge',
            'cation_mass',
            'groups', 
        ]
        
    def test_all_oxide_cation_mappings_complete(self):
        """Verify all oxides have cations and all mappings are consistent."""
                
        # Test 1: All oxides and cations in conversion dicts
        for oxide in core.oxides:
            with self.subTest(oxide=oxide):
                # Check oxide->cation mapping exists
                self.assertIn(oxide, core.oxides_to_cations,
                            f"{oxide} missing from oxides_to_cations")
        
        for cation in core.cations:
            with self.subTest(cation=cation):
                self.assertIn(cation, core.cations_to_oxides,
                              f"{cation} missing from cations_to_oxides")
        
        # Test 2: All oxides have all values defined in oxide_data dict
        for oxide in core.oxides:
            with self.subTest(oxide=oxide):
                for k in core.oxide_data[oxide].keys():
                    self.assertIn(k,
                              self.oxide_data_dict_keys,
                              f"{oxide} has incomplete data in oxide_data dict")
        
        # Test 3: All oxides in property dictionaries 
        # Check oxide->cation mapping exists
        for oxide in core.oxides:
            with self.subTest(oxide=oxide):
                self.assertIn(oxide, core.oxideMass,
                            f"{oxide} missing from oxideMass")
                
                self.assertIn(oxide, core.CationNum,
                            f"{oxide} missing from CationNum")
                
                self.assertIn(oxide, core.OxygenNum,
                            f"{oxide} missing from OxygenNum")
                
                self.assertIn(oxide, core.CationCharge,
                            f"{oxide} missing from CationCharge")
                
                self.assertIn(oxide, core.CationMass,
                            f"{oxide} missing from CationMass")