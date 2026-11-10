"""Independent stoichiometric fixtures for the public molar-input paths."""

from math import fsum

import numpy as np
import pandas as pd


# Project molecular masses, with cation counts stated independently of Sample.
INPUT_SPECIES = {
    'mol_oxides': [('SiO2', 'SiO2', 60.083, 1),
                   ('MgO', 'MgO', 40.304, 1)],
    'mol_cations': [('Si', 'SiO2', 60.083, 1),
                    ('Na', 'Na2O', 61.979, 2),
                    ('Al', 'Al2O3', 101.961, 2)],
}


def molar_frame(units, rows=2, scale=1.0):
    if units == 'mol_oxides':
        values = {'SiO2': [.5, .25], 'MgO': [.5, .75]}
    else:
        values = {'Si': [.25, .5], 'Na': [.25, .25], 'Al': [.5, .25]}
    data = {'Label': ['A', 'B'][:rows],
            'Notes': ['first', 'second'][:rows],
            'Replicate': [2, 5][:rows]}
    data.update({name: [v * scale for v in values[:rows]]
                 for name, values in values.items()})
    return pd.DataFrame(data, index=[19, 3][:rows])


def expected_wtpercent(frame, units):
    """Use Python arithmetic, never Sample or a BatchFile converter."""
    converted = []
    for _, row in frame.iterrows():
        masses = {oxide: float(row[name]) * mass / count
                  for name, oxide, mass, count in INPUT_SPECIES[units]}
        total = fsum(masses.values())
        converted.append({oxide: 100 * value / total
                          for oxide, value in masses.items()})
    return pd.DataFrame(converted, index=frame['Label'].astype(str).tolist())


def assert_molar_batch(case, batch, frame, units):
    expected = expected_wtpercent(frame, units)
    composition = batch.get_composition()
    case.assertEqual(composition.index.tolist(), expected.index.tolist())
    np.testing.assert_allclose(composition[expected.columns], expected,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(composition.sum(axis=1), 100,
                               rtol=0, atol=1e-9)
    absent = composition.drop(columns=expected.columns)
    np.testing.assert_array_equal(absent.to_numpy(), 0)

    data = batch.get_data()
    case.assertTrue(data.columns.is_unique)
    for name in ['Notes', 'Replicate']:
        case.assertEqual(data[name].tolist(), frame[name].tolist())

    molar_data = batch.get_data(units=units)
    case.assertTrue(molar_data.columns.is_unique)
    fields = [name for name, _, _, _ in INPUT_SPECIES[units]]
    amounts = frame[fields].to_numpy(dtype=float)
    expected_moles = amounts / amounts.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(molar_data[fields], expected_moles,
                               rtol=0, atol=1e-9)
    for name in ['Notes', 'Replicate']:
        case.assertEqual(molar_data[name].tolist(), frame[name].tolist())
    if units == 'mol_cations':
        for name in fields:
            case.assertNotIn(name, batch.data.columns)
