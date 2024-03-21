from VESIcal import core

import numpy as np
import warnings as w

giordano_oxide_mass = {'SiO2': 60.0843,
                       'TiO2': 79.8658,
                       'Al2O3': 101.961276,
                       'FeO': 71.8444,
                       'MnO': 70.937449,
                       'MgO': 40.3044,
                       'CaO': 56.0774,
                       'Na2O': 61.97894,
                       'K2O': 94.1960,
                       'P2O5': 141.9446,
                       'H2O': 18.01528,
                       'F2O': 18.9984*2}


def _normalize_Giordano(sample):
    """
    Converts an input, not necessarily normalized composition to a
    composition in mol% oxides, using the normalization routine in the Giordano
    et al. (2008) excel spreadsheet.

    Parameters
    ----------
    sample  Sample class
        Magma major element composition.

    Returns
    -------
    dict
        Normalized major element composition in mol% oxides.
    """

    _sample = sample.get_composition(units='wtpt_oxides',
                                     oxide_masses=giordano_oxide_mass,
                                     asSampleClass=True)

    # convert FeO and Fe2O3 to FeOT as FeO
    with w.catch_warnings():
        w.filterwarnings("ignore", message="Species Fe2O3 not found in composition, " +
                                           "assigning value of 0.0")
        _sample_FeOT = (_sample.get_composition(units='wtpt_oxides',
                                                oxide_masses=giordano_oxide_mass,
                                                species='FeO') +
                        0.8998*_sample.get_composition(units='wtpt_oxides',
                                                       species='Fe2O3'))
    _sample.change_composition({'FeO': _sample_FeOT, 'Fe2O3': 0.0})

    # if passed, convert F to F2O
    _sample_F2O = _sample.get_composition(units='wtpt_oxides',
                                          oxide_masses=giordano_oxide_mass,
                                          species='F2O')
    if 'F' in _sample.get_composition().keys():
        _sample_F2O += (_sample.get_composition()['F'] *
                        giordano_oxide_mass['F2O']/core.oxideMass['F'])
    _sample.change_composition({'F2O': _sample_F2O})

    # ensure all necessary oxides are in sample
    oxides_giordano = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO',
                       'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'F2O']
    sample_wtper = {}
    for oxide in oxides_giordano:
        if oxide in _sample.get_composition().keys():
            sample_wtper[oxide] = _sample.get_composition(units='wtpt_oxides',
                                                          oxide_masses=giordano_oxide_mass,
                                                          species=oxide)
        else:
            sample_wtper[oxide] = 0.0

    # special normalize wt%
    orig_sample_sum_anhy = (sum([v for v in sample_wtper.values()]) -
                            sample_wtper['H2O'])

    sample_wtper_norm = {}
    for oxide in sample_wtper.keys():
        if oxide == 'H2O':
            sample_wtper_norm[oxide] = sample_wtper['H2O']
        else:
            sample_wtper_norm[oxide] = ((100.0-sample_wtper['H2O']) *
                                        sample_wtper[oxide]/(orig_sample_sum_anhy))

    GFW = (100.0 /
           (sum([sample_wtper_norm[oxide]/core.oxideMass[oxide] for oxide
                in oxides_giordano])))

    sample_moloxides_giordano = {oxide: GFW*sample_wtper_norm[oxide] /
                                 giordano_oxide_mass[oxide] for oxide in
                                 oxides_giordano}

    return sample_moloxides_giordano


def calculate_liquid_viscosity(sample, temperature, **kwargs):
    """
    Calculates the viscosity of a liquid given the chemical composition and
    temperature.

    Parameters
    ----------
    sample  Sample class
        Magma major element composition.

    temperature     float or int
        Temperature in degrees C.

    Returns
    -------
    float
        Log viscosity of the liquid in Pa*s, rounded to 4 dp.
    """

    # normalize using Giordano et al. (2008) spreadsheet routine
    comp_molpercent = _normalize_Giordano(sample)

    # known A constant value (Table 1, Giordano et al., 2008)
    A_value = -4.55

    def _calculate_partial_B_value(bulk_comp, b_number):
        """
        Calculate values for B (computed values in Giordano
        spreadsheet).

        B = sum(i=1-7)[b_i*M_i] + sum(j=1-3)[b_1j(M1_1j*M2_1j)]

        Parameters
        ----------
        bulk_comp   dict
            Bulk composition of the major oxides in mol percent.

        b_number    int
            Number of B value to calculate.

        Returns
        -------
        float
            Partial B value for given b_number. E.g., if given b_number=1, will
            return value of B1. B, then, is the sum of all partial B's.
        """
        comp = bulk_comp.copy()

        # Convert FeO, Fe2O3 to FeOT
        comp['FeOT'] = comp['FeO']

        # map b numbers to oxides
        bs_to_oxides = {1: ['SiO2', 'TiO2'],
                        2: ['Al2O3'],
                        3: ['FeOT', 'MnO', 'P2O5'],
                        4: ['MgO'],
                        5: ['CaO'],
                        6: ['Na2O', 'H2O', 'F2O']}

        # known constant values (Table 1, Giordano et al., 2008)
        b_number_to_constants = {1: 159.56,
                                 2: -173.34,
                                 3: 72.13,
                                 4: 75.69,
                                 5: -38.98,
                                 6: -84.08,
                                 7: 141.54,
                                 11: -2.43,
                                 12: -0.91,
                                 13: 17.62}

        if b_number <= 6:
            partial_B_value = (b_number_to_constants[b_number] *
                               sum([comp[bs_to_oxides[b_number][x]] for x in
                                    range(len(bs_to_oxides[b_number]))]))
        if b_number == 7:
            partial_B_value = (b_number_to_constants[b_number] *
                               (comp['H2O'] + comp['F2O'] +
                               np.log(1+comp['H2O'])))

        if b_number == 11:
            partial_B_value = (b_number_to_constants[b_number] *
                               (comp['SiO2']+comp['TiO2']) *
                               (comp['FeOT']+comp['MnO']+comp['MgO']))

        if b_number == 12:
            partial_B_value = (b_number_to_constants[b_number] *
                               (comp['SiO2']+comp['TiO2']+comp['Al2O3'] +
                               comp['P2O5']) *
                               (comp['Na2O']+comp['K2O']+comp['H2O']))

        if b_number == 13:
            partial_B_value = (b_number_to_constants[b_number] *
                               (comp['Al2O3']) *
                               (comp['Na2O']+comp['K2O']))

        return partial_B_value

    def _calculate_partial_C_value(bulk_comp, c_number):
        """
        Calculate values for C (computed values in Giordano
        spreadsheet).

        C = sum(i=1-6)[c_i*N_i] + [c_11*(N1_11*N2_11)]

        Parameters
        ----------
        bulk_comp   dict
            Bulk composition of the major oxides in mol percent.

        c_number    int
            Number of C value to calculate.

        Returns
        -------
        float
            Partial C value for given c_number. E.g., if given c_number=1, will
            return value of C1. C, then, is defined by the equation above.
        """
        comp = bulk_comp.copy()

        # Convert FeO, Fe2O3 to FeOT
        comp['FeOT'] = comp['FeO']

        # map b numbers to oxides
        cs_to_oxides = {1: ['SiO2'],
                        2: ['TiO2', 'Al2O3'],
                        3: ['FeOT', 'MnO', 'MgO'],
                        4: ['CaO'],
                        5: ['Na2O', 'K2O']}

        # known constant values (Table 1, Giordano et al., 2008)
        c_number_to_constants = {1: 2.75,
                                 2: 15.72,
                                 3: 8.32,
                                 4: 10.20,
                                 5: -12.29,
                                 6: -99.54,
                                 11: 0.30}

        if c_number <= 5:
            partial_C_value = (c_number_to_constants[c_number] *
                               sum([comp[cs_to_oxides[c_number][x]] for x in
                                    range(len(cs_to_oxides[c_number]))]))
        if c_number == 6:
            partial_C_value = (c_number_to_constants[c_number] *
                               (np.log(1+comp['H2O'] + comp['F2O'])))

        if c_number == 11:
            partial_C_value = (c_number_to_constants[c_number] *
                               (comp['Al2O3'] + comp['FeOT'] + comp['MnO'] +
                                comp['MgO'] + comp['CaO'] - comp['P2O5']) *
                               (comp['Na2O'] + comp['K2O'] + comp['H2O'] +
                                comp['F2O']))

        return partial_C_value

    b_numbers = [1, 2, 3, 4, 5, 6, 7, 11, 12, 13]
    B_value = sum([_calculate_partial_B_value(comp_molpercent, b) for b in b_numbers])

    c_numbers = [1, 2, 3, 4, 5, 6, 11]
    C_value = sum([_calculate_partial_C_value(comp_molpercent, c) for c in c_numbers])

    log_nu_Pas = A_value + B_value/((temperature+273.15)-C_value)

    return round(log_nu_Pas, 4)
