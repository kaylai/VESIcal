# ---------DEFINE SOME CONSTANTS---------- #
# Partial Molar Volumes
# Volumes for SiO2, Al2O3, MgO, CaO, Na2O, K2O at Tref=1773 K (Lange, 1997; CMP)
# Volume for H2O at Tref=1273 K (Ochs and Lange, 1999)
# Volume for FeO at Tref=1723 K (Guo et al., 2014)
# Volume for Fe2O3 at Tref=1723 K (Liu and Lange, 2006)
# Volume for TiO2 at Tref=1773 K (Lange and Carmichael, 1987)
MolarVolumes = {'SiO2': 26.86, 'TiO2': 28.32, 'Al2O3': 37.42, 'Fe2O3': 41.50,
                'FeO': 12.68, 'MgO': 12.02, 'CaO': 16.90, 'Na2O': 29.65,
                'K2O': 47.28, 'H2O': 22.9}

# Partial Molar Volume uncertainties
# value = 0 if not reported
MolarVolumeErrors = {'SiO2': 0.03, 'TiO2': 0.0, 'Al2O3': 0.09, 'Fe2O3': 0.0,
                     'FeO': 0.0, 'MgO': 0.07, 'CaO': 0.06, 'Na2O': 0.07,
                     'K2O': 0.10, 'H2O': 0.60}

# dV/dT values
# MgO, CaO, Na2O, K2O Table 4 (Lange, 1997)
# SiO2, TiO2, Al2O3 Table 9 (Lange and Carmichael, 1987)
# H2O from Ochs & Lange (1999)
# Fe2O3 from Liu & Lange (2006)
# FeO from Guo et al (2014)
dVdTs = {'SiO2': 0.0, 'TiO2': 0.00724, 'Al2O3': 0.00262, 'Fe2O3': 0.0,
         'FeO': 0.00369, 'MgO': 0.00327, 'CaO': 0.00374, 'Na2O': 0.00768,
         'K2O': 0.01208, 'H2O': 0.0095}

# dV/dT uncertainties
# value = 0 if not reported
dVdTErrors = {'SiO2': 0.0, 'TiO2': 0.0, 'Al2O3': 0.0, 'Fe2O3': 0.0, 'FeO': 0.0,
              'MgO': 0.0, 'CaO': 0.0, 'Na2O': 0.0, 'K2O': 0.0, 'H2O': 0.00080}

# dV/dP values
# Anhydrous component data from Kess and Carmichael (1991)
# H2O data from Ochs & Lange (1999)
dVdPs = {'SiO2': -0.000189, 'TiO2': -0.000231, 'Al2O3': -0.000226,
         'Fe2O3': -0.000253, 'FeO': -0.000045, 'MgO': 0.000027, 'CaO': 0.000034,
         'Na2O': -0.00024, 'K2O': -0.000675, 'H2O': -0.00032}

# dV/dP uncertainties
dVdPErros = {'SiO2': 0.000002, 'TiO2': 0.000006, 'Al2O3': 0.000009,
             'Fe2O3': 0.000009, 'FeO': 0.000003, 'MgO': 0.000007,
             'CaO': 0.000005, 'Na2O': 0.000005, 'K2O': 0.000014,
             'H2O': 0.000060}

# Tref values
Trefs = {'SiO2': 1773, 'TiO2': 1773, 'Al2O3': 1773, 'Fe2O3': 1723, 'FeO': 1723,
         'MgO': 1773, 'CaO': 1773, 'Na2O': 1773, 'K2O': 1773, 'H2O': 1273}

# species used for density calculation
densitySpecies = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'MgO', 'CaO', 'Na2O',
                  'K2O', 'H2O']

# species used for viscosity calculation
viscositySpecies = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'Fe2O3', 'MnO', 'MgO',
                    'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'F2O']
