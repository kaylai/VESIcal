from VESIcal import core
from VESIcal.thermo import thermo_core as tc


def calculate_liquid_density(sample, pressure, temperature, **kwargs):
    """ Calculates the density of the given liquid using the DensityX model of
    Iacovino and Till (2019).

    Parameters
    ----------
    sample  Sample class
        Magma major element composition.
    pressure    float
        Total pressure in bars.
    temperature     float
        Temperature in degrees C

    Returns
    -------
    float
        The density of the liquid in g/L, rounded to 3 dp.
    """

    comp_molfrac = sample.get_composition(units='mol_oxides')
    # get only species considered in density calculation
    comp_molfrac = {k: v for k, v in comp_molfrac.items() if
                    k in tc.densitySpecies}
    tempK = temperature + 273.15

    indiv_Vliqs = []
    for k, v in comp_molfrac.items():
        indiv_Vliq = ((tc.MolarVolumes[k] + (tc.dVdTs[k] * (tempK-tc.Trefs[k])) +
                      (tc.dVdPs[k] * (pressure-1))) * v)
        indiv_Vliqs.append(indiv_Vliq)

    Vliq_sum = sum(indiv_Vliqs)
    Sum_X_MW = sum([v*core.oxideMass[k] for k, v in comp_molfrac.items()])

    density_g_per_cm3 = Sum_X_MW/Vliq_sum
    density_g_per_L = density_g_per_cm3*1000

    return round(density_g_per_L, 3)
