# -*- coding: utf-8 -*-
# Script written by Simon Matthews (simonmatthews@jhu.edu).
# All the code contained in these scripts will ultimately be published
# alongside an updated interface to MagmaSat by Kayla Iacovino and colleagues.
# Before publication some code needs to be tidied up, but this current
# implementation should be fine for calculating saturation pressures in basalts
# with low H2O concentrations. One particular caveat is that the EOS implemented
# for the Dixon (1997) solubility model is not the same as the one they used
# when calibrating. However, they encourage use of a different EOS at high
# pressures anyway, and there should be little practical difference in the results
# using the one implemented here (Tim Holland's CHO model).

# VERSION 2.2- SEPTEMBER 2019

import numpy as np
#import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.interpolate import RectBivariateSpline
#from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
#sns.set_style('ticks')
params = {'text.usetex': False, 'mathtext.fontset': 'stixsans',
          'xtick.labelsize':10, 'ytick.labelsize':10, 'axes.labelsize':10,
          'xtick.direction':'in','ytick.direction':'in','figure.figsize':(6,4),
          'xtick.top':True,'ytick.right':True}
plt.rcParams.update(params)


# TO DO:
# Harmonise temperature units between different parameterisation inputs
# Use the correct CO2 and H2O EOS with Dixon
# Update the H2O wrapper Function
# Allow parameterisation specific variable propagation through the wrapper functions

#P0 = 1 # bar
#T0 = 1473.15 # K
#R = 83.15 # cm3 bar mol-1 K-1


oxideMass = {'SiO2':28.085+32,'MgO':24.305+16,'FeO':55.845+16,'CaO':40.078+16,'Al2O3':2*26.982+16*3,'Na2O':22.99*2+16,
             'K2O':39.098*2+16,'MnO':54.938+16,'TiO2':47.867+32,'P2O5':2*30.974+5*16,'Cr2O3':51.996*2+3*16,
             'NiO':58.693+16,'SO2':32.06+32,'F':18.998,'Cl':35.45,'Fe2O3':55.845*2+16*3,'H2O':18}
CationNum = {'SiO2':1,'MgO':1,'FeO':1,'CaO':1,'Al2O3':2,'Na2O':2,
             'K2O':2,'MnO':1,'TiO2':1,'P2O5':2,'Cr2O3':2,
             'NiO':1,'SO2':1,'F':1,'Cl':1,'Fe2O3':2,'H2O':2}
CationCharge = {'SiO2':4,'MgO':2,'FeO':2,'CaO':2,'Al2O3':3,'Na2O':1,
             'K2O':1,'MnO':2,'TiO2':4,'P2O5':5,'Cr2O3':3,
             'NiO':2,'SO2':2,'F':-1,'Cl':-1,'Fe2O3':3,'H2O':1}
CationMass = {'SiO2':28.085,'MgO':24.305,'FeO':55.845,'CaO':40.078,'Al2O3':26.982,'Na2O':22.990,
             'K2O':39.098,'MnO':54.938,'TiO2':47.867,'P2O5':30.974,'Cr2O3':51.996,
             'NiO':58.693,'SO2':32.06,'F':18.998,'Cl':35.45,'Fe2O3':55.845,'H2O':2}

# MORB
Shishkina_ME_MORB = pd.Series({'SiO2':50.69,
                             'TiO2':1.46,
                             'Al2O3':16.95,
                             'FeO':8.51,
                             'MnO':0.13,
                             'MgO':7.48,
                             'CaO':12.03,
                             'Na2O':2.52,
                             'K2O':0.22,
                             'P2O5':0.0,})

# Alkali Basalt
Shishkina_ME_AB2518 = pd.Series({'SiO2':46.21,
                             'TiO2':2.70,
                             'Al2O3':14.65,
                             'FeO':11.74,
                             'MnO':0.16,
                             'MgO':8.71,
                             'CaO':10.63,
                             'Na2O':3.51,
                             'K2O':1.06,
                             'P2O5':0.56,
                             'Cr2O3':0.06})

# Island Arc Tholeiite
Shishkina_ME_ArcTholeiite = pd.Series({'SiO2':50.17,
                             'TiO2':0.92,
                             'Al2O3':18.28,
                             'FeO':9.37,
                             'MnO':0.17,
                             'MgO':7.00,
                             'CaO':11.37,
                             'Na2O':2.33,
                             'K2O':0.23,
                             'P2O5':0.15})

# Nephelinite
Shishkina_ME_Nephelinite = pd.Series({'SiO2':42.32,
                             'TiO2':2.26,
                             'Al2O3':11.80,
                             'FeO':11.00,
                             'MnO':0.19,
                             'MgO':13.31,
                             'CaO':13.23,
                             'Na2O':3.72,
                             'K2O':0.96,
                             'P2O5':1.13,
                             'Cr2O3':0.08})

# Nephelinite
Shishkina_ME_Basanite = pd.Series({'SiO2':43.64,
                             'TiO2':2.64,
                             'Al2O3':12.65,
                             'FeO':11.54,
                             'MnO':0.19,
                             'MgO':12.07,
                             'CaO':11.82,
                             'Na2O':3.68,
                             'K2O':1.01,
                             'P2O5':0.7,
                             'Cr2O3':0.06})

# Etna
IM_Etna = pd.Series({'SiO2':47.95,
                     'TiO2':1.67,
                     'Al2O3':17.32,
                     'FeO':10.24,
                     'MgO':5.76,
                     'CaO':10.93,
                     'Na2O':3.45,
                     'K2O':1.99,
                     'P2O5':0.51})


RTlnf_CO2_table = pd.read_csv('CO2_RTlnf_PitzerSterner.csv',index_col=0)
RTlnf_H2O_table = pd.read_csv('H2O_RTlnf_PitzerSterner.csv',index_col=0)
RTlnf_CO2 = RectBivariateSpline(RTlnf_CO2_table.index,RTlnf_CO2_table.columns.astype(np.float),RTlnf_CO2_table)
RTlnf_H2O = RectBivariateSpline(RTlnf_H2O_table.index,RTlnf_H2O_table.columns.astype(np.float),RTlnf_H2O_table)


# Dixon (1997) Eq (1)
def Dixon_XCO3_pure(P,XCO3Std=3.8e-7,fCO2=0.0,MajorElements='default'):
    """Dixon (1997) Eqn (1).

    Calculates the mole fraction of carbonate ion dissolved in a tholeiitic basalt in
    equilibrium with pure CO2 fluid at 1200C, at a chosen pressure.

    Parameters
    ----------
    P:  float
        Pressure in bar
    XCO3Std:    float
        Mole fraction of carbonate dissolved in tholeiitic basalt in equilibrium with
        pure CO2 fluid at 1200C and 1 bar. Default value is 3.8e-7, as suggested in paper.
        If MajorElements variable not set to string, this will be replaced.
    fCO2:   float
        Fugacity of CO2 at the pressure of interest. If value set to 0.0, fCO2 will be
        set to that calculated using the Pitzer & Sterner equations by the CHO software.
        This is the default option.
    MajorElements:  series
        If set to string XH2OStd variable passed will be used (or default value).
        Series of major element oxide concentrations in wt%, requiring SiO2 as minimum.
        Used to calculated XCO3std from compositional parametrisation."""

    DeltaVr = 23 #cm3 mole-1
    P0 = 1
    R = 83.15
    T = 1473
    T0 = 1473.15

    if isinstance(MajorElements,str) == False:
        XCO3Std = Dixon_XCO3_Std(MajorElements)
    if fCO2 == 0.0:
        fCO2=PitzerSterner_fCO2(P,T)

    return XCO3Std * fCO2 * np.exp(-DeltaVr * (P-P0)/(R*T0))

# Dixon (1997) Eq (2)
def Dixon_XH2O_pure(P,XH2OStd=3.28e-5,fH2O=0.0,MajorElements='default'):
    """Dixon (1997) Eqn (2).

    Calculates the mole fraction of molecular H2O dissolved in a tholeiitic basalt in
    equilibrium with pure H2O fluid at 1200C, at a chosen pressure.

    Parameters
    ----------
    P: float
        Pressure in bar
    XH2OStd: float
        Mole fraction of cmolecular H2O in tholeiitic basalt in equilibrium with
        pure H2O fluid at 1200C and 1 bar. Default value is 3.28e-5, as suggested in paper.
        If MajorElements variable not set to string, this will be replaced.
    fH2O: float
        Fugacity of H2O at the pressure of interest. If value set to 0.0, fH2O will be
        set to that calculated using the Pitzer & Sterner equations by the CHO software.
        This is the default option.
    MajorElements:  series
        If set to string XH2OStd variable passed will be used (or default value).
        Series of major element oxide concentrations in wt%, requiring SiO2 as minimum.
        Used to calculated XH2Ostd from compositional parametrisation."""
    VH2O = 12 #cm3 mole-1
    P0 = 1
    R = 83.15
    T = 1473
    T0 = 1473.15

    if isinstance(MajorElements,str) == False:
        XH2OStd = Dixon_XH2O_Std(MajorElements)
    if fH2O == 0.0:
        fH2O=PitzerSterner_fH2O(P,T)

    return XH2OStd * fH2O * np.exp(-VH2O * (P-P0)/(R*T0))

# Dixon (1997) Eq (3)
def Dixon_CO2_ppm(XCO3):
    """Dixon (1997) Eqn (3).

    Calculates the concentration (in ppm) of CO2 in the melt given the mole
    fraction of carbonate ion.

    Parameters
    ----------
    XCO3: float
        Mole fraction of carbonate ion in the melt."""
    return 1e4 * (4400 * XCO3) / (36.6 - 44*XCO3)

# Dixon (1997) Eq (4)
def Dixon_XOH_root(XOH,XH2O):
    """Dixon (1997) Eq (4).

    Returns the difference between the two sides of equation (4) which relates
    the mole fraction of hydroxyl groups to the mole fraction of molecular water.
    For use in finding the root of the equation.

    It is a regular solution model (Silver and Stolper, 1989).

    Parameters
    ----------
    XOH: float
        Mole fraction of hydroxyl groups dissolved in melt.
    XH2O: float
        Mole fraction of molecular water dissolved in melt."""
    A = 0.403
    B = 15.333
    C = 10.894

    lhs = - np.log(XOH**2.0/(XH2O*(1.0-XH2O)))
    rhs = A + B*XOH + C*XH2O

    return rhs - lhs

# Find the root of Dixon (1997) Eq (4)
def Dixon_XOH_solve(XH2O):
    """Solves Dixon (1997) Eq (4).

    Finds the mole fraction of hydroxyl groups dissolved in the melt, given the
    mole fraction of molecular water. Uses the Newton method to find the root of
    Eq (4).

    Often struggles to find root when XH2O -> 0.

    Parameters
    ----------
    XH2O: float
        Mole fraction of molecular water dissolved in melt."""
    results = list()
    if isinstance(XH2O,np.ndarray)==False and isinstance(XH2O,list) == False:
        XH2O=[XH2O]
    for i in range(np.shape(XH2O)[0]):
        if XH2O[i] == 0:
            results.append(0)
        else:
            results.append(newton(Dixon_XOH_root,0.01,args=(XH2O[i],)))
    return np.array(results)

# Dixon (1997) Eq (5) and Eq(6)
def Dixon_H2O_wt(XH2O,XOH):
    """Dixon (1997) Eq (5) and Eq (6).

    Calculates total H2O (wt%) dissolved in the melt, given the molar fractions
    of molecular water and hydroxyl groups dissolved.

    Parameters
    ----------
    XH2O: float
        Mole fraction of molecular water dissolved in the melt.
    XOH: float
        Mole fraction of hydroxyl groups dissolved in the melt."""
    XB = XH2O + 0.5*XOH
    return 1801.5*XB/(36.6-18.6*XB)

# Dixon (1997) Eq (8)
def Dixon_XCO3_Std(MajorElements='default'):
    """Dixon (1997) Eq (8).

    The compositional parameterisation for the mole fraction of carbonate ions
    dissolved in a melt in equilibrium with CO2 vapour at 1200C and 1 bar.

    Parameters
    ----------
    MajorElements:    dict or series (or str)
        A dictionary or series containing 'SiO2' as a label for the SiO2 content
        of the melt. If the variable is a str, it will return the default XCO3_std
        value."""
    if isinstance(MajorElements,str):
        return 3.8e-7
    else:
        return 8.7e-6 - 1.7e-7*MajorElements['SiO2']

# Dixon (1997) Eq (8)
def Dixon_XH2O_Std(MajorElements='default'):
    """Dixon (1997) Eq (9).

    The compositional parameterisation for the mole fraction of molecular water
    dissolved in a melt in equilibrium with H2O vapour at 1200C and 1 bar.

    Parameters
    ----------
    MajorElements:    dict or series (or str)
        A dictionary or series containing 'SiO2' as a label for the SiO2 content
        of the melt. If the variable is a str, it will return the default XH2O_std
        value."""
    if isinstance(MajorElements,str):
        return 3.28e-5
    else:
        return -3.04e-5 + 1.29e-6*MajorElements['SiO2']


def PitzerSterner_fCO2(P_list,T):
    """Returns the fugacity of pure CO2 at selected P and T using a lookup table of
    values calculated from the Pitzer and Sterner equations, calculated using the
    CHO software by Tim Holland. Requires the lookup table to be saved in the same
    folder.

    Parameters
    ----------
    P_list:     float or list
        Pressure (bar), or list of pressures to lookup.
    T:  float
        Temperature (K) to lookup"""
    # Output table in units of kbar and degC
    if isinstance(P_list,list)==False and isinstance(P_list,np.ndarray)==False:
        P_list = [P_list]
    f = list()
#    RTlnf_CO2 = RectBivariateSpline(RTlnf_CO2_table.index,RTlnf_CO2_table.columns.astype(np.float),RTlnf_CO2_table)

    for P in P_list:
        if P <=10:
            f.append(P)
        else:
            P=P/1000

            out = RTlnf_CO2(P,T-273.15)[0,0]

            f.append(np.exp(out/(8.314e-3*T)))
    if np.shape(f)[0] == 1:
        f = f[0]
    return np.array(f)

def PitzerSterner_fH2O(P_list,T):
    """Returns the fugacity of pure H2O at selected P and T using a lookup table of
    values calculated from the Pitzer and Sterner equations, calculated using the
    CHO software by Tim Holland. Requires the lookup table to be saved in the same
    folder.

    Parameters
    ----------
    P_list:     float or list
        Pressure (bar), or list of pressures to lookup.
    T:  float
        Temperature (C) to lookup"""
    # Output table in units of kbar and degC
    if isinstance(P_list,list)==False and isinstance(P_list,np.ndarray)==False:
        P_list = [P_list]
    f = list()
#    RTlnf_H2O = RectBivariateSpline(RTlnf_H2O_table.index,RTlnf_H2O_table.columns.astype(np.float),RTlnf_H2O_table)

    for P in P_list:
        if P <= 10:
            f.append(P)
        else:
            P=P/1000

            out = RTlnf_H2O(P,T-273.15)[0,0]

            f.append(np.exp(out/(8.314e-3*T)))
    if np.shape(f)[0] == 1:
        f = f[0]
    return np.array(f)

# Dixon (1997) Eq (11)
def XCO3_mixed(XH2O,XH2O_pure,XCO3_pure):
    """Dixon (1997) Eq (11).

    Returns the molar proportion of carbonate ions dissolved in a melt in the presence
    of a mixed H2O and CO2 fluid, assuming ideal mixing.

    Parameters
    ----------
    XH2O:   float
        The molar proportion of molecular water dissolved in the melt.
    XH2O_pure:  float
        The molar proportion of molecular water dissolved in the melt in the
        presence of a pure H2O vapour.
    XCO3_pure:  float
        The molar proportion of carbonate ion dissolved in the melt in the presence
        or pure CO2 vapour."""
    return XCO3_pure * (1-XH2O/XH2O_pure)

# Dixon (1997) Eq (11)
def XH2O_mixed(XCO3,XH2O_pure,XCO3_pure):
    """Dixon (1997) Eq (11).

    Returns the molar proportion of molecular water dissolved in a melt in the presence
    of a mixed H2O and CO2 fluid, assuming ideal mixing.

    Parameters
    ----------
    XCO3:   float
        The molar proportion of carbonate ion dissolved in the melt.
    XH2O_pure:  float
        The molar proportion of molecular water dissolved in the melt in the
        presence of a pure H2O vapour.
    XCO3_pure:  float
        The molar proportion of carbonate ion dissolved in the melt in the presence
        or pure CO2 vapour."""
    return XH2O_pure * (1-XCO3/XCO3_pure)

def Dixon_CO2(P,MajorElements='default'):
    """Calculates CO2 solubility in equilibrium with a pure CO2 vapour,
    according to the Dixon parameterisation.

    Parameters
    ----------
    P:  float
        Pressure
    MajorElements: series
        Major element oxides in wt%, only SiO2 required. If set to string, it will
        use the non-compositional parameterisation."""
    XCO3= Dixon_XCO3_pure(P,MajorElements=MajorElements)
    return Dixon_CO2_ppm(XCO3)

def Dixon_H2O(P,MajorElements='default'):
    """Calculates H2O solubility in equilibrium with a pure H2O vapour,
    according to the Dixon parameterisation.

    Parameters
    ----------
    P:  float
        Pressure
    MajorElements: series
        Major element oxides in wt%, only SiO2 required. If set to string, it will
        use the non-compositional parameterisation."""
    XH2O = Dixon_XH2O_pure(P,MajorElements=MajorElements)
    XOH = Dixon_XOH_solve(XH2O)
    return Dixon_H2O_wt(XH2O,XOH)

def Dixon_Degassing_beta(P,MajorElements='default',PiStar=0.0,method='Dixon'):
    """Dixon (1995) Eq (A10).

    Returns the equilibrium fractionation factor (beta) for CO2-H2O between melt
    and vapour.

    Parameters
    ----------
    P:  float
        Pressure at which fractionation is allowed to happen.
    MajorElements:  series
        Major element oxide wt%. If set to a string and method is set to Dixon,
        no compositional dependence will be applied. If PiStar is non-zero and
        method is Shishkina, it will be ignored.
    PiStar:     float
        Pi* compositional parameter, used by the Shishkina parameterisation. If
        set to zero, it will be calculated from MajorElements.
    method:     string
        Parameterisation (Dixon or Shishkina) to use in calculating CO2 and H2O
        solubility."""
    fH2O = PitzerSterner_fH2O(P,1473)
    fCO2 = PitzerSterner_fCO2(P,1473)

    if method == 'Dixon':
        XCO3 = Dixon_XCO3_pure(P,MajorElements=MajorElements)
        XH2O = Dixon_XH2O_pure(P,MajorElements=MajorElements)
    elif method == 'Shishkina':
        if PiStar == 0.0:
            PiStar = Shishkina_PiStar(MajorElements)
        XCO3 = Dixon_XCO3(Shishkina_CO2(P,PiStar=PiStar,MajorElements=MajorElements))
        XH2O = Dixon_XH2O(Shishkina_H2O(P,MajorElements=MajorElements))

    return (fCO2/fH2O)/(XCO3/XH2O)

def Dixon_Degassing_dOHdy(y):
    """Dixon (1995) Eq (A22).

    Derivative of mole fraction of hydroxyl ion with respect to the parameter y.

    Paramters
    ---------
    y:  float
        y parameter, as defined by equation A17. The fraction of total water
        present as OH groups."""
    return 0.618104 + 5.45815*y + 17.76879*y**2 - 30.71796*y**3 + 26.84815*y**4 - 9.43803*y**5

def Dixon_Degassing_y(OH_m,H2O_m):
    """Dixon (1995) Eq (A17).

    The fraction of total water present as OH groups.

    Parameters
    ----------
    OH_m:   float
        Mole fraction of OH ions in the melt.
    H2O_m:  float
        Mole fraction of molecular water in the melt."""
    return OH_m/(OH_m+2*H2O_m)

def Dixon_Degassing_Derivatives(P,beta,lamda,CO3_m,CO3_0,H2Otot_0,OH_m,H2O_m):
    """Dixon (1995) Eq (A26).

    Solves the linear set of equations corresponding to the conservation laws
    for degassing. Returns the values of the three derivatives which satisy the
    equations, as an array in the order dOH_m/df, dH2O_m/df, dCO3_m/df.

    Parameters
    ----------
    P:  float
        Pressure at which degassing is occuring.
    beta:   float
        Fraction factor of CO2-H2O between vapour and melt.
    lamda:  float
        Fraction of vapour remaining in the system. Closed system is 1.0, Open
        system is 0.0.
    CO3_m:  float
        Mole fraction carbonate in the melt.
    CO3_0:  float
        Mole fraction carbonate in the melt, prior to degassing.
    H2Otot_0: float
        Total mole fraction of water dissolved in the melt, prior to degassing.
    OH_m:   float
        Mole fraction of hydroxyl ion in the melt.
    H2O_m:  float
        Mole fraction of molecular water dissolved in the melt.
    """
    y = Dixon_Degassing_y(OH_m,H2O_m)
    dOHdy = Dixon_Degassing_dOHdy(y)

    a = np.zeros([3,3])

    a[0,0] = 0.5
    a[0,1] = 1
    a[0,2] = 1
    a[1,0] = -0.5*beta*CO3_m
    a[1,1] = (lamda-beta)*CO3_m - lamda*CO3_0
    a[1,2] = beta*lamda*(H2Otot_0-0.5*OH_m-H2O_m) + H2O_m
    a[2,0] = (y-1)*dOHdy + OH_m + 2*H2O_m
    a[2,1] = 2*y*dOHdy
    a[2,2] = 0

    b = np.zeros([3])

    b[0] = H2Otot_0 + CO3_0
    b[1] = 0
    b[2] = 0

    return np.linalg.solve(a,b)

def Dixon_Degassing(CO2_ppm_0,H2O_wt_0,P,lamda,MajorElements='Default',PiStar=0.0,method='Dixon',df=0.01):
    """
    Calculates a degassing trajectory according to the equations set out in the
    Dixon (1995) appendix. Returns a dataframe of melt compositions and related
    parameters.

    Parameters
    ----------
    CO2_ppm_0:   float
          Concentration of CO2 dissolved in the original melt, in ppm.
    H2O_wt_0:    float
       Concentration of H2O dissolved in the original melt, in wt%.
    P:   float
       Pressure at which degassing is taking place.
    lamda:  float
        Fraction of vapour retained. Closed system is 1.0, open system is 0.0.
    MajorElements:  series
        Major element oxide wt%. For the Dixon parameterisation, a string will yield
        the non-compositional parameterisation.
    PiStar:     Float
        Pi* compositional parameter used in the Shishkina method. If non-zero it will
        be used in preference to that calculated from MajorElements.
    method:     string
        Parameterisation (Dixon or Shishkina) to use for calculting CO2 and H2O
        solubilities.
    df:     float
        Discretization of fraction volatiles remaining in the melt.

   """
    # Calculate initial conditions:
    CO3_0 = Dixon_XCO3(CO2_ppm_0)
    H2Otot_0 = Dixon_XB(H2O_wt_0)
    H2O_0 = Dixon_XH2O(H2O_wt_0)
    OH_0 = (H2Otot_0 - H2O_0)*2
    SatP_0 = SaturationP(CO2_ppm_0,H2O_wt_0,MajorElements,PiStar,method)

    # Initiate Variables for loop
    CO2 = list([CO2_ppm_0])
    H2O = list([H2O_wt_0])
    CO3_m = list([CO3_0])
    H2O_m = list([H2O_0])
    OH_m = list([OH_0])
    f = list([1.0])
    SatP = list([SatP_0])
    beta = Dixon_Degassing_beta(P,MajorElements,PiStar,method)

    while f[-1]>0.0 and SatP[-1]>P:
        derivatives = Dixon_Degassing_Derivatives(P,beta,lamda,CO3_m[-1],CO3_0,H2Otot_0,OH_m[-1],H2O_m[-1])
        OH_m.append(OH_m[-1]-derivatives[0]*df)
        H2O_m.append(H2O_m[-1]-derivatives[1]*df)
        if CO3_m[-1] < 0:
            CO3_m.append(0.0)
        else:
            CO3_m.append(CO3_m[-1]-derivatives[2]*df)
        f.append(f[-1]-df)
        CO2.append(Dixon_CO2_ppm(CO3_m[-1]))
        H2O.append(Dixon_H2O_wt(H2O_m[-1],OH_m[-1]))
        SatP.append(SaturationP(CO2[-1],H2O[-1],MajorElements,PiStar,method,SatP[-1]))

    results = pd.DataFrame(
        {'f': f,
         'Saturation P': SatP,
         'CO2': CO2,
         'H2O': H2O,
         'CO3_m': CO3_m,
         'H2O_m': H2O_m,
         'OH_m' : OH_m
        })

    return results


def Shishkina_PiStar(MajorElements):
    """Shishkina et al. (2014) Eq (11)

    Calculates the Pi* parameter for use in calculating CO2 solubility.

    Parameters
    ----------
    MajorElements:  Series
        Series containing major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions.
    """
    MolarProps = MajorElements.copy()
    for ox in MolarProps.index:
        MolarProps[ox] = CationNum[ox]*MolarProps[ox]/oxideMass[ox]

    # Normalise
    MolarProps = MolarProps/MolarProps.sum()

    pi = ((MolarProps['CaO'] + 0.8*MolarProps['K2O'] + 0.7*MolarProps['Na2O'] +
        0.4*MolarProps['MgO'] + 0.4*MolarProps['FeO'])/(MolarProps['SiO2']+MolarProps['Al2O3']))
    return pi

def Shishkina_CO2(P,PiStar=0.0,MajorElements='default'):
    """Shishkina et al. (2014) Eq (13).

    Calculates the CO2 concentration (ppm) dissolved in a melt in equilibrium
    with a pure CO2 vapour, at a given pressure, for the given melt composition,
    or pi* value.

    Parameters
    ----------
    P: float
        Pressure in bar
    PiStar:  float
        PiStar compositional parameter. Ignored if set to 0.0.
    MajorElements:  series
        Series labelled with major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions. Ignored unless PiStar is set to
        0.0."""

    if PiStar == 0.0:
        PiStar = Shishkina_PiStar(MajorElements)
    P=P/10
    A = 1.150
    B = 6.71
    C= -1.345

    return np.exp(A*np.log(P)+B*PiStar+C)


def molar_props(MajorElements,oxides=True,count_H2O=True,include_O=False):
    """
    Returns the molar proportions for given wt% oxides.

    Parameters
    ----------
    MajorElements: series
        Series labelled with major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions. Ignored unless PiStar is set to
        0.0.
    oxides:     bool
        If true, returns mole fractions of the oxides, rather than the cations.
    count_H2O:  bool
        If true, counts the mole fraction of H2O.
    include_O (bool): Include O(2-) in the sum, calculated by charge balance. Only if
                      oxides = False.
    """
    MolarProps = MajorElements.copy()

    if count_H2O == False:
        MolarProps['H2O'] = 0.0


    if oxides == True:
        for ox in MolarProps.index:
            MolarProps[ox] = MolarProps[ox]/oxideMass[ox]
    else:
        if include_O == True:
            O_mols = 0.0
            for ox in MolarProps.index:
                MolarProps[ox] = CationNum[ox]*MolarProps[ox]/oxideMass[ox]
                O_mols = O_mols + CationCharge[ox]*MolarProps[ox]/2

        else:
            for ox in MolarProps.index:
                MolarProps[ox] = CationNum[ox]*MolarProps[ox]/oxideMass[ox]

    # Normalise
    if include_O == True:
        cation_sum = MolarProps.sum()
        MolarProps = MolarProps/cation_sum
        O_mols = O_mols*(MolarProps.sum()/cation_sum)
        return MolarProps, O_mols
    else:
        MolarProps = MolarProps/MolarProps.sum()
        return MolarProps

def IM_NBO_O(molarProps,hydrous=False):
    X = molarProps
    NBO = 2*(X['K2O']+X['Na2O']+X['CaO']+X['MgO']+X['FeO']-X['Al2O3'])
    O = 2*X['SiO2']+2*X['TiO2']+3*X['Al2O3']+X['MgO']+X['FeO']+X['CaO']+X['Na2O']+X['K2O']

    if hydrous == True:
        NBO = NBO + 2*X['H2O']
        O = O + X['H2O']

    return NBO/O


def IM_CO2(P,T=273+1200,MajorElements='default',P_CO2='pure',H2O=0.0,hydrous=False):
    """
    Iacono-Marziano (2012) Eq (7). Calculates the concentration of CO2 in
    ppm dissolved in the melt at vapour saturation, for specified pressure,
    temperature and partial pressure of CO2. Modifies the output of the equation
    to give CO2 rather than CO3.

    Parameters
    ----------
    P: float
        Pressure in bar
    T: float
        Temperature in K
    MajorElements: series
        Series labelled with major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions.
    P_CO2: float
        Partial pressure of CO2 in bar. If set to 'pure', function will set it as P.
    H2O: float
        Concentration of H2O in the melt in wt%.
    hydrous: bool
        Use the hydrous or anhydrous calibration
    """

    if hydrous == True:
        d = np.array([-16.4,4.4,-17.1,22.8])
        a = 1.0
        b = 17.3
        B = -6.0
        C = 0.12
    else:
        d = np.array([2.3,3.8,-16.3,20.1])
        a = 1.0
        b = 15.8
        B = -5.3
        C = 0.14

    if type(P_CO2) == str:
        P_CO2 = P

    if type(MajorElements) == str:
        MajorElements = Shishkina_ME_MORB

    MajorElements['H2O'] = H2O

    molarProps = molar_props(MajorElements)

    x = list()
    x.append(molarProps['H2O'])
    x.append(molarProps['Al2O3']/(molarProps['CaO']+molarProps['K2O']+molarProps['Na2O']))
    x.append(molarProps['FeO']+molarProps['MgO'])
    x.append(molarProps['Na2O']+molarProps['K2O'])
    x = np.array(x)


    CO3 = np.exp(np.sum(x*d) + a*np.log(P_CO2) +
                 b*IM_NBO_O(molarProps,hydrous)+
                 B+C*P/T)

    CO2 = CO3/(12+16*3)*(12+16*2)


    return CO2


def IM_H2O(P,T=273+1200,MajorElements='default',P_H2O='pure',hydrous=False):
    """
    Iacono-Marziano (2012) Eq (8). Calculates the concentration of CO2 in
    ppm dissolved in the melt at vapour saturation, for specified pressure,
    temperature and partial pressure of CO2. Modifies the output of the equation
    to give CO2 rather than CO3.

    Parameters
    ----------
    P: float
        Pressure in bar
    T: float
        Temperature in K
    MajorElements: series
        Series labelled with major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions.
    P_CO2: float
        Partial pressure of CO2 in bar
    hydrous: bool
        Use the hydrous or anhydrous calibration
    """

    if hydrous == True:
        a = 0.53
        b = 2.35
        B = -3.37
        C = -0.02
    else:
        a = 0.54
        b = 1.24
        B = -2.95
        C = 0.02

    if type(P_H2O) == str:
        P_H2O = P

    if type(MajorElements) == str:
        MajorElements = Shishkina_ME_MORB


    molarProps = molar_props(MajorElements)


    H2O = np.exp(a*np.log(P_H2O) +b*IM_NBO_O(molarProps,hydrous)+
                 B+C*P/T)

    return H2O

def IM_mixed_root(x,CO2,H2O,T=273+1200,MajorElements='default',hydrous=False):
    """
    Function to run root finder on for finding Saturation P for H2O + CO2
    saturation.

    Parameters
    ----------
    x:  list of floats
        [Pressure in bar, XCO2 as a fraction (volume/mole fraction CO2 in vapour phase.)]
    CO2:    float
        CO2 concentration in melt.
    H2O:    float
        H2O concentration in melt.
    T:  float
        Temperature in K.
    MajorElements:  series
        Series labelled with major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions.
    hydrous:    bool
        Use the hydrous or anhydrous calibration.
    """

    XCO2 = x[1]
    P = x[0]

    PCO2 = XCO2*P
    PH2O = (1-XCO2)*P

    if H2O != 0:
        H2O_calc = IM_H2O(P,T,MajorElements,PH2O,hydrous)
    else:
        H2O_calc = 0
    CO2_calc = IM_CO2(P,T,MajorElements,PCO2,H2O,hydrous)

    return [CO2-CO2_calc,H2O-H2O_calc]

def IM_SaturationP(CO2,H2O,T=273+1200,MajorElements='default',hydrous=False,
                   x0=np.array([1000,0.99])):
    """
    Returns the Saturation Pressure for a given CO2 and H2O concentration in melt,
    and a vapour composition. Solves Iacono-Marziano (2012) Eqs 7 and 8.

    Parameters
    ----------
    CO2:    float
        Concentration of CO2 in melt in ppm.
    H2O:    float
        Concentration of H2O in melt in wt%.
    MajorElements:  series
        Series labelled with major element oxides in wt%. All oxides passed will be
        included in calculating molar proportions.
    hydrous:    bool
        Use the hydrous or anhydrous calibration.

    Returns
    -------
    x:  np.array
        Array containing the solved Pressure and XCO2 in the coexisting vapour.

    """
    return root(IM_mixed_root,x0,args=(CO2,H2O,T,MajorElements,hydrous))['x']


def IM_contourPlot(P=4000,T=273+1200,MajorElements=IM_Etna):

    XCO2 = np.linspace(0.01,0.99,100)
    PH2O = (1-XCO2)*P
    PCO2 = XCO2*P

    CO2 = list()
    H2O = list()
    for i in range(len(XCO2)):
        H2O.append(IM_H2O(P,T,MajorElements,P_H2O=PH2O[i]))
        CO2.append(IM_CO2(P,T,MajorElements,P_CO2=PCO2[i],H2O=H2O[i]))


    return CO2,H2O


# Dixon (1997) Eq (3) reversed, in order to get XCO3 from CO2 ppm, as given by Shishkina
def Dixon_XCO3(CO2):
    """Dixon (1997) Eq (3).
    Calculates mole fraction carbonate ions dissolved in melt, given ppm CO2. Use for
    mixed fluid calculations when model outputs CO2 (ppm) directly, e.g. Shishkina (2014).

    Parameters
    ----------
    CO2:    float
        Concentration of CO2 in melt (in ppm)."""
    return (36.6e-4*CO2) / (4400-44e-4*CO2)

def Dixon_XB(H2O):
    """Dixon (1997) Eq (5).
    Calculates the mole fraction of total H2O dissolved in the melt, given H2O
    concentration in wt%.

    Parameters
    ----------
    H2O:    float
        Concentration of H2O in the melt (wt%)."""
    return (36.6*H2O)/(1801.5+18.6*H2O)


def Shishkina_H2O(P,MajorElements):
    """Shishkina et al. (2014) Eq (9).
    Calculates H2O concentration (wt%) in a melt in equilibrium with pure H2O
    vapour.

    Parameters
    ----------
    P:  float
        Pressure in bar.
    MajorElements:  Series
        Series of major elements with oxide labels and values in wt%. All oxides
        passed will be included in calculating molar proportions."""
    P=P/10
    MolarProps = MajorElements.copy()
    for ox in MolarProps.index:
        MolarProps[ox] = CationNum[ox]*MolarProps[ox]/oxideMass[ox]

    # Normalise
    MolarProps = MolarProps/MolarProps.sum()

    a = 3.36e-7 * P**3 - 2.33e-4*P**2 + 0.0711*P - 1.1309
    b = -1.2e-5*P**2 + 0.0196*P+1.1297

    return a*(MolarProps['Na2O']+MolarProps['K2O']) + b

def Dixon_XH2O_root(XH2O,XB):
    """Dixon (1997) Eq (4), substituted XOH for XB.

    Returns the difference between the two sides of equation (4) which relates
    the mole fraction of hydroxyl groups to the mole fraction of molecular water.
    For use in finding the root of the equation. Since both XH2O and XOH are unknown
    when starting from H2O wt% (or XB), XB is substituted for XOH.

    It is a regular solution model (Silver and Stolper, 1989).

    Parameters
    ----------
    XH2O: float
        Mole fraction of molecular water dissolved in melt.
    XB:     float
        XB mole fraction parameter, calculated from total H2O wt%."""
    A = 0.403
    B = 15.333
    C = 10.894

    lhs = ((2*(XB-XH2O))**2/(XH2O*(1-XH2O)))**(-1)
    rhs = np.exp(A + 2*B*(XB-XH2O)+C*XH2O)

    return rhs-lhs

def Dixon_XH2O(H2O):
    """Solves Dixon (1997) modified Eq (4).

    Finds the mole fraction of molecular water dissolved in the melt, given the
    mole fraction parameter XB. Uses the Newton method to find the root of
    modified Eq (4).

    Often struggles to find route when H2O -> 0.

    Parameters
    ----------
    H2O: float or list
        Concentration of water dissolved in melt (wt%)."""
    if isinstance(H2O,float):
        XB = (36.6*H2O)/(1801.5+18.6*H2O)
        return newton(Dixon_XH2O_root,XB/2,args=(XB,))
    else:
        results = list()
        for i in range(np.shape(H2O)[0]):
            XB = (36.6*H2O[i])/(1801.5+18.6*H2O[i])
            results.append(newton(Dixon_XH2O_root,XB/2,args=(XB,)))
        return np.array(results)


def Dixon_SaturationP_root(P,XCO3,XH2O,MajorElements='default'):
    """Based on Dixon (1997) Eq (11).

    Returns the difference between the lhs and rhs. Used for finding the saturation
    pressure of a given combination of CO2 and H2O when combined with a root finder.
    Calculates pressure dependence using the Dixon (1997) equations.

    Parameters
    ----------
    P:  float
        Pressure in bar
    XCO3:   float
        Molar proportion of carbonate molecule observed in the melt.
    XH2O:   float
        Molar proportion of water molecules observed in the melt.
    MajorElements:  series or str
        If a string, no compositional parameterisation is applied. Otherwise it must
        be a series of major element oxides in wt%, containing at least SiO2."""
    a = XH2O/Dixon_XH2O_pure(P,MajorElements=MajorElements)
    b = XCO3/Dixon_XCO3_pure(P,MajorElements=MajorElements)
    return a + b - 1

def Shishkina_SaturationP_root(P,XCO3,XH2O,MajorElements,PiStar=0.0):
    """Based on Dixon (1997) Eq (11).

    Returns the difference between the lhs and rhs. Used for finding the saturation
    pressure of a given combination of CO2 and H2O when combined with a root finder.
    Calculates pressure dependence using the Shishkina et al. (2014) equations.

    Parameters
    ----------
    P:  float
        Pressure in bar
    XCO3:   float
        Molar proportion of carbonate molecule observed in the melt.
    XH2O:   float
        Molar proportion of water molecules observed in the melt.
    PiStar:     float
        Ignored if set to 0.0. Pi* compositional parameter to use when calculating
        CO2 solubility. If ignored, Pi* is calculated using the passed MajorElements.
    MajorElements:  series
        A series of major element oxides in wt%, all oxides contained within it will
        be used for calculating total cations. """
    H2O_pure = Shishkina_H2O(P,MajorElements)
    XH2O_pure = Dixon_XH2O(H2O_pure)
    CO2_pure = Shishkina_CO2(P,PiStar,MajorElements)
    XCO3_pure = Dixon_XCO3(CO2_pure)

    return XH2O/XH2O_pure + XCO3/XCO3_pure - 1


def eguchi_Xi_melt(P,T=1473.0,MajorElements='default',species='CO3'):
    """
    Equation (9) in Eguchi and Dasgupta (2018). Calculates the mole fraction
    of dissolved molecular CO2 or carbonate CO3(2-).

    Parameters
    ----------
    P (float):  Pressure (bar)
    T (float):  Temperature (K)
    MajorElements (pandas Series):  Major Element Oxides.
    species (str):  Which species to calculate, molecular CO2 or carbonate ion.

    Returns
    -------
    Xi (float):     Mole fraction of selected species in the melt

    """

    if species == 'CO3':
        DH = -1.65e5
        DV = 2.38e-5
        DS = -43.64
        B = 1.47e3
        yNBO = 3.29
        A_CaO = 1.68e5
        A_Na2O = 1.76e5
        A_K2O = 2.11e5
    if species == 'CO2':
        DH = -9.02e4
        DV = 1.92e-5
        DS = -43.08
        B = 1.12e3
        yNBO = -7.09
        A_CaO = 0
        A_Na2O = 0
        A_K2O = 0
    R = 8.314

    if isinstance(MajorElements,str):
            MajorElements = Shishkina_ME_MORB

    # Calculate NBO term
    cations, O = molar_props(MajorElements,oxides=False,count_H2O=False,include_O=True)
    oxideprops =molar_props(MajorElements,oxides=True,count_H2O=False,include_O=False)

    NM = (cations['MgO'] + cations['CaO'] + cations['FeO'] + cations['Na2O'] +
          cations['K2O'] + cations['MnO'])
    Al = cations['Al2O3'] - NM
    if Al > 0:
        Al = NM
    else:
        Al = cations['Al2O3']
    Fe = cations['Fe2O3'] + Al
    if Al > 0:
        Fe = 0
    if Al < 0 and Fe > 0:
        Fe = - Al
    if Al < 0 and Fe < 0:
        Fe = cations['Fe2O3']
    Tet = cations['SiO2'] + cations['TiO2'] + cations['P2O5'] + Al + Fe
    NBO = 2*O - 4*Tet

#    lnfCO2 = np.log(float(PitzerSterner_fCO2(P,T)))
    lnfCO2 = np.log(zhang_fCO2(P,T))

    # Convert P to Pa
    P = P*1e5

    lnXi = ((DH/(R*T)-(P*DV)/(R*T)+DS/R) +
            (A_CaO*oxideprops['CaO']+A_Na2O*oxideprops['Na2O']+A_K2O*oxideprops['K2O'])/(R*T) +
            (B*lnfCO2/T) + yNBO*NBO
            )

    return np.exp(lnXi)

def eguchi_wt_convert(XCO2,XCO3,MajorElements='default'):
    """
    Eguchi and Dasgupta (2018) Eq (10). Converts the mole fraction molecular and
    carbonate ions into wt% CO2.

    Parameters
    ----------
    XCO2 (float): Mole fraction molecular CO2 dissolved in melt
    XCO3 (float): Mole fraction carbonate ion CO3 dissolved in melt
    MajorElements (pandas Series): Major element oxides.

    Returns
    -------
    Total CO2 (float):  Total CO2 dissolved in melt (wt%).
    Molecular CO2 (float):  CO2 dissolved as molecular CO2 (wt%).
    Carbonate CO2 (float):  CO2 dissolved as carbonate ion in melt (wt%).
    """
    if isinstance(MajorElements,str):
            MajorElements = Shishkina_ME_MORB

    FW_one = 15.999
    cations, O = molar_props(MajorElements,oxides=False,count_H2O=False,include_O=True)
    cations_one = cations/O
    for ox in cations_one.index:
        FW_one = FW_one + CationMass[ox]*cations_one[ox]

    CO2_CO2 = ((44.01*XCO2)/(44.01*XCO2+(1-(XCO2+XCO3))*FW_one))*100
    CO2_CO3 = ((44.01*XCO3)/(44.01*XCO3+(1-(XCO2+XCO3))*FW_one))*100

    return CO2_CO2+CO2_CO3, CO2_CO2, CO2_CO3

def eguchi_CO2(P,MajorElements='default',T=1473.0):
    XCO3 = eguchi_Xi_melt(P,T,MajorElements,species='CO3')
    XCO2 = eguchi_Xi_melt(P,T,MajorElements,species='CO2')
    CO2 = eguchi_wt_convert(XCO2,XCO3,MajorElements)[0]

    return CO2

def eguchi_CO2_tosolve(P,CO2,MajorElements,T):
    return eguchi_CO2(P,MajorElements,T)-CO2

def eguchi_SatP(CO2,MajorElements='default',T=1473.0):
    P = fsolve(eguchi_CO2_tosolve,1000,args=(CO2,MajorElements,T))
    return P

def convert_iron(MajorElements='default',ferric_total=0.15):
    if isinstance(MajorElements,str):
            MajorElements = Shishkina_ME_MORB
    majors = MajorElements.copy()
    Fe_t = majors['FeO']/oxideMass['FeO']
    Fe3 = ferric_total*Fe_t
    Fe2 = Fe_t - Fe3
    majors['FeO'] = Fe2*oxideMass['FeO']
    majors['Fe2O3'] = Fe3*oxideMass['Fe2O3']/2
    majors = majors/majors.sum()*100
    return majors

def zhang_Vm(Vm,P,T):
    Pm = 3.0636*P*3.79**3/235.0
    Tm = 154*T/235.0
    a = np.array([0.0,
                  2.95177298930e-2,
                  -6.33756452413e3,
                  -2.75265428882e5,
                  1.29128089283e-3,
                  -1.45797416153e2,
                  7.65938947237e4,
                  2.58661493537e-6,
                  0.52126532146,
                  -1.39839523753e2,
                  -2.36335007175e-8,
                  5.35026383543e-3,
                  -0.27110649951,
                  2.50387836486e4,
                  0.73226726041,
                  1.5483335997e-2])

    return ((1+(a[1]+a[2]/Tm**2+a[3]/Tm**3)/Vm+
             (a[4]+a[5]/Tm**2+a[6]/Tm**3)/Vm**2+
             (a[7]+a[8]/Tm**2+a[9]/Tm**3)/Vm**4)*0.08314*Tm/Pm - Vm
            )

def zhang_fCO2(P,T):
    a = np.array([0.0,
                  2.95177298930e-2,
                  -6.33756452413e3,
                  -2.75265428882e5,
                  1.29128089283e-3,
                  -1.45797416153e2,
                  7.65938947237e4,
                  2.58661493537e-6,
                  0.52126532146,
                  -1.39839523753e2,
                  -2.36335007175e-8,
                  5.35026383543e-3,
                  -0.27110649951,
                  2.50387836486e4,
                  0.73226726041,
                  1.5483335997e-2])
    e = 235.0
    s = 3.79

    Pm = 3.0636*P*s**3/e
    Tm = 154*T/e
    Vm = fsolve(zhang_Vm,200,args=(P,T))

    S1 = ((a[1]+a[2]/Tm**2+a[3]/Tm**3)/Vm+
          (a[4]+a[5]/Tm**2+a[6]/Tm**3)/(2*Vm**2)+
          (a[7]+a[8]/Tm**2+a[9]/Tm**3)/(4*Vm**4)+
          (a[10]+a[11]/Tm**2+a[12]/Tm**3)/(5*Vm**5)+
          (a[13]/(2*a[15]*Tm**3)*(a[14]+1-(a[14]+1+a[15]/Vm**2)*
           np.exp(-a[15]/Vm**2)))
         )
#    S2 = ((2*a[2]/Tm**2+3*a[3]/Tm**3)/Vm+
#          (2*a[5]/Tm**2+3*a[6]/Tm**3)/(2*Vm**2)+
#          (2*a[8]/Tm**2+3*a[9]/Tm**3)/(4*Vm**4)+
#          (2*a[11]/Tm**2+3*a[12]/Tm**3)/(5*Vm**5)+
#          (3*a[13]/(2*a[15]*Tm**3)*(a[14]+1-
#           (a[14]+1+a[15]/Vm**2)*np.exp(-a[15]/Vm**2)))
#            )

    Z = Pm*Vm/(0.08314*Tm)

    lnfc = Z - 1 - np.log(Z) + S1

    return P*np.exp(lnfc)


def eguchi_test():
    P = np.linspace(1000,30000,100)
    T = np.linspace(800+273,1600+273,100)
    ferric = np.linspace(0.01,0.3,100)

    alkali_basalt = pd.Series({'SiO2':42.4,'TiO2':4.1,'Al2O3':14.1,'Fe2O3':5.8,
                               'FeO':8.5,'MnO':0.2,'MgO':6.7,'CaO':11.9,'Na2O':2.8,
                               'K2O':2.0,'P2O5':0.6})
    tholeiite_basalt = pd.Series({'SiO2':53.8,'TiO2':2.0,'Al2O3':13.9,'Fe2O3':2.6,
                               'FeO':9.3,'MnO':0.2,'MgO':4.1,'CaO':7.9,'Na2O':3.0,
                               'K2O':1.5,'P2O5':1.4})
    ferr_bas = pd.Series({'SiO2':53.8,'TiO2':2.0,'Al2O3':13.9,
                               'FeO':11.5,'MnO':0.2,'MgO':4.1,'CaO':7.9,'Na2O':3.0,
                               'K2O':1.5,'P2O5':1.4})

    fCO2 = list()
    thol_CO2 = list()
    alk_CO2 = list()
    T_CO2 = list()
    ferr_CO2 = list()

    for i in range(len(P)):
        fCO2.append(zhang_fCO2(P[i],1473))
        thol_CO2.append(eguchi_CO2(P[i],T=1673,MajorElements=tholeiite_basalt)[0])
        alk_CO2.append(eguchi_CO2(P[i],T=1673,MajorElements=alkali_basalt)[0])

    for i in range(len(T)):
        T_CO2.append(eguchi_CO2(5000,T=T[i],MajorElements=tholeiite_basalt)[0])

    for i in range(len(ferric)):
        majors = convert_iron(ferr_bas,ferric[i])
        ferr_CO2.append(eguchi_CO2(5000,T=1673,MajorElements=majors)[0])

    f,a = plt.subplots(2,2)
    a= np.ravel(a).tolist()

    a[0].plot(P,fCO2)
    a[1].plot(P,thol_CO2,label='tholeiite basalt')
    a[1].plot(P,alk_CO2,label='alkali basalt')
    a[1].legend()
    a[2].plot(T,T_CO2)
    a[3].plot(ferric,ferr_CO2)

    plt.show()


def Allison_CO2Sat(P,T=1200,FWone=36.594,loc='sunset',fit='power'):
    if fit == 'thermodynamic':
        P0 = 1000 # bar
        params = dict({'sunset':[16.4,-14.67],
                       'sfvf':[15.02,-14.87],
                       'erebus':[15.83,-14.65],
                       'vesuvius':[24.42,-14.04],
                       'etna':[21.59,-14.28],
                       'stromboli':[14.93,-14.68]})
        DV = params[loc][0]
        lnK0 = params[loc][1]

        lnK = lnK0 - (P-P0)*DV/(8.3141*(T+273.15))
        print(lnK)
        fCO2 = zhang_fCO2(P,T+273.15)
        print(fCO2)
        Kf = np.exp(lnK)*fCO2
        print(Kf)
        XCO3 = Kf/(1-Kf)
        print(XCO3)
        wtCO2 = (44.01*XCO3)/((44.01*XCO3)+(1-XCO3)*FWone)*100

        return wtCO2
    if fit == 'power':
        params = dict({'stromboli':[1.05,0.883],
                       'etna':[2.831,0.797],
                       'vesuvius':[4.796,0.754],
                       'sfvf':[3.273,0.74],
                       'sunset':[4.32,0.728],
                       'erebus':[5.145,0.713]})

        fCO2 = zhang_fCO2(P,T+273.15)

        return params[loc][0]*fCO2**params[loc][1]

def Allison_CO2Sat_tosolve(P,CO2,T=1200,loc='sunset',fit='power'):
    predicted = Allison_CO2Sat(P,T,loc=loc,fit=fit)
    return (CO2-predicted)**2

def Allison_SatP(CO2,T=1200,loc='sunset',fit='power'):
    P = fsolve(Allison_CO2Sat_tosolve,1000,args=(CO2,T,loc,fit))
    return P



def SaturationP(CO2,H2O=0.0,MajorElements='default',PiStar=0.0,method='Dixon',x0=1000,
                Allison_loc='sunset',suppress_warnings=False):
    """Finds the saturation pressure (bar), given a CO2 and H2O melt content.

    Parameters
    ----------
    CO2:    float
        CO2 concentration (ppm) in the melt.
    H2O:    float
        H2O concentration (wt%) in the melt.
    MajorElements:  series or str
        For Dixon method, if set to str, no compositional parameters will be used.
        For Shishkina method, it must not be a str. It should be a series of major
        element oxides in wt%, including the minimum terms required for the relevant
        parameterisation.
    PiStar:     float
        Ignored if set to 0.0. If non-zero it will set the Pi* parameter when calculating
        CO2 solubility, rather than calculating it from the MajorElements object.
    Method:     float
        Dixon or Shishkina. Sets which set of solubility equations to use. Shishkina
        uses the Dixon equations to obtain molar proportions of carbonate ion and molecular
        water dissolved in the melt, from the H2O and CO2 concentrations.
    x0:     float
        Starting guess for root-finder."""

    if method == 'Dixon':
        XCO3 = Dixon_XCO3(CO2)
        if H2O > 0.0:
            XH2O = Dixon_XH2O(H2O)
        else:
            XH2O = 0.0
        try:
            return fsolve(Dixon_SaturationP_root,x0,args=(XCO3,XH2O,MajorElements))[0]
        except:
            return 0.0
    elif method == 'Shishkina':
        if isinstance(MajorElements,str):
            MajorElements = Shishkina_ME_MORB
        if H2O == 0:
            if PiStar == 0:
                PiStar = Shishkina_PiStar(MajorElements)
            return np.exp((np.log(CO2)-6.71*PiStar+1.345)/1.15)*10
        XCO3 = Dixon_XCO3(CO2)
        XH2O = Dixon_XH2O(H2O)
        try:
            return fsolve(Shishkina_SaturationP_root,x0,args=(XCO3,XH2O,MajorElements,PiStar))[0]
        except:
            try:
                fsolve(Shishkina_SaturationP_root,x0,args=(XCO3,0.0,MajorElements,PiStar))[0]
            except:
                return 0.0
    elif method =='Eguchi':
        if H2O != 0.0 and suppress_warnings==False:
            print('Warning: H2O not implemented in Eguchi parameterisation')
        return eguchi_SatP(CO2/1e4,MajorElements)[0]
    elif method =='Iacono-Marziano':
        return IM_SaturationP(CO2,H2O,MajorElements=MajorElements)[0]
    elif method =='Allison':
        if H2O != 0.0 and suppress_warnings==False:
            print('Warning: H2O not implemented in Allison parameterisation')
        return Allison_SatP(CO2,loc=Allison_loc)[0]




def MixedFluidsSaturation(P_list,T=1473,method='Dixon',MajorElements='Default'):
    '''Calculates H2O and CO2 concentration in melts at saturation. Returns a list of
    CO2 concentrations and a list of H2O concentrations for each pressure.

    Parameters
    ----------
    P_list: list of floats
        List of pressures at which to calculate CO2 and H2O saturation.

    T: float
        Temperature in K
    method:     float
        Either Shishkina or Dixon. Determines which parametrization of solubility
        to use. To calculate mixed fluids using the Shishkina equations, the equations
        relating CO2 ppm and XCO3, and H2O wt%, XOH and XH2O, by Dixon (1997) are
        used.
    composition:    Series
        Series of major element oxide concentrations in wt% with their oxides as
        labels. If set as a string the composition will be set to the MORB composition
        used by Shishkina (2014), or default if using Dixon (1997).
        '''
    results_H2O = list()
    results_CO2 = list()
    for P in P_list:
        if method == 'Dixon':
            XH2O_max = Dixon_XH2O_pure(P)
            XH2O = np.linspace(1e-5, XH2O_max,1000)
            XCO3_max = Dixon_XCO3_pure(P,MajorElements=MajorElements)
            XCO3 = XCO3_mixed(XH2O,XH2O_max,XCO3_max)
            XOH = Dixon_XOH_solve(XH2O)

            results_H2O.append(Dixon_H2O_wt(XH2O,XOH))
            results_CO2.append(Dixon_CO2_ppm(XCO3))
        elif method == 'Shishkina':
            if isinstance(MajorElements,str):
                MajorElements = Shishkina_ME_MORB
            H2O_max = Shishkina_H2O(P,MajorElements)
            XH2O_max = Dixon_XH2O([H2O_max])
            CO2_max = Shishkina_CO2(P,MajorElements=MajorElements)
            H2O = np.linspace(1e-2, H2O_max,1000)
            XH2O = Dixon_XH2O(H2O)
            CO2 = Dixon_CO2_ppm(XCO3_mixed(XH2O,XH2O_max,Dixon_XCO3(CO2_max)))

            results_H2O.append(H2O)
            results_CO2.append(CO2)

        else:
            print('Method not recognised')

    return results_H2O,results_CO2

def CO2(P,MajorElements='default',PiStar=0.0,method='Dixon',Allison_loc='sunset'):
    """Wrapper function for calculating CO2 content of a melt in equilibrium with
    a pure CO2 vapour.

    Parameters
    ----------
    P:  float
        Pressure (bar)
    MajorElements:  series
        Major element oxides in wt%. If using the Dixon method, a string will
        result in the non-compositional dependent parameterisation.
    PiStar:     float
        Pi* compositional paramter. If zero, it will be ignored, and instead
        calculated from MajorElements. Used only in the Shishkina model.
    method:     string
        Parameterization (Dixon, Shishkina, Eguchi, Iacono-Marziano or Allison)
        to use when calculating solubility.
    """
    if method == 'Dixon':
        return Dixon_CO2(P,MajorElements)
    elif method == 'Shishkina':
        return Shishkina_CO2(P,PiStar,MajorElements)
    elif method == 'Eguchi':
        return eguchi_CO2(P,MajorElements)
    elif method == 'Iacono-Marziano':
        return IM_CO2(P,MajorElements)
    elif method == 'Allison':
        return Allison_CO2Sat(P,loc=Allison_loc)
    else:
        "Method variable not recognised. Use Dixon, Shishkina, Eguchi, \
         Iacono-Marziano or Allison"

def H2O(P,MajorElements='default',method='Dixon'):
    """Wrapper function for calculating H2O content of a melt in equilibrium with
    a pure H2O vapour.

    Parameters
    ----------
    P:  float
        Pressure (bar)
    MajorElements:  series
        Major element oxides in wt%. If using the Dixon method, a string will
        result in the non-compositional dependent parameterisation.
    method:     string
        Parameterization (Dixon or Shishkina) to use when calculating solubility.
    """
    if method == 'Dixon':
        return Dixon_H2O(P,MajorElements)
    elif method == 'Shishkina':
        return Shishkina_H2O(P,MajorElements)




#H2O_list,CO2_list = MixedFluidsSaturation([1000.0,5000.0],1473,'Dixon')

#fig = plt.figure()
#ax = plt.gca()
#
#for H2O,CO2 in zip(H2O_list,CO2_list):
#    ax.plot(H2O,CO2)
#ax.set_xlabel('H2O (wt%)')
##ax.set_ylim(0,2000)
##ax.set_xlim([0,3.0])
#ax.set_ylabel('CO2 (ppm)')
#plt.show()
