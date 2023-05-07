from VESIcal import calibration_checks
from VESIcal import core

from scipy.optimize import root_scalar
from abc import abstractmethod
import numpy as np


class FugacityModel(object):
    """ The fugacity model object is for implementations of fugacity models
    for individual volatile species, though it may depend on the mole
    fraction of other volatile species. It contains all the methods required
    to calculate the fugacity at a given pressure and mole fraction.
    """

    def __init__(self):
        self.set_calibration_ranges([])

    def set_calibration_ranges(self, calibration_ranges):
        self.calibration_ranges = calibration_ranges

    @abstractmethod
    def fugacity(self, pressure, **kwargs):
        """
        """

    # @abstractmethod
    def check_calibration_range(self, parameters, report_nonexistance=True):
        s = ''
        for cr in self.calibration_ranges:
            if cr.check(parameters) is False:
                s += cr.string(parameters, report_nonexistance)
        return s


# ------------- FUGACITY MODELS -------------------------------- #

class fugacity_idealgas(FugacityModel):
    """ An instance of FugacityModel for an ideal gas.
    """

    def fugacity(self, pressure, X_fluid=1.0, **kwargs):
        """ Returns the fugacity of an ideal gas, i.e., the partial pressure.

        Parameters
        ----------
        pressure    float
            Total pressure of the system, in bars.
        X_fluid     float
            The mole fraction of the species in the vapour phase.

        Returns
        -------
        float
            Fugacity (partial pressure) in bars
        """
        return pressure*X_fluid


class fugacity_KJ81_co2(FugacityModel):
    """ Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids. This class
    will return the properties of the CO2 component of the mixed fluid.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', 20000.0,  calibration_checks.crf_LessThan, 'bar',
                'Kerrick and Jacobs (1981) EOS', fail_msg=calibration_checks.crmsg_LessThan_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'temperature', 1050, calibration_checks.crf_LessThan, 'oC',
                'Kerrick and Jacobs (1981) EOS', fail_msg=calibration_checks.crmsg_LessThan_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description)])

    def fugacity(self, pressure, temperature, X_fluid, **kwargs):
        """ Calculates the fugacity of CO2 in a mixed CO2-H2O fluid. Above 1050C, it assumes H2O
        and CO2 do not interact, as the equations are not defined beyond this point.

        Parameters
        ----------
        pressure    float
            Total pressure of the system in bars.
        temperature     float
            Temperature in degC
        X_fluid     float
            Mole fraction of CO2 in the fluid.

        Returns
        -------
        float
            fugacity of CO2 in bars
        """
        if X_fluid == 0:
            return 0
        elif temperature >= 1050.0:
            return pressure*np.exp(self.lnPhi_mix(pressure, temperature, 1.0))*X_fluid
        else:
            return pressure*np.exp(self.lnPhi_mix(pressure, temperature, X_fluid))*X_fluid

    def volume(self, P, T, X_fluid):
        """ Calculates the volume of the mixed fluid, by solving Eq (28) of Kerrick and
        Jacobs (1981) using scipy.root_scalar.

        Parameters
        ----------
        P   float
            Total pressure of the system, in bars.
        T   float
            Temperature in degC
        X_fluid     float
            Mole fraction of CO2 in the fluid

        Returns
        -------
        float
            Volume of the mixed fluid.
        """
        if X_fluid != 1.0:
            # x0 = self.volume(P,T,1.0)*X_fluid + self.volume_h(P,T)*(1-X_fluid)
            # print(x0)
            if P >= 20000 and T < 800-273.15:
                x0 = (X_fluid*25+(1-X_fluid)*15)
            else:
                x0 = (X_fluid*35+(1-X_fluid)*15)

        else:
            if P >= 20000 and T < 800-273.15:
                x0 = 25
            else:
                x0 = 35
        return root_scalar(self.root_volume, x0=x0, x1=x0*0.9, args=(P, T, X_fluid)).root

    def root_volume(self, v, P, T, X_fluid):
        """ Returns the difference between the lhs and rhs of Eq (28) of Kerrick and Jacobs (1981).
        For use with a root finder to obtain the volume of the mixed fluid.

        Parameters
        ----------
        v   float
            Guess for the volume
        P   float
            Total system pressure in bars.
        T   float
            Temperature in degC
        X_fluid     float
            Mole fraction of CO2 in the fluid.

        Returns
        -------
        float
            Difference between lhs and rhs of Eq (28) of Kerrick and Jacobs (1981), in bars.
        """
        T = T + 273.15
        c = {}
        h = {}

        c['b'] = 58.0
        c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
        c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
        c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
        h['b'] = 29.0
        h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6
        h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
        h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

        if X_fluid == 1:
            bm = c['b']
            cm = c['c']
            c12 = c['c']
            dm = c['d']
            d12 = c['d']
            em = c['e']
            e12 = c['e']
        else:
            bm = X_fluid*c['b'] + (1-X_fluid)*h['b']
            c12 = (c['c']*h['c'])**0.5
            cm = c['c']*X_fluid**2 + h['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
            d12 = (c['d']*h['d'])**0.5
            dm = c['d']*X_fluid**2 + h['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
            e12 = (c['e']*h['e'])**0.5
            em = c['e']*X_fluid**2 + h['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

        am = cm + dm/v + em/v**2

        y = bm/(4*v)

        pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
        pt2 = - am / (T**0.5 * v * (v+bm))

        return -(P - pt1 - pt2)

    def volume_h(self, P, T):
        """ Calculates the volume of a pure H2O fluid, by solving Eq (14) of
        Kerrick and Jacobs (1981).

        Parameters
        ----------
        P   float
            Total pressure in bars.
        T   float
            Temperature in degC.

        Returns
        -------
        Difference between lhs and rhs of Eq (14) of Kerrick and Jacobs (1981), in bars.
        """
        return root_scalar(self.root_volume_h, x0=15, x1=35, args=(P, T)).root

    def root_volume_h(self, v, P, T):
        """ Returns the difference between the lhs and rhs of Eq (14) of
        Kerrick and Jacobs (1981). For use with a root solver to identify the
        volume of a pure H2O fluid.

        Parameters
        ----------
        v   float
            Guess for the volume
        P   float
            Total pressure in bars.
        T   float
            Temperature in degC.

        Returns
        -------
        float
            The difference between the lhs and rhs of Eq (14) of Kerrick and Jacobs (1981),
            in bars.
        """
        T = T + 273.15
        h = {}
        h['b'] = 29.0
        h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6
        h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
        h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6
        h['a'] = h['c'] + h['d']/v + h['e']/v**2

        y = h['b']/(4*v)

        pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
        pt2 = - h['a'] / (T**0.5 * v * (v+h['b']))

        return -(P - pt1 - pt2)

    def lnPhi_mix(self, P, T, X_fluid):
        """ Calculates the natural log of the fugacity coefficient for CO2 in a
        mixed CO2-H2O fluid. Uses Eq (27) of Kerrick and Jacobs (1981).

        Parameters
        ----------
        P   float
            Total pressure in bars.
        T   float
            Temperature in degC
        X_fluid     float
            The mole fraction of CO2 in the fluid.

        Returns
        -------
        float
            The natural log of the fugacity coefficient for CO2 in a mixed fluid.
        """
        T = T + 273.15
        v = self.volume(P, T-273.15, X_fluid)

        c = {}
        h = {}

        c['b'] = 58.0
        c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
        c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
        c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
        h['b'] = 29.0
        h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6
        h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
        h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

        if X_fluid == 1:
            bm = c['b']
            cm = c['c']
            c12 = c['c']
            dm = c['d']
            d12 = c['d']
            em = c['e']
            e12 = c['e']
        else:
            bm = X_fluid*c['b'] + (1-X_fluid)*h['b']
            c12 = (c['c']*h['c'])**0.5
            cm = c['c']*X_fluid**2 + h['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
            d12 = (c['d']*h['d'])**0.5
            dm = c['d']*X_fluid**2 + h['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
            e12 = (c['e']*h['e'])**0.5
            em = c['e']*X_fluid**2 + h['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

        y = bm/(4*v)

        Z = v*P/(83.14*T)

        lnPhi = 0

        lnPhi += (4*y-3*y**2)/(1-y)**2 + (c['b']/bm * (4*y-2*y**2)/(1-y)**3)
        lnPhi += - (2*c['c']*X_fluid+2*(1-X_fluid)*c12)/(83.14*T**1.5*bm)*np.log((v+bm)/v)
        lnPhi += - cm*c['b']/(83.14*T**1.5*bm*(v+bm))
        lnPhi += cm*c['b']/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
        lnPhi += - (2*c['d']*X_fluid+2*d12*(1-X_fluid)+dm)/(83.14*T**1.5*bm*v)
        lnPhi += (2*c['d']*X_fluid+2*(1-X_fluid)*d12+dm)/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
        lnPhi += c['b']*dm/(83.14*T**1.5*v*bm*(v+bm)) + 2*c['b']*dm/(83.14*T**1.5*bm**2*(v+bm))
        lnPhi += - 2*c['b']*dm/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
        lnPhi += - (2*c['e']*X_fluid + 2*(1-X_fluid)*e12+2*em)/(83.14*T**1.5*2*bm*v**2)
        lnPhi += (2*c['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**2*v)
        lnPhi += - (2*c['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
        lnPhi += (em*c['b']/(83.14*T**1.5*2*bm*v**2*(v+bm)) -
                  3*em*c['b']/(83.14*T**1.5*2*bm**2*v*(v+bm)))
        lnPhi += (3*em*c['b']/(83.14*T**1.5*bm**4)*np.log((v+bm)/v) -
                  3*em*c['b']/(83.14*T**1.5*bm**3*(v+bm)))
        lnPhi += - np.log(Z)

        return lnPhi


class fugacity_KJ81_h2o(FugacityModel):
    """Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids. This class
    will return the properties of the H2O component of the mixed fluid.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', 20000.0, calibration_checks.crf_LessThan, 'bar',
                'Kerrick and Jacobs (1981) EOS',
                fail_msg=calibration_checks.crmsg_LessThan_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'temperature', 1050, calibration_checks.crf_LessThan, 'oC',
                'Kerrick and Jacobs (1981) EOS',
                fail_msg=calibration_checks.crmsg_LessThan_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description)])

    def fugacity(self, pressure, temperature, X_fluid, **kwargs):
        """ Calculates the fugacity of H2O in a mixed CO2-H2O fluid. Above 1050C,
        it assumes H2O and CO2 do not interact, as the equations are not defined
        beyond this point.

        Parameters
        ----------
        pressure    float
            Total pressure of the system in bars.
        temperature     float
            Temperature in degC
        X_fluid     float
            Mole fraction of H2O in the fluid.

        Returns
        -------
        float
            fugacity of H2O in bars
        """
        if X_fluid == 0:
            return 0
        elif temperature >= 1050:
            return pressure*np.exp(self.lnPhi_mix(pressure, temperature, 1.0))*X_fluid
        else:
            return pressure*np.exp(self.lnPhi_mix(pressure, temperature, X_fluid))*X_fluid

    def volume(self, P, T, X_fluid):
        """ Calculates the volume of the mixed fluid, by solving Eq (28) of Kerrick and
        Jacobs (1981) using scipy.root_scalar.

        Parameters
        ----------
        P   float
            Total pressure of the system, in bars.
        T   float
            Temperature in degC
        X_fluid     float
            Mole fraction of H2O in the fluid

        Returns
        -------
        float
            Volume of the mixed fluid.
        """
        if X_fluid != 1.0:
            if P >= 20000 and T < 800-273.15:
                x0 = ((1-X_fluid)*25+X_fluid*15)
            else:
                x0 = ((1-X_fluid)*35+X_fluid*15)

        else:
            if P >= 20000 and T < 800-273.15:
                x0 = 10
            else:
                x0 = 15
        return root_scalar(self.root_volume, x0=x0, x1=x0*0.9, args=(P, T, X_fluid)).root

    def root_volume(self, v, P, T, X_fluid):
        """ Returns the difference between the lhs and rhs of Eq (28) of Kerrick and Jacobs (1981).
        For use with a root finder to obtain the volume of the mixed fluid.

        Parameters
        ----------
        v   float
            Guess for the volume
        P   float
            Total system pressure in bars.
        T   float
            Temperature in degC
        X_fluid     float
            Mole fraction of H2O in the fluid.

        Returns
        -------
        float
            Difference between lhs and rhs of Eq (28) of Kerrick and Jacobs (1981), in bars.
        """
        T = T + 273.15
        c = {}
        h = {}

        c['b'] = 58.0
        c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
        c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
        c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
        h['b'] = 29.0
        h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6
        h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
        h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

        if X_fluid == 1:
            bm = h['b']
            cm = h['c']
            dm = h['d']
            em = h['e']
            c12 = h['c']
            d12 = h['d']
            e12 = h['e']
        else:
            bm = X_fluid*h['b'] + (1-X_fluid)*c['b']
            c12 = (c['c']*h['c'])**0.5
            cm = h['c']*X_fluid**2 + c['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
            d12 = (c['d']*h['d'])**0.5
            dm = h['d']*X_fluid**2 + c['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
            e12 = (c['e']*h['e'])**0.5
            em = h['e']*X_fluid**2 + c['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

        am = cm + dm/v + em/v**2

        y = bm/(4*v)

        pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
        pt2 = - am / (T**0.5 * v * (v+bm))

        return -(P - pt1 - pt2)

    def volume_c(self, P, T):
        """ Calculates the volume of a pure CO2 fluid, by solving Eq (14) of
        Kerrick and Jacobs (1981).

        Parameters
        ----------
        P   float
            Total pressure in bars.
        T   float
            Temperature in degC.

        Returns
        -------
        Difference between lhs and rhs of Eq (14) of Kerrick and Jacobs (1981), in bars.
        """
        return root_scalar(self.root_volume_c, x0=15, x1=35, args=(P, T)).root

    def root_volume_c(self, v, P, T):
        """ Returns the difference between the lhs and rhs of Eq (14) of
        Kerrick and Jacobs (1981). For use with a root solver to identify the
        volume of a pure H2O fluid.

        Parameters
        ----------
        v   float
            Guess for the volume
        P   float
            Total pressure in bars.
        T   float
            Temperature in degC.

        Returns
        -------
        float
            The difference between the lhs and rhs of Eq (14) of Kerrick and Jacobs (1981),
            in bars.
        """
        T = T + 273.15
        c = {}
        c['b'] = 58.0
        c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
        c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
        c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
        c['a'] = c['c'] + c['d']/v + c['e']/v**2

        y = c['b']/(4*v)

        pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
        pt2 = - c['a'] / (T**0.5 * v * (v+c['b']))

        return -(P - pt1 - pt2)

    def lnPhi_mix(self, P, T, X_fluid):
        """ Calculates the natural log of the fugacity coefficient for H2O in a
        mixed CO2-H2O fluid. Uses Eq (27) of Kerrick and Jacobs (1981).

        Parameters
        ----------
        P   float
            Total pressure in bars.
        T   float
            Temperature in degC
        X_fluid     float
            The mole fraction of H2O in the fluid.

        Returns
        -------
        float
            The natural log of the fugacity coefficient for H2O in a mixed fluid.
        """
        T = T + 273.15
        v = self.volume(P, T-273.15, X_fluid)

        c = {}
        h = {}

        c['b'] = 58.0
        c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
        c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
        c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
        h['b'] = 29.0
        h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6
        h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
        h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

        if X_fluid == 1:
            bm = h['b']
            cm = h['c']
            dm = h['d']
            em = h['e']
            c12 = h['c']
            d12 = h['d']
            e12 = h['e']
        else:
            bm = X_fluid*h['b'] + (1-X_fluid)*c['b']
            c12 = (c['c']*h['c'])**0.5
            cm = h['c']*X_fluid**2 + c['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
            d12 = (c['d']*h['d'])**0.5
            dm = h['d']*X_fluid**2 + c['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
            e12 = (c['e']*h['e'])**0.5
            em = h['e']*X_fluid**2 + c['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

        y = bm/(4*v)

        # Z = (1+y+y**2-y**3)/(1-y)**2 - am/(83.14*T**1.5*(v+bm))
        Z = v*P/(83.14*T)

        lnPhi = 0

        lnPhi += (4*y-3*y**2)/(1-y)**2 + (h['b']/bm * (4*y-2*y**2)/(1-y)**3)
        lnPhi += - (2*h['c']*X_fluid+2*(1-X_fluid)*c12)/(83.14*T**1.5*bm)*np.log((v+bm)/v)
        lnPhi += - cm*h['b']/(83.14*T**1.5*bm*(v+bm))
        lnPhi += cm*h['b']/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
        lnPhi += - (2*h['d']*X_fluid+2*d12*(1-X_fluid)+dm)/(83.14*T**1.5*bm*v)
        lnPhi += (2*h['d']*X_fluid+2*(1-X_fluid)*d12+dm)/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
        lnPhi += h['b']*dm/(83.14*T**1.5*v*bm*(v+bm)) + 2*h['b']*dm/(83.14*T**1.5*bm**2*(v+bm))
        lnPhi += - 2*h['b']*dm/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
        lnPhi += - (2*h['e']*X_fluid + 2*(1-X_fluid)*e12+2*em)/(83.14*T**1.5*2*bm*v**2)
        lnPhi += (2*h['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**2*v)
        lnPhi += - (2*h['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
        lnPhi += (em*h['b']/(83.14*T**1.5*2*bm*v**2*(v+bm)) -
                  3*em*h['b']/(83.14*T**1.5*2*bm**2*v*(v+bm)))
        lnPhi += (3*em*h['b']/(83.14*T**1.5*bm**4)*np.log((v+bm)/v) -
                  3*em*h['b']/(83.14*T**1.5*bm**3*(v+bm)))
        lnPhi += - np.log(Z)

        return lnPhi


class fugacity_ZD09_co2(FugacityModel):
    """ Implementation of the Zhang and Duan (2009) fugacity model for pure CO2
    fluids."""
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar',
                'Zhang and Duan (2009) EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', [200, 2300], calibration_checks.crf_Between, 'oC',
                'Zhang and Duan (2009) EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description)])

    def fugacity(self, pressure, temperature, X_fluid=1.0, **kwargs):
        """ Calculates the fugacity of a pure CO2 fluid, or a mixed fluid assuming
        ideal mixing. Implements eqn (14) of Zhang and Duan (2009).

        Parameters
        ---------
        pressure     float
            Pressure in bars
        temperature     float
            Temperature in degC
        X_fluid     float
            Mole fraction of CO2 in the fluid. Default is 1.0.

        Returns
        -------
        float
            Fugacity of CO2, standard state 1 bar.
        """

        P = pressure/10
        T = temperature + 273.15

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
        Vm = root_scalar(self.Vm, x0=200, x1=100, args=(P, T)).root

        S1 = ((a[1]+a[2]/Tm**2+a[3]/Tm**3)/Vm +
              (a[4]+a[5]/Tm**2+a[6]/Tm**3)/(2*Vm**2) +
              (a[7]+a[8]/Tm**2+a[9]/Tm**3)/(4*Vm**4) +
              (a[10]+a[11]/Tm**2+a[12]/Tm**3)/(5*Vm**5) +
              (a[13]/(2*a[15]*Tm**3)*(a[14]+1-(a[14]+1+a[15]/Vm**2) *
               np.exp(-a[15]/Vm**2)))
              )

        Z = Pm*Vm/(8.314*Tm)

        lnfc = Z - 1 - np.log(Z) + S1

        return P*np.exp(lnfc)*10

    def Vm(self, Vm, P, T):
        """ Function to use for solving for the parameter Vm, defined by eqn (8) of
        Zhang and Duan (2009). Called by scipy.fsolve in the fugacity method.

        Parameters
        ----------
        Vm     float
            Guessed value of Vm
        P     float
            Pressure in MPa
        T     float
            Temperature in K

        Returns
        -------
        float
            Difference between (rearranged) LHS and RHS of eqn (8) of Zhang and Duan (2009).
        """
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

        return ((1+(a[1]+a[2]/Tm**2+a[3]/Tm**3)/Vm +
                 (a[4]+a[5]/Tm**2+a[6]/Tm**3)/Vm**2 +
                 (a[7]+a[8]/Tm**2+a[9]/Tm**3)/Vm**4)*0.08314*Tm/Pm - Vm
                )


class fugacity_MRK_co2(FugacityModel):
    """ Modified Redlick Kwong fugacity model as used by VolatileCalc. Python implementation by
    D. J. Rasmussen (github.com/DJRgeoscience/VolatileCalcForPython), based on VB code by Newman &
    Lowenstern.
    """
    def __init__(self):
        self.set_calibration_ranges([])

    def fugacity(self, pressure, temperature, X_fluid=1.0, **kwargs):
        """ Calculates the fugacity of CO2 in a pure or mixed H2O-CO2 fluid (assuming ideal
        mixing).

        Parameters
        ----------
        pressure    float
            Total pressure of the system in bars.
        temperature     float
            Temperature in degC
        X_fluid     float
            Mole fraction of CO2 in the fluid.

        Returns
        -------
        float
            fugacity of CO2 in bars
        """
        fug = self.MRK(pressure, temperature+273.15)
        return fug*X_fluid

    def FNA(self, TK):
        return ((166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2
                - 0.071288 * ((TK - 273.15)**3)) * 1.01325)

    def FNB(self, TK):
        return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)

    def FNC(self, TK):
        R = 83.14321
        return (1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3)
                * 0.5 * R * R * TK**2.5 / 1.02668 + 40123800))

    def FNF(self, V, TK, A, B, P):
        R = 83.14321
        return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

    def MRK(self, P, TK):  # Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
        R = 83.14321
        B_1 = 14.6
        B_2 = 29.7

        for X_1 in [0, 1]:
            B = X_1 * B_1 + (1 - X_1) * B_2
            A = (X_1**2 * self.FNA(TK) + 2 * X_1 * (1 - X_1) * self.FNC(TK) +
                 (1 - X_1)**2 * self.FNB(TK))
            Temp2 = B + 5
            Q = 1
            Temp1 = 0
            while abs(Temp2 - Temp1) >= 0.00001:
                Temp1 = Temp2
                F_1 = (self.FNF(Temp1 + 0.01, TK, A, B, P) - self.FNF(Temp1, TK, A, B, P)) / 0.01
                Temp2 = Temp1 - Q * self.FNF(Temp1, TK, A, B, P) / F_1
                F_2 = (self.FNF(Temp2 + 0.01, TK, A, B, P) - self.FNF(Temp2, TK, A, B, P)) / 0.01
                if F_2 * F_1 <= 0:
                    Q = Q / 2.
                if abs(Temp2 - Temp1) > 0.00001:
                    F_1 = F_2
            V = Temp2
            G_1 = (np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * self.FNA(TK) +
                   (1 - X_1) * self.FNC(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B))
            G_1 = (G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) -
                   np.log(P * V / (R * TK)))
            G_1 = np.exp(G_1)
            G_2 = (np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * self.FNC(TK) +
                   (1 - X_1) * self.FNB(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B))
            G_2 = (G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) -
                   np.log(P * V / (R * TK)))
            G_2 = np.exp(G_2)
            if X_1 == 0:
                fCO2o = G_2 * P  # The fugacity of CO2
        # return fCO2o
        return fCO2o


class fugacity_MRK_h2o(FugacityModel):
    """ Modified Redlick Kwong fugacity model as used by VolatileCalc. Python implementation by
    D. J. Rasmussen (github.com/DJRgeoscience/VolatileCalcForPython), based on VB code by Newman &
    Lowenstern.
    """
    def __init__(self):
        self.set_calibration_ranges([])

    def fugacity(self, pressure, temperature, X_fluid=1.0, **kwargs):
        """ Calculates the fugacity of H2O in a pure or mixed H2O-CO2 fluid (assuming ideal
        mixing).

        Parameters
        ----------
        pressure    float
            Total pressure of the system in bars.
        temperature     float
            Temperature in degC
        X_fluid     float
            Mole fraction of H2O in the fluid.

        Returns
        -------
        float
            fugacity of H2O in bars
        """
        fug = self.MRK(pressure, temperature+273.15)
        return fug*X_fluid

    def FNA(self, TK):
        return ((166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2 -
                0.071288 * ((TK - 273.15)**3)) * 1.01325)

    def FNB(self, TK):
        return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)

    def FNC(self, TK):
        R = 83.14321
        return (1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3) *
                0.5 * R * R * TK**2.5 / 1.02668 + 40123800))

    def FNF(self, V, TK, A, B, P):
        R = 83.14321
        return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

    def MRK(self, P, TK):  # Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
        R = 83.14321
        B_1 = 14.6
        B_2 = 29.7

        # X_1 = 1
        for X_1 in [0, 1]:
            B = X_1 * B_1 + (1 - X_1) * B_2
            A = (X_1**2 * self.FNA(TK) + 2 * X_1 * (1 - X_1) * self.FNC(TK) +
                 (1 - X_1)**2 * self.FNB(TK))
            Temp2 = B + 5
            Q = 1
            Temp1 = 0
            while abs(Temp2 - Temp1) >= 0.00001:
                Temp1 = Temp2
                F_1 = (self.FNF(Temp1 + 0.01, TK, A, B, P) - self.FNF(Temp1, TK, A, B, P)) / 0.01
                Temp2 = Temp1 - Q * self.FNF(Temp1, TK, A, B, P) / F_1
                F_2 = (self.FNF(Temp2 + 0.01, TK, A, B, P) - self.FNF(Temp2, TK, A, B, P)) / 0.01
                if F_2 * F_1 <= 0:
                    Q = Q / 2.
                if abs(Temp2 - Temp1) > 0.00001:
                    F_1 = F_2
            V = Temp2
            G_1 = (np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * self.FNA(TK) +
                   (1 - X_1) * self.FNC(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B))
            G_1 = (G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) -
                   np.log(P * V / (R * TK)))
            G_1 = np.exp(G_1)
            G_2 = (np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * self.FNC(TK) + (1 - X_1)
                   * self.FNB(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B))
            G_2 = (G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) -
                   np.log(P * V / (R * TK)))
            G_2 = np.exp(G_2)
            if X_1 == 1:
                fH2Oo = G_1 * P  # The fugacity of H2O
                # return fH2Oo
        return fH2Oo


class fugacity_HB_co2(FugacityModel):
    """
    Implementation of the Holloway and Blank (1994) Modified Redlich Kwong EoS for CO2.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', 500.0, calibration_checks.crf_GreaterThan, 'oC',
                'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_GreaterThan_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])
        self.HBmodel = fugacity_HollowayBlank()

    def fugacity(self, pressure, temperature, X_fluid=1.0, **kwargs):
        pure_f = self.HBmodel.fugacity(pressure=pressure, temperature=temperature, species='CO2')
        return pure_f * X_fluid


class fugacity_HB_h2o(FugacityModel):
    """
    Implementation of the Holloway and Blank (1994) Modified Redlich Kwong EoS for H2O.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', 500.0, calibration_checks.crf_GreaterThan, 'oC',
                'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_GreaterThan_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])
        self.HBmodel = fugacity_HollowayBlank()

    def fugacity(self, pressure, temperature, X_fluid=1.0, **kwargs):
        pure_f = self.HBmodel.fugacity(pressure=pressure, temperature=temperature, species='H2O')
        return pure_f * X_fluid


class fugacity_HollowayBlank(FugacityModel):
    """
    Implementation of the Modified Redlich Kwong presented in Holloway and Blank (1994) Reviews
    in Mineralogy and Geochemistry vol. 30. Originally written in Quickbasic. CO2 calculations
    translated to Matlab by Chelsea Allison and translated to python by K. Iacovino for VESIcal.
    H2O calculations translated to VisualBasic by Gordon M. Moore and translated to python by
    K. Iacovino for VESIcal.

    """

    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar',
                'MRK EOS (Holloway and Blank, 1994)',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', 500, calibration_checks.crf_GreaterThan, 'oC',
                'MRK EOS (Holloway and Blank, 1994)',
                fail_msg=calibration_checks.crmsg_GreaterThan_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])

    def REDKW(self, BP, A2B):
        """
        The RK routine. A routine to calculate compressibility factor and fugacity coefficient
        with the Redlich-Kwong equation following Edmister (1968). This solution for supercritical
        fluid.

        Parameters
        ----------
        BP: float
            B parameter sum from RKCALC

        A2B: float
            A parameter sum from RKCALC

        Returns
        -------
        float
            XLNFP (fugacity coefficient?)
        """
        if A2B < 1*10**(-10):
            A2B = 0.001

        # Define constants
        TH = 0.333333
        RR = -A2B*BP**2
        QQ = BP*(A2B-BP-1)
        XN = QQ*TH+RR-0.074074
        XM = QQ-TH
        XNN = XN*XN*0.25
        XMM = XM**3 / 27.0
        ARG = XNN+XMM

        if ARG > 0:
            X = np.sqrt(ARG)
            F = 1
            XN2 = -XN*0.5
            iXMM = XN2+X

            if iXMM < 0:
                F = -1
            XMM = F*((F*iXMM)**TH)
            F = 1
            iXNN = XN2 - X

            if iXNN < 0:
                F = -1
            XNN = F*((F*iXNN)**TH)
            Z = XMM+XNN+TH
            ZBP = Z-BP
            if ZBP < 0.000001:
                ZBP = 0.000001

            BPZ = 1+BP/Z
            FP = Z-1-np.log(ZBP)-A2B*np.log(BPZ)

            if FP < -37 or FP > 37:
                FP = 0.000001

        elif ARG < 0:
            COSPHI = np.sqrt(-XNN/XMM)
            if XN > 0:
                COSPHI = -COSPHI

            TANPHI = np.sqrt(1-COSPHI**2)/COSPHI
            PHI = np.arctan(TANPHI)*TH
            FAC = 2*np.sqrt(-XM*TH)

            # sort for largest root
            R1 = np.cos(PHI)
            R2 = np.cos(PHI+2.0944)
            R3 = np.cos(PHI+4.18879)
            RH = R2

            if R1 > R2:
                RH = R1
            if R3 > RH:
                RH = R3

            Z = RH*FAC+TH
            ZBP = Z-BP
            if ZBP < 0.000001:
                ZBP = 0.000001
            BPZ = 1+BP/Z
            FP = Z-1-np.log(ZBP)-A2B*np.log(BPZ)
            if FP < -37 or FP > 37:
                FP = 0.000001
        else:
            FP = 1
            Z = 1
        XLNFP = FP

        return XLNFP

    def Saxena(self, TK, pb):
        """
        High pressure corresponding states routines from Saxena and Fei (1987) GCA
        vol. 51, 783-791.

        Parameters
        ----------
        TK: float
            Temperature in K.

        pb: float
            Pressure in bars.

        Returns
        -------
        float
            XLNF, Natural log of the ratio F(P)/F(4000 bar)
        """

        # Define integration limit
        PO = 4000

        # Critical temperatures and pressures for CO2
        TR = TK/304.2
        PC = 73.9

        # Virial coeficients
        A = 2.0614-2.2351/TR**2 - 0.39411*np.log(TR)
        B = 0.055125/TR + 0.039344/TR**2
        C = -1.8935*10**(-6)/TR - 1.1092*10**(-5)/TR**2 - 2.1892*10**(-5)/TR**3
        D = 5.0527*10**(-11)/TR - 6.3033*10**(-21)/TR**3

        # integrate from PO (4000 bars) to P to calculate ln fugacity
        LNF = A*np.log(pb/PO)+(B/PC)*(pb-PO)+(C/(2*PC**2))*(pb**2-PO**2)
        LNF = LNF+(D/(3*PC**3))*(pb**3-PO**3)
        XLNF = LNF

        return XLNF

    def RKCALC(self, temperature, pressure, species):
        """
        Calculation of pure gas MRK properties following Holloway 1981, 1987

        Parameters
        ----------
        temperature: float
            Temperature in degrees K.

        pressure: float
            Pressure in atmospheres.

        Returns
        -------
        float
            Natural log of the fugacity of a pure gas.
        """
        # Define constants
        R = 82.05736
        pb = 1.013*pressure
        PBLN = np.log(pb)
        TCEL = temperature-273.15
        RXT = R*temperature
        RT = R*temperature**1.5 * 10**(-6)

        if species == 'CO2':
            # Calculate T-dependent MRK A parameter CO2
            ACO2M = 73.03 - 0.0714*TCEL + 2.157*10**(-5)*TCEL**2

            # Define MRK B parameter for CO2
            BSUM = 29.7

            ASUM = ACO2M / (BSUM*RT)

        elif species == 'H2O':
            # Calculate T-dependent MRK A parameter H2O
            AH2OM = 115.98 - np.double(0.0016295)*temperature - 1.4984*10**(-5)*temperature**2

            # Define MRK B parameter for H2O
            BSUM = 14.5

            ASUM = AH2OM / (BSUM*RT)

        BSUM = pressure*BSUM/RXT
        XLNFP = self.REDKW(BSUM, ASUM)

        # Convert to ln(fugacity)
        PUREG = XLNFP + PBLN
        return PUREG

    def fugacity(self, pressure, temperature, species, **kwargs):
        """
        Calculates fugacity.

        Parameters
        ----------
        temperature: float
            Temperature in degrees C.

        pressure: float
            Pressure in bars.

        species: str
            Choose which species to calculate. Options are 'H2O' and 'CO2'.

        Returns
        -------
        float
            Fugacity coefficient for passed species
        """

        # convert temp and press to atmospheres and Kelvin
        pressureAtmo = pressure/1.013
        temperatureK = temperature + 273.15
        PO = 4000/1.013

        # Use the MRK below 4,000 bars, Saxena above 4,000 bars
        if pressure > 4000 and species == 'CO2':
            iPUREG = self.RKCALC(temperatureK, PO, species)
            XLNF = self.Saxena(temperatureK, pressure)
            PUREG = iPUREG + XLNF
        else:
            PUREG = self.RKCALC(temperatureK, pressureAtmo, species)

        # Convert from ln(fugacity) to fugacity
        stdf = np.exp(PUREG)
        return stdf


class fugacity_RK_co2(FugacityModel):
    """
    Implementation of the Redlich Kwong EoS for CO2.
    Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30
    October 2003.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', [500], calibration_checks.crf_GreaterThan, 'oC',
                'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_GreaterThan_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])

        self.RKmodel = fugacity_RedlichKwong()

    def fugacity(self, pressure, temperature, X_fluid, **kwargs):
        return self.RKmodel.fugacity(pressure, temperature, X_fluid, 'CO2')


class fugacity_RK_h2o(FugacityModel):
    """
    Implementation of the Redlich Kwong EoS for H2O.
    Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30
    October 2003.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', 500, calibration_checks.crf_GreaterThan, 'oC', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_GreaterThan_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])
        self.RKmodel = fugacity_RedlichKwong()

    def fugacity(self, pressure, temperature, X_fluid, **kwargs):
        return self.RKmodel.fugacity(pressure, temperature, X_fluid, 'H2O')


class fugacity_RedlichKwong(FugacityModel):
    """
    Implementation of the Redlich Kwong EoS
    Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30
    October 2003.
    """
    def __init__(self):
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [1, 1e5], calibration_checks.crf_Between, 'bar', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', 500, calibration_checks.crf_GreaterThan, 'oC', 'Redlich Kwong EOS',
                fail_msg=calibration_checks.crmsg_GreaterThan_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])

    def gamma(self, pressure, temperature, species):
        """
        Calculates fugacity coefficients.

        Parameters
        ----------
        temperature: fload
            Temperature in degrees C.

        pressure: float
            Pressure in bars.

        species: str
            Choose which species to calculate. Options are 'H2O' and 'CO2'.

        Returns
        -------
        float
            Fugacity coefficient for passed species.
        """

        temperatureK = temperature + 273.15
        R = 8.3145

        critical_params = {'CO2': {"cT":   304.15,
                                   "cP":   73.8659,
                                   "o":    0.225
                                   },
                           'H2O': {"cT":   647.25,
                                   "cP":   221.1925,
                                   "o":    0.334
                                   }
                           }

        # Calculate a and b parameters (depend only on critical parameters)...
        a = (0.42748 * R**2.0 * critical_params[species]["cT"]**(2.5) /
             (critical_params[species]["cP"] * 10.0**5))
        b = (0.08664 * R * critical_params[species]["cT"] /
             (critical_params[species]["cP"] * 10.0**5))

        # Calculate coefficients in the cubic equation of state...
        # coeffs: (C0, C1, C2, A, B)
        A = a * pressure * 10.0**5 / (np.sqrt(temperatureK) * (R * temperatureK)**2.0)
        B = b * pressure * 10.0**5 / (R * temperatureK)
        C2 = -1.0
        C1 = A - B - B * B
        C0 = -A * B

        # Solve the cubic equation for Z0 - Z2, D...
        Q1 = C2 * C1 / 6.0 - C0 / 2.0 - C2**3.0 / 27.0
        P1 = C2**2.0 / 9.0 - C1 / 3.0
        D = Q1**2.0 - P1**3.0

        if D >= 0:
            kOneThird = 1.0 / 3.0

            absQ1PSqrtD = np.fabs(Q1 + np.sqrt(D))
            temp1 = absQ1PSqrtD**kOneThird
            temp1 *= (Q1 + np.sqrt(D)) / absQ1PSqrtD

            absQ1MSqrtD = np.fabs(Q1 - np.sqrt(D))
            temp2 = absQ1MSqrtD**kOneThird
            temp2 *= (Q1 - np.sqrt(D)) / absQ1MSqrtD

            Z0 = temp1 + temp2 - C2 / 3.0
        else:
            temp1 = Q1**2.0 / (P1**3.0)
            temp2 = np.sqrt(1.0 - temp1) / np.sqrt(temp1)
            temp2 *= Q1 / np.fabs(Q1)

            gamma = np.arctan(temp2)

            if gamma < 0:
                gamma = gamma + np.pi

            Z0 = 2.0 * np.sqrt(P1) * np.cos(gamma/3.0) - C2 / 3.0
            Z1 = 2.0 * np.sqrt(P1) * np.cos((gamma + 2.0 * np.pi) / 3.0) - C2/3.0
            Z2 = 2.0 * np.sqrt(P1) * np.cos((gamma + 4.0 * np.pi) / 3.0) - C2/3.0

            if Z0 < Z1:
                temp0 = Z0
                Z0 = Z1
                Z1 = temp0

            if Z1 < Z2:
                temp0 = Z1
                Z1 = Z2
                Z2 = temp0

            if Z0 < Z1:
                temp0 = Z0
                Z0 = Z1
                Z1 = temp0

        # Calculate Departure Functions
        gamma = np.exp(Z0 - 1.0 - np.log(Z0-B) - A * np.log(1.0+B/Z0)/B)

        return gamma

    def fugacity(self, pressure, temperature, X_fluid=1.0, species='H2O', **kwargs):
        """
        Calculates the fugacity of H2O in a mixed H2O-CO2 fluid using the universal relationships:
        P_i = f_i/gamma_i = (fpure_i * Xfluid_i) / gamma_i
        See Iacovino (2015) EPSL for further explanation.
        """

        gammaH2O = self.gamma(pressure, temperature, 'H2O')
        gammaCO2 = self.gamma(pressure, temperature, 'CO2')

        fugacityH2Opure = pressure * gammaH2O
        fugacityCO2pure = pressure * gammaCO2

        if species == 'H2O':
            return fugacityH2Opure * X_fluid
        elif species == 'CO2':
            return fugacityCO2pure * X_fluid
        else:
            raise core.InputError("Species must be H2O or CO2.")
