from copy import copy
import numpy as np


class CalibrationRange(object):
    """ The CalibrationRange object allows the range of allowable parameters
    to be specified and used in checking and reporting of the results.
    """
    def __init__(self, parameter_name, value, checkfunction=None, units='',
                 model_name='', fail_msg='', fail_dict={}, pass_msg='',
                 pass_dict={}, description_msg='', description_dict={}):
        self.parameter_name = parameter_name
        self.value = value
        self.checkfunction = checkfunction
        self.units = units
        self.model_name = model_name
        self.fail_msg = (copy(fail_msg), copy(fail_dict))
        self.pass_msg = (copy(pass_msg), copy(pass_dict))
        self.description_msg = (copy(description_msg), copy(description_dict))

    def check(self, parameters):
        """Method for checking whether parameters satisfy the calibration
        range.
        """
        if self.parameter_name in parameters:
            return self.checkfunction(self.value,
                                      parameters[self.parameter_name])
        else:
            return None

    def string(self, parameters, report_nonexistance=True):
        """Returns a string statement of the calibration check"""
        if parameters is None:
            msgdict = self.description_msg[1]
            if isinstance(self.value, float) or isinstance(self.value, int):
                msgdict['calib_val'] = self.value
            elif (isinstance(self.value, list) or isinstance(self.value, tuple) or
                  isinstance(self.value, np.ndarray)):
                for i in range(len(self.value)):
                    msgdict['calib_val'+str(i)] = self.value[i]
            if 'param_name' not in msgdict:
                msgdict['param_name'] = self.parameter_name
            if 'units' not in msgdict:
                msgdict['units'] = self.units
            if 'model_name' not in msgdict:
                msgdict['model_name'] = self.model_name
            return self.description_msg[0].format(**msgdict)
        else:
            check = self.check(parameters)
            if check:
                msgdict = self.pass_msg[1]
                msgdict['param_val'] = parameters[self.parameter_name]
                if isinstance(self.value, float) or isinstance(self.value, int):
                    msgdict['calib_val'] = self.value
                elif (isinstance(self.value, list) or isinstance(self.value, tuple) or
                      isinstance(self.value, np.ndarray)):
                    for i in range(len(self.value)):
                        msgdict['calib_val'+str(i)] = self.value[i]
                if 'param_name' not in msgdict:
                    msgdict['param_name'] = self.parameter_name
                if 'units' not in msgdict:
                    msgdict['units'] = self.units
                if 'model_name' not in msgdict:
                    msgdict['model_name'] = self.model_name
                return self.pass_msg[0].format(**msgdict)
            elif check is False:
                msgdict = self.fail_msg[1]
                msgdict['param_val'] = parameters[self.parameter_name]
                if isinstance(self.value, float) or isinstance(self.value, int):
                    msgdict['calib_val'] = self.value
                elif (isinstance(self.value, list) or isinstance(self.value, tuple) or
                      isinstance(self.value, np.ndarray)):
                    for i in range(len(self.value)):
                        msgdict['calib_val'+str(i)] = self.value[i]
                if 'param_name' not in msgdict:
                    msgdict['param_name'] = self.parameter_name
                if 'units' not in msgdict:
                    msgdict['units'] = self.units
                if 'model_name' not in msgdict:
                    msgdict['model_name'] = self.model_name
                return self.fail_msg[0].format(**msgdict)
            else:
                if report_nonexistance:
                    return "A value for {} was not provided.".format(
                                                           self.parameter_name)
                else:
                    return ''

# ------------- DEFAULT CALIBRATIONRANGE OBJECTS --------------- #


def crf_EqualTo(calibval, paramval):
    return calibval == paramval


crmsg_EqualTo_pass = ("The {param_name} ({param_val:.1f} {units}) is equal to "
                      "{calib_val:.1f} {units} as required by the calibration "
                      "range of the {model_name} model. ")
crmsg_EqualTo_fail = ("{param_name} ({param_val:.1f} {units}) is outside the "
                      "calibration range of the {model_name} model "
                      "({calib_val:.1f} {units}). ")
crmsg_EqualTo_description = ("The {model_name} model is calibrated for "
                             "{param_name} equal to {calib_val:.1f} {units}. ")


def crf_MixedFluidWarning(calibval, paramval):
    if isinstance(paramval, float) or isinstance(paramval, int) or len(paramval) == 1:
        return paramval == 0 or paramval == 1
    elif len(paramval) == 2:
        return list(paramval) == [0, 1] or list(paramval) == [1, 0]
    else:
        return True


crmsg_MixedFluidWarning_pass = "The pure fluid model is being used. "
crmsg_MixedFluidWarning_fail = ("The {model_name} model should be used to "
                                "model mixed fluids with caution. The mixed "
                                "fluid model does not recreate the "
                                "experimental observations. ")
crmsg_MixedFluidWarning_description = ("The {model_name} shoule be used to "
                                       "model mixed fluids with caution. ")


def crf_GreaterThan(calibval, paramval):
    return paramval >= calibval


crmsg_GreaterThan_pass = ("The {param_name} ({param_val:.1f} {units}) is "
                          "greater than {calib_val:.1f} {units} as required "
                          "by the calibration range of the {model_name} "
                          "model. ")
crmsg_GreaterThan_fail = ("The {param_name} is outside the calibration range "
                          "of {model_name} ({param_val:.1f}<{calib_val:.1f} "
                          "{units}. ")
crmsg_GreaterThan_description = ("The {model_name} model is calibrated for "
                                 "{param_name} greater than {calib_val:.1f} "
                                 "{units}. ")


def crf_LessThan(calibval, paramval):
    return paramval <= calibval


crmsg_LessThan_pass = ("The {param_name} ({param_val:.1f} {units}) is less "
                       "than {calib_val:.1f} {units} as required by the "
                       "calibration range of the {model_name} model. ")
crmsg_LessThan_fail = ("The {param_name} is outside the calibration range "
                       "of {model_name} ({param_val:.1f}>{calib_val:.1f} "
                       "{units}. ")
crmsg_LessThan_description = ("The {model_name} model is calibrated for "
                              "{param_name} less than {calib_val:.1f} "
                              "{units}. ")


def crf_Between(calibval, paramval):
    if isinstance(paramval, np.ndarray):
        return paramval.any() >= calibval[0] and paramval.any() <= calibval[1]
    else:
        return paramval >= calibval[0] and paramval <= calibval[1]


crmsg_Between_pass = ("The {param_name} ({param_val:.1f} {units}) is between "
                      "{calib_val0:.1f} and {calib_val1:.1f} {units} as "
                      "required by the calibration range of the {model_name} "
                      "model. ")
crmsg_Between_fail = ("{param_name} ({param_val:.1f} {units}) is outside the "
                      "calibration range of the {model_name} model "
                      "({calib_val0:.1f}-{calib_val1:.1f} {units}). ")
crmsg_Between_description = ("The {model_name} model is calibrated for "
                             "{param_name} between {calib_val0:.1f} and "
                             "{calib_val1:.1f} {units}. ")
# Different wording for compositional checks (loosing "the")
crmsg_BC_pass = ("{param_name} ({param_val:.1f} {units}) is between "
                 "{calib_val0:.1f} and {calib_val1:.1f} {units} as required "
                 "by the calibration range of the {model_name} model. ")
crmsg_BC_fail = ("{param_name} ({param_val:.1f} {units}) is outside the "
                 "calibration range ({calib_val0:.1f}-{calib_val1:.1f} "
                 "{units}). ")
