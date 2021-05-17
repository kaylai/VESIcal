from abc import abstractmethod


class activity_model(object):
    """ The activity model object is for implementing activity models
    for volatile species in melts. It contains all the methods required to
    evaluate the activity.
    """
    def __init__(self):
        self.set_calibration_ranges([])

    def set_calibration_ranges(self, calibration_ranges):
        self.calibration_ranges = calibration_ranges

    @abstractmethod
    def activity(self, X, **kwargs):
        """
        """

    # @abstractmethod
    def check_calibration_range(self, parameters, report_nonexistance=True):
        s = ''
        for cr in self.calibration_ranges:
            if cr.check(parameters) is False:
                s += cr.string(parameters, report_nonexistance)
        return s


# --------------- ACTVITY MODELS ------------------------------- #

class activity_idealsolution(activity_model):
    """ Implements an ideal solution activity model, i.e. it
    will always return the mole fraction.
    """

    def activity(self, X):
        """ The activity of the component in an ideal solution, i.e., it
        will return the mole fraction.

        Parameters
        ----------
        X   float
            The mole fraction of the species in the solution.

        Returns
        -------
        float
            The activity of the species in the solution, i.e.,
            the mole fraction.
        """
        return X
