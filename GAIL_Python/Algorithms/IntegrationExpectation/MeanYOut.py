#!/usr/bin/env python
"""Definition of MeanYOut object"""

__author__ = ["Anil Simon", "Divya Vasireddy"]

import numbers
from operator import attrgetter
import inspect

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.helper_functions import set_named_args
except:
    from helper_functions import set_named_args

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.meanYParam import MeanYParam
except:
    from meanYParam import MeanYParam




class MeanYOut(MeanYParam):
    """
    MeanYOut is a class containing the parameters related to the
    outputs from the algorithms that find the mean of a random variable.
    This class includes the time and sample size required for the
    computation
    """
    stddev = property(attrgetter('_stddev'))

    @stddev.setter
    def stddev(self, stddev_value):  # {'numeric'}, {'nonnegative'})
        msg = "stddev should be positive number"
        if stddev_value is None:
            raise Exception("stddev is empty. " + msg)
        if not isinstance(stddev_value, numbers.Number):
            raise Exception("stddev value should be a number. " + msg)
        else:
            if stddev_value < 0:
                raise Exception("stddev is negative. " + msg)

        self._stddev = stddev_value

    sol = property(attrgetter('_sol'))

    @sol.setter
    def sol(self, sol_value):
        msg = "sol should be positive number"
        if sol_value is None:
            raise Exception("sol is empty. " + msg)
        if not isinstance(sol_value, numbers.Number):
            raise Exception("sol value should be a number. " + msg)
        self._sol = sol_value

    nSample = property(attrgetter('_nSample'))

    @nSample.setter
    def nSample(self, nSample_value):
        msg = "nSample should be positive number"
        if nSample_value is None:
            raise Exception("nSample is empty. " + msg)
        if not isinstance(nSample_value, numbers.Number):
            raise Exception("nSample value should be a number. " + msg)
        else:
            if nSample_value < 0:
                raise Exception("nSample is negative. " + msg)

        self._nSample = nSample_value

    errBd = property(attrgetter('_errBd'))

    @errBd.setter
    def errBd(self, errBd_value):
        msg = "errBd should be positive number"
        if errBd_value is None:
            raise Exception("errBd is empty. " + msg)
        if not isinstance(errBd_value, numbers.Number):
            raise Exception("errBd value should be a number. " + msg)
        else:
            if errBd_value < 0:
                raise Exception("errBd is negative. " + msg)

        self._errBd = errBd_value

    def __init__(self, Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, **kwargs):
        super().__init__(Y, absTol, relTol, alpha, nSig, **kwargs)
        meanYOutParams = {'stddev', 'errBd', 'nSample', 'sol'}
        for k, v in kwargs.items():
            if k in meanYOutParams:
                setattr(self, k, v)
        self.time = 0 # time in seconds

        #
        # # set mis-spelled input arguments
        # params_names = ['Y', 'absTol', 'relTol', 'alpha', 'nSig']
        # set_named_args(self, kwargs, params_names)


    def __str__(self):
        try:
            function_source = inspect.getsource(self.Y)
        except:
            function_source = str(self.Y)

        return 'meanYOut with properties:\n\tY\t: {}\n\tabsTol\t: {}\n\trelTol\t: {}\n\talpha\t: {}\n\tnSig\t: {}' \
            .format(function_source, self.absTol, self.relTol, self.alpha, self.nSig)