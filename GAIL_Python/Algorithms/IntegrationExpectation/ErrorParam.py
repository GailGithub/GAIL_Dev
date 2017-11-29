#!/usr/bin/env python
"""Definition of ErrorParam object"""

__author__ = ["Anil Simon", "Divya Vasireddy"]

from operator import attrgetter
import numbers

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.helper_functions import identity_function, ci_calculator, \
        set_named_args, default_random_generator
except ModuleNotFoundError:
    from helper_functions import identity_function, ci_calculator, set_named_args, default_random_generator


class ErrorParam(object):
    """ Definition of ErrorParam Object
    This class contains the error tolerances, solution function, and
    related parameters determining the error criterion.
    """
    solFun = property(attrgetter('_solFun'))

    @solFun.setter
    def solFun(self, sol_fun_func):
        """To judge if input is a function handle or not"""
        if not callable(sol_fun_func):
            raise Exception("solFun has to be a function")
        self._solFun = sol_fun_func

    solBdFun = property(attrgetter('_solBdFun'))

    @solBdFun.setter
    def solBdFun(self, solbd_fun_func):
        """To judge if input is a function handle or not"""
        if not callable(solbd_fun_func):
            raise Exception("solBdFun has to be a function")
        self._solBdFun = solbd_fun_func

    absTol = property(attrgetter('_absTol'))

    @absTol.setter
    def absTol(self, absTol_value):  # {'numeric'}, {'scalar','nonnegative', '<', 1})
        msg = "absTol should be positive number"
        if absTol_value is None:
            raise Exception("value cannot be empty. " + msg)
        if not isinstance(absTol_value, numbers.Number):
            raise Exception("absTol value should be a number. " + msg)
        else:
            if absTol_value < 0:
                raise Exception("absTol value cannot be negative. " + msg)
        self._absTol = absTol_value

    relTol = property(attrgetter('_relTol'))

    @relTol.setter
    def relTol(self, relTol_value):  # {'numeric'}, {'scalar','nonnegative'})
        msg = "relTol should be a positive number"
        if relTol_value is None:
            raise Exception("value cannot be empty. " + msg)
        if not isinstance(relTol_value, numbers.Number):
            raise Exception("relTol value should be a number. " + msg)
        else:
            if relTol_value < 0:
                raise Exception("relTol value cannot be negative. " + msg)
        self._relTol = relTol_value

    def __init__(self, absTol=1e-2, relTol=0, solFun=None, solBdFun=None, **kwargs):
        self.absTol = absTol
        self.relTol = relTol
        if solFun is None:
            self.solFun = identity_function
        else:
            self.solFun = solFun
        if solBdFun is None:
            self.solBdFun = ci_calculator
        else:
            self.solBdFun = solBdFun

        # set mis-spelled input arguments
        params_names = ['absTol', 'relTol', 'solFun', 'solBdFun']
        set_named_args(self, kwargs, params_names)

