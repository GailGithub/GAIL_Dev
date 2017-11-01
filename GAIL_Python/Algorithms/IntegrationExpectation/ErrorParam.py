from operator import attrgetter
import numbers
import numpy as np


class ErrorParam(object):
    """
    This class contains the error tolerances, solution function, and
    related parameters determining the error criterion.
    """
    sol_fun = property(attrgetter('_sol_fun'))

    @sol_fun.setter
    def sol_fun(self, sol_fun_func):
        # ISFCN To judge if input is a function handle or not
        if not callable(sol_fun_func):
            raise Exception("sol_fun has to be a function")
        self._sol_fun = sol_fun_func

    solbd_fun = property(attrgetter('_solbd_fun'))

    @solbd_fun.setter
    def solbd_fun(self, solbd_fun_func):
        # ISFCN To judge if input is a function handle or not
        if not callable(solbd_fun_func):
            raise Exception("solbd_fun has to be a function")
        self._solbd_fun = solbd_fun_func

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

    def __init__(self, absTol=1e-2, relTol=0, solFun=None, solbdFun=None, **kwargs):
        self.absTol = absTol
        self.relTol = relTol
        if solFun is None:
            self.solFun = lambda mu: mu
        else:
            self.solFun = solFun
        if solbdFun is None:
            self.solbdFun = lambda muhat, errbd: (muhat - errbd, muhat + errbd)
        else:
            self.solbdFun = solbdFun
        meanYparams = {'sol_fun': 'solFun', 'solbd_fun': 'solbdFun', 'solfun': 'solFun'}
        for k, v in kwargs.items():
            if k in meanYparams:
                setattr(self, meanYparams[k], v)
