from operator import attrgetter
import numbers
import numpy as np


class ErrorParam(object):
    """
    This class contains the error tolerances, solution function, and
    related parameters determining the error criterion.
    """

    def __init__(self, abs_tol=1e-2, rel_tol=0):
        self.abs_tol = abs_tol
        self.rel_tol = rel_tol
        self.sol_fun = lambda mu: mu
        self.solbd_fun = lambda muhat, errbd: (muhat - errbd, muhat + errbd)

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

    abs_tol = property(attrgetter('_abs_tol'))

    @abs_tol.setter
    def abs_tol(self, abs_tol_value):  # {'numeric'}, {'scalar','nonnegative', '<', 1})
        msg = "abs_tol should be positive number"
        if not abs_tol_value:
            raise Exception("value cannot be empty. " + msg)
        if not isinstance(abs_tol_value, numbers.Number):
            raise Exception("abs_tol value should be a number. " + msg)
        else:
            if abs_tol_value < 0:
                raise Exception("abs_tol value cannot be negative. " + msg)
        self._abs_tol = abs_tol_value

    rel_tol = property(attrgetter('_rel_tol'))

    @rel_tol.setter
    def rel_tol(self, rel_tol_value):  # {'numeric'}, {'scalar','nonnegative'})
        msg = "rel_tol should be a positive number"
        if not rel_tol_value:
            raise Exception("value cannot be empty. " + msg)
        if not isinstance(rel_tol_value, numbers.Number):
            raise Exception("rel_tol value should be a number. " + msg)
        else:
            if rel_tol_value < 0:
                raise Exception("rel_tol value cannot be negative. " + msg)
        self._rel_tol = rel_tol_value
