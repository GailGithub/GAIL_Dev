import inspect
import numbers
from operator import attrgetter

import numpy as np
try:
    from GAIL_Python.Algorithms.IntegrationExpectation.CubMeanParam import CubMeanParam
except:
    from CubMeanParam import CubMeanParam
try:
    from GAIL_Python.Algorithms.IntegrationExpectation.ErrorParam import ErrorParam
except:
    from ErrorParam import ErrorParam


def default_random_generator(n):
    """
    X_i^2 where X_i comes from a U[0,1]
    :param n: n is length of random numbers needed
    :return: vector of random numbers
    """
    return np.random.uniform(size=n) ** 2


class MeanYParam(object):
    """
    MeanYParam is a class containing the parameters related to
    algorithms that find the mean of a random variable
    This class contains the random number generator, the uncertainty
    """

    @property
    def Yout(self):
        if self._Yout is None:
            self._Yout = len(self.y(1))
        return self._Yout

    @Yout.setter
    def Yout(self, Yout_value):
        self._Yout = Yout_value

    Y = property(attrgetter('_Y'))

    @Y.setter
    def Y(self, Y_func):
        # ISFCN To judge if input is a function handle or not
        if not callable(Y_func):
            raise Exception("y has to be a function")
        self._Y = Y_func

    alpha = property(attrgetter('_alpha'))

    @alpha.setter
    def alpha(self, alpha_value):  # {'numeric'}, {'scalar','nonnegative', '<', 1})
        msg = "alpha should be positive and less than 1"
        if alpha_value is None:
            raise Exception("alpha is empty. " + msg)
        if not isinstance(alpha_value, numbers.Number):
            raise Exception("alpha value should be a number. " + msg)
        else:
            if alpha_value < 0:
                raise Exception("alpha is negative. " + msg)
            elif alpha_value > 1:
                raise Exception("alpha is greater than 1. " + msg)
        self._alpha = alpha_value

    nSig = property(attrgetter('_nSig'))

    @nSig.setter
    def nSig(self, nSig_value):  # {'numeric'}, {'scalar','nonnegative'})
        msg = "nSig should be a positive number"
        if nSig_value is None:
            raise Exception("nSig is empty. " + msg)
        if not isinstance(nSig_value, numbers.Number):
            raise Exception("nSig is not a number. " + msg)
        else:
            if nSig_value < 0:
                raise Exception("nSig is negative. " + msg)
        self._nSig = nSig_value

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

    def __init__(self, Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, **kwargs):
        if isinstance(Y, MeanYParam):
            # copy the MeanYParam object
            self.make_copy(Y)
        else:
            if Y is None:
                self.Y = default_random_generator
            else:
                self.Y = Y  # random number generator

            self.absTol = absTol
            self.relTol = relTol

            self.alpha = alpha  # uncertainty
            self.nSig = nSig  # sample size to estimate variance

            self.err = ErrorParam(absTol=absTol, relTol=relTol, **kwargs)  # an errorParam object
            self.CM = CubMeanParam(**kwargs)  # a cubMeanParam object
            try:
                self.nY = kwargs['nY']
            except:
                self.nY = 1  # number of Y for each mean

    def make_copy(self, Y):
        self.Y = Y.Y
        self.absTol = Y.absTol
        self.relTol = Y.relTol
        self.alpha = Y.alpha  # uncertainty
        self.nSig = Y.nSig  # sample size to estimate variance
        self.err = ErrorParam(absTol=Y.err.absTol, relTol=Y.err.relTol)
        self.CM = CubMeanParam(inflate=Y.CM.inflate, inflateFun=Y.CM.inflateFun, nInit=Y.CM.nInit, nMax=Y.CM.nMax,
                               nMu=Y.CM.nMu, trueMuCV=Y.CM.trueMuCV)
        self.nY = Y.nY

    def __str__(self):

        try:
            function_source = inspect.getsource(self.Y)
        except:
            function_source = str(self.Y)

        return 'meanYParam with properties:\n\tY\t: {}\n\tabsTol\t: {}\n\trelTol\t: {}\n\talpha\t: {}\n\tnSig\t: {}' \
            .format(function_source, self.absTol, self.relTol, self.alpha, self.nSig)