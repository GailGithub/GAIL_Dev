import numbers
from operator import attrgetter

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
        self._sol = sol_value

    nSample = property(attrgetter('_nSample'))

    @nSample.setter
    def nSample(self, nSample_value):
        self._nSample = nSample_value

    errBd = property(attrgetter('_errBd'))

    @errBd.setter
    def errBd(self, errBd_value):
        self._errBd = errBd_value

    def __init__(self, Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, **kwargs):
        super().__init__(Y, absTol, relTol, alpha, nSig, **kwargs)
        meanYOutParams = {'stddev', 'errBd', 'nSample', 'sol'}
        for k, v in kwargs.items():
            if k in meanYOutParams:
                setattr(self, k, v)
