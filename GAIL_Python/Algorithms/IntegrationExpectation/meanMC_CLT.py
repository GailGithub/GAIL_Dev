# MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
import time
import numpy as np
from scipy.special import erfcinv
from copy import copy, deepcopy

class CubMeanParam(object):
    """
    CubMeanParam is a class containing the parameters related to
    algorithms that act on functions of x
    This class contains the function, its domain, etc.
    """
    pass

class ErrorParam(object):
    """
    ErrorParam is a class containing the parameters related to the
    error tolerance
    This class contains the error tolerances, solution function, and
    related parameters determining the error criterion
    """
    pass



class MeanYParam(object):
    """
    MeanYParam is a class containing the parameters related to
    algorithms that find the mean of a random variable
    This class contains the random number generator, the uncertainty
    """
    def __init__(self, y=None, abs_tol=1e-2, rel_tol=0, alpha=0.01, n_sig=1000):
        self.y = y  # random number generator

        self.abs_tol = abs_tol
        self.rel_tol = rel_tol
        self.alpha = alpha  # uncertainty
        self.nSig = n_sig  # sample size to estimate variance

        def default_random_generator(n=self.nSig):
            # X_i^2 where X_i comes from a U[0,1]
            return np.random.uniform(size=n)**2

        if self.y is None:
            self.y = default_random_generator

        self.err = {}  # an errorParam object #TODO
        self.CM = {}  # a cubMeanParam object #TODO
        self.nY = 1  # number of Y for each mean



    def check_parameters(self):
        pass


    # def __copy__(self):
    #     cls =self.__class__
    #     result = cls.__new__(cls)
    #     result.__dict__.update(self.__dict__)
    #     return result
    #
    # def __deepcopy__(self, memodict={}):
    #     cls = self.__class__
    #     result = cls.__new__(cls)
    #     memodict[id(self)] = result
    #     for k, v in self.__dict__
    #     result.__dict__.update(self.__dict__)
    #     return result


class MeanYOut(MeanYParam):
    """
    MeanYOut is a class containing the parameters related to the
    outputs from the algorithms that find the mean of a random variable.
    This class includes the time and sample size required for the
    computation
    """
    def __init__(self, y=None, abs_tol=1e-2, rel_tol=0, alpha=0.01, n_sig=1000, inflate=1.2):
        super().__init__(y, abs_tol, rel_tol, alpha, n_sig, inflate)
        self.inflate = inflate
        self.nSig = n_sig
        self.timekeeper = {}
        self.out = {}


class meanMC_CLT(MeanYOut):
    def __init__(self, y=None, abs_tol=1e-2, rel_tol=0, alpha=0.01, n_sig=1000, inflate=1.2):
        super().__init__(y, abs_tol, rel_tol, alpha, n_sig, inflate)
        self.out = {}


    def generate_y_values(self):
        self.yy = np.empty(self.n_sig)
        for i, el in enumerate(self.y): self.yy[i] = el

def check_meanMC_CLT_params(**kwargs):
    checked_params = kwargs

    #TODO meanYOut and meanYParam functions
    # if not (0 <= relTol <= 1):
    #     logging.error("Relative error tolerance should be in the range [0,1]")
    return checked_params


def stdnorminv(p):
    """
    This function is the inverse function of CDF of standard normal distribution
    :param p:
    :return:
    """
    erfinv_val = erfcinv(np.multiply(2, p))
    return np.multiply(-1 * np.sqrt(2), erfinv_val)


def meanMC_CLT(Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, inflate=1.2):
    start_time = time.time()

    out = check_meanMC_CLT_params(**{'Y': Y, 'absTol': absTol, 'relTol': relTol,
                                     'alpha': alpha, 'nSig': nSig, 'inflate': inflate})

    Yrand = out['Y']
    p = 0  # out['CM']['nCV']
    q = 1  # out['nY']
    val = (y for y in Yrand(out['nSig']))

    if p == 0 and q == 1:
        YY = list(val)
    else:
        # TODO replace with control variate procedures
        YY = list(val)

    out['stddev'] = np.std(YY)
    sig0up = np.multiply(out['inflate'], out['stddev'])  # TODO change inflate factor out.CM.inflate
    hmu0 = np.mean(YY)

    nmu = int(max(1, np.power(np.ceil(
        np.multiply(
            np.multiply(-1, stdnorminv(out['alpha'] / 2)), sig0up)
        / max(out['absTol'], out['relTol'] * abs(hmu0)) #TODO get tolerances from out[err]
    ), 2)))
    # TODO
    # if nmu > out['CM']['nMax']: # don't exceed sample budget
    # # warning(['The algorithm wants to use nmu = ' int2str(nmu)...
    # #          ', which is too big. Using ' int2str(out.CM.nMax) ' instead.'])
    #     nmu = out['CM']['nMax'] # revise nmu

    YY = (y for y in Yrand(nmu)) #get samples for computing the mean

    out['sol'] = np.mean(list(YY))

    out['nSample'] = out['nSig'] + nmu

    out['errBd'] = np.multiply(
        np.multiply(-1, stdnorminv(out['alpha'] / 2)), sig0up / np.sqrt(nmu))

    out['time'] = time.time() - start_time

    return out
