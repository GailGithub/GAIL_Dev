# MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
import time
import warnings

import numpy as np
from scipy.special import erfcinv

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.MeanYOut import MeanYOut, stdnorminv
except:
    from MeanYOut import MeanYOut


def check_meanMC_CLT_params(**kwargs):
    checked_params = kwargs

    # TODO meanYOut and meanYParam functions
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


def meanMC_CLT(Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, inflate=1.2, **kwargs):
    start_time = time.time()
    input_kwargs = {'Y': Y, 'absTol': absTol, 'relTol': relTol,
                    'alpha': alpha, 'nSig': nSig, 'inflate': inflate}
    input_kwargs.update(kwargs)
    out = MeanYOut(**input_kwargs)

    Yrand = out.Y  # the random number generator

    p = out.CM.nCV  # the number of control variates
    q = out.nY  # the number of target random variable
    val = np.array(Yrand(out.nSig))  # get samples to estimate variance

    if p == 0 and q == 1:
        YY = val
    else:
        # if there is control variate, construct a new random variable that has the
        # same expected value and smaller variance
        meanVal = np.mean(val, axis=0)  # the mean of each column
        A = val - meanVal  # covariance matrix of the samples

        U, Sdiag, V = np.linalg.svd(
            np.concatenate([A, np.concatenate([np.ones(shape=[1, q]), np.zeros(shape=[1, p])], axis=1)], axis=0),
            full_matrices=False)  # use SVD to solve a constrained least square problem

        U2 = U[-1]  # last row of U

        beta = np.dot(V.T,
                      np.divide(np.divide(U2.T, np.dot(U2, U2.T)), Sdiag))  # get the coefficient for control variates

        # TODO ask V.T, V which one to use for calculating beta.
        YY = np.dot(np.concatenate([val[:, :q], A[:, q:]], axis=1), beta)  # get samples of the new random variable

    out.stddev = np.std(YY)  # standard deviation of the new samples

    sig0up = np.multiply(out.CM.inflate, out.stddev)  # upper bound on the standard deviation

    hmu0 = np.mean(YY)  # mean of the samples

    nmu = int(max(1, np.power(np.ceil(
        np.multiply(
            np.multiply(-1, stdnorminv(out.alpha / 2)), sig0up)
        / max(out.err.absTol, out.err.relTol * abs(hmu0))  # TODO get tolerances from out[err]
    ), 2)))

    if nmu > out.CM.nMax:  # don't exceed sample budget
        warnings.warn(' '.join(['The algorithm wants to use nmu =', str(nmu),
                                ', which is too big. Using', str(out.CM.nMax), 'instead.']))
        nmu = out.CM.nMax  # revise nmu

    W = np.array(Yrand(nmu))  # get samples for computing the mean

    if p == 0 and q == 1:  # no control variates
        YY = W
    else:  # samples of the new random variable
        if out.CM.trueMuCV:
            W[:, q:] = W[:, q:] - out.CM.trueMuCV  # subtract true mean from control variates

        YY = np.dot(W, beta)  # incorporate the control variates and multiple Y's

    sol = np.mean(YY)

    out.sol = sol

    out.nSample = out.nSig + nmu

    out.errBd = np.multiply(
        np.multiply(-1, stdnorminv(out.alpha / 2)), sig0up / np.sqrt(nmu))

    out.time = time.time() - start_time

    return sol, out
