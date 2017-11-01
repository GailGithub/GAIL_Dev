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


def meanMC_CLT(Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, inflate=1.2):
    start_time = time.time()

    out = MeanYOut(**{'Y': Y, 'absTol': absTol, 'relTol': relTol,
                      'alpha': alpha, 'nSig': nSig, 'inflate': inflate})
    # out = check_meanMC_CLT_params(**{'Y': Y, 'absTol': absTol, 'relTol': relTol,
    #                                  'alpha': alpha, 'nSig': nSig, 'inflate': inflate})

    Yrand = out.Y  # the random number generator
    # Yrand = out['Y']
    p = out.CM.nCV  # the number of control variates
    q = out.nY  # the number of target random variable
    val = Yrand(out.nSig)  # get samples to estimate variance

    if p == 0 and q == 1:
        YY = val
    else:
        # if there is control variate, construct a new random variable that has the
        # same expected value and smaller variance
        # TODO replace with control variate procedures
        #    meanVal = mean(val); %the mean of each column
        #    A = bsxfun(@minus, val, meanVal); %covariance matrix of the samples
        #    [U, S, V] = svd([A; [ones(1,q) zeros(1,p)] ],0); %use SVD to solve a constrained least square problem
        #    Sdiag = diag(S); %the vector of the single values
        #    U2 = U(end,:); %last row of U
        #    beta = V*(U2'/(U2*U2')./Sdiag); %get the coefficient for control variates
        #    YY = [val(:,1:q) A(:,q+1:end)] * beta; %get samples of the new random variable
        meanVal = np.mean(val, axis=0)
        A = np.diff(val, meanVal, axis=1)
        U, S, V = np.linalg.svd(A, full_matrices=True)
        Sdiag = np.diag(S)
        nrow = U.shape[0]
        U2 = U[-1]
        U2_t = np.transpose(U2)
        divide = np.divide(U2_t, np.multiply(U2, U2_t))
        multiply = np.multiply(V, divide)
        beta = np.divide(multiply, Sdiag)

        YY = val

    out.stddev = np.std(YY) # standard deviation of the new samples
    sig0up = np.multiply(out.CM.inflate, out.stddev) # upper bound on the standard deviation
    hmu0 = np.mean(YY) # mean of the samples

    nmu = int(max(1, np.power(np.ceil(
        np.multiply(
            np.multiply(-1, stdnorminv(out.alpha / 2)), sig0up)
        / max(out.err.absTol, out.err.relTol * abs(hmu0))  # TODO get tolerances from out[err]
    ), 2)))

    if nmu > out.CM.nMax:  # don't exceed sample budget
        warnings.warn(' '.join(['The algorithm wants to use nmu =', str(nmu),
                                ', which is too big. Using', str(out.CM.nMax), 'instead.']))
        nmu = out.CM.nMax  # revise nmu

    YY = (y for y in Yrand(nmu))  # get samples for computing the mean

    sol = np.mean(list(YY))
    out.sol = sol

    out.nSample = out.nSig + nmu

    out.errBd = np.multiply(
        np.multiply(-1, stdnorminv(out.alpha / 2)), sig0up / np.sqrt(nmu))

    out.time = time.time() - start_time

    return sol, out
