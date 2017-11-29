#!/usr/bin/env python

""" Definition of meanMC_CLT function
    Definition of a Monte Carlo method to calculate the mean of a random variable
    using the central limit theorem
    :author
"""

__author__ = ["Anil Simon", "Divya Vasireddy"]

import time
import warnings

import numpy as np

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.helper_functions import stdnorminv
except ModuleNotFoundError:
    from helper_functions import stdnorminv

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.MeanYOut import MeanYOut
except ModuleNotFoundError:
    from MeanYOut import MeanYOut


def meanMC_CLT(Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, inflate=1.2, **kwargs):
    """Returns an estimate of the mean of a random variable
    meanMC_CLT is a Monte Carlo method to estimate the mean of a random variable
 
    sol = meanMC_CLT(Y,absTol,relTol,alpha,nSig,inflate) estimates the
    mean, mu, of a random variable to within a specified error tolerance,
    i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
    1-alpha, where abstol is the absolute error tolerance.  The default
    values are abstol=1e-2 and alpha=1 . Input Y is a function handle that
    accepts a positive integer input n and returns an n x 1 vector of IID
    instances of the random variable.
 
    This is a heuristic algorithm based on a Central Limit Theorem
    approximation.


    :param Y: The function or structure for generating n IID instances of a
            random variable Y whose mean we want to estimate. Y is often defined
            as a function of some random variable X with a simple distribution.
            The input of Yrand should be the number of random variables n, the
            output of Yrand should be n function values. For example, if Y = X.^2
            where X is a standard uniform random variable, then one may define
            Yrand = @(n) rand(n,1).^2.
        
    :param absTol: the absolute error tolerance, which should be
            non-negative --- default = 1e-2
     
    :param relTol: the relative error tolerance, which should be
            non-negative and no greater than 1 --- default = 0

    :param alpha: the uncertainty, which should be a small positive
            percentage --- default = 1%

    :param nSig: the number of samples used to compute the sample variance
            --- default = 1000

    :param inflate: the standard deviation inflation factor --- default = 1.2

    :param kwargs: Any additional named arguments that are input arguments
            to MeanYParam, ErrorParam or CubMeanParam



    :return:

    """
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

        # TODO ask V.T, V which one to use for calculating beta.
        beta = np.dot(V.T,
                      np.divide(np.divide(U2.T, np.dot(U2, U2.T)), Sdiag))  # get the coefficient for control variates

        YY = np.dot(np.concatenate([val[:, :q], A[:, q:]], axis=1), beta)  # get samples of the new random variable

    out.stddev = np.std(YY)  # standard deviation of the new samples

    sig0up = np.multiply(out.CM.inflate, out.stddev)  # upper bound on the standard deviation

    hmu0 = np.mean(YY)  # mean of the samples

    nmu = int(max(1, np.power(np.ceil(
        np.multiply(
            np.multiply(-1, stdnorminv(out.alpha / 2)), sig0up)
        / max(out.err.absTol, out.err.relTol * abs(hmu0))
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
