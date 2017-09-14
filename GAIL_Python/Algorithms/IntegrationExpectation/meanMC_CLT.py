# MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
import time

import numpy as np
from scipy.special import erfcinv


def check_meanMC_CLT_params(**kwargs):
    checked_params = kwargs

    # if not (0 <= relTol <= 1):
    #     logging.error("Relative error tolerance should be in the range [0,1]")
    return checked_params


def stdnorminv(p):
    # this function is the inverse function of CDF of standard normal distribution
    erfinv_val = erfcinv(np.multiply(2, p))
    return np.multiply(-1 * np.sqrt(2), erfinv_val)


def meanMC_CLT(Y=None, absTol=1e-2, relTol=0, alpha=0.01, nSig=1000, inflate=1.2):
    start_time = time.time()

    out = check_meanMC_CLT_params(**{'Y': Y, 'absTol': absTol, 'relTol': relTol,
                                     'alpha': alpha, 'nSig': nSig, 'inflate': inflate})

    Yrand = out['Y']
    p = 0  # out['CM']['nCV']
    q = 1  # out['nY']
    val = (for y in Yrand(out['nSig']))

    if p == 0 and q == 1:
        YY = list(val)
    else:
        # TODO replace with control variate procedures
        YY = list(val)

    out['stddev'] = np.std(YY)
    sig0up = np.multiply(out['inflate'], out['stddev'])  # TODO change inflate factor out.CM.inflate
    hmu0 = np.mean(YY)

    nmu = max(1, np.power(np.ceil(
        np.multiply(
            np.multiply(-1, stdnorminv(out['alpha'] / 2)), sig0up)
        / max(out['absTol'], out['relTol'] * abs(hmu0))
    ), 2))

    out['sol'] = np.mean(YY)

    out['nSample'] = out['nSig'] + nmu

    out['errBd'] = np.multiply(
        np.multiply(-1, stdnorminv(out['alpha'] / 2)), sig0up / np.sqrt(nmu))

    out['time'] = time.time() - start_time

    return out
