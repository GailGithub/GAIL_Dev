import numpy as np
import math
import time
from GAIL_Python.Algorithms.IntegrationExpectation.meanMC_CLT import meanMC_CLT,stdnorminv
# from GAIL_Python.Algorithms.IntegrationExpectation.TestCases import normsqd_fun
def normsqd_fun(t, a, d):
    normt = sum(np.multiply(t,t),axis = 0)
    cos_value = np.multiply(np.power((2 * math.pi * a ^ 2), d / 2), math.cos(np.multiply(a, math.sqrt(normt))))
    exp_value = math.exp(np.multiply(np.diff(1 / 2, np.power(a, 2))),normt)
    final_value = np.multiply(cos_value, exp_value)
    return final_value

d=4
abstol = 0.01
reltol = 0.0
avec=[1/math.sqrt(2),1,1.2]
for a in avec:
    start_time = time.time()
    mu, sigma = 0,1
    normal = lambda n: np.random.normal(mu, sigma, n)
    f = normsqd_fun(normal, a,d)
    print(meanMC_CLT(f, abstol, reltol))



