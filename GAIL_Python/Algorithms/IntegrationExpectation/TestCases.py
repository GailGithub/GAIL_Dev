import numpy as np
import math
import time
from GAIL_Python.Algorithms.IntegrationExpectation.meanMC_CLT import meanMC_CLT,stdnorminv
def f1(normt, a, d):
    cos_value = np.multiply(np.power((2 * np.pi * np.power(a, 2)), d / 2), np.cos(np.multiply(a, np.sqrt(normt))))
    exp_value = np.exp(np.multiply(0.5 - np.power(a, 2), normt))
    final_value = np.multiply(cos_value, exp_value)
    return final_value
d_val=4
abstol = 0
reltol = 0.01
avec=[1/math.sqrt(2),1,1.2]
mu, sigma = 0, 1
norm_sqd = lambda t: np.sum(np.multiply(t, t), axis=1)
f = lambda t, a, d: f1(norm_sqd(t), a, d=4)
for a_val in avec:
    print(meanMC_CLT(lambda n: f(np.random.standard_normal(size=(n,d_val)),a_val,d_val), abstol, reltol))

