import numpy as np
import math
import time
from GAIL_Python.Algorithms.IntegrationExpectation.meanMC_CLT import meanMC_CLT
d=4
abstol = 0.01
reltol = 0.0
avec=[1/math.sqrt(2),1,1.2]
for a in avec:
    start_time = time.time()
    mu, sigma = 0,1
    f = lambda n: np.random.normal(mu, sigma, n)
    print(meanMC_CLT(f, a, abstol, reltol))