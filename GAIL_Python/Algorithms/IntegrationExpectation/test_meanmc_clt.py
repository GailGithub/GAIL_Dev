import numpy as np
import math
import time
import pytest

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.meanMC_CLT import meanMC_CLT, stdnorminv
except:
    from meanMC_CLT import meanMC_CLT, stdnorminv


# Test case1
def f1(normt, a, d):
    cos_value = np.multiply(np.power((2 * np.pi * np.power(a, 2)), d / 2), np.cos(np.multiply(a, np.sqrt(normt))))
    exp_value = np.exp(np.multiply(0.5 - np.power(a, 2), normt))
    final_value = np.multiply(cos_value, exp_value)
    return final_value


@pytest.mark.parametrize(
    'a_val, expected', [
        (1 / math.sqrt(2), 2.1701),
        (1, 2.1601),
        (1.2, 2.1649)
    ]
)
def test_complex_f1(a_val, expected):
    d_val = 4
    abstol = 0
    reltol = 0.01
    norm_sqd = lambda t: np.sum(np.multiply(t, t), axis=1)
    f = lambda t, a, d: f1(norm_sqd(t), a, d=4)
    out = meanMC_CLT(lambda n: f(np.random.standard_normal(size=(n, d_val)), a_val, d_val),
                     abstol, reltol)
    assert abs(expected - out['sol']) <= expected * reltol


# Test case2
default_random_generator = lambda n: np.random.uniform(size=n) ** 2


def test_default_uniform_sq():
    abstol = 0
    reltol = 0.01
    out = meanMC_CLT(default_random_generator, abstol, reltol)
    exact = 1 / 3;
    assert abs(exact - out['sol']) < 2e-3

# test Case3
