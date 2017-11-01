import numpy as np
import math
import pytest

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.meanMC_CLT import meanMC_CLT, stdnorminv
except:
    from meanMC_CLT import meanMC_CLT, stdnorminv

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.MeanYOut import MeanYOut
except:
    from MeanYOut import MeanYOut


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
    sol, out = meanMC_CLT(lambda n: f(np.random.standard_normal(size=(n, d_val)), a_val, d_val),
                          abstol, reltol)
    assert abs(expected - sol) <= expected * reltol


# Test case2
default_random_generator = lambda n: np.random.uniform(size=n) ** 2


def test_default_uniform_sq():
    abstol = 0
    reltol = 0.01
    sol, out = meanMC_CLT(default_random_generator, abstol, reltol)
    exact = 1 / 3
    assert abs(exact - sol) < exact * reltol


# Test Case3

def test_CV_1():
    f = lambda x: np.hstack([np.exp(-1 * (np.power(x, 2))), x])
    YXn = lambda n: f(np.random.uniform(size=n))
    inp = dict(Y=YXn, nY=1, trueMuCV=1 / 2, absTol=0, relTol=1e-3)
    hmu, out = meanMC_CLT(**inp)
    exact = math.erf(1) * math.sqrt(math.pi) / 2
    assert abs(exact - hmu) < max(0, 1e-3 * abs(exact))


# Test Case 4
# Estimate the Keister's integration in dimension 1 with a=1, 1/sqrt(2) and
# using cos(x) as a control variate:
def test_cv_2():
    a = 1
    d = 1
    normsqd = lambda x: np.sum(np.dot(x, x), axis=1)
    f = lambda normt, a, d: np.power((2 * math.pi * a ** 2), (d / 2)) * np.dot(math.cos(a * math.sqrt(normt)),
                                                                               math.exp((1 / 2 - a ** 2) * normt))
    f1 = lambda x, a, d: f(normsqd(x), a, d)
    f2 = lambda x: [f1(x, 1, 1), f1(x, 1 / math.sqrt(2), 1), math.cos(x)]
    YXn = lambda n: f2(np.random.uniform(size=n))
    s = dict(Y=YXn, nY=2, trueMuCV=1 / math.sqrt(math.exp(1)), absTol=0, relTol=1e-3)
    hmu, out = meanMC_CLT(**s)
    exact = 1.380388447043143
    assert abs(exact - hmu) < max(0, 1e-3 * abs(exact))


# Test Case 5
# Estimate the integral with integrand f(x) = x1.^3.*x2.^3.*x3.^3 in the
# interval [0,1]^3 with pure absolute error 1e-3 using x1.*x2.*x3 as a
# control variate:

def test_cv_3():
    f = lambda x: [np.power(x[:, 1], 3) * np.power(x[:, 2], 3) * np.power(x[:, 3], 3), x[:, 1] * x[:, 2] * x[:, 3]]
    s = dict(Y=lambda n: f(np.random.uniform(size=n)), nY=1, trueMuCV=1 / 8, absTol=1e-3, relTol=0)
    hmu, out = meanMC_CLT(*s)
    exact = 1 / 64
    assert abs(exact - hmu) < max(1e-3, 1e-3 * abs(exact))


# test_case 6
# Estimate the integrals with integrands f1(x) = x1.^3.*x2.^3.*x3.^3 and
# f2(x)= x1.^2.*x2.^2.*x3.^2-1/27+1/64 in the interval [0,1]^3 using
# x1.*x2.*x3 and x1+x2+x3 as control variates:

def test_cv_4():
    f = lambda x: [np.power(x[:, 1], 3) * np.power(x[:, 2], 3) * np.power(x[:, 3], 3),
                   np.power(x[:, 1], 2) * np.power(x[:, 2], 2) * np.power(x[:, 3], 2) - 1 / 27 + 1 / 64,
                   x[:, 1] * x[:, 2] * x[:, 3], x[:, 1] + x[:, 2] + x[:, 3]];
    s = dict(Y=lambda n: f(np.random.uniform(size=(n, 3)), nY=2, trueMuCV=[1 / 8, 1.5], absTol=1e-4, relTol=1e-3))
    hmu, out = meanMC_CLT(**s)
    exact = 1 / 64
    assert abs(exact - hmu) < max(1e-4, 1e-3 * abs(exact))
