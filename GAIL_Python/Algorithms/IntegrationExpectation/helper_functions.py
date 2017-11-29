#!/usr/bin/env python
"""Helper functions used throughout the GAIL project

"""

__author__ = ["Anil Simon", "Divya Vasireddy"]

import numpy as np
from scipy.special import erfcinv
import re


def stdnorminv(p):
    """Inverse function of CDF of standard normal distribution
    :param p:
    :return:
    """
    erfinv_val = erfcinv(np.multiply(2, p))
    return np.multiply(-1 * np.sqrt(2), erfinv_val)


def default_random_generator(n):
    """
    X_i^2 where X_i comes from a U[0,1]
    :param n: n is length of random numbers needed
    :return: vector of random numbers
    """
    return np.random.uniform(size=n) ** 2


def identity_function(mu):
    """Identity function that returns what ever is passed to it"""
    return mu


def ci_calculator(mean_estimate, error_bound):
    return mean_estimate - error_bound, mean_estimate + error_bound


def get_param_variants(actual_param):
    """Generates variants of a camel case string"""
    add_underscore = re.compile(r'(.)([A-Z][a-z]+)')
    lowercase_special = re.compile('([a-z0-9])([A-Z])')

    with_underscore = add_underscore.sub(r'\1_\2', actual_param)
    lowercased_underscore = lowercase_special.sub(r'\1_\2', with_underscore)
    return {actual_param, actual_param.lower(), actual_param.upper(),
            with_underscore, with_underscore.lower(), with_underscore.upper(),
            lowercased_underscore, lowercased_underscore.lower(), lowercased_underscore.upper()}


def alias_param_generator(actual_params):
    """Generates a lookup dictionary for parameter variants"""
    result = {}
    for actual_param in actual_params:
        result.update({variant: actual_param for variant in get_param_variants(actual_param)})
    return result

def set_named_args(some_object, kwargs, params_names):
    alias_params = alias_param_generator(params_names)
    for k, v in kwargs.items():
        if k in alias_params:
            setattr(some_object, alias_params[k], v)



def simple_test(property_name, tests=None):
    """Constructs parameter lists for use in parametrised test"""
    if tests is None:
        tests = ['non-numeric', 'negative', 'empty', 'null']

    test_value = {'non-numeric': 'text',
                  'negative': -0.1,
                  'negative large': -0.1,
                  'greater than 1': 2,
                  'empty': '',
                  'null': None,
                  'list': [4, 5, 6],
                  'numeric': 5,
                  'positive fraction': 0.3,
                  'callable': lambda x: x,
                  'float': 2.5,
                  'float large': 2000.332
                  }

    test_parameters_list = []
    for test in tests:
        test_parameters_list.append(
            (property_name, test_value[test],
             '   Test for {} {} = {} '.format(test, str(property_name), str(test_value[test]))))

    return test_parameters_list


def defaultInflateFun(m):
    return np.multiply((16 / 3), np.power(2., (-1 * m)))
