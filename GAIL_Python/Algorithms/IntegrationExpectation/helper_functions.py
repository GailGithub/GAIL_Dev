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
