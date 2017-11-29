#!/usr/bin/env python
"""Tests for validation of setter functions of ErrorParam object"""

__author__ = ["Anil Simon", "Divya Vasireddy"]

import pytest
try:
    from GAIL_Python.Algorithms.IntegrationExpectation.ErrorParam import ErrorParam, ci_calculator, identity_function, default_random_generator
except:
    from ErrorParam import ErrorParam, ci_calculator, identity_function, default_random_generator


try:
    from GAIL_Python.Algorithms.IntegrationExpectation.helper_functions import simple_test
except ModuleNotFoundError:
    from helper_functions import simple_test

error_param_default_tests = [
    ('absTol', 0.01, '  Comment: default absTol = 0.01'),
    ('relTol', 0, '  Comment: default relTol = 0'),
    ('solFun', identity_function, '  Comment: default solFun = identity Function'),
    ('solBdFun', ci_calculator, '  Comment: default solBdFun = CI calculator')]


@pytest.mark.parametrize(
    'property, default_value, comment', error_param_default_tests)
def test_ErrorParam_defaults_(property, default_value, comment):
    myp = ErrorParam()
    assert getattr(myp, property) == default_value



error_exception_tests = simple_test('solFun', ['non-numeric', 'numeric', 'list']) + \
                  simple_test('solBdFun', ['non-numeric', 'numeric', 'list']) + \
                  simple_test('absTol') + \
                  simple_test('relTol')


@pytest.mark.parametrize('property_name, property_value,  comment', error_exception_tests)
def test_errorparam_exceptions(property_name, property_value, comment):
    with pytest.raises(Exception):
        myp = ErrorParam(**{property_name: property_value})




valid_tests = simple_test('solFun', ['callable']) + \
              simple_test('solBdFun', ['callable']) + \
              simple_test('absTol', ['numeric']) + \
              simple_test('relTol', ['numeric'])

@pytest.mark.parametrize('property_name, valid_value, comment', valid_tests)
def test_errorparam_valid(property_name, valid_value, comment):
    myp = ErrorParam(**{property_name: valid_value})
    assert getattr(myp, property_name) == valid_value
