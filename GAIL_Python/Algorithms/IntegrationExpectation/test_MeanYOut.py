#!/usr/bin/env python
"""Tests for validation of MeanYOut input parameters
"""
__author__ = ["Anil Simon", "Divya Vasireddy"]

import pytest
from MeanYOut import MeanYOut

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.helper_functions import default_random_generator, simple_test
except ModuleNotFoundError:
    from helper_functions import default_random_generator, simple_test

default_tests = [
    ('alpha', 0.01, '  Test for default alpha = 0.01'),
    ('absTol', 0.01, '  Test for default absTol = 0.01'),
    ('nSig', 1000, '  Test for default nSig = 1000'),
    ('relTol', 0, '  Test for default relTol = 0'),
    ('Y', default_random_generator, '  Test for default Y = Uniform Square Random generator')
]


@pytest.mark.parametrize('property_name, default_value, comment', default_tests)
def test_meanyout_defaults(property_name, default_value, comment):
    myp = MeanYOut()
    assert getattr(myp, property_name) == default_value


exception_tests = simple_test('Y', ['non-numeric', 'numeric', 'list']) + \
                  simple_test('absTol') + \
                  simple_test('relTol') + \
                  simple_test('relTol', ['negative large']) + \
                  simple_test('alpha') + simple_test('alpha', ['greater than 1']) + \
                  simple_test('nSig')


@pytest.mark.parametrize('property_name, property_value,  comment', exception_tests)
def test_meanyout_exceptions(property_name, property_value, comment):
    with pytest.raises(Exception):
        myp = MeanYOut(**{property_name: property_value})


valid_tests = simple_test('Y', ['callable']) + \
              simple_test('absTol', ['numeric']) + \
              simple_test('relTol', ['numeric']) + \
              simple_test('alpha', ['positive fraction']) + \
              simple_test('nSig', ['numeric'])


@pytest.mark.parametrize('property_name, valid_value, comment', valid_tests)
def test_meanyout_valid(property_name, valid_value, comment):
    myp = MeanYOut(**{property_name: valid_value})
    assert getattr(myp, property_name) == valid_value


meanyout_setter_exceptions_tests = simple_test('stddev', ['non-numeric', 'empty', 'negative']) + \
                                   simple_test('errBd', ['non-numeric', 'empty', 'negative']) + \
                                   simple_test('nSample', ['non-numeric', 'empty', 'negative']) + \
                                   simple_test('sol', ['non-numeric', 'empty'])
@pytest.mark.parametrize(
    'property,  property_value,  comment', meanyout_setter_exceptions_tests)
def test_meanyout_setter_Exceptions(property, property_value, comment):
    with pytest.raises(Exception):
        myo = MeanYOut()
        setattr(myo, property, property_value)
        assert getattr(myo, property) == property_value


@pytest.mark.parametrize(
    'property,  property_value, comment', [
        ('stddev', 2, '  Comment: stddev is numeric'),
        ('errBd', 1.5, '  Comment: errBd is numeric'),
        ('nSample', 1000, '  Comment: nSample is numeric values'),
        ('sol', 1.2, '  Comment: sol is float value')
    ])
def test_meanYOut_setter_valid(property, property_value, comment):
    myo = MeanYOut()
    setattr(myo, property, property_value)
    assert getattr(myo, property) == property_value
