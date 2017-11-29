#!/usr/bin/env python
"""Tests for validation of CubMeanParam input parameters
"""
__author__ = ["Anil Simon", "Divya Vasireddy"]

import pytest
from CubMeanParam import CubMeanParam

try:
    from GAIL_Python.Algorithms.IntegrationExpectation.helper_functions import defaultInflateFun, simple_test
except ModuleNotFoundError:
    from helper_functions import defaultInflateFun, simple_test

cubmeanparam_defaults_tests = [
    ('inflate', 1.2, '  Comment: default inflate = 1.2'),
    ('nInit', 1024, '  Comment: default nInit = 1024'),
    ('nMax', 2 ** 24, '  Comment: default nMax = 2**24'),
    ('nMu', 1, '  Comment: default nMu = 1'),
    ('trueMuCV', [], '  Comment: default trueMuCV = []'),
    ('inflateFun', defaultInflateFun,
     '  Comment: default inflateFun = np.multiply((16 / 3), np.power(2., (-1 * m)))')
]


@pytest.mark.parametrize(
    'property, default_value, comment', cubmeanparam_defaults_tests)
def test_CubMeanParams_defaults(property, default_value, comment):
    cmp = CubMeanParam()
    assert getattr(cmp, property) == default_value


@pytest.mark.parametrize(
    'inflate,  comment', [
        ('5', '  Comment: inflate is non-numeric'),
        ('', '  Comment: empty inflate'),
    ])
def test_inflate_Exceptions(inflate, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(inflate=inflate)


exception_tests = simple_test('nInit', ['non-numeric', 'empty', 'null', 'negative', 'negative large', 'float']) + \
                  simple_test('nmax', ['non-numeric', 'empty', 'null', 'negative', 'negative large', 'float']) + \
                  simple_test('nMu', ['non-numeric', 'empty', 'null', 'negative', 'negative large', 'float']) + \
                  simple_test('inflateFun', ['non-numeric', 'list']) + \
                  simple_test('trueMuCV', ['non-numeric'])


@pytest.mark.parametrize('property_name, property_value,  comment', exception_tests)
def test_cubmeanparam_exceptions(property_name, property_value, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(**{property_name: property_value})


cubmeanparam_valid_values = [
    ('inflate', 1.3, 1.3, '  Comment: valid inflate = 1.3'),
    ('nInit', 10000, 10000, '  Comment: valid nInit = 10000'),
    ('nMax', 2 ** 15, 2 ** 15, '  Comment: valid nMax = 2**15'),
    ('trueMuCV', 5, [5], '  Comment: valid trueMuCV = 5'),
    ('trueMuCV', [2, 3], [2, 3], '  Comment: valid trueMuCV = [2,3]')
]


@pytest.mark.parametrize(
    'property, valid_input, test_value, comment', cubmeanparam_valid_values)
def test_CubMeanParams_valid(property, valid_input, test_value, comment):
    cmp = CubMeanParam()
    setattr(cmp, property, valid_input)
    assert getattr(cmp, property) == test_value
