import pytest
try:
    from GAIL_Python.Algorithms.IntegrationExpectation.meanYParam import MeanYParam, default_random_generator
except:
    from meanYParam import MeanYParam, default_random_generator


def simple_test(property_name, tests=None):
    """Constructs parameter lists for use in parametrised test"""
    if tests is None:
        tests = ['non-numeric', 'negative', 'empty', 'null']

    test_value = {'non-numeric': 'text',
                  'negative': -0.1,
                  'greater than 1': 2,
                  'empty': '',
                  'null': None,
                  'list': [4, 5, 6],
                  'numeric': 5,
                  'positive fraction': 0.3,
                  'callable': lambda x: x
                  }

    test_parameters_list = []
    for test in tests:
        test_parameters_list.append(
            (property_name, test_value[test],
             '   Test for {} {} = {} '.format(test, str(property_name), str(test_value[test]))))

    return test_parameters_list


default_tests = [
    ('alpha', 0.01, '  Test for default alpha = 0.01'),
    ('absTol', 0.01, '  Test for default absTol = 0.01'),
    ('nSig', 1000, '  Test for default nSig = 1000'),
    ('relTol', 0, '  Test for default relTol = 0'),
    ('Y', default_random_generator, '  Test for default Y = Uniform Square Random generator')
]


@pytest.mark.parametrize('property_name, default_value, comment', default_tests)
def test_meanyparam_defaults(property_name, default_value, comment):
    myp = MeanYParam()
    assert getattr(myp, property_name) == default_value


exception_tests = simple_test('Y', ['non-numeric', 'numeric', 'list']) + \
                  simple_test('absTol') + \
                  simple_test('relTol') + \
                  simple_test('alpha') + simple_test('alpha', ['greater than 1']) + \
                  simple_test('nSig')


@pytest.mark.parametrize('property_name, property_value,  comment', exception_tests)
def test_meanyparam_exceptions(property_name, property_value, comment):
    with pytest.raises(Exception):
        myp = MeanYParam(**{property_name: property_value})


valid_tests = simple_test('Y', ['callable']) + \
              simple_test('absTol', ['numeric']) + \
              simple_test('relTol', ['numeric']) + \
              simple_test('alpha', ['positive fraction']) + \
              simple_test('nSig', ['numeric'])


@pytest.mark.parametrize('property_name, valid_value, comment', valid_tests)
def test_meanyparam_valid(property_name, valid_value, comment):
    myp = MeanYParam(**{property_name: valid_value})
    assert getattr(myp, property_name) == valid_value
