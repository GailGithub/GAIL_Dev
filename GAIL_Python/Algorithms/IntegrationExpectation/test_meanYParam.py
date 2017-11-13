from meanYParam import MeanYParam, default_random_generator
import inspect
import pytest


@pytest.mark.parametrize(
    'property, default_value, comment', [
        ('alpha', 0.01, '  Comment: default alpha = 0.01'),
        ('absTol', 0.01, '  Comment: default absTol = 0.01'),
        ('nSig', 1000, '  Comment: default nSig = 1000'),
        ('relTol', 0, '  Comment: default relTol = 0'),
        ('Y', default_random_generator, '  Comment: default Y = Uniform Square Random generator')
    ])
def test_defaults(property, default_value, comment):
    myp = MeanYParam()
    assert getattr(myp, property) == default_value


@pytest.mark.parametrize(
    'alpha,  comment', [
        (-0.1, '  Comment: negative alpha'),
        ('0.1', '  Comment: non-numeric alpha'),
        (2, '  Comment: greater than 1'),
        ('', '  Comment: empty alpha'),
        (None, '  Comment: null alpha '),
    ])
def test_alpha_excepetions(alpha, comment):
    with pytest.raises(Exception):
        myp = MeanYParam(alpha=alpha)


@pytest.mark.parametrize(
    'absTol,  comment', [
        (-0.1, '  Comment: negative absTol'),
        ('0.1', '  Comment: non-numeric absTol'),
        ('', '  Comment: empty absTol'),
        (None, '  Comment: null absTol '),
    ])
def test_absTol_excepetions(absTol, comment):
    with pytest.raises(Exception):
        myp = MeanYParam(absTol=absTol)


@pytest.mark.parametrize(
    'nSig,  comment', [
        (-0.1, '  Comment: negative nSig'),
        ('0.1', '  Comment: non-numeric nSig'),
        ('', '  Comment: empty nSig'),
        (None, '  Comment: null nSig '),
    ])
def test_nSig_excepetions(nSig, comment):
    with pytest.raises(Exception):
        myp = MeanYParam(nSig=nSig)


@pytest.mark.parametrize(
    'relTol,  comment', [
        (-0.1, '  Comment: negative relTol'),
        ('0.1', '  Comment: non-numeric relTol'),
        ('', '  Comment: empty relTol'),
        (None, '  Comment: null relTol '),
    ])
def test_relTol_excepetions(relTol, comment):
    with pytest.raises(Exception):
        myp = MeanYParam(relTol=relTol)


@pytest.mark.parametrize(
    'Y,  comment', [
        ('5', '  Comment: Y not String'),
        ([5, 6, 6], '  Comment: Y not list'),
        (45, '  Comment: Y not numeric'),
    ])
def test_y(Y, comment):
    with pytest.raises(Exception):
        myp = MeanYParam(Y=Y)
