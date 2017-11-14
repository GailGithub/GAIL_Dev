from CubMeanParam import CubMeanParam, defaultInflateFun
import pytest


@pytest.mark.parametrize(
    'property, default_value, comment', [
        ('inflate', 1.2, '  Comment: default inflate = 1.2'),
        ('nInit', 1024, '  Comment: default inflate = 1024'),
        ('nMax', 2 ** 24, '  Comment: default inflate = 2**24'),
        ('nMu', 1, '  Comment: default nMu = 1'),
        ('trueMuCV', [], '  Comment: default trueMuCV = []'),
        ('inflateFun', defaultInflateFun,
         '  Comment: default inflateFun = np.multiply((16 / 3), np.power(2., (-1 * m)))')
    ])
def test_CubMeanParams_defaults(property, default_value, comment):
    cmp = CubMeanParam()
    assert getattr(cmp, property) == default_value


# CubMeanParam functions

@pytest.mark.parametrize(
    'inflate,  comment', [
        ('5', '  Comment: inflate is non-numeric'),
        ('', '  Comment: empty inflate'),
    ])
def test_inflate_Exceptions(inflate, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(inflate=inflate)


@pytest.mark.parametrize(
    'nInit,  comment', [
        ('5', '  Comment: nInit is non-numeric'),
        ('', '  Comment: empty nInit'),
        (-1.1, '  Comment: nInit negative nSig'),
        (2.1, '  Comment: nInit is float')
    ])
def test_nInit_Exceptions(nInit, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(nInit=nInit)


@pytest.mark.parametrize(
    'nMax,  comment', [
        ('5', '  Comment: nMax is non-numeric'),
        ('', '  Comment: empty nMax'),
        (-1.1, '  Comment: nMax negative nSig'),
        (2.1, '  Comment: nMax is float')
    ])
def test_nMax_Exceptions(nMax, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(nMax=nMax)


@pytest.mark.parametrize(
    'nMu,  comment', [
        ('5', '  Comment: nMu is non-numeric'),
        ('', '  Comment: empty nMu'),
        (-1.1, '  Comment: nMu negative nSig'),
        (2.1, '  Comment: nMu is float')
    ])
def test_nMu_Exceptions(nMu, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(nMu=nMu)


@pytest.mark.parametrize(
    'inflateFun,  comment', [
        ('5', '  Comment: inflateFun is non-numeric'),
        ('', '  Comment: inflateFun not empty'),
    ])
def test_inflateFun_Exceptions(inflateFun, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(inflateFun=inflateFun)


@pytest.mark.parametrize(
    'trueMuCV,  comment', [
        ('5', '  Comment: trueMuCV is non-numeric'),
        (['2', '3'], '  Comment: list of non-numeric values')
    ])
def test_trueMuCV_Exceptions(trueMuCV, comment):
    with pytest.raises(Exception):
        cmp = CubMeanParam(trueMuCV=trueMuCV)
