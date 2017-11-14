from MeanYOut import MeanYOut
import pytest


# MeanYOut functions

@pytest.mark.parametrize(
    'stddev,  comment', [
        ('5', '  Comment: stddev is non-numeric'),
        (None, '  Comment: empty stddev'),
        (-2, '  Comment: negative stddev')
    ])
def test_stddev_Exceptions(stddev, comment):
    with pytest.raises(Exception):
        myo = MeanYOut(stddev=stddev)


@pytest.mark.parametrize(
    'property,  property_value, comment', [
        ('stddev', 2, '  Comment: stddev is numeric'),
        ('errBd', 1.5, '  Comment: errBd is numeric'),
        ('nSample', 1000, '  Comment: nSample is numeric values'),
        ('sol', 1.2, '  Comment: sol is float value')
    ])
def test_meanYOut_attribute_setter(property, property_value, comment):
    myo = MeanYOut()
    setattr(myo, property, property_value)
    assert getattr(myo, property) == property_value
