from meanYParam import MeanYParam, default_random_generator
import inspect


def test_absTol():
    myp = MeanYParam()
    assert myp.absTol == 0.01
