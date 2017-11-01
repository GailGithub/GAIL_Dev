import numbers
from operator import attrgetter

import numpy as np


class CubMeanParam(object):
    """
    CubMeanParam is a class containing the parameters related to
    algorithms that act on functions of x
    This class contains the function, its domain, etc.
    """

    @property
    def nCV(self):
        if self._nCV is None:
            self._nCV = len(self.trueMuCV)
        return self._nCV

    @nCV.setter
    def nCV(self, nCV_value):
        self._nCV = nCV_value

    inflate = property(attrgetter('_inflate'))

    @inflate.setter
    def inflate(self, inflate_value):  # {'function_handle', 'numeric'})
        if inflate_value is None:
            raise Exception("inflate is empty.")
        if not (callable(inflate_value) or isinstance(inflate_value, numbers.Number)):
            raise Exception("inflate should be a function or numeric")

        self._inflate = inflate_value

    inflateFun = property(attrgetter('_inflateFun'))

    @inflateFun.setter
    def inflateFun(self, inflateFun_value):  # {'function_handle', 'numeric'})
        if inflateFun_value is None:
            raise Exception("inflateFun is empty.")
        if not (callable(inflateFun_value) or isinstance(inflateFun_value, numbers.Number)):
            raise Exception("inflateFun should be a function or numeric")

        self._inflateFun = inflateFun_value

        nInit = property(attrgetter('_nInit'))

        @nInit.setter
        def nInit(self, nInit_value):  # {'numeric'}, {'scalar', 'positive', 'integer'})
            msg = "nInit should be positive integer"
            if nInit_value is None:
                raise Exception("nInit value is empty. " + msg)
            if not isinstance(nInit_value, numbers.Number):
                raise Exception("nInit value is not number. " + msg)
            else:
                if not isinstance(nInit_value, int):
                    raise Exception("nInit value is not an integer. " + msg)
                if nInit_value < 0:
                    raise Exception("nInit value is negative. " + msg)

            self._nInit = nInit_value

        nMax = property(attrgetter('_nMax'))

        @nMax.setter
        def nMax(self, nMax_value):  # {'numeric'}, {'scalar', 'positive', 'integer'})
            msg = "nMax should be positive integer"
            if nMax_value is None:
                raise Exception("nMax value is empty. " + msg)
            if not isinstance(nMax_value, numbers.Number):
                raise Exception("nMax value is not number. " + msg)
            else:
                if not isinstance(nMax_value, int):
                    raise Exception("nMax value is not an integer. " + msg)
                if nMax_value < 0:
                    raise Exception("nMax value is negative. " + msg)

            self._nMax = nMax_value

        nMu = property(attrgetter('_nMu'))

        @nMu.setter
        def nMu(self, nMu_value):  # {'numeric'}, {'positive', 'integer'})
            msg = "nMu should be positive integer"
            if nMu_value is None:
                raise Exception("nMu value is empty. " + msg)
            if not isinstance(nMu_value, numbers.Number):
                raise Exception("nMu value is not number. " + msg)
            else:
                if not isinstance(nMu_value, int):
                    raise Exception("nMu value is not an integer. " + msg)
                if nMu_value < 0:
                    raise Exception("nMu value is negative. " + msg)

            self._nMu = nMu_value

        trueMuCV = property(attrgetter('_trueMuCV'))

        @trueMuCV.setter
        def trueMuCV(self, trueMuCV_value):  # {'numeric'}
            if not isinstance(trueMuCV_value, numbers.Number):
                raise Exception("trueMuCV value is not number.")

            self._trueMuCV = trueMuCV_value

    def __init__(self, inflate=1.2, inflateFun=None, nInit=1024, nMax=2 ** 24, nMu=1, trueMuCV=None, **kwargs):
        if inflate is not None: self.inflate = inflate  # inflation factor for bounding the error
        if inflateFun is None:
            self.inflateFun = lambda m: np.multiply((16 / 3), np.power(2., (-1 * m)))
        else:
            self.inflateFun = inflateFun  # inflateFun factor

        self.nInit = nInit  # initial sample size
        self.nMax = nMax  # maximum sample size
        self.nMu = nMu  # number of integrals for solution function
        if trueMuCV is None:
            self.trueMuCV = []
        else:
            self.trueMuCV = trueMuCV  # true integral for control variates

        self._nCV = None  # number of control variates


