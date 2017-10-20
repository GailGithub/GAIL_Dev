from operator import attrgetter
import numbers
import numpy as np
import inspect


def default_random_generator(n):
    """
    X_i^2 where X_i comes from a U[0,1]
    :param n: n is length of random numbers needed
    :return: vector of random numbers
    """
    return np.random.uniform(size=n) ** 2


class MeanYParam(object):
    """
    MeanYParam is a class containing the parameters related to
    algorithms that find the mean of a random variable
    This class contains the random number generator, the uncertainty
    """

    def __init__(self, y=None, abs_tol=1e-2, rel_tol=0, alpha=0.01, n_sig=1000):
        if y is None:
            self.y = default_random_generator
        else:
            self.y = y  # random number generator

        self.abs_tol = abs_tol
        self.rel_tol = rel_tol

        self.alpha = alpha  # uncertainty
        self.n_sig = n_sig  # sample size to estimate variance

        self.err = {}  # an errorParam object #TODO
        self.CM = {}  # a cubMeanParam object #TODO
        self.nY = 1  # number of Y for each mean#TODO
        self.Yout = None  #

    # def __copy__(self):
    #     cls =self.__class__
    #     result = cls.__new__(cls)
    #     result.__dict__.update(self.__dict__)
    #     return result
    #
    # def __deepcopy__(self, memodict={}):
    #     cls = self.__class__
    #     result = cls.__new__(cls)
    #     memodict[id(self)] = result
    #     for k, v in self.__dict__
    #     result.__dict__.update(self.__dict__)
    #     return result



    y = property(attrgetter('_y'))

    @y.setter
    def y(self, y_func):
        # ISFCN To judge if input is a function handle or not
        if not callable(y_func):
            raise Exception("y has to be a function")
        self._y = y_func

    alpha = property(attrgetter('_alpha'))

    @alpha.setter
    def alpha(self, alpha_value):  # {'numeric'}, {'scalar','nonnegative', '<', 1})
        msg = "alpha should be positive and less than 1"
        if not alpha_value:
            raise Exception("value cannot be empty. " + msg)
        if not isinstance(alpha_value, numbers.Number):
            raise Exception("alpha value should be a number. " + msg)
        else:
            if alpha_value < 0:
                raise Exception("alpha value cannot be negative. " + msg)
            elif alpha_value > 1:
                raise Exception("alpha value cannot be greater than 1. " + msg)
        self._alpha = alpha_value

    n_sig = property(attrgetter('_n_sig'))

    @n_sig.setter
    def n_sig(self, n_sig_value):  # {'numeric'}, {'scalar','nonnegative'})
        msg = "n_sig should be a positive number"
        if not n_sig_value:
            raise Exception("value cannot be empty. " + msg)
        if not isinstance(n_sig_value, numbers.Number):
            raise Exception("n_sig value should be a number. " + msg)
        else:
            if n_sig_value < 0:
                raise Exception("n_sig value cannot be negative. " + msg)
        self._n_sig = n_sig_value

    def __str__(self):
        return 'meanYParam with properties:\n\tY\t: {}\n\tabs_tol\t: {}\n\trel_tol\t: {}\n\talpha\t: {}\n\tn_sig\t: {}' \
            .format(inspect.getsource(self.y), self.abs_tol, self.rel_tol, self.alpha, self.n_sig)


if __name__ == '__main__':
    def f(n): return range(n)


    mean_y_parm_inst = MeanYParam(None)

    # mean_y_parm_inst.alpha = -1
    # mean_y_parm_inst.alpha = 'dfs'
    # mean_y_parm_inst.alpha = 5
    mean_y_parm_inst.alpha = 0.25
    # mean_y_parm_inst.alpha = None
    print(mean_y_parm_inst)
