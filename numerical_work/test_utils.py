# test for returning to position of SHO
# test for amplitude insensitivity of SHO

import unittest
import numpy

import utils

class AdaptiveRK4Test(unittest.TestCase):
    def test_meta(self):
        self.assertTrue(True)
    def test_period_check(self):
        stateI = numpy.array([0.,0.,1.])
        print('pre starting')
        path_out = utils.RK4_adapt(base_SHO, stateI, 2*numpy.pi, max_steps=500, precision=10**-6)
        print(path_out[-1][0]-2*numpy.pi, path_out[-1][1])
        self.assertTrue((path_out[-1][0]-2*numpy.pi) < 10**-5)
        self.assertTrue(abs(path_out[-1][1]) < 10**-5)

class ReturnTimeTest(unittest.TestCase):
    def test_meta(self):
        self.assertTrue(True)


def base_SHO(state):
    deltas = numpy.zeros(len(state))
    deltas[0] = 1.
    deltas[1] = state[2]
    deltas[2] = -state[1]
    return deltas

