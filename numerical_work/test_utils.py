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
    def setUp(self):
        self.stateI = numpy.array([0.,0.,1.])
        self.path_out = utils.RK4_adapt(base_SHO, self.stateI, 2*numpy.pi, max_steps=500, precision=10**-6)


    def test_meta(self):
        self.assertTrue(True)

    def test_returns_1_time_per_degree_of_freedom(self):
        gam0 = self.stateI
        gam1 = self.path_out[1]
        gam2 = self.path_out[2]
        return_times = utils.return_time(gam0, gam1, gam2)
        self.assertTrue(len(return_times)==(len(gam0)-1))

    def test_finds_previous_times(self):
        gam0 = self.stateI
        gam1 = self.path_out[1]
        gam2 = self.path_out[2]
        return_times = utils.return_time(gam0, gam1, gam2)
        self.assertTrue(return_times[0]<gam1[0])
        self.assertTrue(return_times[1]<gam1[0])

    def test_finds_future_times(self):
        gam0 = self.stateI
        gam1 = self.path_out[-3]
        gam2 = self.path_out[-2]
        return_times = utils.return_time(gam0, gam1, gam2)
        self.assertTrue(return_times[0]>gam1[0])
        self.assertTrue(return_times[1]>gam1[0])

    def test_depends_on_more_than_spatial(self):
        gam0 = self.stateI
        mid_point = len(self.path_out)/2
        gam1 = self.path_out[mid_point-1]
        gam2 = self.path_out[mid_point]
        return_times = utils.return_time(gam0, gam1, gam2)
        time_diff = abs(return_times[0]-return_times[1])
        end_time = self.path_out[-1][0]
        self.assertTrue(
            time_diff>(end_time)/100.
            )

    def test_agrees_on_period(self):
        gam0 = self.stateI
        gam1 = self.path_out[-2]
        gam2 = self.path_out[-1]
        return_times = utils.return_time(gam0, gam1, gam2)
        time_diff = abs(return_times[0]-return_times[1])
        end_time = self.path_out[-1][0]
        avg_time =(return_times[0]+return_times[1])/2.
        self.assertTrue(
            time_diff<10**-6
            )
        self.assertAlmostEqual(
            avg_time,2*numpy.pi
            )



def base_SHO(state):
    deltas = numpy.zeros(len(state))
    deltas[0] = 1.
    deltas[1] = state[2]
    deltas[2] = -state[1]
    return deltas

