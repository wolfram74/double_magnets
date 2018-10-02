import numpy
import utils
from matplotlib import pyplot
import time

def double_dipole_eqs(state):
    #sliding equations of motion for double dipoles
    deltas = numpy.zeros(len(state))
    del_phi = state[1]-state[2]
    tot_phi = state[1]+state[2]
    deltas[0] = 1. #dt
    deltas[1] = 10.*state[4] #dot phi1
    deltas[2] = 10.*state[5] #dot phi2
    deltas[3] = 2.*state[6] #dot tht
    deltas[4] = -(
        numpy.sin(del_phi)
        +3*numpy.sin(tot_phi-2*state[3])
    )/12. #dot p_phi1
    deltas[5] = (
        numpy.sin(del_phi)
        -3*numpy.sin(tot_phi-2*state[3])
    )/12. #dot p_phi2
    deltas[6] = numpy.sin(tot_phi-2*state[3])/2 #dot p_tht
    return deltas

def reduced_double_dipole(state):
    #sliding equations of motion for double dipoles with relative angles
    deltas = numpy.zeros(len(state))
    deltas[0] = 1. #dt
    deltas[1] = 20.*state[4] #dot phid
    deltas[2] = 20.*state[5] #dot phit
    deltas[3] = 2.*state[6] #dot tht
    deltas[4] = -(
        numpy.sin(state[1])
    )/12. #dot p_phid
    deltas[5] = -(
        numpy.sin(state[2]-2.*state[3])
    )/4. #dot p_phit
    deltas[6] = numpy.sin(state[2]-2*state[3])/2 #dot p_tht
    return deltas


def mag_oscillation(KE = .001):
    # print(KE)
    # p_phi1 = (KE/10.)**.5 #give half the KE to p_phi1
    # p_phi2 = -p_phi1
    # p_tht = 0.
    # T_long = 2*numpy.pi/((5./3.)**.5)
    p_phi1 = (KE/18.)**.5
    p_phi2 = p_phi1
    p_tht = -2*p_phi1
    T_long = 2*numpy.pi/((7.)**.5)
    calced_KE = p_phi1**2*5+p_phi2**2*5+2*p_tht**2
    print(calced_KE, KE)
    print(p_tht**2/2)
    print(T_long)
    init_state = numpy.array([
        0,
        0,0,0,
        p_phi1,p_phi2,p_tht
        ])
    sim_path = utils.RK4_adapt(
        double_dipole_eqs, init_state, T_long*2,
        max_steps=10**6,precision=10**-6
        )
    t_vals = utils.read_column(sim_path, 0)
    phi1_vals = numpy.array(utils.read_column(sim_path, 1))
    phi2_vals = numpy.array(utils.read_column(sim_path, 2))
    tht_vals = numpy.array(utils.read_column(sim_path, 3))
    print(phi1_vals[-1]-tht_vals[-1])
    pyplot.plot(t_vals ,phi1_vals)
    pyplot.plot(t_vals ,phi2_vals)
    pyplot.plot(t_vals ,tht_vals)
    pyplot.show()

def rf_red(state):
    return 2.*state[6]**2-(
        numpy.cos(state[1])
        +numpy.cos(state[2]-2.*state[3])
        )/4.

def random_point_examine():
    #phd, pht, tht, pd, pt, ptht
    state0p, state0m = utils.reduced_state_gen()
    figure, subplots = pyplot.subplots(2,2)
    states = [state0p, state0m]
    print(state0p)
    for state0_ind in range(2):
        start = numpy.array(states[state0_ind])
        print(start[-2:])
        # print(start[-2:])
        path = utils.RK4_adapt(
            reduced_double_dipole, start, 20.,
            max_steps=10**6, precision=10**-7
            )
        t_vals = utils.read_column(path, 0)
        phid_vals = numpy.array(utils.read_column(path, 1))
        phit_vals = numpy.array(utils.read_column(path, 2))
        tht_vals = numpy.array(utils.read_column(path, 3))
        pphid_vals = numpy.array(utils.read_column(path, 4))
        pphit_vals = numpy.array(utils.read_column(path, 5))
        ptht_vals = numpy.array(utils.read_column(path, 6))
        force_r = numpy.array(map(rf_red, path))
        subplots[state0_ind,0].plot(t_vals ,phid_vals)
        subplots[state0_ind,0].plot(t_vals ,phit_vals)
        subplots[state0_ind,0].plot(t_vals ,tht_vals)
        subplots[state0_ind,0].plot(t_vals ,force_r)
        subplots[state0_ind,1].plot(t_vals ,pphid_vals)
        subplots[state0_ind,1].plot(t_vals ,pphit_vals)
        subplots[state0_ind,1].plot(t_vals ,ptht_vals)
    pyplot.show()

def random_point_period():
    state0p, state0m = utils.reduced_state_gen()
    # state0p = [0.0, 0.0, 0.5166312725363087, 0.33848195776981216, 0.0, 0.10794186242929231, 0.10008255831561577]
    states = [state0p, state0m]
    for state0_ind in range(2):
        start = numpy.array(states[state0_ind])
        print('initial momenta')
        print(start[-2:])
        return_times = utils.return_times_finder(
            reduced_double_dipole, start,
            max_steps=10**6, precision=10**-8, max_time=110
            )
        filtered_times = []
        for index in range(len(return_times)):
            if index%2 ==1:
                filtered_times.append(return_times[index])
        # print(filtered_times[0:4])
        # print([val[0]/filtered_times[0][0] for val in filtered_times][0:4])
        period_chisqr = utils.fit_slope_and_chisqr(filtered_times)
        print('best fit period %f' % period_chisqr[0])
        # print('chi sqr %f' % period_chisqr[1])
        print('')
    print('square root of 7 is %f' % 7**.5)

def sample_space():
    start = time.time()
    data_out = open('./%d.txt' % time.time(), 'w')
    for i in range(100):
        state0p = utils.reduced_state_gen()[0]
        start = numpy.array(state0p)
        return_times = utils.return_times_finder(
            reduced_double_dipole, start,
            max_steps=10**6, precision=10**-8, max_time=50.0
            )
        filtered_times = []
        for index in range(len(return_times)):
            if index%2 ==1:
                filtered_times.append(return_times[index])
        period = utils.fit_slope_and_chisqr(filtered_times)
        state0p.append(period[0])
        data_out.write(('%f '*7) % tuple(state0p[1:]))
        data_out.write('\n')
    data_out.close()
    print('took %d seconds' % (time.time()-start))

if __name__ == '__main__':
    # sho_plots()
    # mag_oscillation(.001)
    # random_point_examine()
    # random_point_period()
    sample_space()
'''
choosing a point in phase space:
    pick random phi_t and theta, set phi_d to 0
    calculate U_0
    pick random T_0 such that E_0 < 0
    pick L_0 such that L_0^2 < -9/2 E_0
    determine P_t and P_theta
'''
