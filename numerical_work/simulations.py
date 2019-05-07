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
        +3*numpy.cos(state[2]-2.*state[3])
        )/4.

def plot_path_from_point(state0):
    start = numpy.array(state0)
    figure, subplots = pyplot.subplots(1,2)
    # print(start[-2:])
    path = utils.RK4_adapt(
        reduced_double_dipole, start, 6.,
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
    ratio = numpy.array(map(lambda arr: arr[2]/arr[3], path))
    # subplots[0].plot(t_vals ,phid_vals)
    # subplots[0].plot(t_vals ,phit_vals)
    # subplots[0].plot(t_vals ,tht_vals)
    subplots[0].plot(t_vals ,force_r)
    # subplots[0].plot(t_vals ,ratio) #between phit and tht
    subplots[1].plot(t_vals ,pphid_vals)
    subplots[1].plot(t_vals ,pphit_vals)
    subplots[1].plot(t_vals ,ptht_vals)
    pyplot.show()

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
    state0p = [0.0, 0.0, 0.5166312725363087, 0.33848195776981216, 0.0, 0.10794186242929231, 0.10008255831561577]
    states = [state0p, state0m]
    for state0_ind in range(1):
        start = numpy.array(states[state0_ind])
        print(start[-2:])
        return_times = utils.return_times_finder(
            reduced_double_dipole, start,
            max_steps=10**6, precision=10**-7, max_time=2*2.8
            )
        print(return_times)
        print([val/return_times[0] for val in return_times])

def orbital_period_sim():
    samples = 100
    # coarse_res = 50
    # coarseT = 0./12.
    maxT = 1./3.
    data_out = open('./%d.txt' % time.time(), 'w')
    start = time.time()
    for frac in range(1, samples):
        # if frac <=coarse_res:
        #     T_0 = (coarseT/coarse_res)*frac
        # else:
        #     T_0 = coarseT + (frac-coarse_res)*(maxT-coarseT)/(samples-coarse_res)
        T_0 = (maxT/samples)*frac
        # T_0 = .5
        print(T_0)
        P_tht = (2.*T_0/7.)**.5
        init = numpy.array([0.0,0.0,0.0,0.0,0.0,-P_tht/2.,P_tht])
        sim_path = utils.RK4_adapt(
            reduced_double_dipole, init, .01,
            max_steps=10**6,precision=10**-10
            )
        S_0 = sim_path[-1]
        S_0[0] = 0.0
        # S_0 = init
        # print(S_0)
        return_times = utils.return_times_finder(
            reduced_double_dipole, S_0,
            max_steps=10**6, precision=10**-9, max_time=75
            )
        # print(return_times)
        # print([val[0]/return_times[0][0] for val in return_times])
        period = utils.fit_slope_and_chisqr(return_times)
        print(period)
        S_0 = list(S_0)
        S_0.append(period[0])
        S_0.append(period[1])
        data_out.write(('%f '*8) % tuple(S_0[1:]))
        data_out.write('\n')
        print(' ')
    data_out.close()
    print('took %d seconds' % (time.time()-start))        # print(len(sim_path), T_0)

def spinning_period_sim():
    samples = 200
    fine_res = 150
    fineT = 2.5/12.
    maxT = 18./12.
    data_out = open('./%d.txt' % time.time(), 'w')
    start = time.time()
    energies = []
    for frac in range(1, samples):
        if frac <=fine_res:
            T_0 = (fineT/fine_res)*frac
        else:
            T_0 = fineT + (frac-fine_res)*(maxT-fineT)/(samples-fine_res)
        energies.append(T_0)
        # T_0 = .5
        print(T_0)
        Pd = (T_0/10.)**.5
        init = numpy.array([0.0,0.0,0.0,0.0,Pd,0.,0.0])
        sim_path = utils.RK4_adapt(
            reduced_double_dipole, init, .05,
            max_steps=10**6,precision=10**-10
            )
        S_0 = sim_path[-1]
        S_0[0] = 0.0
        # S_0 = init
        # print(S_0)
        return_times = utils.return_times_finder(
            reduced_double_dipole, S_0,
            max_steps=10**6, precision=10**-10, max_time=75
            )
        # print(return_times)
        # print([val[0]/return_times[0][0] for val in return_times])
        period = utils.fit_slope_and_chisqr(return_times)
        print(period)
        S_0 = list(S_0)
        S_0.append(period[0])
        S_0.append(period[1])
        data_out.write(('%f '*8) % tuple(S_0[1:]))
        data_out.write('\n')
        print(' ')
    data_out.close()
    # periods = utils.read_column(S_0, -2)
    # matplotlib.plot(energies, periods)
    # pyplot.show()
    print('took %d seconds' % (time.time()-start))        # print(len(sim_path), T_0)



if __name__ == '__main__':
    # sho_plots()
    # mag_oscillation(.001)
    # random_point_examine()
    # random_point_period()
    orbital_period_sim()
    # spinning_period_sim()
    # T_0 = 1./3.-.0001
    # T_0 = 1./2.-.0001 # still maintains contact
    # # T_0 = 7./4.
    # P_tht = (2.*T_0/7.)**.5
    # # init = [0.0,0.0,0.0,0.0,0.0,-P_tht/2.,P_tht]
    # init = [0.0,0.0,0.0,0.0,1.1,0.,0.]
    # plot_path_from_point(init)

'''
choosing a point in phase space:
    pick random phi_t and theta, set phi_d to 0
    calculate U_0
    pick random T_0 such that E_0 < 0
    pick L_0 such that L_0^2 < -9/2 E_0
    determine P_t and P_theta
    ---
    L_0 = P_o + 2P_t
    T = P_o^2 + 10P_t^2
    L_0 = 0
        -> P_t = -P_o/2
        -> T = P_o^2+5/2P_o^2= 7/2P_o^2
        -> (2/7T)^.5 = P_0
period notes
w = 7^.5
period = 2*pi/w = 2.37482
w = (5/3)^.5
period = 2*pi/w = 4.86693

'''
