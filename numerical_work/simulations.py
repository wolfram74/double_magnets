import numpy
import utils
from matplotlib import pyplot

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

if __name__ == '__main__':
    # sho_plots()
    mag_oscillation(.001)
