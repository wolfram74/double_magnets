import numpy
import utils
import matplotlib.animation as animation
import matplotlib.pyplot as pyplot

def double_dipole_eqs(state):
    #sliding equations of motion for double dipoles
    deltas = numpy.zeros(len(state))
    del_phi = state[1]-state[2]
    tot_phi = state[1]+state[2]
    deltas[0] = 1. #dt
    deltas[1] = 10.*state[4] #dot phi1
    deltas[2] = 10.*state[5] #dot phi2
    deltas[3] = state[6] #dot tht
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

def frame_gen(energy):
    # del x_vals[:]
    # del y_vals[:]
    axes.clear()
    sim_path = simulate_energy(energy)
    t_vals = numpy.array(utils.read_column(sim_path, 0))
    phi1_vals = numpy.array(utils.read_column(sim_path, 1))
    phi2_vals = numpy.array(utils.read_column(sim_path, 2))
    tht_vals = numpy.array(utils.read_column(sim_path, 3))
    pyplot.plot(t_vals ,phi1_vals)
    pyplot.plot(t_vals ,phi2_vals)
    pyplot.plot(t_vals ,tht_vals)
    axes.set_ylim(-2., 2.)
    axes.set_title('initial KE=%f' % energy)
    # axes.set_xlim(0, 170)

def simulate_energy(energy):
    p_phi1 = (energy/10.)**.5 #give half the KE to p_phi1
    p_phi2 = -p_phi1
    p_tht = 0.
    T_long = 2*numpy.pi/((5./3.)**.5)
    init_state = numpy.array([
        0,
        0,0,0,
        p_phi1,p_phi2,p_tht
        ])
    return utils.RK4_adapt(
        double_dipole_eqs, init_state, T_long*2,
        max_steps=10**6,precision=10**-6
        )

if __name__=='__main__':
    # data = load_data('../produced_data/tht_0.6666_w_137/t_0.txt')
    fig, axes = pyplot.subplots()
    x_vals = [.5]
    y_vals = [.2]
    field = axes.scatter(x_vals, y_vals)
    axes.set_ylim(-1, 1)
    axes.set_xlim(0, 1)
    energies = [1./3*(i)/200 for i in range(1, 101)]
    print(energies[-1])
    movie = animation.FuncAnimation(
        fig,
        frame_gen,
        frames=energies
        )
    movie.save('increasing_energy.mp4')
