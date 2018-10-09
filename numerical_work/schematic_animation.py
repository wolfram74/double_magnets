import numpy
import utils
import matplotlib.animation as animation
import matplotlib.pyplot as pyplot
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

cos, sin = numpy.cos, numpy.sin

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

def frame_gen(state):
    axes.clear()
    phi1 = (state[2]+state[1])/2.
    phi2 = (state[2]-state[1])/2.
    tht = state[3]
    c1 = numpy.array([.5*numpy.cos(tht), .5*numpy.sin(tht)])
    c2 = -c1
    mu1 = c1+numpy.array([numpy.cos(phi1), numpy.sin(phi1)])*.5
    mu2 = c2+numpy.array([numpy.cos(phi2), numpy.sin(phi2)])*.5
    colors = []
    mu1_arrow = mpatches.Arrow(
        c1[0],c1[1], cos(phi1)/2., sin(phi1)/2.
        , width=.2, zorder=1.
        )
    colors.append((.1,.5,.1))
    circ_1 = mpatches.Circle(c1, radius=.5, alpha=.1, zorder=.5)
    colors.append((.5,.1,.1,.1))
    mu2_arrow = mpatches.Arrow(
        c2[0],c2[1], cos(phi2)/2., sin(phi2)/2., width=.2, zorder=1.
        )
    colors.append((.1,.5,.1))
    circ_2 = mpatches.Circle(c2, radius=.5, alpha=.1, zorder=.5)
    colors.append((.5,.1,.1,.1))
    shapes = PatchCollection(
        [mu1_arrow, circ_1, mu2_arrow, circ_2],
        facecolors = colors
        )
    axes.add_collection(shapes)
    t_T0 = state[0]/(2*numpy.pi/7**.5)
    axes.set_title('t=%f, t/T0=%f' % (state[0], t_T0))

    # pyplot.show()
    # pyplot.plot(t_vals ,phi1_vals)
    # pyplot.plot(t_vals ,phi2_vals)
    # pyplot.plot(t_vals ,tht_vals)

def simulate_initial(init_state):
    # p_phi1 = (energy/10.)**.5 #give half the KE to p_phi1
    # p_phi2 = -p_phi1
    # p_tht = 0.
    # T_long = 2*numpy.pi/((5./3.)**.5)
    # p_phi1 = (energy/12.)**.5
    # p_phi2 = p_phi1
    # p_tht = -2*p_phi1
    T_long = 12*numpy.pi/((7.)**.5)
    # init_state = numpy.array([
    #     0,
    #     0,0,0,
    #     p_phi1,p_phi2,p_tht
    #     ])
    return utils.RK4_adapt(
        reduced_double_dipole, init_state, T_long*2,
        max_steps=10**6,precision=10**-6
        )

if __name__=='__main__':
    # data = load_data('../produced_data/tht_0.6666_w_137/t_0.txt')
    fig, axes = pyplot.subplots()
    axes.set_ylim(-1, 1)
    axes.set_xlim(-1, 1)
    axes.set_aspect(aspect='equal')
    # frame_gen([0.,0.,.5,-0.25,0.,0.,0.])
    # pt = (3./84-.002)**.5
    pt = (1./28.-.001)**.5
    state = numpy.array([0.,0.,0.,0.,0.,pt,-2.*pt])
    path = simulate_initial(state)
    every4th = path[::4]
    print('s0=%f_%f_%f_%f_%f_%f.mp4' % tuple(state[1:]))
    movie = animation.FuncAnimation(
        fig,
        frame_gen,
        frames=every4th
        )
    movie.save('s0=%f_%f_%f_%f_%f_%f.mp4' % tuple(state[1:]))

'''
L = pth + 2pt
L=0 -> pth = -2pt
K = pth**2 + 10pt**2
L=0 -> K = 14pt**2
1/3 = 14pt**2
'''
