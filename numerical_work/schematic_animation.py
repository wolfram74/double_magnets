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

def frame_gen(state, axes):
    axes.clear()
    phi1 = (state[2]+state[1])/2.
    w1 = 10*(state[5]+state[4])
    phi2 = (state[2]-state[1])/2.
    w2 = 10*(state[5]-state[4])
    tht = state[3]
    wtht = 2*state[6]
    # print(w1, w2, wtht)
    c1 = numpy.array([.5*numpy.cos(tht), .5*numpy.sin(tht)])
    c2 = -c1
    mu1 = c1+numpy.array([numpy.cos(phi1), numpy.sin(phi1)])*.5
    mu2 = c2+numpy.array([numpy.cos(phi2), numpy.sin(phi2)])*.5
    colors = []
    mu1_arrow = mpatches.Arrow(
        c1[0],c1[1], cos(phi1)/2., sin(phi1)/2.
        , width=.2, zorder=1.
        )
    mu1_terminus = [c1[0]+cos(phi1)/2, c1[1]+sin(phi1)/2]
    w1_terminus = [c1[0]+cos(phi1+w1)/2, c1[1]+sin(phi1+w1)/2]
    colors.append((.1,.5,.1))
    circ_1 = mpatches.Circle(c1, radius=.5, alpha=.1, zorder=.5)
    colors.append((.5,.1,.1,.1))
    mu2_arrow = mpatches.Arrow(
        c2[0],c2[1], cos(phi2)/2., sin(phi2)/2., width=.2, zorder=1.
        )
    mu2_terminus = [c2[0]+cos(phi2)/2, c2[1]+sin(phi2)/2]
    w2_terminus = [c2[0]+cos(phi2+w2)/2, c2[1]+sin(phi2+w2)/2]
    colors.append((.1,.5,.1))
    circ_2 = mpatches.Circle(c2, radius=.5, alpha=.1, zorder=.5)
    colors.append((.5,.1,.1,.1))
    style="Simple,tail_width=0.02,head_width=.1,head_length=.1"
    kw = dict(arrowstyle=style, color="k")
    vel_color = (.1, .1, .1, 1)
    vel_scale = .28
    del_tht = [-wtht*vel_scale*sin(tht), wtht*vel_scale*cos(tht) ]

    wtht_arrow1 = mpatches.Arrow(
        c1[0],c1[1], del_tht[0], del_tht[1],
        width=.2,zorder=1.
        )
    colors.append(vel_color)
    wtht_arrow2 = mpatches.Arrow(
        c2[0],c2[1], -del_tht[0], -del_tht[1],
        width=.2,zorder=1.
        )
    colors.append(vel_color)
    w1_arrow = mpatches.FancyArrowPatch(
        mu1_terminus,
        w1_terminus,
        connectionstyle='arc3,rad=%f'%(vel_scale*w1), **kw
        )
    colors.append(vel_color)
    w2_arrow = mpatches.FancyArrowPatch(
        mu2_terminus,
        w2_terminus,
        connectionstyle='arc3,rad=%f'%(vel_scale*w2), **kw
        )
    colors.append(vel_color)
    shapes = PatchCollection(
        [mu1_arrow, circ_1, mu2_arrow, circ_2,
        wtht_arrow1,wtht_arrow2, w1_arrow, w2_arrow
        ],
        facecolors = colors
        )
    axes.add_collection(shapes)
    t_T0 = state[0]/(2*numpy.pi/7**.5)
    # axes.set_title('t=%.3f, t/T0=%.3f' % (state[0], t_T0))
    axes.set_title('t=%.3f' % (state[0]))

    # pyplot.show()

def simulate_initial(init_state, T_long=(24*numpy.pi/((7.)**.5))):
    return utils.RK4_adapt(
        reduced_double_dipole, init_state, T_long,
        max_steps=10**6,precision=10**-6
        )

def custom_animation():
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

def animatic(frames=5):
    fig, axes = pyplot.subplots(frames,1)
    #orbital mode settings
    # T_0 = 4.0 #orbital period
    # E_0 = .123 #orbital energy
    # pt = ((1./3.+E_0)/14.)**.5 #orbital momentum
    # state = numpy.array([0.,0.,0.,0.,0.,pt,-2.*pt]) # orbital mode

    #spinning mode settings
    T_0 = 8.0 #spinning period
    E_0 = .1505 #spinning kinetic energy
    pd = (E_0/10.)**.5 #spinning momentum
    state = numpy.array([0.,0.,0.,0.,pd,0.,0.]) # spinning mode

    path = simulate_initial(state,T_0)
    total_frames = len(path)
    tique_count = 0
    frame_num = 0
    while tique_count<frames:
        state_i = path[frame_num]
        if(state_i[0] < tique_count*(T_0/frames)):
            frame_num+=1
            continue
        i = tique_count
        axes[i].set_ylim(-1, 1)
        axes[i].set_xlim(-1, 1)
        axes[i].get_xaxis().set_visible(False)
        axes[i].get_yaxis().set_visible(False)
        axes[i].set_axis_off()
        axes[i].set_aspect(aspect='equal')
        print(tique_count, frame_num, state_i)
        frame_gen(state_i, axes[i])
        tique_count+=1
    # every4th = path[::4]
    pyplot.show()
    print('s0=%f_%f_%f_%f_%f_%f.mp4' % tuple(state[1:]))

if __name__=='__main__':
    fig, axes = pyplot.subplots()

    # custom_animation()
    animatic(4)

'''
L = pth + 2pt
L=0 -> pth = -2pt
K = pth**2 + 10pt**2
L=0 -> K = 14pt**2
1/3 = 14pt**2

orbital mode with T=4
0.000000 -0.036117 0.007223
-(1+3*cos(-0.0361-2*0.0072))/12
-.333
0.000000 -0.180544 0.361088
K =.456
E0 =.123
'''
