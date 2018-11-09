import matplotlib.pyplot as pyplot
import matplotlib.text as text_drawer
import utils
import numpy


def parse_file(target):
    data_in = open(target, 'r')
    output = []
    for line in data_in:
        vals = [float(num) for num in line.rstrip().split(' ')]
        state = vals[:-2]
        period = vals[-2]
        unc = vals[-1]
        L0 = 2*state[4]+state[5]
        T0 = (20*state[3]**2+20*state[4]**2+2*state[5]**2)/2.
        U0 = -(numpy.cos(state[0])+3*numpy.cos(state[1]-2*state[2]))/12.
        E0 = T0+U0
        output.append((period, unc, E0, L0, T0, U0))
    return output

def plot_period_vs_energy():
    # plot spinning periods
    data_spinning = parse_file('./clean_spinnning_periods.txt')
    p_vals = utils.read_column(data_spinning, 0)
    sig_vals = utils.read_column(data_spinning, 1)
    E_vals = utils.read_column(data_spinning, 2)
    pyplot.plot(E_vals, p_vals,

        label='Large Amplitude Spinning Mode')
    T_s = numpy.pi*2*(5./3.)**(-.5)
    E_vals = numpy.linspace(-1./3., -1/6., 20)
    ones = numpy.ones(20)
    pyplot.plot(E_vals, T_s*ones, 'b--', label='Small Amplitude spinning mode')
    #plot orbital periods
    data_orbital = parse_file('./clean_orbital_periods.txt')
    p_vals = utils.read_column(data_orbital, 0)
    sig_vals = utils.read_column(data_orbital, 1)
    E_vals = utils.read_column(data_orbital, 2)
    p2_vals = map(lambda x: x*2, p_vals)
    pyplot.plot(E_vals, p_vals
        ,
        label='Large Amplitude Orbital Mode')
    # pyplot.plot(E_vals, p2_vals)

    # plot pendulum periods
    # E_p_vals_spin = spin_periods()
    # pyplot.plot(E_p_vals_spin[0], E_p_vals_spin[1])

    #plot analytic low-amplitude
    T_o = numpy.pi*2*(7.)**(-.5)
    E_vals = numpy.linspace(-1./3., 0., 20)
    ones = numpy.ones(20)
    pyplot.plot(E_vals, T_o*ones, 'g--',label='Small Amplitude orbital mode')
    pyplot.xlabel('Energy', fontsize=16)
    pyplot.ylabel('Period', fontsize=16)
    pyplot.ylim(0, 12)
    pyplot.xlim(-1./3., 1./6.)
    pyplot.text(
        (-0.3), 8, 'Spinning Mode'
        )
    pyplot.text(
        (.05), 5, 'Orbital Mode'
        )
    # pyplot.text(
    #     x=(-1/3.+.005), y=6.,
    #     text='S', size=15
    #     )
    # pyplot.legend(loc='upper right')
    # pyplot.legend(handles = legend_elements, loc='upper right')
    pyplot.title('Large Amplitude Period vs System Energy', loc='left', fontsize=18)
    pyplot.savefig('plot.png')
    pyplot.show()

def pendu_integrand_gen(tht_m):
    max_term = numpy.sin(tht_m/2)**2
    def pendu_integrand(tht):
        #from rubin chp 5, pdf page 131
        return (
            max_term-numpy.sin(tht/2)**2
            )**(-.5)
    return pendu_integrand

# def pendu_integrand_gen(tht_m):
#     #from https://arxiv.org/pdf/physics/0510206.pdf
#     max_term = numpy.cos(tht_m)
#     def pendu_integrand(tht):
#         #from rubin chp 5, pdf page 131
#         return (
#             numpy.cos(tht)-max_term
#             )**(-.5)
#     return pendu_integrand


def spin_periods():
    samples = 150
    Emax = 2./12.
    T0 = numpy.pi*2*(5./3.)**(-.5)
    E_vals = []
    T_vals = []
    for i in range(1, samples):
        E_val = i*Emax/samples
        tht_max = numpy.arccos(-12.*(E_val-1./12.))
        # print(E_val, tht_max)
        # integrand = pendu_integrand_gen(tht_max)
        # T_val = T0*utils.gaussian_leg(integrand, 10, [0, tht_max])/numpy.pi
        # integral = utils.gaussian_leg(integrand, 16, [0, tht_max])
        # T_val = integral*2**.5*T0/numpy.pi
        T_val = T0*(
            1
            + tht_max**2/16.
            + 11.*tht_max**4/3072.
            # + 9.*25.*tht_max**6/(2.**12*36.)
            )
        E_vals.append((E_val-Emax-2./12.))
        T_vals.append(T_val)
    return [E_vals, T_vals]


if __name__ =='__main__':
    plot_period_vs_energy()
