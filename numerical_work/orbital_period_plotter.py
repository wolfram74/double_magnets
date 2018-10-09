import matplotlib.pyplot as pyplot
import utils
import numpy


def parse_file():
    data_in = open('./clean_orbital_periods.txt', 'r')
    output = []
    for line in data_in:
        vals = [float(num) for num in line.rstrip().split(' ')]
        state = vals[:-2]
        period = vals[-2]
        unc = vals[-1]
        L0 = 2*state[4]+state[5]
        T0 = (20*state[4]**2+2*state[5]**2)/2.
        U0 = -(1+3*numpy.cos(state[1]-2*state[2]))/12.
        E0 = T0+U0
        output.append((period, unc, E0, L0, T0, U0))
    return output

def plot_period_vs_energy():
    data = parse_file()
    p_vals = utils.read_column(data, 0)
    sig_vals = utils.read_column(data, 1)
    E_vals = utils.read_column(data, 2)
    pyplot.plot(E_vals, p_vals)
    pyplot.errorbar(E_vals, p_vals, yerr= sig_vals)
    E_p_vals_spin = spin_periods()
    print(E_p_vals_spin[0][:5])
    print(E_p_vals_spin[1][:5])
    pyplot.plot(E_p_vals_spin[0], E_p_vals_spin[1])
    pyplot.show()

def pendu_integrand_gen(tht_m):
    max_term = numpy.sin(tht_m/2)**2
    def pendu_integrand(tht):
        #from rubin chp 5, pdf page 131
        return (
            max_term-numpy.sin(tht/2)**2
            )**(-.5)
    return pendu_integrand

def spin_periods():
    samples = 150
    Emax = 1./12.
    T0 = numpy.pi*2*(5./3.)**(-.5)
    E_vals = []
    T_vals = []
    for i in range(1, samples):
        E_val = i*Emax/samples
        tht_max = numpy.arccos(-12.*E_val)
        integrand = pendu_integrand_gen(tht_max)
        T_val = T0*utils.gaussian_leg(integrand, 10, [0, tht_max])/numpy.pi
        E_vals.append((E_val-Emax))
        T_vals.append(T_val)
    return [E_vals, T_vals]


if __name__ =='__main__':
    plot_period_vs_energy()
