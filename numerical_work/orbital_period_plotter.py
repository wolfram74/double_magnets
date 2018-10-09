import matplotlib.pyplot as pyplot
import utils
import numpy

def parse_file():
    data_in = open('./1538506548.txt', 'r')
    output = []
    for line in data_in:
        vals = [float(num) for num in line.rstrip().split(' ')]
        state = vals[:-1]
        period = vals[-1]
        L0 = 2*state[4]+state[5]
        T0 = (20*state[4]**2+2*state[5]**2)/2.
        U0 = -(1+3*numpy.cos(state[1]-2*state[2]))/12.
        E0 = T0+U0
        output.append((period, E0, L0, T0, U0))
    return output
