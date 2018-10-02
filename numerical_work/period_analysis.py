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

def sort_by_index(data, index):
    def interested_element(arr):
        return arr[index]
    return sorted(data ,
        key=interested_element
        )

def pretty_print(data):
    for line in data:
        print(line)

def ordering():
    data_points = parse_file()
    # print(data_points[0])
    #T, E, L, KE, PE
    ordered_data = sort_by_index(data_points, 2)
    # print(data_points[0])
    pretty_print(ordered_data)
    # print(data_points[0])
if __name__ == '__main__':
    ordering()
