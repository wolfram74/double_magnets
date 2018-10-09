import numpy
import random

def read_column(rank_2_tensor, col_num):
    return [row[col_num] for row in rank_2_tensor]

def rk4_step(gradient_function, state, step_size):
    kernel1 = step_size*gradient_function(state + 0)
    kernel2 = step_size*gradient_function(state + kernel1*.5)
    kernel3 = step_size*gradient_function(state + kernel2*.5)
    kernel4 = step_size*gradient_function(state + kernel3)
    delta = (6.**-1)*(
        kernel1 +
        2.*kernel2 +
        2.*kernel3 +
        kernel4
        )
    return delta

def RK4_adapt(
        gradient_function, state, end_time,
        precision=10.**-4, step_size=10.**-3, max_steps=10**4
        ):
    path = [state]
    running = True
    last_loop = False
    time_left = end_time - state[0]
    print('starting')
    while running:
        last_state = path[-1]
        if last_loop:
            running = False
        if len(path) > max_steps:
            running = False
            print('did not finish, remaining time')
            print(time_left)
            return path
        double_step = last_state + rk4_step(
            gradient_function,last_state, 2.*step_size
            )
        mid_step = last_state + rk4_step(
            gradient_function,last_state, step_size
            )
        two_single_step = mid_step + rk4_step(
            gradient_function, mid_step, step_size
            )
        max_disagreement = max_relative_error(double_step, two_single_step)
        if max_disagreement > precision:
            step_size = step_size/1.25
            continue

        path.append(mid_step)
        path.append(two_single_step)

        if max_disagreement < precision*(10.**-3):
            step_size *= 1.25

        time_left = end_time - path[-1][0]
        if time_left < 2.*step_size and running:
            print('finishing early', len(path)+2)
            step_size = time_left/2
            last_loop = True
    return path

def max_relative_error(vec1, vec2):
    deltas = vec1 - vec2
    return numpy.amax(numpy.absolute(deltas))

def return_time(state0, stateN, stateNp1):
    times = [0.0 for ti in range(len(state0)-1)]
    for i in range(1,len(times)+1):
        delT = stateNp1[0]-stateN[0]
        en, enp1 = stateN[i], stateNp1[i]
        e0 = state0[i]
        pi2 = numpy.pi*2
        if i in [1,2,3]: #cyclic coordinates
            if  en > e0 and enp1 > e0: #possibly gained more than a full rotation
                while en > e0 and enp1 > e0:
                    en  -= pi2
                    enp1-= pi2
            elif en < e0 and enp1 < e0: # possibly lost more than a full rotation
                while en < e0 and enp1 < e0:
                    en  += pi2
                    enp1+= pi2
        delQ = enp1-en
        if delQ == 0:
            times[i-1] = None
            continue
        times[i-1] = (
            stateN[0]
            +delT*(state0[i]-en)/delQ
            )
    return times

def parse_returned_times(time_vals, tN, tNp1):
    clean_vals = filter(lambda x: not x is None, time_vals)
    # print(time_vals)
    # print(clean_vals)
    # print('')
    for val in clean_vals:
        if not(val>tN and val <tNp1):
            return False
    average = sum(clean_vals)/len(clean_vals)
    std_dev = 0
    for val in clean_vals:
        std_dev += (average-val)**2/(len(clean_vals)-1)
    return (average, std_dev**.5)

def return_times_finder(
    gradient_function,state0,
    precision=10.**-4.,step_size=10.**-3,
    max_time=10**2, max_steps=10**4
    ):
    returned_to_state0 = []
    running = True
    last_loop = False
    time_left = max_time - state0[0]
    steps = 0
    last_state = state0
    next_state = state0
    while running:
        if last_loop:
            running = False
        double_step = next_state + rk4_step(
            gradient_function,next_state, 2.*step_size
            )
        mid_step = next_state + rk4_step(
            gradient_function,next_state, step_size
            )
        two_single_step = mid_step + rk4_step(
            gradient_function, mid_step, step_size
            )
        max_disagreement = max_relative_error(double_step, two_single_step)
        if max_disagreement > precision:
            step_size = step_size/1.25
            continue

        last_state = next_state
        next_state = two_single_step
        steps+=2

        time_guesses = return_time(state0, last_state, next_state)
        # per1 = (last_state[0]>2.62 and last_state[0] < 2.63)
        # per2 = (last_state[0]>5.05 and last_state[0] < 5.08)
        # # if last_state[0]>5.05 and last_state[0] < 5.1:
        # if per1 or per2:
        #     print(last_state[0],time_guesses, next_state[0])
        #     print(last_state[:4])
        #     print(last_state[:4]-numpy.pi)
        #     print(state0[:4])
        #     print(next_state[:4])
        #     print(next_state[:4]-numpy.pi)
        #     # print(numpy.cos(last_state[2]-2*last_state[3]))
        #     # print(numpy.cos(state0[2]-2*state0[3]))
        #     # print(numpy.cos(next_state[2]-2*next_state[3]))
        #     print(last_state[2]-2*last_state[3])
        #     print(state0[2]-2*state0[3])
        #     print(next_state[2]-2*next_state[3])
        #     print('')
        linear_results = parse_returned_times(
            time_guesses, last_state[0], next_state[0])
        if linear_results:
            returned_to_state0.append(linear_results)

        if max_disagreement < precision*(10.**-3):
            step_size *= 1.25

        time_left = max_time - next_state[0]
        if time_left < 2.*step_size and running:
            print('finishing early', steps+2)
            step_size = time_left/2
            last_loop = True
        if steps == max_steps:
            running = False
            print('did not finish, remaining time')
            print(time_left)
            break

    return returned_to_state0

def fit_slope_and_chisqr(data):
    denom = 0
    numer = 0
    print('in fitter')
    print(len(data))
    if len(data) == 0:
        return -1, -1
    for i in range(len(data)):
        # print(data[i])
        denom += (i+1)**2
        numer += data[i][0]*(i+1)
    T_fit = float(numer)/float(denom)
    chi_sqr = 0
    for i in range(len(data)):
        chi_sqr += (T_fit*(i+1)-data[i][0])**2/data[i][0]**2
    return T_fit, chi_sqr

def gaussian_leg(func, samples, interval):
    '''
    unlike other methods doesn't use evenly spaced inputs
    for n evaluation points, you would use the n solutions for the nth order legendre polynomial
    '''
    poly_coeffs = numpy.zeros(samples+1)
    poly_coeffs[-1]=1
    legendre = numpy.polynomial.legendre
    roots = legendre.legroots(poly_coeffs)
    poly_der_coeffs = legendre.legder(poly_coeffs)
    poly_der_vals = legendre.legval(roots, poly_der_coeffs)
    weights = []
    for index in range(len(roots)):
        weights.append(
            2.0/((1-roots[index]**2)*(poly_der_vals[index]**2))
            )
    slope = (interval[1]-interval[0])/2.0
    const = (interval[1]+interval[0])/2.0
    sample_points = map(lambda y: slope*y+const, roots)
    total = 0
    for index in range(len(sample_points)):
        total += func(sample_points[index])*slope*weights[index]
    return total

def reduced_state_gen():
    # t, phd, pht, tht, pd, pt, pth
    state_outp = [0.,0.,0.,0.,0.,0.,0.]
    state_outm = [0.,0.,0.,0.,0.,0.,0.]
    bad = 0
    good = 0
    while True:
        pht = (random.random()-.5)*numpy.pi
        tht = (random.random()-.5)*numpy.pi/2.+pht/2.
        U_0 = -(1+3*numpy.cos(pht-2*tht))/12.
        T_0 = -random.random()*U_0
        E_0 = T_0+U_0
        # print(U_0, T_0, E_0)
        L_0 = (-random.random()*4.5*E_0)**.5
        # print(L_0)
        radical = 56*T_0 - 40*L_0**2
        if radical < 0:
            continue
        ptp = (4.*L_0+(56.*T_0-40*L_0**2)**.5)/28.
        ptm = (4.*L_0-(56.*T_0-40*L_0**2)**.5)/28.
        pthp = L_0-2*ptp
        pthm = L_0-2*ptm
        Tc = pthp**2 + 10*ptp**2
        state_outp = [0.,0.,pht,tht,0.,ptp,pthp]
        state_outm = [0.,0.,pht,tht,0.,ptm,pthm]
        break
    return state_outp, state_outm
'''
2.62786

5.076827


'''
