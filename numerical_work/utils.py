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
        delQ = stateNp1[i]-stateN[i]
        if delQ == 0:
            times[i-1] = None
            continue
        times[i-1] = (
            stateN[0]
            +delT*(state0[i]-stateN[i])/delQ
            )
    return times

def parse_returned_times(time_vals, tN, tNp1):
    clean_vals = filter(lambda x: not x is None, time_vals)
    for val in clean_vals:
        if not(val>tN and val <tNp1):
            return False
    return sum(clean_vals)/len(clean_vals)

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
        # print(last_state[0],time_guesses, next_state[0])
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
