import numpy

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
        times[i-1] = (
            stateN[0]
            +delT*(state0[i]-stateN[i])/delQ
            )
    return times

