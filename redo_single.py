import sympy

def return_T_U_config_velocities():
    r, phi, tht = sympy.symbols('r phi theta', real=True)
    vr, vphi, vtht = sympy.symbols('v_r v_phi v_theta', real=True)
    gamma, m, mu, I = sympy.symbols('gamma m mu I', real=True, positive=True)
    # location = [phi, tht, r]
    # velocity = [vphi, vtht, vr]
    # T = (vr**2+r**2*vtht**2+vphi**2/10)/2
    # U = -r**(-3)*(
    #     sympy.cos(phi)
    #     +3*sympy.cos(phi-2*tht)
    #     )/12
    location = [phi, tht]
    velocity = [vphi, vtht]
    T = (vtht**2+vphi**2/10)/2
    U = -(
        sympy.cos(phi)
        +3*sympy.cos(phi-2*tht)
        )/12

    return [T, U, location, velocity]

def K_mat_gen(L, q_vars, equi_vals):
    subs = []
    for index in range(len(q_vars)):
        subs.append((
            q_vars[index],
            equi_vals[index]
            ))
    return sympy.Matrix(
        [[
            -L.diff(row_ind, col_ind).subs(subs) for col_ind in q_vars
        ] for row_ind in q_vars]
        )

def normal_modes():
    T, U, xs, vs = return_T_U_config_velocities()
    w = sympy.symbols('omega')
    L = (T-U)
    pi = sympy.pi
    for vel in vs:
        sympy.pprint(L.diff(vel))
    masses = sympy.Matrix([
        [1,0],
        [0, 10]
        ])/10
    k_vals2 = K_mat_gen(L, xs, [0, 0])
    # sympy.pprint(masses)
    # sympy.pprint(k_vals2)
    # sympy.pprint(characteristic)
    characteristic = k_vals2 - w**2*masses
    normal_modes = characteristic.eigenvects()
    for mode in normal_modes:
        # sympy.pprint(mode)
        # sympy.pprint(mode)
        # sympy.pprint(sympy.simplify(mode[2]))
        # sympy.pprint(sympy.simplify(mode[2][0][0]))
        sympy.pprint(sympy.solve(mode[0], w**2))
    k_vals2 = K_mat_gen(L, xs, [pi, -pi/2])
    characteristic = k_vals2 - w**2*masses
    normal_modes = characteristic.eigenvects()
    for mode in normal_modes:
        # sympy.pprint(mode)
        # sympy.pprint(mode)
        # sympy.pprint(sympy.simplify(mode[2]))
        # sympy.pprint(sympy.simplify(mode[2][0][0]))
        sympy.pprint(sympy.solve(mode[0], w**2))

if __name__ == '__main__':
    normal_modes()
