import sympy

def return_T_U_config_velocities():
    phi1, phi2, tht = sympy.symbols('phi_1 phi_2 theta', real=True)
    vphi1, vphi2, vtht = sympy.symbols('v_phi1 v_phi2 v_theta', real=True)
    location = [phi1, phi2, tht]
    velocity = [vphi1, vphi2, vtht]
    T = (vtht**2)/4 + (vphi1**2+vphi2**2)/10
    U = -(
        sympy.cos(phi1 - phi2)
        +3*sympy.cos(phi1 + phi2 - 2*tht)
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

def equi_from_jk(j, k):
    return [
        (j+k)*sympy.pi/2,
        (j-k)*sympy.pi/2,
        0
    ]

def normal_modes():
    T, U, xs, vs = return_T_U_config_velocities()
    w = sympy.symbols('omega')
    L = (T-U)
    pi = sympy.pi
    sympy.pprint(xs)
    sympy.pprint(vs)
    for vel in vs:
        sympy.pprint(L.diff(vel))
    masses = sympy.Matrix([
        [1,0,0],
        [0,1,0],
        [0,0,10]
        ])/10
    equi_curves = [
        [0,0],
        [1,0],
        [0,1],
        [1,1],
    ]
    for curve in equi_curves:
        gamma0 = equi_from_jk(curve[0], curve[1])
        k_vals2 = K_mat_gen(L, xs, gamma0)
        characteristic = k_vals2 - w**2*masses
        normal_modes = characteristic.eigenvects()
        print('mode analysis for ', curve)
        sympy.pprint(gamma0)
        for mode in normal_modes:
            # sympy.pprint(mode)
            # sympy.pprint(mode)
            # sympy.pprint(sympy.simplify(mode[2][0][0]))
            freq_sqr = sympy.solve(mode[0], w**2)[0]
            sympy.pprint(freq_sqr)
            eig_vecs = [vi.subs(w**2, freq_sqr) for vi in mode[2][0]]
            # eig_vecs = [vi for vi in mode[2][0]]
            # sympy.pprint(eig_vecs[0].subs(w**2, freq_sqr))
            sympy.pprint(
                eig_vecs
                )

if __name__ == '__main__':
    normal_modes()
