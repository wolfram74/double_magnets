import sympy

def dot_prod(vec1, vec2):
    total = 0
    for index in range(len(vec1)):
        total += vec1[index]*vec2[index]
    return total

def math_jaxify(string):
    return string.replace('\\', '\\\\')

def ham_equations(state_vec, hamiltonian):
    variables = len(state_vec)/2
    print('force-equivalents')
    for i in range(variables):
        sympy.pprint(state_vec[i+variables])
        force = -hamiltonian.diff(state_vec[i])
        sympy.pprint(-hamiltonian.diff(state_vec[i]))
        print(math_jaxify(sympy.latex(force)))
    print('velocity equivalents')
    for i in range(variables, variables*2):
        sympy.pprint(state_vec[i-variables])
        sympy.pprint(hamiltonian.diff(state_vec[i]))


def cartesian():
    cos, sin = (sympy.cos, sympy.sin)
    x1, x2, y1, y2, x, y = sympy.symbols('x1 x2 y1 y2 x y', real=True)
    vx1, vx2, vy1, vy2 = sympy.symbols('vx1 vx2 vy1 vy2', real=True)
    phi1, phi2 = sympy.symbols('phi1 phi2', real=True)
    omega1, omega2 = sympy.symbols('omega1 omega2', real=True)
    mu, m, I = sympy.symbols('mu m I', real=True, positive=True)
    T = (
        (m/2)*(vx1**2+vy1**2+vx2**2+vy2**2)
        +
        (I/2)*(omega1**2+omega2**2)
        )
    sympy.pprint(T)

def com_polar():
    cos, sin = (sympy.cos, sympy.sin)
    r, tht, phi1, phi2 = sympy.symbols('r theta phi1 phi2', real=True)
    vr, vtht, vphi1, vphi2 = sympy.symbols('r dot_theta dot_phi1 dot_phi2', real=True)
    pr, ptht, pphi1, pphi2 = sympy.symbols('p_r p_theta p_phi1 p_phi2', real=True)
    mt, mr, mu1, mu2, I1, I2, gamma= sympy.symbols('m_t m_r mu_1 mu_2 I_1 I_2 gamma', real=True, positive=True)
    state_variables = [r, tht, phi1, phi2, pr, ptht, pphi1, pphi2]
    Ua = mu1*mu2*(cos(phi1)*cos(phi2)+sin(phi1)*sin(phi2))
    Ub = mu1*(cos(phi1)*cos(tht)+sin(phi1)*sin(tht))
    Uc = mu2*(cos(tht)*cos(phi2)+sin(tht)*sin(phi2))
    Ut = gamma*r**-3*(Ua-3*Ub*Uc)
    sympy.pprint(Ut)
    U_clean = Ut.simplify()
    sympy.pprint(Ut.simplify())
    print(sympy.printing.latex(Ut.simplify()))
    T = (
        pr**2/mr
        +ptht**2/(mr*r**2)
        +pphi1**2/I1
        +pphi2**2/I2
        )/2
    sympy.pprint(T)
    sympy.pprint(U_clean)
    H = T+U_clean
    sympy.pprint(H)
    ham_equations(state_variables, H)
    # L1 = H.diff(tht)
    # L2 = H.diff(phi1)
    # L3 = H.diff(phi2)
    # print('total angular momentum acceleration')
    # sympy.pprint(L1+L2+L3)
    # sympy.pprint((L1+L2+L3).simplify())

def equilibrium_solutions():
    tht, phi1, phi2 = sympy.symbols('theta phi1 phi2', real=True)
    n, m = sympy.symbols('n m', integer=True)
    constraints = sympy.Matrix([
        [-2, 1, 1],
        [0, 1, 1]
        ])
    variables = sympy.Matrix([tht, phi1, phi2])
    out = sympy.Matrix([[n*sympy.pi],[m*sympy.pi]])
    sympy.pprint(constraints)
    sympy.pprint(constraints*variables)
    sympy.pprint(constraints.solve(variables))
    # sympy.pprint(out)
    # sympy.pprint((constraints*variables).solve(out))

def enumerate_equilibria():
    for n in range(6):
        # print [ ((n,m), ((n+m)/2., (n-m)/2.)) for m in range(6)]
        print [ ((n+m)/2., (n-m)/2.) for m in range(6)]
            # print()

if __name__ == '__main__':
    # cartesian()
    # com_polar()
    # equilibrium_solutions()
    enumerate_equilibria()
