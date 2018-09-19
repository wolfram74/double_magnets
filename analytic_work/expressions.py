import sympy

pphi1, pphi2, ptht, pr = sympy.symbols('p_phi1 p_phi2 p_theta p_r')
vphi1, vphi2, vtht, vr = sympy.symbols('v_phi1 v_phi2 v_theta v_r')
phi1, phi2, tht, r = sympy.symbols('phi_1 phi_2 theta r')

positions = [phi1, phi2, tht, r]
velcoties = [vphi1, vphi2, vtht, vr]
momenta = [pphi1, pphi2, ptht, pr]

Tp = (2*pr**2+2*ptht**2*r**(-2) + 10*pphi1**2 + 10*pphi2**2)/2
Tl = (vr**2/2+vtht**2*r**(-2)/2 + vphi1**2/2 + vphi2**2/10)/2
U = -r**(-3)*(
    sympy.cos(phi1-phi2)
    +3*sympy.cos(phi1+phi2-2*tht)
    )/12
Ham = Tp+U
Lag = Tl-U
sympy.pprint(Tp)
sympy.pprint(U)
sympy.pprint(Ham)
