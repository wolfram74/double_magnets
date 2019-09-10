import numpy
import sympy
cos = sympy.cos
sin = sympy.sin
mu, mu0, pi = sympy.symbols('mu mu0 pi')
phi1, phi2, tht, r = sympy.symbols('phi_1 phi_2 theta r')

mu1hat = mu*numpy.array([cos(phi1), sin(phi1)])
mu2hat = mu*numpy.array([cos(phi2), sin(phi2)])
rhat = numpy.array([cos(tht), sin(tht)])
U1 = (mu0/(4*pi*r**3))*(
    numpy.dot(mu1hat, mu2hat)
    - 3*numpy.dot(mu2hat,rhat)*numpy.dot(mu1hat,rhat)
    )
sympy.pprint(numpy.dot(mu1hat, mu2hat))
sympy.pprint(U1)
sympy.pprint(U1.expand().collect(mu**2*mu0))
simp_U1 = U1.expand().simplify()
sympy.pprint(simp_U1)
sympy.pprint(simp_U1*2)
