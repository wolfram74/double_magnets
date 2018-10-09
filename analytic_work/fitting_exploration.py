import sympy

def fits():
    e0, enm1, en, enp1 = sympy.symbols('e0 en1 en2 en3')
    sympy.pprint(e0+enm1+en+enp1)

    t0, tnm1, tn, tnp1, t = sympy.symbols('t0 tn1 tn2 tn3 t')
    sympy.pprint(t0+tnm1+tn+tnp1)
    curve = (
        enm1*(t-tn)*(t-tnp1)/((tnm1-tn)*(tnm1-tnp1))
        +en*(t-tnm1)*(t-tnp1)/((tn-tnm1)*(tn-tnp1))
        +enp1*(t-tnm1)*(t-tn)/((tnp1-tnm1)*(tnp1-tn))
        )
    sympy.pprint(curve)
    sympy.pprint(curve.expand())
    sympy.pprint(curve.expand().simplify())

if __name__=='__main__':
    fits()
