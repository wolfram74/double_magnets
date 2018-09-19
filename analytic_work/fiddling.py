import sympy
from expressions import *
''' includes
positions, velcoties, momenta
Tp,Tl,U,
Ham, Lag
'''

def ham_equations():
    for qi in positions:
        sympy.pprint(qi)
        sympy.pprint(-Ham.diff(qi))
        print(sympy.latex(-Ham.diff(qi).subs(positions[3], 1)))

if __name__ == '__main__':
    ham_equations()
