from sympy import *
import numpy as np
from sympy.utilities.autowrap import ufuncify
import matplotlib.pyplot as plt

plot = False
libceed = True

# From Wheeler2009, the BDM1 spaces span
#  a1 x + b1 y + g1 + rx**2 + 2sxy
#  a2 x + b2 y + g2 - 2rxy - sy**2

# spatial coordinates in the reference element x,y \in [0,+1]
x,y = symbols('x y')

# coefficients of the space
a1,a2,b1,b2,g1,g2,r,s = symbols('a1 a2 b1 b2 g1 g2 r s')
vx  = a1*x + b1*y + g1 + r*x**2 + 2*s*x*y
vy  = a2*x + b2*y + g2 - 2*r*x*y - s*y**2
"""
3---4
|   |
1---2
"""
nb = [0,-1.]
nt = [0,1.]
nl = [-1.,0]
nr = [1.,0]
normals = [nb, nb, nr, nr, nt, nt, nl, nl]
n1 = [-1.,-1.]
n2 = [1.,-1.]
n3 = [-1.,1.]
n4 = [1.,1.]
nodes = [n1, n2, n2, n4, n3, n4, n1, n3]

if plot: fig,ax = plt.subplots(nrows=2,ncols=4)
for i in range(4):
    for j in range(2):
        k = 2*i+j
        eqs = []
        for n in range(8):
            eqs.append(np.dot([vx,vy],normals[n]).subs({x:nodes[n][0],y:nodes[n][1]}))

        eqs[k] -= 1
        sol = solve(eqs)
        ux = vx.subs(sol)
        uy = vy.subs(sol)
        
        if libceed:
            def _f(fcn):
                fcn = fcn.replace("x**2","x*x")
                fcn = fcn.replace("y**2","y*y")
                fcn = fcn.replace("x","x[0]")
                fcn = fcn.replace("y","x[1]")
                if "/4" in fcn: fcn = "(%s)*0.25;" % (fcn.replace("/4",""))
                if "/8" in fcn: fcn = "(%s)*0.125;" % (fcn.replace("/8",""))
                return fcn
            print("Bx[%d] = " % (k) + _f("%s" % (ux)),";")
            print("By[%d] = " % (k) + _f("%s" % (uy)),";")
            
        if plot:
            X, Y = np.meshgrid(np.linspace(-1,1,9), np.linspace(-1,1,9))
            uxy = ufuncify((x, y), ux)
            vxy = ufuncify((x, y), uy)
            ax[j,i].quiver(X, Y, uxy(X, Y), vxy(X, Y))
            ax[j,i].set_title("v%d%d = [%s, %s]" % (i+1,j+1,ux,uy))
            ax[j,i].plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],'-k')
            ax[j,i].set_xlim([-1.5,1.5])
            ax[j,i].set_ylim([-1.5,1.5])
            print("b%d%d = [%s, %s]" % (i+1,j+1,ux,uy))
if plot: plt.show()



