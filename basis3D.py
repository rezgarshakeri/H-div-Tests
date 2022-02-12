from sympy import *
import numpy as np

libceed = True

def div(v):
    d  = diff(v[0],x)
    d += diff(v[1],y)
    d += diff(v[2],z)
    return d
    
def curl(v):
    c = []
    c.append(  diff(v[2],y)-diff(v[1],z) )
    c.append(-(diff(v[2],x)-diff(v[0],z)))
    c.append(  diff(v[1],x)-diff(v[0],y) )
    return np.asarray(c)
             
x,y,z = symbols('x y z')
a0,b0,c0,d0,a1,b1,c1,d1,a2,b2,c2,d2 = symbols("a0 b0 c0 d0 a1 b1 c1 d1 a2 b2 c2 d2")
r0,r1,r2,r3,s0,s1,s2,s3,t0,t1,t2,t3 = symbols("r0 r1 r2 r3 s0 s1 s2 s3 t0 t1 t2 t3")

BDDF1 = np.asarray([a0+b0*x+c0*y+d0*z,
                    a1+b1*x+c1*y+d1*z,
                    a2+b2*x+c2*y+d2*z])
BDDF1 += r0*curl([0       ,0       ,x*y*z   ])
BDDF1 += r1*curl([0       ,0       ,x*y**2  ])
BDDF1 += r2*curl([0       ,0       ,x**2*z  ])
BDDF1 += r3*curl([0       ,0       ,x**2*y*z])
BDDF1 += s0*curl([x*y*z   ,0       ,0       ])
BDDF1 += s1*curl([y*z**2  ,0       ,0       ])
BDDF1 += s2*curl([x*y**2  ,0       ,0       ])
BDDF1 += s3*curl([x*y**2*z,0       ,0       ])
BDDF1 += t0*curl([0       ,x*y*z   ,0       ])
BDDF1 += t1*curl([0       ,x**2*z  ,0       ])
BDDF1 += t2*curl([0       ,y*z**2  ,0       ])
BDDF1 += t3*curl([0       ,x*y*z**2,0       ])

"""
     local numbering of Hex
 
     5--------7           z
   / |      / |            |
 6 --|---- 8  |            |
 |   |     |  |            /----y
 |   |     |  |           /
 |   1 -------3         x
 | /       | /
 2 --------4
"""
# normals
nl = [0.0,-1.0, 0.0]
nr = [0.0, 1.0, 0.0]
nbt = [0.0, 0.0, -1.0]
nt = [0.0, 0.0, 1.0]
nf = [1.0, 0.0, 0.0]
nbk = [-1.0, 0.0, 0.0]
# nodes
n1 = [-1., -1., -1.]
n2 = [1., -1., -1]
n4 = [1., 1., -1.]
n3 = [-1., 1., -1.]
n5 = [-1., -1., 1.]
n6 = [1., -1., 1]
n8 = [1., 1., 1.]
n7 = [-1., 1., 1.]

nodes = [n1, n2, n3, n4, n5, n6, n7, n8, n2, n1, n6, n5, n4, n3, n8, n7, n2, n4, n6, n8, n1, n3, n5, n7]
normals = [nbt, nbt, nbt, nbt, nt, nt, nt, nt, nl, nl, nl, nl, nr, nr, nr, nr, nf, nf, nf, nf, nbk, nbk, nbk, nbk]

for i in range(8): # for each vertex
    for j in range(3): # for each direction
        k = 3*i+j

        eqs = []
        for n in range(24):
            eqs.append(np.dot(BDDF1,normals[n]).subs({x:nodes[n][0],y:nodes[n][1],z:nodes[n][2]}))
        
        eqs[k] -= 1 # the k^th functions should be a 1, rest are 0
        sol = solve(eqs)
        ux = BDDF1[0].subs(sol)
        uy = BDDF1[1].subs(sol)
        uz = BDDF1[2].subs(sol)
        
        if libceed:
            def _f(fcn):
                fcn = fcn.replace("x**2","x*x")
                fcn = fcn.replace("y**2","y*y")
                fcn = fcn.replace("z**2","z*z")
                fcn = fcn.replace("x","x[0]")
                fcn = fcn.replace("y","x[1]")
                fcn = fcn.replace("z","x[2]")
                if "/8" in fcn: fcn = "(%s)*0.125;" % (fcn.replace("/8",""))
                if "/16" in fcn: fcn = "(%s)*0.0625;" % (fcn.replace("/16",""))
                return fcn
            print("Bx[%2d] = " % (k) + _f("%s" % (ux)),";")
            print("By[%2d] = " % (k) + _f("%s" % (uy)),";")
            print("Bz[%2d] = " % (k) + _f("%s" % (uz)),";")
