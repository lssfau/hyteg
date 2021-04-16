#!/usr/bin/env python3

# This script generates exact solutions for HyTeG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(k*grad(u)) = f

import numpy as np
from sympy import *
import random

#------- coordinates ---------#
x, y, z = symbols('x y z')

def torusCoords(R, x, y, z):
    ψ = atan2(y, x)  # angle anlong torus
    # polar coordinates (φ,r) on crosssection
    φ = atan2(z, sqrt(x*x + y*y) - R)
    r = sqrt(R*R + x*x + y*y + z*z - 2*R*sqrt(x*x + y*y))

    return (ψ,φ,r)


#------ sample solutions -------#

def tokamak_jump():
    R_torus = symbols('R_torus')
    r_torus = symbols('r_torus')
    ψ,φ,r = torusCoords(R_torus, x,y,z)
    # parameters to adjust position, width and limits of jump
    d_jump,r_jump,k_min,k_max = symbols('d_jump r_jump k_min k_max')

    u = cos(pi/2*r/r_torus)*sin(pi*z/r_torus)
    k = k_min + 0.5 * (k_max - k_min) * (1 + tanh(7 / d_jump * (r - r_jump)))

    return u,k


def example():
    u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    k = 1
    return (u,k)


#------------ functions to compute rhs -------------#

def gradient(u):
    return np.array([diff(u, x), diff(u, y), diff(u, z)])

def divergence(u):
    return diff(u[0],x) + diff(u[1],y) + diff(u[2],z)



if __name__ == '__main__':

    # exact solution and diffusion coefficient
    u,k = tokamak_jump()

    print('u =', u)
    u = simplify(u)
    print(' = ', u)

    print('k =', k)
    k = simplify(k)
    print(' = ', k)

    # Right hand side
    f = -divergence(k * gradient(u))

    print('∇⋅(k∇u) = f =', f)
    f = simplify(f)
    print(' = ', f)

    print('=========================================')



    # Generate code
    x.name = 'x[0]'
    y.name = 'x[1]'
    z.name = 'x[2]'

    def function_template(name, body):
        return 'std::function<real_t(const hyteg::Point3D&)> %s = [](const hyteg::Point3D& x) { return %s; };' % (name, body)

    print(function_template('coeff', ccode(simplify(k))))
    print(function_template('exact', ccode(simplify(u))))
    print(function_template('rhs', ccode(simplify(f))))

    # tokamak jump c++
    #
    # parameters for analytic solution and diffusion coefficient:
    # double R_torus;       // torus radius
    # double r_torus;       // radius of torus crosssection
    # double r_jump;        // location of the jump (r_jump in [0,r_torus])
    # double k_min, k_max;  // min/max value of coefficient
    # double d_jump;        // smoothed width of jump, i.e. |coeff(r) - k_min| ≈ 0 for r < (r_jump - d_jump/2) and |coeff(r) - k_max| ≈ 0 for r > (r_jump + d_jump/2)