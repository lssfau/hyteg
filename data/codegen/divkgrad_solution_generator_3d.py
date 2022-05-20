#!/usr/bin/env python3

# This script generates exact solutions for HyTeG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(k*grad(u)) = f

import numpy as np
from sympy import *
import random

#------- coordinates ---------#
x, y, z = symbols('x y z')

def torusCoords(R0,R1,R2, x, y, z):
    φ = atan2(y, x)  # angle anlong torus

    # helper variables
    ξ = (sqrt(x*x + y*y) - R0)/R1
    ζ = z/R2
    # normalized polar coordinates (θ,ϱ) on crosssection
    θ = atan2(ζ,ξ)
    ϱ = sqrt(ξ*ξ + ζ*ζ)
    # ϱ = sqrt((R0**2 + x**2 + y**2 - 2*R0*sqrt(x**2 + y**2))/(R1**2) + (z**2)/(R2**2))

    return (φ,θ,ϱ)

def φ(x, y):
    return atan2(y, x)


def r(x, y, z):
    return sqrt(x*x + y*y + z*z)

#------ sample solutions -------#

def tokamak_jump():
    # tokamak parameters
    R0 = symbols('R0')
    R1 = symbols('R1')
    R2 = symbols('R2')
    δ  = symbols('delta')

    # parameters to adjust position, width and limits of jump
    d_jump,r_jump,k_min,k_max = symbols('d_jump r_jump k_min k_max')

    # helper variables
    r = sqrt(x**2 + y**2) - R0
    rR = r/R1
    zR = z/R2

    φ, θ, ϱ = torusCoords(R0, R1, R2, x, y, z)

    # boundary points r_min,r_max s.th. (x,y,z) with sqrt(x^2 + y^2) - R0 ∈ {r_min, r_max} are on ∂Ω
    r_min = -cos(asin(zR) - asin(δ)*zR)
    r_max = cos(asin(zR) + asin(δ)*zR)

    u = sin(rR-δ)*(rR - r_min)*(rR - r_max) * zR*zR*sin(pi*zR)
    # u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    # u = (rR-δ)*sin(pi*(rR - r_min)/(r_max-r_min)) * sin(pi*zR)
    # u = (rR-δ)*cos(pi/2*ϱ)*sin(pi*zR)
    # u = (rR-δ)*cos(pi/2*ϱ)*sin(pi*zR**3)
    k = k_min + 0.5 * (k_max - k_min) * (1 + tanh(3.5 / d_jump * (ϱ - r_jump)))

    print("// parameters for analytic solution and diffusion coefficient:")
    print("double %s;\t\t// torus radius" %R0);
    print("double %s;\t\t// semi-minor half axis of crossection" %R1)
    print("double %s;\t\t// semi-major half axis of crossection" %R2)
    print("double %s;\t\t// triangularity parameter" %δ)
    print("double %s;\t\t// triangularity parameter" %δ)
    print("double %s;\t\t// location of the jump, relative to radius of crosssection" %r_jump)
    print("double %s;\t\t// smoothed width of the jump, relative to diameter of crosssection" %d_jump)
    print("double %s, %s;\t// min/max value of jump coefficient" %(k_min,k_max))
    print("\n// TODO initialize above variables!\n")

    return u,k,ϱ*ϱ


def spherical_shell_example():
    u = sin(2*r(x, y, z))
    return (u, 1)


def spherical_shell_example2():
    u = sin(1/r(x,y,z)**2)
    return (u, 1)


def example():
    u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    k = x*x + 1
    return (u,k)


#------------ functions to compute rhs -------------#

def gradient(u):
    return np.array([diff(u, x), diff(u, y), diff(u, z)])

def divergence(u):
    return diff(u[0],x) + diff(u[1],y) + diff(u[2],z)


def function_template(name, body):
    return 'std::function<real_t(const hyteg::Point3D&)> %s = [=](const hyteg::Point3D& x) { return %s; };' % (name, body)

def to_cpp_function(name, expr):
    cs = cse(simplify(expr))
    # print(cs)
    print('std::function<real_t(const hyteg::Point3D&)> %s = [=](const hyteg::Point3D& x) \n{'%name)
    for (var,val) in cs[0]:
        print('   auto %s = %s;'%(ccode(var),ccode(val)))
    print('   return %s;\n};'%ccode(cs[1][0]))


# ---------- tokamak ---------------
def tokamak_solution(foo):

    # exact solution and diffusion coefficient
    u,k,ϱ2 = foo()

    # Right hand side
    f = -divergence(k * gradient(u))

    r2 = symbols('r2_xyz')
    r2.name = 'r2(x,y,z)'

    u = u.subs(ϱ2,r2)
    k = k.subs(ϱ2,r2)
    f = f.subs(ϱ2,r2)

    # print('r2 = @(x,y,z) max(1e-100,', octave_code(ϱ2), ');')
    # print('u = @(x,y,z)', octave_code(u), ';')
    # print('k = @(x,y,z)', octave_code(k), ';')
    # print('f = @(x,y,z)', octave_code(f), ';')


    # Generate code
    x.name = 'x[0]'
    y.name = 'x[1]'
    z.name = 'x[2]'
    r2.name = 'r2(x)'

    print(function_template('r2', 'std::max(1e-100, %s)' %ccode(ϱ2)))
    # to_cpp_function('r2', max(1e-100,ϱ2))
    to_cpp_function('coeff', k)
    to_cpp_function('exact', u)
    to_cpp_function('rhs', f)
    # print(function_template('r2', 'std::max(1e-100,  %s)' %ccode(ϱ2)))
    # print(function_template('coeff', ccode(k)))
    # print(function_template('exact', ccode(u)))
    # print(function_template('rhs', ccode(f)))


def generic_solution(foo):
    # exact solution and diffusion coefficient
    u,k = foo()

    # Right hand side
    f = -divergence(k * gradient(u))

    # r2 = symbols('r2_xyz')
    # r2.name = 'r2(x,y,z)'

    # u = u.subs(ϱ2,r2)
    # k = k.subs(ϱ2,r2)
    # f = f.subs(ϱ2,r2)

    # print('r2 = @(x,y,z) max(1e-100,', octave_code(ϱ2), ');')
    # print('u = @(x,y,z)', octave_code(u), ';')
    # print('k = @(x,y,z)', octave_code(k), ';')
    # print('f = @(x,y,z)', octave_code(f), ';')


    # Generate code
    x.name = 'x[0]'
    y.name = 'x[1]'
    z.name = 'x[2]'
    # r2.name = 'r2(x)'

    # print(function_template('r2', 'std::max(1e-100, %s)' %ccode(ϱ2)))
    # to_cpp_function('r2', max(1e-100,ϱ2))
    to_cpp_function('coeff', k)
    to_cpp_function('exact', u)
    to_cpp_function('rhs', f)


if __name__ == '__main__':
    # tokamak_solution(tokamak_jump)
    generic_solution(spherical_shell_example2)
    # generic_solution(example)
