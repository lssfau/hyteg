#!/usr/bin/env python3

# This script generates exact solutions for HyTeG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(2*eps(u) - p*Id) = f
#                    div(u) = 0

import math
import numpy as np
from sympy import *
import random

x, y = symbols('x y')

p = x**2 * y**3

Id = Matrix([[1,0], [0,1]])

# Divergence of the curl is always zero
U = sin(x)*cos(y) / (x+1)

def divergence(u,v):
    return diff(u,x) + diff(v,y)

def constructDivergenceFreeSolution(U):
    
    U_x = simplify(diff(U,x))
    V = integrate(-U_x, y)

    print("divergence(U,V) = {}".format(simplify(divergence(U,V))))

    return (U,V)

u, v = constructDivergenceFreeSolution(U)

def gradient(u, v):
    return Matrix([[diff(u,x), diff(u,y)], [diff(v,x), diff(v,y)]])

def epsilon(u, v):
    return 0.5 * (gradient(u,v) + gradient(u,v).T)

def divergence(u,v):
    return diff(u,x) + diff(v,y)

def divergence_mat(eps):
    return Matrix([[diff(eps[(0,0)],x) + diff(eps[(0,1)],y)], [diff(eps[(1,0)],x) + diff(eps[(1,1)],y)]])

# Right hand side
f = -divergence_mat(2 * epsilon(u,v) - p * Id)

# Generate code
x.name = 'x[0]'
y.name = 'x[1]'

def function_template(name, body):
    return 'std::function<real_t(const hyteg::Point3D&)> %s = [](const hyteg::Point3D& x) { return %s; };' % (name, body)

print(function_template('exact_u', ccode(simplify(u))))
print(function_template('exact_v', ccode(simplify(v))))
print(function_template('exact_p', ccode(simplify(p))))

print(function_template('rhs_u', ccode(simplify(f[0]))))
print(function_template('rhs_v', ccode(simplify(f[1]))))
