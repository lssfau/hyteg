#!/usr/bin/env python3

# This script generates exact solutions for HyTeG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(k*grad(u)) = f

import numpy as np
from sympy import *
import random

x, y = symbols('x y')

def basic_example():
    k_11 = 1 + 2*x + 3*y + x**2
    k_12 = 0
    k_22 = 1 + 2*x + 3*y + x**2

    K = Matrix([[k_11, k_12], [k_12, k_22]])

    u = sin(x)*sinh(y)

    return (K,u)

K, u = basic_example()

def gradient(u):
    return Matrix([diff(u,x), diff(u,y)])

def divergence(u):
    return diff(u[0],x) + diff(u[1],y)

# Right hand side
f = -divergence(K * gradient(u))


# Generate code
x.name = 'x[0]'
y.name = 'x[1]'

def function_template(name, body):
    return 'std::function<real_t(const hhg::Point3D&)> %s = [](const hhg::Point3D& x) { return %s; };' % (name, body)

print(function_template('coeff_11', ccode(simplify(K[0]))))
print(function_template('coeff_12', ccode(simplify(K[1]))))
print(function_template('coeff_22', ccode(simplify(K[3]))))
print(function_template('exact', ccode(simplify(u))))
print(function_template('rhs', ccode(simplify(f))))
