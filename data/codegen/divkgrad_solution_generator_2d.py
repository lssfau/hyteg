#!/usr/bin/env python3

# This script generates exact solutions for TinyHHG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(k*grad(u)) = f

import numpy as np
from sympy import *

x, y = symbols('x y')

# Coefficient
k = x + 1

# Solution
u = sin(x)*sinh(y)

def gradient(u):
    return np.array([diff(u,x), diff(u,y)])

def divergence(u):
    return diff(u[0],x) + diff(u[1],y)

# Right hand side
f = -divergence(k * gradient(u))


# Generate code
x.name = 'x[0]'
y.name = 'x[1]'

def function_template(name, body):
    return 'std::function<real_t(const hhg::Point3D&)> %s = [](const hhg::Point3D& x) { return %s; };' % (name, body)

print(function_template('coeff', ccode(k)))
print(function_template('exact', ccode(u)))
print(function_template('rhs', ccode(f)))
