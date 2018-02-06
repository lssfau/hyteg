#!/usr/bin/env python3

# This script generates exact solutions for TinyHHG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(K*grad(u)) = f

import math
import numpy as np
from sympy import *
import random

x, y = symbols('x y')


b = 1
h = 1
r_1 = 1.0
r_2 = 2.0
alpha = math.pi / 4.0

def radius(x,y):
  return (r_2 - r_1) * y / h + r_1

def phi(x,y):
  return (-2.0 * alpha) / b * x + alpha + math.pi / 2.0;

F = Matrix([radius(x,y) * cos(phi(x,y)), radius(x,y) * sin(phi(x,y))])
DF = F.jacobian([x,y])

K = simplify(DF * DF.T / DF.det())
u = sin(x)*sinh(y)

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

print(function_template('map_x', ccode(simplify(F[0]))))
print(function_template('map_y', ccode(simplify(F[1]))))
print(function_template('coeff_11', ccode(K[0])))
print(function_template('coeff_12', ccode(K[1])))
print(function_template('coeff_22', ccode(K[3])))
print(function_template('exact', ccode(simplify(u))))
print(function_template('rhs', ccode(simplify(f))))
