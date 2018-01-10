#!/usr/bin/env python3

# This script generates exact solutions for TinyHHG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(k*grad(u)) = f

import numpy as np
from sympy import *

x, y = symbols('x y')

def basic_example():
    k = 1 + 3*x + 4*y + 7*x**2
    u = sin(x)*sinh(y)

    return (k,u)

def plume_example():

    alpha = 30
    beta = 10

    R = sqrt((x-0.5)**2)

    trunk = exp(-10*beta*R**2)
    top = exp(-beta*R**2 - alpha*(y-1.5*(R+0.5)*(1.0-R-0.5)-0.3)**2)

    def sigmf(var, a, c):
        return 1/(1+exp(-a*(var-c)))

    k = sigmf(1-y, 17, 0.465) * trunk * 0.98 + 0.89*top + 1

    u = sin(x)*sinh(y)

    return (k,u)


k, u = basic_example()

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

print(function_template('coeff', ccode(simplify(k))))
print(function_template('exact', ccode(simplify(u))))
print(function_template('rhs', ccode(simplify(f))))
