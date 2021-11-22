#!/usr/bin/env python3

# This script generates exact solutions for HyTeG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#     -div(k*grad(u)) = f

import numpy as np
from sympy import *
import random

x, y = symbols('x y')

def φ(x,y):
    return atan2(y,x)

def r(x,y):
    return sqrt(x*x + y*y)

def sigmoidK():
    α = symbols('alpha')
    k = 1 + 1/(1+exp(-(α*(x-0.5))))
    u = log(2*exp(α*(x-0.5)) + 1)/α - 2*x
    return (k,u)

def annulus():
    u = sin(4*φ(x,y))*sin(2*r(x,y))
    return (1,u)

def smooth_jump():
    α,φ = symbols('alpha phi')
    k = 2+tanh(α*(x-0.5));
    u = sin(φ*pi*x)*sinh(pi*y)
    return (k,u)


def basic_example():
    u = sin(pi*x)*sin(pi*y)
    k = x*x + 1

    return (k,u)

def poly_5():
    random.seed(5618947216424)
    k = 1 + random.random()*x + random.random()*y + random.random()*x**2 + random.random() * x * y + random.random() * y**2
    k = k + random.random()*x**3 + random.random() * x**2 * y + random.random() * x * y**2 + random.random()*y**3
    k = k + random.random()*x**4 + random.random() * x**3 * y + random.random() * x**2 * y**2 + random.random()* x * y**3 + random.random() * y**4
    k = k + random.random()*x**5 + random.random() * x**4 * y + random.random() * x**3 * y**2 + random.random() * x**2 * y**3 + random.random() * x * y**4 + random.random() * y**5
    u = x**3 * y**2 / (x*y + 1)

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


k, u = sigmoidK()

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
    return 'std::function<real_t(const hyteg::Point3D&)> %s = [](const hyteg::Point3D& x) { return %s; };' % (name, body)

print(function_template('coeff', ccode(simplify(k))))
print(function_template('exact', ccode(simplify(u))))
print(function_template('rhs', ccode(simplify(f))))
