#!/usr/bin/env python3

from sympy import *

xi, eta = symbols('xi eta')
x1, y1, x2, y2, x3, y3 = symbols('x1 y1 x2 y2 x3 y3')
s1, s3bar = symbols('s1 s3bar')
x_c, y_c = symbols('x_c y_c')
radius = symbols('radius')

x2bar = x2 - x1
y2bar = y2 - y1
x3bar = x3 - x1
y3bar = y3 - y1

A = Matrix([[x2bar, x3bar], [y2bar, y3bar]])

def phi(s):
  return radius * cos(s) + x_c;

def psi(s):
  return radius * sin(s) + y_c;


def Phi(eta):
  return (phi(s1 + s3bar * eta) - x1 - x3bar * eta) / (1-eta)

def Psi(eta):
  return (psi(s1 + s3bar * eta) - y1 - y3bar * eta) / (1-eta)

T = A * Matrix([[xi], [eta]]) + Matrix([[x1], [y1]])
x = T[0]
y = T[1]

F = simplify(T + (1 - xi - eta) * Matrix([[Phi(eta)], [Psi(eta)]]))
DF = simplify(F.jacobian([xi,eta]))


xi.name = 'xRef[0]'
eta.name = 'xRef[1]'
x1.name = 'x1_[0]'
y1.name = 'x1_[1]'
x2.name = 'x2_[0]'
y2.name = 'x2_[1]'
x3.name = 'x3_[0]'
y3.name = 'x3_[1]'
s1.name = 's1_'
s3bar.name = 's3bar_'
x_c.name = 'center_[0]'
y_c.name = 'center_[1]'
radius.name = 'radius_'

print('Fx[0] =  {};'.format(ccode(F[0])))
print('Fx[1] =  {};'.format(ccode(F[1])))

print('DFx(0,0) =  {};'.format(ccode(DF[(0,0)])))
print('DFx(0,1) =  {};'.format(ccode(DF[(0,1)])))
print('DFx(1,0) =  {};'.format(ccode(DF[(1,0)])))
print('DFx(1,1) =  {};'.format(ccode(DF[(1,1)])))
