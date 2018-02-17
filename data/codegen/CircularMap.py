#!/usr/bin/env python3

from sympy import *

x, y = symbols('x y')
x1, y1, x2bar, y2bar, x3bar, y3bar = symbols('x1 y1 x2bar y2bar x3bar y3bar')
s1, s3bar = symbols('s1 s3bar')
x_c, y_c = symbols('x_c y_c')
radius = symbols('radius')

A = Matrix([[x2bar, x3bar], [y2bar, y3bar]])
Ainv = A**(-1)

Tinv = Ainv * Matrix([[x - x1],[y - y1]])

xi = Tinv[0]
eta = Tinv[1]

#xi, eta = symbols('xi eta')

def phi(s):
  return radius * cos(s) + x_c;

def psi(s):
  return radius * sin(s) + y_c;


def Phi(eta):
  return (phi(s1 + s3bar * eta) - x1 - x3bar * eta) / (1-eta)

def Psi(eta):
  return (psi(s1 + s3bar * eta) - y1 - y3bar * eta) / (1-eta)


F = A * Tinv + Matrix([[x1], [y1]]) + (1 - xi - eta) * Matrix([[Phi(eta)], [Psi(eta)]])
DF = F.jacobian([x,y])

det = DF.det()


x.name = 'x[0]'
y.name = 'x[1]'
x1.name = 'x1_[0]'
y1.name = 'x1_[1]'
x2bar.name = 'x2bar_[0]'
y2bar.name = 'x2bar_[1]'
x3bar.name = 'x3bar_[0]'
y3bar.name = 'x3bar_[1]'
s1.name = 's1_'
s3bar.name = 's3bar_'
x_c.name = 'center_[0]'
y_c.name = 'center_[1]'
radius.name = 'radius_'

# print('real_t xi = {};'.format(ccode(simplify(xi_))))
# print('real_t eta = {};'.format(ccode(simplify(eta_))))

print('Fx[0] =  {};'.format(ccode(F[0])))
print('Fx[1] =  {};'.format(ccode(F[1])))

print('DFx(0,0) =  {};'.format(ccode(DF[(0,0)])))
print('DFx(0,1) =  {};'.format(ccode(DF[(0,1)])))
print('DFx(1,0) =  {};'.format(ccode(DF[(1,0)])))
print('DFx(1,1) =  {};'.format(ccode(DF[(1,1)])))

# print('DFxInv(0,0) =  {};'.format(ccode(DFinv[(0,0)])))
# print('DFxInv(0,1) =  {};'.format(ccode(DFinv[(0,1)])))
# print('DFxInv(1,0) =  {};'.format(ccode(DFinv[(1,0)])))
# print('DFxInv(1,1) =  {};'.format(ccode(DFinv[(1,1)])))

print('det = {};'.format(ccode(det)))
