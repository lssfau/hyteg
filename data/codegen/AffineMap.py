#!/usr/bin/env python3

from sympy import *

x, y = symbols('x y')
x1, y1, x2, y2, x3, y3 = symbols('x1 y1 x2 y2 x3 y3')
xb1, yb1, xb2, yb2, xb3, yb3 = symbols('xb1 yb1 xb2 yb2 xb3 yb3')

x2bar = x2 - x1
y2bar = y2 - y1
x3bar = x3 - x1
y3bar = y3 - y1
xb2bar = xb2 - xb1
yb2bar = yb2 - yb1
xb3bar = xb3 - xb1
yb3bar = yb3 - yb1

A = Matrix([[x2bar, x3bar], [y2bar, y3bar]])
Ainv = A**(-1)

Tinv = Ainv * Matrix([[x - x1],[y - y1]])

Ab = Matrix([[xb2bar, xb3bar], [yb2bar, yb3bar]])

F = simplify(Ab * Tinv + Matrix([[xb1], [yb1]]))
DF = simplify(F.jacobian([x,y]))
DFinv = simplify(DF**(-1))
det = DF.det()


x.name = 'x[0]'
y.name = 'x[1]'
x1.name = 'x1_[0]'
y1.name = 'x1_[1]'
x2.name = 'x2_[0]'
y2.name = 'x2_[1]'
x3.name = 'x3_[0]'
y3.name = 'x3_[1]'
xb1.name = 'xb1_[0]'
yb1.name = 'xb1_[1]'
xb2.name = 'xb2_[0]'
yb2.name = 'xb2_[1]'
xb3.name = 'xb3_[0]'
yb3.name = 'xb3_[1]'

# print('real_t xi = {};'.format(ccode(simplify(xi_))))
# print('real_t eta = {};'.format(ccode(simplify(eta_))))

print('Fx[0] =  {};'.format(ccode(F[0])))
print('Fx[1] =  {};'.format(ccode(F[1])))

print('DFx(0,0) =  {};'.format(ccode(DF[(0,0)])))
print('DFx(0,1) =  {};'.format(ccode(DF[(0,1)])))
print('DFx(1,0) =  {};'.format(ccode(DF[(1,0)])))
print('DFx(1,1) =  {};'.format(ccode(DF[(1,1)])))

print('DFxInv(0,0) =  {};'.format(ccode(DFinv[(0,0)])))
print('DFxInv(0,1) =  {};'.format(ccode(DFinv[(0,1)])))
print('DFxInv(1,0) =  {};'.format(ccode(DFinv[(1,0)])))
print('DFxInv(1,1) =  {};'.format(ccode(DFinv[(1,1)])))

print('det = {};'.format(ccode(det)))
