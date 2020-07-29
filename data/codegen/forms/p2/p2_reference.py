import numpy as np
from sympy import *

def create_2d():
  # Reference triangle
  reference = np.array([[0,0], [1,0], [0,1], [0.5,0.5], [0,0.5], [0.5,0]])

  # Symbolic coordinate
  eta, xi = symbols('eta xi')

  # Construct local nodal basis functions
  basis = [1, eta, xi, eta*xi, eta**2, xi**2]
  phi = []

  for i in range(len(basis)):
    a, b, c, d, e, f = symbols('a b c d e f')
    phi_l = a*basis[0] + b*basis[1] + c*basis[2] + d*basis[3] + e*basis[4] + f*basis[5]
    system = []
    for j in range(6):
        phi_l_t = phi_l.subs([(eta,reference[j][0]), (xi,reference[j][1])])
        system.append(Eq(phi_l_t, float(i==j)))

    sol = solve(system, [a,b,c,d,e,f])
    phi.append(simplify(phi_l.subs(sol)))

  # Construct local nodal gradient basis functions
  def gradient(f):
    return simplify(Matrix([[diff(f,eta)], [diff(f,xi)]]))

  grad_phi = []

  for p in phi:
    grad_phi.append(gradient(p))

  return phi, grad_phi

def create_3d():
  # Reference tet
  reference = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0]])

  # Symbolic coordinate
  eta, xi, nu = symbols('eta xi nu')

  # Construct local nodal basis functions
  basis = [1, eta, xi, nu, eta*xi, eta*nu, xi*nu, eta**2, xi**2, nu**2]
  phi = []

  for i in range(len(basis)):
    a, b, c, d, e, f, g, h, k, l = symbols('a b c d e f g h k l')
    phi_l = a*basis[0] + b*basis[1] + c*basis[2] + d*basis[3] + e*basis[4] + f*basis[5]+ g*basis[6] + h*basis[7] + k*basis[8] + l*basis[9]
    system = []
    for j in range(10):
        phi_l_t = phi_l.subs([(eta,reference[j][0]), (xi,reference[j][1]), (nu,reference[j][2])])
        system.append(Eq(phi_l_t, float(i==j)))

    sol = solve(system, [a,b,c,d,e,f,g,h,k,l])
    phi.append(simplify(phi_l.subs(sol)))

  # Construct local nodal gradient basis functions
  def gradient(f):
    return simplify(Matrix([[diff(f,eta)], [diff(f,xi)], [diff(f,nu)]]))

  grad_phi = []

  for p in phi:
    grad_phi.append(gradient(p))

  return phi, grad_phi
