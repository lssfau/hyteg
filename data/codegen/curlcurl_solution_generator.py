# This script generates exact solutions for HyTeG in order to analyze convergence of errors
# The generated solutions are for following PDE:
#   curl curl u + u = f  in Ω
#             u × n = 0  on ∂Ω

import sympy as sp
from sympy.abc import x, y, z, r, R

def curl(v):
  return sp.Matrix([ v[2].diff(y) - v[1].diff(z)
                   , v[0].diff(z) - v[2].diff(x)
                   , v[1].diff(x) - v[0].diff(y)
                   ])

def solid_torus():
  # f and a can be choosen freely
  f = sp.Matrix([sp.sin(x), sp.cos(z), sp.sin(y)])
  a = 8 * sp.cos(x) * sp.sin(z) * sp.sin(y)

  rXY = sp.sqrt(x*x + y*y)
  u = f * ((rXY - R)**2 + z**2 - r**2) + a / r * (sp.Matrix([x, y, z]) - R / rXY * sp.Matrix([x, y, 0]))

  rhs = curl(curl(u)) + u

  return (u, rhs)

def solid_torus_slice():
  # a can be choosen freely
  a = 1

  rXY = sp.sqrt(x*x + y*y)
  rCenter = (rXY - R)**2 + z**2 - r**2
  u = a * rCenter / rXY * sp.Matrix([y, -x, 0])

  rhs = curl(curl(u)) + u

  return (u, rhs)

def code_gen(u):
  code = []

  code.append("const real_t x = xVec[0];")
  code.append("const real_t y = xVec[1];")
  code.append("const real_t z = xVec[2];")

  tmp_symbols = sp.numbered_symbols("tmp")
  tmp_assignments, u = sp.cse(u, symbols=tmp_symbols)

  for a in tmp_assignments:
    code.append(f"const real_t {a[0]} = {sp.cxxcode(a[1])};")
  code.append(f"const real_t u0 = {sp.cxxcode(u[0][0])};")
  code.append(f"const real_t u1 = {sp.cxxcode(u[0][1])};")
  code.append(f"const real_t u2 = {sp.cxxcode(u[0][2])};")
  code.append("return Point3D{ u0, u1, u2 };")

  code = "\n".join(code)

  return code

u, rhs = solid_torus_slice()

print("#########")
print("### u ###")
print("#########")
print()
print(code_gen(u))
print()

print("###########")
print("### rhs ###")
print("###########")
print()
print(code_gen(rhs))
print()
