import sympy as sp
from sympy.utilities.codegen import codegen, CCodeGen
from sympy.codegen.ast import Assignment, CodeBlock

x, y, u, v, a, b = sp.symbols("x y u v a b")

r, phi, k, T, p, r0, r_, phi_ = sp.symbols("r phi k T p r0 r_ phi_", positive=True)

eta_0, eta_a, eta_r, x_a = sp.symbols("eta_0 eta_a eta_r x_a")

# r = sp.sqrt(x**2 + y**2)
# phi = sp.atan(y/x)

# ur = (1/r)*k*sp.sin(r)*sp.cos(k*phi)
# uphi = -1*sp.cos(r)*sp.sin(k*phi)

# eta_0 = 1.0
# eta_a = 1.0
# eta_r = 100.0
# x_a = 0.5

# eta = 1.0
eta = eta_0 + eta_a * sp.exp(-eta_r * (r - x_a)**2)

rho = 1.0
# rho = a * sp.exp(b * (r0 - r))

ur = (k / r) * sp.sin(r) * sp.cos(k * phi)
uphi = -sp.cos(r) * sp.sin(k * phi)

# p = -((k**2 + r**2)/k)*sp.cos(r)*sp.cos(k*phi) + ((2.0/r) - (3*r/(k**2)))*ur
p = (
    ((2 / r**2) * k * sp.sin(r) * sp.cos(k * phi))
    - ((1.0 / k) * sp.sin(r) * sp.cos(k * phi))
    - ((k**2 + r**2 + 1) / (k * r)) * sp.cos(r) * sp.cos(k * phi)
)

u = sp.cos(phi) * ur - sp.sin(phi) * uphi
v = sp.sin(phi) * ur + sp.cos(phi) * uphi

u = u / rho
v = v / rho

# u = sp.sin(phi)/rho
# v = -sp.cos(phi)/rho

divu = sp.diff(u * rho, x) + sp.diff(v * rho, y)

print(sp.simplify(divu))


drdx = x / r
drdy = y / r

dphidx = -y / (r**2)
dphidy = x / (r**2)

dudx = sp.diff(u, r) * drdx + sp.diff(u, phi) * dphidx
dudy = sp.diff(u, r) * drdy + sp.diff(u, phi) * dphidy

dvdx = sp.diff(v, r) * drdx + sp.diff(v, phi) * dphidx
dvdy = sp.diff(v, r) * drdy + sp.diff(v, phi) * dphidy

dudx = dudx * eta
dudy = dudy * eta

dvdx = dvdx * eta
dvdy = dvdy * eta

d2udx2 = sp.diff(dudx, x) + sp.diff(dudx, r) * drdx + sp.diff(dudx, phi) * dphidx
d2udy2 = sp.diff(dudy, y) + sp.diff(dudy, r) * drdy + sp.diff(dudy, phi) * dphidy

d2vdx2 = sp.diff(dvdx, x) + sp.diff(dvdx, r) * drdx + sp.diff(dvdx, phi) * dphidx
d2vdy2 = sp.diff(dvdy, y) + sp.diff(dvdy, r) * drdy + sp.diff(dvdy, phi) * dphidy

dpdx = sp.diff(p, r) * drdx + sp.diff(p, phi) * dphidx
dpdy = sp.diff(p, r) * drdy + sp.diff(p, phi) * dphidy

rhs_x = -d2udx2 - d2udy2 + dpdx
rhs_y = -d2vdx2 - d2vdy2 + dpdy

def printWithCse(exprToPrint):
    rhs_assignments = []
    rhs_cse = sp.cse(exprToPrint)
    for (sym, expr) in rhs_cse[0]:
        rhs_assignments.append(Assignment(sp.Symbol('real_t ' + sym.name), expr))

    rhs_assignments.append(Assignment(sp.Symbol('rhs_x'), rhs_cse[1][0]))

    rhs_x_codeblock = CodeBlock(*rhs_assignments)

    print(sp.ccode(rhs_x_codeblock))
