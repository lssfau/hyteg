#!/usr/bin/env python3

import numpy as np
from sympy import *
import p2_reference as p2ref
import quadpy

eta, xi, nu = symbols('eta xi nu')

x1, x2, x3, x4 = symbols('coords[0][0] coords[1][0] coords[2][0] coords[3][0]')
y1, y2, y3, y4 = symbols('coords[0][1] coords[1][1] coords[2][1] coords[3][1]')
z1, z2, z3, z4 = symbols('coords[0][2] coords[1][2] coords[2][2] coords[3][2]')

DFinv_11, DFinv_12, DFinv_13, DFinv_21, DFinv_22, DFinv_23, DFinv_31, DFinv_32, DFinv_33 = symbols('DFinv_11, DFinv_12, DFinv_13, DFinv_21, DFinv_22, DFinv_23, DFinv_31, DFinv_32, DFinv_33')
DFinv_11.name = 'DFinv(0,0)'
DFinv_12.name = 'DFinv(0,1)'
DFinv_13.name = 'DFinv(0,2)'
DFinv_21.name = 'DFinv(1,0)'
DFinv_22.name = 'DFinv(1,1)'
DFinv_23.name = 'DFinv(1,2)'
DFinv_31.name = 'DFinv(2,0)'
DFinv_32.name = 'DFinv(2,1)'
DFinv_33.name = 'DFinv(2,2)'

DF_11, DF_12, DF_13, DF_21, DF_22, DF_23, DF_31, DF_32, DF_33 = symbols('DF_11, DF_12, DF_13, DF_21, DF_22, DF_23, DF_31, DF_32, DF_33')
DF_11.name = 'DF(0,0)'
DF_12.name = 'DF(0,1)'
DF_13.name = 'DF(0,2)'
DF_21.name = 'DF(1,0)'
DF_22.name = 'DF(1,1)'
DF_23.name = 'DF(1,2)'
DF_31.name = 'DF(2,0)'
DF_32.name = 'DF(2,1)'
DF_33.name = 'DF(2,2)'

x_hat, y_hat, z_hat = symbols('x_hat[0], x_hat[1], x_hat[2]')
x_tilde, y_tilde, z_tilde = symbols('x_tilde[0], x_tilde[1], x_tilde[2]')

K = symbols('K')
K.name = 'callback( x_tilde, geometryMap_ )'

N_2d = 6 # 6
N_3d = 10 # 10

def generate_2d_integration_rule():
    scheme_tri = quadpy.triangle.strang_fix_cowper_02()
    reference_tri = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    vol = quadpy.nsimplex.get_vol(reference_tri)
    x_q_tri = quadpy.nsimplex.transform(scheme_tri.points.T, reference_tri.T).T
    w_tri = vol * scheme_tri.weights
    return (w_tri, x_q_tri)

def generate_2d_forms():
    R = simplify(Matrix([[x2-x1, x3-x1], [y2-y1, y3-y1]]) * Matrix([[eta], [xi]]) + Matrix([[x1], [y1]]))
    DR = R.jacobian([eta, xi])
    DRinv = DR**(-1)

    phi_hat, grad_phi_hat = p2ref.create_2d()

    DFinv = Matrix([[DFinv_11, DFinv_12], [DFinv_21, DFinv_22]])

    grad_phi = []
    for i in range(N_2d):
      grad_phi.append(simplify(DFinv.T * DRinv.T * grad_phi_hat[i]))

    dx_tilde = simplify(abs(det(DR)))
    dx = simplify(1 / abs(det(DFinv)) * dx_tilde)

    forms_2d = []
    for i in range(N_2d):
        for j in range(N_2d):
            print('Processing 2D form {},{}'.format(i,j))
            form = grad_phi[i].T * K * grad_phi[j] * dx
            form = form[0]
            form = form.subs([(eta, x_hat), (xi, y_hat)])
            # form = simplify(form)
            forms_2d.append(form)

    return forms_2d, R

def generate_3d_integration_rule():
    scheme_tet = quadpy.tetrahedron.yu_1()
    reference_tet = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    vol = quadpy.nsimplex.get_vol(reference_tet)
    x_q_tet = quadpy.nsimplex.transform(scheme_tet.points.T, reference_tet.T).T
    w_tet = vol * scheme_tet.weights
    return (w_tet, x_q_tet)

def generate_3d_forms():
    R = simplify(Matrix([[x2-x1, x3-x1, x4-x1], [y2-y1, y3-y1, y4-y1], [z2-z1, z3-z1, z4-z1]]) * Matrix([[eta], [xi], [nu]]) + Matrix([[x1], [y1], [z1]]))
    DR = R.jacobian([eta, xi, nu])
    DRinv = DR**(-1)

    print('generating 3d basis functions')
    phi_hat, grad_phi_hat = p2ref.create_3d()

    DF = Matrix([[DF_11, DF_12, DF_13], [DF_21, DF_22, DF_23], [DF_31, DF_32, DF_33]])
    DFinv = DF**(-1)

    # DFinv = Matrix([[DFinv_11, DFinv_12, DFinv_13], [DFinv_21, DFinv_22, DFinv_23], [DFinv_31, DFinv_32, DFinv_33]])

    grad_phi = []
    for i in range(N_3d):
      print('computing gradient of {} basis function'.format(i))
      grad_phi.append(simplify(DFinv.T * DRinv.T * grad_phi_hat[i]))

    print('computing dx_tilde')
    # dx_tilde = simplify(1 / abs(det(DRinv)))
    # dx_tilde = 1 / abs(det(DRinv))
    dx_tilde = simplify(abs(det(DR)))

    print('computing dx')
    # dx = simplify(1 / abs(det(DFinv)) * dx_tilde)
    dx = simplify(1 / abs(det(DFinv)) * dx_tilde)

    forms_3d = []
    for i in range(N_3d):
        for j in range(N_3d):
            print('Processing 3D form {},{}'.format(i,j))
            form = grad_phi[i].T * K * grad_phi[j] * dx
            form = form[0]
            form = form.subs([(eta, x_hat), (xi, y_hat), (nu, z_hat)])
            # form = simplify(form)
            forms_3d.append(form)

    return forms_3d, R

def sympyToC(classname):
    c_code = "class {} : public P2FormHyTeG {{\n".format(classname)
    c_code += "public:\n"

    ### 2D

    tmpsyms = numbered_symbols("tmp")
    forms_2d, R = generate_2d_forms()
    symbols, simple = cse(forms_2d, symbols=tmpsyms)
    w, x_q = generate_2d_integration_rule()

    c_code += "  void evalQuadraturePoint2D(const Point3D& x_hat, const Point3D& x_tilde, const std::array<Point3D,3>& coords, const Matrix2r& DFinv, real_t w, Matrix6r& elMat) const\n"
    c_code += "  {\n"
    for s in symbols:
        #print s
        c_code +=  "    real_t " +cxxcode(s[0]) + " = " + cxxcode(s[1]) + ";\n"

    for i in range(N_2d):
        for j in range(N_2d):
            c_code +=  "    elMat({},{}) += w * ".format(i,j) + cxxcode(simple[N_2d*i + j])+";\n"
    c_code += "  }\n\n"
    c_code += "  void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final\n"
    c_code += "  {\n"

    c_code += "    Point3D x_hat;\n"
    c_code += "    Point3D x_tilde;\n"
    c_code += "    Matrix2r DFinv;\n"
    for i in range(N_2d):
        for j in range(N_2d):
            c_code += "    elMat({},{}) = 0;\n".format(i,j)

    for q in range(len(w)):
        for j in range(2):
            c_code += "    x_hat[{}] = {};\n".format(j, x_q[q][j])

        x_tilde = R.subs([(eta, x_q[q][0]), (xi, x_q[q][1])])
        for j in range(2):
            c_code += "    x_tilde[{}] = {};\n".format(j, x_tilde[j])
        c_code += "    geometryMap_->evalDFinv(x_tilde, DFinv);\n"
        c_code += "    evalQuadraturePoint2D(x_hat, x_tilde, coords, DFinv, {}, elMat);\n".format(w[q])

    c_code += "  }\n\n"

    ### 3D

    forms_3d, R = generate_3d_forms()
    tmpsyms = numbered_symbols("tmp")
    symbols, simple = cse(forms_3d, symbols=tmpsyms)
    w, x_q = generate_3d_integration_rule()

    c_code += "  void evalQuadraturePoint3D(const Point3D& x_hat, const Point3D& x_tilde, const std::array<Point3D,4>& coords, const Matrix3r& DF, real_t w, Matrix10r& elMat) const\n"
    c_code += "  {\n"
    for s in symbols:
        #print s
        c_code +=  "    real_t " +cxxcode(s[0]) + " = " + cxxcode(s[1]) + ";\n"

    for i in range(N_3d):
        for j in range(N_3d):
            c_code +=  "    elMat({},{}) += w * ".format(i,j) + cxxcode(simple[N_3d*i + j])+";\n"
    c_code += "  }\n\n"
    c_code += "  void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final\n"
    c_code += "  {\n"

    c_code += "    Point3D x_hat;\n"
    c_code += "    Point3D x_tilde;\n"
    c_code += "    Matrix3r DF;\n"
    for i in range(N_3d):
        for j in range(N_3d):
            c_code += "    elMat({},{}) = 0;\n".format(i,j)

    for q in range(len(w)):
        for j in range(3):
            c_code += "    x_hat[{}] = {};\n".format(j, x_q[q][j])

        x_tilde = R.subs([(eta, x_q[q][0]), (xi, x_q[q][1]), (nu, x_q[q][2])])
        for j in range(3):
            c_code += "    x_tilde[{}] = {};\n".format(j, x_tilde[j])
        c_code += "    geometryMap_->evalDF(x_tilde, DF);\n"
        c_code += "    evalQuadraturePoint3D(x_hat, x_tilde, coords, DF, {}, elMat);\n".format(w[q])

    # c_code += "    }\n"
    c_code += "  }\n\n"

    c_code += "  static std::function< real_t( const Point3D&, const std::shared_ptr< GeometryMap >& ) > callback;\n\n"
    c_code += "};\n\n"
    return c_code

c_code = "#pragma once\n\n"
c_code += "#include \"hyteg/forms/form_hyteg_base/P2FormHyTeG.hpp\"\n\n"
c_code += "using walberla::real_c;\n\n"
c_code += "namespace hyteg {\n\n"

c_code += sympyToC('P2Form_divKgradBlending')

c_code += "}\n"

print(c_code)
