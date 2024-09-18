#!/usr/bin/env python3

import sympy as sp

# A: outer ray vertex
A = sp.Matrix([sp.Symbol("A_x", real=True), sp.Symbol("A_y", real=True)])
# B: inner ray vertex
B = sp.Matrix([sp.Symbol("B_x", real=True), sp.Symbol("B_y", real=True)])
# C: the other vertex
C = sp.Matrix([sp.Symbol("C_x", real=True), sp.Symbol("C_y", real=True)])

# D: vertex != A and vertex != B that spans the outer edge
#    it may be C, may be the other vertex in the prism
D = sp.Matrix([sp.Symbol("D_x", real=True), sp.Symbol("D_y", real=True)])

# AC "cuts" through the prism - we need the plane (line here in 2D) equation
# eq: n . x + d = 0
n_cut = sp.Matrix([-(C[1] - A[1]), (C[0] - A[0])])
d_cut = -n_cut.dot(A)

# same for line at the outer boundary
n_outer = sp.Matrix([-(D[1] - A[1]), (D[0] - A[0])])
d_outer = -n_outer.dot(A)

one_over_A_norm = sp.Symbol("one_over_A_norm")

MtM_inv_Mt = sp.Matrix(
    [sp.Symbol("MtM_inv_Mt_0", real=True), sp.Symbol("MtM_inv_Mt_1", real=True)]
).reshape(1, 2)

NtN_inv_Nt = sp.Matrix(
    [sp.Symbol("NtN_inv_Nt_0", real=True), sp.Symbol("NtN_inv_Nt_1", real=True)]
).reshape(1, 2)


def eval_f():
    # the original point
    P = sp.Matrix([sp.Symbol("P_x", real=True), sp.Symbol("P_y", real=True)])

    # ray parallel to AB through p:
    # r_bar = P + beta * (B - A)
    #
    # intersection with cutting plane AC:
    beta = -(n_cut.dot(P) + d_cut) / (n_cut.dot(B - A))

    # find w in (0, 1) as the portion of the cut plane where the intersection is located
    b = P + beta * (B - A) - A
    omega = MtM_inv_Mt * b

    # plane parallel to outer plane - going through p - we move p along that
    p_plane_n = n_outer
    p_plane_d = -n_outer.dot(P)

    s = -p_plane_d / (p_plane_n.dot(A + omega[0] * (D - A)))

    P_tilde = s * (A + omega[0] * (D - A))

    # finding the radius
    q = -(p_plane_d + p_plane_n.dot(A)) / p_plane_n.dot(B - A)

    radius = (A + q * (B - A)).norm()

    P_tilde = (P_tilde / P_tilde.norm()) * radius

    tmp_syms = sp.numbered_symbols("tmp")
    symbols, xnew = sp.cse(P_tilde, symbols=tmp_syms)

    latex_algo = []
    c_code = ["void evalF( const Point3D& x, Point3D& Fx ) const {"]

    for i, s in enumerate(symbols):
        c_code += ["   real_t " + sp.ccode(s[0]) + " = " + sp.ccode(s[1]) + ";"]
        latex_algo += [f"$t_{i} \\gets {sp.latex(s[1])}$"]

    c_code.append(f"   xnew[0] = {sp.ccode(xnew[0][0])};")
    c_code.append(f"   xnew[1] = {sp.ccode(xnew[0][1])};")

    latex_algo += [f"$\\tilde{{P}}_x \\gets {sp.latex(xnew[0][0])}$"]
    latex_algo += [f"$\\tilde{{P}}_y \\gets {sp.latex(xnew[0][1])}$"]

    c_code += ["}"]

    print()
    print("\n".join(c_code))
    print()
    # print("\n".join(latex_algo))

    return P, P_tilde


def eval_f_inv():
    # the original point
    P_tilde = sp.Matrix(
        [sp.Symbol("P_tilde_x", real=True), sp.Symbol("P_tilde_y", real=True)]
    )

    Q = A * (P_tilde.norm() * one_over_A_norm)
    d_parallel = -n_outer.dot(Q)

    alpha = -d_outer / (n_outer.dot(P_tilde))
    R = alpha * P_tilde
    # omega = (A - R).norm() / (A - D).norm()  # only in 2D

    # N = D - A
    c = alpha * P_tilde - A

    omega = NtN_inv_Nt * c

    E = A + omega[0] * (C - A)

    s = -(n_outer.dot(E) + d_parallel) / n_outer.dot(B - A)

    P = E + s * (B - A)

    tmpsyms = sp.numbered_symbols("tmp")
    symbols, xnew = sp.cse(P, symbols=tmpsyms)

    latex_algo = []
    c_code = ["void evalFinv( const Point3D& Fx, Point3D& x ) const {"]

    for i, s in enumerate(symbols):
        c_code += ["   real_t " + sp.ccode(s[0]) + " = " + sp.ccode(s[1]) + ";"]
        latex_algo += [f"$t_{i} \\gets {sp.latex(s[1])}$"]

    c_code.append(f"   xnew[0] = {sp.ccode(xnew[0][0])};")
    c_code.append(f"   xnew[1] = {sp.ccode(xnew[0][1])};")

    latex_algo += [f"$\\tilde{{P}}_x \\gets {sp.latex(xnew[0][0])}$"]
    latex_algo += [f"$\\tilde{{P}}_y \\gets {sp.latex(xnew[0][1])}$"]

    c_code += ["}"]

    print()
    print("\n".join(c_code))
    print()
    # print("\n".join(latex_algo))

    return P_tilde, P


def eval_df(x_old, x_new):

    code = []

    f = sp.Matrix([x_new[0], x_new[1]])

    df = f.jacobian([x_old[0], x_old[1]])

    rows, cols = df.shape

    entries = []
    for row in range(rows):
        for col in range(cols):
            entries.append(df[row, col])

    tmp_symbols = sp.numbered_symbols("tmp")
    tmp_assignments, entries = sp.cse(entries, symbols=tmp_symbols)

    code.append("void evalDF( const Point3D& x, Matrix2r& DFx ) const {")

    code.append("    auto xold_0 = xold[0];")
    code.append("    auto xold_1 = xold[1];")

    for a in tmp_assignments:
        code.append(f"    auto {a[0]} = {sp.ccode(a[1])};")

    for row in range(rows):
        for col in range(cols):
            code.append(
                f"    DF({row}, {col}) = {sp.ccode(entries[row * cols + col])};"
            )

    code.append("    return DF.det();")
    code.append("}")

    code = "\n".join(code)

    print()
    print(code)
    print()

    return df


if __name__ == "__main__":
    xold, xnew = eval_f()
    eval_f_inv()
    eval_df(xold, xnew)
