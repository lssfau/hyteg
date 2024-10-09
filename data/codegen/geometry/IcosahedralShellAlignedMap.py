#!/usr/bin/env python3

import sympy as sp

# A: outer ray vertex
A = sp.Matrix(["A_x", "A_y", "A_z"])
# B: inner ray vertex
B = sp.Matrix(["B_x", "B_y", "B_z"])
# C1, C2: the other vertices
C1 = sp.Matrix(["C1_x", "C1_y", "C1_z"])
C2 = sp.Matrix(["C2_x", "C2_y", "C2_z"])

# D1, D2: vertices != A and vertices != B that spans the outer triangle
#         it may be that they are the same as C1, C2 or different -- i.e. the other vertices in the prism
D1 = sp.Matrix(["D1_x", "D1_y", "D1_z"])
D2 = sp.Matrix(["D2_x", "D2_y", "D2_z"])

# A-C1-C2 "cuts" through the prism - we need the plane equation
# eq: n . x + d = 0
n_cut = (C1 - A).cross(C2 - A)
d_cut = -n_cut.dot(A)

# same for triangle at the outer boundary
n_outer = (D1 - A).cross(D2 - A)
d_outer = -n_outer.dot(A)

one_over_A_norm = sp.Symbol("one_over_A_norm")

MtM_inv_Mt = sp.Matrix(
    [
        [
            sp.Symbol("MtM_inv_Mt_0_0", real=True),
            sp.Symbol("MtM_inv_Mt_0_1", real=True),
            sp.Symbol("MtM_inv_Mt_0_2", real=True),
        ],
        [
            sp.Symbol("MtM_inv_Mt_1_0", real=True),
            sp.Symbol("MtM_inv_Mt_1_1", real=True),
            sp.Symbol("MtM_inv_Mt_1_2", real=True),
        ],
    ]
)

NtN_inv_Nt = sp.Matrix(
    [
        [
            sp.Symbol("NtN_inv_Nt_0_0", real=True),
            sp.Symbol("NtN_inv_Nt_0_1", real=True),
            sp.Symbol("NtN_inv_Nt_0_2", real=True),
        ],
        [
            sp.Symbol("NtN_inv_Nt_1_0", real=True),
            sp.Symbol("NtN_inv_Nt_1_1", real=True),
            sp.Symbol("NtN_inv_Nt_1_2", real=True),
        ],
    ]
)


def eval_f_lateral(P: sp.Matrix):

    # plane parallel to outer plane - going through p - we move p along that
    p_plane_n = n_outer
    p_plane_d = -n_outer.dot(P)

    # ray parallel to AB through P:
    # r_bar = P + beta * (B - A)
    #
    # intersection with cutting plane AC:
    beta = -(n_cut.dot(P) + d_cut) / (n_cut.dot(B - A))

    # find w1, w2 in (0, 1) as the portion of the cut plane where the intersection is located
    omega = MtM_inv_Mt * (P + beta * (B - A) - A)

    # intersect that plane with our target ray
    alpha = -p_plane_d / (p_plane_n.dot(A + omega[0] * (D1 - A) + omega[1] * (D2 - A)))

    P_tilde = alpha * (A + omega[0] * (D1 - A) + omega[1] * (D2 - A))

    return P_tilde


def eval_f_inv():
    # the original point
    P_tilde = sp.Matrix(
        [
            sp.Symbol("P_tilde_x", real=True),
            sp.Symbol("P_tilde_y", real=True),
            sp.Symbol("P_tilde_z", real=True),
        ]
    )

    Q = A * (P_tilde.norm() * one_over_A_norm)
    d_parallel = -n_outer.dot(Q)

    alpha = -d_outer / (n_outer.dot(P_tilde))

    c = alpha * P_tilde - A

    omega = NtN_inv_Nt * c

    E = A + omega[0] * (C1 - A) + omega[1] * (C2 - A)

    s = -(n_outer.dot(E) + d_parallel) / n_outer.dot(B - A)

    P = E + s * (B - A)

    return P


def eval_f_radial(P: sp.Matrix):

    return (n_outer.dot(P) / n_outer.dot(A)) * (P / P.norm())


def code_f(P_tilde: sp.Matrix, function_name: str):

    tmp_syms = sp.numbered_symbols("tmp")
    symbols, xnew = sp.cse(P_tilde, symbols=tmp_syms)

    c_code = [f"void {function_name}( const Point3D& x, Point3D& Fx ) const {{"]

    for s in symbols:
        c_code += [f"   real_t {sp.ccode(s[0])} = {sp.ccode(s[1])};"]

    c_code.append(f"   xnew[0] = {sp.ccode(xnew[0][0])};")
    c_code.append(f"   xnew[1] = {sp.ccode(xnew[0][1])};")
    c_code.append(f"   xnew[2] = {sp.ccode(xnew[0][2])};")

    c_code += ["}"]

    print()
    print("\n".join(c_code))
    print()


def code_jacobian(df):
    code = []

    rows, cols = df.shape

    entries = []
    for row in range(rows):
        for col in range(cols):
            entries.append(df[row, col])

    tmp_symbols = sp.numbered_symbols("tmp")
    tmp_assignments, entries = sp.cse(entries, symbols=tmp_symbols)

    code.append("void evalDF( const Point3D& x, Matrix3r& DFx ) const {")

    code.append("    real_t xold_0 = xold[0];")
    code.append("    real_t xold_1 = xold[1];")
    code.append("    real_t xold_2 = xold[2];")

    for a in tmp_assignments:
        code.append(f"    real_t {a[0]} = {sp.ccode(a[1])};")

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


if __name__ == "__main__":

    # The original point
    P = sp.Matrix(
        [
            sp.Symbol("P_x", real=True),
            sp.Symbol("P_y", real=True),
            sp.Symbol("P_z", real=True),
        ]
    )

    P_tilde_only_lateral = eval_f_lateral(P)
    P_tilde_only_radial = eval_f_radial(P)
    P_tilde = eval_f_radial(P_tilde_only_lateral)

    # Forward map
    code_f(P_tilde, "evalF")

    # Inverse map
    P_inv = eval_f_inv()
    code_f(P_inv, "evalFinv")

    J_lateral = P_tilde_only_lateral.jacobian([P[0], P[1], P[2]])
    J_radial = P_tilde_only_radial.jacobian([P[0], P[1], P[2]])

    # Unfortunately sympy is not able to generate the Jacobian directly.
    # So we simply compute the product of the two Jacobians to obtain the Jacobian of the composite function.
    J = (
            J_radial.subs(
                {
                    P[0]: P_tilde_only_lateral[0],
                    P[1]: P_tilde_only_lateral[1],
                    P[2]: P_tilde_only_lateral[2],
                }
            )
            * J_lateral
    )

    # Jacobian of the forward map
    code_jacobian(J)
