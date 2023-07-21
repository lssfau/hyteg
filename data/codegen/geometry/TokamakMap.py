import sympy as sp
import numpy as np


def torus_coordinates(
    torusRadiusToTubeCenter, tubeRadius, toroidalAngle, poloidalAngle
):
    return (
        (torusRadiusToTubeCenter + tubeRadius * sp.cos(poloidalAngle))
        * sp.cos(toroidalAngle),
        (torusRadiusToTubeCenter + tubeRadius * sp.cos(poloidalAngle))
        * sp.sin(toroidalAngle),
        tubeRadius * sp.sin(poloidalAngle),
    )


def evalF():

    code = []

    code.append("auto xold_0 = xold[0];")
    code.append("auto xold_1 = xold[1];")
    code.append("auto xold_2 = xold[2];")

    xold_0, xold_1, xold_2 = sp.symbols("xold_0 xold_1 xold_2")
    toroidalStartAngle_ = sp.symbols("toroidalStartAngle_")
    toroidalPrism_ = sp.symbols("toroidalPrism_")
    toroidalAngleIncrement_ = sp.symbols("toroidalAngleIncrement_")

    radiusOriginToCenterOfTube_ = sp.symbols("radiusOriginToCenterOfTube_")
    tubeLayerRadiiBack_ = sp.symbols("tubeLayerRadiiBack_")

    poloidalAngleIncrement_ = sp.symbols("poloidalAngleIncrement_")
    poloidalPrism_ = sp.symbols("poloidalPrism_")
    poloidalStartAngle_ = sp.symbols("poloidalStartAngle_")

    # 1. Blend onto torus (toroidal)

    toroidalAngle = sp.atan2(xold_1, xold_0)

    alpha = (
        toroidalAngle - toroidalStartAngle_ - toroidalPrism_ * toroidalAngleIncrement_
    )
    beta = 0.5 * (np.pi - toroidalAngleIncrement_)
    gamma = np.pi - alpha - beta

    toroidalRadiusNew = sp.sin(gamma) * (
        sp.sqrt(xold_0 * xold_0 + xold_1 * xold_1) / sp.sin(beta)
    )

    xnew_0 = toroidalRadiusNew * sp.cos(toroidalAngle)
    xnew_1 = toroidalRadiusNew * sp.sin(toroidalAngle)
    xnew_2 = xold_2

    # 2. Blend onto torus (poloidal)

    c_0, c_1, c_2 = torus_coordinates(radiusOriginToCenterOfTube_, 0, toroidalAngle, 0)

    xTrafoToOrigin_0 = xnew_0 - c_0
    xTrafoToOrigin_1 = xnew_1 - c_1
    xTrafoToOrigin_2 = xnew_2 - c_2

    xTrafoToOrigin_0 = (
        sp.cos(-toroidalAngle) * xTrafoToOrigin_0
        - sp.sin(-toroidalAngle) * xTrafoToOrigin_1
    )
    xTrafoToOrigin_1 = (
        sp.sin(-toroidalAngle) * xTrafoToOrigin_0
        + sp.cos(-toroidalAngle) * xTrafoToOrigin_1
    )
    xTrafoToOrigin_2 = xTrafoToOrigin_2

    poloidalAngle = sp.atan2(xTrafoToOrigin_2, xTrafoToOrigin_0)
    poloidalAngle -= poloidalStartAngle_
    poloidalAngle = sp.Piecewise(
        (poloidalAngle + 2 * np.pi, poloidalAngle < 0), (poloidalAngle, True)
    )

    poloidalAnglePrism = poloidalPrism_ * poloidalAngleIncrement_
    poloidalAnglePrism = sp.Piecewise(
        (poloidalAnglePrism + 2 * np.pi, poloidalAnglePrism < 0),
        (poloidalAnglePrism, True),
    )
    a_0, a_1, a_2 = torus_coordinates(
        radiusOriginToCenterOfTube_,
        tubeLayerRadiiBack_,
        toroidalAngle,
        poloidalAnglePrism,
    )

    pLocal_0 = xnew_0 - c_0
    pLocal_1 = xnew_1 - c_1
    pLocal_2 = xnew_2 - c_2

    alpha = poloidalAngle - poloidalAnglePrism
    beta = 0.5 * (np.pi - poloidalAngleIncrement_)
    gamma = np.pi - alpha - beta
    poloidalRadiusNew = sp.sin(gamma) * (
        (sp.sqrt(pLocal_0 ** 2 + pLocal_1 ** 2 + pLocal_2 ** 2) / sp.sin(beta))
    )

    xnew_poloidal_0, xnew_poloidal_1, xnew_poloidal_2 = torus_coordinates(
        radiusOriginToCenterOfTube_,
        poloidalRadiusNew,
        toroidalAngle,
        poloidalAngle + poloidalStartAngle_,
    )

    # There is an issue if the poloidal radius (the distance to the center of the tube) is zero.
    # In that case, the poloidal angle cannot be defined properly.
    # We fix this here by just projecting points that are near the poloidal origin to the poloidal origin, and by
    # setting the angle to some random value (it is zero here).
    # The eps could probably be set to just trigger when we have _exact_ floating point zero, but maybe(?) it's a good
    # idea to prevent stuff blowing up by making it a little larger.
    # Keeping the xnew_ variables. If the Tokamak part was removed, those are the correct coordinates to return.
    eps = 1e-14
    distToTubeOrigin = sp.sqrt(xTrafoToOrigin_0**2 + xTrafoToOrigin_2**2)
    xnew_0 = sp.Piecewise((xnew_0, distToTubeOrigin < eps), (xnew_poloidal_0, True))
    xnew_1 = sp.Piecewise((xnew_1, distToTubeOrigin < eps), (xnew_poloidal_1, True))
    xnew_2 = sp.Piecewise((xnew_2, distToTubeOrigin < eps), (xnew_poloidal_2, True))
    poloidalAngle = sp.Piecewise((0, distToTubeOrigin < eps), (poloidalAngle, True))
    poloidalRadiusNew = sp.Piecewise((0, distToTubeOrigin < eps), (poloidalRadiusNew, True))

    # # 3. It's tokamak time

    r1_, r2_, delta_ = sp.symbols("r1_ r2_ delta_")

    poloidalAngle += poloidalStartAngle_

    xnew_0 = (
        radiusOriginToCenterOfTube_
        + (poloidalRadiusNew / tubeLayerRadiiBack_)
        * r1_
        * sp.cos(poloidalAngle + sp.asin(delta_) * sp.sin(poloidalAngle))
    ) * sp.cos(toroidalAngle)

    xnew_1 = (
        radiusOriginToCenterOfTube_
        + (poloidalRadiusNew / tubeLayerRadiiBack_)
        * r1_
        * sp.cos(poloidalAngle + sp.asin(delta_) * sp.sin(poloidalAngle))
    ) * sp.sin(toroidalAngle)

    xnew_2 = (poloidalRadiusNew / tubeLayerRadiiBack_) * r2_ * sp.sin(poloidalAngle)

    tmp_symbols = sp.numbered_symbols("tmp")
    tmp_assignments, xnew = sp.cse([xnew_0, xnew_1, xnew_2], symbols=tmp_symbols)

    for a in tmp_assignments:
        code.append(f"auto {a[0]} = {sp.ccode(a[1])};")
    code.append(f"xnew[0] = {sp.ccode(xnew[0])};")
    code.append(f"xnew[1] = {sp.ccode(xnew[1])};")
    code.append(f"xnew[2] = {sp.ccode(xnew[2])};")

    code = "\n".join(code)

    return code, (xold_0, xold_1, xold_2), (xnew_0, xnew_1, xnew_2)


def evalDF():

    code = []

    code.append("auto xold_0 = xold[0];")
    code.append("auto xold_1 = xold[1];")
    code.append("auto xold_2 = xold[2];")

    _, xold, xnew = evalF()

    f = sp.Matrix([xnew[0], xnew[1], xnew[2]])

    df = f.jacobian([xold[0], xold[1], xold[2]])

    rows, cols = df.shape

    entries = []
    for row in range(rows):
        for col in range(cols):
            entries.append(df[row, col])

    tmp_symbols = sp.numbered_symbols("tmp")
    tmp_assignments, entries = sp.cse(entries, symbols=tmp_symbols)

    for a in tmp_assignments:
        code.append(f"auto {a[0]} = {sp.ccode(a[1])};")

    for row in range(rows):
        for col in range(cols):
            code.append(f"DF({row}, {col}) = {entries[row * cols + col]};")

    code.append("return DF.det();")

    code = "\n".join(code)

    return code


if __name__ == "__main__":

    code, _, _ = evalF()

    print("###############")
    print("### evalF() ###")
    print("###############")
    print()
    print(code)
    print()

    code = evalDF()

    print("################")
    print("### evalDF() ###")
    print("################")
    print()
    print(code)
    print()
