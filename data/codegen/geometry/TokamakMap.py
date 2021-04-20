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

    """

      //      auto toroidalAngle = atan2( xold[1], xold[0] );
      //
      //      {
      //         auto alpha             = toroidalAngle - toroidalStartAngle_ - real_c( toroidalPrism_ ) * toroidalAngleIncrement_;
      //         auto beta              = 0.5 * ( pi - toroidalAngleIncrement_ );
      //         auto gamma             = pi - alpha - beta;
      //         auto toroidalRadiusNew = ( std::sin( gamma ) * ( std::sqrt( xold[0] * xold[0] + xold[1] * xold[1] ) / std::sin( beta ) ) );
      //
      //         xnew[0] = toroidalRadiusNew * std::cos( toroidalAngle );
      //         xnew[1] = toroidalRadiusNew * std::sin( toroidalAngle );
      //         xnew[2] = xold[2];
      //      }

    """

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

    """
    // then we rotate the mapped centroid around the z-axis and translate it to the origin
      // this way we can find the angle and therefore the prism ID via polar coordinates in the x-z-plane

      auto C = torusCoordinates( radiusOriginToCenterOfTube_, 0, toroidalAngle,  0 );
      auto centroidTrafoToOrigin = centroid - C;
      centroidTrafoToOrigin = Point3D( {
                                       std::cos( -toroidalAngle ) * centroidTrafoToOrigin[0] - std::sin( -toroidalAngle ) * centroidTrafoToOrigin[1],
                                       std::sin( -toroidalAngle ) * centroidTrafoToOrigin[0] + std::cos( -toroidalAngle ) * centroidTrafoToOrigin[1],
                                       centroidTrafoToOrigin[2]
                                   } );

      auto poloidalAngle = std::atan2( centroidTrafoToOrigin[2], centroidTrafoToOrigin[0] );
      if ( poloidalAngle < 0 )
      {
         poloidalAngle += 2 * pi;
      }
      poloidalPrism_ = uint_c((poloidalAngle - poloidalStartAngle_) /  poloidalAngleIncrement_);
    """

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

    """

         auto poloidalAnglePrism = poloidalStartAngle_ + real_c( poloidalPrism_ ) * poloidalAngleIncrement_;
         auto C = torusCoordinates( radiusOriginToCenterOfTube_, 0, toroidalAngle,  0 );
         auto A = torusCoordinates( radiusOriginToCenterOfTube_, tubeLayerRadii_.back(), toroidalAngle,  poloidalAnglePrism );

         auto pLocal = xnew - C;
         auto aLocal = A - C;

         auto poloidalAngle = angle( aLocal, pLocal );
         auto alpha = poloidalAngle - poloidalAnglePrism;
         auto beta = 0.5 * ( pi - poloidalAngleIncrement_ );
         auto gamma = pi - alpha - beta;

         auto poloidalRadiusNew = ( std::sin( gamma ) * ( pLocal.norm() / std::sin( beta ) ) );

         xnew = torusCoordinates( radiusOriginToCenterOfTube_, poloidalRadiusNew, toroidalAngle, poloidalAngle );


    """

    poloidalAnglePrism = poloidalStartAngle_ + poloidalPrism_ * poloidalAngleIncrement_
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
        sp.sqrt(pLocal_0 ** 2 + pLocal_1 ** 2 + pLocal_2 ** 2) / sp.sin(beta)
    )

    xnew_0, xnew_1, xnew_2 = torus_coordinates(
        radiusOriginToCenterOfTube_, poloidalRadiusNew, toroidalAngle, poloidalAngle
    )

    # 3. It's tokamak time

    r0_, r1_, r2_, delta_ = sp.symbols("r0_ r1_ r2_ delta_")

    xnew_0 = (
        1
        + (poloidalRadiusNew / tubeLayerRadiiBack_)
        * (r1_ / r0_)
        * sp.cos(poloidalAngle + sp.asin(delta_) * sp.sin(poloidalAngle))
    ) * sp.cos(toroidalAngle)

    xnew_1 = (
        1
        + (poloidalRadiusNew / tubeLayerRadiiBack_)
        * (r1_ / r0_)
        * sp.cos(poloidalAngle + sp.asin(delta_) * sp.sin(poloidalAngle))
    ) * sp.sin(toroidalAngle)

    xnew_2 = (
        (poloidalRadiusNew / tubeLayerRadiiBack_) * (r2_ / r0_) * sp.sin(poloidalAngle)
    )

    tmp_symbols = sp.numbered_symbols("tmp")
    tmp_assignments, xnew = sp.cse([xnew_0, xnew_1, xnew_2], symbols=tmp_symbols)

    for a in tmp_assignments:
        code.append(f"auto {a[0]} = {sp.ccode(a[1])};")
    code.append(f"xnew[0] = {xnew[0]};")
    code.append(f"xnew[1] = {xnew[1]};")
    code.append(f"xnew[2] = {xnew[2]};")

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

    code = "\n".join(code)

    return code


if __name__ == "__main__":
    code, _, _ = evalF()
    print(code)

    print()

    code = evalDF()
    print(code)
