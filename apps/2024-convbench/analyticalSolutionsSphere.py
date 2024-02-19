import assess
import numpy as np
from math import pi

mtRes = lambda nTan, maxLevel: 2 * (nTan - 1) * 2 ** (maxLevel)

# rDash = 3 * pi / 2
rDash = 1.72
# rDashdp = 1.97
# rDashdm = 1.47

solutionDirichletSmooth2d = assess.CylindricalStokesSolutionSmoothZeroSlip(2, 2)
solutionFreeslipSmooth2d = assess.CylindricalStokesSolutionSmoothFreeSlip(2, 2)
solutionFreeZeroslipSmooth2d = assess.CylindricalStokesSolutionSmoothFreeZeroSlip(2, 2)

solutionAboveDirichletDelta2d = assess.CylindricalStokesSolutionDeltaZeroSlip(2, 1)
solutionBelowDirichletDelta2d = assess.CylindricalStokesSolutionDeltaZeroSlip(2, -1)

solutionAboveFreeslipDelta2d = assess.CylindricalStokesSolutionDeltaFreeSlip(2, 1)
solutionBelowFreeslipDelta2d = assess.CylindricalStokesSolutionDeltaFreeSlip(2, -1)

solutionAboveFreeZeroslipDelta2d = assess.CylindricalStokesSolutionDeltaFreeZeroSlip(
    2, 1
)
solutionBelowFreeZeroslipDelta2d = assess.CylindricalStokesSolutionDeltaFreeZeroSlip(
    2, -1
)

#################### 3D ###################

l = 2
m = 2
k = 2

solutionDirichletSmooth3d = assess.SphericalStokesSolutionSmoothZeroSlip(l, m, k)
solutionFreeZeroslipSmooth3d = assess.SphericalStokesSolutionSmoothFreeZeroSlip(l, m, k)
solutionFreeslipSmooth3d = assess.SphericalStokesSolutionSmoothFreeSlip(l, m, k)

solutionAboveDirichletDelta3d = assess.SphericalStokesSolutionDeltaZeroSlip(l, m, 1)
solutionBelowDirichletDelta3d = assess.SphericalStokesSolutionDeltaZeroSlip(l, m, -1)

solutionAboveFreeslipDelta3d = assess.SphericalStokesSolutionDeltaFreeSlip(l, m, 1)
solutionBelowFreeslipDelta3d = assess.SphericalStokesSolutionDeltaFreeSlip(l, m, -1)

solutionAboveFreeZeroslipDelta3d = assess.SphericalStokesSolutionDeltaFreeZeroSlip(
    l, m, 1
)
solutionBelowFreeZeroslipDelta3d = assess.SphericalStokesSolutionDeltaFreeZeroSlip(
    l, m, -1
)

###########################################################################
######################## General utils ####################################
###########################################################################


def dYdtheta(x):
    return [assess.spherical.dYdtheta(2, 2, x[1], x[2])]


def dYdphi(x):
    return [assess.spherical.dYdphi(2, 2, x[1], x[2])]


###########################################################################
################################### 2D ####################################
###########################################################################


def getDirichletPressureSmooth2d(x):
    return [solutionDirichletSmooth2d.pressure_cartesian(x)]


def getDirichletPressureDelta2d(x):
    r = np.linalg.norm(x)
    if r > rDash:
        return [solutionAboveDirichletDelta2d.pressure_cartesian(x)]
    else:
        return [solutionBelowDirichletDelta2d.pressure_cartesian(x)]


def getFreeslipPressureSmooth2d(x):
    return [solutionFreeslipSmooth2d.pressure_cartesian(x)]


def getFreeslipPressureDelta2d(x):
    r = np.linalg.norm(x)
    if r > rDash:
        return [solutionAboveFreeslipDelta2d.pressure_cartesian(x)]
    else:
        return [solutionBelowFreeslipDelta2d.pressure_cartesian(x)]

    # if r > rDashdp:
    #     return [solutionAboveAboveFreeslipDelta2d.pressure_cartesian(x)]
    # elif r < rDashdm:
    #     return [solutionBelowBelowFreeslipDelta2d.pressure_cartesian(x)]
    # else:
    #     return [solutionAboveBelowFreeslipDelta2d.pressure_cartesian(x)]


def getFreeZeroslipPressureSmooth2d(x):
    return [solutionFreeZeroslipSmooth2d.pressure_cartesian(x)]


def getFreeZeroslipPressureDelta2d(x):
    r = np.linalg.norm(x)
    if r > rDash:
        return [solutionAboveFreeZeroslipDelta2d.pressure_cartesian(x)]
    else:
        return [solutionBelowFreeZeroslipDelta2d.pressure_cartesian(x)]


def getDirichletVelocitySmooth2d(x):
    # print(x)
    # print(solution2d.velocity_cartesian(x))

    return list(solutionDirichletSmooth2d.velocity_cartesian(x))


def getDirichletVelocityDelta2d(x):
    # print(x)
    # print(solution2d.velocity_cartesian(x))

    r = np.linalg.norm(x)
    if r > rDash:
        return list(solutionAboveDirichletDelta2d.velocity_cartesian(x))
    else:
        return list(solutionBelowDirichletDelta2d.velocity_cartesian(x))


def getFreeslipVelocitySmooth2d(x):
    # print(x)
    # print(solution2d.velocity_cartesian(x))
    # print("Freeslip from Python")
    return list(solutionFreeslipSmooth2d.velocity_cartesian(x))


def getFreeslipVelocityDelta2d(x):
    # print(x)
    # print(solution2d.velocity_cartesian(x))

    r = np.linalg.norm(x)

    if r > rDash:
        return list(solutionAboveFreeslipDelta2d.velocity_cartesian(x))
    else:
        return list(solutionBelowFreeslipDelta2d.velocity_cartesian(x))

    # if r > rDashdp:
    #     return list(solutionAboveAboveFreeslipDelta2d.velocity_cartesian(x))
    # elif r < rDashdm:
    #     return list(solutionBelowBelowFreeslipDelta2d.velocity_cartesian(x))
    # else:
    #     return list(solutionAboveBelowFreeslipDelta2d.velocity_cartesian(x))


def getFreeZeroslipVelocitySmooth2d(x):
    # print(x)
    # print(solution2d.velocity_cartesian(x))
    # print("Freeslip from Python")
    return list(solutionFreeZeroslipSmooth2d.velocity_cartesian(x))


def getFreeZeroslipVelocityDelta2d(x):
    # print(x)
    # print(solution2d.velocity_cartesian(x))

    r = np.linalg.norm(x)
    if r > rDash:
        return list(solutionAboveFreeZeroslipDelta2d.velocity_cartesian(x))
    else:
        return list(solutionBelowFreeZeroslipDelta2d.velocity_cartesian(x))


def getDeltaRho2d(x):
    return [solutionDirichletSmooth2d.delta_rho_cartesian(x)]


###########################################################################
################################### 3D ####################################
###########################################################################


def getDirichletPressureSmooth3d(x):
    return [solutionDirichletSmooth3d.pressure_cartesian(x)]


def getDirichletPressureDelta3d(x):
    r = np.linalg.norm(x)
    if r > rDash:
        return [solutionAboveDirichletDelta3d.pressure_cartesian(x)]
    else:
        return [solutionBelowDirichletDelta3d.pressure_cartesian(x)]


def getDirichletVelocitySmooth3d(x):
    return list(solutionDirichletSmooth3d.velocity_cartesian(x))


def getDirichletVelocityDelta3d(x):
    # print("x from python = ", x)
    # print("velocity from python = ", solution.velocity_cartesian(x))
    r = np.linalg.norm(x)
    if r > rDash:
        return list(solutionAboveDirichletDelta3d.velocity_cartesian(x))
    else:
        return list(solutionBelowDirichletDelta3d.velocity_cartesian(x))


# Freeslip


def getFreeslipPressureSmooth3d(x):
    return [solutionFreeslipSmooth3d.pressure_cartesian(x)]


def getFreeslipPressureDelta3d(x):
    r = np.linalg.norm(x)
    if r > rDash:
        return [solutionAboveFreeslipDelta3d.pressure_cartesian(x)]
    else:
        return [solutionBelowFreeslipDelta3d.pressure_cartesian(x)]


def getFreeslipVelocitySmooth3d(x):
    return list(solutionFreeslipSmooth3d.velocity_cartesian(x))


def getFreeslipVelocityDelta3d(x):
    # print("x from python = ", x)
    # print("velocity from python = ", solution.velocity_cartesian(x))
    r = np.linalg.norm(x)
    if r > rDash:
        return list(solutionAboveFreeslipDelta3d.velocity_cartesian(x))
    else:
        return list(solutionBelowFreeslipDelta3d.velocity_cartesian(x))


# Freeslip and Zeroslip


def getFreeZeroslipPressureSmooth3d(x):
    return [solutionFreeZeroslipSmooth3d.pressure_cartesian(x)]


def getFreeZeroslipPressureDelta3d(x):
    r = np.linalg.norm(x)
    if r > rDash:
        return [solutionAboveFreeZeroslipDelta3d.pressure_cartesian(x)]
    else:
        return [solutionBelowFreeZeroslipDelta3d.pressure_cartesian(x)]


def getFreeZeroslipVelocitySmooth3d(x):
    return list(solutionFreeZeroslipSmooth3d.velocity_cartesian(x))


def getFreeZeroslipVelocityDelta3d(x):
    # print("x from python = ", x)
    # print("velocity from python = ", solution.velocity_cartesian(x))
    r = np.linalg.norm(x)
    if r > rDash:
        return list(solutionAboveFreeZeroslipDelta3d.velocity_cartesian(x))
    else:
        return list(solutionBelowFreeZeroslipDelta3d.velocity_cartesian(x))


def getSPH(x):
    # print("x from python = ", x)
    _, theta, phi = assess.to_spherical(x)
    # print("SPH from python = ", assess.Y(l, m, theta, phi))
    return [assess.Y(l, m, theta, phi)]
