/*
* Copyright (c) 2017-2022 Nils Kohl.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/p0functionspace/P0P0MassForm.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/egfunctionspace/EGConvTestUtils.hpp"

using walberla::real_t;
using walberla::uint_t;


using hyteg::Point3D;
using hyteg::dg::eg::EGSIPGLaplaceOperator;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;

namespace hyteg {
    namespace dg {
        namespace eg {

            void
            Epsilon2D(const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr<PrimitiveStorage> &storage);

            void SmoothViscosityTest2D(const uint_t minLevel, const uint_t maxLevel,
                                       const std::shared_ptr<PrimitiveStorage> &storage);
            void
            Epsilon3D(const uint_t minLevel, const uint_t maxLevel);

        } // namespace eg
    } // namespace dg
} // namespace hyteg

int main(int argc, char *argv[]) {
    walberla::MPIManager::instance()->initializeMPI(&argc, &argv);
    walberla::MPIManager::instance()->useWorldComm();
    hyteg::PETScManager petscManager(&argc, &argv);

    /* commandline arguments for petsc solver:
   -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space
   */

    uint_t minLevel = 3;

    if(false)
    {
        WALBERLA_LOG_INFO_ON_ROOT("### Testing const Epsilon 2D ###")
        auto meshInfo = hyteg::MeshInfo::meshRectangle(
                hyteg::Point2D({-1, -1}), hyteg::Point2D({1, 1}), hyteg::MeshInfo::CRISSCROSS, 1, 1);
        hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                  walberla::uint_c(
                                                          walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

        uint_t maxLevel = 6;

        hyteg::dg::eg::Epsilon2D( minLevel, maxLevel, storage );
    }

    if(false)
    {
        WALBERLA_LOG_INFO_ON_ROOT("### Testing variable Epsilon 2D ###")
        auto meshInfo = hyteg::MeshInfo::meshRectangle(
                hyteg::Point2D({-1, -1}), hyteg::Point2D({1, 1}), hyteg::MeshInfo::CRISSCROSS, 1, 1);
        hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                  walberla::uint_c(
                                                          walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

        uint_t maxLevel = 5;

        hyteg::dg::eg::SmoothViscosityTest2D( minLevel, maxLevel, storage );
    }
    if(true)
    {
        WALBERLA_LOG_INFO_ON_ROOT("### Testing const. Epsilon 3D ###")

        uint_t maxLevel = 5;

        hyteg::dg::eg::Epsilon3D( minLevel, maxLevel );
    }

    return 0;
}

namespace hyteg {
    namespace dg {
        namespace eg {

            void
            Epsilon2D(const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr<PrimitiveStorage> &storage) {
                EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_1(storage, minLevel, maxLevel,
                                                             [](const hyteg::Point3D &) { return 1; });

                auto dummyLambda = [](const Point3D &) -> real_t { return 0; };
                StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                        "EGP0ConstEpsilonOp2D_hom_asym",
                        std::make_tuple(
                                [](const Point3D &xx) -> real_t {
                                    return std::sin(M_PI * (xx[0] + 1.0) / 2.0) * std::sin(M_PI * (xx[1] + 1.0) / 2.0) *
                                           std::exp(xx[1]);
                                },
                                [](const Point3D &xx) -> real_t {
                                    return std::sin(M_PI * (xx[0] + 1.0) / 2.0) * std::sin(M_PI * (xx[1] + 1.0) / 2.0) *
                                           std::exp(xx[1]);
                                },
                                dummyLambda,
                                [](const Point3D &xx) -> real_t { return xx[0] + xx[1]; }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x1 = std::sin(x0);
                                    const real_t x2 = std::exp(y);
                                    const real_t x3 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x4 = std::sin(x3);
                                    const real_t x5 = x2 * x4;
                                    const real_t x6 = x2 * std::cos(x3);
                                    const real_t x7 = std::cos(x0);
                                    const real_t x8 = std::pow(M_PI, 2);
                                    return 0.75 * x1 * x2 * x4 * x8 - 1.0 * x1 * x5 - 1.0 * M_PI * x1 * x6 -
                                           0.5 * M_PI * x5 * x7 -
                                           0.25 * x6 * x7 * x8 + 1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x1 = std::sin(x0);
                                    const real_t x2 = std::exp(y);
                                    const real_t x3 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x4 = std::sin(x3);
                                    const real_t x5 = x2 * x4;
                                    const real_t x6 = x2 * std::cos(x3);
                                    const real_t x7 = std::cos(x0);
                                    const real_t x8 = std::pow(M_PI, 2);
                                    return 0.75 * x1 * x2 * x4 * x8 - 2.0 * x1 * x5 - 2.0 * M_PI * x1 * x6 -
                                           0.5 * M_PI * x5 * x7 -
                                           0.25 * x6 * x7 * x8 + 1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x1 = std::sin(x0);
                                    const real_t x2 = std::exp(y);
                                    const real_t x3 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x4 = x2 * std::sin(x3);
                                    const real_t x5 = M_PI_2;
                                    return -x1 * x2 * x5 * std::cos(x3) - x1 * x4 - x4 * x5 * std::cos(x0);
                                }),
                    EGP0EpsilonOp_mu_1,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                        "EGP0ConstEpsilonOp2D_hom_sym",
                        std::make_tuple(
                                [](const Point3D &xx) -> real_t {
                                    return std::sin(M_PI * (xx[0] + 1.0) / 2.0) * std::sin(M_PI * (xx[1] + 1.0) / 2.0);
                                },
                                [](const Point3D &xx) -> real_t {
                                    return std::sin(M_PI * (xx[0] + 1.0) / 2.0) * std::sin(M_PI * (xx[1] + 1.0) / 2.0);
                                },
                                dummyLambda,
                                [](const Point3D &xx) -> real_t { return xx[0] + xx[1]; }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    const real_t x1 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x2 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    return 0.75 * x0 * std::sin(x1) * std::sin(x2) -
                                           0.25 * x0 * std::cos(x1) * std::cos(x2) + 1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    const real_t x1 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x2 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    return 0.75 * x0 * std::sin(x1) * std::sin(x2) -
                                           0.25 * x0 * std::cos(x1) * std::cos(x2) + 1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x1 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x2 = M_PI_2;
                                    return -(x2 * std::sin(x0) * std::cos(x1) + x2 * std::sin(x1) * std::cos(x0));
                                }),
                    EGP0EpsilonOp_mu_1,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                        "EGP0ConstEpsilonOp2D_inhom_sym",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];

                                    return std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];

                                    return std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                },
                                dummyLambda,
                                [](const Point3D &xx) -> real_t { return xx[0] + xx[1]; }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    return 0.5 * x0 * std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           0.25 * x0 * std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) + 1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    return 0.25 * x0 * std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           0.5 * x0 * std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) + 1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI_2;
                                    return -x0 * std::cos(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) -
                                           x0 * std::cos(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                }),
                    EGP0EpsilonOp_mu_1,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                        "EGP0EpsilonOp2D_inhom_constvisc",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                },
                                dummyLambda,
                                [](const Point3D &xx) -> real_t { return xx[0] + xx[1]; }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    return 0.5 * x0 * std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           0.25 * x0 * std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) + 1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    return 0.25 * x0 * std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           0.5 * x0 * std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) + 1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI_2;
                                    return -x0 * std::cos(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) -
                                           x0 * std::cos(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                }),
                        EGP0EpsilonOp_mu_1,
                        storage,
                        minLevel,
                        maxLevel);
       }
/*
            void IncreasingSteepnessTest(const uint_t minLevel, const uint_t maxLevel,
                                         const std::shared_ptr<PrimitiveStorage> &storage) {
                real_t alpha = 1;
                auto dummyLambda = [](const Point3D &) -> real_t { return 0; };

                auto viscosity = [&alpha](const hyteg::Point3D &p) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    return 500.0 * std::tanh(alpha * M_PI * (x + y)) + 501.0;
                };
                auto u_x = [&alpha](const Point3D &p) -> real_t {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    return (0.5 * std::tanh(alpha * M_PI * (x + y)) + 0.5) * std::sin(M_PI * (2 * x + 2 * y));
                };
                auto u_y = [&alpha](const Point3D &p) -> real_t {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    return -(0.5 * std::tanh(alpha * M_PI * (x + y)) + 0.5) * std::sin(M_PI * (2 * x + 2 * y));
                };
                auto pressure = [](const Point3D &p) -> real_t {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    return 2.5 * M_PI * cos(2 * M_PI * x) * cos(M_PI * y);
                };
                auto f_x = [&alpha](const Point3D &p) -> real_t {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t x0 = std::pow(M_PI, 2);
                    const real_t x1 = 2 * x;
                    const real_t x2 = M_PI * alpha;
                    const real_t x3 = std::tanh(x2 * (x + y));
                    const real_t x4 = 0.5 * x3 + 0.5;
                    const real_t x5 = -x4;
                    const real_t x6 = M_PI * (x1 + 2 * y);
                    const real_t x7 = std::sin(x6);
                    const real_t x8 = x0 * x7;
                    const real_t x9 = 2.0 * x4;
                    const real_t x10 = 1000.0 * x3 + 1002.0;
                    const real_t x11 = std::cos(x6);
                    const real_t x12 = M_PI * x11;
                    const real_t x13 = 1.0 * x12;
                    const real_t x14 = 1 - std::pow(x3, 2);
                    const real_t x15 = x14 * x2;
                    const real_t x16 = 1000.0 * x15;
                    return -5.0 * x0 * std::sin(M_PI * x1) * std::cos(M_PI * y) - x10 * (-2.0 * x5 * x8 - x8 * x9) -
                           x10 *
                           (-1.0 * std::pow(alpha, 2) * x14 * x3 * x8 + 2.0 * alpha * x0 * x11 * x14 - 4.0 * x4 * x8) -
                           x16 * (x12 * x9 + 0.5 * x15 * x7) - x16 * (x13 * x4 + x13 * x5);
                };
                auto f_y = [&alpha](const Point3D &p) -> real_t {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t x0 = std::pow(M_PI, 2);
                    const real_t x1 = 2 * x;
                    const real_t x2 = M_PI * alpha;
                    const real_t x3 = std::tanh(x2 * (x + y));
                    const real_t x4 = 0.5 * x3 + 0.5;
                    const real_t x5 = -x4;
                    const real_t x6 = M_PI * (x1 + 2 * y);
                    const real_t x7 = std::sin(x6);
                    const real_t x8 = x0 * x7;
                    const real_t x9 = 2.0 * x8;
                    const real_t x10 = 1000.0 * x3 + 1002.0;
                    const real_t x11 = std::cos(x6);
                    const real_t x12 = 1.0 * M_PI * x11;
                    const real_t x13 = 1 - std::pow(x3, 2);
                    const real_t x14 = x13 * x2;
                    const real_t x15 = 1000.0 * x14;
                    return -2.5 * x0 * std::sin(M_PI * y) * std::cos(M_PI * x1) - x10 * (-x4 * x9 - x5 * x9) -
                           x10 * (1.0 * std::pow(alpha, 2) * x0 * x13 * x3 * x7 - 2.0 * alpha * x0 * x11 * x13 -
                                  4.0 * x5 * x8) -
                           x15 * (x12 * x4 + x12 * x5) - x15 * (2.0 * M_PI * x11 * x5 - 0.5 * x14 * x7);
                };

                {
                    hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_alpha_1_smooth(
                            storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<hyteg::P2P1ElementwiseAffineEpsilonStokesOperator>(
                            "P2P1EpsilonOp_mu_alpha1_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            P2P1ElementwiseEpsilonOp_mu_alpha_1_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);

                    EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_alpha_1_smooth(storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                            "EGP0EpsilonOp_mu_alpha1_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            EGP0EpsilonOp_mu_alpha_1_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);
                }

                {
                    alpha = 10;
                    hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_alpha_10_smooth(
                            storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<hyteg::P2P1ElementwiseAffineEpsilonStokesOperator>(
                            "P2P1EpsilonOp_mu_alpha10_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            P2P1ElementwiseEpsilonOp_mu_alpha_10_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);

                    EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_alpha_10_smooth(storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                            "EGP0EpsilonOp_mu_alpha10_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            EGP0EpsilonOp_mu_alpha_10_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);
                }

                {
                    alpha = 50;
                    hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_alpha_50_smooth(
                            storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<hyteg::P2P1ElementwiseAffineEpsilonStokesOperator>(
                            "P2P1EpsilonOp_mu_alpha50_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            P2P1ElementwiseEpsilonOp_mu_alpha_50_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);

                    EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_alpha_50_smooth(storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                            "EGP0EpsilonOp_mu_alpha50_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            EGP0EpsilonOp_mu_alpha_50_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);
                }
                {
                    alpha = 1000;
                    hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_alpha_1000_smooth(
                            storage, minLevel, maxLevel, viscosity);
                    StokesConvergenceOrderTest<hyteg::P2P1ElementwiseAffineEpsilonStokesOperator>(
                            "P2P1EpsilonOp_mu_alpha1000_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            P2P1ElementwiseEpsilonOp_mu_alpha_1000_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);

                    EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_alpha_1000_smooth(storage, minLevel, maxLevel,
                                                                                 viscosity);
                    StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                            "EGP0EpsilonOp_mu_alpha1000_smooth",
                            std::make_tuple(u_x, u_y, dummyLambda, pressure),
                            std::make_tuple(f_x, f_y, dummyLambda, [](const Point3D &) -> real_t { return 0; }),
                            EGP0EpsilonOp_mu_alpha_1000_smooth,
                            storage,
                            minLevel,
                            maxLevel,
                            false);
                }
            }
*/
            void SmoothViscosityTest2D(const uint_t minLevel, const uint_t maxLevel,
                                       const std::shared_ptr<PrimitiveStorage> &storage) {
                auto dummyLambda = [](const Point3D &) -> real_t { return 0; };

                EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_smooth(storage, minLevel, maxLevel,
                                                                  [](const hyteg::Point3D &p) {
                                                                      const real_t x = p[0];
                                                                      const real_t y = p[1];
                                                                      return std::exp(x) * std::sin(
                                                                              M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) *
                                                                             std::sin(M_PI *
                                                                                      ((1.0 / 2.0) * y + 1.0 / 2.0)) +
                                                                             1;
                                                                  });

                StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                        "EGP0EpsilonOp2D_divFree_smoothVisc",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return y + 2 * std::sin(M_PI * (x + y)) + 4;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return -x - 2 * std::sin(M_PI * (x + y)) + 3;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return 2 * x - y + 1;
                                }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * (x + y);
                                    const real_t x1 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x2 = std::exp(x) * std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                    const real_t x3 = x2 * std::sin(x1);
                                    return 4.0 * std::pow(M_PI, 2) * (x3 + 1) * std::sin(x0) -
                                           4.0 * M_PI * ((1.0 / 2.0) * M_PI * x2 * std::cos(x1) + x3) * std::cos(x0) +
                                           2;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::exp(x);
                                    const real_t x1 = std::pow(M_PI, 2);
                                    const real_t x2 = M_PI * (x + y);
                                    const real_t x3 = std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0));
                                    const real_t x4 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    return 2.0 * x0 * x1 * x3 * std::cos(x2) * std::cos(x4) -
                                           4.0 * x1 * (x0 * x3 * std::sin(x4) + 1) * std::sin(x2) - 1;
                                },
                                dummyLambda,
                                [](const Point3D &) -> real_t { return 0; }),
                        EGP0EpsilonOp_mu_smooth,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                        "EGP0EpsilonOp2D_asym_smoothVisc",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::exp(y) * std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::exp(y) * std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) +
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0));
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];

                                    return x + y;
                                }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];

                                    const real_t x0 = std::exp(y);
                                    const real_t x1 = std::pow(M_PI, 2);
                                    const real_t x2 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x3 = std::sin(x2);
                                    const real_t x4 = std::exp(x);
                                    const real_t x5 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x6 = std::sin(x5);
                                    const real_t x7 = x4 * x6;
                                    const real_t x8 = x3 * x7;
                                    const real_t x9 = x8 + 1;
                                    const real_t x10 = std::cos(x2);
                                    const real_t x11 = x0 * x10;
                                    const real_t x12 = std::cos(x5);
                                    const real_t x13 = 0.25 * M_PI;
                                    const real_t x14 = 0.5 * x0 * x3 + x11 * x13;
                                    return 0.5 * x0 * x1 * x3 * x9 -
                                           1.0 * M_PI * x11 * ((1.0 / 2.0) * M_PI * x10 * x7 + x8) -
                                           M_PI * x12 * x3 * x4 * (x12 * x13 + x14) -
                                           2 * x9 * (-0.125 * x1 * x6 + x14) + 1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::exp(x);
                                    const real_t x1 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x2 = std::sin(x1);
                                    const real_t x3 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x4 = M_PI * std::cos(x3);
                                    const real_t x5 = std::exp(y);
                                    const real_t x6 = x2 * x5;
                                    const real_t x7 = 1.0 * x6;
                                    const real_t x8 = std::pow(M_PI, 2);
                                    const real_t x9 = std::sin(x3);
                                    const real_t x10 = x0 * x9;
                                    const real_t x11 = x10 * x2;
                                    const real_t x12 = 2 * x11 + 2;
                                    const real_t x13 = std::cos(x1);
                                    const real_t x14 = M_PI * x13;
                                    return -x0 * x2 * x4 * (0.5 * x4 + x7) - x12 * (x7 - 0.25 * x8 * x9) -
                                           x12 * (0.25 * M_PI * x13 * x5 - 0.125 * x6 * x8) -
                                           2 * ((1.0 / 2.0) * x10 * x14 + x11) *
                                           (0.25 * x14 * x5 + 0.25 * x4 + 0.5 * x6) + 1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI_2;
                                    const real_t x1 = std::exp(y);
                                    const real_t x2 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    return -x0 * x1 * std::cos(x2) -
                                           x0 * std::cos(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) - x1 * std::sin(x2);
                                }),
                        EGP0EpsilonOp_mu_smooth,
                        storage,
                        minLevel,
                        maxLevel);
            }

            void Epsilon3D(const uint_t minLevel, const uint_t maxLevel) {

                {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/tet_1el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    EGP0EpsilonStokesOperator EGP0ConstEpsilonOp(storage, minLevel, maxLevel,
                                                                  [](const hyteg::Point3D &) { return 1; });


                    // tet_1el, inhom.
                    if (true) {
                        WALBERLA_LOG_INFO_ON_ROOT("### tet_1el, inhom. ###");

                        StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                                "EGP0EpsilonOp3D_tet_1el_inhom",
                                std::make_tuple(
                                        [](const hyteg::Point3D &xx) { return -real_c(4) * std::cos(real_c(4) * xx[2]); },
                                        [](const hyteg::Point3D &xx) { return real_c(8) * std::cos(real_c(8) * xx[0]); },
                                        [](const hyteg::Point3D &xx) { return -real_c(2) * std::cos(real_c(2) * xx[1]); },
                                        [](const hyteg::Point3D &xx) {
                                            return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                        }),
                                std::make_tuple(
                                        [](const hyteg::Point3D &p) {
                                            const real_t x = p[0];
                                            const real_t y = p[1];
                                            const real_t z = p[2];
                                            return 4*std::sin(8*y)*std::sin(2*z)*std::cos(4*x) - 64.0*std::cos(4*z);
                                        },
                                        [](const hyteg::Point3D &p) {
                                            const real_t x = p[0];
                                            const real_t y = p[1];
                                            const real_t z = p[2];
                                            return 8*std::sin(4*x)*std::sin(2*z)*std::cos(8*y) + 512.0*std::cos(8*x);
                                        },
                                        [](const hyteg::Point3D &p) {
                                            const real_t x = p[0];
                                            const real_t y = p[1];
                                            const real_t z = p[2];
                                            return 2*std::sin(4*x)*std::sin(8*y)*std::cos(2*z) - 8.0*std::cos(2*y);
                                        },
                                        [](const hyteg::Point3D &) { return 0; }),
                                EGP0ConstEpsilonOp,
                                storage,
                                minLevel,
                                maxLevel, 5, true);


                    }
                }
                    // cube_6el, inhom. solution
                if (true) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_6el.msh");

                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);


                       EGP0EpsilonStokesOperator EGP0ConstEpsilonOp(storage, minLevel, maxLevel,
                                                                     [](const hyteg::Point3D &) { return 1; });

                    // cube, hom. solution
                    WALBERLA_LOG_INFO_ON_ROOT("### cube_6el, inhom. solution ###");
                    StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                            "EGP0ConstEpsilonOp3D_cube_6el_inhom",
                            std::make_tuple(
                                    [](const hyteg::Point3D &xx) { return -real_c(4) * std::cos(real_c(4) * xx[2]); },
                                    [](const hyteg::Point3D &xx) { return real_c(8) * std::cos(real_c(8) * xx[0]); },
                                    [](const hyteg::Point3D &xx) { return -real_c(2) * std::cos(real_c(2) * xx[1]); },
                                    [](const hyteg::Point3D &xx) {
                                        return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                    }),
                            std::make_tuple(
                                    [](const hyteg::Point3D &p) {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        const real_t z = p[2];
                                        return 4*std::sin(8*y)*std::sin(2*z)*std::cos(4*x) - 64.0*std::cos(4*z);
                                    },
                                    [](const hyteg::Point3D &p) {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        const real_t z = p[2];
                                        return 8*std::sin(4*x)*std::sin(2*z)*std::cos(8*y) + 512.0*std::cos(8*x);
                                    },
                                    [](const hyteg::Point3D &p) {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        const real_t z = p[2];
                                        return 2*std::sin(4*x)*std::sin(8*y)*std::cos(2*z) - 8.0*std::cos(2*y);
                                    },
                                    [](const hyteg::Point3D &) { return 0; }),
                            EGP0ConstEpsilonOp,
                            storage,
                            minLevel,
                            maxLevel,
                            5,
                            true);
                }

                // cube_24el, inhom. solution
                if (true) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_24el.msh");

                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    EGP0EpsilonStokesOperator EGP0ConstEpsilonOp(storage, minLevel, maxLevel,
                                                                  [](const hyteg::Point3D &) { return 1; });

                    // cube, hom. solution
                    WALBERLA_LOG_INFO_ON_ROOT("### cube_24el, inhom. solution ###");
                    StokesConvergenceOrderTest<EGP0EpsilonStokesOperator>(
                            "EGP0ConstEpsilonOp3D_cube_24el_inhom",
                            std::make_tuple(
                                    [](const hyteg::Point3D &xx) { return -real_c(4) * std::cos(real_c(4) * xx[2]); },
                                    [](const hyteg::Point3D &xx) { return real_c(8) * std::cos(real_c(8) * xx[0]); },
                                    [](const hyteg::Point3D &xx) { return -real_c(2) * std::cos(real_c(2) * xx[1]); },
                                    [](const hyteg::Point3D &xx) {
                                        return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                    }),
                            std::make_tuple(
                                    [](const hyteg::Point3D &p) {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        const real_t z = p[2];
                                        return 4*std::sin(8*y)*std::sin(2*z)*std::cos(4*x) - 64.0*std::cos(4*z);
                                    },
                                    [](const hyteg::Point3D &p) {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        const real_t z = p[2];
                                        return 8*std::sin(4*x)*std::sin(2*z)*std::cos(8*y) + 512.0*std::cos(8*x);
                                    },
                                    [](const hyteg::Point3D &p) {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        const real_t z = p[2];
                                        return 2*std::sin(4*x)*std::sin(8*y)*std::cos(2*z) - 8.0*std::cos(2*y);
                                    },
                                    [](const hyteg::Point3D &) { return 0; }),
                            EGP0ConstEpsilonOp,
                            storage,
                            minLevel,
                            maxLevel,
                            5,
                            true);
                }
            }
        } // namespace eg
    } // namespace dg
} // namespace hyteg
