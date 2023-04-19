/*
* Copyright (c) 2017-2022 Nils Kohl, Andreas Wagner, Fabian Böhm.
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

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P1P0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/egfunctionspace/EGConvTestUtils.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

using walberla::real_t;
using walberla::uint_t;

using hyteg::Point3D;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;
using hyteg::dg::eg::EGSIPGLaplaceOperator;

namespace hyteg {
    namespace dg {
        namespace eg {
            void Stokes2D(const uint_t minLevel, const uint_t maxLevel);

            void HytegSolverCheck2D(const uint_t minLevel, const uint_t maxLevel,
                                    const std::shared_ptr<PrimitiveStorage> &storage);

            void Stokes3D(const uint_t minLevel, const uint_t maxLevel);
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


    if constexpr (true) {
        WALBERLA_LOG_INFO_ON_ROOT("### Testing 2D ###")

        hyteg::dg::eg::Stokes2D(4, 5);
    }

    if constexpr (true) {
        WALBERLA_LOG_INFO_ON_ROOT("### Testing 3D ###")

        hyteg::dg::eg::Stokes3D(3, 4);
    }

    return 0;
}

namespace hyteg {
    namespace dg {
        namespace eg {

            void Stokes2D(const uint_t minLevel, const uint_t maxLevel) {
                auto dummyLambda = [](const Point3D &) -> real_t { return 0; };
                if constexpr (false) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);


                  }

                // quad_4el, inhom
                if constexpr (true) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);




                    if (false) {
                        StokesConvergenceOrderTest<P2P1TaylorHoodStokesOperator>(
                                "P2P1StokesOp2D_quad_4el_inhom",
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
                                            return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                                        },
                                        [](const Point3D &p) -> real_t {
                                            const real_t x = p[0];
                                            const real_t y = p[1];
                                            return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                                        },
                                        dummyLambda,
                                        [](const Point3D &) -> real_t { return 0; }),
                                std::make_shared<P2P1TaylorHoodStokesOperator>(storage, minLevel, maxLevel),
                                storage,
                                minLevel,
                                maxLevel,
                                2);
                    }

                    if (true) {
                        StokesConvergenceOrderTest<EGP0StokesOperatorNitscheBC>(
                                "EGP0StokesOpNitscheBC2D_quad_4el_inhom",
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
                                            return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                                        },
                                        [](const Point3D &p) -> real_t {
                                            const real_t x = p[0];
                                            const real_t y = p[1];
                                            return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                                        },
                                        dummyLambda,
                                        [](const Point3D &) -> real_t { return 0; }),
                                std::make_shared<EGP0StokesOperatorNitscheBC>(storage, minLevel, maxLevel),
                                storage,
                                minLevel,
                                maxLevel,
                                2);
                    }
                }

                // IIPG: quad_4el, inhom
                if constexpr (false) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    EGP0IIPGStokesOperator EGP0StokesOp(storage, minLevel, maxLevel);
                    StokesConvergenceOrderTest<EGP0IIPGStokesOperator>(
                            "EGP0IIPGStokesOp2D_quad_4el_inhom",
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
                                        return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                                    },
                                    [](const Point3D &p) -> real_t {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                                    },
                                    dummyLambda,
                                    [](const Point3D &) -> real_t { return 0; }),
                            std::make_shared<EGP0IIPGStokesOperator>(storage, minLevel, maxLevel),
                            storage,
                            minLevel,
                            maxLevel);
                }

                // quad_16el, inhom
                if constexpr (false) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_16el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    StokesConvergenceOrderTest<EGP0StokesOperator>("EGP0StokesOp2D_quad_16el_inhom",
                                                                   std::make_tuple(
                                                                           [](const Point3D &p) -> real_t {
                                                                               const real_t x = p[0];
                                                                               const real_t y = p[1];
                                                                               return y + 2 * std::sin(M_PI * (x + y)) +
                                                                                      4;
                                                                           },
                                                                           [](const Point3D &p) -> real_t {
                                                                               const real_t x = p[0];
                                                                               const real_t y = p[1];
                                                                               return -x -
                                                                                      2 * std::sin(M_PI * (x + y)) + 3;
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
                                                                               return 4 * std::pow(M_PI, 2) *
                                                                                      std::sin(M_PI * (x + y)) + 2;
                                                                           },
                                                                           [](const Point3D &p) -> real_t {
                                                                               const real_t x = p[0];
                                                                               const real_t y = p[1];
                                                                               return -4 * std::pow(M_PI, 2) *
                                                                                      std::sin(M_PI * (x + y)) - 1;
                                                                           },
                                                                           dummyLambda,
                                                                           [](const Point3D &) -> real_t { return 0; }),
                                                                   std::make_shared<EGP0StokesOperator>(storage,
                                                                                                        minLevel,
                                                                                                        maxLevel),
                                                                   storage,
                                                                   minLevel,
                                                                   maxLevel);
                }

                // IIPG: quad_16el, inhom
                if constexpr (false) {
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_16el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    EGP0IIPGStokesOperator EGP0StokesOp(storage, minLevel, maxLevel);
                    StokesConvergenceOrderTest<EGP0IIPGStokesOperator>(
                            "EGP0IIPGStokesOp2D_quad_16el_inhom",
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
                                        return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                                    },
                                    [](const Point3D &p) -> real_t {
                                        const real_t x = p[0];
                                        const real_t y = p[1];
                                        return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                                    },
                                    dummyLambda,
                                    [](const Point3D &) -> real_t { return 0; }),
                            std::make_shared<EGP0IIPGStokesOperator>(storage, minLevel, maxLevel),
                            storage,
                            minLevel,
                            maxLevel);
                }
            }

            void Stokes3D(const uint_t minLevel, const uint_t maxLevel) {
               // cube_6el, inhom., Nitsche BCs
                if (true) {
                    //MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), 4, 4, 4);
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_6el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);


                    // EG
                    if constexpr (true) {
                        EGP0StokesOperatorNitscheBC EGP0StokesOp(storage, minLevel, maxLevel);
                        WALBERLA_LOG_INFO_ON_ROOT("### EG cube_6el, inhom.,  Nitsche BCs and rhs-integration ###");
                        StokesConvergenceOrderTest<EGP0StokesOperatorNitscheBC>(
                                "EGP0StokesOperatorNitscheBC3D_cube_6el_inhom",
                                std::make_tuple([](const hyteg::Point3D &xx) {
                                                    return -real_c(4) * std::cos(real_c(4) * xx[2]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return real_c(8) * std::cos(real_c(8) * xx[0]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) *
                                                           std::sin(2 * xx[2]);
                                                }),
                                std::make_tuple(
                                        [](const hyteg::Point3D &xx) {
                                            return 4 * std::sin(8 * xx[1]) * std::sin(2 * xx[2]) * std::cos(4 * xx[0]) -
                                                   64 * std::cos(4 * xx[2]);
                                        },
                                        [](const hyteg::Point3D &xx) {
                                            return 8 * std::sin(4 * xx[0]) * std::sin(2 * xx[2]) * std::cos(8 * xx[1]) +
                                                   512 * std::cos(8 * xx[0]);
                                        },
                                        [](const hyteg::Point3D &xx) {
                                            return 2 * std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::cos(2 * xx[2]) -
                                                   8 * std::cos(2 * xx[1]);
                                        },
                                        [](const hyteg::Point3D &) { return 0; }),
                                std::make_shared<EGP0StokesOperatorNitscheBC>(storage, minLevel, maxLevel),
                                storage,
                                minLevel,
                                maxLevel,
                                2);
                    }

                    // P1P0
                    if constexpr (false) {
                        WALBERLA_LOG_INFO_ON_ROOT("### P1P0 cube_6el, inhom. ###");
                        StokesConvergenceOrderTest<hyteg::P1P0StokesOperator>(
                                "P1P0StokesOperator3D_cube_6el_inhom",
                                std::make_tuple([](const hyteg::Point3D &xx) {
                                                    return -real_c(4) * std::cos(real_c(4) * xx[2]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return real_c(8) * std::cos(real_c(8) * xx[0]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) *
                                                           std::sin(2 * xx[2]);
                                                }),
                                std::make_tuple(
                                        [](const hyteg::Point3D &xx) {
                                            return 4 * std::sin(8 * xx[1]) * std::sin(2 * xx[2]) * std::cos(4 * xx[0]) -
                                                   64 * std::cos(4 * xx[2]);
                                        },
                                        [](const hyteg::Point3D &xx) {
                                            return 8 * std::sin(4 * xx[0]) * std::sin(2 * xx[2]) * std::cos(8 * xx[1]) +
                                                   512 * std::cos(8 * xx[0]);
                                        },
                                        [](const hyteg::Point3D &xx) {
                                            return 2 * std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::cos(2 * xx[2]) -
                                                   8 * std::cos(2 * xx[1]);
                                        },
                                        [](const hyteg::Point3D &) { return 0; }),
                                std::make_shared<hyteg::P1P0StokesOperator>(storage, minLevel, maxLevel, 0.1),
                                storage,
                                minLevel,
                                maxLevel,
                                5);
                    }

                    // P2P1
                    if constexpr (false) {
                        WALBERLA_LOG_INFO_ON_ROOT("### P2P1 cube_6el, inhom. ###");
                        StokesConvergenceOrderTest<hyteg::P2P1TaylorHoodStokesOperator>(
                                "P2P1StokesOperator3D_cube_6el_inhom",
                                std::make_tuple([](const hyteg::Point3D &xx) {
                                                    return -real_c(4) * std::cos(real_c(4) * xx[2]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return real_c(8) * std::cos(real_c(8) * xx[0]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                                },
                                                [](const hyteg::Point3D &xx) {
                                                    return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) *
                                                           std::sin(2 * xx[2]);
                                                }),
                                std::make_tuple(
                                        [](const hyteg::Point3D &xx) {
                                            return 4 * std::sin(8 * xx[1]) * std::sin(2 * xx[2]) * std::cos(4 * xx[0]) -
                                                   64 * std::cos(4 * xx[2]);
                                        },
                                        [](const hyteg::Point3D &xx) {
                                            return 8 * std::sin(4 * xx[0]) * std::sin(2 * xx[2]) * std::cos(8 * xx[1]) +
                                                   512 * std::cos(8 * xx[0]);
                                        },
                                        [](const hyteg::Point3D &xx) {
                                            return 2 * std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::cos(2 * xx[2]) -
                                                   8 * std::cos(2 * xx[1]);
                                        },
                                        [](const hyteg::Point3D &) { return 0; }),
                                std::make_shared<hyteg::P2P1TaylorHoodStokesOperator>(storage, minLevel, maxLevel),
                                storage,
                                minLevel,
                                maxLevel,
                                2);
                    }
                }

                // cube_24el, hom. solution
                if (false) {
                    // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_24el.msh");

                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    // cube, hom. solution
                    WALBERLA_LOG_INFO_ON_ROOT("### cube_24el, hom. solution ###");
                    StokesConvergenceOrderTest<EGP0StokesOperator>(
                            "EGP0StokesOp3D_cube_24el_hom",
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        return std::sin(x * x0) * std::sin(x0 * y) * std::sin(x0 * z);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        return std::sin(x * x0) * std::sin(x0 * y) * std::sin(x0 * z);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        return std::sin(x * x0) * std::sin(x0 * y) * std::sin(x0 * z);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                    }),
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * z;
                                        const real_t x1 = 2 * M_PI;
                                        return 4 * std::sin(x0) * std::sin(8 * y) * std::cos(4 * x) +
                                               12 * std::pow(M_PI, 2) * std::sin(M_PI * x0) * std::sin(x * x1) *
                                               std::sin(x1 * y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * z;
                                        const real_t x1 = 2 * M_PI;
                                        return 8 * std::sin(4 * x) * std::sin(x0) * std::cos(8 * y) +
                                               12 * std::pow(M_PI, 2) * std::sin(M_PI * x0) * std::sin(x * x1) *
                                               std::sin(x1 * y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * z;
                                        const real_t x1 = 2 * M_PI;
                                        return 2 * std::sin(4 * x) * std::sin(8 * y) * std::cos(x0) +
                                               12 * std::pow(M_PI, 2) * std::sin(M_PI * x0) * std::sin(x * x1) *
                                               std::sin(x1 * y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        const real_t x1 = x0 * y;
                                        const real_t x2 = std::sin(x1);
                                        const real_t x3 = x * x0;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = x0 * z;
                                        const real_t x6 = x0 * std::sin(x5);
                                        return -x0 * x2 * x4 * std::cos(x5) - x2 * x6 * std::cos(x3) -
                                               x4 * x6 * std::cos(x1);
                                    }),
                            std::make_shared<EGP0StokesOperator>(storage, minLevel, maxLevel),
                            storage,
                            minLevel,
                            maxLevel - 1,
                            2);
                }

                // IIPG: cube_24el, hom. solution
                if (false) {
                    // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_24el.msh");

                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    EGP0IIPGStokesOperator EGP0StokesOp(storage, minLevel - 1, maxLevel);

                    // cube, hom. solution
                    WALBERLA_LOG_INFO_ON_ROOT("### cube_24el, hom. solution ###");
                    StokesConvergenceOrderTest<EGP0IIPGStokesOperator>(
                            "EGP0IIPGStokesOp3D_cube_24el_hom",
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        return std::sin(x * x0) * std::sin(x0 * y) * std::sin(x0 * z);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        return std::sin(x * x0) * std::sin(x0 * y) * std::sin(x0 * z);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        return std::sin(x * x0) * std::sin(x0 * y) * std::sin(x0 * z);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                    }),
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * z;
                                        const real_t x1 = 2 * M_PI;
                                        return 4 * std::sin(x0) * std::sin(8 * y) * std::cos(4 * x) +
                                               12 * std::pow(M_PI, 2) * std::sin(M_PI * x0) * std::sin(x * x1) *
                                               std::sin(x1 * y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * z;
                                        const real_t x1 = 2 * M_PI;
                                        return 8 * std::sin(4 * x) * std::sin(x0) * std::cos(8 * y) +
                                               12 * std::pow(M_PI, 2) * std::sin(M_PI * x0) * std::sin(x * x1) *
                                               std::sin(x1 * y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * z;
                                        const real_t x1 = 2 * M_PI;
                                        return 2 * std::sin(4 * x) * std::sin(8 * y) * std::cos(x0) +
                                               12 * std::pow(M_PI, 2) * std::sin(M_PI * x0) * std::sin(x * x1) *
                                               std::sin(x1 * y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2 * M_PI;
                                        const real_t x1 = x0 * y;
                                        const real_t x2 = std::sin(x1);
                                        const real_t x3 = x * x0;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = x0 * z;
                                        const real_t x6 = x0 * std::sin(x5);
                                        return -x0 * x2 * x4 * std::cos(x5) - x2 * x6 * std::cos(x3) -
                                               x4 * x6 * std::cos(x1);
                                    }),
                            std::make_shared<EGP0IIPGStokesOperator>(storage, minLevel, maxLevel),

                            storage,
                            minLevel,
                            maxLevel - 1,
                            2,
                            true);
                }

      }

        } // namespace eg
    } // namespace dg
} // namespace hyteg
