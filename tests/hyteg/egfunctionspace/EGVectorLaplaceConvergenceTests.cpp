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

#include <random>

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "mixedOperator//P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"
#include "mixedOperator/EGOperators.hpp"
#include "constantStencilOperator/P1ConstantOperator.cpp"

using hyteg::MeshInfo;
using hyteg::Point2D;
using hyteg::Point3D;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using hyteg::dg::eg::EGSIPGLaplaceOperator;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0StokesOperator;

// scalar lambda for one component of analytical solution and rhs
typedef std::function<real_t(const hyteg::PointND<real_t, 3> &p)> ScalarLambda;

// tuple of function for solution (u,p) and rhs of vector values stokes equation
typedef std::tuple<ScalarLambda, ScalarLambda, ScalarLambda> LambdaTuple;

namespace hyteg {
    real_t VectorLaplace(const std::string &name,
                         MeshInfo meshInfo,
                         const uint_t &level,
                         LambdaTuple sol_tuple,
                         LambdaTuple rhs_tuple,
                         uint_t solverType,
                         bool writeVTK = false) {
        // MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
        SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, 1);

        EGFunction<idx_t> numerator("numerator", storage, level, level);
        numerator.enumerate(level);

        EGSIPGLaplaceOperator L(storage, level, level);
        EGMassOperator M_EG(storage, level, level);

        EGFunction<real_t> u("u", storage, level, level);
        EGFunction<real_t> f("f", storage, level, level);
        EGFunction<real_t> rhs("rhs", storage, level, level);
        EGFunction<real_t> sol("sol", storage, level, level);
        EGFunction<real_t> err("err", storage, level, level);
        EGFunction<real_t> tmp("tmp", storage, level, level);

        EGFunction<real_t> tmp2("tmp2", storage, level, level);
        EGFunction<real_t> Merr("Merr", storage, level, level);

        //auto [u_x_expr, u_y_expr, u_z_expr] = sol_tuple;
        auto u_x_expr = std::get<0>(sol_tuple);
        auto u_y_expr = std::get<1>(sol_tuple);
        auto u_z_expr = std::get<2>(sol_tuple);
        //auto [f_x_expr, f_y_expr, f_z_expr] = rhs_tuple;
        auto f_x_expr = std::get<0>(rhs_tuple);
        auto f_y_expr = std::get<1>(rhs_tuple);
        auto f_z_expr = std::get<2>(rhs_tuple);

        if (storage->hasGlobalCells()) {
            sol.interpolate({u_x_expr, u_y_expr, u_z_expr}, level, All);
            u.getConformingPart()->interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
            f.interpolate({f_x_expr, f_y_expr, f_z_expr}, level, All);
            rhs.getConformingPart()->interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
        } else {
            sol.interpolate({u_x_expr, u_y_expr}, level, All);
            u.getConformingPart()->interpolate({u_x_expr, u_y_expr}, level, DirichletBoundary);
            f.interpolate({f_x_expr, f_y_expr}, level, All);
            rhs.getConformingPart()->interpolate({u_x_expr, u_y_expr}, level, DirichletBoundary);
        }
        M_EG.apply(f, rhs, level, All, Replace);


        switch (solverType) {
            case 0: {
                PETScCGSolver<EGSIPGLaplaceOperator> solver(storage, level, numerator, 1e-15, 1e-15);
                solver.solve(L, u, rhs, level);
                break;
            }
            case 1: {
                CGSolver<EGSIPGLaplaceOperator> solver(storage, level, level, 10000, 1e-15);
                solver.solve(L, u, rhs, level);
                break;
            }
            default: WALBERLA_ABORT("No solver chosen.");
        }

        // u.getDiscontinuousPart()->interpolate( 0, level, DirichletBoundary );
        err.assign({1.0, -1.0}, {sol, u}, level, All);

        // calculate the error in the L2 norm
        M_EG.apply(err, Merr, level, Inner, Replace);
        auto discrL2 = sqrt(err.dotGlobal(Merr, level, Inner));
        //auto discrL2 = sqrt( err.dotGlobal( err, level, Inner ) / real_c( numberOfGlobalDoFs( u, level ) ) );

        if (writeVTK) {
            VTKOutput vtk("/mnt/c/Users/Fabia/OneDrive/Desktop/hyteg_premerge/hyteg-build/output", name, storage);
            vtk.add(u);
            vtk.add(*u.getConformingPart());
            vtk.add(*u.getDiscontinuousPart());
            vtk.add(sol);
            vtk.add(*sol.getConformingPart());
            vtk.add(*sol.getDiscontinuousPart());
            vtk.add(f);
            vtk.add(*f.getConformingPart());
            vtk.add(*f.getDiscontinuousPart());
            vtk.add(rhs);
            vtk.add(*rhs.getConformingPart());
            vtk.add(*rhs.getDiscontinuousPart());
            vtk.add(err);
            vtk.add(*err.getConformingPart());
            vtk.add(*err.getDiscontinuousPart());
            vtk.write(level);
        }

        return discrL2;
    }

    void runTestcase(const std::string &name,
                     const uint_t &minLevel,
                     const uint_t &maxLevel,
                     MeshInfo meshInfo,
                     LambdaTuple sol_tuple,
                     LambdaTuple rhs_tuple,
                     uint_t solverType,
                     bool writeVTK = false) {
        auto l2ConvRate = std::pow(2, -(int(1) + 1));
        auto convRateEps = l2ConvRate * 0.1;
        auto err = VectorLaplace(name, meshInfo, minLevel, sol_tuple, rhs_tuple, solverType, writeVTK);
        WALBERLA_LOG_INFO_ON_ROOT("error level " << minLevel << ": " << err);
        real_t computedRate = 0.0;
        uint_t l = 0;
        for (l = minLevel + 1; l <= maxLevel; l++) {
            auto errFiner = VectorLaplace(name, meshInfo, l, sol_tuple, rhs_tuple, solverType, writeVTK);
            computedRate = errFiner / err;

            WALBERLA_LOG_INFO_ON_ROOT("error level " << l << ": " << errFiner);
            WALBERLA_LOG_INFO_ON_ROOT(
                    "Convergence L2 rate level " << l << " vs level " << l - 1 << "; computed: " << computedRate);

            err = errFiner;
        }
        WALBERLA_CHECK_LESS_EQUAL(computedRate,
                                  l2ConvRate + convRateEps,
                                  "Convergence L2 rate level " << l << " vs level " << l - 1
                                                               << " not sufficiently small (computed: " << computedRate
                                                               << ", estimated + eps: " << l2ConvRate + convRateEps
                                                               << ")");
    }

    void testEGVectorLaplace2D(uint_t solverType, bool writeVTK, uint_t minLevel, uint_t maxLevel);

    void testEGVectorLaplace3D(uint_t solverType, bool writeVTK, uint_t minLevel, uint_t maxLevel);

} // namespace hyteg

int main(int argc, char **argv) {
    walberla::mpi::Environment MPIenv(argc, argv);
    walberla::MPIManager::instance()->useWorldComm();

    hyteg::PETScManager petscManager(&argc, &argv);

    const bool writeVTK = false;

    for (uint_t solverType = 0; solverType < 1; solverType++) {
        WALBERLA_LOG_INFO_ON_ROOT("### " << (solverType == 0 ? "PETScCG: " : "HytegCG: ") << " ###");

        hyteg::testEGVectorLaplace2D(solverType, writeVTK, 3, 4);

        hyteg::testEGVectorLaplace3D(solverType, writeVTK, 3, 4);

    }
    return EXIT_SUCCESS;
}

namespace hyteg {


    void testEGVectorLaplace2D(uint_t solverType, bool writeVTK, uint_t minLevel, uint_t maxLevel) {

          WALBERLA_LOG_INFO_ON_ROOT("### Test on single triangle, inhom. BC, rhs != 0 ###");
          {
              std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) { return sin(x[0]) * sin(x[1]); };
              std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &x) { return 2 * sin(x[0]) * sin(x[1]); };
              hyteg::runTestcase("EGVectorLaplaceConvergence_tri_1el",
                                 minLevel,
                                 maxLevel,
                                 MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh"),
                                 std::make_tuple(solFunc, solFunc, solFunc),
                                 std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                                 solverType,
                                 writeVTK);
          }



          WALBERLA_LOG_INFO_ON_ROOT("### Test on single triangle, inhom. BC, rhs = 0 ###");
          {
              std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) { return sin(x[0]) * sinh(x[1]); };
              std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &) { return 0; };

              hyteg::runTestcase("EGVectorLaplaceConvergence_tri_1el",
                                 minLevel,
                                 maxLevel,
                                 MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh"),
                                 std::make_tuple(solFunc, solFunc, solFunc),
                                 std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                                 solverType,
                                 writeVTK);
          }

          WALBERLA_LOG_INFO_ON_ROOT("### Test on 2 triangles, inhom. BC, rhs = 0 ###");
          {
              MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/tri_2el.msh");

              std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                  return sin(pi * x[0]) * sin(pi * x[1]);
              };

              std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &x) {
                  return 2 * pi * pi * sin(pi * x[0]) * sin(pi * x[1]);
              };
              hyteg::runTestcase("EGVectorLaplaceConvergence_tri_2el",
                                 minLevel,
                                 maxLevel,
                                 meshInfo,
                                 std::make_tuple(solFunc, solFunc, solFunc),
                                 std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                                 solverType,
                                 writeVTK);
          }

          WALBERLA_LOG_INFO_ON_ROOT("### Test on quad 2D, inhom. BC, rhs = 0 ###");
          {
              MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");

              std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                  return sin(pi * x[0]) * sin(pi * x[1]);
              };

              std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &x) {
                  return 2 * pi * pi * sin(pi * x[0]) * sin(pi * x[1]);
              };
              hyteg::runTestcase("EGVectorLaplaceConvergence_quad_4el",
                                 minLevel,
                                 maxLevel,
                                 meshInfo,
                                 std::make_tuple(solFunc, solFunc, solFunc),
                                 std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                                 solverType,
                                 writeVTK);
          }


          WALBERLA_LOG_INFO_ON_ROOT("### Test on square, multiple macros, inhom. BC, rhs != 0 ###");
          {
              MeshInfo meshInfo =
                      hyteg::MeshInfo::meshRectangle(Point2D({-1, -1}), Point2D({1, 1}), hyteg::MeshInfo::CRISS, 2, 2);

              std::function<real_t(const hyteg::Point3D &)> solFunc = [](const hyteg::Point3D &x) {
                  return std::exp(-x[0] - (x[1] * x[1]));
              };
              std::function<real_t(const hyteg::Point3D &)> rhsFunc = [](const hyteg::Point3D &x) {
                  return -(4 * x[1] * x[1] - 1) * std::exp(-x[0] - (x[1] * x[1]));
              };
              hyteg::runTestcase("EGVectorLaplaceConvergence_square_expSol",
                                 minLevel,
                                 maxLevel,
                                 meshInfo,
                                 std::make_tuple(solFunc, solFunc, solFunc),
                                 std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                                 solverType,
                                 writeVTK);
          }

    }

    void testEGVectorLaplace3D(uint_t solverType, bool writeVTK, uint_t minLevel, uint_t maxLevel) {




        WALBERLA_LOG_INFO_ON_ROOT("### Test pyramid_2el, first ###");
        {
            MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/3D/pyramid_2el.msh");

            std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                return sin(x[0]) * sin(x[1]) * sin(x[2]);
            };

            std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &x) {
                return 3 * sin(x[0]) * sin(x[1]) * sin(x[2]);
            };

            hyteg::runTestcase("EGVectorLaplaceConvergence_pyramid_2el",
                               minLevel,
                               maxLevel,
                               meshInfo,
                               std::make_tuple(solFunc, solFunc, solFunc),
                               std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                               solverType,
                               writeVTK);
        }

        WALBERLA_LOG_INFO_ON_ROOT("### Test pyramid_2el, second ###");
        {
            MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/3D/pyramid_2el.msh");

            std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                return sin(x[0]) * sinh(x[1]) * x[2];
            };
            std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &) { return 0; };

            hyteg::runTestcase("EGVectorLaplaceConvergence_pyramid_2el",
                               minLevel,
                               maxLevel,
                               meshInfo,
                               std::make_tuple(solFunc, solFunc, solFunc),
                               std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                               solverType,
                               writeVTK);
        }

        WALBERLA_LOG_INFO_ON_ROOT("### Test pyramid_4el, first ###");
        {
            MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/3D/pyramid_4el.msh");

            std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                return sin(pi * x[0]) * sin(pi * x[1]) * sin(pi * x[2]);
            };

            std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &x) {
                return 3 * pi * pi * sin(pi * x[0]) * sin(pi * x[1]) * sin(pi * x[2]);
            };

            hyteg::runTestcase("EGVectorLaplaceConvergence_pyramid_4el",
                               minLevel,
                               maxLevel,
                               meshInfo,
                               std::make_tuple(solFunc, solFunc, solFunc),
                               std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                               solverType,
                               writeVTK);
        }

        WALBERLA_LOG_INFO_ON_ROOT("### Test pyramid_4el, second ###");
        {
            MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/3D/pyramid_4el.msh");

            std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                return sin(x[0]) * sinh(x[1]) * x[2];
            };
            std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &) { return 0; };

            hyteg::runTestcase("EGVectorLaplaceConvergence_pyramid_4el",
                               minLevel,
                               maxLevel,
                               meshInfo,
                               std::make_tuple(solFunc, solFunc, solFunc),
                               std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                               solverType,
                               writeVTK);
        }

        WALBERLA_LOG_INFO_ON_ROOT("### Test on cube, first ###");
        {
            MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), 1, 1, 1);

            std::function<real_t(const Point3D &)> solFunc = [](const Point3D &x) {
                return sin(x[0]) * sinh(x[1]) * x[2];
            };
            std::function<real_t(const Point3D &)> rhsFunc = [](const Point3D &) { return 0; };

            hyteg::runTestcase("EGVectorLaplaceConvergence_cube",
                               minLevel,
                               maxLevel - 1,
                               meshInfo,
                               std::make_tuple(solFunc, solFunc, solFunc),
                               std::make_tuple(rhsFunc, rhsFunc, rhsFunc),
                               solverType,
                               writeVTK);
        }



        WALBERLA_LOG_INFO_ON_ROOT("### Test on cube, split solution ###");
        {
            MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), 1, 1, 1);

            std::function<real_t(const hyteg::Point3D &)> solFunc_u = [](const hyteg::Point3D &xx) {
                return -real_c(4) * std::cos(real_c(4) * xx[2]);
            };
            std::function<real_t(const hyteg::Point3D &)> solFunc_v = [](const hyteg::Point3D &xx) {
                return real_c(8) * std::cos(real_c(8) * xx[0]);
            };
            std::function<real_t(const hyteg::Point3D &)> solFunc_w = [](const hyteg::Point3D &xx) {
                return -real_c(2) * std::cos(real_c(2) * xx[1]);
            };

            std::function<real_t(const Point3D &)> rhsFunc_u = [](const Point3D &xx) {
                return real_c(-64) * std::cos(4 * xx[2]);
            };
            std::function<real_t(const Point3D &)> rhsFunc_v = [](const Point3D &xx) {
                return real_c(512) * std::cos(8 * xx[0]);
            };
            std::function<real_t(const Point3D &)> rhsFunc_w = [](const Point3D &xx) {
                return real_c(-8) * std::cos(2 * xx[1]);
            };

            hyteg::runTestcase("EGVectorLaplaceConvergence_cube",
                               minLevel,
                               maxLevel,
                               meshInfo,
                               std::make_tuple(solFunc_u, solFunc_v, solFunc_w),
                               std::make_tuple(rhsFunc_u, rhsFunc_v, rhsFunc_w),
                               solverType,
                               writeVTK);
        }
    }
} // namespace hyteg