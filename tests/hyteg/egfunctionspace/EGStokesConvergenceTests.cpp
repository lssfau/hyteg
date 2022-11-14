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

using walberla::real_t;
using walberla::uint_t;

// scalar lambda for one component of analytical solution and rhs
typedef std::function<real_t(const hyteg::PointND<real_t, 3> &p)> ScalarLambda;

// tuple of function for solution (u,p) and rhs of vector values stokes equation
typedef std::tuple<ScalarLambda, ScalarLambda, ScalarLambda, ScalarLambda> LambdaTuple;

using hyteg::Point3D;
using hyteg::dg::eg::EGLaplaceOperator;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0ConstEpsilonStokesOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;

namespace hyteg {
    namespace dg {
        namespace eg {
            auto copyBdry = [](EGP0StokesFunction<real_t> fun) {
                fun.p().setBoundaryCondition(fun.uvw().getBoundaryCondition());
            };

            template<typename StokesOperatorType>
            constexpr bool isEGP0Discr() {
                return std::is_same<StokesOperatorType, EGP0EpsilonStokesOperator>::value ||
                       std::is_same<StokesOperatorType, EGP0StokesOperator>::value ||
                       std::is_same<StokesOperatorType, EGP0ConstEpsilonStokesOperator>::value;
            }

            template<typename StokesOperatorType>
            constexpr bool isP2P1Discr() {
                return std::is_same<StokesOperatorType, hyteg::P2P1ElementwiseAffineEpsilonStokesOperator>::value ||
                       std::is_same<StokesOperatorType, hyteg::P2P1TaylorHoodStokesOperator>::value;
            }

            template<typename StokesOperatorType>
            class StokesConvergenceOrderTest {
            public:
                using StokesFunctionType = typename StokesOperatorType::srcType;
                using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType<idx_t>;

                StokesConvergenceOrderTest(const std::string &testName,
                                           LambdaTuple solTuple,
                                           LambdaTuple rhsTuple,
                                           StokesOperatorType &Op,
                                           const std::shared_ptr<PrimitiveStorage> &storage,
                                           const uint_t minLevel,
                                           const uint_t maxLevel,
                                           const uint_t &solverType = 5,
                                           bool writeVTK = false) : testName_(testName), solTuple_(solTuple),
                                                                    rhsTuple_(rhsTuple), Op_(Op), storage_(storage),
                                                                    solverType_(solverType), writeVTK_(writeVTK) {
                    // std::vector<std::tuple<real_t, real_t, real_t> > errors_per_h;
                    WALBERLA_LOG_INFO_ON_ROOT("Running " << testName);
                    auto ratesCounter = EGConvRatesCounter();
                    ratesCounter.printHeader();

                    for (uint_t level = minLevel; level <= maxLevel; level++) {
                        ratesCounter.update(RunStokesTestOnLevel(level));
                        ratesCounter.printCurrentRates();
                    }

                    //const real_t expectedRate = 4.;
                    //WALBERLA_CHECK_LESS( 0.9 * expectedRate, currentRate, "unexpected rate!" );
                    //WALBERLA_CHECK_GREATER( 1.1 * expectedRate, currentRate, "unexpected rate!" );
                    //WALBERLA_LOG_INFO_ON_ROOT( "Test " << testName << " converged correctly." );

                    // write to plot file
                    /*
                      std::ofstream err_file;
                      auto err_file_name = "../../../hyteg-plots/EG_ConvOrders/" + testName;
                      err_file.open(err_file_name);
                      for (auto err: errors_per_h) {
                          err_file << std::get<0>(err) << ", " << std::get<1>(err) << ", " << std::get<2>(err) << "\n";
                      }
                      err_file.close();
                      */
                }

            private:
                class EGConvRatesCounter {
                public:
                    EGConvRatesCounter() : errors_({std::nan(""), std::nan(""), std::nan(""), std::nan("")}),
                                           rates_({std::nan(""), std::nan(""), std::nan(""), std::nan("")}) {
                    }

                    void update(const std::array<real_t, 4> &newErrors) {
                        std::transform(errors_.begin(), errors_.end(),
                                       newErrors.begin(), rates_.begin(),
                                       std::divides<real_t>());
                        errors_ = newErrors;
                    }

                    void printHeader() {
                        WALBERLA_LOG_INFO_ON_ROOT(
                                walberla::format("%6s|%15s|%15s|%15s|%15s|%15s|%15s|%15s|%15s", "level",
                                                 "e_v", "e_v_conf", "e_v_disc", "e_p",
                                                 "rate_v", "rate_v_conf", "rate_v_disc", "rate_p"));

                    }

                    void printCurrentRates() {
                        WALBERLA_LOG_INFO_ON_ROOT(walberla::format(
                                "%6d|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e", level_, errors_[0],
                                errors_[1], errors_[2], errors_[3],
                                rates_[0], rates_[1], rates_[2], rates_[3]));
                    }

                private:
                    // e_v e_v_conf e_v_disc e_p
                    std::array<real_t, 4> errors_;
                    std::array<real_t, 4> rates_;
                    uint_t level_;
                };

                std::array<real_t, 4>
                computeErrorNorms(uint_t level, StokesFunctionType &err, StokesFunctionType &Merr) {
                    real_t e_v, e_v_conf, e_v_disc, e_p;
                    e_v_conf = 0.;
                    e_v_disc = 0.;
                    if constexpr (isEGP0Discr<StokesOperatorType>()) {
                        EGMassOperator M_vel(storage_, level, level);
                        M_vel.apply(err.uvw(), Merr.uvw(), level, Inner, Replace);
                        if (!storage_->hasGlobalCells()) {
                            auto mass_form = std::make_shared<dg::DGMassFormP0P0>();
                            dg::DGOperator M_pressure(storage_, level, level, mass_form);
                            M_pressure.apply(*err.p().getDGFunction(), *Merr.p().getDGFunction(), level, All, Replace);
                            e_p = sqrt(err.p().dotGlobal(Merr.p(), level, Inner));
                        }

                        e_v = sqrt(err.uvw().dotGlobal(Merr.uvw(), level, Inner));
                    } else if constexpr (isP2P1Discr<StokesOperatorType>()) {
                        WALBERLA_UNUSED(e_v_conf);
                        WALBERLA_UNUSED(e_v_disc);
                        P2ConstantMassOperator M_vel(storage_, level, level);
                        if (!storage_->hasGlobalCells()) {
                            M_vel.apply(err.uvw()[0], Merr.uvw()[0], level, Inner, Replace);
                            M_vel.apply(err.uvw()[1], Merr.uvw()[1], level, Inner, Replace);
                        } else {
                            M_vel.apply(err.uvw()[0], Merr.uvw()[0], level, Inner, Replace);
                            M_vel.apply(err.uvw()[1], Merr.uvw()[1], level, Inner, Replace);
                            M_vel.apply(err.uvw()[2], Merr.uvw()[2], level, Inner, Replace);
                        }
                        P1ConstantMassOperator M_pressure(storage_, level, level);
                        M_pressure.apply(err.p(), Merr.p(), level, All, Replace);

                        /*
                        P2toP2QuadraticProlongation P2P2ProlongationOp;

                        P1toP1LinearProlongation P1P1ProlongationOp;
                        for (uint_t k = 0; k < err.uvw().getDimension(); k++) {
                            P2P2ProlongationOp.prolongate(err.uvw()[k], level, Inner);
                        }
                        P1P1ProlongationOp.prolongate(err.p(), level, Inner);
                        //     P2ConstantMassOperator M_vel( storage, level+1, level+1 );
                        //    M_vel.apply( err.uvw()[0], Merr.uvw()[0], level+1, Inner, Replace );
                        //    M_vel.apply( err.uvw()[1], Merr.uvw()[1], level+1, Inner, Replace );
                        discrL2_velocity_err =
                                sqrt(err.uvw().dotGlobal(err.uvw(), level + 1, Inner) /
                                     real_c(numberOfGlobalDoFs(u.uvw(), level + 1)));
*/
                        e_v = sqrt(err.uvw().dotGlobal(Merr.uvw(), level, Inner));
                        e_p = sqrt(err.p().dotGlobal(Merr.p(), level, Inner));
                    }
                    if constexpr (isEGP0Discr<StokesOperatorType>()) {
                        e_v_disc = sqrt(err.uvw().getDiscontinuousPart()->dotGlobal(
                                *err.uvw().getDiscontinuousPart(), level, All) /
                                             real_c(numberOfGlobalDoFs(*err.uvw().getDiscontinuousPart(),
                                                                       level)));
                        e_v_conf = sqrt(err.uvw().getConformingPart()->dotGlobal(
                                *err.uvw().getConformingPart(), level, All) /
                                             real_c(numberOfGlobalDoFs(*err.uvw().getConformingPart(),
                                                                       level)));
                        WALBERLA_LOG_INFO_ON_ROOT(
                                "||e_v_disc|| = " << e_v_disc
                                                  << ", ||e_v_conf|| = "
                                                  << e_v_conf);

                        return {e_v, e_v_conf, e_v_disc, e_p};
                    } else {
                        return {e_v, 0, 0, e_p};
                    }
                    // discrL2_velocity_err = sqrt( err.uvw().dotGlobal( Merr.uvw(), level, Inner ) );
                    //      discrL2_pressure_err = sqrt( err.p().dotGlobal( Merr.p(), level, Inner ) );
                    //  discrL2_velocity_err = sqrt( err.uvw().dotGlobal( err.uvw(), level, All )/ real_c( numberOfGlobalDoFs( u.uvw(), level ) ) );
                    // discrL2_pressure_err = sqrt( err.p().dotGlobal( err.p(), level, All ) / real_c( numberOfGlobalDoFs( u.p(), level ) ));


                }

                void setupRHSandBC(uint_t level, const StokesFunctionType &f, StokesFunctionType &rhs,
                                   StokesFunctionType &u) {

                    // solution, rhs as a lambda function
                    auto [u_x_expr, u_y_expr, u_z_expr, p_expr] = solTuple_;
                    auto [f_x_expr, f_y_expr, f_z_expr, g_expr] = rhsTuple_;

                    if constexpr (isP2P1Discr<StokesOperatorType>()) {
                        P2ConstantMassOperator M_vel(storage_, level, level);
                        if (!storage_->hasGlobalCells()) {
                            M_vel.apply(f.uvw()[0], rhs.uvw()[0], level, All);
                            M_vel.apply(f.uvw()[1], rhs.uvw()[1], level, All);
                            u.uvw().interpolate({u_x_expr, u_y_expr}, level, hyteg::DirichletBoundary);
                        } else {
                            M_vel.apply(f.uvw()[0], rhs.uvw()[0], level, All);
                            M_vel.apply(f.uvw()[1], rhs.uvw()[1], level, All);
                            M_vel.apply(f.uvw()[2], rhs.uvw()[2], level, All);
                            u.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, hyteg::DirichletBoundary);

                            rhs.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
                        }

                        P1ConstantMassOperator M_pressure(storage_, level, level);
                        M_pressure.apply(f.p(), rhs.p(), level, All, Replace);
                    } else {
                        if constexpr (isEGP0Discr<StokesOperatorType>()) {
                            EGMassOperator M_vel(storage_, level, level);
                            M_vel.apply(f.uvw(), rhs.uvw(), level, All, Replace);

                            if (!storage_->hasGlobalCells()) {
                                u.uvw().getConformingPart()->interpolate({u_x_expr, u_y_expr}, level,
                                                                         DirichletBoundary);
                                auto mass_form = std::make_shared<dg::DGMassFormP0P0>();
                                dg::DGOperator M_pressure(storage_, level, level, mass_form);
                                M_pressure.apply(*f.p().getDGFunction(), *rhs.p().getDGFunction(), level, All, Replace);
                            } else {
                                u.uvw().getConformingPart()->interpolate({u_x_expr, u_y_expr, u_z_expr}, level,
                                                                         DirichletBoundary);
                                rhs.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);

                            }
                        } else {
                            WALBERLA_ABORT("Benchmark not implemented for other discretizations!");
                        }
                    }
                }

                std::array<real_t, 4> RunStokesTestOnLevel(
                        const uint_t &level
                ) {
                    StokesFunctionNumeratorType numerator("numerator", storage_, level, level);
                    numerator.enumerate(level);
                    WALBERLA_LOG_INFO_ON_ROOT("Global DoFs: " << numberOfGlobalDoFs(numerator, level) << ", u DoFs: "
                                                              << numberOfGlobalDoFs(numerator.uvw(), level)
                                                              << ", p DoFs: "
                                                              << numberOfGlobalDoFs(numerator.p(), level));

                    // solution, rhs as a lambda function
                    auto [u_x_expr, u_y_expr, u_z_expr, p_expr] = solTuple_;
                    auto [f_x_expr, f_y_expr, f_z_expr, g_expr] = rhsTuple_;

                    StokesFunctionType u("u", storage_, level, level);
                    StokesFunctionType f("f", storage_, level, level);
                    StokesFunctionType rhs("rhs", storage_, level, level);
                    StokesFunctionType sol("sol", storage_, level, level);
                    StokesFunctionType err("err", storage_, level, level + 1);
                    StokesFunctionType Merr("Merr", storage_, level, level + 1);
                    if constexpr (isEGP0Discr<StokesOperatorType>()) {
                        copyBdry(u);
                        copyBdry(f);
                        copyBdry(rhs);
                        copyBdry(sol);
                        copyBdry(err);
                        copyBdry(Merr);
                    }

                    // interpolate analytical solution and rhs
                    if (!storage_->hasGlobalCells()) {
                        sol.uvw().interpolate({u_x_expr, u_y_expr}, level, All);
                        f.uvw().interpolate({f_x_expr, f_y_expr}, level, All);
                    } else {
                        sol.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, All);
                        f.uvw().interpolate({f_x_expr, f_y_expr, f_z_expr}, level, Inner);
                    }
                    sol.p().interpolate(p_expr, level, All);
                    f.p().interpolate(g_expr, level, All);


                    setupRHSandBC(level, f, rhs, u);

                    // solve
                    switch (solverType_) {
                        case 0: {
                            MinResSolver<StokesOperatorType> solver(storage_, level, level);
                            solver.setPrintInfo(true);
                            solver.solve(Op_, u, rhs, level);
                            break;
                        }
                        default: {
                            PETScMinResSolver<StokesOperatorType> solver(storage_, level);
                            solver.setFromOptions(true);
                            /*StokesFunctionType nullSpace("ns", storage, level, level);
                            nullSpace.uvw().interpolate(0, level, All);
                            nullSpace.p().interpolate(1, level, All);
                            solver.setNullSpace(nullSpace);*/
                            solver.solve(Op_, u, rhs, level);
                            break;
                        }
                    }

                    // pressure projection to space of mean-value-0-functions
                    if constexpr (isEGP0Discr<StokesOperatorType>()) {
                        hyteg::dg::projectMean(u.p(), level);
                        hyteg::dg::projectMean(sol.p(), level);
                    } else if constexpr (isP2P1Discr<StokesOperatorType>()) {
                        hyteg::vertexdof::projectMean(u.p(), level);
                        hyteg::vertexdof::projectMean(sol.p(), level);
                    }

                    err.assign({1.0, -1.0}, {u, sol}, level, All);

                    if (writeVTK_) {
                        if constexpr (isEGP0Discr<StokesOperatorType>()) {
                            VTKOutput vtk("/mnt/c/Users/Fabia/OneDrive/Desktop/hyteg_premerge/hyteg-build/output",
                                          testName_, storage_);
                            vtk.add(u);
                            vtk.add(*u.uvw().getConformingPart());
                            vtk.add(*u.uvw().getDiscontinuousPart());

                            vtk.add(sol);
                            vtk.add(*sol.uvw().getConformingPart());
                            vtk.add(*sol.uvw().getDiscontinuousPart());

                            vtk.add(err);
                            vtk.add(*err.uvw().getConformingPart());
                            vtk.add(*err.uvw().getDiscontinuousPart());

                            vtk.add(f);

                            vtk.write(level, 1);
                        }
                    }

                    //real_t h = MeshQuality::getMaximalEdgeLength(storage, level);
                    return computeErrorNorms(level, err, Merr);
                }

                std::string testName_;
                LambdaTuple solTuple_;
                LambdaTuple rhsTuple_;
                StokesOperatorType Op_;
                std::shared_ptr<PrimitiveStorage> storage_;
                uint_t solverType_;
                bool writeVTK_;


            };

            void
            Stokes2D(const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr<PrimitiveStorage> &storage);

            void
            Epsilon2D(const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr<PrimitiveStorage> &storage);

            void SmoothViscosityTest2D(const uint_t minLevel, const uint_t maxLevel,
                                       const std::shared_ptr<PrimitiveStorage> &storage);

            void HytegSolverCheck2D(const uint_t minLevel, const uint_t maxLevel,
                                    const std::shared_ptr<PrimitiveStorage> &storage);

            void
            Stokes3D(const uint_t minLevel, const uint_t maxLevel);
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
    uint_t maxLevel = 4;

    {
        WALBERLA_LOG_INFO_ON_ROOT("### Testing 2D ###")
        auto meshInfo = hyteg::MeshInfo::meshRectangle(
                hyteg::Point2D({-1, -1}), hyteg::Point2D({1, 1}), hyteg::MeshInfo::CRISSCROSS, 1, 1);
        hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                  walberla::uint_c(
                                                          walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);


        hyteg::dg::eg::Stokes2D(minLevel, maxLevel, storage);
        // hyteg::dg::eg::Epsilon2D( minLevel, maxLevel, storage );
        // hyteg::dg::eg::SmoothViscosityTest2D( minLevel, maxLevel, storage );
    }


    hyteg::dg::eg::Stokes3D(minLevel, maxLevel);


    return 0;
}

namespace hyteg {
    namespace dg {
        namespace eg {

            void
            Stokes2D(const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr<PrimitiveStorage> &storage) {
                WALBERLA_LOG_INFO_ON_ROOT("### Stokes2D ###")
                EGP0StokesOperator EGP0StokesOp(storage, minLevel, maxLevel);
                auto dummyLambda = [](const Point3D &) -> real_t { return 0; };

                StokesConvergenceOrderTest<EGP0StokesOperator>(
                        "EGP0StokesOp2D_hom_asym",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * (x + 1) / 2) * std::sin(M_PI * (y + 1) / 2) * std::exp(y);
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * (x + 1) / 2) * std::sin(M_PI * (y + 1) / 2) * std::exp(y);
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

                                    const real_t x0 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x1 = std::sin(x0);
                                    const real_t x2 = std::exp(y);
                                    const real_t x3 = std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0));
                                    const real_t x4 = x2 * x3;
                                    return (1.0 / 2.0) * std::pow(M_PI, 2) * x1 * x2 * x3 - x1 * x4 -
                                           M_PI * x4 * std::cos(x0) + 1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x1 = std::sin(x0);
                                    const real_t x2 = std::exp(y);
                                    const real_t x3 = std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0));
                                    const real_t x4 = x2 * x3;
                                    return (1.0 / 2.0) * std::pow(M_PI, 2) * x1 * x2 * x3 - x1 * x4 -
                                           M_PI * x4 * std::cos(x0) + 1;
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
                        EGP0StokesOp,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0StokesOperator>(
                        "EGP0StokesOp2D_hom_pointsym",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * (x + 1) / 2) * std::sin(M_PI * (y + 1) / 2);
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * (x + 1) / 2) * std::sin(M_PI * (y + 1) / 2);
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
                                    return (1.0 / 2.0) * std::pow(M_PI, 2) *
                                           std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) *
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) +
                                           1;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return (1.0 / 2.0) * std::pow(M_PI, 2) *
                                           std::sin(M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0)) *
                                           std::sin(M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0)) +
                                           1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * ((1.0 / 2.0) * y + 1.0 / 2.0);
                                    const real_t x1 = M_PI * ((1.0 / 2.0) * x + 1.0 / 2.0);
                                    const real_t x2 = M_PI_2;
                                    return -x2 * std::sin(x0) * std::cos(x1) - x2 * std::sin(x1) * std::cos(x0);
                                }),
                        EGP0StokesOp,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0StokesOperator>(
                        "EGP0StokesOp2D_hom_asym_p",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * (x + y));
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * (x + y));
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return 2 * x - y + 4;
                                }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * x;
                                    const real_t x1 = std::sin(x0);
                                    const real_t x2 = M_PI * (x + y);
                                    const real_t x3 = std::pow(M_PI, 2);
                                    const real_t x4 = M_PI * y;
                                    const real_t x5 = x3 * std::sin(x4);
                                    const real_t x6 = 2 * std::cos(x2);
                                    return -x1 * x3 * x6 * std::cos(x4) + 4 * x1 * x5 * std::sin(x2) -
                                           x5 * x6 * std::cos(x0) + 2;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = std::pow(M_PI, 2);
                                    const real_t x1 = M_PI * x;
                                    const real_t x2 = std::sin(x1);
                                    const real_t x3 = M_PI * y;
                                    const real_t x4 = std::sin(x3);
                                    const real_t x5 = M_PI * (x + y);
                                    const real_t x6 = 2 * std::cos(x5);
                                    return 4 * x0 * x2 * x4 * std::sin(x5) - x0 * x2 * x6 * std::cos(x3) -
                                           x0 * x4 * x6 * std::cos(x1) - 1;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    const real_t x0 = M_PI * y;
                                    const real_t x1 = M_PI * x;
                                    const real_t x2 = std::sin(x1);
                                    const real_t x3 = M_PI * (x + y);
                                    const real_t x4 = M_PI * std::sin(x3);
                                    const real_t x5 = std::sin(x0);
                                    return -x2 * x4 * std::cos(x0) - 2 * M_PI * x2 * x5 * std::cos(x3) -
                                           x4 * x5 * std::cos(x1);
                                }),
                        EGP0StokesOp,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0StokesOperator>("EGP0StokesOp2D_inhom_asym",
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
                                                               EGP0StokesOp,
                                                               storage,
                                                               minLevel,
                                                               maxLevel);
            }

            void
            Epsilon2D(const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr<PrimitiveStorage> &storage) {
                WALBERLA_LOG_INFO_ON_ROOT("### Epsilon2D ###")
                EGP0ConstEpsilonStokesOperator EGP0ConstantEpsilonOp(storage, minLevel, maxLevel);
                EGP0EpsilonStokesOperator EGP0EpsilonOp_mu_1(storage, minLevel, maxLevel,
                                                             [](const hyteg::Point3D &) { return 1; });

                auto dummyLambda = [](const Point3D &) -> real_t { return 0; };
                StokesConvergenceOrderTest<EGP0ConstEpsilonStokesOperator>(
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
                        EGP0ConstantEpsilonOp,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0ConstEpsilonStokesOperator>(
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
                        EGP0ConstantEpsilonOp,
                        storage,
                        minLevel,
                        maxLevel);

                StokesConvergenceOrderTest<EGP0ConstEpsilonStokesOperator>(
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
                        EGP0ConstantEpsilonOp,
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

                /*
   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_1(
       storage, minLevel, maxLevel, []( const hyteg::Point3D& ) { return 1; } );
   hyteg::P2P1TaylorHoodStokesOperator P2P1StokesOp( storage, minLevel, maxLevel );

   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_smooth(
       storage, minLevel, maxLevel, []( const hyteg::Point3D& p ) {
          const real_t x = p[0];
          const real_t y = p[1];
          return std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) +
                 1;
       } );
   hyteg::StokesConvergenceOrderTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >(
       "THEpsilonInhomogeneousDirichlet:smoothViscosity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( y ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( y ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0  = std::exp( y );
              const real_t x1  = std::pow( M_PI, 2 );
              const real_t x2  = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x3  = std::sin( x2 );
              const real_t x4  = std::exp( x );
              const real_t x5  = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x6  = std::sin( x5 );
              const real_t x7  = x4 * x6;
              const real_t x8  = x3 * x7;
              const real_t x9  = x8 + 1;
              const real_t x10 = std::cos( x2 );
              const real_t x11 = x0 * x10;
              const real_t x12 = std::cos( x5 );
              const real_t x13 = 0.25 * M_PI;
              const real_t x14 = 0.5 * x0 * x3 + x11 * x13;
              return 0.5 * x0 * x1 * x3 * x9 - 1.0 * M_PI * x11 * ( ( 1.0 / 2.0 ) * M_PI * x10 * x7 + x8 ) -
                     M_PI * x12 * x3 * x4 * ( x12 * x13 + x14 ) - 2 * x9 * ( -0.125 * x1 * x6 + x14 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x   = p[0];
              const real_t y   = p[1];
              const real_t x0  = std::exp( x );
              const real_t x1  = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2  = std::sin( x1 );
              const real_t x3  = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4  = M_PI * std::cos( x3 );
              const real_t x5  = std::exp( y );
              const real_t x6  = x2 * x5;
              const real_t x7  = 1.0 * x6;
              const real_t x8  = std::pow( M_PI, 2 );
              const real_t x9  = std::sin( x3 );
              const real_t x10 = x0 * x9;
              const real_t x11 = x10 * x2;
              const real_t x12 = 2 * x11 + 2;
              const real_t x13 = std::cos( x1 );
              const real_t x14 = M_PI * x13;
              return -x0 * x2 * x4 * ( 0.5 * x4 + x7 ) - x12 * ( x7 - 0.25 * x8 * x9 ) -
                     x12 * ( 0.25 * M_PI * x13 * x5 - 0.125 * x6 * x8 ) -
                     2 * ( ( 1.0 / 2.0 ) * x10 * x14 + x11 ) * ( 0.25 * x14 * x5 + 0.25 * x4 + 0.5 * x6 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI_2;
              const real_t x1 = std::exp( y );
              const real_t x2 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              return -x0 * x1 * std::cos( x2 ) - x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) - x1 * std::sin( x2 );
           } ),
       P2P1ElementwiseEpsilonOp_mu_smooth,
       storage,
       3,
       5,
       false );

   hyteg::StokesConvergenceOrderTest< hyteg::P2P1TaylorHoodStokesOperator >(
       "P2P1StokesHomogeneousDirichlet_asymmetric_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = x2 * std::sin( x3 );
              const real_t x5 = M_PI_2;
              return -x1 * x2 * x5 * std::cos( x3 ) - x1 * x4 - x4 * x5 * std::cos( x0 );
           } ),
       P2P1StokesOp,
       storage,
       3,
       5,
       false );

   hyteg::StokesConvergenceOrderTest< hyteg::P2P1TaylorHoodStokesOperator >(
       "P2P1StokesHomogeneousDirichlet:simple_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * x );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::sin( M_PI * y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::pow( M_PI, 2 ) * std::sin( M_PI * x ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::pow( M_PI, 2 ) * std::sin( M_PI * y ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -M_PI * std::cos( M_PI * x ) - M_PI * std::cos( M_PI * y );
           } ),
       P2P1StokesOp,
       storage,
       3,
       5,
       false );

   hyteg::StokesConvergenceOrderTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >(
       "P2P1EpsilonHomogeneousDirichlet_asymmetric_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = x2 * std::sin( x3 );
              const real_t x5 = M_PI_2;
              return -x1 * x2 * x5 * std::cos( x3 ) - x1 * x4 - x4 * x5 * std::cos( x0 );
           } ),
       P2P1ElementwiseEpsilonOp_mu_1,
       storage,
       3,
       5,
       false );
   */
            }

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

            void SmoothViscosityTest2D(const uint_t minLevel, const uint_t maxLevel,
                                       const std::shared_ptr<PrimitiveStorage> &storage) {
                WALBERLA_LOG_INFO_ON_ROOT("### Epsilon2D with varying but smooth viscosity ###")
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

            void HytegSolverCheck2D(const uint_t minLevel, const uint_t maxLevel,
                                    const std::shared_ptr<PrimitiveStorage> &storage) {
                WALBERLA_LOG_INFO_ON_ROOT("### Using Hyteg MINRES ###")
                EGP0StokesOperator EGP0StokesOp(storage, minLevel, maxLevel);
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
                        "EGP0EpsilonOp_divFree_smoothVisc",
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
                        maxLevel,
                        1);

                StokesConvergenceOrderTest<EGP0StokesOperator>("EGP0StokesOp_inhom_asym",
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
                                                               EGP0StokesOp,
                                                               storage,
                                                               minLevel,
                                                               maxLevel,
                                                               1);
            }

            void Stokes3D(const uint_t minLevel, const uint_t maxLevel) {
                WALBERLA_LOG_INFO_ON_ROOT("### Stokes3D ###")
                {

                    //WALBERLA_LOG_INFO_ON_ROOT("### tet_1el, hom. ###");
                    auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/tet_1el.msh");
                    hyteg::SetupPrimitiveStorage setupStorage(meshInfo,
                                                              walberla::uint_c(
                                                                      walberla::mpi::MPIManager::instance()->numProcesses()));
                    setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
                    auto storage = std::make_shared<hyteg::PrimitiveStorage>(setupStorage, 1);

                    EGP0StokesOperator EGP0StokesOp(storage, minLevel, maxLevel);
                    hyteg::P2P1TaylorHoodStokesOperator P2P1StokesOp(storage, minLevel, maxLevel);

                    WALBERLA_LOG_INFO_ON_ROOT("### tet_1el, inhom. ###");
                    StokesConvergenceOrderTest<EGP0StokesOperator>(
                            "EGP0StokesOp3D",
                            std::make_tuple([](const Point3D &xx) -> real_t {
                                                return -real_c(4) * std::cos(real_c(4) * xx[2]);
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                return real_c(8) * std::cos(real_c(8) * xx[0]);
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                            }),
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        return 4 * std::sin(8 * xx[1]) * std::sin(2 * xx[2]) * std::cos(4 * xx[0]) -
                                               64 * std::cos(4 * xx[2]);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        return 8 * std::sin(4 * xx[0]) * std::sin(2 * xx[2]) * std::cos(8 * xx[1]) +
                                               512 * std::cos(8 * xx[0]);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        return 2 * std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::cos(2 * xx[2]) -
                                               8 * std::cos(2 * xx[1]);
                                    },
                                    [](const Point3D &) -> real_t { return 0; }),
                            EGP0StokesOp,
                            storage,
                            minLevel,
                            maxLevel,
                            2,
                            true);

                    StokesConvergenceOrderTest<hyteg::P2P1TaylorHoodStokesOperator>(
                            "P2P1StokesOp3D",
                            std::make_tuple([](const Point3D &xx) -> real_t {
                                                return -real_c(4) * std::cos(real_c(4) * xx[2]);
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                return real_c(8) * std::cos(real_c(8) * xx[0]);
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
                                            }),
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        return 4 * std::sin(8 * xx[1]) * std::sin(2 * xx[2]) * std::cos(4 * xx[0]) -
                                               64 * std::cos(4 * xx[2]);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        return 8 * std::sin(4 * xx[0]) * std::sin(2 * xx[2]) * std::cos(8 * xx[1]) +
                                               512 * std::cos(8 * xx[0]);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        return 2 * std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::cos(2 * xx[2]) -
                                               8 * std::cos(2 * xx[1]);
                                    },
                                    [](const Point3D &) -> real_t { return 0; }),
                            P2P1StokesOp,
                            storage,
                            minLevel,
                            maxLevel,
                            2);
                    /*
                    WALBERLA_LOG_INFO_ON_ROOT("### tet_1el, hom. ###");
                    StokesConvergenceOrderTest<EGP0StokesOperator>(
                            "EGP0StokesOp3D",
                            std::make_tuple([](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                const real_t x0 = 2 * x;
                                                const real_t x1 = 2 * y;
                                                const real_t x2 = 2 * z;
                                                return std::sin(M_PI * x0) * std::sin(M_PI * x1) * std::sin(M_PI * x2) *
                                                       std::sin(M_PI * (x0 + x1 + x2 - 2));
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                const real_t x0 = 2 * x;
                                                const real_t x1 = 2 * y;
                                                const real_t x2 = 2 * z;
                                                return std::sin(M_PI * x0) * std::sin(M_PI * x1) * std::sin(M_PI * x2) *
                                                       std::sin(M_PI * (x0 + x1 + x2 - 2));
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                const real_t x0 = 2 * x;
                                                const real_t x1 = 2 * y;
                                                const real_t x2 = 2 * z;
                                                return std::sin(M_PI * x0) * std::sin(M_PI * x1) * std::sin(M_PI * x2) *
                                                       std::sin(M_PI * (x0 + x1 + x2 - 2));
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                return std::sin(4 * x) * std::sin(8 * y) * std::sin(2 * z);
                                            }),
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*z;
                                        const real_t x1 = std::pow(M_PI, 2);
                                        const real_t x2 = 2*x;
                                        const real_t x3 = M_PI*x2;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = 2*y;
                                        const real_t x6 = M_PI*x5;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x0;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*(x0 + x2 + x5 - 2);
                                        const real_t x11 = 8*std::cos(x10);
                                        const real_t x12 = x1*x11*x4;
                                        return -x1*x11*x7*x9*std::cos(x3) + 24*x1*x4*x7*x9*std::sin(x10) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6) + 4*std::sin(x0)*std::sin(8*y)*std::cos(4*x);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*z;
                                        const real_t x1 = std::pow(M_PI, 2);
                                        const real_t x2 = 2*x;
                                        const real_t x3 = M_PI*x2;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = 2*y;
                                        const real_t x6 = M_PI*x5;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x0;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*(x0 + x2 + x5 - 2);
                                        const real_t x11 = 8*std::cos(x10);
                                        const real_t x12 = x1*x11*x4;
                                        return -x1*x11*x7*x9*std::cos(x3) + 24*x1*x4*x7*x9*std::sin(x10) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6) + 8*std::sin(4*x)*std::sin(x0)*std::cos(8*y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*z;
                                        const real_t x1 = std::pow(M_PI, 2);
                                        const real_t x2 = 2*x;
                                        const real_t x3 = M_PI*x2;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = 2*y;
                                        const real_t x6 = M_PI*x5;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x0;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*(x0 + x2 + x5 - 2);
                                        const real_t x11 = 8*std::cos(x10);
                                        const real_t x12 = x1*x11*x4;
                                        return -x1*x11*x7*x9*std::cos(x3) + 24*x1*x4*x7*x9*std::sin(x10) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6) + 2*std::sin(4*x)*std::sin(8*y)*std::cos(x0);
                                    },
                                    [](const Point3D & xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*x;
                                        const real_t x1 = 2*y;
                                        const real_t x2 = 2*z;
                                        const real_t x3 = M_PI*(x0 + x1 + x2 - 2);
                                        const real_t x4 = M_PI*x0;
                                        const real_t x5 = std::sin(x4);
                                        const real_t x6 = M_PI*x1;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x2;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*x7*x9;
                                        const real_t x11 = 2*std::sin(x3);
                                        const real_t x12 = M_PI*x11*x5;
                                        return -x10*x11*std::cos(x4) - 6*x10*x5*std::cos(x3) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6);
                                    }),
                            EGP0StokesOp,
                            storage,
                            minLevel,
                            maxLevel,
                            2,
                            true);

                    StokesConvergenceOrderTest<hyteg::P2P1TaylorHoodStokesOperator>(
                            "P2P1StokesOp3D",
                            std::make_tuple([](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                const real_t x0 = 2 * x;
                                                const real_t x1 = 2 * y;
                                                const real_t x2 = 2 * z;
                                                return std::sin(M_PI * x0) * std::sin(M_PI * x1) * std::sin(M_PI * x2) *
                                                       std::sin(M_PI * (x0 + x1 + x2 - 2));
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                const real_t x0 = 2 * x;
                                                const real_t x1 = 2 * y;
                                                const real_t x2 = 2 * z;
                                                return std::sin(M_PI * x0) * std::sin(M_PI * x1) * std::sin(M_PI * x2) *
                                                       std::sin(M_PI * (x0 + x1 + x2 - 2));
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                const real_t x0 = 2 * x;
                                                const real_t x1 = 2 * y;
                                                const real_t x2 = 2 * z;
                                                return std::sin(M_PI * x0) * std::sin(M_PI * x1) * std::sin(M_PI * x2) *
                                                       std::sin(M_PI * (x0 + x1 + x2 - 2));
                                            },
                                            [](const Point3D &xx) -> real_t {
                                                const real_t x = xx[0];
                                                const real_t y = xx[1];
                                                const real_t z = xx[2];
                                                return std::sin(4 * x) * std::sin(8 * y) * std::sin(2 * z);
                                            }),
                            std::make_tuple(
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*z;
                                        const real_t x1 = std::pow(M_PI, 2);
                                        const real_t x2 = 2*x;
                                        const real_t x3 = M_PI*x2;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = 2*y;
                                        const real_t x6 = M_PI*x5;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x0;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*(x0 + x2 + x5 - 2);
                                        const real_t x11 = 8*std::cos(x10);
                                        const real_t x12 = x1*x11*x4;
                                        return -x1*x11*x7*x9*std::cos(x3) + 24*x1*x4*x7*x9*std::sin(x10) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6) + 4*std::sin(x0)*std::sin(8*y)*std::cos(4*x);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*z;
                                        const real_t x1 = std::pow(M_PI, 2);
                                        const real_t x2 = 2*x;
                                        const real_t x3 = M_PI*x2;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = 2*y;
                                        const real_t x6 = M_PI*x5;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x0;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*(x0 + x2 + x5 - 2);
                                        const real_t x11 = 8*std::cos(x10);
                                        const real_t x12 = x1*x11*x4;
                                        return -x1*x11*x7*x9*std::cos(x3) + 24*x1*x4*x7*x9*std::sin(x10) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6) + 8*std::sin(4*x)*std::sin(x0)*std::cos(8*y);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*z;
                                        const real_t x1 = std::pow(M_PI, 2);
                                        const real_t x2 = 2*x;
                                        const real_t x3 = M_PI*x2;
                                        const real_t x4 = std::sin(x3);
                                        const real_t x5 = 2*y;
                                        const real_t x6 = M_PI*x5;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x0;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*(x0 + x2 + x5 - 2);
                                        const real_t x11 = 8*std::cos(x10);
                                        const real_t x12 = x1*x11*x4;
                                        return -x1*x11*x7*x9*std::cos(x3) + 24*x1*x4*x7*x9*std::sin(x10) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6) + 2*std::sin(4*x)*std::sin(8*y)*std::cos(x0);
                                    },
                                    [](const Point3D &xx) -> real_t {
                                        const real_t x = xx[0];
                                        const real_t y = xx[1];
                                        const real_t z = xx[2];
                                        const real_t x0 = 2*x;
                                        const real_t x1 = 2*y;
                                        const real_t x2 = 2*z;
                                        const real_t x3 = M_PI*(x0 + x1 + x2 - 2);
                                        const real_t x4 = M_PI*x0;
                                        const real_t x5 = std::sin(x4);
                                        const real_t x6 = M_PI*x1;
                                        const real_t x7 = std::sin(x6);
                                        const real_t x8 = M_PI*x2;
                                        const real_t x9 = std::sin(x8);
                                        const real_t x10 = M_PI*x7*x9;
                                        const real_t x11 = 2*std::sin(x3);
                                        const real_t x12 = M_PI*x11*x5;
                                        return -x10*x11*std::cos(x4) - 6*x10*x5*std::cos(x3) - x12*x7*std::cos(x8) - x12*x9*std::cos(x6);
                                    }),
                            P2P1StokesOp,
                            storage,
                            minLevel,
                            maxLevel,
                            2);
                    */
/*
     StokesConvergenceOrderTest< EGP0StokesOperator >(
         "EGP0StokesOp3D",
         std::make_tuple( []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                          []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                          []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                          []( const Point3D& xx ) -> real_t {
                             return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                          } ),
         std::make_tuple(
             []( const Point3D& xx ) -> real_t {
                return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
             },
             []( const Point3D& xx ) -> real_t {
                return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
             },
             []( const Point3D& xx ) -> real_t {
                return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
             },
             []( const Point3D& ) -> real_t { return 0; } ),
         EGP0StokesOp,
         storage,
         minLevel,
         maxLevel,
         2,
         true );

     StokesConvergenceOrderTest< hyteg::P2P1TaylorHoodStokesOperator >(
         "P2P1StokesOp3D",
         std::make_tuple( []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                          []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                          []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                          []( const Point3D& xx ) -> real_t {
                             return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                          } ),
         std::make_tuple(
             []( const Point3D& xx ) -> real_t {
                return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
             },
             []( const Point3D& xx ) -> real_t {
                return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
             },
             []( const Point3D& xx ) -> real_t {
                return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
             },
             []( const Point3D& ) -> real_t { return 0; } ),
         P2P1StokesOp,
         storage,
         minLevel,
         maxLevel,
         2 );*/
                }
            }

        } // namespace eg
    } // namespace dg
} // namespace hyteg
