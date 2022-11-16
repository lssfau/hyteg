
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

#pragma once

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/composites/P1P0StokesOperator.hpp"
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

namespace hyteg {
    namespace dg {
        namespace eg {
            auto copyBdry = [](EGP0StokesFunction<real_t> fun) {
                fun.p().setBoundaryCondition(fun.uvw().getBoundaryCondition());
            };

            // scalar lambda for one component of analytical solution and rhs
            typedef std::function<real_t(const hyteg::PointND<real_t, 3> &p)> ScalarLambda;

            // tuple of function for solution (u,p) and rhs of vector values stokes equation
            typedef std::tuple<ScalarLambda, ScalarLambda, ScalarLambda, ScalarLambda> LambdaTuple;

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
            constexpr bool isP1P0Discr() {
                return std::is_same<StokesOperatorType, hyteg::P1P0StokesOperator>::value;
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
                        ratesCounter.printCurrentRates(level);
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

                    void printCurrentRates(uint_t level) {
                        WALBERLA_LOG_INFO_ON_ROOT(walberla::format(
                                "%6d|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e|%15.2e", level, errors_[0],
                                errors_[1], errors_[2], errors_[3],
                                rates_[0], rates_[1], rates_[2], rates_[3]));
                    }

                private:
                    // e_v e_v_conf e_v_disc e_p
                    std::array<real_t, 4> errors_;
                    std::array<real_t, 4> rates_;
                };

                std::array<real_t, 4>
                computeErrorNorms(uint_t level, StokesFunctionType &err, StokesFunctionType &Merr) {
                    real_t e_v, e_v_conf, e_v_disc, e_p;
                    e_v_conf = 0.;
                    e_v_disc = 0.;
                    if constexpr (isEGP0Discr<StokesOperatorType>()) {
                        EGMassOperator M_vel(storage_, level, level);
                        M_vel.apply(err.uvw(), Merr.uvw(), level, All, Replace);
                        auto mass_form = std::make_shared<dg::P0P0MassForm>();
                        dg::DGOperator M_pressure(storage_, level, level, mass_form);
                        M_pressure.apply(*err.p().getDGFunction(), *Merr.p().getDGFunction(), level, All, Replace);
                        e_p = sqrt(err.p().dotGlobal(Merr.p(), level, All));
                        e_v = sqrt(err.uvw().dotGlobal(Merr.uvw(), level, All));
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
                        e_v = sqrt(err.uvw().dotGlobal(Merr.uvw(), level, All));
                        e_p = sqrt(err.p().dotGlobal(Merr.p(), level, All));
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
                        return {e_v, e_v_conf, e_v_disc, e_p};
                    } else {
                        return {e_v, -1, -1, e_p};
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
                            u.uvw().interpolate({u_x_expr, u_y_expr}, level, DirichletBoundary);
                            rhs.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
                        } else {
                            M_vel.apply(f.uvw()[0], rhs.uvw()[0], level, All);
                            M_vel.apply(f.uvw()[1], rhs.uvw()[1], level, All);
                            M_vel.apply(f.uvw()[2], rhs.uvw()[2], level, All);
                            u.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
                            rhs.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
                        }

                        P1ConstantMassOperator M_pressure(storage_, level, level);
                        M_pressure.apply(f.p(), rhs.p(), level, All, Replace);
                    } else if constexpr (isEGP0Discr<StokesOperatorType>()) {
                        EGMassOperator M_vel(storage_, level, level);
                        M_vel.apply(f.uvw(), rhs.uvw(), level, All, Replace);

                        if (!storage_->hasGlobalCells()) {
                            u.uvw().getConformingPart()->interpolate({u_x_expr, u_y_expr}, level,
                                                                     DirichletBoundary);
                            rhs.uvw().interpolate({u_x_expr, u_y_expr}, level, DirichletBoundary);
                        } else {
                            u.uvw().getConformingPart()->interpolate({u_x_expr, u_y_expr, u_z_expr}, level,
                                                                     DirichletBoundary);
                            rhs.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
                        }

                        auto mass_form = std::make_shared<dg::P0P0MassForm>();
                        dg::DGOperator M_pressure(storage_, level, level, mass_form);
                        M_pressure.apply(*f.p().getDGFunction(), *rhs.p().getDGFunction(), level, All, Replace);
                    } else if constexpr (isP1P0Discr<StokesOperatorType>()) {
                        if (!storage_->hasGlobalCells()) {
                            u.uvw().interpolate({u_x_expr, u_y_expr}, level, DirichletBoundary);
                        } else {
                            u.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, DirichletBoundary);
                        }
                    } else {
                        WALBERLA_ABORT("Benchmark not implemented for other discretizations!");
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
                    /*
                     if constexpr (isEGP0Discr<StokesOperatorType>()) {
                         copyBdry(u);
                         copyBdry(f);
                         copyBdry(rhs);
                         copyBdry(sol);
                         copyBdry(err);
                         copyBdry(Merr);
                     }
                     */

                    // interpolate analytical solution and rhs
                    if (!storage_->hasGlobalCells()) {
                        sol.uvw().interpolate({u_x_expr, u_y_expr}, level, All);
                        f.uvw().interpolate({f_x_expr, f_y_expr}, level, All);
                    } else {
                        sol.uvw().interpolate({u_x_expr, u_y_expr, u_z_expr}, level, All);
                        f.uvw().interpolate({f_x_expr, f_y_expr, f_z_expr}, level, All);
                    }
                    sol.p().interpolate(p_expr, level, All);
                    f.p().interpolate(g_expr, level, All);

                    setupRHSandBC(level, f, rhs, u);

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


        } // namespace eg
    } // namespace dg
} // namespace hyteg


