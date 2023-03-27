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

#include <memory>

#include "hyteg/functions/FunctionIterator.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"
#include "hyteg/volumedofspace/VolumeDoFPackInfo.hpp"

#include "PETScSparseMatrix.hpp"
#include "PETScVector.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

    // check for data member blockPrec
    template<typename, typename = void>
    constexpr bool hasBlockPrec = false;
    template<typename T>
    constexpr bool hasBlockPrec<T, std::void_t<decltype(std::declval<T>().blockPrec)> > = true;

    template<class OperatorType>
    class PETScBlockPreconditionedStokesSolver : public Solver<OperatorType> {
    public:
        typedef typename OperatorType::srcType FunctionType;
        typedef typename OperatorType::BlockPreconditioner_T BlockPreconditioner_T;


        /// \brief PETSc-based block preconditioned Krylov solver for the Stokes problem.
        ///
        /// \param velocityPreconditionerType choose from different velocity preconditioners:
        ///                                   - 0: PCGAMG
        ///                                   - 1: PCJACOBI
        ///                                   - 2: Schur complement
        ///                                   - 3: Hypre (BoomerAMG)
        ///                                   - 4: none
        ///                                   - 5: GKB
        /// \param pressurePreconditionerType choose from different pressure preconditioners:
        ///                                   - 0: none
        ///                                   - 1: PCJACOBI (lumped mass)
        /// \param krylovSolverType choose from different solvers:
        ///                                   - 0: MINRES
        ///                                   - 1: FGMRES
        PETScBlockPreconditionedStokesSolver(const std::shared_ptr<PrimitiveStorage> &storage,
                                             const uint_t &level,
                                             const real_t tolerance = 1e-12,
                                             const PetscInt maxIterations = std::numeric_limits<PetscInt>::max(),
                                             const uint_t &velocityPreconditionerType = 1,
                                             const uint_t &pressurePreconditionerType = 1,
                                             const uint_t &krylovSolverType = 0)
                : allocatedLevel_(level), petscCommunicator_(storage->getSplitCommunicatorByPrimitiveDistribution()),
                  num("numerator", storage, level, level), Amat("Amat", petscCommunicator_),
                  AmatNonEliminatedBC("AmatNonEliminatedBC", petscCommunicator_), Pmat("Pmat", petscCommunicator_),
                  xVec("xVec", petscCommunicator_), bVec("bVec", petscCommunicator_),
                  nullspaceVec_("nullspaceVec", petscCommunicator_), storage_(storage), tolerance_(tolerance),
                  maxIterations_(maxIterations), flag_(hyteg::All), nullSpaceSet_(false),
                  velocityPreconditionerType_(velocityPreconditionerType),
                  pressurePreconditionerType_(pressurePreconditionerType), krylovSolverType_(krylovSolverType),
                  verbose_(false), reassembleMatrix_(true), matrixWasAssembledOnce_(false),
                  disableApplicationBC_(false), setFromOptions_(false) {
            //num.enumerate(level);
        }

        ~PETScBlockPreconditionedStokesSolver() = default;

        /// \brief If set to true, the operator is reassembled for every solve / manual assembly call.
        ///        Default is true.
        void reassembleMatrix(bool reassembleMatrix) { reassembleMatrix_ = reassembleMatrix; }

        void disableApplicationBC(bool dis) { disableApplicationBC_ = dis; }

        void setNullSpace(const FunctionType &nullspace) {
            nullSpaceSet_ = true;
            nullspaceVec_.createVectorFromFunction(nullspace, num, allocatedLevel_);
            MatNullSpaceCreate(petscCommunicator_, PETSC_FALSE, 1, &nullspaceVec_.get(), &nullspace_);
        }

        void setTolerance(real_t tol) { tolerance_ = tol; }

        void setMaxIt(uint_t maxIt) { maxIterations_ = maxIt; }

        void setVerbose(bool verbose) { verbose_ = verbose; }

        void setFromOptions(bool doIt) { setFromOptions_ = doIt; }

        void solve(const OperatorType &A, const FunctionType &x, const FunctionType &b, const uint_t level) {
            WALBERLA_CHECK_EQUAL(level, allocatedLevel_)

            walberla::WcTimer timer;

            x.getStorage()->getTimingTree()->start("PETSc block prec MinRes Solver");

            x.getStorage()->getTimingTree()->start("Setup");
            timer.start();

            num.copyBoundaryConditionFromFunction(x);
            num.enumerate(level);

            x.getStorage()->getTimingTree()->start("Index set setup");

            // gather index sets to split matrix into block matrix
            // therefore we need the row indices of the velocity and pressure
            std::vector<PetscInt> velocityIndices;
            std::vector<PetscInt> pressureIndices;

            gatherIndices(velocityIndices, pressureIndices, *storage_, level, num);

            std::sort(velocityIndices.begin(), velocityIndices.end());
            std::sort(pressureIndices.begin(), pressureIndices.end());
            ISCreateGeneral(petscCommunicator_,
                            static_cast< PetscInt >( velocityIndices.size()),
                            velocityIndices.data(),
                            PETSC_COPY_VALUES,
                            &is_[0]);
            ISCreateGeneral(petscCommunicator_,
                            static_cast< PetscInt >( pressureIndices.size()),
                            pressureIndices.data(),
                            PETSC_COPY_VALUES,
                            &is_[1]);

            x.getStorage()->getTimingTree()->stop("Index set setup");

            KSPCreate(petscCommunicator_, &ksp);

            switch (krylovSolverType_) {
                case 0:
                    KSPSetType(ksp, KSPMINRES);
                    break;
                case 1:
                    KSPSetType(ksp, KSPFGMRES);
                    break;
                case 2:
                    KSPSetType(ksp, KSPPREONLY);
                    break;
                default: WALBERLA_ABORT("Invalid solver type for PETSc block prec MinRes solver.")
            }

            KSPSetTolerances(ksp, 1e-30, tolerance_, PETSC_DEFAULT, maxIterations_);
            KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
            KSPSetFromOptions(ksp);

            x.getStorage()->getTimingTree()->start("Vector copy");
            xVec.createVectorFromFunction(x, num, level);
            bVec.createVectorFromFunction(b, num, level, All);
            x.getStorage()->getTimingTree()->stop("Vector copy");

            if (reassembleMatrix_ || !matrixWasAssembledOnce_) {
                x.getStorage()->getTimingTree()->start("Matrix assembly");
                Amat.zeroEntries();
                Pmat.zeroEntries();
                AmatNonEliminatedBC.zeroEntries();
                Amat.createMatrixFromOperator(A, level, num, All);
                AmatNonEliminatedBC.createMatrixFromOperator(A, level, num, All);

                // some operators offer their own blockpreconditioner as member:
                // if the operator has a variable viscosity, the blockPrec cannot be constructed here
                // but must be constructed with the construction of the operator, lablese only then the viscosity function is available
                if constexpr (hasBlockPrec<OperatorType>) {
                    Pmat.createMatrixFromOperator(A.blockPrec, level, num, All);
                } else {
                    BlockPreconditioner_T blockPrec(storage_, level, level);
                    Pmat.createMatrixFromOperator( blockPrec, level, num, All );
                }

                x.getStorage()->getTimingTree()->stop("Matrix assembly");

                x.getStorage()->getTimingTree()->start("Dirichlet BCs");
                if (!disableApplicationBC_) {
                    Amat.applyDirichletBCSymmetrically(x, num, bVec, level);
                    Pmat.applyDirichletBCSymmetrically(num, level);
                }
                x.getStorage()->getTimingTree()->stop("Dirichlet BCs");

                matrixWasAssembledOnce_ = true;
            } else {
                MatCopy(AmatNonEliminatedBC.get(), Amat.get(), SAME_NONZERO_PATTERN);
                x.getStorage()->getTimingTree()->start("Dirichlet BCs");
                if (!disableApplicationBC_) {
                    Amat.applyDirichletBCSymmetrically(x, num, bVec, level);
                }
                x.getStorage()->getTimingTree()->stop("Dirichlet BCs");
            }

            if (nullSpaceSet_) {
                MatSetNullSpace(Amat.get(), nullspace_);
            }


/*
       std::string Pname;

      Pname.append("FOR_MATLAB_PETSCBLOCKPREC_ORDERS_P");
      Pname.append(std::to_string(level));
      Pname.append(".m");
      Pmat.print(Pname, false,PETSC_VIEWER_ASCII_MATLAB);
      WALBERLA_ABORT("byse.")*/
            KSPSetOperators(ksp, Amat.get(), Pmat.get());

            if (velocityPreconditionerType_ == 2) {
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCFIELDSPLIT);
                PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
                //  PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG);
                PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);
                PCFieldSplitSetIS(pc, "0", is_[0]);
                PCFieldSplitSetIS(pc, "1", is_[1]);

                PetscInt numSubKsps;

                PCSetUp(pc);
                PCFieldSplitSchurGetSubKSP(pc, &numSubKsps, &sub_ksps_);

                KSPSetType(sub_ksps_[0], KSPPREONLY);
                KSPSetType(sub_ksps_[1], KSPPREONLY);
                PC pc_u, pc_p;
                KSPGetPC(sub_ksps_[0], &pc_u);
                KSPGetPC(sub_ksps_[1], &pc_p);
                PCSetType(pc_u, KSPCG);
                PCSetType(pc_p, KSPCG);
                KSPSetTolerances(sub_ksps_[0], 1e-15, 1e-15, PETSC_DEFAULT, maxIterations_);
                KSPSetTolerances(sub_ksps_[1], 1e-15, 1e-15, PETSC_DEFAULT, maxIterations_);
            } else if (velocityPreconditionerType_ == 5) {
                // Original system matrix A is used for the GKB preconditioner
                KSPSetOperators(ksp, Amat.get(), Amat.get());

                // preconditioner setup
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCFIELDSPLIT);
                PCFieldSplitSetType(pc, PC_COMPOSITE_GKB);
                PCFieldSplitSetIS(pc, "0", is_[0]);
                PCFieldSplitSetIS(pc, "1", is_[1]);

                // parameters of GKB
                PCFieldSplitSetGKBDelay(pc, 5);
                PCFieldSplitSetGKBMaxit(pc, maxIterations_);
                PCFieldSplitSetGKBNu(pc, 0);
                PCFieldSplitSetGKBTol(pc, tolerance_);
                PCSetUp(pc);

                // one SubKsp: solver for M
                PetscInt numSubKsps;
                PCFieldSplitGetSubKSP(pc, &numSubKsps, &sub_ksps_);

                // CG for M system
                KSPSetType(sub_ksps_[0], KSPCG);
                KSPSetTolerances(sub_ksps_[0], tolerance_ / 10, tolerance_ / 10, PETSC_DEFAULT, maxIterations_);
                PC H_pc;
                KSPGetPC(sub_ksps_[0], &H_pc);
                PCSetType(H_pc, PCNONE);
            } else {
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCFIELDSPLIT);
                PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
                PCFieldSplitSetIS(pc, "u", is_[0]);
                PCFieldSplitSetIS(pc, "p", is_[1]);

                PetscInt numSubKsps;
                PC pc_u, pc_p;
                PCFieldSplitGetSubKSP(pc, &numSubKsps, &sub_ksps_);

                KSPGetPC(sub_ksps_[0], &pc_u);
                KSPGetPC(sub_ksps_[1], &pc_p);

                switch (velocityPreconditionerType_) {
                    case 0:
                        PCSetType(pc_u, PCGAMG);
                        PCGAMGSetType(pc_u, PCGAMGAGG);
                        PCGAMGSetNSmooths(pc_u, 1);
                        break;
                    case 1:
                        PCSetType(pc_u, PCJACOBI);
                        break;
                    case 3: {
                        PCSetType(pc_u, PCHYPRE);
                        break;
                    }
                    case 4:
                        PCSetType(pc_u, PCNONE);
                        break;
                    case 6:
                        KSPSetType(sub_ksps_[0], KSPCG);
                        KSPSetTolerances(sub_ksps_[0], 1e-8, 1e-8, PETSC_DEFAULT, std::numeric_limits<PetscInt>::max());
                        PCSetType(pc_u, PCGAMG);
                        break;
                    default: WALBERLA_ABORT("Invalid velocity preconditioner for PETSc block prec MinRes solver.")
                }

                switch (pressurePreconditionerType_) {
                    case 0:
                        PCSetType(pc_p, PCNONE);
                        break;
                    case 1:
                        // inv. lumped mass
                        PCSetType(pc_p, PCJACOBI);
                        break;
                    case 2:
                        KSPSetType(sub_ksps_[1], KSPPREONLY);
                        KSPSetTolerances(sub_ksps_[0], 1e-6, 1e-6, PETSC_DEFAULT, std::numeric_limits<PetscInt>::max());
                        PCSetType(pc_p, PCLU);
                        break;
                    default: WALBERLA_ABORT("Invalid pressure preconditioner for PETSc block prec MinRes solver.")
                }
            }

            timer.end();
            const double hytegToPetscSetup = timer.last();
            x.getStorage()->getTimingTree()->stop("Setup");
            x.getStorage()->getTimingTree()->start("Solve");

            timer.start();
            KSPSolve(ksp, bVec.get(), xVec.get());
            timer.end();
            const double petscKSPTimer = timer.last();

            x.getStorage()->getTimingTree()->stop("Solve");

            if (verbose_) {
                PetscInt numKSPIterations;
                KSPGetIterationNumber(ksp, &numKSPIterations);
                WALBERLA_LOG_INFO_ON_ROOT("[PETScBlockPreconditionedStokesSolver] num KSP iterations: "
                                                  << numKSPIterations << " | "
                                                  << "PETSc KSPSolver time: " << petscKSPTimer << " | "
                                                  << "HyTeG to PETSc setup: " << hytegToPetscSetup)
            }

            xVec.createFunctionFromVector(x, num, level, flag_);

            x.getStorage()->getTimingTree()->stop("PETSc block prec MinRes Solver");

            PetscFree(sub_ksps_);

            ISDestroy(&is_[0]);
            ISDestroy(&is_[1]);
            KSPDestroy(&ksp);
            if (nullSpaceSet_)
                MatNullSpaceDestroy(&nullspace_);
        }

    private:
        void gatherIndices(std::vector<PetscInt> &velocityIndices,
                           std::vector<PetscInt> &pressureIndices,
                           const PrimitiveStorage &storage,
                           uint_t level,
                           const P1StokesFunction<idx_t> &numerator) {
            for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.uvw()[0], level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }

            for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.uvw()[1], level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }

            if (storage.hasGlobalCells()) {
                for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.uvw()[2], level)) {
                    velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
                }
            }
            for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.p(), level)) {
                pressureIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
        }

        void gatherIndices(std::vector<PetscInt> &velocityIndices,
                           std::vector<PetscInt> &pressureIndices,
                           const PrimitiveStorage &storage,
                           uint_t level,
                           const P2P1TaylorHoodFunction<idx_t> &numerator) {
            for (auto dof:
                    FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.uvw()[0].getVertexDoFFunction(),
                                                                           level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
            for (auto dof: FunctionIterator<EdgeDoFFunction<idx_t> >(numerator.uvw()[0].getEdgeDoFFunction(), level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
            for (auto dof:
                    FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.uvw()[1].getVertexDoFFunction(),
                                                                           level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
            for (auto dof: FunctionIterator<EdgeDoFFunction<idx_t> >(numerator.uvw()[1].getEdgeDoFFunction(), level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
            if (storage.hasGlobalCells()) {
                for (auto dof:
                        FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(
                                numerator.uvw()[2].getVertexDoFFunction(), level)) {
                    velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
                }
                for (auto dof: FunctionIterator<EdgeDoFFunction<idx_t> >(numerator.uvw()[2].getEdgeDoFFunction(),
                                                                         level)) {
                    velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
                }
            }
            for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(numerator.p(), level)) {
                pressureIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
        }

        void gatherIndices(std::vector<PetscInt> &velocityIndices,
                           std::vector<PetscInt> &pressureIndices,
                           const PrimitiveStorage &storage,
                           uint_t level,
                           const EGP0StokesFunction<idx_t> &numerator) {
            for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(
                    numerator.uvw().getConformingPart()->component(0), level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
            for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(
                    numerator.uvw().getConformingPart()->component(1), level)) {
                velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
            }
            if (storage.hasGlobalCells()) {
                for (auto dof: FunctionIterator<vertexdof::VertexDoFFunction<idx_t> >(
                        numerator.uvw().getConformingPart()->component(2), level)) {
                    velocityIndices.push_back(static_cast< PetscInt >( dof.value()));
                }
            }

            gatherP0Indices(velocityIndices, storage, level, *(numerator.uvw().getDiscontinuousPart()));
            gatherP0Indices(pressureIndices, storage, level, numerator.p());

            /*
            for ( auto dof : FunctionIterator< P0Function< idx_t > >( *(numerator.uvw().getDiscontinuousPart()), level ) )
            {
               velocityIndices.push_back( static_cast< PetscInt >( dof.value() ) );
            }
            for ( auto dof : FunctionIterator< P0Function< idx_t > >( numerator.p(), level ) )
            {
               pressureIndices.push_back( static_cast< PetscInt >( dof.value() ) );
            }
             */
        }

        void gatherP0Indices(std::vector<PetscInt> &indices,
                             const PrimitiveStorage &storage,
                             uint_t level,
                             const P0Function<idx_t> &numerator) {
            if (storage.hasGlobalCells()) {
                for (const auto &cellIt: storage.getCells()) {
                    const auto cellId = cellIt.first;
                    const auto cell = *cellIt.second;

                    const auto mem = numerator.getDGFunction()->volumeDoFFunction()->dofMemory(cellId, level);
                    const auto layout = numerator.getDGFunction()->volumeDoFFunction()->memoryLayout();

                    for (auto cellType: celldof::allCellTypes) {
                        for (auto elementIdx: celldof::macrocell::Iterator(level, cellType)) {
                            const auto idx = index(elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, 0, 1,
                                                   level, layout);
                            indices.push_back(static_cast< PetscInt >( mem[idx] ));

                        }
                    }
                }
            } else {
                for (const auto &faceIt: storage.getFaces()) {
                    const auto faceId = faceIt.first;
                    const auto face = *faceIt.second;

                    const auto mem = numerator.getDGFunction()->volumeDoFFunction()->dofMemory(faceId, level);
                    const auto layout = numerator.getDGFunction()->volumeDoFFunction()->memoryLayout();

                    for (auto faceType: facedof::allFaceTypes) {
                        for (auto idxIt: facedof::macroface::Iterator(level, faceType)) {
                            const auto idx = volumedofspace::indexing::index(idxIt.x(), idxIt.y(), faceType, 0, 1,
                                                                             level,
                                                                             layout);//index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, layout );
                            indices.push_back(static_cast< PetscInt >( mem[idx] ));

                        }

                    }
                }
            }
        }

        uint_t allocatedLevel_;
        MPI_Comm petscCommunicator_;
        typename OperatorType::srcType::template FunctionType<idx_t> num;
        PETScSparseMatrix<OperatorType> Amat;
        PETScSparseMatrix<OperatorType> AmatNonEliminatedBC;
        PETScSparseMatrix<BlockPreconditioner_T> Pmat;
        PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> xVec;
        PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> bVec;
        PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> nullspaceVec_;

        std::shared_ptr<PrimitiveStorage> storage_;

        real_t tolerance_;
        idx_t maxIterations_;


        KSP ksp;
        KSP *sub_ksps_;
        PC pc;
        IS is_[2];
        MatNullSpace nullspace_;
        DoFType flag_;
        bool nullSpaceSet_;

        uint_t velocityPreconditionerType_;
        uint_t pressurePreconditionerType_;
        uint_t krylovSolverType_;

        bool disableApplicationBC_;

        bool verbose_;
        bool reassembleMatrix_;
        bool matrixWasAssembledOnce_;
        bool setFromOptions_;
    };

} // namespace hyteg

#endif
