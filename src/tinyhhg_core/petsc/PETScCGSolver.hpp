#pragma once

#include <memory>

#include "PETScVector.hpp"
#include "PETScSparseMatrix.hpp"
#include "tinyhhg_core/solvers/Solver.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#if (PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 9)
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverType
#else
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

namespace hhg {


    template<class OperatorType>
    class PETScCGSolver : public Solver<OperatorType> {
    public:

        typedef typename OperatorType::srcType FunctionType;

        PETScCGSolver(const std::shared_ptr<PrimitiveStorage> &storage, const uint_t &level,
                      uint_t iterations = PETSC_DEFAULT, real_t tol = 1e-15)
                : allocatedLevel_(level), num("numerator", storage, level, level),
                  Amat(numberOfLocalDoFs<typename FunctionType::Tag>(*storage, level),
                       numberOfGlobalDoFs<typename FunctionType::Tag>(*storage, level)),
                  CopyMat(numberOfLocalDoFs<typename FunctionType::Tag>(*storage, level),
                          numberOfGlobalDoFs<typename FunctionType::Tag>(*storage, level)),
                  xVec(numberOfLocalDoFs<typename FunctionType::Tag>(*storage, level)),
                  bVec(numberOfLocalDoFs<typename FunctionType::Tag>(*storage, level))
#if 0
                , inKernel( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
#endif
                , flag_(hhg::All) {
            num.enumerate(level);

            maxIter = iterations;
            tolerance = tol;

            KSPCreate(walberla::MPIManager::instance()->comm(), &ksp);
            KSPSetType(ksp, KSPCG);
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCNONE);
            HHG_PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
            KSPSetTolerances(ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxIter);
        }

        ~PETScCGSolver() {
            KSPDestroy(&ksp);
        }

#if 0
        void setNullSpace( FunctionType & inKernel, const uint_t & level )
  {
    inKernel.createVectorFromFunction( inKernel, *num, level, All );
    VecNormalize(inKernel.get(), NULL);
    MatNullSpace nullspace;
    MatNullSpaceCreate( walberla::MPIManager::instance()->comm(), PETSC_FALSE, 1, &(inKernel.get()), &nullspace );
    MatSetNullSpace( Amat.get(), nullspace );
  }
#endif

        void setConstantNullSpace() {
            MatNullSpace nullspace;
            MatNullSpaceCreate(walberla::MPIManager::instance()->comm(), PETSC_TRUE, 0, NULL, &nullspace);
            MatSetNullSpace(Amat.get(), nullspace);
        }

        void solve(const OperatorType &A, const FunctionType &x, const FunctionType &b, const uint_t level) {
            WALBERLA_CHECK_EQUAL(level, allocatedLevel_);

            x.getStorage()->getTimingTree()->start("PETSc CG Solver");

            bVec.createVectorFromFunction(b, num, level, All);

            const bool firstAssembly = Amat.createMatrixFromFunctionOnce(A, level, num, All);

            CopyMat.copyMat(Amat);
            CopyMat.applyDirichletBCSymmetrically(x, num, bVec, level);

            if (firstAssembly) {
                KSPSetOperators(ksp, CopyMat.get(), CopyMat.get());
                KSPSetFromOptions(ksp);
                KSPSetUp(ksp);
            }

            KSPSolve(ksp, bVec.get(), xVec.get());

            xVec.createFunctionFromVector(x, num, level, flag_);

            x.getStorage()->getTimingTree()->stop("PETSc CG Solver");
        }


    private:
        uint_t allocatedLevel_;
        typename OperatorType::srcType::template FunctionType<PetscInt> num;
        PETScSparseMatrix<OperatorType, OperatorType::srcType::template FunctionType> Amat, CopyMat;
        PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> xVec;
        PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> bVec;
#if 0
        PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> inKernel;
#endif

        KSP ksp;
        PC pc;
        hhg::DoFType flag_;
        uint_t maxIter;
        real_t tolerance;
    };
}
#endif
