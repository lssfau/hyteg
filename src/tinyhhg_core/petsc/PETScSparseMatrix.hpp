#pragma once

#include "PETScWrapper.hpp"
#include "tinyhhg_core/types/flags.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"

#include "tinyhhg_core/composites/petsc/P1StokesPetsc.hpp"

#include "tinyhhg_core/p2functionspace/P2Petsc.hpp"
#include "tinyhhg_core/composites/petsc/P2P1TaylorHoodPetsc.hpp"

#include "tinyhhg_core/petsc/PETScVector.hpp"

namespace hhg {

template <class OperatorType, template <class> class FunctionType>
class PETScSparseMatrix {
protected:
  Mat mat;
  bool assembled;

public:
  PETScSparseMatrix() = delete;

  PETScSparseMatrix(uint_t localSize, uint_t globalSize, const char name[] = "Mat") {
    MatCreate(walberla::MPIManager::instance()->comm(),&mat);
    MatSetType(mat,MATMPIAIJ);
    MatSetSizes(mat,(PetscInt)localSize,(PetscInt)localSize,(PetscInt)globalSize,(PetscInt)globalSize);
    // Roughly overestimate number of non-zero entries for faster assembly of matrix
    MatMPIAIJSetPreallocation(mat, 1000, NULL, 1000, NULL);
    //MatSetUp(mat);
    //MatCreateAIJ(walberla::MPIManager::instance()->comm(),(PetscInt)localSize,(PetscInt)localSize,(PetscInt)globalSize,(PetscInt)globalSize,7, NULL, 0, NULL,&mat);
    setName(name);
    reset();
  }



  virtual ~PETScSparseMatrix() {
    MatDestroy(&mat);
  }

  inline void createMatrixFromFunction(const OperatorType& op, uint_t level,const FunctionType<PetscInt>& numerator,DoFType flag = All){
    //WALBERLA_LOG_INFO_ON_ROOT("Creating PETSc Matrix")
    hhg::petsc::createMatrix<OperatorType>(op, numerator, numerator, mat, level, flag);

    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

  inline bool createMatrixFromFunctionOnce(const OperatorType& op, uint_t level,const FunctionType<PetscInt>& numerator,DoFType flag = All){
    if(assembled)
      return false;
    createMatrixFromFunction(op,level,numerator,flag);
    assembled = true;
    return true;
  }

  inline void print(const std::string& name) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,name.c_str(),&viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB );
    //PetscViewerMatlabOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);
    MatView(mat,viewer);
    PetscViewerDestroy(&viewer);
  }

  void applyDirichletBC(const FunctionType<PetscInt>& numerator, uint_t level){
    //WALBERLA_LOG_INFO_ON_ROOT("")
    std::vector<PetscInt> ind;
    hhg::petsc::applyDirichletBC(numerator,ind,level);

    MatZeroRows(mat,ind.size(),ind.data(),1.0,0,0);


    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);

  }

    /// \brief Applies Dirichlet BCs to a linear system without losing symmetry.
    ///
    /// Uses the PETSc function MatZeroRowsColumns() which does that automatically.
    /// Still, we need to think how we can easily integrate this to use more efficient
    /// solvers in HyTeG, because the RHS is modified depending on the original system.
    ///
    /// So far I do not know any solution to this without re-assembling the system every time
    /// we solve it since we need to also rebuild the RHS.
    /// It should be possible to store a copy of the original system and circumvent re-assembling by
    /// copying it and applying only MatZeroRowsColumns() (without re-assembly) before calling the solver.
    /// If PETSc is only used as a coarse grid solver this might be a good solution.
    ///
    /// \param dirichletSolution a function that has the respective values interpolated on the Dirichlet boundary
    /// \param numerator an enumerated function
    /// \param rhsVec RHS of the system as PETSc vector - NOTE THAT THIS IS MODIFIED IN PLACE
    /// \param level the refinement level
    ///
    void applyDirichletBCSymmetrically( const FunctionType< real_t >&        dirichletSolution,
                                        const FunctionType< PetscInt >&      numerator,
                                        PETScVector< real_t, FunctionType >& rhsVec,
                                        const uint_t&                        level )
    {
       std::vector< PetscInt > bcIndices;
       hhg::petsc::applyDirichletBC( numerator, bcIndices, level );

       PETScVector< real_t, FunctionType > dirichletSolutionVec( dirichletSolution, numerator, level );

       WALBERLA_ASSERT(
           isSymmetric(),
           "PETSc: Dirichlet boundary conditions can only be applied symmetrically if the original system is symmetric." );

       // This is required as the implementation of MatZeroRowsColumns() checks (for performance reasons?!)
       // if there are zero diagonals in the matrix. If there are, the function halts.
       // To disable that check, we need to allow setting MAT_NEW_NONZERO_LOCATIONS to true.
       MatSetOption( mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE );

       MatZeroRowsColumns( mat, bcIndices.size(), bcIndices.data(), 1.0, dirichletSolutionVec.get(), rhsVec.get() );

       WALBERLA_ASSERT( isSymmetric() );
    }

  inline void reset()  { assembled = false; }

  inline void setName(const char name[]){ PetscObjectSetName((PetscObject)mat,name); }

  inline Mat& get() { return mat; }

  bool isSymmetric(real_t tol = real_c(1e-13)) {
    Mat B;
    PetscReal norm;
    MatTranspose(mat, MAT_INITIAL_MATRIX, &B);
    MatAYPX(B, -1.0, mat, DIFFERENT_NONZERO_PATTERN);
    MatNorm(B, NORM_INFINITY, &norm);
    // WALBERLA_LOG_DEVEL_ON_ROOT( "PETSC_NORM = " << norm );
    MatDestroy( &B );
    return norm < tol;
  }

};

}

#endif
