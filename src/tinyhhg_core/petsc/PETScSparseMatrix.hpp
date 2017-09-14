#pragma once

#include "PETScWrapper.hpp"
#include "tinyhhg_core/types/flags.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"

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
    // Usually, we have about 7 nz entries in each matrix row, except in macro-vertex rows.
    // Therefore, we can assume following preallocation for a huge performance gain.
    MatMPIAIJSetPreallocation(mat, 7, NULL, 6, NULL);
    setName(name);
    reset();
  }



  virtual ~PETScSparseMatrix() {
    MatDestroy(&mat);
  }

  inline void createMatrixFromFunction(OperatorType& op, uint_t level,FunctionType<PetscInt>& numerator,DoFType flag = All){
    //WALBERLA_LOG_INFO_ON_ROOT("Creating PETSc Matrix")
    hhg::petsc::createMatrix<OperatorType>(op, numerator, numerator, mat, level, flag);

    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

  inline bool createMatrixFromFunctionOnce(OperatorType& op, uint_t level,FunctionType<PetscInt>& numerator,DoFType flag = All){
    if(assembled)
      return false;
    createMatrixFromFunction(op,level,numerator,flag);
    assembled = true;
    return true;
  }

  inline void print(const char name[]) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB );
    //PetscViewerMatlabOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);
    MatView(mat,viewer);
    PetscViewerDestroy(&viewer);
  }

  void applyDirichletBC(FunctionType<PetscInt>& numerator, uint_t level){
    //WALBERLA_LOG_INFO_ON_ROOT("")
    std::vector<PetscInt> ind;
    hhg::petsc::applyDirichletBC(numerator,ind,level);

    MatZeroRows(mat,ind.size(),ind.data(),1.0,0,0);


    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);

  }

  inline void reset()  { assembled = false; }

  inline void setName(const char name[]){ PetscObjectSetName((PetscObject)mat,name); }

  inline Mat& get() { return mat; }





};

}

#endif
