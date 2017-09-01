#pragma once

#include "PETScWrapper.hpp"
#include "tinyhhg_core/types/flags.hpp"

#ifdef HHG_BUILD_WITH_PETSC

namespace hhg {

template <typename OperatorType,typename FunctionType>
class PETScSparseMatrix {
public:
  Mat mat;

public:
  PETScSparseMatrix() = delete;

  PETScSparseMatrix(const char name[],OperatorType& op, uint_t level,FunctionType& numerator,uint_t num,DoFType flag = All) {
    MatCreate(walberla::MPIManager::instance()->comm(),&mat);
    MatSetType(mat,MATAIJ);
    MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,(PetscInt)num,(PetscInt)num);
    MatSetUp(mat);
    op.createMatrix(numerator, numerator, mat, level, flag);

    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    PetscObjectSetName((PetscObject)mat,name);
  }

  virtual ~PETScSparseMatrix() {
    MatDestroy(&mat);
  }

  void print(const char name[]) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB );
    //PetscViewerMatlabOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);
    MatView(mat,viewer);
    PetscViewerDestroy(&viewer);
  }

  void applyDirichletBC(FunctionType& numerator, uint_t level){
    //WALBERLA_LOG_INFO_ON_ROOT("")
    std::vector<PetscInt> ind;
    numerator.applyDirichletBC(ind,level);

    MatZeroRows(mat,ind.size(),ind.data(),1.0,0,0);


    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);

  }





};

}

#endif
