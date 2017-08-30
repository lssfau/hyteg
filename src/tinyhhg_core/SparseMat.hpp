#pragma once


#include <petscmat.h>
#include "tinyhhg_core/types/flags.hpp"



namespace hhg {

template <typename Operator,typename Function>
class SparseMat {
public:
  Mat mat;

public:
  SparseMat() = delete;

  SparseMat(Operator& op, uint_t level,Function& numerator,uint_t num,DoFType flag = All) {
    MatCreate(walberla::MPIManager::instance()->comm(),&mat);
    MatSetType(mat,MATAIJ);
    MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,num,num);
    MatSetUp(mat);

    op.createMatrix(numerator, numerator, mat, level, flag);

    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

  }

  virtual ~SparseMat() {
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

};


}
