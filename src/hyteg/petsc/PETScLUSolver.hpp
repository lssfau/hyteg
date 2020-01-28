/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "PETScVector.hpp"
#include "PETScSparseMatrix.hpp"
#include "hyteg/solvers/Solver.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

#if (PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 9)
#define HYTEG_PCFactorSetMatSolverType PCFactorSetMatSolverType
#else
#define HYTEG_PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

namespace hyteg {


template < class OperatorType >
class PETScLUSolver : public Solver< OperatorType >
{
public:

  typedef typename OperatorType::srcType FunctionType;

  PETScLUSolver( const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level )
  : allocatedLevel_( level )
  , num( "numerator", storage, level, level )
  , Amat( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ), numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level ) )
  , xVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
  , bVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
#if 0
  , inKernel( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
#endif
  , flag_( hyteg::All )
  {
    num.enumerate( level );
    KSPCreate( walberla::MPIManager::instance()->comm(), &ksp );
    KSPSetType( ksp, KSPPREONLY );
    KSPSetFromOptions( ksp );
  }

  ~PETScLUSolver(){
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

  void setConstantNullSpace()
  {
    MatNullSpace nullspace;
    MatNullSpaceCreate( walberla::MPIManager::instance()->comm(), PETSC_TRUE, 0, NULL, &nullspace );
    MatSetNullSpace( Amat.get(), nullspace );
  }

  void solve(const OperatorType& A,const FunctionType& x, const FunctionType& b, const uint_t level)
  {
    WALBERLA_CHECK_EQUAL( level, allocatedLevel_ );

    x.getStorage()->getTimingTree()->start( "PETSc LU Solver" );

    b.assign({1.0}, {x}, level, DirichletBoundary);

    bVec.createVectorFromFunction(b,num,level,All);

    x.getStorage()->getTimingTree()->start( "Matrix assembly and solver setup" );
    const bool matrixAssembledForTheFirstTime = Amat.createMatrixFromOperatorOnce(A, level, num, All);

    if( matrixAssembledForTheFirstTime )
    {
      Amat.applyDirichletBC(num,level);
      KSPSetOperators(ksp,Amat.get(),Amat.get());

      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      HYTEG_PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
      //PCFactorSetUpMatSolverPackage(pc); /* call MatGetFactor() to create F */
      //PCFactorGetMatrix(pc,&F);
    }
    x.getStorage()->getTimingTree()->stop( "Matrix assembly and solver setup" );

    x.getStorage()->getTimingTree()->start( "Solver" );
    KSPSolve(ksp,bVec.get(),xVec.get());
    x.getStorage()->getTimingTree()->stop( "Solver" );

    xVec.createFunctionFromVector(x,num,level,flag_);

    x.getStorage()->getTimingTree()->stop( "PETSc LU Solver" );
  }




private:
  uint_t allocatedLevel_;
  typename OperatorType::srcType::template FunctionType<PetscInt> num;
  PETScSparseMatrix<OperatorType,OperatorType::srcType::template FunctionType> Amat;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> xVec;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> bVec;
#if 0
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> inKernel;
#endif

  KSP ksp;
  PC pc;
  hyteg::DoFType flag_;
  //Mat F; //factored Matrix



};



}

#endif
