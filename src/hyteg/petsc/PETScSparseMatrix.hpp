/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/types/flags.hpp"

#include "PETScWrapper.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/composites/petsc/P1StokesPetsc.hpp"
#include "hyteg/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2Petsc.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScSparseMatrixInfo.hpp"
#include "hyteg/petsc/PETScSparseMatrixProxy.hpp"
#include "hyteg/petsc/PETScVector.hpp"

namespace hyteg {

/// Wrapper class for PETSc sparse matrix usage
template < class OperatorType >
class PETScSparseMatrix
{
 public:
   template < typename ValueType >
   using FunctionTypeSrc = typename OperatorType::srcType::template FunctionType< ValueType >;

   template < typename ValueType >
   using FunctionTypeDst = typename OperatorType::dstType::template FunctionType< ValueType >;

   PETScSparseMatrix() = delete;

   PETScSparseMatrix( uint_t          localRows,
                      uint_t          globalRows,
                      uint_t          globalCols,
                      const char      name[]            = "Mat",
                      const MPI_Comm& petscCommunicator = walberla::mpi::MPIManager::instance()->comm() )
   : petscCommunicator_( petscCommunicator )
   , assembled_( false )
   {
      MatCreate( petscCommunicator, &mat );
      MatSetType( mat, MATMPIAIJ );
      MatSetSizes( mat, (PetscInt) localRows, (PetscInt) globalCols, (PetscInt) globalRows, (PetscInt) globalCols );
      // Roughly overestimate number of non-zero entries for faster assembly of matrix
      MatMPIAIJSetPreallocation( mat, 500, NULL, 500, NULL );
      setName( name );
      reset();
   }

   PETScSparseMatrix( uint_t          localSize,
                      uint_t          globalSize,
                      const char      name[]            = "Mat",
                      const MPI_Comm& petscCommunicator = walberla::mpi::MPIManager::instance()->comm() )
   : PETScSparseMatrix( localSize, globalSize, globalSize, name, petscCommunicator )
   {}

   PETScSparseMatrix( const std::shared_ptr< PrimitiveStorage >& storage,
                      const uint_t&                              level,
                      const char                                 name[] = "Mat",
                      const MPI_Comm& petscCommunicator                 = walberla::mpi::MPIManager::instance()->comm() )
   : PETScSparseMatrix( numberOfLocalDoFs< typename OperatorType::dstType::Tag >( *storage, level ),
                        numberOfGlobalDoFs< typename OperatorType::dstType::Tag >( *storage, level, petscCommunicator ),
                        name,
                        petscCommunicator )
   {}

   PETScSparseMatrix( const FunctionTypeSrc< PetscInt >& enumerator,
                      const uint_t&                      level,
                      const char                         name[]            = "Mat",
                      const MPI_Comm&                    petscCommunicator = walberla::mpi::MPIManager::instance()->comm() )
   : PETScSparseMatrix( numberOfLocalDoFs( enumerator, level ),
                        numberOfGlobalDoFs( enumerator, level, petscCommunicator ),
                        name,
                        petscCommunicator )
   {}

   PETScSparseMatrix( const FunctionTypeSrc< PetscInt >& enumeratorSrc,
                      const FunctionTypeDst< PetscInt >& enumeratorDst,
                      const uint_t&                      level,
                      const char                         name[]            = "Mat",
                      const MPI_Comm&                    petscCommunicator = walberla::mpi::MPIManager::instance()->comm() )
   : PETScSparseMatrix( numberOfLocalDoFs( enumeratorDst, level ),
                        numberOfGlobalDoFs( enumeratorDst, level, petscCommunicator ),
                        numberOfGlobalDoFs( enumeratorSrc, level, petscCommunicator ),
                        name,
                        petscCommunicator )
   {}

   virtual ~PETScSparseMatrix() { MatDestroy( &mat ); }

   inline void createMatrixFromOperator( const OperatorType&                op,
                                         uint_t                             level,
                                         const FunctionTypeSrc< matIdx_t >& numerator,
                                         DoFType                            flag = All )
   {
      auto proxy = std::make_shared< PETScSparseMatrixProxy >( mat );
      op.toMatrix( proxy, numerator, numerator, level, flag );

      MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY );
      assembled_ = true;
   }

   inline void createMatrixFromOperator( const OperatorType&                op,
                                         uint_t                             level,
                                         const FunctionTypeSrc< matIdx_t >& numeratorSrc,
                                         const FunctionTypeDst< matIdx_t >& numeratorDst,
                                         DoFType                            flag = All )
   {
      auto proxy = std::make_shared< PETScSparseMatrixProxy >( mat );
      op.toMatrix( proxy, numeratorSrc, numeratorDst, level, flag );

      MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY );
      assembled_ = true;
   }

   inline bool createMatrixFromOperatorOnce( const OperatorType&                op,
                                             uint_t                             level,
                                             const FunctionTypeSrc< matIdx_t >& numerator,
                                             DoFType                            flag = All )
   {
      if ( assembled_ )
         return false;
      createMatrixFromOperator( op, level, numerator, flag );
      return true;
   }

   inline void print( const std::string& name, bool binary = false, PetscViewerFormat format = PETSC_VIEWER_ASCII_MATRIXMARKET )
   {
      PetscViewer viewer;
      if ( binary )
      {
         PetscViewerBinaryOpen( petscCommunicator_, name.c_str(), FILE_MODE_WRITE, &viewer );
      }
      else
      {
         PetscViewerASCIIOpen( petscCommunicator_, name.c_str(), &viewer );
         PetscViewerPushFormat( viewer, format );
      }
      MatView( mat, viewer );
      PetscViewerDestroy( &viewer );
   }

   void applyDirichletBC( const FunctionTypeSrc< PetscInt >& numerator, uint_t level )
   {
      std::vector< PetscInt > ind;
      hyteg::petsc::applyDirichletBC( numerator, ind, level );

      // This is required as the implementation of MatZeroRows() checks (for performance reasons?!)
      // if there are zero diagonals in the matrix. If there are, the function halts.
      // To disable that check, we need to allow setting MAT_NEW_NONZERO_LOCATIONS to true.
      MatSetOption( mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE );

      MatZeroRows( mat, ind.size(), ind.data(), 1.0, 0, 0 );

      MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY );
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
   void applyDirichletBCSymmetrically( const FunctionTypeSrc< real_t >&                                     dirichletSolution,
                                       const FunctionTypeSrc< PetscInt >&                                   numerator,
                                       PETScVector< real_t, OperatorType::dstType::template FunctionType >& rhsVec,
                                       const uint_t&                                                        level )
   {
      std::vector< PetscInt > bcIndices;
      hyteg::petsc::applyDirichletBC( numerator, bcIndices, level );

      PETScVector< real_t, FunctionTypeSrc > dirichletSolutionVec(
          dirichletSolution, numerator, level, All, "dirichletSolutionVec", rhsVec.getCommunicator() );

      // This is required as the implementation of MatZeroRowsColumns() checks (for performance reasons?!)
      // if there are zero diagonals in the matrix. If there are, the function halts.
      // To disable that check, we need to allow setting MAT_NEW_NONZERO_LOCATIONS to true.
      MatSetOption( mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE );

      MatZeroRowsColumns( mat, bcIndices.size(), bcIndices.data(), 1.0, dirichletSolutionVec.get(), rhsVec.get() );
   }

   /// \brief Variant of applyDirichletBCSymmetrically() that only modifies the matrix itself
   ///
   /// \return Vector with global indices of the Dirichlet DoFs
   std::vector< PetscInt > applyDirichletBCSymmetrically( const FunctionTypeSrc< PetscInt >& numerator, const uint_t& level )
   {
      std::vector< PetscInt > bcIndices;
      hyteg::petsc::applyDirichletBC( numerator, bcIndices, level );

      // This is required as the implementation of MatZeroRowsColumns() checks (for performance reasons?!)
      // if there are zero diagonals in the matrix. If there are, the function halts.
      // To disable that check, we need to allow setting MAT_NEW_NONZERO_LOCATIONS to true.
      MatSetOption( mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE );

      MatZeroRowsColumns( mat, bcIndices.size(), bcIndices.data(), 1.0, NULL, NULL );

      return bcIndices;
   }

   inline void reset() { assembled_ = false; }

   /// \brief Sets all entries of the matrix to zero.
   inline void zeroEntries() { MatZeroEntries( mat ); }

   inline void setName( const char name[] ) { PetscObjectSetName( (PetscObject) mat, name ); }

   inline Mat& get() { return mat; }

   bool isSymmetric( real_t tol = real_c( 1e-13 ) )
   {
      Mat       B;
      PetscReal norm;
      MatTranspose( mat, MAT_INITIAL_MATRIX, &B );
      MatAYPX( B, -1.0, mat, DIFFERENT_NONZERO_PATTERN );
      MatNorm( B, NORM_INFINITY, &norm );
      // WALBERLA_LOG_DEVEL_ON_ROOT( "PETSC_NORM = " << norm );
      MatDestroy( &B );
      return norm < tol;
   }

   bool isDiagonal( real_t tol = real_c( 1e-13 ) )
   {
      Mat       B;
      PetscReal norm;
      Vec       diag;
      MatCreate( petscCommunicator_, &B );
      MatSetType( B, MATMPIAIJ );
      PetscInt localSize, globalSize;
      MatGetSize( mat, &localSize, &globalSize );
      MatSetSizes( B, localSize, localSize, globalSize, globalSize );
      MatSetUp( B );
      MatAssemblyBegin( B, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( B, MAT_FINAL_ASSEMBLY );
      VecCreate( petscCommunicator_, &diag );
      VecSetType( diag, VECMPI );
      VecSetSizes( diag, localSize, globalSize );
      VecSetUp( diag );
      MatCopy( mat, B, DIFFERENT_NONZERO_PATTERN );
      MatGetDiagonal( B, diag );
      VecScale( diag, -1.0 );
      MatDiagonalSet( B, diag, ADD_VALUES );
      MatNorm( B, NORM_INFINITY, &norm );
      // WALBERLA_LOG_DEVEL_ON_ROOT( "PETSC_NORM = " << norm );
      MatDestroy( &B );
      VecDestroy( &diag );
      return norm < tol;
   }

   PETScSparseMatrixInfo getInfo()
   {
      if ( !assembled_ )
      {
         WALBERLA_ABORT( "Matrix assembly must be complete before calling getInfo()!" );
      }
      MatInfo info;
      MatGetInfo( mat, MAT_GLOBAL_SUM, &info );
      PETScSparseMatrixInfo matInfo( info );
      return matInfo;
   };

 protected:
   MPI_Comm petscCommunicator_;
   Mat      mat;
   bool     assembled_;
};

} // namespace hyteg

#endif
