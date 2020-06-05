/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/composites/petsc/P1StokesPetsc.hpp"
#include "hyteg/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "hyteg/elementwiseoperators/ElementwiseOperatorPetsc.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2Petsc.hpp"
#include "hyteg/trilinos/KokkosWrapper.hpp"
#include "hyteg/trilinos/TeuchosWrapper.hpp"
#include "hyteg/trilinos/TpetraWrapper.hpp"
#include "hyteg/trilinos/TrilinosSparseMatrixProxy.hpp"
#include "hyteg/trilinos/TrilinosVector.hpp"

namespace hyteg {
namespace trilinos {

using walberla::real_t;
using walberla::uint_t;

using Teuchos::RCP;
using Teuchos::rcp;

template < typename OperatorType, template < class > class FunctionType, typename MatrixScalarType = real_t >
class TrilinosSparseMatrix
{
 public:
   typedef Tpetra::Map<>       MapType;
   typedef Tpetra::CrsMatrix<> MatrixType;

   /// \brief Allocates the parallel sparse matrix data structure.
   ///
   /// This constructor may only pre-allocate memory, not even the non-zero
   /// structure is determined. A subsequent call to assembleSystem() is necessary.
   ///
   TrilinosSparseMatrix( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t& level )
   : level_( level )
   {
      trilinosCommunicatorRaw_ = walberla::mpi::MPIManager::instance()->comm();
      trilinosCommunicator_    = rcp( new Teuchos::MpiComm< int >( trilinosCommunicatorRaw_ ) );

      const uint_t numGlobalUnknowns =
          numberOfGlobalDoFs< typename OperatorType::dstType::Tag >( *storage, level, trilinosCommunicatorRaw_ );
      const size_t maxEntriesPerRow = 100; // only a hint, not restriction

      rowMap_    = rcp( new MapType( Tpetra::global_size_t( numGlobalUnknowns ), 0, trilinosCommunicator_ ) );
      crsMatrix_ = rcp( new MatrixType( rowMap_, maxEntriesPerRow ) );
   }

   /// \brief Assembles the sparse matrix from the passed operator.
   ///
   /// This routine assembles the "Neumann" system. To incorporate Dirichlet boundary
   /// conditions, subsequent calls to applyDirichletBoundaryConditions*() must follow.
   void assembleSystem( const OperatorType& op, const FunctionType< PetscInt >& numerator )
   {
      if ( crsMatrix_->isFillComplete() )
      {
         crsMatrix_->resumeFill();
      }
      crsMatrix_->setAllToScalar( 0 );
      auto proxy = std::make_shared< TrilinosSparseMatrixProxy >( crsMatrix_ );
      hyteg::petsc::createMatrix( op, numerator, numerator, proxy, level_, All );
      crsMatrix_->fillComplete();
   }

   /// \brief Modifies the system matrix to conform to the underlying Dirichlet problem.
   ///
   /// Sets the diagonal entry of all rows that correspond to DoFs on a Dirichlet boundary to 1.
   void applyDirichletBoundaryConditions( const FunctionType< PetscInt >& numerator )
   {
      if ( crsMatrix_->isFillComplete() )
      {
         crsMatrix_->resumeFill();
      }
      std::vector< PetscInt > dirichletRowIndices;
      hyteg::petsc::applyDirichletBC( numerator, dirichletRowIndices, level_ );

      for ( auto row : dirichletRowIndices )
      {
         Teuchos::ArrayView< const MatrixType::global_ordinal_type > columnIndices;
         crsMatrix_->getCrsGraph()->getGlobalRowView( row, columnIndices );
         for ( auto col : columnIndices )
         {
            real_t val = 0;
            if ( row == col )
            {
               val = 1.0;
            }
            crsMatrix_->replaceGlobalValues( row,
                                             Teuchos::tuple< Tpetra::Vector<>::global_ordinal_type >( col ),
                                             Teuchos::tuple< Tpetra::Vector<>::scalar_type >( val ) );
         }
      }
      crsMatrix_->fillComplete();
   }

   /// \brief Performs a matrix-vector multiplication in Trilinos.
   void apply( const TrilinosVector< FunctionType, MatrixScalarType >& src,
               TrilinosVector< FunctionType, MatrixScalarType >&       dst )
   {
      crsMatrix_->apply( *( src.getTpetraVector() ), *( dst.getTpetraVector() ) );
   }

   /// \brief Returns a string representation of this matrix.
   ///
   /// Must be called collectively by all processes.
   std::string to_string() const
   {
      std::stringstream ss;
      crsMatrix_->print( ss );
      return ss.str();
   }

   RCP< MatrixType > getTpetraMatrix() const { return crsMatrix_; }

 private:
   MPI_Comm                          trilinosCommunicatorRaw_;
   RCP< const Teuchos::Comm< int > > trilinosCommunicator_;
   RCP< const MapType >              rowMap_;
   RCP< MatrixType >                 crsMatrix_;
   uint_t                            level_;
};

} // namespace trilinos
} // namespace hyteg