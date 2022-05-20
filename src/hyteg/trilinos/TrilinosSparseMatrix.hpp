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
#include "core/Format.hpp"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/sparseassembly/DirichletBCs.hpp"
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
   typedef Tpetra::Map<>                                 MapType;
   typedef Tpetra::Map<>::local_ordinal_type             LO;
   typedef Tpetra::Map<>::global_ordinal_type            GO;
   typedef Tpetra::CrsMatrix< MatrixScalarType, LO, GO > MatrixType;

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

      numGlobalUnknowns_ = numberOfGlobalDoFs< typename OperatorType::dstType::Tag >( *storage, level, trilinosCommunicatorRaw_ );
      const size_t maxEntriesPerRow = 100; // only a hint, not restriction

      rowMap_    = rcp( new MapType( Tpetra::global_size_t( numGlobalUnknowns_ ), 0, trilinosCommunicator_ ) );
      crsMatrix_ = rcp(
          new MatrixType( rowMap_, maxEntriesPerRow ) );
   }

   /// \brief Assembles the sparse matrix from the passed operator.
   ///
   /// This routine assembles the "Neumann" system. To incorporate Dirichlet boundary
   /// conditions, subsequent calls to applyDirichletBoundaryConditions*() must follow.
   void assembleSystem( const OperatorType& op, const FunctionType< idx_t >& numerator )
   {
      if ( crsMatrix_->isFillComplete() )
      {
         crsMatrix_->resumeFill();
      }

      auto proxy = std::make_shared< TrilinosSparseMatrixProxy >( crsMatrix_ );
      // We insert zeros on the diagonal to be able to later modify the diagonal values in case
      // there is a zero in the original system. It seems like we cannot insert values later on
      // after calling fillComplete() (?).
      for ( uint_t idx = 0; idx < numGlobalUnknowns_; idx++ )
      {
         proxy->addValue( idx, idx, 0 );
      }
      op.toMatrix( proxy, numerator, numerator, level_, All );
      crsMatrix_->fillComplete();
   }

   /// \brief Modifies the system matrix to conform to the underlying Dirichlet problem.
   ///
   /// Sets the diagonal entry of all rows that correspond to DoFs on a Dirichlet boundary to 1.
   void applyDirichletBoundaryConditions( const FunctionType< idx_t >& numerator )
   {
      std::vector< idx_t > dirichletRowIndices;
      hyteg::applyDirichletBC( numerator, dirichletRowIndices, level_ );

      zeroRows( dirichletRowIndices, 1.0 );
   }

   /// \brief Zeroes all entries of the specified matrix row but the diagonal value, which is set to the passed value.
   ///
   /// Can be used for asymmetric pinning if the nullspace of the operator is not zero.
   ///
   /// \param row the global row index of the row that shall be zeroed
   /// \param diagonalValue the value to be put on the diagonal (set to 0.0 to zero entire row)
   void zeroRow( const idx_t& row, const MatrixScalarType& diagonalValue = 1.0 ) { zeroRows( { row }, diagonalValue ); }

   /// \brief Zeroes all entries of the specified matrix rows but the diagonal values, which are set to the passed value.
   ///
   /// Can be used for asymmetric pinning if the nullspace of the operator is not zero.
   ///
   /// \param rows the global row indices of the rows that shall be zeroed
   /// \param diagonalValue the value to be put on the diagonal (set to 0.0 to zero entire row)
   void zeroRows( const std::vector< idx_t >& rows, const MatrixScalarType& diagonalValue = 1.0 )
   {
      WALBERLA_CHECK( crsMatrix_->hasColMap(), "Trilinos matrix was not assembled correctly." );
      RCP< const MapType > rowMap    = crsMatrix_->getRowMap();
      RCP< const MapType > columnMap = crsMatrix_->getColMap();

      if ( crsMatrix_->isFillComplete() )
      {
         crsMatrix_->resumeFill();
      }

      for ( auto row : rows )
      {
         if ( !rowMap->isNodeGlobalElement( row ) )
            continue;

         auto                                         localRow = rowMap->getLocalElement( GO( row ) );
         Teuchos::ArrayView< const LO >               columnIndicesView;
         Teuchos::ArrayView< const MatrixScalarType > columnValuesView;
         crsMatrix_->getLocalRowView( localRow, columnIndicesView, columnValuesView );
         for ( const auto& localCol : columnIndicesView )
         {
            auto col = columnMap->getGlobalElement( localCol );
            if ( row != col )
               crsMatrix_->replaceGlobalValues( row, Teuchos::tuple< GO >( col ), Teuchos::tuple< MatrixScalarType >( 0 ) );
         }
         crsMatrix_->replaceGlobalValues( row, Teuchos::tuple< GO >( row ), Teuchos::tuple< MatrixScalarType >( diagonalValue ) );
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
      ss << Teuchos::describe( *crsMatrix_, Teuchos::EVerbosityLevel::VERB_EXTREME );
      return ss.str();
   }

   /// \brief Exports the matrix into a Matlab file. When resulting script is executed, the matrix appears
   ///        in the Matlab workspace.
   void exportToMatlabFormat( const std::string& filePath, const std::string& matlabVariableName = "Mat" )
   {
      WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(),
                            1,
                            "Trilinos sparse matrix Matlab output currently only supported for serial runs." )

      std::ofstream out( filePath );
      out << matlabVariableName << " = zeros(" << crsMatrix_->getGlobalNumEntries() << ",3);" << std::endl;
      out << matlabVariableName << " = [" << std::endl;

      for ( size_t r = 0; r < crsMatrix_->getNodeNumRows(); ++r )
      {
         const size_t nE  = crsMatrix_->getNumEntriesInLocalRow( r );
         GO           gid = crsMatrix_->getRowMap()->getGlobalElement( r );
         if ( crsMatrix_->isGloballyIndexed() )
         {
            Teuchos::ArrayView< const GO >               rowinds;
            Teuchos::ArrayView< const MatrixScalarType > rowvals;
            crsMatrix_->getGlobalRowView( gid, rowinds, rowvals );
            for ( size_t j = 0; j < nE; ++j )
            {
               out << gid + 1 << ", " << rowinds[j] + 1 << ", " << walberla::format( "%23.16e", rowvals[j] ) << std::endl;
            }
         }
         else if ( crsMatrix_->isLocallyIndexed() )
         {
            Teuchos::ArrayView< const LO >               rowinds;
            Teuchos::ArrayView< const MatrixScalarType > rowvals;
            crsMatrix_->getLocalRowView( r, rowinds, rowvals );
            for ( size_t j = 0; j < nE; ++j )
            {
               out << gid + 1 << ", " << crsMatrix_->getColMap()->getGlobalElement( rowinds[j] ) + 1 << ", "
                   << walberla::format( "%23.16e", rowvals[j] ) << std::endl;
            }
         }
      }

      out << "];" << std::endl;
      out << matlabVariableName << " = spconvert(" << matlabVariableName << ");";
   }

   RCP< MatrixType > getTpetraMatrix() const { return crsMatrix_; }

 private:
   MPI_Comm                          trilinosCommunicatorRaw_;
   RCP< const Teuchos::Comm< int > > trilinosCommunicator_;
   RCP< const MapType >              rowMap_;
   RCP< MatrixType >                 crsMatrix_;
   uint_t                            level_;
   uint_t                            numGlobalUnknowns_;
};

} // namespace trilinos
} // namespace hyteg
