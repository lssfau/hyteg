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
   TrilinosSparseMatrix( const OperatorType&                        op,
                         const std::shared_ptr< PrimitiveStorage >& storage,
                         const uint_t&                              level,
                         const FunctionType< PetscInt >&            numerator,
                         DoFType                                    flag = All )
   {
      trilinosCommunicatorRaw_ = walberla::mpi::MPIManager::instance()->comm();
      trilinosCommunicator_    = rcp( new Teuchos::MpiComm< int >( trilinosCommunicatorRaw_ ) );

      const uint_t numGlobalUnknowns =
          numberOfGlobalDoFs< typename OperatorType::dstType::Tag >( *storage, level, trilinosCommunicatorRaw_ );
      const size_t maxEntriesPerRow = 100; // only a hint, not restriction

      rowMap_    = rcp( new MapType( Tpetra::global_size_t( numGlobalUnknowns ), 0, trilinosCommunicator_ ) );
      crsMatrix_ = rcp( new MatrixType( rowMap_, maxEntriesPerRow ) );

      auto proxy = std::make_shared< TrilinosSparseMatrixProxy >( crsMatrix_ );
      hyteg::petsc::createMatrix( op, numerator, numerator, proxy, level, flag );
      crsMatrix_->fillComplete();
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

 private:
   typedef MatrixScalarType                            ST;
   typedef unsigned                                    LO;
   typedef unsigned                                    GO;
   typedef KokkosClassic::DefaultNode::DefaultNodeType NT;
   typedef Tpetra::Map<>                               MapType;
   typedef Tpetra::CrsMatrix<>                         MatrixType;

   MPI_Comm                          trilinosCommunicatorRaw_;
   RCP< const Teuchos::Comm< int > > trilinosCommunicator_;
   RCP< const MapType >              rowMap_;
   RCP< MatrixType >                 crsMatrix_;
};

} // namespace trilinos
} // namespace hyteg