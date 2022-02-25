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

#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/trilinos/KokkosWrapper.hpp"
#include "hyteg/trilinos/TeuchosWrapper.hpp"
#include "hyteg/trilinos/TpetraWrapper.hpp"
#include "hyteg/trilinos/TrilinosVectorProxy.hpp"

namespace hyteg {
namespace trilinos {

using walberla::real_t;
using walberla::uint_t;

using Teuchos::RCP;
using Teuchos::rcp;

template < template < class > class FunctionType, typename FunctionScalarType = real_t >
class TrilinosVector
{
 public:
   typedef Tpetra::Map<>                                     MapType;
   typedef Tpetra::Map<>::local_ordinal_type                 LO;
   typedef Tpetra::Map<>::global_ordinal_type                GO;
   typedef Tpetra::MultiVector< FunctionScalarType, LO, GO > VectorType;

   /// \brief Allocates the Trilinos vector data structure.
   ///
   /// Meaningful values have to be written to the vector through the assemble() call.
   TrilinosVector( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t& level )
   : level_( level )
   {
      trilinosCommunicatorRaw_ = walberla::mpi::MPIManager::instance()->comm();
      trilinosCommunicator_    = rcp( new Teuchos::MpiComm< int >( trilinosCommunicatorRaw_ ) );

      const uint_t numGlobalUnknowns =
          numberOfGlobalDoFs< typename FunctionType< FunctionScalarType >::Tag >( *storage, level, trilinosCommunicatorRaw_ );

      rowMap_ = rcp( new MapType( Tpetra::global_size_t( numGlobalUnknowns ), 0, trilinosCommunicator_ ) );
      vector_ = rcp( new VectorType( rowMap_, 1 ) );
   }

   void fillFromFunction( const FunctionType< FunctionScalarType >& function,
                          const FunctionType< idx_t >&              numerator,
                          DoFType                                   flag = All )
   {
      auto proxy = std::make_shared< TrilinosVectorProxy< VectorType > >( vector_ );
      function.toVector( numerator, proxy, level_, flag );
   }

   RCP< VectorType > getTpetraVector() const { return vector_; }

   void writeToFunction( const FunctionType< FunctionScalarType >& function,
                         const FunctionType< idx_t >&              numerator,
                         DoFType                                   flag = All )
   {
      auto proxy = std::make_shared< TrilinosVectorProxy< VectorType > >( vector_ );
      function.fromVector( numerator, proxy, level_, flag );
   }

   /// \brief Returns a string representation of this vector.
   ///
   /// Must be called collectively by all processes.
   std::string to_string() const
   {
      std::stringstream ss;
      vector_->print( ss );
      return ss.str();
   }

 private:
   MPI_Comm                          trilinosCommunicatorRaw_;
   RCP< const Teuchos::Comm< int > > trilinosCommunicator_;
   RCP< const MapType >              rowMap_;
   RCP< VectorType >                 vector_;
   uint_t                            level_;
};

} // namespace trilinos
} // namespace hyteg
