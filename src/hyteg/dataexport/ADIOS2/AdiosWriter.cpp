/*
 * Copyright (c) 2023 Marcus Mohr, Roman Freissler.
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

#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"

#include <adios2.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP1.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

void AdiosWriter::write( const uint_t level, const uint_t timestep )
{
   communication::syncRegisteredFunctions( feFunctionRegistry_, level, communication::syncDirection_t::LOW2HIGH );

   // for each registered function type check whether a writer for the given level
   // already exists

   // -------------------
   //  P1 Type Functions
   // -------------------
   uint_t numP1TypeFunctions = feFunctionRegistry_.getP1Functions().size() + feFunctionRegistry_.getP1VectorFunctions().size();
   // WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter: " << numP1TypeFunctions << " functions of P1 type are registered!" );
   if ( numP1TypeFunctions > 0 )
   {
      if ( p1Writers_.count( level ) == 0 )
      {
         p1Writers_[level] =
             std::make_unique< AdiosWriterForP1 >( adios_, filePath_, fileBaseName_, engineType_, level, storage_ );
      }

      p1Writers_[level]->write( feFunctionRegistry_, timestep, userProvidedParameters_, userDefinedAttributes_ );
   }

   // -------------------
   //  P2 Type Functions
   // -------------------
   uint_t numP2TypeFunctions = feFunctionRegistry_.getP2Functions().size() + feFunctionRegistry_.getP2VectorFunctions().size();
   // WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter: " << numP2TypeFunctions << " functions of P2 type are registered!" );
   if ( numP2TypeFunctions > 0 )
   {
      if ( p2Writers_.count( level ) == 0 )
      {
         p2Writers_[level] =
             std::make_unique< AdiosWriterForP2 >( adios_, filePath_, fileBaseName_, engineType_, level, storage_ );
      }

      p2Writers_[level]->write( feFunctionRegistry_, timestep, userProvidedParameters_, userDefinedAttributes_ );
   }

   // remember that we had our first write() episode
   firstWriteDidHappen_ = true;
}

} // namespace hyteg
