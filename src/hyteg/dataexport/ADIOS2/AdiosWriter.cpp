/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include <adios2.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP1.hpp"
#include "hyteg/dataexport/FEFunctionRegistry.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

void AdiosWriter::write( const uint_t level, const uint_t timestep )
{
   communication::syncRegisteredFunctions( feFunctionRegistry_, level );

   // for each registered function type check whether a writer for the given level
   // already exists
   uint_t numP1TypeFunctions = feFunctionRegistry_.getP1Functions().size();
   WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter: " << numP1TypeFunctions << " functions of P1 type are registered!" );
   if ( numP1TypeFunctions > 0 )
   {
      if ( p1Writers_.count( level ) == 0 )
      {
        p1Writers_[level] = std::make_unique< AdiosWriterForP1 >( adios_, filePath_, fileBaseName_, level, storage_ );
      }

      p1Writers_[level]->write( feFunctionRegistry_, timestep );
   }

   WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter::write() not fully functional, yet!" );
}

} // namespace hyteg
