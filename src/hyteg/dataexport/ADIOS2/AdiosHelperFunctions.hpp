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
#pragma once

#include <adios2.h>

#include "hyteg/functions/FunctionMultiStore.hpp"

namespace hyteg::adiosHelpers {

std::string generateVTKMetaInfo( const std::vector< std::string >& namesOfPointDataFunctions,
                                 const std::vector< std::string >& namesOfCellDataFunctions );

/// Name of scalar variable to be used in exporting time-step information
///
/// We use this variable to ensure consistency between the different places where this
/// variable name is needed, e.g in genewrateVTKMetaInfo and putTimeStepInfo. Note, though,
/// that currently TIME is the only allowed value.
extern const std::string nameOfTimeStepVariable;

/// Schedule information on current time step to be exported
inline void putTimeStepInfo( adios2::IO& io, adios2::Engine& engine, uint_t timestep )
{
   adios2::Variable< real_t > varTimeStep = io.InquireVariable< real_t >( nameOfTimeStepVariable );
   if ( !varTimeStep )
   {
      varTimeStep = io.DefineVariable< real_t >( nameOfTimeStepVariable );
   }
   engine.Put( varTimeStep, real_c( timestep ) );
}

/// Check whether the current process owns primitives with data to export (faces, or cell depending on dimension)
inline bool mpiProcessHasMacrosOfHighestDimension( const std::shared_ptr< PrimitiveStorage >& storage )
{
   bool weStoreRelevantData = false;
   if ( storage->hasGlobalCells() && storage->getCells().size() > 0 )
   {
      return true;
   }
   else if ( storage->getFaces().size() > 0 )
   {
      return true;
   }
   return false;
}

} // namespace hyteg::adiosHelpers
