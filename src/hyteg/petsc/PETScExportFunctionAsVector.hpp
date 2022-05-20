/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <fstream>

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

/// \brief Exports finite-element function as vector to an ASCII file
///
/// Uses PETSc functionality to generate the vector associated with the given function
/// and export if to a file in ASCII forma which can be executed as a Matlab function.
///
/// \param function         the HyTeG FE function to export
/// \param fileName         name of file to write the vector to
/// \param vectorName       name of the vector, must be different from base part of fileName;
///                         when the script is executed in Matlab this will be the final
///                         name of the vector
/// \param storage          primitive storage
/// \param level            refinement level on which the operator matrix should be generated
/// \param beVerbose        should function be talkative or not
///
template < template < class > class FunctionType, class FunctionTag >
void exportFunction( const FunctionType< real_t >&       function,
                     std::string                         fileName,
                     std::string                         vectorName,
                     std::shared_ptr< PrimitiveStorage > storage,
                     uint_t                              level,
                     bool                                beVerbose = false,
                     bool                                binary    = false,
                     PetscViewerFormat                   format    = PETSC_VIEWER_ASCII_MATRIXMARKET )
{
   // Get dimension of function space
   uint_t localDoFs  = hyteg::numberOfLocalDoFs< FunctionTag >( *storage, level );
   uint_t globalDoFs = hyteg::numberOfGlobalDoFs< FunctionTag >( *storage, level );
   if ( localDoFs != globalDoFs )
   {
      WALBERLA_ABORT( "localDoFs and globalDoFs must agree for this app!" );
   }
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Dimension of function space is " << globalDoFs );
   }
   FunctionType< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   // Fire up PETSc
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Firing up PETSc" );
   }
   PETScManager pmgr;

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Converting function to PETSc vector" );
   }
   PETScVector< real_t, FunctionType > petscVector( function, numerator, level, All, vectorName );

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Exporting vector to file '" << fileName << "'" );
   }
   petscVector.print( fileName.c_str(), binary, format );
}

} // namespace hyteg

#endif
