/*
 * Copyright (c) 2024 Michael Zikeli.
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

#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/LaTeX/Table.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {
using namespace walberla;

template < bool is3D, bool isP2, class Function_t >
void innerFkt( latex::Table< 9 >& benchmarkTable )
{
   std::string domain;
   uint_t      minLevel = 0;
   uint_t      maxLevel;
   if constexpr ( is3D )
   {
      domain = "../../../data/meshes/3D/cube_24el.msh";
      if constexpr ( isP2 )
      {
         maxLevel = 8;
      }
      else
      {
         maxLevel = 10;
      }
   }
   else
   {
      domain = "../../../data/meshes/quad_8el.msh";
      if constexpr ( isP2 )
      {
         maxLevel = 12;
      }
      else
      {
         maxLevel = 15;
      }
   }

   SetupPrimitiveStorage setupStorage( MeshInfo::fromGmshFile( domain ), uint_c( mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   Function_t function( "P1Fct", storage, minLevel, maxLevel );

   for ( auto lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      uint_t colIdx = 0;
      if constexpr ( is3D )
      {
         colIdx += 2;
      }
      else
      {
         colIdx += 0;
      }
      if constexpr ( isP2 )
      {
         colIdx += 6;
      }
      else
      {
         colIdx += 2;
      }

      if constexpr ( !isP2 && !is3D )
      {
         benchmarkTable.addElement( lvl, 0, lvl );
         benchmarkTable.addElement( lvl, 1, '&' );
         benchmarkTable.addElement( lvl, 3, '&' );
         benchmarkTable.addElement( lvl, 5, '&' );
         benchmarkTable.addElement( lvl, 7, '&' );
      }

      const auto dofs = real_c( function.getNumberOfGlobalDoFs( lvl ) );
      benchmarkTable.addElement( lvl, colIdx, dofs );
   }
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   using namespace walberla;
   using namespace hyteg;

   // +++ Setup Environment +++
   Environment env( argc, argv );
   logging::Logging::instance()->setLogLevel( logging::Logging::INFO );
   MPIManager::instance()->useWorldComm();

   latex::Table< 9 > benchmarkTable( { "Level",
                                       "Sep1",
                                       "DoFs-2D-P1",
                                       "Sep2"
                                       "DoFs-3D-P1",
                                       "Sep3"
                                       "DoFs-2D-P2",
                                       "Sep4"
                                       "DoFs-3D-P2" } );

   innerFkt< false, false, P1Function< real_t > >( benchmarkTable );
   innerFkt< true, false, P1Function< real_t > >( benchmarkTable );
   innerFkt< false, true, P2Function< real_t > >( benchmarkTable );
   innerFkt< true, true, P2Function< real_t > >( benchmarkTable );

   benchmarkTable.write( "./", "memoryConsumption" );

   return EXIT_SUCCESS;
}