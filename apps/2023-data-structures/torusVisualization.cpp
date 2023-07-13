/*
* Copyright (c) 2023 Daniel Bauer.
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

#include "core/DataTypes.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/TokamakMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t level      = 3;
   const uint_t resolution = uint_c( std::pow( level, 2 ) ) + 1;

   const uint_t                toroidalResolution         = 34;
   const uint_t                poloidalResolution         = 6;
   const real_t                radiusOriginToCenterOfTube = 2;
   const std::vector< real_t > tubeLayerRadii             = { 0.4 };
   const real_t                torodialStartAngle         = 0.0;
   const real_t                polodialStartAngle         = 0.0;
   const real_t                delta                      = 0;
   const real_t                r1                         = tubeLayerRadii.back();
   const real_t                r2                         = tubeLayerRadii.back();

   const MeshInfo        torus = MeshInfo::meshTorus( toroidalResolution,
                                               poloidalResolution,
                                               radiusOriginToCenterOfTube,
                                               tubeLayerRadii,
                                               torodialStartAngle,
                                               polodialStartAngle );
   SetupPrimitiveStorage setupStorage( torus, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   TokamakMap::setMap( setupStorage,
                       toroidalResolution,
                       poloidalResolution,
                       radiusOriginToCenterOfTube,
                       tubeLayerRadii,
                       torodialStartAngle,
                       polodialStartAngle,
                       delta,
                       r1,
                       r2 );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // edges
   writeBlendedCoarseMeshVTK( storage, "output", "torusVisualization", resolution );

   // faces
   n1e1::N1E1VectorFunction< real_t > u( "u", storage, level, level );
   VTKOutput                          vtk( "output", "torusVisualization", storage );
   vtk.add( u );
   vtk.write( level );

   return EXIT_SUCCESS;
}
