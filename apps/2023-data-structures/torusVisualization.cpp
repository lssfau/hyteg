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
#include "hyteg/geometry/TorusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

std::shared_ptr< PrimitiveStorage >
    createStorage( const uint_t toroidalResolution, const uint_t poloidalResolution, const std::vector< real_t > tubeLayerRadii )
{
   const real_t radiusOriginToCenterOfTube = 2;
   const real_t torodialStartAngle         = 0.0;
   const real_t polodialStartAngle         = 0.0;

   const MeshInfo torus = MeshInfo::meshTorus( toroidalResolution,
                                               poloidalResolution,
                                               radiusOriginToCenterOfTube,
                                               tubeLayerRadii,
                                               torodialStartAngle,
                                               polodialStartAngle );

   SetupPrimitiveStorage setupStorage( torus, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   TorusMap::setMap( setupStorage,
                     toroidalResolution,
                     poloidalResolution,
                     radiusOriginToCenterOfTube,
                     tubeLayerRadii,
                     torodialStartAngle,
                     polodialStartAngle );
   return std::make_shared< PrimitiveStorage >( setupStorage );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // edges
   {
      std::shared_ptr< PrimitiveStorage > storage = createStorage( 136, 10, { 0.1, 0.2, 0.3, 0.4 } );

      const uint_t resolution = 3;
      writeBlendedCoarseMeshVTK( storage, "output", "torusVisualization", resolution );
   }

   // faces
   {
      std::shared_ptr< PrimitiveStorage > storage = createStorage( 34, 6, { 0.4 } );
      const uint_t                        level   = 5;

      const real_t R = 2;
      const real_t r = real_c( 0.4 );

      const auto analyticalSol = [R, r]( const Point3D& xVec ) {
         const real_t x    = xVec[0];
         const real_t y    = xVec[1];
         const real_t z    = xVec[2];
         const real_t tmp0 = real_c( std::sqrt( std::pow( x, 2 ) + std::pow( y, 2 ) ) );
         const real_t tmp1 = real_c( -std::pow( r, 2 ) + std::pow( z, 2 ) + std::pow( -R + tmp0, 2 ) ) / tmp0;
         const real_t u0   = tmp1 * y;
         const real_t u1   = -tmp1 * x;
         const real_t u2   = 0;
         return Point3D{ u0, u1, u2 };
      };
      n1e1::N1E1VectorFunction< real_t > sol( "sol", storage, level, level );
      sol.interpolate( analyticalSol, level );

      VTKOutput vtk( "output", "torusVisualization", storage );
      vtk.add( sol );
      vtk.write( level );
   }

   return EXIT_SUCCESS;
}
