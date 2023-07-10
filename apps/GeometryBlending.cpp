/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <core/Environment.h>
#include <core/config/Config.h>

#include "core/timing/Timer.h"

#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   /// create enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const size_t level = 2;

   /// read mesh file and create storage
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../data/meshes/unitsquare_with_circular_hole.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{0.5, 0.5, 0}};
   real_t  circleRadius = 0.25;

   for( auto it : setupStorage.getFaces() )
   {
      Face& face = *(it.second);

      std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
      neighborEdgesOnBoundary.erase(
         std::remove_if( neighborEdgesOnBoundary.begin(), neighborEdgesOnBoundary.end(),
                      [ &setupStorage ]( const PrimitiveID & id ){ return !setupStorage.onBoundary( id ); } )
         , neighborEdgesOnBoundary.end());

      if( neighborEdgesOnBoundary.size() > 0 )
      {
         Edge& edge = *setupStorage.getEdge( neighborEdgesOnBoundary[0] );

         if( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.4 )
         {
            setupStorage.setGeometryMap( edge.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
            setupStorage.setGeometryMap( face.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
         }
      }
   }

   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto x = std::make_shared< hyteg::P2Function< real_t > >( "x", storage, level, level );
   auto y = std::make_shared< hyteg::P2Function< real_t > >( "y", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > tmp_x = [&]( const hyteg::Point3D& x_ ) { return x_[0]; };

   std::function< real_t( const hyteg::Point3D& ) > tmp_y = [&]( const hyteg::Point3D& x_ ) { return x_[1]; };

   x->interpolate( tmp_x, level, hyteg::All );
   y->interpolate( tmp_y, level, hyteg::All );

   VTKOutput vtkOutput("../output", "GeometryBlending", storage);
   vtkOutput.add( *x );
   vtkOutput.add( *y );
   vtkOutput.write( level );

   return 0;
}
