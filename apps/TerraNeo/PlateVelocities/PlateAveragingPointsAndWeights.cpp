/*
* Copyright (c) 2025 Nils Kohl
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

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKPointCloudOutput.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "terraneo/plates/LocalAveragingPointWeightProvider.hpp"
#include "terraneo/plates/PlateVelocityProvider.hpp"
#include "vtk/VTKOutput.h"

using walberla::int_c;
using walberla::real_c;
using walberla::real_t;
using namespace hyteg;
using namespace terraneo;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   // plates::UniformCirclesPointWeightProvider provider( { { 1.0 / 200.0, 6 }, { 1.0 / 100.0, 12 } }, 1e-2 );
   plates::FibonacciLatticePointWeightProvider provider( 100, 1 / 100., 1e-2 );

   const uint_t level  = 4;
   const uint_t nTan   = 3;
   const real_t radius = real_c( 1 );

   MeshInfo              meshInfo = hyteg::MeshInfo::meshThinSphericalShell( nTan, radius );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   ThinShellMap::setMap( setupStorage, radius );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   P2VectorFunction< real_t > u( "u", storage, level, level );

   VTKOutput vtkFunction( "./output", "vtk_function", storage );
   vtkFunction.add( u );

   VTKPointCloudOutput vtkPoints( "./output", "vtk_points" );

   std::vector< vec3D > centersLonLat;
   centersLonLat.push_back( vec3D( 36.0, 36.0, 1.0 ) );

   for ( auto centerLonLat : centersLonLat )
   {
      const auto samples = provider.samplePointsAndWeightsLonLat( centerLonLat );
      for ( const auto& [p, w] : samples )
      {
         const auto pcart = conversions::sph2cart( { p.x(), p.y() }, p.z() );
         vtkPoints.addPoints( { pcart } );
         vtkPoints.addValues( "weights", { w } );
      }
   }

   vtkFunction.write( level );
   vtkPoints.write();
}