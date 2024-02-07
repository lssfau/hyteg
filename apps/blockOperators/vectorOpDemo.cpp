/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "mixed_operator/VectorLaplaceOperator.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using terraneo::SphericalHarmonicsTool;
using walberla::real_c;
using walberla::real_t;
using walberla::math::pi;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ============
   //  Parameters
   // ============

   // check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      // auto defaultFile = "./SPHdemo.prm";
      // WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      // cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   // const walberla::Config::BlockHandle params = cfg->getBlock( "Parameters" );

   // =========
   //  Meshing
   // =========
   const uint_t level = 5; // params.getParameter< uint_t >( "level" );

   WALBERLA_LOG_INFO_ON_ROOT( "Generating criss mesh on unit square" );
   Point2D  cornerLL( 0.0, 0.0 );
   Point2D  cornerUR( 1.0, 1.0 );
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   meshInfo          = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CRISS, 1, 1 );

   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   // ============
   //  Testing P2
   // ============
   P2VectorFunction< real_t > input( "input", storage, level, level );
   P2VectorFunction< real_t > output( "output", storage, level, level );

   uint_t freq1 = 1;
   uint_t freq2 = 3;

   std::function< real_t( const Point3D& ) > inFunc1 = [freq1]( const Point3D& x ) {
      real_t m = real_c( freq1 );
      return std::sin( m * pi * x[0] ) * std::sin( m * pi * x[1] );
   };

   std::function< real_t( const Point3D& ) > inFunc2 = [freq2]( const Point3D& x ) {
      real_t m = real_c( freq2 );
      return std::sin( m * pi * x[0] ) * std::sin( m * pi * x[1] );
   };

   input.interpolate( {inFunc1, inFunc2}, level );

   P2ConstantVectorLaplaceOperator vecLap( storage, level, level );
   vecLap.apply( input, output, level, Inner );

   // ============
   //  Testing P1
   // ============
   P1VectorFunction< real_t > inputP1( "input", storage, level, level );
   P1VectorFunction< real_t > outputP1( "output", storage, level, level );

   inputP1.interpolate( {inFunc1, inFunc2}, level );

   P1ElementwiseVectorLaplaceOperator vecLapP1( storage, level, level );
   vecLapP1.apply( inputP1, outputP1, level, Inner );

   // output data for visualisation
   bool outputVTK = true;
   if ( outputVTK )
   {
      VTKOutput vtkOutput( "./output", "vectorOpDemo", storage );
      vtkOutput.add( input );
      vtkOutput.add( output );
      vtkOutput.add( inputP1 );
      vtkOutput.add( outputP1 );
      vtkOutput.write( level );
   }

   return EXIT_SUCCESS;
}
