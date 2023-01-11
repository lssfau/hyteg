/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl.
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
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/geometry/TokamakMap.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

void testTokamak( uint_t level, real_t errorL2Max )
{
   const uint_t maxLevel = level;

   const uint_t minLevel = maxLevel;
   const bool   writeVTK = false;

   // ITER configuration

   const uint_t                numToroidalSlices          = 8;
   const uint_t                numPoloidalSlices          = 6;
   const real_t                radiusOriginToCenterOfTube = 6.2;
   const std::vector< real_t > tubeLayerRadii             = { 3 };
   const real_t                torodialStartAngle         = 0.0;
   const real_t                polodialStartAngle         = 2.0 * pi / real_c( 2 * numPoloidalSlices );

   real_t delta = sin( 0.33 );
   real_t r1    = 2.0;
   real_t r2    = 3.7;

   const auto meshInfo = MeshInfo::meshTorus(
       numToroidalSlices, numPoloidalSlices, radiusOriginToCenterOfTube, tubeLayerRadii, torodialStartAngle, polodialStartAngle );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   TokamakMap::setMap( setupStorage,
                       numToroidalSlices,
                       numPoloidalSlices,
                       radiusOriginToCenterOfTube,
                       tubeLayerRadii,
                       torodialStartAngle,
                       polodialStartAngle,
                       delta,
                       r1,
                       r2 );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK( storage->hasGlobalCells() );

   writeDomainPartitioningVTK( storage, "../../output", "TokamakLaplaceTest" );

   typedef P1ElementwiseBlendingLaplaceOperator LaplaceOperator_T;

   LaplaceOperator_T laplaceOperator( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > exact = [&]( const Point3D& p ) -> real_t {
      return sin( p[0] / radiusOriginToCenterOfTube ) * sinh( p[1] / radiusOriginToCenterOfTube ) *
             ( p[2] / tubeLayerRadii.back() );
   };

   std::function< real_t( const Point3D& ) > zero = []( const Point3D& ) -> real_t { return 0.0; };

   std::function< real_t( const Point3D& ) > one = []( const Point3D& ) -> real_t { return 1.0; };

   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) -> real_t {
      return real_c( walberla::math::realRandom( 0.0, 1.0 ) );
   };

   P1Function< real_t > res( "r", storage, minLevel, maxLevel );
   P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   P1Function< real_t > uExact( "u_exact", storage, minLevel, maxLevel );
   P1Function< real_t > err( "err", storage, minLevel, maxLevel );

   u.interpolate( rand, maxLevel, DoFType::Inner );
   u.interpolate( exact, maxLevel, DoFType::DirichletBoundary );

   uExact.interpolate( exact, maxLevel, DoFType::All );

   auto solver = CGSolver< LaplaceOperator_T >( storage, minLevel, maxLevel );
   solver.setPrintInfo( false );

   VTKOutput vtkOutput( "../../output", "TokamakLaplaceTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( uExact );
   vtkOutput.add( err );

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   solver.solve( laplaceOperator, u, f, maxLevel );

   err.assign( { 1.0, -1.0 }, { u, uExact }, maxLevel );
   laplaceOperator.apply( u, res, minLevel, DoFType::Inner );

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, 1 );
   }

   auto unknowns = numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel );

   auto discrL2Residual = std::sqrt( res.dotGlobal( res, maxLevel, DoFType::Inner ) / real_c( unknowns ) );
   auto discrL2Error    = std::sqrt( err.dotGlobal( err, maxLevel, DoFType::Inner ) / real_c( unknowns ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residual L2 on level " << maxLevel << ": " << std::scientific << discrL2Residual
                                                      << " | Error L2: " << discrL2Error );

   WALBERLA_CHECK_LESS( discrL2Error, errorL2Max );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   bool longrun = false;
   for ( int i = 0; i < argc; i++ )
   {
      auto arg = std::string( argv[i] );
      if ( arg == "--longrun" )
      {
         longrun = true;
      }
   }

   testTokamak( 2, 8.2e-03 );
   if ( longrun )
   {
      testTokamak( 3, 2.9e-03 );
   }

   return 0;
}
