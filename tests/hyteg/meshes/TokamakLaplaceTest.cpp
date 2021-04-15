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

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t minLevel = 3;
   const uint_t maxLevel = 3;
   const bool   writeVTK = true;

   const auto            meshInfo = MeshInfo::meshTokamak( 24, 3, 0.7, 1.6, 0.7, true, true );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK( storage->hasGlobalCells() );

   writeDomainPartitioningVTK( storage, "../../output", "TokamakLaplaceTest" );

   P1ConstantLaplaceOperator laplaceOperator( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > exact = []( const Point3D& p ) -> real_t {
      return sin( p[0] ) * sinh( p[1] ) * p[2];
   };

   std::function< real_t( const Point3D& ) > zero = []( const Point3D& ) -> real_t { return 0.0; };

   std::function< real_t( const Point3D& ) > one = []( const Point3D& ) -> real_t { return 1.0; };

   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) -> real_t {
      return walberla::math::realRandom( 0.0, 1.0 );
   };

   P1Function< real_t > res( "r", storage, minLevel, maxLevel );
   P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   P1Function< real_t > uExact( "u_exact", storage, minLevel, maxLevel );
   P1Function< real_t > err( "err", storage, minLevel, maxLevel );

   u.interpolate( rand, maxLevel, DoFType::Inner );
   u.interpolate( exact, maxLevel, DoFType::DirichletBoundary );

   uExact.interpolate( exact, maxLevel, DoFType::All );

   auto solver = CGSolver< P1ConstantLaplaceOperator >( storage, minLevel, maxLevel );
   solver.setPrintInfo( true );

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

   WALBERLA_LOG_INFO( "Residual L2 on level " << maxLevel << ": " << std::scientific << discrL2Residual
                                              << " | Error L2: " << discrL2Error );

   return 0;
}
