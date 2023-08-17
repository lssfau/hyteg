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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"

namespace hyteg {
real_t runTest( uint_t coarseRefinements, uint_t level, hyteg::MeshInfo meshInfo )
{
   if ( coarseRefinements > 0 )
   {
      meshInfo = MeshInfo::refinedCoarseMesh( meshInfo, coarseRefinements );
   }
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1ConstantLaplaceOperator laplaceOperator3D( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& p ) -> real_t {
      return sin( p[0] ) * sinh( p[1] ) * p[2];
   };

   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) -> real_t { return 0.0; };

   std::function< real_t( const hyteg::Point3D& ) > one = []( const hyteg::Point3D& ) -> real_t { return 1.0; };

   std::function< real_t( const hyteg::Point3D& ) > rand = []( const hyteg::Point3D& ) -> real_t {
      return real_c( walberla::math::realRandom( 0.0, 1.0 ) );
   };

   hyteg::P1Function< real_t > res( "r", storage, level, level );
   hyteg::P1Function< real_t > f( "f", storage, level, level );
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > uExact( "u_exact", storage, level, level );
   hyteg::P1Function< real_t > err( "err", storage, level, level );
   hyteg::P1Function< real_t > oneFunction( "oneFunction", storage, level, level );

   auto numUnknowns = numberOfGlobalDoFs< P1FunctionTag >( *storage, level );

   u.interpolate( rand, level, DoFType::Inner );
   u.interpolate( exact, level, DoFType::DirichletBoundary );
   f.interpolate( zero, level, DoFType::All );
   res.interpolate( zero, level, DoFType::All );
   uExact.interpolate( exact, level, DoFType::All );
   oneFunction.interpolate( one, level, DoFType::All );

   auto solver = CGSolver< P1ConstantLaplaceOperator >( storage, level, level );
   //solver.setPrintInfo(true);
   solver.solve( laplaceOperator3D, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, uExact }, level );
   laplaceOperator3D.apply( u, res, level, DoFType::Inner );

   auto discrL2Error = std::sqrt( err.dotGlobal( err, level, DoFType::Inner ) / real_c( numUnknowns ) );

   return discrL2Error;
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   WALBERLA_LOG_DEVEL( "Start" )
   hyteg::MeshInfo  meshInfo      = hyteg::MeshInfo::meshCuboid( { -0.03, 0.05, -0.07 }, { 1.01, 0.99, 1.03 }, 1, 1, 1 );
   walberla::real_t normalLevel3  = hyteg::runTest( 0, 3, meshInfo );
   walberla::real_t coarsenLevel1 = hyteg::runTest( 1, 1, meshInfo );
   walberla::real_t coarsenLevel2 = hyteg::runTest( 1, 2, meshInfo );
   WALBERLA_CHECK_FLOAT_EQUAL( normalLevel3, coarsenLevel2 )
   WALBERLA_CHECK_LESS( coarsenLevel2 / coarsenLevel1, 0.3 )

   WALBERLA_LOG_DEVEL( "discrL2Error            level 3 : " << normalLevel3 )

   WALBERLA_LOG_DEVEL( "discrL2Error Refined 1, level 2 : " << coarsenLevel1 )
   WALBERLA_LOG_DEVEL( "discrL2Error Refined 1, level 3 : " << coarsenLevel2 << " rate: " << coarsenLevel2 / coarsenLevel1 )
}