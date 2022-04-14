/*
 * Copyright (c) 2020 Daniel Drzisga.
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
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_c;
using walberla::real_t;

using namespace hyteg;


template < typename StokesFunctionType, typename ProjectNormalOperatorType, bool use3D >
static void testProjectNormal( )
{
   const int level = use3D ? 2 : 3;

   auto meshInfo = MeshInfo::emptyMeshInfo();
   if ( use3D )
   {
      // this spheres tets are very distorted, but it is good (and fast) enough for this test:
      meshInfo = MeshInfo::meshSphericalShell( 2, 2, 0.5, 1.0 );
   }
   else
   {
      meshInfo = MeshInfo::meshAnnulus( 0.5, 1.0, MeshInfo::CRISS, 6, 6 );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 3, 0, true );

   if ( use3D )
   {
      IcosahedralShellMap::setMap( setupStorage );
   }
   else
   {
      AnnulusMap::setMap( setupStorage );
   }

   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto normalInterpolant = [] ( const Point3D & p ) {
     real_t norm = p.norm();
     real_t sign = (norm > 0.75) ? 1.0 : -1.0;
     return sign/norm * p;
   };

   auto normalFunction = [=]( const Point3D& p, Point3D& n ) -> void {
      n = normalInterpolant( p );
   };

   ProjectNormalOperatorType projectNormalOperator( storage, level, level, normalFunction );

   StokesFunctionType u( "u", storage, level, level );

   // we check if a radial function gets set to zero on the free slip boundary:
   u.uvw()[0].interpolate( [=](auto & p){ return normalInterpolant(p)[0]; }, level );
   u.uvw()[1].interpolate( [=](auto & p){ return normalInterpolant(p)[1]; }, level );
   if (storage->hasGlobalCells())
      u.uvw()[2].interpolate( [=](auto & p){ return normalInterpolant(p)[2]; }, level );
   WALBERLA_CHECK_GREATER( u.dotGlobal(u, level, FreeslipBoundary), 1 );
   projectNormalOperator.project( u, level, FreeslipBoundary );
   WALBERLA_CHECK_LESS( u.dotGlobal(u, level, FreeslipBoundary), 1e-14 );

   // we check if a tangential function is not changed by the projection operator:
   StokesFunctionType uTan( "uTan", storage, level, level );
   uTan.uvw()[0].interpolate( [=](auto & p){ return -p[1]; }, level );
   uTan.uvw()[1].interpolate( [=](auto & p){ return p[0]; }, level );
   if (storage->hasGlobalCells())
      uTan.uvw()[2].interpolate(0, level );
   u.assign({1}, {uTan}, level, All);
   WALBERLA_CHECK_GREATER( u.dotGlobal(u, level, FreeslipBoundary), 1 );
   projectNormalOperator.project( u, level, FreeslipBoundary );
   StokesFunctionType diff( "diff", storage, level, level );
   diff.assign( {1, -1}, {u, uTan}, level, All );
   WALBERLA_CHECK_LESS( diff.dotGlobal(diff, level, All), 1e-14 );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT("normal projection P1-P1 in 2D");
   testProjectNormal< P1StokesFunction< real_t >, P1ProjectNormalOperator, false >( );
   WALBERLA_LOG_INFO_ON_ROOT("normal projection P2-P1-TH in 2D");
   testProjectNormal< P2P1TaylorHoodFunction< real_t >, P2ProjectNormalOperator, false >( );
   WALBERLA_LOG_INFO_ON_ROOT("normal projection P1-P1 in 3D");
   testProjectNormal< P1StokesFunction< real_t >, P1ProjectNormalOperator, true >( );
   WALBERLA_LOG_INFO_ON_ROOT("normal projection P2-P1-TH in 3D");
   testProjectNormal< P2P1TaylorHoodFunction< real_t >, P2ProjectNormalOperator, true >( );
}
