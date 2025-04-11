/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P
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
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./SphericalShellBenchRotationMinimalMinimal.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle parameterConfig = cfg->getBlock( "Parameters" );
   const walberla::Config::BlockHandle setupConfig     = cfg->getBlock( "Setup" );

   WALBERLA_ROOT_SECTION()
   {
      setupConfig.listParameters();
      parameterConfig.listParameters();
   }

   const uint_t nTan = parameterConfig.getParameter< uint_t >( "nTan" );
   const uint_t nRad = parameterConfig.getParameter< uint_t >( "nRad" );

   const real_t rMin = parameterConfig.getParameter< real_t >( "rMin" );
   const real_t rMax = parameterConfig.getParameter< real_t >( "rMax" );

   auto meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( setupConfig.getParameter< bool >( "blending" ) )
   {
      hyteg::IcosahedralShellMap::setMap( *setupStorage );
   }

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   const uint_t minLevel = parameterConfig.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = parameterConfig.getParameter< uint_t >( "maxLevel" );

   const uint_t seed = parameterConfig.getParameter< uint_t >( "seed" );

   uint_t nFunctions = parameterConfig.getParameter< uint_t >( "nFunctions" );

   std::vector< std::shared_ptr< P2P1TaylorHoodFunction< real_t > > > testFunctions;

   for ( uint_t i = 0u; i < nFunctions; i++ )
   {
      testFunctions.emplace_back( std::make_shared< P2P1TaylorHoodFunction< real_t > >(
          std::string( "test" ) + std::to_string( i ), storage, minLevel, maxLevel ) );
   }

   walberla::math::seedRandomGenerator( seed );
   std::function< real_t( const Point3D& ) > randFunc = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   for ( uint_t i = 0u; i < nFunctions; i++ )
   {
      testFunctions[i]->interpolate( randFunc, maxLevel, All );
      real_t testDot  = testFunctions[i]->dotGlobal( *(testFunctions[i]), maxLevel, All );
      WALBERLA_LOG_INFO_ON_ROOT( "testDot" << i << " = " << testDot );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   return 0;
}
