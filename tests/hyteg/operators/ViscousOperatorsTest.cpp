/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

// Basic testing of the P2EpsilonOperator and P2FullViscousOperator
// TODO: extend

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/p2functionspace/P2EpsilonOperator.hpp"
#include "hyteg/p2functionspace/P2FullViscousOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

template < typename oper_t >
void checkObjectGeneration( std::string                                label,
                            std::shared_ptr< PrimitiveStorage >        primStore,
                            uint_t                                     minLevel,
                            uint_t                                     maxLevel,
                            std::function< real_t( const Point3D& ) >& callback )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Generating object of type '" << label << "'" );
   oper_t op( primStore, minLevel, maxLevel, callback );
}

template < typename oper_t >
void checkObjectGeneration( std::string                                label,
                            std::shared_ptr< PrimitiveStorage >        primStore,
                            uint_t                                     minLevel,
                            uint_t                                     maxLevel )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Generating object of type '" << label << "'" );
   oper_t op( primStore, minLevel, maxLevel );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::unique_ptr< SetupPrimitiveStorage > setStore;
   std::shared_ptr< PrimitiveStorage >      primStore;

   uint_t minLevel = 2;
   uint_t maxLevel = 3;

   Matrix2r mat;
   Point2D  vec;

   // isoviscous setting
   std::function< real_t( const Point3D& ) > viscosity = []( const Point3D& ) { return real_c( 1 ); };

   // ----------
   //  2D Tests
   // ----------
   logSectionHeader( "Testing 2D with BFS" );
   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkObjectGeneration< P2ConstantEpsilonOperator >( "P2ConstantEpsilonOperator", primStore, minLevel, maxLevel );
   checkObjectGeneration< P2ElementwiseAffineEpsilonOperator >(
       "P2ElementwiseAffineEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ElementwiseBlendingEpsilonOperator >(
       "P2ElementwiseBlendingEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ConstantFullViscousOperator >( "P2ConstantFullViscousOperator", primStore, minLevel, maxLevel );

   // ----------
   //  3D Tests
   // ----------
   logSectionHeader( "Testing 3D with pyramid_2el" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkObjectGeneration< P2ConstantEpsilonOperator >( "P2ConstantEpsilonOperator", primStore, minLevel, maxLevel );
   checkObjectGeneration< P2ElementwiseAffineEpsilonOperator >(
       "P2ElementwiseAffineEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ElementwiseBlendingEpsilonOperator >(
       "P2ElementwiseBlendingEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ConstantFullViscousOperator >( "P2ConstantFullViscousOperator", primStore, minLevel, maxLevel );

   return EXIT_SUCCESS;
}
