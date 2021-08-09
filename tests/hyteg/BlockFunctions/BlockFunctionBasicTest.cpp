/*
 * Copyright (c) 2021 Marcus Mohr.
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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodBlockFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// Perform some basic test to check that BlockFunctions can be instantiated,
// called, executed and exported

using namespace hyteg;

template < typename value_t >
class P1P1TaylorHoodStokesBlockFunction : public BlockFunction< value_t >
{
 public:
   P1P1TaylorHoodStokesBlockFunction( const std::string&                         name,
                                      const std::shared_ptr< PrimitiveStorage >& storage,
                                      size_t                                     minLevel,
                                      size_t                                     maxLevel )
   : BlockFunction< value_t >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1VectorFunction< value_t > > >( name + "_uvw", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_p", storage, minLevel, maxLevel ) );
   };
};

void logCall( const std::string& msg )
{
   size_t      len = msg.length();
   size_t      fw  = 30;
   std::string dots( len < fw ? fw - len : 3, '.' );
   WALBERLA_LOG_INFO_ON_ROOT( " -> " << msg << "() " << dots << " called" );
}

template < typename func_t >
void runTest( const std::string& kind )
{
   WALBERLA_LOG_INFO_ON_ROOT( "RUNNING WITH '" << kind << "'" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   BlockFunction< double > emptyFunc1( "Dummy BlockFunction<double>" );
   BlockFunction< float >  emptyFunc2( "Dummy BlockFunction<float>" );
   BlockFunction< uint_t > emptyFunc3( "Dummy BlockFunction<uint_t>" );

   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << emptyFunc1.getFunctionName() << "'" );
   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << emptyFunc2.getFunctionName() << "'" );
   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << emptyFunc3.getFunctionName() << "'" );

   uint_t minLevel = 1;
   uint_t maxLevel = 3;

   func_t stokes1( "Stokes #1", storage, minLevel, maxLevel );
   func_t stokes2( "Stokes #2", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << stokes1.getFunctionName() << "'" );
   WALBERLA_LOG_INFO_ON_ROOT( " -> dimension = '" << stokes1.getDimension() << "'" );
   WALBERLA_LOG_INFO_ON_ROOT( " -> # blocks  = '" << stokes1.getNumberOfBlocks() << "'" );

   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& x ) { return real_c( 2 ) * x[0]; };

   // -----------------------------------------
   //  Check whether we can call class methods
   // -----------------------------------------

#ifdef TIMING_TREE
   // Fails for P2P1TaylorHoodBlockFunction! Why?
   stokes1.enableTiming( storage->getTimingTree() );
   logCall( "enableTiming" );
#endif

   stokes1.interpolate( real_c( 5 ), maxLevel );
   logCall( "interpolate[constant]" );

   stokes1.interpolate( expr, maxLevel );
   logCall( "interpolate[expression]" );

   stokes1.interpolate( {expr}, maxLevel );
   logCall( "interpolate[vector]" );

   stokes1.add( real_c( 5 ), maxLevel );
   logCall( "add[scalar]" );

   stokes1.add( {real_c( 2 )}, {stokes2}, maxLevel );
   logCall( "add[functions]" );

   stokes1.swap( stokes2, maxLevel );
   logCall( "swap" );

   stokes1.assign( {2.0, -3.0}, {stokes1, stokes2}, maxLevel );
   logCall( "assign" );

   real_t aux = stokes1.dotGlobal( stokes2, maxLevel );
   WALBERLA_UNUSED( aux );
   logCall( "dotGlobal" );

   stokes1.multElementwise( {stokes1, stokes2}, maxLevel );
   logCall( "multElementwise" );

   uint_t nDoFs = stokes1.getNumberOfLocalDoFs( maxLevel );
   logCall( "getNumberOfLocalDoFs" );
   if constexpr ( std::is_same< P2P1TaylorHoodBlockFunction< real_t >, func_t >::value )
   {
      WALBERLA_CHECK_EQUAL( nDoFs, numberOfLocalDoFs< P2P1TaylorHoodBlockFunctionTag >( *storage, maxLevel ) );
   }
   else
   {
      WALBERLA_UNUSED( nDoFs );
   }
   // WALBERLA_LOG_INFO( " -> # of local DoFs = " << nDoFs );

   // -----------------------------------------------
   //  Check whether we can export the BlockFunction
   // -----------------------------------------------
   VTKOutput vtkOutput( "../../output", "BlockFunctionBasicTest", storage );
   vtkOutput.add( stokes1 );
   vtkOutput.write( maxLevel );
   logCall( "VTKOutput::add" );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing BlockFunction" );
   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );

   runTest< P1P1TaylorHoodStokesBlockFunction< real_t > >( "P1P1TaylorHoodStokesBlockFunction" );
   runTest< P2P1TaylorHoodBlockFunction< real_t > >( "P2P1TaylorHoodBlockFunction" );

   return EXIT_SUCCESS;
}
