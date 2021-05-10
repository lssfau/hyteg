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

// #include "hyteg/communication/Syncing.hpp"
// #include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
// #include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
// #include "hyteg/dataexport/VTKOutput.hpp"
// #include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
// #include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// Perform some basic test on the FunctionWrapper class

#include "hyteg/functions/FunctionWrapper.hpp"

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing FunctionWrapper" );
   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t minLevel = 1;
   uint_t maxLevel = 1;

   FunctionWrapper< P1Function< real_t > >       p1Wrap( "P1func", storage, minLevel, maxLevel );
   FunctionWrapper< P2Function< real_t > >       p2Wrap( "P2func", storage, minLevel, maxLevel );
   FunctionWrapper< P1VectorFunction< real_t > > p1vecWrap( "P1VecFunc", storage, minLevel, maxLevel );
   FunctionWrapper< P2VectorFunction< real_t > > p2vecWrap( "P2VecFunc", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "P1func -> dimension = " << p1Wrap.getDimension() );
   WALBERLA_LOG_INFO_ON_ROOT( "P2func -> dimension = " << p2Wrap.getDimension() );
   WALBERLA_LOG_INFO_ON_ROOT( "P1VecFunc -> dimension = " << p1vecWrap.getDimension() );
   WALBERLA_LOG_INFO_ON_ROOT( "P2VecFunc -> dimension = " << p2vecWrap.getDimension() );

   // access storage
   WALBERLA_CHECK_EQUAL( p1Wrap.getStorage()->hasGlobalCells(), false );
   WALBERLA_CHECK_EQUAL( p2Wrap.getStorage()->hasGlobalCells(), false );
   WALBERLA_CHECK_EQUAL( p1vecWrap.getStorage()->hasGlobalCells(), false );
   WALBERLA_CHECK_EQUAL( p2vecWrap.getStorage()->hasGlobalCells(), false );

   // test unwrapping
   GenericFunction< real_t >* ptr    = &p1Wrap;
   P1Function< real_t >&      p1Func = ptr->template unwrap< P1Function< real_t > >();
   WALBERLA_CHECK_EQUAL( p1Func.getFunctionName(), p1Wrap.getFunctionName() );
   WALBERLA_LOG_INFO_ON_ROOT( "Successfully unwrapped '" << p1Func.getFunctionName() << "'" );

   // detection of incorrect unwrapping
   // (only works w/o further debugging, as then assertions kill the process)
#ifdef NDEBUG
   uint_t failed = 0;
   try
   {
      real_t aux = p2Wrap.dotGlobal( p1Wrap, maxLevel, Inner );
      WALBERLA_UNUSED( aux );
   } catch ( const std::exception& ex )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Exception encountered is '" << ex.what() << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( "Unwrapping failed, as it should ;-)" );
      failed = 1;
   }
   WALBERLA_CHECK( failed );
#endif

   // check assign
   FunctionWrapper< P1Function< real_t > > p1WrapOther( "Another P1func", storage, minLevel, maxLevel );
   p1Wrap.assign( {1.0, -2.0}, {p1Wrap, p1WrapOther}, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "assign() -> check" );

   // check inner product
   p2Wrap.interpolate( real_t( 2 ), maxLevel, All );
   p2Wrap.interpolate( real_t( 0 ), maxLevel, Boundary );
   real_t aux   = p2Wrap.dotGlobal( p2Wrap, maxLevel, Inner );
   uint_t nDoFs = numberOfGlobalInnerDoFs< FunctionTrait< P2Function< real_t > >::Tag >( *storage, maxLevel );
   WALBERLA_CHECK_FLOAT_EQUAL( aux, real_c( nDoFs * 4 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "dotGlobal() -> check" );

   // difficult check
   p2vecWrap.interpolate( real_t( 2 ), maxLevel, All );
   P2VectorFunction< real_t > p2vec = p2vecWrap.unwrap();
   WALBERLA_CHECK_FLOAT_EQUAL( p2vec[0].getMaxMagnitude( maxLevel ), real_c( 2 ) );

   p2vecWrap.multElementwise( {p2vecWrap, p2vecWrap}, maxLevel, All );
   WALBERLA_CHECK_FLOAT_EQUAL( p2vec[0].getMaxMagnitude( maxLevel ), real_c( 4 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "P2VecFunc.interpolate() -> check" );

   return EXIT_SUCCESS;
}
