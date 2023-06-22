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
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionWrapper.hpp"
#include "hyteg/operators/GenericOperator.hpp"
#include "hyteg/operators/OperatorWrapper.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// Perform some basic construction and apply test for OperatorWrapper class

using namespace hyteg;

template < typename oper_t >
static void runTest( std::string opName, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with --> " << opName << " <--" );

   typedef typename oper_t::srcType src_t;
   typedef typename oper_t::dstType dst_t;

   WALBERLA_LOG_INFO_ON_ROOT( "* Constructing OperatorWrapper<" << opName << ">" );
   OperatorWrapper< oper_t > wrapped( storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "* Constructing GenericOperator pointer" );
   std::shared_ptr< GenericOperator< real_t > > generic =
       std::make_shared< OperatorWrapper< oper_t > >( storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "* Unwrapping GenericOperator" );
   const oper_t& unwrapped = generic->unwrap< oper_t >();

   WALBERLA_LOG_INFO_ON_ROOT( "* Testing apply() #1" );
   src_t srcFunc( "source", storage, minLevel, maxLevel );
   dst_t dstFunc( "destination", storage, minLevel, maxLevel );
   unwrapped.apply( srcFunc, dstFunc, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT( "* Testing apply() #2" );
   std::shared_ptr< GenericFunction< real_t > > srcGeneric =
       std::make_shared< FunctionWrapper< src_t > >( "gSrc", storage, minLevel, maxLevel );
   std::shared_ptr< GenericFunction< real_t > > dstGeneric =
       std::make_shared< FunctionWrapper< dst_t > >( "gDst", storage, minLevel, maxLevel );
   wrapped.apply( *srcGeneric, *dstGeneric, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT( "done\n" );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing OperatorWrapper" );
   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t minLevel = 1;
   uint_t maxLevel = 1;

   runTest< P2ConstantLaplaceOperator >( "P2ConstantLaplaceOperator", storage, minLevel, maxLevel );
   runTest< P1ElementwiseBlendingVectorMassOperator >( "P1ElementwiseBlendingVectorMassOperator", storage, minLevel, maxLevel );

   return EXIT_SUCCESS;
}
