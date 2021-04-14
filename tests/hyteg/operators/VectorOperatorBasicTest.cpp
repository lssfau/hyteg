/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Marcus Mohr.
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
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// Perform a bacis compile and apply test for some VectorToVectorOperators

namespace hyteg {

template < typename vfType, typename opType >
static void runTest( bool beVerbose, std::string tag, std::string opName )
{

   WALBERLA_LOG_INFO_ON_ROOT( "checking " << opName );

   const uint_t minLevel = 3;
   const uint_t maxLevel = 3;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // setup two vector functions
   vfType srcVec( "srcVecFunc", storage, minLevel, maxLevel );
   vfType dstVec( "dstVecFunc", storage, minLevel, maxLevel );

   // setup our operator
   opType testOp( storage, minLevel, maxLevel );

   // check that we can apply it
   testOp.apply( srcVec, dstVec, maxLevel, All );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing VectorLaplaceOperators" );
   WALBERLA_LOG_INFO_ON_ROOT( "================================" );
   hyteg::runTest< hyteg::P1VectorFunction< walberla::real_t >, hyteg::P1ConstantVectorLaplaceOperator >(
       true, "P1", "P1ConstantVectorLaplaceOperator" );
   hyteg::runTest< hyteg::P1VectorFunction< walberla::real_t >, hyteg::P1ElementwiseVectorLaplaceOperator >(
       true, "P1", "P1ElementwiseVectorLaplaceOperator" );
   hyteg::runTest< hyteg::P1VectorFunction< walberla::real_t >, hyteg::P1ElementwiseBlendingVectorLaplaceOperator >(
       true, "P1", "P1ElementwiseBlendingVectorLaplaceOperator" );

   hyteg::runTest< hyteg::P2VectorFunction< walberla::real_t >, hyteg::P2ConstantVectorLaplaceOperator >(
       true, "P2", "P2ConstantVectorLaplaceOperator" );
   hyteg::runTest< hyteg::P2VectorFunction< walberla::real_t >, hyteg::P2ElementwiseVectorLaplaceOperator >(
       true, "P2", "P2ElementwiseVectorLaplaceOperator" );
   hyteg::runTest< hyteg::P2VectorFunction< walberla::real_t >, hyteg::P2ElementwiseBlendingVectorLaplaceOperator >(
       true, "P2", "P2ElementwiseBlendingVectorLaplaceOperator" );

   WALBERLA_LOG_INFO_ON_ROOT( "=============================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing VectorMassOperators" );
   WALBERLA_LOG_INFO_ON_ROOT( "=============================" );
   hyteg::runTest< hyteg::P1VectorFunction< walberla::real_t >, hyteg::P1ConstantVectorMassOperator >(
       true, "P1", "P1ConstantVectorMassOperator" );
   hyteg::runTest< hyteg::P1VectorFunction< walberla::real_t >, hyteg::P1ElementwiseVectorMassOperator >(
       true, "P1", "P1ElementwiseVectorMassOperator" );
   hyteg::runTest< hyteg::P1VectorFunction< walberla::real_t >, hyteg::P1ElementwiseBlendingVectorMassOperator >(
       true, "P1", "P1ElementwiseBlendingVectorMassOperator" );

   hyteg::runTest< hyteg::P2VectorFunction< walberla::real_t >, hyteg::P2ConstantVectorMassOperator >(
       true, "P2", "P2ConstantVectorMassOperator" );
   hyteg::runTest< hyteg::P2VectorFunction< walberla::real_t >, hyteg::P2ElementwiseVectorMassOperator >(
       true, "P2", "P2ElementwiseVectorMassOperator" );
   hyteg::runTest< hyteg::P2VectorFunction< walberla::real_t >, hyteg::P2ElementwiseBlendingVectorMassOperator >(
       true, "P2", "P2ElementwiseBlendingVectorMassOperator" );

   return EXIT_SUCCESS;
}
