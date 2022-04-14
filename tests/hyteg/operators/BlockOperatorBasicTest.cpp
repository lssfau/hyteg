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

#include "hyteg/composites/P2P1TaylorHoodBlockFunction.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// Perform a bacis compile and apply test for BlockOperator class

using namespace hyteg;

using thType = P2P1TaylorHoodBlockFunction< real_t >;

std::shared_ptr< OperatorWrapper< VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction > > >
    createP1EpsilonOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
{
   std::shared_ptr< OperatorWrapper< VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction > > > wrapped =
       std::make_shared< OperatorWrapper< VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction > > >(
           storage, minLevel, maxLevel );

   VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction > oper = wrapped->unwrap();

   oper.setSubOperator( 0, 0, std::make_shared< P1ConstantEpsilonOperator_11 >( storage, minLevel, maxLevel ) );
   oper.setSubOperator( 1, 0, std::make_shared< P1ConstantEpsilonOperator_12 >( storage, minLevel, maxLevel ) );
   oper.setSubOperator( 0, 1, std::make_shared< P1ConstantEpsilonOperator_21 >( storage, minLevel, maxLevel ) );
   oper.setSubOperator( 1, 1, std::make_shared< P1ConstantEpsilonOperator_22 >( storage, minLevel, maxLevel ) );

   return wrapped;
}

// NOTE: This is, of course, bad coding, since for kind = 2 our operator does not work on a function of thType.
//       It is funny that it works at all :) Need to change that by adding templates and a definition of a
//       P1P1StokesBlockFunction
std::shared_ptr< BlockOperator< thType, thType > > createOperator( uint_t                                     kind,
                                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                                   size_t                                     minLevel,
                                                                   size_t                                     maxLevel,
                                                                   uint_t                                     nRows,
                                                                   uint_t                                     nCols )
{
   // setup empty block operator for Taylor Hood functions
   auto oper = std::make_shared< BlockOperator< thType, thType > >( storage, minLevel, maxLevel, nRows, nCols );

   // fill operator with sub-blocks, depending on kind value
   switch ( kind )
   {
   case 1:
      oper->setSubOperator(
          0, 0, std::make_shared< OperatorWrapper< P2ConstantVectorLaplaceOperator > >( storage, minLevel, maxLevel ) );
      oper->setSubOperator(
          0, 1, std::make_shared< OperatorWrapper< P1ToP2ConstantDivTOperator > >( storage, minLevel, maxLevel ) );
      oper->setSubOperator(
          1, 0, std::make_shared< OperatorWrapper< P2ToP1ConstantDivOperator > >( storage, minLevel, maxLevel ) );
      oper->setSubOperator( 1, 1, nullptr );
      break;

   case 2:
      oper->setSubOperator( 0, 0, createP1EpsilonOperator( storage, minLevel, maxLevel ) );
      oper->setSubOperator(
          0, 1, std::make_shared< OperatorWrapper< P1ToP2ConstantDivTOperator > >( storage, minLevel, maxLevel ) );
      oper->setSubOperator(
          1, 0, std::make_shared< OperatorWrapper< P2ToP1ConstantDivOperator > >( storage, minLevel, maxLevel ) );
      oper->setSubOperator( 1, 1, std::make_shared< OperatorWrapper< P1PSPGOperator > >( storage, minLevel, maxLevel ) );
      break;

   default:
      WALBERLA_LOG_INFO_ON_ROOT( "Unsupported kind value '" << kind << "'" );
   }
   return oper;
}

void runTest( uint_t kind, std::string tag, const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
{
   WALBERLA_LOG_INFO_ON_ROOT( "RUNNING with " << tag );

   // need functions to work on
   thType src( "src", storage, minLevel, maxLevel );
   thType dst( "dst", storage, minLevel, maxLevel );

   // let factory generate an operator
   auto stokesOper = createOperator( kind, storage, minLevel, maxLevel, 2, 2 );
   WALBERLA_LOG_INFO_ON_ROOT( "Factory generation ............. worked" );

   // apply block operator to block function
   stokesOper->apply( src, dst, minLevel, All );
   WALBERLA_LOG_INFO_ON_ROOT( "BlockOperator application ...... worked" );

   // extract and apply "A" block separately to velocity
   const std::shared_ptr< GenericOperator< real_t > > aBlock = stokesOper->getSubOperator( 0, 0 );
   aBlock->apply( src.getSubFunction( 0 ), dst.getSubFunction( 0 ), maxLevel, Inner, Replace );
   WALBERLA_LOG_INFO_ON_ROOT( "'A-block' application .......... worked" );

   // apply stabilisation operator, if present
   auto stabOper = stokesOper->getSubOperator( 1, 1 );
   if ( stabOper != nullptr )
   {
      stabOper->apply( src.getSubFunction( 1 ), dst.getSubFunction( 1 ), maxLevel, All );
      WALBERLA_LOG_INFO_ON_ROOT( "PSPG block application ......... worked" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "=======================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing BlockOperator" );
   WALBERLA_LOG_INFO_ON_ROOT( "=======================" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t minLevel = 1;
   uint_t maxLevel = 1;

   BlockOperator< thType, thType > blockOper( storage, minLevel, maxLevel, 2, 2 );

   thType src( "src", storage, minLevel, maxLevel );
   thType dst( "dst", storage, minLevel, maxLevel );

   // checking no-op
   blockOper.apply( src, dst, maxLevel, All );

   // use factory
   runTest( 1, "P2P1ConstantStokesOperator", storage, minLevel, maxLevel );
   runTest( 2, "P1P1ConstantEpsilonOperator", storage, minLevel, maxLevel );

   return EXIT_SUCCESS;
}
