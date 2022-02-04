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

#include "hyteg/functions/FunctionWrapper.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/facedofspace/FaceDoFFunction.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// Perform some basic test on the FunctionWrapper class

using namespace hyteg;

template < typename func_t >
void wrapFunction( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, uint_t dim )

{
   std::string               kind = FunctionTrait< func_t >::getTypeName();
   FunctionWrapper< func_t > wrappedFunc( kind, storage, minLevel, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "-> wrapped a '" << kind << "'" );
   WALBERLA_CHECK_EQUAL( dim, wrappedFunc.getDimension() );
}

#ifdef HYTEG_BUILD_WITH_PETSC
template < template < typename > class FunctionKind >
void testPETScConversion( const std::shared_ptr< PrimitiveStorage >& storage )
{
   WALBERLA_LOG_INFO_ON_ROOT( "-> Running test for " << FunctionTrait< FunctionKind< real_t > >::getTypeName() );

   uint_t level = 4;

   FunctionWrapper< FunctionKind< real_t > > src( "testing", storage, level, level );
   FunctionWrapper< FunctionKind< real_t > > dst( "testing", storage, level, level );
   FunctionWrapper< FunctionKind< idx_t > >  numerator( "numerator", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > expression = []( const hyteg::Point3D& x ) {
      real_t value;
      value = std::sin( real_c( 3 ) * x[0] ) + real_c( 0.5 ) * x[1] * x[1];
      return value;
   };

   numerator.enumerate( level );
   src.interpolate( expression, level );

   PETScVector< real_t, GenericFunction > vector( src, numerator, level, All );
   vector.createFunctionFromVector( dst, numerator, level, All );
   dst.assign( { real_c( 1 ), real_c( -1 ) }, { dst, src }, level, All );
   // real_t diff = dst.getMaxComponentMagnitude( level, All );
   // WALBERLA_CHECK_FLOAT_EQUAL( diff, real_c( 0 ) );
   // WALBERLA_LOG_INFO_ON_ROOT( "   max. difference = " << diff );
}
#endif

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

   // WALBERLA_LOG_INFO_ON_ROOT( "P1func -> dimension = " << p1Wrap.getDimension() );
   // WALBERLA_LOG_INFO_ON_ROOT( "P2func -> dimension = " << p2Wrap.getDimension() );
   // WALBERLA_LOG_INFO_ON_ROOT( "P1VecFunc -> dimension = " << p1vecWrap.getDimension() );
   // WALBERLA_LOG_INFO_ON_ROOT( "P2VecFunc -> dimension = " << p2vecWrap.getDimension() );

   // check wrapping for different function classes
   WALBERLA_LOG_INFO_ON_ROOT( "Testing wrapping functions:" );
   wrapFunction< vertexdof::VertexDoFFunction< real_t > >( storage, minLevel, maxLevel, 1 );
   wrapFunction< EdgeDoFFunction< real_t > >( storage, minLevel, maxLevel, 1 );
   wrapFunction< FaceDoFFunction< real_t > >( storage, minLevel, maxLevel, 1 );
   wrapFunction< P1Function< real_t > >( storage, minLevel, maxLevel, 1 );
   wrapFunction< P2Function< real_t > >( storage, minLevel, maxLevel, 1 );
   wrapFunction< P1VectorFunction< real_t > >( storage, minLevel, maxLevel, 2 );
   wrapFunction< P2VectorFunction< real_t > >( storage, minLevel, maxLevel, 2 );

   // some wraps for eating later ;-)
   FunctionWrapper< P1Function< real_t > >       p1Wrap( "P1func", storage, minLevel, maxLevel );
   FunctionWrapper< P2Function< real_t > >       p2Wrap( "P2func", storage, minLevel, maxLevel );
   FunctionWrapper< P1VectorFunction< real_t > > p1vecWrap( "P1VecFunc", storage, minLevel, maxLevel );
   FunctionWrapper< P2VectorFunction< real_t > > p2vecWrap( "P2VecFunc", storage, minLevel, maxLevel );

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

   // test unwrapping with free-function
   const P2Function< real_t >& foo = unwrap( p2Wrap );
   WALBERLA_LOG_INFO_ON_ROOT( "Successfully unwrapped '" << foo.getFunctionName() << "'" );

   // detection of incorrect unwrapping
   // (only works w/o further debugging, as then assertions kill the process)
#ifdef NDEBUG
   uint_t failed = 0;
   if ( walberla::mpi::MPIManager::instance()->numProcesses() == 1 )
   {
      try
      {
         WALBERLA_LOG_INFO( "Trying to unwrap p1Wrap to P2Function (bad idea)" );
         real_t aux = p2Wrap.dotGlobal( p1Wrap, maxLevel, Inner );
         WALBERLA_UNUSED( aux );
      } catch ( const std::exception& ex )
      {
         WALBERLA_LOG_INFO( "Unwrapping failed, as it should ;-)" );
         WALBERLA_LOG_INFO( "Exception encountered is '" << ex.what() << "'" );
         failed = 1;
      } catch ( ... )
      {
         WALBERLA_LOG_INFO( "Unwrapping failed, as it should ;-)" );
         WALBERLA_LOG_INFO( "Precise type of exception unknown" );
         failed = 1;
      }
      if ( !failed )
      {
         WALBERLA_LOG_INFO( "Unwrapping worked, but it should not :-(" );
      }
      WALBERLA_CHECK( failed );
   }
#endif

   // check getting info from a GenericFunction on its type
   GenericFunction< real_t >*   ptr1 = &p1Wrap;
   GenericFunction< real_t >*   ptr2 = &p2Wrap;
   functionTraits::FunctionKind fk1  = ptr1->getFunctionKind();
   functionTraits::FunctionKind fk2  = ptr2->getFunctionKind();
   WALBERLA_CHECK_EQUAL( fk1, functionTraits::P1_FUNCTION );
   WALBERLA_CHECK_EQUAL( fk2, functionTraits::P2_FUNCTION );
   WALBERLA_LOG_INFO_ON_ROOT( "getFunctionKind() -> check" );

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

#ifdef HYTEG_BUILD_WITH_PETSC
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " Function <-> Vector Conversion" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------" );
   PETScManager petscManager( &argc, &argv );
   testPETScConversion< P1Function >( storage );
   testPETScConversion< P2Function >( storage );
   testPETScConversion< EdgeDoFFunction >( storage );
   testPETScConversion< P1VectorFunction >( storage );
   testPETScConversion< P2VectorFunction >( storage );
#endif

   return EXIT_SUCCESS;
}
