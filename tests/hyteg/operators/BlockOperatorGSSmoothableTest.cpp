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

#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"
#include "mixedOperator/VectorToVectorOperator.hpp"

// Test, that the Gauss-Seidel-Smoother of the block operator, finds in one iteration the correct solution.

using namespace hyteg;

// A simple P1xP1 test block function
template < typename value_t >
class TestBlockFunction : public BlockFunction< value_t >
{
 public:

   template < typename VType >
   using FunctionType = TestBlockFunction< VType >;

   TestBlockFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      size_t                                     minLevel,
                      size_t                                     maxLevel )
   : BlockFunction< value_t >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_p1", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_p1", storage, minLevel, maxLevel ) );
   };
};

// generate meta information
class TestBlockFunctionTag
{};

namespace hyteg {

template < typename VType >
struct FunctionTrait< TestBlockFunction< VType > >
{
   typedef VType                ValueType;
   typedef TestBlockFunctionTag Tag;

   static std::string getTypeName() { return "TestBlockFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

} // namespace hyteg

typedef TestBlockFunction< real_t > fType;

class DiagonalOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                         public GSSmoothable< P1Function< real_t > >
{
 public:
   using FType = P1Function< real_t >;

   DiagonalOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , diagonal_( "diagonal", storage, minLevel, maxLevel )
   , invDiagonal_( "invDiagonal", storage, minLevel, maxLevel )
   , tmp_( "tmp", storage, minLevel, maxLevel )
   {}

   using Operator< P1Function< real_t >, P1Function< real_t > >::getMinLevel;
   using Operator< P1Function< real_t >, P1Function< real_t > >::getMaxLevel;

   void setDiagonal( const std::function< real_t( Point3D ) >& f )
   {
      for ( uint_t level = getMinLevel(); level <= getMaxLevel(); level += 1 )
      {
         diagonal_.interpolate( f, level, All );
         invDiagonal_.assign( { 1 }, { diagonal_ }, level, All );
         invDiagonal_.invertElementwise( level, All, false );
      }
   }

   void apply( const FType& src,
               const FType& dst,
               size_t       level,
               DoFType      flag,
               UpdateType   updateType = Replace ) const override final
   {
      tmp_.multElementwise( { src, diagonal_ }, level, flag );
      if ( updateType == Replace )
         dst.assign( { 1 }, { tmp_ }, level, flag );
      else
         dst.assign( { 1, 1 }, { tmp_, dst }, level, flag );
   }

   void smooth_gs( const FType& dst, const FType& rhs, size_t level, DoFType flag ) const override
   {
      dst.multElementwise( { invDiagonal_, rhs }, level, flag );
   }

 protected:
   P1Function< real_t > diagonal_;
   P1Function< real_t > invDiagonal_;
   P1Function< real_t > tmp_;
};

std::shared_ptr< BlockOperator< fType, fType > >
    createOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
{
   // setup empty block operator for Taylor Hood functions
   auto oper = std::make_shared< BlockOperator< fType, fType > >( storage, minLevel, maxLevel, 2, 2 );

   auto op00 = std::make_shared< OperatorWrapper< DiagonalOperator > >( storage, minLevel, maxLevel );
   auto op01 = std::make_shared< OperatorWrapper< DiagonalOperator > >( storage, minLevel, maxLevel );
   auto op10 = std::make_shared< OperatorWrapper< DiagonalOperator > >( storage, minLevel, maxLevel );
   auto op11 = std::make_shared< OperatorWrapper< DiagonalOperator > >( storage, minLevel, maxLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level += 1 )
   {
      op00->unwrap().setDiagonal( []( auto p ) { return 2 * p[0] + p[1] + 1; } );
      op10->unwrap().setDiagonal( []( auto p ) { return p[0] - p[1]; } );
      op11->unwrap().setDiagonal( []( auto p ) { return 5 * p[0] - 0.5 * p[1] + 0.75; } );
      op01->unwrap().setDiagonal( []( auto ) { return 0; } );
   }

   oper->setSubOperator( 0, 0, op00 );
   oper->setSubOperator( 0, 1, op01 );
   oper->setSubOperator( 1, 0, op10 );
   oper->setSubOperator( 1, 1, op11 );

   return oper;
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing BlockGSSmootable " );
   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t minLevel = 3;
   uint_t maxLevel = 3;

   auto op = createOperator( storage, minLevel, maxLevel );

   fType x( "x", storage, minLevel, maxLevel );
   fType x_exact( "x_exact", storage, minLevel, maxLevel );
   fType b( "b", storage, minLevel, maxLevel );
   fType error( "error", storage, minLevel, maxLevel );

   auto dirichlet0 = []( auto p ) { return std::sin( 5 * p[0] ) * p[1]; };
   auto dirichlet1 = []( auto p ) { return std::cos( 5 * p[1] ) * p[0]; };

   for ( uint_t level = minLevel; level <= maxLevel; level += 1 )
   {
      // manufactured solution
      x_exact[0].interpolate( []( auto p ) { return p[0] + p[1]; }, level, All );
      x_exact[1].interpolate( []( auto p ) { return p[0] - p[1]; }, level, All );
      x_exact[0].interpolate( dirichlet0, level, DirichletBoundary );
      x_exact[1].interpolate( dirichlet1, level, DirichletBoundary );
      // random stuff
      x[0].interpolate( dirichlet0, level, All );
      x[1].interpolate( dirichlet1, level, All );
      // righthand side
      b[0].interpolate( []( auto p ) { return ( 2 * p[0] + p[1] + 1 ) * ( p[0] + p[1] ); }, level, All );
      b[1].interpolate(
          []( auto p ) { return ( p[0] - p[1] ) * ( p[0] + p[1] ) + ( 5. * p[0] - 0.5 * p[1] + 0.75 ) * ( p[0] - p[1] ); },
          level,
          All );
   }

   op->smooth_gs( x, b, maxLevel, Inner | NeumannBoundary );

   error.assign( { 1, -1 }, { x, x_exact }, maxLevel, All );

   const auto errorValue = error.dotGlobal( error, maxLevel, All );

   WALBERLA_CHECK_FLOAT_EQUAL( errorValue, 0. );

   return EXIT_SUCCESS;
}
