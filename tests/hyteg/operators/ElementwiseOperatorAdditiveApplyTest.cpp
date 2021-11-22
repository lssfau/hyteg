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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

// This test checks if the additive application of the elementwise
// operators works correctly.
//
// "Additive application" here means that
//     u := Ax + By
// can be performed without any auxiliary function.

using walberla::real_t;
using namespace hyteg;

template < typename OpTypeConst, typename OpTypeElem, typename FuncType >
void additiveApplyTest( std::shared_ptr< PrimitiveStorage > storage, const uint_t level )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Additive apply test" )

   const real_t epsilon = 1e-12;

   // functions
   FuncType srcA( "srcA", storage, level, level );
   FuncType srcB( "srcB", storage, level, level );

   FuncType dstManualConst( "dstManualConst", storage, level, level );
   FuncType dstManualElementwise( "dstManualElementwise", storage, level, level );
   FuncType dstAdditiveConst( "dstAdditiveConst", storage, level, level );
   FuncType dstAdditiveElementwise( "dstAdditiveElementwise", storage, level, level );
   FuncType dstConstReplaceElementwiseAdd( "dstConstReplaceElementwiseAdd", storage, level, level );

   FuncType tmpA( "tmpA", storage, level, level );
   FuncType tmpB( "tmpB", storage, level, level );

   FuncType error( "error", storage, level, level );

   // setup operators
   OpTypeElem  elemWiseOp( storage, level, level );
   OpTypeConst constantOp( storage, level, level );

   auto functionA = []( const Point3D& x ) { return std::sin( x[0] ) + 6.0 * std::sin( x[1] * x[1] * x[1] ) + x[2] * x[2] * x[2]; };

   auto functionB = []( const Point3D& x ) { return 2 * std::sin( 2.0 * x[0] ) + 7.0 * std::cos( x[1] ) + std::cosh( x[2] * x[2] ); };

   srcA.interpolate( functionA, level );
   srcB.interpolate( functionB, level );

   // manual add const
   constantOp.apply( srcA, tmpA, level, All, Replace );
   constantOp.apply( srcB, tmpB, level, All, Replace );
   dstManualConst.assign( {1.0, 1.0}, {tmpA, tmpB}, level, All );

   // manual add elementwise
   elemWiseOp.apply( srcA, tmpA, level, All, Replace );
   elemWiseOp.apply( srcB, tmpB, level, All, Replace );
   dstManualElementwise.assign( {1.0, 1.0}, {tmpA, tmpB}, level, All );

   // const replace, const add
   constantOp.apply( srcA, dstAdditiveConst, level, All, Replace );
   constantOp.apply( srcB, dstAdditiveConst, level, All, Add );

   // elementwise replace, elementwise add
   elemWiseOp.apply( srcA, dstAdditiveElementwise, level, All, Replace );
   elemWiseOp.apply( srcB, dstAdditiveElementwise, level, All, Add );

   // const replace, elementwise add
   //
   // It's not obvious, but this is the most important test case.
   // Due to the ghost-layer handling, if the elementwise operator is applied twice
   // (1. replace, 2. add) the result is added up in the ghost-layers of the highest-dimensional
   // primitives. Therefore, it looks like the implementation is correct although the additive
   // communication still overwrites the data on the lower-dimensional primitives.
   //
   // In this test, the application of the constant (stencil-based) operator followed by
   // the additive elementwise operator triggers that exact issue.
   //
   constantOp.apply( srcA, dstConstReplaceElementwiseAdd, level, All, Replace );
   elemWiseOp.apply( srcB, dstConstReplaceElementwiseAdd, level, All, Add );

   error.assign( {1.0, -1.0}, {dstManualConst, dstManualElementwise}, level, All );
   auto errorMax = error.getMaxMagnitude( level );
   WALBERLA_CHECK_LESS( errorMax, epsilon, "manual elementwise" );

   error.assign( {1.0, -1.0}, {dstManualConst, dstAdditiveConst}, level, All );
   errorMax = error.getMaxMagnitude( level );
   WALBERLA_CHECK_LESS( errorMax, epsilon, "const replace, const add" );

   error.assign( {1.0, -1.0}, {dstManualConst, dstAdditiveElementwise}, level, All );
   errorMax = error.getMaxMagnitude( level );
   WALBERLA_CHECK_LESS( errorMax, epsilon, "elementwise replace, elementwise add" );

   error.assign( {1.0, -1.0}, {dstManualConst, dstConstReplaceElementwiseAdd}, level, All );
   errorMax = error.getMaxMagnitude( level );
   WALBERLA_CHECK_LESS( errorMax, epsilon, "const replace, elementwise add" );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   std::string           meshFileName = "../../data/meshes/quad_16el.msh";
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t level = 4;

   // -------------------
   //  Run some 2D tests
   // -------------------
   WALBERLA_LOG_INFO_ON_ROOT( "P1, Laplace, 2D" )
   additiveApplyTest< P1ConstantLaplaceOperator, P1ElementwiseLaplaceOperator, P1Function< real_t > >( storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "P1, Mass, 2D" )
   additiveApplyTest< P1ConstantMassOperator, P1ElementwiseMassOperator, P1Function< real_t > >( storage, level );

   level = 3;
   WALBERLA_LOG_INFO_ON_ROOT( "P2, Laplace, 2D" )
   additiveApplyTest< P2ConstantLaplaceOperator, P2ElementwiseLaplaceOperator, P2Function< real_t > >( storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "P2, Mass, 2D" )
   additiveApplyTest< P2ConstantMassOperator, P2ElementwiseMassOperator, P2Function< real_t > >( storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "P2, DivKGrad, 2D" )
   additiveApplyTest< P2ConstantLaplaceOperator, P2ElementwiseDivKGradOperator, P2Function< real_t > >( storage, level );

   // ----------------------------
   //  Prepare setup for 3D tests
   // ----------------------------
   meshFileName                     = "../../data/meshes/3D/pyramid_tilted_4el.msh";
   MeshInfo              meshInfo3D = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage3D( meshInfo3D, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage3D );
   std::shared_ptr< PrimitiveStorage > storage3D = std::make_shared< PrimitiveStorage >( setupStorage3D );

   // -------------------
   //  Run some 3D tests
   // -------------------
   level = 3;
   WALBERLA_LOG_INFO_ON_ROOT( "P1, Laplace, 3D" )
   additiveApplyTest< P1ConstantLaplaceOperator, P1ElementwiseLaplaceOperator, P1Function< real_t > >( storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "P1, Mass, 3D" )
   additiveApplyTest< P1ConstantMassOperator, P1ElementwiseMassOperator, P1Function< real_t > >( storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "P2, Laplace, 3D" )
   additiveApplyTest< P2ConstantLaplaceOperator, P2ElementwiseLaplaceOperator, P2Function< real_t > >( storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "P2, Mass, 3D" )
   additiveApplyTest< P2ConstantMassOperator, P2ElementwiseMassOperator, P2Function< real_t > >( storage, level );

   return 0;
}
