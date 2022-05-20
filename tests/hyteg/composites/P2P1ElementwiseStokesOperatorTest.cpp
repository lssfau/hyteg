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

#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hyteg;

void compareApplyCC( const MeshInfo & meshInfo, const uint_t level, bool precomputeElementMatrices )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Test compare apply CC" )

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const real_t epsilon = 1e-13;

   // functions
   P2P1TaylorHoodFunction< real_t > src( "src", storage, level, level );

   P2P1TaylorHoodFunction< real_t > dstStencilCC( "dstStencilCC", storage, level, level );
   P2P1TaylorHoodFunction< real_t > dstElementwiseCC( "dstElementwiseCC", storage, level, level );

   P2P1TaylorHoodFunction< real_t > error( "error", storage, level, level );

   // setup operators
   P2P1TaylorHoodStokesOperator stencilOp( storage, level, level );
   P2P1ElementwiseConstantCoefficientStokesOperator elemWiseOp( storage, level, level );

   if ( precomputeElementMatrices )
   {
      elemWiseOp.computeAndStoreLocalElementMatrices();
   }

   // interpolate something on src function
   auto vel_x = []( const Point3D & p ){
      return std::sin( p[0] ) + std::sin( 1.3 * p[1] ) + std::sin( 1.7 * p[2] );
   };

   auto vel_y = []( const Point3D & p ){
     return std::sin( p[0] ) + std::sin( 2.3 * p[1] ) + std::sin( 2.7 * p[2] );
   };

   auto vel_z = []( const Point3D & p ){
     return std::sin( p[0] ) + std::sin( 3.3 * p[1] ) + std::sin( 3.7 * p[2] );
   };

   auto pressure = []( const Point3D & p ){
      return p[0] * p[0] * p[0] + std::cos( p[1] * p[2] );
   };

   src.uvw().interpolate( { vel_x, vel_y, vel_z }, level, All );
   src.p().interpolate( pressure, level, All );

   stencilOp.apply( src, dstStencilCC, level, Inner | NeumannBoundary );
   elemWiseOp.apply( src, dstElementwiseCC, level, Inner | NeumannBoundary );

   error.assign( {1.0, -1.0}, {dstStencilCC, dstElementwiseCC}, level, All );

   auto errorMaxMagnitudeU = error.uvw()[0].getMaxMagnitude( level );
   auto errorMaxMagnitudeV = error.uvw()[1].getMaxMagnitude( level );
   auto errorMaxMagnitudeW = error.uvw().getDimension() == 3 ? error.uvw()[2].getMaxMagnitude( level ) : real_c( 0 );
   auto errorMaxMagnitudeP = error.p().getMaxMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "Error max magnitude u: " << errorMaxMagnitudeU );
   WALBERLA_LOG_INFO_ON_ROOT( "Error max magnitude v: " << errorMaxMagnitudeV );
   WALBERLA_LOG_INFO_ON_ROOT( "Error max magnitude w: " << errorMaxMagnitudeW );
   WALBERLA_LOG_INFO_ON_ROOT( "Error max magnitude p: " << errorMaxMagnitudeP );

   WALBERLA_CHECK_LESS( errorMaxMagnitudeU, epsilon );
   WALBERLA_CHECK_LESS( errorMaxMagnitudeV, epsilon );
   WALBERLA_CHECK_LESS( errorMaxMagnitudeW, epsilon );
   WALBERLA_CHECK_LESS( errorMaxMagnitudeP, epsilon );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" ), 3, false );
   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" ), 4, false );

   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_tilted_4el.msh" ), 3, false );
   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_tilted_4el.msh" ), 4, false );

   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" ), 3, true );
   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" ), 4, true );

   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_tilted_4el.msh" ), 3, true );
   compareApplyCC( MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_tilted_4el.msh" ), 4, true );

   return 0;
}
