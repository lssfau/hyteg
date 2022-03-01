/*
 * Copyright (c) 2020 Daniel Drzisga, Andreas Wagner, Nils Kohl.
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

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1BlendingStokesOperator.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

std::shared_ptr< SetupPrimitiveStorage >
    setupStorageRectangle( const double channelLength, const double channelHeight, const uint_t ny )
{
   Point2D left( { -channelLength / 2, 0 } );
   Point2D right( { channelLength / 2, channelHeight } );

   const uint_t    nx           = ny * static_cast< uint_t >( channelLength / channelHeight );
   hyteg::MeshInfo meshInfo     = hyteg::MeshInfo::meshRectangle( left, right, MeshInfo::CROSS, nx, ny );
   auto            setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( *setupStorage );

   // Boundaries
   auto inflow   = [=]( auto p ) { return p[0] <= -channelLength / 2 + 1e-14; };
   auto outflow  = [=]( auto p ) { return p[0] >= +channelLength / 2 - 1e-14; };
   auto noslip   = [=]( auto p ) { return p[1] >= +channelHeight - 1e-14; };
   auto freeslip = [=]( auto p ) { return p[1] <= 1e-14; };

   setupStorage->setMeshBoundaryFlagsByVertexLocation( 2, outflow );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 3, freeslip );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 1, noslip );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 1, inflow );

   return setupStorage;
}

template < typename StokesFunctionType, typename StokesOperatorType, typename ProjectNormalOperatorType >
void run( const real_t absErrorTolerance )
{
   // solver parameters
   const uint_t minLevel = 2;
   const uint_t maxLevel = 2;

   // geometry rectangle
   real_t       channelLength = 0.5;
   real_t       channelHeight = 0.5;
   const uint_t ny            = 1;

   auto setupStorage = setupStorageRectangle( channelLength, channelHeight, ny );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage );

   StokesFunctionType                                          u_src( "u_src", storage, minLevel, maxLevel );
   StokesFunctionType                                          u_dst_hyteg( "u_dst_hyteg", storage, minLevel, maxLevel );
   StokesFunctionType                                          u_dst_petsc( "u_dst_petsc", storage, minLevel, maxLevel );
   StokesFunctionType                                          diff( "diff", storage, minLevel, maxLevel );
   typename StokesFunctionType::template FunctionType< idx_t > numerator( "numerator", storage, minLevel, maxLevel );

   numerator.enumerate( maxLevel );

   walberla::math::seedRandomGenerator( 1234 );
   auto rand = []( const Point3D& ) { return walberla::math::realRandom(); };

   u_src.uvw().interpolate( { rand, rand, rand }, maxLevel, All );
   u_src.p().interpolate( rand, maxLevel, All );

   using StokesOperatorFS = hyteg::StrongFreeSlipWrapper< StokesOperatorType, ProjectNormalOperatorType >;
   auto stokes            = std::make_shared< StokesOperatorType >( storage, minLevel, maxLevel );
   auto normalsRect       = []( auto, Point3D& n ) { n = Point3D( { 0, -1 } ); };

   auto projection = std::make_shared< ProjectNormalOperatorType >( storage, minLevel, maxLevel, normalsRect );

   StokesOperatorFS L( stokes, projection, FreeslipBoundary );

   L.apply( u_src, u_dst_hyteg, maxLevel, All );

   PETScSparseMatrix< StokesOperatorType > petscMatStokes( "stokes_pure" );
   petscMatStokes.createMatrixFromOperator( *stokes, maxLevel, numerator );
   // petscMatStokes.print( "/tmp/stokes.m" );

   PETScSparseMatrix< StokesOperatorFS > petscMatFS( "stokes_fs" );
   petscMatFS.createMatrixFromOperator( L, maxLevel, numerator );
   // petscMatFS.print( "/tmp/free_slip.m" );

   PETScVector< real_t, StokesFunctionType::template FunctionType > srcVectorPetsc( u_src, numerator, maxLevel );
   PETScVector< real_t, StokesFunctionType::template FunctionType > dstVectorPetsc( u_dst_petsc, numerator, maxLevel );

   MatMult( petscMatFS.get(), srcVectorPetsc.get(), dstVectorPetsc.get() );

   dstVectorPetsc.createFunctionFromVector( u_dst_petsc, numerator, maxLevel );

   diff.assign( { 1, -1 }, { u_dst_hyteg, u_dst_petsc }, maxLevel, All );
   auto norm = sqrt( diff.dotGlobal( diff, maxLevel, All ) );

   const bool outputVTK = false;

   if ( outputVTK )
   {
      VTKOutput vtk( "../../output", "FreeslipPetscApplyTest", storage );
      vtk.add( u_src );
      vtk.add( diff );
      vtk.add( u_dst_hyteg );
      vtk.add( u_dst_petsc );
      vtk.write( maxLevel );
   }

   WALBERLA_CHECK_LESS( norm, absErrorTolerance );
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager manager( &argc, &argv );

   WALBERLA_LOG_INFO_ON_ROOT( "free-slip PETSc assembly P2-P1-TH test" );
   run< P2P1TaylorHoodFunction< real_t >, // function type
        P2P1TaylorHoodStokesOperator,     // operator
        P2ProjectNormalOperator           // projection
        >( 1e-13 );
}
