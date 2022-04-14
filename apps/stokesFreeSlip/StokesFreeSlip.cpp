/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <hyteg/p2functionspace/P2ProjectNormalOperator.hpp>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1BlendingStokesOperator.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
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
   Point2D left( {-channelLength / 2, 0} );
   Point2D right( {channelLength / 2, channelHeight} );

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

template < typename StokesFunction >
void applyDirichletBCRectangle( const double channelLength, const double channelHeight, const uint_t level, StokesFunction& u )
{
   auto inflow                = [=]( auto p ) { return p[0] <= -channelLength / 2 + 1e-14; };
   auto dirichletInterpolantX = [=]( auto p ) {
      return inflow( p ) ? ( channelHeight - p[1] ) * ( channelHeight + p[1] ) : 0;
      // return inflow(p) ? (channelHeight-p[1])*p[1] : 0;
   };

   u.uvw()[0].interpolate( dirichletInterpolantX, level, DirichletBoundary );
   u.uvw()[1].interpolate( 0, level, DirichletBoundary );
}

std::shared_ptr< SetupPrimitiveStorage >
    setupStorageAnnulus( const double rmin, const double rmax, const uint_t nTan, const uint_t nRad )
{
   hyteg::MeshInfo meshInfo = hyteg::MeshInfo::meshAnnulus( rmin, rmax, MeshInfo::CRISS, nTan, nRad );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( *setupStorage );

   auto rminSq = rmin * rmin;
   auto rmaxSq = rmax * rmax;

   // Boundaries
   auto inflow = [=]( auto p ) {
      auto normSq = p.normSq();
      return rminSq - 1e-12 < normSq && normSq < rminSq + 1e-12;
   };

   auto freeslip = [=]( auto p ) {
      auto normSq = p.normSq();
      return rmaxSq - 1e-12 < normSq && normSq < rmaxSq + 1e-12;
   };

   setupStorage->setMeshBoundaryFlagsOnBoundary( 0, 0, true );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 3, freeslip );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 1, inflow );

   AnnulusMap::setMap( *setupStorage );

   return setupStorage;
}

template < typename StokesFunction >
void applyDirichletBCAnnulus( const double rmin, const double, const uint_t level, StokesFunction& u )
{
   auto rminSq = rmin * rmin;

   // Boundaries
   auto inflow = [=]( auto p ) {
      auto normSq = p.normSq();
      return rminSq - 1e-12 < normSq && normSq < rminSq + 1e-12;
   };

   auto dirichletInterpolantX = [=]( auto p ) { return inflow( p ) ? -p[1] / rmin : 0; };
   auto dirichletInterpolantY = [=]( auto p ) { return inflow( p ) ? +p[0] / rmin : 0; };

   u.uvw().interpolate( { dirichletInterpolantX, dirichletInterpolantY }, level, DirichletBoundary );
}

template < typename StokesFunctionType,
           typename StokesFunctionTag,
           typename StokesOperatorType,
           typename ProjectNormalOperatorType >
void run( std::shared_ptr< walberla::config::Config > cfg )
{
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   if ( mainConf.getParameter< bool >( "printParameters" ) )
      mainConf.listParameters();

   // solver parameters
   const uint_t minLevel               = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel               = mainConf.getParameter< uint_t >( "maxLevel" );
   const real_t absoluteTargetResidual = mainConf.getParameter< real_t >( "absoluteTargetResidual" );
   const uint_t maxIterations          = mainConf.getParameter< uint_t >( "maxIterations" );

   // geometry rectangle
   real_t       channelLength = mainConf.getParameter< real_t >( "channelLength" );
   real_t       channelHeight = mainConf.getParameter< real_t >( "channelHeight" );
   const uint_t ny            = mainConf.getParameter< uint_t >( "ny" );

   // geometry annulus
   real_t rmin = mainConf.getParameter< real_t >( "rmin" );
   real_t rmax = mainConf.getParameter< real_t >( "rmax" );
   uint_t nTan = mainConf.getParameter< uint_t >( "nTan" );
   uint_t nRad = mainConf.getParameter< uint_t >( "nRad" );

   uint_t solverType = mainConf.getParameter< uint_t >( "solverType" );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   bool useRectangleScenario = mainConf.getParameter< bool >( "useRectangleScenario" );

   std::shared_ptr< SetupPrimitiveStorage > setupStorage = nullptr;
   if ( useRectangleScenario )
   {
      setupStorage = setupStorageRectangle( channelLength, channelHeight, ny );
   }
   else
   {
      setupStorage = setupStorageAnnulus( rmin, rmax, nTan, nRad );
   }

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, timingTree );

   if ( mainConf.getParameter< bool >( "printGlobalStorageInfo" ) )
   {
      auto globalInfo = storage->getGlobalInfo();
      WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
   }

   if ( mainConf.getParameter< bool >( "writeDomainVTK" ) )
   {
      hyteg::writeDomainPartitioningVTK( storage, "./output", "StokesFreeSlip_domain" );
   }

   StokesFunctionType f( "f", storage, minLevel, maxLevel );

   f.interpolate( 0, maxLevel, All );

   StokesFunctionType u( "u", storage, minLevel, maxLevel );

   if ( useRectangleScenario )
      applyDirichletBCRectangle( channelLength, channelHeight, maxLevel, u );
   else
      applyDirichletBCAnnulus( rmin, rmax, maxLevel, u );

   if ( mainConf.getParameter< bool >( "printDoFCount" ) )
   {
      uint_t totalGlobalDofsStokes = 0;
      for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         uint_t tmpDofStokes = numberOfGlobalDoFs< StokesFunctionTag >( *storage, lvl );
         WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
         totalGlobalDofsStokes += tmpDofStokes;
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );
   }

   hyteg::VTKOutput vtkOutput( "./output", "StokesFreeSlip", storage );
   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.add( u );
   }

   // using StokesOperator = hyteg::StrongFreeSlipWrapper< hyteg::P1P1StokesOperator, hyteg::P1ProjectNormalOperator >;
   using StokesOperatorFS = hyteg::StrongFreeSlipWrapper< StokesOperatorType, ProjectNormalOperatorType >;
   auto stokes            = std::make_shared< StokesOperatorType >( storage, minLevel, maxLevel );
   auto normalsRect       = []( auto, Point3D& n ) { n = Point3D( {0, -1} ); };
   auto normalsAnn        = [=]( Point3D p, Point3D& n ) { n = Point3D( {p[0] / rmax, p[1] / rmax} ); };

   std::shared_ptr< ProjectNormalOperatorType > projection = nullptr;
   if ( useRectangleScenario )
   {
      projection = std::make_shared< ProjectNormalOperatorType >( storage, minLevel, maxLevel, normalsRect );
   }
   else
   {
      projection = std::make_shared< ProjectNormalOperatorType >( storage, minLevel, maxLevel, normalsAnn );
   }

   StokesOperatorFS L( stokes, projection, FreeslipBoundary );

   std::shared_ptr< Solver< StokesOperatorFS > > solver;

   if ( solverType == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: MINRES + pressure prec. (HyTeG)" );
      solver = hyteg::solvertemplates::stokesMinResSolver< StokesOperatorFS >(
          storage, maxLevel, absoluteTargetResidual, maxIterations );
   }
   else if ( solverType == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: MINRES + no prec. (PETSc)" );
#ifdef HYTEG_BUILD_WITH_PETSC
      solver =
          std::make_shared< PETScMinResSolver< StokesOperatorFS > >( storage, maxLevel, 1e-30, absoluteTargetResidual, maxIterations );
#else
      WALBERLA_ABORT( "PETSc not activated." );
#endif
   }

   // print info about convergence
   if ( auto s = std::dynamic_pointer_cast< MinResSolver< StokesOperatorFS > >( solver ) )
      s->setPrintInfo( true );

   solver->solve( L, u, f, maxLevel );

   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write( maxLevel );
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager manager( &argc, &argv );
#endif

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./StokesFreeSlip.prm" );
   }
   else
   {
      cfg = env.config();
   }

   if ( cfg->getBlock( "Parameters" ).getParameter< bool >( "runP1P1" ) )
   {
      run< P1StokesFunction< real_t >,      // function type
           P1StokesFunction< real_t >::Tag, // tag
           P1BlendingStokesOperator,        // operator
           P1ProjectNormalOperator          // projection
           >( cfg );
   }
   else
   {
      run< P2P1TaylorHoodFunction< real_t >,      // function type
           P2P1TaylorHoodFunction< real_t >::Tag, // tag
           P2P1TaylorHoodStokesOperator,          // operator
           P2ProjectNormalOperator                // projection
           >( cfg );
   }
}
