/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/composites/MMOCTransport.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"

using walberla::real_c;
using walberla::real_t;

namespace hyteg {

void simulate( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petScManager;

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./StokesSphereTransport.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   /////////////// Parameters ///////////////
   const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );
   const walberla::Config::BlockHandle layersParam = cfg->getBlock( "Layers" );

   const uint_t          ntan = mainConf.getParameter< uint_t >( "ntan" );
   std::vector< double > layers;
   for ( const auto& it : layersParam )
   {
      layers.push_back( layersParam.getParameter< double >( it.first ) );
   }

   const double rmin = layers.front();
   const double rmax = layers.back();

   const uint_t minLevel            = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel            = mainConf.getParameter< uint_t >( "maxLevel" );
   const real_t stokesResidual      = mainConf.getParameter< real_t >( "stokesResidual" );
   const uint_t stokesMaxNumVCycles = mainConf.getParameter< uint_t >( "stokesMaxNumVCycles" );

   const uint_t timeSteps = mainConf.getParameter< uint_t >( "timeSteps" );

   const real_t diffusivity = mainConf.getParameter< real_t >( "diffusivity" );
   const real_t dt          = mainConf.getParameter< real_t >( "dt" );

   const real_t rhsScaleFactor = mainConf.getParameter< real_t >( "rhsScaleFactor" );

   const uint_t VTKOutputFrequency = mainConf.getParameter< uint_t >( "VTKFrequency" );

   /////////////////// Mesh / Domain ///////////////////////

   MeshInfo meshInfo     = MeshInfo::meshSphericalShell( ntan, layers );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( *setupStorage, timingTree );

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   if ( mainConf.getParameter< bool >( "writeDomainVTK" ) )
   {
      writeDomainPartitioningVTK( storage, "./output", "StokesSphereTransport_domain" );
   }

   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2Function< real_t >             temp( "temperature", storage, minLevel, maxLevel );
   P2Function< real_t >             normalX( "normalX", storage, minLevel, maxLevel );
   P2Function< real_t >             normalY( "normalY", storage, minLevel, maxLevel );
   P2Function< real_t >             normalZ( "normalZ", storage, minLevel, maxLevel );

   uint_t totalGlobalDofsStokes = 0;
   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      uint_t tmpDofStokes = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
      totalGlobalDofsStokes += tmpDofStokes;
   }

   VTKOutput vtkOutput( "./output", "StokesSphereTransport", storage, VTKOutputFrequency );
   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.add( u.u );
      vtkOutput.add( u.v );
      vtkOutput.add( u.w );
      // vtkOutput.add( &u.p );
      // vtkOutput.add( &f.u );
      // vtkOutput.add( &f.v );
      // vtkOutput.add( &f.w );
      // vtkOutput.add( &f.p );
      vtkOutput.add( temp );
   }

   P2P1TaylorHoodStokesOperator L( storage, minLevel, maxLevel );
   P2ConstantMassOperator       M( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > temperature = [rmin, rmax]( const Point3D& x ) {
      return std::pow( ( rmax - x.norm() ) / ( rmax - rmin ), 3.0 );
   };

   temp.interpolate( temperature, maxLevel );

   std::function< real_t( const Point3D& ) > zero = []( const Point3D& ) { return 0.0; };
   std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return 1.0; };

   std::function< real_t( const Point3D& ) > nX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > nY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > nZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   normalX.interpolate( nX, maxLevel );
   normalY.interpolate( nY, maxLevel );
   normalZ.interpolate( nZ, maxLevel );

   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   auto coarseGridSolver = std::make_shared< PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > >(
       storage, minLevel, 1e-12, std::numeric_limits< PetscInt >::max(), 1 );
   auto stokesRestriction  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );
   auto stokesProlongation = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P2P1TaylorHoodStokesOperator > >( storage, minLevel, maxLevel, 0.3 );
   uzawaSmoother->estimateAndSetRelaxationParameter( 2, 20 );
   auto gmgSolver = std::make_shared< GeometricMultigridSolver< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, stokesRestriction, stokesProlongation, minLevel, maxLevel, 3, 3, 2 );

   MMOCTransport< P2Function< real_t > > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );
   real_t                                time = 0.0;

   auto calculateResidual = [&]() {
      L.apply( u, r, maxLevel, Inner | NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, r}, maxLevel, Inner | NeumannBoundary );

      real_t res = sqrt( r.dotGlobal( r, maxLevel, All ) ) /
                   real_c( numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel ) );
      return res;
   };

   walberla::WcTimer timer;

   for ( uint_t step = 0; step < timeSteps; ++step )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Time step " << step )

      // Updating right-hand side (Boussinesq approximation)

      M.apply( temp, f.u, maxLevel, All );
      M.apply( temp, f.v, maxLevel, All );
      M.apply( temp, f.w, maxLevel, All );

      f.u.multElementwise( {f.u, normalX}, maxLevel, All );
      f.v.multElementwise( {f.v, normalY}, maxLevel, All );
      f.w.multElementwise( {f.w, normalZ}, maxLevel, All );

      f.u.assign( {rhsScaleFactor}, {f.u}, maxLevel, All );
      f.v.assign( {rhsScaleFactor}, {f.v}, maxLevel, All );
      f.w.assign( {rhsScaleFactor}, {f.w}, maxLevel, All );

      // Stokes solver

      // Check residual, if less than tolerance, quit v-cycle loop early

      timer.start();

      real_t currentResidual = calculateResidual();
      real_t lastResidual    = currentResidual;
      WALBERLA_LOG_INFO_ON_ROOT( " iteration | residual (l2) | convergence rate " );
      WALBERLA_LOG_INFO_ON_ROOT( "-----------+---------------+------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "   initial | %13.3e |                - ", currentResidual ) )

      for ( uint_t i = 0; i < stokesMaxNumVCycles; i++ )
      {
         if ( currentResidual < stokesResidual )
            break;

         gmgSolver->solve( L, u, f, maxLevel );

         currentResidual      = calculateResidual();
         auto convergenceRate = currentResidual / lastResidual;
         lastResidual         = currentResidual;
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9d | %13.3e | %16.3e ", i, currentResidual, convergenceRate ) )
      }

      timer.end();
      WALBERLA_LOG_INFO_ON_ROOT("")
      WALBERLA_LOG_INFO_ON_ROOT( "... " << timer.last() << " seconds" )
      WALBERLA_LOG_INFO_ON_ROOT("")

      // Advection-diffusion

      WALBERLA_LOG_INFO_ON_ROOT( "Advection time step ..." )

      timer.start();

      time += dt;

      WALBERLA_UNUSED( diffusivity );
      transport.step( setupStorage, temp, u.u, u.v, u.w, maxLevel, All, dt, 1, true );

      timer.end();
      WALBERLA_LOG_INFO_ON_ROOT("")
      WALBERLA_LOG_INFO_ON_ROOT( "... " << timer.last() << " seconds" )
      WALBERLA_LOG_INFO_ON_ROOT( "" )


      if ( mainConf.getParameter< bool >( "VTKOutput" ) )
      {
         vtkOutput.write( maxLevel, step + 1 );
      }
   }

   if ( mainConf.getParameter< bool >( "PrintTiming" ) )
   {
      auto tt = timingTree->getReduced();
      tt      = timingTree->getCopyWithRemainder();
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::simulate( argc, argv );
}
