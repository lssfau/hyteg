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
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./ConvectionCell2D.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   /////////////// Parameters ///////////////
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   if ( mainConf.getParameter< bool >( "printParameters" ) )
   {
      mainConf.listParameters();
   }

   const uint_t minLevel   = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel   = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t numVCycles = mainConf.getParameter< uint_t >( "numVCycles" );
   const real_t dt         = mainConf.getParameter< real_t >( "dt" );
   const uint_t stepsTotal = mainConf.getParameter< uint_t >( "timesteps" );

   //////////////////////////////////////////

   MeshInfo meshInfo     = hyteg::MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 2, 2 );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 3 );

   storage->getTimingTree()->start( "Total" );

   if ( mainConf.getParameter< bool >( "printGlobalStorageInfo" ) )
   {
      auto globalInfo = storage->getGlobalInfo();
      WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
   }

   if ( mainConf.getParameter< bool >( "writeDomainVTK" ) )
   {
      hyteg::writeDomainPartitioningVTK( storage, "./output", "ConvectionCell2D_Domain" );
   }

   hyteg::P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > uLast( "uLast", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );

   P2Function< real_t > c( "c", storage, minLevel, maxLevel );
   P2Function< real_t > cMass( "cError", storage, minLevel, maxLevel );

   if ( mainConf.getParameter< bool >( "printDoFCount" ) )
   {
      uint_t totalGlobalDofsStokes = 0;
      for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         uint_t tmpDofStokes = numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, lvl );
         WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
         totalGlobalDofsStokes += tmpDofStokes;
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );
   }

   hyteg::VTKOutput vtkOutput( "./output", "ConvectionCell2D", storage );
   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.add( u );
      vtkOutput.add( f );
      vtkOutput.add( c );
   }

   std::function< real_t( const hyteg::Point3D& ) > temperature = []( const hyteg::Point3D& x ) {
      real_t temp_ = 1.0 - std::pow( x[1], 0.5 );
      return temp_ + 0.1 * ( 1.0 - x[1] ) * ( std::sin( walberla::math::pi * x[0] ) );
   };

   const real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " - dt:           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - h (min, max): " << hMin << ", " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - level:        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - time steps:   " << stepsTotal )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::P2P1TaylorHoodStokesOperator   L( storage, minLevel, maxLevel );
   P2ConstantMassOperator                M( storage, minLevel, maxLevel );
   MMOCTransport< P2Function< real_t > > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   auto pressurePreconditioner = std::make_shared<
       hyteg::StokesPressureBlockPreconditioner< hyteg::P2P1TaylorHoodStokesOperator, hyteg::P1LumpedInvMassOperator > >(
       storage, minLevel, maxLevel );
   auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< P2P1TaylorHoodStokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< P2P1TaylorHoodStokesOperator > >( storage,
                                                                                                                  gaussSeidel );
   auto smoother = std::make_shared< hyteg::UzawaSmoother< hyteg::P2P1TaylorHoodStokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, 0.37 );
   auto coarseGridSolver = std::make_shared< hyteg::MinResSolver< hyteg::P2P1TaylorHoodStokesOperator > >(
       storage, minLevel, minLevel, 1000, 1e-12, pressurePreconditioner );
   auto restrictionOperator  = std::make_shared< hyteg::P2P1StokesToP2P1StokesRestriction >();
   auto prolongationOperator = std::make_shared< hyteg::P2P1StokesToP2P1StokesProlongation >();

   auto gmgSolver = std::make_shared< hyteg::GeometricMultigridSolver< hyteg::P2P1TaylorHoodStokesOperator > >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );
   auto gmgLoop = hyteg::SolverLoop< hyteg::P2P1TaylorHoodStokesOperator >( gmgSolver, numVCycles );

   c.interpolate( temperature, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " outer step | timestep | max temperature | total mass | mass lost since last outer step " )
   WALBERLA_LOG_INFO_ON_ROOT( "------------+----------+-----------------+------------+---------------------------------" )

   auto max_temp = c.getMaxMagnitude( maxLevel, All );
   M.apply( c, cMass, maxLevel, All );
   auto total_mass = cMass.sumGlobal( maxLevel );

   vtkOutput.write( maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %10d | %8d | %15.3e | %10.3e | %30.2f%% ", 0, 0, max_temp, total_mass, 0, 0. ) )

   for ( uint_t i = 1; i <= stepsTotal; i++ )
   {
      uLast.assign( { 1.0 }, { u }, maxLevel, All );
      M.apply( c, f.uvw()[1], maxLevel, All );
      f.uvw()[1].assign( { mainConf.getParameter< real_t >( "convectivity" ) }, { f.uvw()[1] }, maxLevel, All );

      gmgLoop.solve( L, u, f, maxLevel );

      transport.step( c, u.uvw(), uLast.uvw(), maxLevel, Inner, dt, 1, true );

      max_temp = c.getMaxMagnitude( maxLevel, All );
      M.apply( c, cMass, maxLevel, All );
      auto total_mass_new  = cMass.sumGlobal( maxLevel );
      auto total_mass_lost = 1.0 - ( total_mass_new / total_mass );

      total_mass = total_mass_new;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( " %10d | %8d | %15.3e | %10.3e | %30.2f%% ", i, i, max_temp, total_mass, total_mass_lost * 100. ) )

      vtkOutput.write( maxLevel, i );
   }

   storage->getTimingTree()->stop( "Total" );
   WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );

   return EXIT_SUCCESS;
}
