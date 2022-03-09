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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P1Transport.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

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
   if( env.config() == nullptr )
   {
      auto defaultFile = "./StokesCubeTransport.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   /////////////// Parameters ///////////////
   const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );

   if( mainConf.getParameter< bool >( "printParameters" ) )
   {
      mainConf.listParameters();
   }

   const uint_t minLevel            = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel            = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t numVCycle           = mainConf.getParameter< uint_t >( "numVCycle" );
   const std::string meshFile       = mainConf.getParameter< std::string >( "meshFile" );

   const real_t uzawaTolerance = mainConf.getParameter< double >( "uzawaTolerance" );
   const uint_t uzawaMaxIter   = mainConf.getParameter< uint_t >( "uzawaMaxIter" );

   //////////////////////////////////////////

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFile );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   if( mainConf.getParameter< bool >( "useParMETIS" ) )
   {
      hyteg::loadbalancing::distributed::parmetis( *storage );
   }

   if( mainConf.getParameter< bool >( "printGlobalStorageInfo" ) )
   {
      auto globalInfo = storage->getGlobalInfo();
      WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
   }

   if( mainConf.getParameter< bool >( "writeDomainVTK" ) )
   {
      hyteg::writeDomainPartitioningVTK( storage, "./output", "StokesCubeTransport_domain" );
   }

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t >       temp( "temperature", storage, minLevel, maxLevel );

   if( mainConf.getParameter< bool >( "printDoFCount" ) )
   {
      uint_t totalGlobalDofsStokes = 0;
      for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         uint_t tmpDofStokes = numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, lvl );
         WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
         totalGlobalDofsStokes += tmpDofStokes;
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );
   }

   hyteg::VTKOutput vtkOutput("./output", "StokesCubeTransport", storage);
   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.add( u.uvw() );
//      vtkOutput.add( &u.p );
//      vtkOutput.add( &f.u );
//      vtkOutput.add( &f.v );
//      vtkOutput.add( &f.w );
//      vtkOutput.add( &f.p );
      vtkOutput.add( temp );
   }

   hyteg::P1P1StokesOperator L( storage, minLevel, maxLevel );
   hyteg::P1ConstantMassOperator   M( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > temperature = []( const hyteg::Point3D& x )
   {
     real_t temp_ = 1.0 - std::pow(x[2], 0.5);
     return temp_ + 0.1 * (1.0 - x[2]) * (std::sin(walberla::math::pi * x[0]) * std::sin(walberla::math::pi * x[1]));
   };

   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   temp.interpolate( temperature, maxLevel );

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   std::string solverType = mainConf.getParameter< std::string >( "solver" );

   ///// MinRes coarse grid solver for UZAWA /////
   typedef StokesPressureBlockPreconditioner< hyteg::P1P1StokesOperator, hyteg::P1LumpedInvMassOperator >
        PressurePreconditioner_T;
   auto pressurePrec = std::make_shared< PressurePreconditioner_T >( storage, minLevel, minLevel );
   typedef hyteg::MinResSolver< hyteg::P1P1StokesOperator > PressurePreconditionedMinRes_T;
   auto pressurePreconditionedMinResSolver = std::make_shared< PressurePreconditionedMinRes_T >(
       storage, minLevel, minLevel, uzawaMaxIter, uzawaTolerance, pressurePrec );

   ///// UZAWA solver /////
   typedef GeometricMultigridSolver< hyteg::P1P1StokesOperator > UzawaSolver_T;
   auto stokesRestriction  = std::make_shared< hyteg::P1P1StokesToP1P1StokesRestriction >();
   auto stokesProlongation = std::make_shared< hyteg::P1P1StokesToP1P1StokesProlongation >();
   auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< P1P1StokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< P1P1StokesOperator > >( storage, gaussSeidel );
   auto uzawaSmoother =
       std::make_shared< hyteg::UzawaSmoother< P1P1StokesOperator > >( storage, uzawaVelocityPreconditioner, minLevel, maxLevel, 0.3 );

   UzawaSolver_T uzawaSolver( storage,
                              uzawaSmoother,
                              pressurePreconditionedMinResSolver,
                              stokesRestriction,
                              stokesProlongation,
                              minLevel,
                              maxLevel,
                              2,
                              2,
                              2 );

   auto count = hyteg::Function< hyteg::vertexdof::VertexDoFFunction< real_t > >::getLevelWiseFunctionCounter();
   if( mainConf.getParameter< bool >( "printFunctionCount" ) ) {
      for (uint_t i = minLevel; i <= maxLevel; ++i) {
         WALBERLA_LOG_INFO_ON_ROOT("Total number of P1 Functions on " << i << " : " << count[i]);
      }
   }

   P1Transport transportOperator(storage, minLevel, maxLevel);
   real_t time = 0.0;
   const real_t dt = mainConf.getParameter< real_t >( "dt" );
   const real_t plotDt = mainConf.getParameter< real_t >( "plotDt" );
   const real_t viscosity = mainConf.getParameter< real_t >( "viscosity" );
   const uint_t steps = mainConf.getParameter< uint_t >( "timesteps" );
   const uint_t plotFrequency = walberla::uint_c(std::ceil(plotDt / dt));
   uint_t plotStep = 1;
   uint_t transportStep = 0;

   for (uint_t step = 0; step < steps; ++step)
   {
      M.apply( temp, f.uvw()[2], maxLevel, All );
      f.uvw()[2].assign({mainConf.getParameter< real_t >( "convectivity" )}, {f.uvw()[2]}, maxLevel, All);

      L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, r}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      real_t currentResidualL2 = sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) ) /
                                 real_c( hyteg::numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, maxLevel ) );
      real_t lastResidualL2 = currentResidualL2;
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] iteration | residual (L2) | convergence rate " );
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] ----------+---------------+------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] "
                                 << std::setw( 9 ) << 0 << " | " << std::setw( 13 ) << std::scientific << currentResidualL2
                                 << " | " << std::setw( 16 ) << std::scientific << currentResidualL2 / lastResidualL2 );
      for( uint_t i = 0; i < numVCycle; i++ )
      {
         uzawaSolver.solve( L, u, f, maxLevel );

         lastResidualL2 = currentResidualL2;
         L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
         r.assign( {1.0, -1.0}, {f, r}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
         currentResidualL2 = sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) ) /
                             real_c( hyteg::numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, maxLevel ) );
         WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] "
                                    << std::setw( 9 ) << i + 1 << " | " << std::setw( 13 ) << std::scientific << currentResidualL2
                                    << " | " << std::setw( 16 ) << std::scientific << currentResidualL2 / lastResidualL2 )
         //WALBERLA_LOG_INFO_ON_ROOT( "after it " << i << ": " << std::scientific << residualMG );
      }

      for (uint_t innerSteps = 0; innerSteps < mainConf.getParameter< uint_t >( "innerTransportSteps" ); ++innerSteps) {
         time += dt;
         WALBERLA_LOG_INFO("time = " << time);

         transportOperator.step( temp, u.uvw(), maxLevel, Inner, dt, viscosity );
         ++transportStep;

         if( transportStep % plotFrequency == 0 && mainConf.getParameter< bool >( "VTKOutput" ) )
         {
            WALBERLA_LOG_INFO("Writing output...");
            vtkOutput.write( maxLevel, plotStep );
            ++plotStep;
         }
      }
   }

   if( mainConf.getParameter< bool >( "PrintTiming" ) ) {
      auto tt = timingTree->getReduced();
      //19.07.2018 this is not in walberla master yet
      //auto tt = timingTree->getCopyWithRemainder();
      WALBERLA_LOG_INFO_ON_ROOT(tt);
   }

   return EXIT_SUCCESS;
}
