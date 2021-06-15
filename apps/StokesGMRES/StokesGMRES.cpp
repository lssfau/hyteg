/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1StokesOperator.hpp"
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
#include "hyteg/solvers/preconditioners/JacobiPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

#include "hyteg/solvers/GMRESSolver.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      auto defaultFile = "./StokesGMRES.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   /////////////// Parameters ///////////////
   const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );
   const walberla::Config::BlockHandle layersParam = cfg->getBlock( "Layers" );

   const uint_t          ntan = mainConf.getParameter< uint_t >( "ntan" );
   std::vector< double > layers;
   for( auto it : layersParam )
   {
      layers.push_back( layersParam.getParameter< double >( it.first ) );
   }

   const double rmin = layers.front();
   const double rmax = layers.back();

   const Point3D sourcePoint  = Point3D( {rmin, 0, 0} ) + 0.5 * Point3D( {rmax - rmin, 0, 0} );
   const real_t  sourceRadius = 0.5;

   const uint_t minLevel            = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel            = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t numVCycle           = mainConf.getParameter< uint_t >( "numVCycle" );
   const uint_t maxKrylowDim        = mainConf.getParameter< uint_t >( "maxKrylowDim" );
   const uint_t restartLength       = mainConf.getParameter< uint_t >( "restartLength" );
   const real_t arnoldiTOL          = mainConf.getParameter< real_t >( "arnoldiTOL" );
   const real_t approxTOL           = mainConf.getParameter< real_t >( "approxTOL" );
   const real_t doubleOrthoTOL      = mainConf.getParameter< real_t >( "doubleOrthoTOL" );

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshSphericalShell( ntan, layers );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );

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

   hyteg::VTKOutput vtkOutput("./output", "StokesGMRES", storage);
   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.add( u.uvw );
      vtkOutput.add( u.p );
      vtkOutput.add( f.uvw );
      vtkOutput.add( f.p );
   }

   hyteg::P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > rhsPlumeX = [sourcePoint, sourceRadius]( const hyteg::Point3D& x ) {
      const real_t distToSourcePoint = ( x - sourcePoint ).norm();
      if( distToSourcePoint < sourceRadius )
         return x[0] * ( sourceRadius - distToSourcePoint );
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > rhsPlumeY = [sourcePoint, sourceRadius]( const hyteg::Point3D& x ) {
      const real_t distToSourcePoint = ( x - sourcePoint ).norm();
      if( distToSourcePoint < sourceRadius )
         return x[1] * ( sourceRadius - distToSourcePoint );
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > rhsPlumeZ = [sourcePoint, sourceRadius]( const hyteg::Point3D& x ) {
      const real_t distToSourcePoint = ( x - sourcePoint ).norm();
      if( distToSourcePoint < sourceRadius )
         return x[2] * ( sourceRadius - distToSourcePoint );
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   f.uvw.interpolate( {rhsPlumeX, rhsPlumeY, rhsPlumeZ}, maxLevel );

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   std::string solverType = mainConf.getParameter< std::string >( "solver" );

   if( solverType == "plain" )
   {
      typedef hyteg::GMRESSolver< hyteg::P1StokesOperator  >GMRESsolv;
      auto KleinerFeinerGMRES = GMRESsolv( storage, minLevel, maxLevel, maxKrylowDim, restartLength, arnoldiTOL, approxTOL, doubleOrthoTOL);
      KleinerFeinerGMRES.solve( L, u, f, maxLevel );

      L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, r}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      real_t currentResidualL2 = sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) );
      WALBERLA_LOG_INFO_ON_ROOT("currentResidualL2 : " << currentResidualL2);
   } else if( solverType == "simplePrec" )
   {
      /// A block Preconditioner for GMRES /////
      typedef StokesBlockDiagonalPreconditioner< hyteg::P1StokesOperator, hyteg::P1LumpedInvMassOperator >Preconditioner_T;

      auto prec = std::make_shared< Preconditioner_T >( storage, minLevel, maxLevel, 5 );

      typedef hyteg::GMRESSolver< hyteg::P1StokesOperator  >GMRESsolv;
      auto KleinerFeinerGMRES = GMRESsolv( storage, minLevel, maxLevel, maxKrylowDim, restartLength, arnoldiTOL, approxTOL, doubleOrthoTOL, prec );
      KleinerFeinerGMRES.solve( L, u, f, maxLevel );

      L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, r}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      real_t currentResidualL2 = sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) );
      WALBERLA_LOG_INFO_ON_ROOT("currentResidualL2 : " << currentResidualL2);
   } else if( solverType == "multigrid") 
   {
      typedef StokesPressureBlockPreconditioner< hyteg::P1StokesOperator, hyteg::P1LumpedInvMassOperator >PressurePreconditioner_T;
      auto pressurePrec = std::make_shared< PressurePreconditioner_T>( storage, minLevel, minLevel );

      typedef hyteg::GMRESSolver< hyteg::P1StokesOperator >PressurePreconditionedGMRES_T;
      auto pressurePreconditionedGMRESSolver = std::make_shared< PressurePreconditionedGMRES_T >(
          storage, minLevel, maxLevel, maxKrylowDim, restartLength, arnoldiTOL, approxTOL, doubleOrthoTOL, pressurePrec );

      typedef GeometricMultigridSolver< hyteg::P1StokesOperator >GMRESMultigrid_T;

      auto stokesRestriction = std::make_shared< hyteg::P1P1StokesToP1P1StokesRestriction>();
      auto stokesProlongation = std::make_shared< hyteg::P1P1StokesToP1P1StokesProlongation >();
      auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1StokesOperator::VelocityOperator_T > >();
      auto multigridVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< hyteg::P1StokesOperator > >( storage, gaussSeidel );

      auto GMRESMultigridSmoother = std::make_shared< hyteg::UzawaSmoother< P1StokesOperator > >(storage, multigridVelocityPreconditioner,minLevel, maxLevel, 0.3);

      GMRESMultigrid_T GMRESsolv(
          storage, GMRESMultigridSmoother, pressurePreconditionedGMRESSolver, stokesRestriction, stokesProlongation, minLevel, maxLevel, 2, 2, 2 );

      auto count = hyteg::Function< hyteg::vertexdof::VertexDoFFunction< real_t > >::getLevelWiseFunctionCounter();
      if( mainConf.getParameter< bool >( "printFunctionCount" ) ) {
         for (uint_t i = minLevel; i <= maxLevel; ++i) {
            WALBERLA_LOG_INFO_ON_ROOT("Total number of P1 Functions on " << i << " : " << count[i]);
         }
      }
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
         GMRESsolv.solve( L, u, f, maxLevel );

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
   } else 
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Unkown solver type" );
   }

#if 0
  auto numerator = std::make_shared< hyteg::P1StokesFunction< PetscInt > >( "numerator", storage, level, level );
   uint_t globalSize = 0;
   const uint_t localSize = numerator->enumerate(level, globalSize);
   PETScManager petscManager( &argc, &argv );
   PETScLUSolver< real_t, hyteg::P1StokesFunction, hyteg::P1StokesOperator > petScLUSolver( numerator, localSize, globalSize );
   f.u.assign( {1.0}, {&u.u}, level, DirichletBoundary );
   f.v.assign( {1.0}, {&u.v}, level, DirichletBoundary );
   f.w.assign( {1.0}, {&u.w}, level, DirichletBoundary );
   petScLUSolver.solve( L, u, f, r, level, uzawaTolerance, maxIterations, Inner | NeumannBoundary );
#endif
   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write( maxLevel, 1 );
   }

   if( mainConf.getParameter< bool >( "PrintTiming" ) ) {
      auto tt = timingTree->getReduced();
      //19.07.2018 this is not in walberla master yet
      //auto tt = timingTree->getCopyWithRemainder();
      WALBERLA_LOG_INFO_ON_ROOT(tt);
   }

   return EXIT_SUCCESS;
}
