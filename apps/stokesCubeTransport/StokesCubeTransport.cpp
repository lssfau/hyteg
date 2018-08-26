#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/composites/P1Transport.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSolver.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;

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
   const uint_t maxMinResIterations = mainConf.getParameter< uint_t >( "maxMinResIterations" );
   const std::string meshFile       = mainConf.getParameter< std::string >( "meshFile" );

   const real_t uzawaTolerance = mainConf.getParameter< double >( "uzawaTolerance" );
   const uint_t uzawaMaxIter   = mainConf.getParameter< uint_t >( "uzawaMaxIter" );

   //////////////////////////////////////////

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFile );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hhg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   if( mainConf.getParameter< bool >( "useParMETIS" ) )
   {
      hhg::loadbalancing::distributed::parmetis( *storage );
   }

   if( mainConf.getParameter< bool >( "printGlobalStorageInfo" ) )
   {
      auto globalInfo = storage->getGlobalInfo();
      WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
   }

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree( timingTree );

   if( mainConf.getParameter< bool >( "writeDomainVTK" ) )
   {
      hhg::writeDomainPartitioningVTK( storage, "./output", "StokesCubeTransport_domain" );
   }

   hhg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P1Function< real_t >       temp( "temperature", storage, minLevel, maxLevel );

   if( mainConf.getParameter< bool >( "printDoFCount" ) )
   {
      uint_t totalGlobalDofsStokes = 0;
      for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         uint_t tmpDofStokes = numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, lvl );
         WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
         totalGlobalDofsStokes += tmpDofStokes;
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );
   }

   hhg::VTKOutput vtkOutput( "./output", "StokesCubeTransport" );
   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.set3D();
      vtkOutput.add( &u.u );
      vtkOutput.add( &u.v );
      vtkOutput.add( &u.w );
//      vtkOutput.add( &u.p );
//      vtkOutput.add( &f.u );
//      vtkOutput.add( &f.v );
//      vtkOutput.add( &f.w );
//      vtkOutput.add( &f.p );
      vtkOutput.add( &temp );
   }

   hhg::P1StokesOperator L( storage, minLevel, maxLevel );
   hhg::P1MassOperator   M( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > temperature = []( const hhg::Point3D& x )
   {
     real_t temp_ = 1.0 - std::pow(x[2], 0.5);
     return temp_ + 0.1 * (1.0 - x[2]) * (std::sin(walberla::math::PI * x[0]) * std::sin(walberla::math::PI * x[1]));
   };

   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   temp.interpolate( temperature, maxLevel );

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   std::string solverType = mainConf.getParameter< std::string >( "solver" );

   ///// MinRes coarse grid solver for UZAWA /////
   typedef StokesPressureBlockPreconditioner< hhg::P1StokesFunction< real_t >, hhg::P1LumpedInvMassOperator >
       PressurePreconditioner_T;

   P1LumpedInvMassOperator  massOperator( storage, minLevel, minLevel );
   PressurePreconditioner_T pressurePrec( massOperator, storage, minLevel, minLevel );

   typedef hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator, PressurePreconditioner_T >
       PressurePreconditionedMinRes_T;

   auto pressurePreconditionedMinResSolver = PressurePreconditionedMinRes_T( storage, minLevel, minLevel, pressurePrec );

   ///// UZAWA solver /////
   typedef UzawaSolver< hhg::P1StokesFunction< real_t >,
                        hhg::P1StokesOperator,
                        PressurePreconditionedMinRes_T,
                        P1P1StokesToP1P1StokesRestriction,
                        P1P1StokesToP1P1StokesProlongation,
                        false >
       UzawaSolver_T;

   P1P1StokesToP1P1StokesRestriction  stokesRestriction{};
   P1P1StokesToP1P1StokesProlongation stokesProlongation{};

   UzawaSolver_T uzawaSolver(
       storage, pressurePreconditionedMinResSolver, stokesRestriction, stokesProlongation, minLevel, maxLevel, 2, 2, 2 );

   auto count = hhg::Function< hhg::vertexdof::VertexDoFFunction< real_t > >::getFunctionCounter();
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
   uint plotStep = 1;
   uint_t transportStep = 0;

   for (uint_t step = 0; step < steps; ++step)
   {
      M.apply( temp, f.w, maxLevel, All );
      f.w.assign({mainConf.getParameter< real_t >( "convectivity" )}, {&f.w}, maxLevel, All);

      L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      r.assign( {1.0, -1.0}, {&f, &r}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      real_t currentResidualL2 = sqrt( r.dotGlobal( r, maxLevel, hhg::Inner ) ) /
                                 real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, maxLevel ) );
      real_t lastResidualL2 = currentResidualL2;
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] iteration | residual (L2) | convergence rate " );
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] ----------+---------------+------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] "
                                 << std::setw( 9 ) << 0 << " | " << std::setw( 13 ) << std::scientific << currentResidualL2
                                 << " | " << std::setw( 16 ) << std::scientific << currentResidualL2 / lastResidualL2 );
      for( uint_t i = 0; i < numVCycle; i++ )
      {
         uzawaSolver.solve( L,
                            u,
                            f,
                            r,
                            maxLevel,
                            uzawaTolerance,
                            uzawaMaxIter,
                            hhg::Inner | hhg::NeumannBoundary,
                            UzawaSolver_T::CycleType::VCYCLE,
                            false );

         lastResidualL2 = currentResidualL2;
         L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
         r.assign( {1.0, -1.0}, {&f, &r}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
         currentResidualL2 = sqrt( r.dotGlobal( r, maxLevel, hhg::Inner ) ) /
                             real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, maxLevel ) );
         WALBERLA_LOG_INFO_ON_ROOT( "[StokesSphere] "
                                    << std::setw( 9 ) << i + 1 << " | " << std::setw( 13 ) << std::scientific << currentResidualL2
                                    << " | " << std::setw( 16 ) << std::scientific << currentResidualL2 / lastResidualL2 )
         //WALBERLA_LOG_INFO_ON_ROOT( "after it " << i << ": " << std::scientific << residualMG );
      }

      for (uint_t innerSteps = 0; innerSteps < mainConf.getParameter< uint_t >( "innerTransportSteps" ); ++innerSteps) {
         time += dt;
         WALBERLA_LOG_INFO("time = " << time);

         transportOperator.step(temp, u.u, u.v, u.w, maxLevel, Inner, dt, viscosity);
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
