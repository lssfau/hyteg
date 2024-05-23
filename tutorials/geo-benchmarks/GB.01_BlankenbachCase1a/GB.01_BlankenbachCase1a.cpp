/*
 * Copyright (c) 2023 Ponsuganth Ilangovan P
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

#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

/**
 * \page GB.01_BlankenbachCase1a Tutorial GB.01 - Blankenbach Benchmark (Case 1a)
 *
 * \dontinclude tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp
 * 
 * \brief We have demonstrated a plume rising on a cube and convection on annulus before. In this tutorial we will implement 
 * a standard benchmark from Geophysical literature, the Blankenbach problem, specially Case 1a, see
 * <a href="https://doi.org/10.1111/j.1365-246X.1989.tb05511.x">A benchmark comparison for mantle convection codes,
 * Blankenbach et al., 1989, GJI</a>.
 * It involves solving the respective equations on an unit square and calculation of Nusselt numbers and velocity RMS values
 * to verify the code.
 * 
 * \section GB01-BlankenbachBenchmark-GoverningEquations Model and Equations
 * 
 * The governing equations are basically the Stokes and energy equations,
 * 
 * \f{align*}{
 *  -\nabla\cdot\tau + \nabla p &= \text{Ra}T \\
 *  \nabla u &= 0 \\
 *  \frac{\partial T}{\partial t} + u \cdot \nabla T - \nabla \cdot \kappa \nabla T &= 0
 * \f}
 * 
 * The initial conditions for temperature is prescribed as,
 * 
 * \f{equation*}{
 *  T(x, y, t = 0) = (1-y) + A\cos{\pi x}\sin{\pi y}
 * \f}
 * 
 * while freeslip boundary conditions is imposed on all four walls of the square.
 * For the temperature field, a Dirichlet of T = 0 and T = 1 is imposed 
 * at the top and bottom and zero flux on the sides. This induces a single convection cell 
 * in the square and eventually reaches a steady state. Once this state is reached, 
 * the Nusselt number values are calculated and verified.
 *
 * \section GB01-BlankenbachBenchmark-Domain Domain
 * 
 * Through the geometry module in HyTeG, a rectangle mesh can be created with the `meshRectangle`, 
 * which is used to create the unit square mesh. In addition to creating the macro mesh with required subdivisions, 
 * the boundary nodes must be marked properly which will be used downstream when defining the boundary conditions
 * for the finite element functions, which is what is done in the following snippet.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp SetupStorageAndMarkings
 * 
 * \section GB01-BlankenbachBenchmark-BCs Boundary Conditions
 * 
 * The boundary conditions for the Temperature and Velocity are defined appropriately. For temperature, a Dirichlet is defined
 * on the top and bottom and zero flux on the side walls.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp BoundaryConditionsTemperature
 * 
 * For the velocity freeslip boundary conditions must be imposed on all four walls. For this, the coordinates can be seperately
 * defined a Dirichlet and Neumann condition, but here we generally impose freeslip boundary condition.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp BoundaryConditionsVelocity
 * 
 * \section GB01-BlankenbachBenchmark-OpStokes Operators -- Stokes
 * 
 * The application of freeslip boundary condition is done with a project normal operator 
 * which needs the normals defined at the boundaries, with which the
 * velocity field is projected and the normal component is set to zero. Here the normals are defined through a lambda
 * function,
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp NormalsFunctionForFS
 * 
 * And then passed on to the project normal operator
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp ProjectNormalForFS
 * 
 * This is then inturn used in defining the freeslip wrapper which wraps the Stokes operator and uses the projection operator
 * to set the normal components to zero.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp StokesFreeslipOperator
 * 
 * \section GB01-BlankenbachBenchmark-OpEnergy Operators -- Energy
 * 
 * For the solution of the energy equations, we use a particle approach based on the modified method of characteristics (MMOC).
 * Hence the advection operator can be defined with
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp MMOCForTransport
 * 
 * As the operator splitting approach is used, this must be taken care in the time stepping algorithm. We use an implicit Euler
 * scheme for timestepping the diffusion and the MMOC is used to step for advection. The `apply` function of the 
 * `P2TransportTimesteppingOperator` applies the time discretized form of the weak form considering only the diffusion.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp TransportOperatorApply
 * 
 * \section GB01-BlankenbachBenchmark-SolverStokes Solver -- Stokes
 * 
 * Minres solver is used to compute the solution to the Stokes system. The RHS is calculated and the projection is applied 
 * before starting the iterative solve. As the Minres solver calls the wrapped Stokes operator, it is ensured that the 
 * projection is done at every step of the iterative solve, hence strongly applying the freeslip boundary condition.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp StokesSolverLambdaFunction
 * 
 * \section GB01-BlankenbachBenchmark-SolverEnergy Solver -- Energy
 * 
 * The timestep for the energy equation solver is calculated according to the CFL condition using the velocity from the most recent
 * Stokes solution. This is then used for the MMOC solver to step the advection, then also used for the diffusion solver.
 * Here we use CG solver for the timestepping solver as the system is nicely symmetric with only diffusion present.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp TransportSolverLambdaFunction
 * 
 * \section GB01-BlankenbachBenchmark-Timestepping Timestepping
 * 
 * First the temperature \f$ T(x, t_0) \f$ is initialized to the prescribed field according to the above equation.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp TemperatureInitialization
 * 
 * Then the Stokes system is solved to get the initial velocity field, else the program can also be started from a checkpoint.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp InitialSolveOrCheckpoint
 * 
 * This is then used to step the transport solver, both advection and the diffusion parts, with 
 * MMOC and implicit Euler respectively. With this we compute the temperature field at time \f$ t_1 \f$. Then the Stokes system is
 * solved again to obtain the corresponding velocity field at time \f$ t_1 \f$ and so on.
 * 
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp TimeStepLoopStart
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp SolveEnergy
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp SolveStokes
 * \snippet tutorials/geo-benchmarks/GB.01_BlankenbachCase1a/GB.01_BlankenbachCase1a.cpp WriteDataOut
 * 
 * \section GB01-BlankenbachBenchmark-Results Results
 * 
 * 
 * <img src="GB.01_BlankenbachSquareUpdated.png" width="50%" />
 * 
 * [Blankenbach] BRENNER, Susanne C.; DIEGEL, Amanda E.; SUNG, Li-Yeng.
 *               <a href="https://arxiv.org/pdf/1811.03430.pdf">A robust solver for a mixed finite element method for the Cahn–Hilliard equation</a>.
 *               Journal of Scientific Computing, 2018, 77. Jg., Nr. 2, S. 1234-1249.
 * 
 */

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

class P2TransportTimesteppingOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2TransportTimesteppingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                    const uint_t                               minLevel,
                                    const uint_t                               maxLevel,
                                    real_t                                     k )
   : Operator( storage, minLevel, maxLevel )
   , diffusionOperator( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   , k_( k )
   {}

   ///[TransportOperatorApply]
   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      diffusionOperator.apply( src, dst, level, flag, updateType );
      dst.assign( { k_ * dt }, { dst }, level, flag );
      massOperator.apply( src, dst, level, flag, Add );
   }
   ///[TransportOperatorApply]

   void setDt( real_t dt_ ) { dt = dt_; }

 private:
   P2ElementwiseBlendingLaplaceOperator diffusionOperator;
   P2ElementwiseBlendingMassOperator    massOperator;

   real_t dt = 0.01, k_ = 1.0;
};

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./GB.01_BlankenbachCase1a.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   const uint_t nx = mainConf.getParameter< uint_t >( "nx" );
   const uint_t ny = mainConf.getParameter< uint_t >( "ny" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   real_t diffusivity  = mainConf.getParameter< real_t >( "diffusivity" );
   uint_t stokesIter   = mainConf.getParameter< uint_t >( "stokesIter" );
   real_t stokesRelTol = mainConf.getParameter< real_t >( "stokesRelTol" );

   uint_t       transportIter   = mainConf.getParameter< uint_t >( "transportIter" );
   real_t       transportRelTol = mainConf.getParameter< real_t >( "transportRelTol" );
   bool         verbose         = mainConf.getParameter< bool >( "verbose" );
   const real_t Ra              = mainConf.getParameter< real_t >( "RayleighNumber" );
   real_t       cflMax          = mainConf.getParameter< real_t >( "cflMax" );
   uint_t       NTimesteps      = mainConf.getParameter< uint_t >( "NTimesteps" );
   real_t       endTime         = mainConf.getParameter< real_t >( "endTime" );
   uint_t       logNsFreq       = mainConf.getParameter< uint_t >( "logNsFreq" );

   std::string outputPath          = mainConf.getParameter< std::string >( "outputPath" );
   std::string outputFilename      = mainConf.getParameter< std::string >( "outputFilename" );
   bool        logFilenameWithTime = mainConf.getParameter< bool >( "logFilenameWithTime" );
   uint_t      dataWriteFrequency  = mainConf.getParameter< uint_t >( "dataWriteFrequency" );
   bool        useAdios2           = mainConf.getParameter< bool >( "useAdios2" );
   std::string adiosXmlConfig      = mainConf.getParameter< std::string >( "adiosXmlConfig" );

   std::string checkpointFilepath = mainConf.getParameter< std::string >( "checkpointFilepath" );
   std::string checkpointFilename = mainConf.getParameter< std::string >( "checkpointFilename" );

   bool storeCheckpoint     = mainConf.getParameter< bool >( "storeCheckpoint" );
   bool startFromCheckpoint = mainConf.getParameter< bool >( "startFromCheckpoint" );

   uint_t checkpointFreq = mainConf.getParameter< uint_t >( "checkpointFreq" );

   enum BoundaryMarkers
   {
      Bottom = 23,
      Right,
      Left,
      Top,
      Corners
   };

   const real_t threshold = 1e-6;

   std::function< bool( const Point3D& ) > bottomMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[1] ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > rightMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[0] - 1.0 ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > leftMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[0] ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > topMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[1] - 1.0 ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > cornerMarker = [&]( const Point3D& x ) {
      if ( ( bottomMarker( x ) && rightMarker( x ) ) || ( bottomMarker( x ) && leftMarker( x ) ) ||
           ( topMarker( x ) && rightMarker( x ) ) || ( topMarker( x ) && leftMarker( x ) ) )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   ///[SetupStorageAndMarkings]
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISSCROSS, nx, ny );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Bottom, bottomMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Left, leftMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Right, rightMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Top, topMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Corners, cornerMarker );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );
   ///[SetupStorageAndMarkings]

   BoundaryCondition bcTemp, bcVelocity;

   ///[BoundaryConditionsTemperature]
   bcTemp.createAllInnerBC();
   bcTemp.createDirichletBC( "DirichletTopAndBottom", { BoundaryMarkers::Bottom, BoundaryMarkers::Top } );
   bcTemp.createNeumannBC( "NeumannSides", { BoundaryMarkers::Left, BoundaryMarkers::Right } );
   ///[BoundaryConditionsTemperature]

   ///[BoundaryConditionsVelocity]
   bcVelocity.createAllInnerBC();
   bcVelocity.createDirichletBC( "DirichletCorners", { BoundaryMarkers::Corners } );
   bcVelocity.createFreeslipBC(
       "FreeslipAll", { BoundaryMarkers::Bottom, BoundaryMarkers::Top, BoundaryMarkers::Left, BoundaryMarkers::Right } );
   ///[BoundaryConditionsVelocity]

   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > uPrev( "uPrev", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > fStrong( "fStrong", storage, minLevel, maxLevel, bcVelocity );

   P2Function< real_t > T( "T", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > TPrev( "TPrev", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > fT( "fT", storage, minLevel, maxLevel, bcTemp );

   auto stokesOperator = std::make_shared< P2P1ElementwiseBlendingStokesOperator >( storage, minLevel, maxLevel );

   ///[TransportOperatorDeclaration]
   P2TransportTimesteppingOperator transportOperator( storage, minLevel, maxLevel, diffusivity );
   ///[TransportOperatorDeclaration]

   P2ElementwiseBlendingMassOperator       massOperator( storage, minLevel, maxLevel );
   P2ElementwiseBlendingVectorMassOperator vectorMassOperator( storage, minLevel, maxLevel );

   ///[NormalsFunctionForFS]
   std::function< void( const Point3D&, Point3D& ) > normalsFS = [&]( const Point3D& x, Point3D& normal ) {
      if ( rightMarker( x ) )
      {
         normal[0] = 1.0;
         normal[1] = 0.0;
      }
      else if ( leftMarker( x ) )
      {
         normal[0] = -1.0;
         normal[1] = 0.0;
      }
      else if ( topMarker( x ) )
      {
         normal[0] = 0.0;
         normal[1] = 1.0;
      }
      else if ( bottomMarker( x ) )
      {
         normal[0] = 0.0;
         normal[1] = -1.0;
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Probably shouldn't be here!" );
      }
   };
   ///[NormalsFunctionForFS]

   ///[ProjectNormalForFS]
   auto projectionNormal = std::make_shared< P2ProjectNormalOperator >( storage, minLevel, maxLevel, normalsFS );
   ///[ProjectNormalForFS]

   ///[MMOCForTransport]
   auto mmocTransport = MMOCTransport< P2Function< real_t > >( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );
   ///[MMOCForTransport]

   ///[StokesFreeslipOperator]
   StrongFreeSlipWrapper< P2P1ElementwiseBlendingStokesOperator, P2ProjectNormalOperator, true > stokesOperatorFS(
       stokesOperator, projectionNormal, FreeslipBoundary );
   ///[StokesFreeslipOperator]

   MinResSolver< StrongFreeSlipWrapper< P2P1ElementwiseBlendingStokesOperator, P2ProjectNormalOperator, true > > minresSolver(
       storage, minLevel, maxLevel, stokesIter, stokesRelTol );
   CGSolver< P2TransportTimesteppingOperator > transportSolver( storage, minLevel, maxLevel, transportIter, transportRelTol );

   minresSolver.setPrintInfo( verbose );

   const real_t AiniPerturb = 0.05;

   ///[TemperatureInitialization]
   std::function< real_t( const Point3D& ) > TIni = [&]( const Point3D& x ) {
      return ( 1 - x[1] ) + AiniPerturb * std::cos( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] );
   };

   T.interpolate( TIni, maxLevel, All );
   ///[TemperatureInitialization]

   real_t dt = 0.0001;

   real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );

   auto clockTime  = std::time( nullptr );
   auto clockTimeM = *std::localtime( &clockTime );

   if ( logFilenameWithTime )
   {
      std::ostringstream ossVtkName;
      ossVtkName << outputFilename << "_" << std::put_time( &clockTimeM, "%d-%m-%Y_%H-%M-%S" );

      outputFilename = ossVtkName.str();
   }

   ///[ADIOSVTKOutputDeclaration]
   VTKOutput vtkOutput( outputPath, outputFilename, storage );

#ifdef HYTEG_BUILD_WITH_ADIOS2
   AdiosWriter adios2Output( outputPath, outputFilename, adiosXmlConfig, storage );
#endif
   ///[ADIOSVTKOutputDeclaration]

   if ( useAdios2 )
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      adios2Output.add( u );
      adios2Output.add( T );
#else
      WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
   }
   else
   {
      vtkOutput.add( u );
      vtkOutput.add( T );
   }

   std::function< void( uint_t ) > writeDataOut = [&]( uint_t timestep ) {
      if ( useAdios2 )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         adios2Output.write( maxLevel, timestep );
#else
         WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
      }
      else
      {
         vtkOutput.write( maxLevel, timestep );
      }
   };

   ///[StokesSolverLambdaFunction]
   std::function< void() > solveU = [&]() {
      u.interpolate( 0.0, maxLevel, DirichletBoundary );
      fStrong.uvw().component( 1U ).interpolate( Ra, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
      fStrong.uvw().component( 1U ).multElementwise(
          { fStrong.uvw().component( 1U ), T }, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );

      vectorMassOperator.apply( fStrong.uvw(), f.uvw(), maxLevel, All );

      projectionNormal->project( f, maxLevel, FreeslipBoundary );
      minresSolver.solve( stokesOperatorFS, u, f, maxLevel );
   };
   ///[StokesSolverLambdaFunction]

   ///[TransportSolverLambdaFunction]
   std::function< void() > solveT = [&]() {
      real_t vMax = u.uvw().getMaxComponentMagnitude( maxLevel, All );

      dt = cflMax * hMin / vMax;

      T.interpolate( TIni, maxLevel, DirichletBoundary );
      mmocTransport.step( T, u.uvw(), uPrev.uvw(), maxLevel, Inner | NeumannBoundary | FreeslipBoundary, dt, 1 );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "dt = %4.4e", dt ) );

      transportOperator.setDt( dt );

      T.interpolate( TIni, maxLevel, DirichletBoundary );
      massOperator.apply( T, fT, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
      transportSolver.solve( transportOperator, T, fT, maxLevel );
   };
   ///[TransportSolverLambdaFunction]

   ///[InitialSolveOrCheckpoint]
   if ( startFromCheckpoint )
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      std::shared_ptr< AdiosCheckpointImporter > adios2Importer =
          std::make_shared< AdiosCheckpointImporter >( checkpointFilepath, checkpointFilename, adiosXmlConfig );

      adios2Importer->restoreFunction( T, minLevel, maxLevel );
#else
      WALBERLA_ABORT( "ADIOS2 checkpoint requested in prm file but ADIOS2 was not compiled!" );
#endif
   }
   else
   {
      solveU();
   }
   ///[InitialSolveOrCheckpoint]

   uPrev.assign( { 1.0 }, { u }, maxLevel, All );

   writeDataOut( 0U );

   real_t simulationTime = 0.0;

   uint_t iTimestep = 1U;

   ///[TimeStepLoopStart]
   for ( iTimestep = 1U; iTimestep <= NTimesteps; iTimestep++ )
   ///[TimeStepLoopStart]
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Started step %d at time %4.7e", iTimestep, simulationTime ) );

      ///[SolveEnergy]
      solveT();
      ///[SolveEnergy]

      simulationTime += dt;

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Transport done", iTimestep ) );

      uPrev.assign( { 1.0 }, { u }, maxLevel, All );

      ///[SolveStokes]
      solveU();
      ///[SolveStokes]

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Stokes done", iTimestep ) );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Ended step %d", iTimestep ) );

      if ( simulationTime > endTime )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Specified endTime reached!" );
         break;
      }

      ///[WriteDataOut]
      if ( iTimestep % dataWriteFrequency == 0 )
      {
         writeDataOut( iTimestep );
      }
      ///[WriteDataOut]

      if ( iTimestep % logNsFreq == 0 )
      {
         uint_t numDoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Number of global DoFs = %d", numDoFs ) );

         real_t nusseltNumber = nusseltcalc::calculateNusseltNumber2D( T, maxLevel, 0.001, 1e-6, 101 );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "nusseltNumber = %4.7e", nusseltNumber ) );
      }

      if ( iTimestep % checkpointFreq == 0 && storeCheckpoint )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         std::shared_ptr< AdiosCheckpointExporter > adios2Exporter =
             std::make_shared< AdiosCheckpointExporter >( adiosXmlConfig );

         adios2Exporter->registerFunction( T, minLevel, maxLevel );

         adios2Exporter->storeCheckpoint( checkpointFilepath, checkpointFilename );
#else
         WALBERLA_ABORT( "ADIOS2 checkpoint requested in prm file but ADIOS2 was not compiled!" );
#endif
      }
   }

   uint_t numDoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Number of global DoFs = %d", numDoFs ) );

   real_t nusseltNumber = nusseltcalc::calculateNusseltNumber2D( T, maxLevel, 0.001, 1e-6, 101 );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "nusseltNumber = %4.7e", nusseltNumber ) );

   return 0;
}
