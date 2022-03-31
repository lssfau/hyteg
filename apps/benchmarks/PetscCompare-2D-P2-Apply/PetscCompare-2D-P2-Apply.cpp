/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Dominik Thoennes, Nils Kohl.
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
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::WcTimingTree wcTimingTreeApp;

   LIKWID_MARKER_THREADINIT;

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./PetscCompare-2D-P2-Apply.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );
   const uint_t                        level    = mainConf.getParameter< uint_t >( "level" );
   const std::string level_string = "-level" + ( level < 10 ? "0" + std::to_string( level ) : std::to_string( level ) );
   const uint_t      iterations   = mainConf.getParameter< uint_t >( "iterations" );

   wcTimingTreeApp.start( "Mesh setup + load balancing" );

   std::shared_ptr< hyteg::MeshInfo > meshInfo;
   if ( mainConf.getParameter< bool >( "useMeshFile" ) )
   {
      std::string meshFileName = mainConf.getParameter< std::string >( "mesh" );
      meshInfo                 = std::make_shared< hyteg::MeshInfo >( hyteg::MeshInfo::fromGmshFile( meshFileName ) );
   }
   else
   {
      uint_t numberOfFaces = mainConf.getParameter< uint_t >( "numberOfFaces" );
      if ( mainConf.getParameter< bool >( "facesTimesProcs" ) )
      {
         meshInfo = std::make_shared< hyteg::MeshInfo >(
             hyteg::MeshInfo::meshFaceChain( numberOfFaces * uint_c( walberla::MPIManager::instance()->numProcesses() ) ) );
      }
   }

   hyteg::SetupPrimitiveStorage setupStorage( *meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   uint_t numberOfFaces = setupStorage.getNumberOfFaces();

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );
   wcTimingTreeApp.stop( "Mesh setup + load balancing" );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      //return 5.0;
      return std::sin( walberla::math::pi * xx[0] ) + std::cos( walberla::math::pi * xx[1] );
      //return ( real_c(std::rand()) / real_c(RAND_MAX));
   };

   wcTimingTreeApp.start( "Function allocation" );
   hyteg::P2Function< double > x( "x", storage, level, level );
   hyteg::P2Function< double > y( "y", storage, level, level );
#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::P2Function< double > z( "z", storage, level, level );
   hyteg::P2Function< idx_t >  numerator( "numerator", storage, level, level );
#endif
   wcTimingTreeApp.stop( "Function allocation" );

   const uint_t totalDoFs = numberOfGlobalDoFs< hyteg::P2FunctionTag >( *storage, level );

   wcTimingTreeApp.start( "Operator assembly" );
   hyteg::P2ConstantLaplaceOperator mass( storage, level, level );
   wcTimingTreeApp.stop( "Operator assembly" );

   wcTimingTreeApp.start( "Interpolation" );
   x.interpolate( exact, level, hyteg::Inner );
   wcTimingTreeApp.stop( "Interpolation" );

#ifdef HYTEG_BUILD_WITH_PETSC
   wcTimingTreeApp.start( "Enumeration" );
   numerator.enumerate( level );
   wcTimingTreeApp.stop( "Enumeration" );

   LIKWID_MARKER_START( "PETSc-setup" );
   wcTimingTreeApp.start( "Petsc setup" );
   hyteg::PETScSparseMatrix< hyteg::P2ConstantLaplaceOperator > matPetsc;
   matPetsc.createMatrixFromOperator( mass, level, numerator, hyteg::Inner );
   hyteg::PETScVector< real_t, hyteg::P2Function > vecPetsc;
   vecPetsc.createVectorFromFunction( x, numerator, level, hyteg::Inner );
   hyteg::PETScVector< real_t, hyteg::P2Function > dstvecPetsc;
   dstvecPetsc.createVectorFromFunction( x, numerator, level, hyteg::Inner );
   wcTimingTreeApp.stop( "Petsc setup" );
   LIKWID_MARKER_STOP( "PETSc-setup" );
#endif
   wcTimingTreeApp.start( "HyTeG apply" );
   LIKWID_MARKER_START( std::string( "HyTeG-apply" + level_string ).c_str() );
   for ( uint_t i = 0; i < iterations; i++ )
   {
      mass.apply( x, y, level, hyteg::Inner );
   }
   LIKWID_MARKER_STOP( std::string( "HyTeG-apply" + level_string ).c_str() );
   wcTimingTreeApp.stop( "HyTeG apply" );

#ifdef HYTEG_BUILD_WITH_PETSC
   //vecPetsc.print("../output/vector0.vec");
   PetscErrorCode ierr;

   wcTimingTreeApp.start( "Petsc apply" );
   LIKWID_MARKER_START( std::string( "Petsc-MatMult" + level_string ).c_str() );
   for ( uint_t i = 0; i < iterations; i++ )
   {
      ierr = MatMult( matPetsc.get(), vecPetsc.get(), dstvecPetsc.get() );
   }
   LIKWID_MARKER_STOP( std::string( "Petsc-MatMult" + level_string ).c_str() );
   wcTimingTreeApp.stop( "Petsc apply" );

   CHKERRQ( ierr );

   dstvecPetsc.createFunctionFromVector( z, numerator, level, hyteg::Inner );
#else
   wcTimingTreeApp.start( "Petsc apply" );
   wcTimingTreeApp.stop( "Petsc apply" );
#endif

   //dstvecPetsc.print("../output/vector1.vec");

   if ( mainConf.getParameter< bool >( "printTiming" ) )
   {
      auto wcTPReduced = wcTimingTreeApp.getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( wcTPReduced );

      walberla::WcTimingTree tt  = timingTree->getReduced();
      auto                   tt2 = tt.getCopyWithRemainder();
      WALBERLA_LOG_INFO_ON_ROOT( tt2 );

      nlohmann::json ttJson;
      walberla::timing::to_json( ttJson, tt2 );
      std::ofstream jsonOutput;
      jsonOutput.open( "TimingTree.json" );
      jsonOutput << ttJson.dump( 4 );
      jsonOutput.close();
   }

   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "writing VTK output" );
      hyteg::VTKOutput vtkOutput( "./output", "PetscCompare-2D-P2-Apply", storage );
      vtkOutput.add( x );
#ifdef HYTEG_BUILD_WITH_PETSC
      vtkOutput.add( z );
#endif
      vtkOutput.add( y );
      vtkOutput.write( level );
   }

#ifdef HYTEG_BUILD_WITH_PETSC
   WALBERLA_CHECK_FLOAT_EQUAL( y.sumGlobal( level, hyteg::Inner ), z.sumGlobal( level, hyteg::Inner ) )
#endif

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "%6s |%6s |%14s |%6s |%6s |%10s |%10s", "faces", "level", "dof", "iter", "procs", "hyteg", "petsc" ) )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6u |%6u |%14u |%6u |%6u |%10.2e |%10.2e",
                                                numberOfFaces,
                                                level,
                                                totalDoFs,
                                                iterations,
                                                walberla::MPIManager::instance()->numProcesses(),
                                                wcTimingTreeApp["HyTeG apply"].last(),
                                                wcTimingTreeApp["Petsc apply"].last() ) )
#ifdef HYTEG_BUILD_WITH_PETSC
   if ( mainConf.getParameter< bool >( "petscMatOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( matPetsc.getInfo() )
   }
   printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );
   WALBERLA_LOG_INFO_ON_ROOT( matPetsc.getInfo() )
#endif
   //printFunctionAllocationInfo( *storage );
   printCurrentMemoryUsage();
   LIKWID_MARKER_CLOSE;
}
