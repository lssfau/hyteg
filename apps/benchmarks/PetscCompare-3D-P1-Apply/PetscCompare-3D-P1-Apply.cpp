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
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"
#include "core/math/Constants.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
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

   PETScManager petscManager( &argc, &argv );

   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      auto defaultFile = "./PetscCompare-3D-P1-Apply.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );
   const uint_t                        level    = mainConf.getParameter< uint_t >( "level" );

   wcTimingTreeApp.start( "Mesh setup + load balancing" );

   std::shared_ptr< hyteg::MeshInfo > meshInfo;
   if( mainConf.getParameter< bool >( "useMeshFile" ) )
   {
      std::string meshFileName = mainConf.getParameter< std::string >( "mesh" );
      meshInfo                 = std::make_shared< hyteg::MeshInfo >( hyteg::MeshInfo::fromGmshFile( meshFileName ) );
   } else
   {
      uint_t numberOfFaces = mainConf.getParameter< uint_t >( "numberOfFaces" );
      if( mainConf.getParameter< bool >( "facesTimesProcs" ) )
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

   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };
   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      //return 5.0;
      return std::sin( walberla::math::pi * xx[0] ) + std::cos( walberla::math::pi * xx[1] );
      //return ( real_c(std::rand()) / real_c(RAND_MAX));
   };

   wcTimingTreeApp.start( "Function allocation" );
   hyteg::P1Function< double >   oneFunc( "x", storage, level, level );
   hyteg::P1Function< double >   x( "x", storage, level, level );
   hyteg::P1Function< double >   y( "y", storage, level, level );
   hyteg::P1Function< double >   z( "z", storage, level, level );
   hyteg::P1Function< double >   diff( "diff", storage, level, level );
   hyteg::P1Function< idx_t >    numerator( "numerator", storage, level, level );
   wcTimingTreeApp.stop( "Function allocation" );

   const uint_t totalDoFs = numberOfGlobalDoFs< hyteg::P1FunctionTag >( *storage, level );

   wcTimingTreeApp.start( "Operator assembly" );
   hyteg::P1ConstantLaplaceOperator mass( storage, level, level );
   wcTimingTreeApp.stop( "Operator assembly" );

   wcTimingTreeApp.start( "Interpolation" );
   x.interpolate( exact, level, hyteg::Inner );
   oneFunc.interpolate( ones, level );
   wcTimingTreeApp.stop( "Interpolation" );

   wcTimingTreeApp.start( "Enumeration" );
   const uint_t globalDoFs = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level );
   const uint_t localDoFs  = hyteg::numberOfLocalDoFs< P1FunctionTag >( *storage, level );
   numerator.enumerate( level );
   wcTimingTreeApp.stop( "Enumeration" );

   LIKWID_MARKER_START( "PETSc-setup" );
   wcTimingTreeApp.start( "Petsc setup" );
   hyteg::PETScSparseMatrix< hyteg::P1ConstantLaplaceOperator, hyteg::P1Function > matPetsc( localDoFs, globalDoFs );
   matPetsc.createMatrixFromOperator( mass, level, numerator, hyteg::Inner );
   hyteg::PETScVector< real_t, hyteg::P1Function > vecPetsc( localDoFs );
   vecPetsc.createVectorFromFunction( x, numerator, level, hyteg::Inner );
   hyteg::PETScVector< real_t, hyteg::P1Function > dstvecPetsc( localDoFs );
   wcTimingTreeApp.stop( "Petsc setup" );
   LIKWID_MARKER_STOP( "PETSc-setup" );

   wcTimingTreeApp.start( "HyTeG apply" );
   LIKWID_MARKER_START( "HyTeG-apply" );
   mass.apply( x, y, level, hyteg::Inner );
   LIKWID_MARKER_STOP( "HyTeG-apply" );
   wcTimingTreeApp.stop( "HyTeG apply" );

   //vecPetsc.print("../output/vector0.vec");
   PetscErrorCode ierr;

   wcTimingTreeApp.start( "Petsc apply" );
   LIKWID_MARKER_START( "Petsc-MatMult" );
   ierr = MatMult( matPetsc.get(), vecPetsc.get(), dstvecPetsc.get() );
   LIKWID_MARKER_STOP( "Petsc-MatMult" );
   wcTimingTreeApp.stop( "Petsc apply" );

   CHKERRQ( ierr );

   dstvecPetsc.createFunctionFromVector( z, numerator, level, hyteg::Inner );

   // WALBERLA_LOG_INFO_ON_ROOT( y.dotGlobal( oneFunc, level, hyteg::Inner ) );
   // WALBERLA_LOG_INFO_ON_ROOT( z.dotGlobal( oneFunc, level, hyteg::Inner ) );

   //dstvecPetsc.print("../output/vector1.vec");

   if( mainConf.getParameter< bool >( "printTiming" ) )
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

   diff.assign( {1.0, -1.0}, {z, y}, level, hyteg::All );

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "writing VTK output" );
      hyteg::VTKOutput vtkOutput("./output", "PetscCompare-2D-P2-Apply", storage);
      vtkOutput.add( x );
      vtkOutput.add( z );
      vtkOutput.add( y );
      vtkOutput.add( diff );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_FLOAT_EQUAL( y.dotGlobal( oneFunc, level, hyteg::Inner ), z.dotGlobal( oneFunc, level, hyteg::Inner ) )

   WALBERLA_LOG_INFO_ON_ROOT( std::scientific << " | " << numberOfFaces << " | " << level << " | " << totalDoFs << " | "
                                              << walberla::MPIManager::instance()->numProcesses() << " | "
                                              << wcTimingTreeApp["HyTeG apply"].last() << " | "
                                              << wcTimingTreeApp["Petsc apply"].last() << " | " );

   LIKWID_MARKER_CLOSE;
}
