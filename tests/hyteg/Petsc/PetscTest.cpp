/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Marcus Mohr.
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
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/CCRStokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2PlusBubbleFunction.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/types/Concepts.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

std::shared_ptr< PrimitiveStorage > generateStorage( bool use3D )
{
   std::string meshFileName = use3D ? prependHyTeGMeshDir( "3D/cube_24el.msh" ) : prependHyTeGMeshDir( "2D/quad_8el.msh" );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   return storage;
}

void solveP1Laplace()
{
   const size_t level = 2;

   std::shared_ptr< PrimitiveStorage > storage = generateStorage( false );

   hyteg::P1Function< real_t > x( "x", storage, level, level );
   hyteg::P1Function< real_t > x_exact( "x_exact", storage, level, level );
   hyteg::P1Function< real_t > err( "err", storage, level, level );

   hyteg::P1ConstantLaplaceOperator A( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1] + 10;
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   x.interpolate( exact, level, hyteg::DirichletBoundary );
   x_exact.interpolate( exact, level );

   uint_t globalDoFs = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "Num dofs = " << uint_c( globalDoFs ) )

   PETScLUSolver< hyteg::P1ConstantLaplaceOperator > solver( storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Solving System" )
   walberla::WcTimer timer;
   solver.solve( A, x, x, level );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_CHECK_LESS( discr_l2_err, 1e-14 );
}

template < template < typename > class func_t >
   requires concepts::fe_function< func_t< real_t > >
void conversionTest( std::string funcDescriptor, bool run3D = true )
{
   WALBERLA_LOG_INFO_ON_ROOT( "* Checking " << funcDescriptor << " [in " << ( run3D ? "3D" : "2D" ) << "]" );

   std::shared_ptr< PrimitiveStorage > storage = generateStorage( false );

   const size_t level = 2;

   // setup FE functions
   func_t< real_t > feFunctionBefore( "before conversion", storage, level, level );
   func_t< real_t > feFunctionAfter( "after conversion", storage, level, level );
   func_t< real_t > difference( "difference", storage, level, level );

   // initialise original function with random values
   WALBERLA_LOG_INFO_ON_ROOT( "  --> preparing FEFunction" );
   walberla::math::seedRandomGenerator( 1234 );
   auto rand = []( const Point3D& ) { return real_c( walberla::math::realRandom() ); };
   feFunctionBefore.interpolate( rand, level, All );

   // prepare vector conversion
   func_t< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   // convert function to vector
   WALBERLA_LOG_INFO_ON_ROOT( "  --> converting to PETScVector" );
   PETScVector< real_t, func_t > petscVector;
   petscVector.createVectorFromFunction( feFunctionBefore, numerator, level, All );

   // convert vector to function
   WALBERLA_LOG_INFO_ON_ROOT( "  --> converting back " << funcDescriptor );
   petscVector.createFunctionFromVector( feFunctionAfter, numerator, level, All );

   // compare the two versions
   WALBERLA_LOG_INFO_ON_ROOT( "  --> comparing functions before and after" );
   difference.assign( { real_c( 1 ), real_c( -1 ) }, { feFunctionBefore, feFunctionAfter }, level, All );
   real_t error = difference.dotGlobal( difference, level, All );
   WALBERLA_LOG_INFO_ON_ROOT( "  --> error measure is " << error );
   WALBERLA_CHECK_FLOAT_EQUAL( error, real_c( 0 ) );
}

template < template < typename > class func_t >
   requires concepts::fe_function< func_t< real_t > >
void saveIdentityOperatorTest( std::string funcDescriptor )
{
   WALBERLA_LOG_INFO_ON_ROOT( "* Checking identity operator for " << funcDescriptor );

   std::shared_ptr< PrimitiveStorage > storage = generateStorage( false );

   const size_t level = 2;

   // Assemble identity matrix
   WALBERLA_LOG_INFO_ON_ROOT( "  --> assembling identity matrix" );

   func_t< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   MPI_Comm petscCommunicator = walberla::mpi::MPIManager::instance()->comm();

   const uint_t localRows  = numberOfLocalDoFs( numerator, level );
   const uint_t localCols  = numberOfLocalDoFs( numerator, level );
   const uint_t globalRows = numberOfGlobalDoFs( numerator, level, petscCommunicator );
   const uint_t globalCols = numberOfGlobalDoFs( numerator, level, petscCommunicator );

   WALBERLA_LOG_INFO_ON_ROOT( "      (numRows x numCols) = (" << globalRows << " x " << globalCols << ")" );

   Mat matrix;
   MatCreate( petscCommunicator, &matrix );
   MatSetType( matrix, MATMPIAIJ );
   MatSetSizes( matrix,
                static_cast< PetscInt >( localRows ),
                static_cast< PetscInt >( localCols ),
                static_cast< PetscInt >( globalRows ),
                static_cast< PetscInt >( globalCols ) );
   PetscObjectSetName( (PetscObject) matrix, "Identity" );
   MatMPIAIJSetPreallocation( matrix, 1, NULL, 1, NULL );
   MatZeroEntries( matrix );

   auto proxy = std::make_shared< PETScSparseMatrixProxy >( matrix );

   saveIdentityOperator( numerator, proxy, level, All );

   MatAssemblyBegin( matrix, MAT_FINAL_ASSEMBLY );
   MatAssemblyEnd( matrix, MAT_FINAL_ASSEMBLY );

   // store matrix to file
   bool storeMatrix = false;
   if ( storeMatrix )
   {
      std::stringstream fname;
      fname << "IdentityMatrix-" << funcDescriptor << ".txt";
      PetscViewer viewer;
      PetscViewerASCIIOpen( petscCommunicator, fname.str().c_str(), &viewer );
      PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATRIXMARKET );
      MatView( matrix, viewer );
      PetscViewerDestroy( &viewer );
   }

   // setup FE functions
   func_t< real_t > srcFunction( "source", storage, level, level );
   func_t< real_t > dstFunction( "destination", storage, level, level );
   func_t< real_t > difference( "difference", storage, level, level );
   func_t< real_t > aux( "auxilliary", storage, level, level );

   // initialise original function with random values
   WALBERLA_LOG_INFO_ON_ROOT( "  --> preparing FEFunctions and Vectors" );

   walberla::math::seedRandomGenerator( 1234 );
   auto rand = []( const Point3D& ) { return real_c( walberla::math::realRandom() ); };
   srcFunction.interpolate( rand, level, All );

   PETScVector< real_t, func_t > srcVector( srcFunction, numerator, level, All, "src" );
   PETScVector< real_t, func_t > dstVector( aux, numerator, level, All, "dst" );

   WALBERLA_LOG_INFO_ON_ROOT( "  --> multiplying vector with identity matrix" );
   MatMult( matrix, srcVector.get(), dstVector.get() );

   // compare the two versions
   WALBERLA_LOG_INFO_ON_ROOT( "  --> comparing functions before and after" );

   dstVector.createFunctionFromVector( dstFunction, numerator, level, All );

   difference.assign( { real_c( 1 ), real_c( -1 ) }, { srcFunction, dstFunction }, level, All );
   real_t error = difference.getMaxDoFMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "  --> error measure is " << std::scientific << error );
   WALBERLA_CHECK_FLOAT_EQUAL( error, real_c( 0 ) );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " TEST: Convert Function to Vector and Back" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );

   conversionTest< P0Function >( "P0Function" );
   conversionTest< P1Function >( "P1Function" );
   conversionTest< P2Function >( "P2Function" );
   conversionTest< DG1Function >( "DG1Function" );
   conversionTest< P2PlusBubbleFunction >( "P2PlusBubbleFunction", false );
   conversionTest< P2P1TaylorHoodFunction >( "P2P1TaylorHoodFunction" );
   conversionTest< CCRStokesFunction >( "CCRStokesFunction", false );

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " TEST: Assemble Identity Matrix for FE Space" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );

   saveIdentityOperatorTest< P1Function >( "P1Function" );
   saveIdentityOperatorTest< DG1Function >( "DG1Function" );

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " TEST: Solve 2D Laplace Problem with P1 via PETSc" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );

   solveP1Laplace();

   return EXIT_SUCCESS;
}
