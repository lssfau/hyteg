/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/trilinos/TrilinosSparseMatrix.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/trilinos/TrilinosVector.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

template < typename OperatorType >
void testSparseMatrix()
{
   WALBERLA_LOG_INFO_ON_ROOT( Tpetra::version() )

   const uint_t level = 3;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, 1, 1 );
   auto     setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   OperatorType                       op( storage, level, level );
   P2P1TaylorHoodFunction< real_t >   src( "src", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   dstTrilinos( "dstTrilinos", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   dstHyteg( "dstHyteg", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   error( "error", storage, level, level );
   P2P1TaylorHoodFunction< idx_t >    numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   auto f = []( const Point3D& p ) -> real_t { return std::sin( p[0] ) + 0.5 * p[1]; };

   src.uvw().interpolate( {f, f}, level, All );
   src.p().interpolate( f, level, All );

   trilinos::TrilinosSparseMatrix< OperatorType, P2P1TaylorHoodFunction > matrix( storage, level );
   matrix.assembleSystem( op, numerator );
   auto matrixString = matrix.to_string();

   trilinos::TrilinosVector< P2P1TaylorHoodFunction > vectorSrc( storage, level );
   vectorSrc.fillFromFunction( src, numerator );
   trilinos::TrilinosVector< P2P1TaylorHoodFunction > vectorDst( storage, level );
   vectorDst.fillFromFunction( dstTrilinos, numerator );
   matrix.apply( vectorSrc, vectorDst );

   vectorDst.writeToFunction( dstTrilinos, numerator );

   op.apply( src, dstHyteg, level, All );

   error.assign( {1.0, -1.0}, {dstTrilinos, dstHyteg}, level, All );

   const auto maxMagnitudeU = error.uvw()[0].getMaxMagnitude( level );
   const auto maxMagnitudeV = error.uvw()[1].getMaxMagnitude( level );
   const auto maxMagnitudeP = error.p().getMaxMagnitude( level );

   WALBERLA_CHECK_LESS( maxMagnitudeU, 1e-14 );
   WALBERLA_CHECK_LESS( maxMagnitudeV, 1e-14 );
   WALBERLA_CHECK_LESS( maxMagnitudeP, 1e-14 );
}

#ifdef HYTEG_BUILD_WITH_PETSC
template < typename OperatorType >
void compareSparseMatrixMatlabOutput()
{
   PETScManager manager;

   WALBERLA_LOG_INFO_ON_ROOT( Tpetra::version() )

   const uint_t level = 3;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, 1, 1 );
   auto     setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   OperatorType                       op( storage, level, level );
   P2P1TaylorHoodFunction< idx_t >    numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   trilinos::TrilinosSparseMatrix< OperatorType, P2P1TaylorHoodFunction > trilinosMatrix( storage, level );
   trilinosMatrix.assembleSystem( op, numerator );
   trilinosMatrix.applyDirichletBoundaryConditions( numerator );
   trilinosMatrix.exportToMatlabFormat( "../../output/TrilinosMatlabExport.m", "MyTrilinosMatrix" );

   // PETScSparseMatrix< OperatorType, P2P1TaylorHoodFunction > petscMatrix( storage, level );
   PETScSparseMatrix< OperatorType > petscMatrix;
   petscMatrix.createMatrixFromOperator( op, level, numerator );
   petscMatrix.applyDirichletBC( numerator, level );
   petscMatrix.print( "../../output/PetscMatlabExport.m" );
}
#endif

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "CC" )
   testSparseMatrix< hyteg::P2P1TaylorHoodStokesOperator >();
   WALBERLA_LOG_INFO_ON_ROOT( "element-wise" )
   testSparseMatrix< hyteg::P2P1ElementwiseBlendingStokesOperator >();

#ifdef HYTEG_BUILD_WITH_PETSC
   compareSparseMatrixMatlabOutput< hyteg::P2P1TaylorHoodStokesOperator >();
#endif

   return EXIT_SUCCESS;
}
