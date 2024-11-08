/*
* Copyright (c) 2023 Nils Kohl.
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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/all.h"
#include "core/timing/all.h"

#include "hyteg/eigen/EigenSparseMatrix.hpp"
#include "hyteg/eigen/EigenVector.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

/// Testing the Eigen sparse matrix assembly proxy here.

using namespace hyteg;

using walberla::int_c;
using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

void testMatVec()
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/cube_6el.msh" ) );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t level = 3;

   P1Function< idx_t > numeratorSrc( "numSrc", storage, level, level );
   P1Function< idx_t > numeratorDst( "numDst", storage, level, level );

   numeratorSrc.enumerate( level );
   numeratorDst.enumerate( level );

   communication::syncFunctionBetweenPrimitives( numeratorSrc, level );
   communication::syncFunctionBetweenPrimitives( numeratorDst, level );

   P1ElementwiseLaplaceOperator op( storage, level, level );

   P1Function< real_t > x( "x", storage, level, level );
   P1Function< real_t > y( "y", storage, level, level );
   P1Function< real_t > yFromEigen( "yFromEigen", storage, level, level );
   P1Function< real_t > err( "err", storage, level, level );

   auto rand = []( const Point3D& ) { return walberla::math::realRandom(); };

   x.interpolate( rand, level );

   op.apply( x, y, level, All );

   Eigen::SparseMatrix< real_t > AEigen = createEigenSparseMatrixFromOperator( op, level, numeratorSrc, numeratorDst );

   VectorXr xEigen = createEigenVectorFromFunction( x, numeratorSrc, level );
   VectorXr yEigen = AEigen * xEigen;

   assignFunctionFromEigenVector( yEigen, yFromEigen, numeratorDst, level );

   err.assign( { 1.0, -1.0 }, { y, yFromEigen }, level );

   auto errNormInf = err.getMaxDoFMagnitude( level );

   WALBERLA_LOG_DEVEL_VAR( errNormInf );

   WALBERLA_CHECK_LESS( errNormInf, real_c( 10 ) * std::numeric_limits< real_t >::epsilon() );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testMatVec();

   return EXIT_SUCCESS;
}
