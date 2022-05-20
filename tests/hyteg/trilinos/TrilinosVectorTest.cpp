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

#include "hyteg/trilinos/TrilinosVector.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/trilinos/TrilinosSparseMatrix.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( Tpetra::version() )

   const uint_t level = 3;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, 1, 1 );
   auto     setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   P2P1TaylorHoodFunction< real_t >   xTrilinos( "xTrilinos", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   xOriginal( "xOriginal", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   error( "error", storage, level, level );
   P2P1TaylorHoodFunction< idx_t >    numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   auto f = []( const Point3D& p ) -> real_t { return std::sin( p[0] ) + 0.5 * p[1]; };

   xOriginal.uvw().interpolate( {f, f}, level, All );
   xOriginal.p().interpolate( f, level, All );

   trilinos::TrilinosVector< P2P1TaylorHoodFunction > vector( storage, level );
   vector.fillFromFunction( xOriginal, numerator );
   vector.writeToFunction( xTrilinos, numerator );

   error.assign( {1.0, -1.0}, {xTrilinos, xOriginal}, level );
   const auto maxMagnitudeU = error.uvw()[0].getMaxMagnitude( level );
   const auto maxMagnitudeV = error.uvw()[1].getMaxMagnitude( level );
   const auto maxMagnitudeP = error.p().getMaxMagnitude( level );

   WALBERLA_CHECK_LESS( maxMagnitudeU, 1e-14 );
   WALBERLA_CHECK_LESS( maxMagnitudeV, 1e-14 );
   WALBERLA_CHECK_LESS( maxMagnitudeP, 1e-14 );

   return EXIT_SUCCESS;
}
