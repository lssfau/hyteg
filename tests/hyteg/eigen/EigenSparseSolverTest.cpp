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

#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"

#include "hyteg/eigen/EigenSparseDirectSolver.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

real_t test( const std::string& meshFile, uint_t level )
{
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1ElementwiseLaplaceOperator L( storage, level, level );

   P1Function< real_t > r( "r", storage, level, level );
   P1Function< real_t > f( "f", storage, level, level );
   P1Function< real_t > u( "u", storage, level, level );
   P1Function< real_t > u_exact( "u_exact", storage, level, level );
   P1Function< real_t > err( "err", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   EigenSparseDirectSolver< P1ElementwiseLaplaceOperator > solver( storage, level );
   solver.solve( L, u, f, level );

   // Solving twice to test whether that breaks anything (has happened in the past).
   // This should produce the same result!
   solver.solve( L, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / real_c( u.getNumberOfGlobalDoFs( level ) ) );

   L.apply( u, r, level, Inner, Replace );

   real_t discr_l2_residual = std::sqrt( r.dotGlobal( r, level ) / real_c( u.getNumberOfGlobalDoFs( level ) ) );

   WALBERLA_CHECK_LESS( discr_l2_residual, 10 * std::numeric_limits< real_t >::epsilon() );

   return discr_l2_err;
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   real_t lastError = 0;

   for ( uint_t level = 2; level < 6; level++ )
   {
      auto err = hyteg::test( hyteg::prependHyTeGMeshDir( "2D/quad_4el.msh" ), level );
      WALBERLA_LOG_INFO_ON_ROOT( "level: " << level << " | err: " << err );

      if ( level > 2 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "rate: " << err / lastError );
         WALBERLA_CHECK_LESS( err / lastError, .275 );
      }

      lastError = err;
   }
}
