/*
* Copyright (c) 2025 Nils Kohl.
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

#include "hyteg/polynomial/PiecewiseLSQPInterpolation.hpp"

#include "core/Environment.h"
#include "core/mpi/MPIIO.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;

template < typename PiecewiseLSQPPolyKDTree >
real_t test( const MeshInfo& meshInfo, uint_t level, uint_t degree, uint_t depth, real_t aabbExtensionFactor )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Degree: " << degree )
   WALBERLA_LOG_INFO_ON_ROOT( "Depth:  " << depth )
   WALBERLA_LOG_INFO_ON_ROOT( "ExtFac: " << aabbExtensionFactor )

   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t dim = PiecewiseLSQPPolyKDTree::NodeType::AABBType::D;
   WALBERLA_CHECK_EQUAL( PiecewiseLSQPPolyKDTree::NodeType::AABBType::D, storage->hasGlobalCells() ? 3 : 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "Dim:    " << dim )

   const auto fileBaseName = "PiecewiseLSQPInterpolationTest_dim" + std::to_string( dim ) + "_deg_" + std::to_string( degree ) +
                             "_depth_" + std::to_string( depth );

   P1Function< real_t > u( "u", storage, level, level );
   P1Function< real_t > uPolyEval( "uPolyEval", storage, level, level );
   P1Function< real_t > err( "err", storage, level, level );

   // Interpolate something oscillatory.
   auto f = []( const Point3D& x ) { return std::sin( 10.0 * x( 0 ) ) + std::cos( 20.0 * x( 1 ) ) + std::sin( 15.0 * x( 2 ) ); };
   u.interpolate( f, level );

   // Compute AABB of domain.
   PointND< real_t, dim > aabbMin;
   aabbMin.fill( std::numeric_limits< real_t >::max() );
   PointND< real_t, dim > aabbMax;
   aabbMax.fill( std::numeric_limits< real_t >::min() );
   for ( auto vertex : meshInfo.getVertices() )
   {
      for ( uint_t d = 0; d < dim; ++d )
      {
         aabbMin[d] = std::min( vertex.second.getCoordinates()( d ), aabbMin[d] );
         aabbMax[d] = std::max( vertex.second.getCoordinates()( d ), aabbMax[d] );
      }
   }
   using AABB = AABB< PiecewiseLSQPPolyKDTree::NodeType::AABBType::D >;
   AABB aabb( aabbMin, aabbMax );

   // Building a KDTree.
   PiecewiseLSQPPolyKDTree kdTree( aabb, depth );

   WALBERLA_LOG_INFO_ON_ROOT( "Collecting points for interpolation and solving LSQP problems." );
   setupPiecewiseLSQPPolyKDTree( kdTree, degree, aabbExtensionFactor, level, u );

   const auto numCoeffs = dim == 2 ? Polynomial2D< MonomialBasis2D >::getNumCoefficients( degree ) :
                                     Polynomial3D< MonomialBasis3D >::getNumCoefficients( degree );
   const auto numLeaves = kdTree.getNumberOfLeaves();
   const auto dofs      = u.getNumberOfGlobalDoFs( level );

   WALBERLA_LOG_INFO_ON_ROOT( "DoFs: " << dofs << "." );
   WALBERLA_LOG_INFO_ON_ROOT( "Built KDTree with " << numLeaves << " leaves." );
   WALBERLA_LOG_INFO_ON_ROOT( "Polynomials with " << numCoeffs << " coefficients." );
   WALBERLA_LOG_INFO_ON_ROOT( "Overall: " << numCoeffs * numLeaves << " coefficients in tree." );
   WALBERLA_LOG_INFO_ON_ROOT( "Compression: " << 100 * ( numCoeffs * numLeaves ) / dofs << " % of DoFs." );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing (de-)serialization." );

   SendBuffer sendBuffer;
   kdTree.serialize( sendBuffer );

   walberla::mpi::writeMPIIO( "../../../output/" + fileBaseName + ".buf", sendBuffer );

   RecvBuffer recvBuffer;

   walberla::mpi::readMPIIO( "../../../output/" + fileBaseName + ".buf", recvBuffer );

   PiecewiseLSQPPolyKDTree kdTreeRead( aabb, 0 );
   kdTreeRead.deserialize( recvBuffer );

   WALBERLA_LOG_INFO_ON_ROOT( "Evaluating polynomial." );
   uPolyEval.interpolate( piecewiseKDTreePolyEvaluator( kdTreeRead ), level );

   WALBERLA_LOG_INFO_ON_ROOT( "Computing error." );
   err.assign( { 1.0, -1.0 }, { u, uPolyEval }, level );

   const auto infErr = err.getMaxDoFMagnitude( level, All );
   const auto maxMag = u.getMaxDoFMagnitude( level, All );
   const auto relErr = infErr / maxMag;
   WALBERLA_LOG_INFO_ON_ROOT( "Max relative error: " << relErr );

   WALBERLA_LOG_INFO_ON_ROOT( "VTK." );
   VTKOutput vtk( "../../../output", fileBaseName, storage );
   vtk.add( u );
   vtk.add( uPolyEval );
   vtk.add( err );
   vtk.write( level );

   for ( uint_t d = 0; d <= kdTreeRead.maxDepth; d++ )
   {
      writeKDTreeAABBsToVTK( kdTreeRead, d, "../../../output", fileBaseName + "_aabbs_depth_" + std::to_string( d ) );
   }

   return relErr;
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   {
      MeshInfo   meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/bfs_12el.msh" ) );
      const auto err      = test< PiecewiseLSQPPolyKDTree2D >( meshInfo, 5, 4, 6, 1.0 );

      WALBERLA_CHECK_LESS( err, 0.044 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      MeshInfo   meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/bfs_12el.msh" ) );
      const auto err      = test< PiecewiseLSQPPolyKDTree2D >( meshInfo, 5, 5, 5, 1.0 );

      if ( std::is_same_v< real_t, double > )
      {
         WALBERLA_CHECK_LESS( err, 0.012 );
      }
      else if ( std::is_same_v< real_t, float > )
      {
         WALBERLA_CHECK_LESS( err, 0.081 );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      MeshInfo   meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/bfs_12el.msh" ) );
      const auto err      = test< PiecewiseLSQPPolyKDTree2D >( meshInfo, 5, 6, 5, 1.0 );

      if ( std::is_same_v< real_t, double > )
      {
         WALBERLA_CHECK_LESS( err, 0.0016 );
      }
      else if ( std::is_same_v< real_t, float > )
      {
         WALBERLA_CHECK_LESS( err, 0.042 );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      MeshInfo   meshInfo = MeshInfo::meshTorus( 5, 5, 1.0, { 0.1, 0.2 } );
      const auto err      = test< PiecewiseLSQPPolyKDTree3D >( meshInfo, 3, 5, 6, 1.0 );

      if ( std::is_same_v< real_t, double > )
      {
         WALBERLA_CHECK_LESS( err, 0.014 );
      }
      else if ( std::is_same_v< real_t, float > )
      {
         WALBERLA_CHECK_LESS( err, 0.086 );
      }
   }
}