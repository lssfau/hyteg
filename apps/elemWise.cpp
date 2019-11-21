/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
// #include "hyteg/solvers/CGSolver.hpp"
// #include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
// #include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
// #include "hyteg/solvers/GeometricMultigridSolver.hpp"
// #include "hyteg/solvers/GaussSeidelSmoother.hpp"

#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::pi;

using namespace hyteg;


void localMatrixVectorMultiply( Face& face, uint_t level, uint_t xIdx, uint_t yIdx,
                                const P2Elements::P2Element& element,
                                real_t* srcVertexData, real_t* srcEdgeData,
                                real_t* dstVertexData, real_t* dstEdgeData,
                                UpdateType update ) {

  WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

  P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > p2mass;

  Matrix6r elMat;
  Point6D elVecOld, elVecNew;
  indexing::Index nodeIdx;
  indexing::IndexIncrement offset;
  Point3D v0, v1, v2;
  std::array<uint_t,6> dofDataIdx;

  // determine element's vertices
  nodeIdx = indexing::Index( xIdx, yIdx, 0 );
  v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
  offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
  v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
  offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
  v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

  // assemble local element matrix
  p2mass.integrateAll( {v0, v1, v2}, elMat );

  // assemble local element vector (note the tweaked ordering to go along with FEniCS indexing)
  dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
  dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
  dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

  dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
  dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
  dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

  elVecOld[0] = srcVertexData[ dofDataIdx[0] ];
  elVecOld[1] = srcVertexData[ dofDataIdx[1] ];
  elVecOld[2] = srcVertexData[ dofDataIdx[2] ];

  elVecOld[3] = srcEdgeData[ dofDataIdx[3] ];
  elVecOld[4] = srcEdgeData[ dofDataIdx[4] ];
  elVecOld[5] = srcEdgeData[ dofDataIdx[5] ];

  // apply matrix (operator locally)
  elVecNew = elMat.mul( elVecOld );

  // redistribute result from "local" to "global vector"
  if ( update == Replace ) {
    dstVertexData[ dofDataIdx[0] ] = elVecNew[0];
    dstVertexData[ dofDataIdx[1] ] = elVecNew[1];
    dstVertexData[ dofDataIdx[2] ] = elVecNew[2];

    dstEdgeData[ dofDataIdx[3] ] = elVecNew[3];
    dstEdgeData[ dofDataIdx[4] ] = elVecNew[4];
    dstEdgeData[ dofDataIdx[5] ] = elVecNew[5];
  }
  else {
    dstVertexData[ dofDataIdx[0] ] += elVecNew[0];
    dstVertexData[ dofDataIdx[1] ] += elVecNew[1];
    dstVertexData[ dofDataIdx[2] ] += elVecNew[2];

    dstEdgeData[ dofDataIdx[3] ] += elVecNew[3];
    dstEdgeData[ dofDataIdx[4] ] += elVecNew[4];
    dstEdgeData[ dofDataIdx[5] ] += elVecNew[5];
  }
}


int main(int argc, char* argv[])
{

  // Setup enviroment
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  // MeshInfo meshInfo = MeshInfo( MeshInfo::fromGmshFile( "../data/meshes/tri_2el.msh" ) );
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {-1.0, -1.0} ), Point2D( {1.0, 1.0} ), MeshInfo::CRISS, 2, 2 );

  // Prepare storage
  SetupPrimitiveStorage setupStorage( meshInfo,
                                      uint_c( walberla::mpi::MPIManager::instance()->numProcesses()));
  hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  uint_t level = 2;
  UpdateType update = Add;
  // UpdateType update = Replace;
  
  // We need two simple functions
  P2Function< real_t > out( "out", storage, level, level );
  P2Function< real_t > in( "vecOfOnes", storage, level, level );
  in.interpolate( real_c(1), level );

  P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > p2mass;

  // We need to initialise the "global output vector" to zero, as the
  // results is computed additively from the local element contributions
  out.interpolate( real_c(0), level );

  // Feels like deja vu
  communication::syncP2FunctionBetweenPrimitives( out, level );
  communication::syncP2FunctionBetweenPrimitives( in, level );

  for( auto& it : storage->getFaces() ) {
    Face& face = *it.second;
    
    Point3D x0( face.coords[0] );
    Point3D x1( face.coords[1] );
    Point3D x2( face.coords[2] );

    WALBERLA_LOG_INFO_ON_ROOT( "Vertex #0:" << x0 );
    WALBERLA_LOG_INFO_ON_ROOT( "Vertex #1:" << x1 );
    WALBERLA_LOG_INFO_ON_ROOT( "Vertex #2:" << x2 );

    WALBERLA_LOG_INFO_ON_ROOT( "For level = " << level << " we have " << levelinfo::num_microfaces_per_face(level)
                               << " microfaces" );

    uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
    uint_t inner_rowsize = rowsize;
    uint_t xIdx, yIdx;
    Point3D v0, v1, v2;
    indexing::Index nodeIdx;
    indexing::IndexIncrement offset;

    // get hold of the actual numerical data in the two functions
    PrimitiveDataID<FunctionMemory< real_t >, Face> dstVertexDoFIdx = out.getVertexDoFFunction().getFaceDataID();
    PrimitiveDataID<FunctionMemory< real_t >, Face> srcVertexDoFIdx = in.getVertexDoFFunction().getFaceDataID();

    PrimitiveDataID<FunctionMemory< real_t >, Face> dstEdgeDoFIdx = out.getEdgeDoFFunction().getFaceDataID();
    PrimitiveDataID<FunctionMemory< real_t >, Face> srcEdgeDoFIdx = in.getEdgeDoFFunction().getFaceDataID();

    real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
    real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

    real_t* srcEdgeData = face.getData( srcEdgeDoFIdx )->getPointer( level );
    real_t* dstEdgeData = face.getData( dstEdgeDoFIdx )->getPointer( level );

    for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx ) {

      WALBERLA_LOG_INFO_ON_ROOT( "yIdx = " << yIdx );

      // loop over vertices in row with two associated triangles
      for( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx ) {

        // element one associated with vertex
        localMatrixVectorMultiply( face, level, xIdx, yIdx, P2Elements::P2Face::elementN,
                                   srcVertexData, srcEdgeData, dstVertexData, dstEdgeData, update );

        // element two associated with vertex
        localMatrixVectorMultiply( face, level, xIdx, yIdx, P2Elements::P2Face::elementNW,
                                   srcVertexData, srcEdgeData, dstVertexData, dstEdgeData, update );
      }
      --inner_rowsize;

      // final vertex in row has only one associated triangle
      localMatrixVectorMultiply( face, level, xIdx, yIdx, P2Elements::P2Face::elementNW,
                                 srcVertexData, srcEdgeData, dstVertexData, dstEdgeData, update );
    }

    // top north-west element not treated, yet
    WALBERLA_LOG_INFO_ON_ROOT( "yIdx = " << yIdx );
    localMatrixVectorMultiply( face, level, 1, yIdx, P2Elements::P2Face::elementNW,
                               srcVertexData, srcEdgeData, dstVertexData, dstEdgeData, update );
  }

  // Push result to lower-dimensional primitives
  out.getVertexDoFFunction().communicateAdditively< Face, Edge >( level );
  out.getVertexDoFFunction().communicateAdditively< Face, Vertex >( level );
  out.getEdgeDoFFunction().communicateAdditively< Face, Edge >( level );

  // Check consistency with P2ElementwiseOperator
  P2ElementwiseMassOperator p2MassOp( storage, level, level, true );
  P2Function< real_t > elOpOut( "P2ElementwiseOperator result", storage, level, level, DoFType::None ); // None because apply( All )!!
  // p2MassOp.apply( in, elOpOut, level, Inner );
  p2MassOp.apply( in, elOpOut, level, All );
  P2Function< real_t > diff2( "difference 2", storage, level, level );
  diff2.assign( {real_c(1), real_c(-1)}, {elOpOut, out}, level );

  // Consistency check with stencil-based operator
  P2ConstantMassOperator massOp( storage, level, level );
  P2Function< real_t > diff( "difference", storage, level, level );
  P2Function< real_t > aux( "out (stencil-based)", storage, level, level );
  massOp.apply( in, aux, level, Inner );
  diff.assign( {real_c(1), real_c(-1)}, {aux, out}, level );

  // ---------------------------------------------------------
  P2ElementwiseLaplaceOperator laplace( storage, level, level, true );
  P2ConstantLaplaceOperator laplaceConst( storage, level, level );

  P2Function< real_t > uOld( "u0", storage, level, level );
  P2Function< real_t > uNewEl( "u smoothed (elem)", storage, level, level );
  P2Function< real_t > uNewSt( "u smoothed (stencil)", storage, level, level );
  P2Function< real_t > rhs( "right-hand side", storage, level, level );
  rhs.interpolate( real_c(0), level );
  std::function<real_t(const Point3D&)> randFunc = []( const Point3D & ) { return walberla::math::realRandom<real_t>(); };
  uOld.interpolate( randFunc, level, Inner );

  // uNewSt.assign( {1.0}, {uOld}, level, All );
  // laplaceConst.smooth_gs( uNewSt, rhs, level, Inner );

  laplaceConst.smooth_jac( uNewSt, rhs, uOld, 1.0, level, Inner );
  laplace.smooth_jac( uNewEl, rhs, uOld, 1.0, level, Inner );
  // communication::syncP2FunctionBetweenPrimitives( uNewSt, level );

  // --------------------------
  //  Diagonal Stencil Entries
  // --------------------------

  P2Function< real_t > diagMass( "mass operator - diagonal entries", storage, level, level );
  P2Function< real_t > diagDiff( "diffusion operator - diagonal entries", storage, level, level );
  diagMass.assign( {1.0}, {*(p2MassOp.diagonalValues_)}, level, All );
  diagDiff.assign( {1.0}, {*(laplace.diagonalValues_)}, level, All );

  // Analytical values for mass operator
  std::function< real_t( const Point3D& ) > analyticFunc = [level]( const Point3D& x ) {

      // determine width of mesh (fixed for current mesh)
      real_t h = 0.25;

      walberla::int32_t xType = (int) std::round( x[0] / h * 2.0 );
      walberla::int32_t yType = (int) std::round( x[1] / h * 2.0 );

      real_t val = real_c( 0 );

      if ( xType % 2 == 0 && yType % 2 == 0 )
      {
        val = real_c( 0.1 ) * h * h; // vertex dof
      }
      else
      {
        val = real_c( 8.0/45.0 ) * h * h; // edge dof
      }

      return val;
  };
  P2Function< real_t > massDiagStencil( "diagonal values mass operator (analyt)", storage, level, level );
  P2Function< real_t > massDiagStencilError( "error in diagonal stencil values mass opertor", storage, level, level );
  massDiagStencil.interpolate( analyticFunc, level, All );
  massDiagStencilError.assign( {1.0,-1.0}, {massDiagStencil,diagMass}, level, All );

  // output data for visualisation
  hyteg::VTKOutput vtkOutput( "../output", "elementwise", storage );
  vtkOutput.add( in   );
  vtkOutput.add( out  );
  vtkOutput.add( aux  );
  vtkOutput.add( diff );
  vtkOutput.add( diff2 );
  vtkOutput.add( elOpOut );
  vtkOutput.add( uOld );
  vtkOutput.add( uNewEl );
  vtkOutput.add( uNewSt );
  vtkOutput.add( diagMass );
  vtkOutput.add( diagDiff );
  vtkOutput.add( massDiagStencil );
  vtkOutput.add( massDiagStencilError );
  vtkOutput.write( level );

  // WALBERLA_LOG_INFO_ON_ROOT( "num_microvertices_per_face( 1 ) = " << levelinfo::num_microvertices_per_face( 1 ) );
  // WALBERLA_LOG_INFO_ON_ROOT( "num_microvertices_per_face( 2 ) = " << levelinfo::num_microvertices_per_face( 2 ) );
  // WALBERLA_LOG_INFO_ON_ROOT( "num_microvertices_per_face( 3 ) = " << levelinfo::num_microvertices_per_face( 3 ) );

  // WALBERLA_LOG_INFO_ON_ROOT( "num_microedges_per_face( 1 ) = " << levelinfo::num_microedges_per_face( 1 ) );
  // WALBERLA_LOG_INFO_ON_ROOT( "num_microedges_per_face( 2 ) = " << levelinfo::num_microedges_per_face( 2 ) );
  // WALBERLA_LOG_INFO_ON_ROOT( "num_microedges_per_face( 3 ) = " << levelinfo::num_microedges_per_face( 3 ) );

  P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > massForm;
  Matrix6r elMat;
  // massForm.integrateAll( { Point3D( {1.0, 1.0, 0.0} ), Point3D( {2.0, 1.0, 0.0} ), Point3D( {1.0, 2.0, 0.0} ) }, elMat );
  massForm.integrateAll( { Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 0.0, 0.0} ), Point3D( {0.0, 1.0, 0.0} ) }, elMat );
  // massForm.integrateAll( { Point3D( {0.0, 0.0, 0.0} ), Point3D( {0.25, 0.0, 0.0} ), Point3D( {0.0, 0.25, 0.0} ) }, elMat );
  // massForm.integrateAll( { Point3D( {0.0, 1.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {1.0, 0.0, 0.0} ) }, elMat );
  elMat *= 360.0;
  WALBERLA_LOG_INFO_ON_ROOT( "Local Mass Matrix for Unit Triangle (scaled by 360):\n" << elMat );

  auto eStencil = massOp.getEdgeToEdgeOpr();
  auto vStencil = massOp.getVertexToVertexOpr();

  for ( auto& it : storage->getFaces() )
    {
      Face& face = *it.second;

      WALBERLA_LOG_INFO_ON_ROOT( "FACE #" << face.getID() );

      real_t* vS = face.getData( vStencil.getFaceStencilID() )->getPointer( level );
      WALBERLA_LOG_INFO_ON_ROOT( "geometric node = " << vS[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] );

      real_t* eS = face.getData( eStencil.getFaceStencilID() )->getPointer( level );
      WALBERLA_LOG_INFO_ON_ROOT( "edge node (horz) = " << eS[ edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C ) ] );
      WALBERLA_LOG_INFO_ON_ROOT( "edge node (vert) = " << eS[ edgedof::stencilIndexFromVerticalEdge( stencilDirection::EDGE_VE_C ) ] );
      WALBERLA_LOG_INFO_ON_ROOT( "edge node (diag) = " << eS[ edgedof::stencilIndexFromDiagonalEdge( stencilDirection::EDGE_DI_C ) ] << "\n" );
    }

  for ( auto& it : storage->getEdges() )
    {
      Edge& edge = *it.second;

      if ( edge.getMeshBoundaryFlag() == 0 )
        {
          WALBERLA_LOG_INFO_ON_ROOT( "EDGE #" << edge.getID() );

          real_t* vS = edge.getData( vStencil.getEdgeStencilID() )->getPointer( level );
          WALBERLA_LOG_INFO_ON_ROOT( "geometric node = " << vS[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] );

          real_t* eS = edge.getData( eStencil.getEdgeStencilID() )->getPointer( level );
          WALBERLA_LOG_INFO_ON_ROOT( "edge node (horz) = " << eS[ edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C ) ] << "\n" );
        }
    }

  for ( auto& it : storage->getVertices() )
    {
      Vertex& vertex = *it.second;

      if ( vertex.getMeshBoundaryFlag() == 0 )
        {
          WALBERLA_LOG_INFO_ON_ROOT( "VERTEX #" << vertex.getID() );

          real_t* vS = vertex.getData( vStencil.getVertexStencilID() )->getPointer( level );
          // WALBERLA_LOG_INFO_ON_ROOT( "geometric node = " << vS[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] );
          WALBERLA_LOG_INFO_ON_ROOT( "geometric node = " << vS[0] );
        }
    }

  return 0;
}
