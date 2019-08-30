
#include "hyteg/mesh/MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/math/Constants.h"

#include <array>
#include <vector>

using walberla::math::pi;

namespace hyteg {

MeshInfo MeshInfo::meshAnnulus( const real_t rhoMin, const real_t rhoMax,
                                const real_t phiLeft, const real_t phiRight,
                                const meshFlavour flavour, uint_t nTan, uint_t nRad )
{

  WALBERLA_ASSERT_LESS( rhoMin, rhoMax );
  WALBERLA_ASSERT_LESS( phiLeft, phiRight );
  WALBERLA_ASSERT_GREATER( rhoMin, 1e-8 );
  WALBERLA_ASSERT_GREATER( nTan, 0 );
  WALBERLA_ASSERT_GREATER( nRad, 0 );

  // mesh partial annulus in polar coordinates
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { rhoMin, phiLeft  } ),
                                               Point2D( { rhoMax, phiRight } ),
                                               flavour, nRad, nTan );

  // map vertex coordinates to cartesian domain
  Point3D node;
  node[2] = 0.0;
  uint_t boundaryFlag;
  real_t rho, phi;
  for ( size_t id = 0; id < meshInfo.vertices_.size(); ++id )
    {
      node = meshInfo.vertices_[id].getCoordinates();
      boundaryFlag = meshInfo.vertices_[id].getBoundaryFlag();
      rho = node[0];
      phi = node[1];
      node[0] = rho*std::cos(phi);
      node[1] = rho*std::sin(phi);
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), boundaryFlag );
    }

  return meshInfo;

}


MeshInfo MeshInfo::meshAnnulus( const real_t rhoMin, const real_t rhoMax, uint_t nTan, uint_t nRad )
{

  WALBERLA_ASSERT_LESS( rhoMin, rhoMax );
  WALBERLA_ASSERT_GREATER( rhoMin, 1e-8 );
  WALBERLA_ASSERT_GREATER( nTan, 0 );
  WALBERLA_ASSERT_GREATER( nRad, 0 );

  MeshInfo meshInfo;

  // -------------------
  //  generate vertices
  // -------------------
  IDType id = 0;
  std::map< std::array<uint_t,2>, IDType > indices2id;

  std::array<real_t,3> node;
  node[2] = (real_t)0.0;

  real_t deltaPhi = 2.0 * pi / (real_t)nTan;
  real_t deltaRho = (rhoMax - rhoMin) / (real_t)nRad;

  for ( uint_t i = 0; i < nTan; ++i )
    {
      // inner boundary
      node[0] = rhoMin * std::cos( (real_t)i * deltaPhi );
      node[1] = rhoMin * std::sin( (real_t)i * deltaPhi );
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), 1 );
      indices2id.insert( { {i,0}, id } );
      ++id;

      // interior nodes
      for ( uint_t j = 1; j < nRad; ++j )
        {
          node[0] = (rhoMin + (real_t)j * deltaRho) * std::cos( (real_t)i * deltaPhi );
          node[1] = (rhoMin + (real_t)j * deltaRho) * std::sin( (real_t)i * deltaPhi );
          meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), 0 );
          indices2id.insert( { {i,j}, id } );
          ++id;
        }

      // outer boundary
      node[0] = rhoMax * std::cos( (real_t)i * deltaPhi );
      node[1] = rhoMax * std::sin( (real_t)i * deltaPhi );
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), 1 );
      indices2id.insert( { {i,nRad}, id } );
      ++id;
    }

  // map node indices to linear id
  auto getIDX = [ indices2id ] ( uint_t ii, uint_t jj ) -> IDType
    {
      std::map< std::array<uint_t,2>, IDType >::const_iterator found = indices2id.find( { ii, jj } );
      WALBERLA_CHECK( found != indices2id.end(), "Could not map tupled to index!" );
      return found->second;
    };

  // --------------------
  //  generate triangles
  // --------------------
  real_t midPhi = 0.0;
  real_t midRho = 0.0;
  for ( uint_t i = 0; i < nTan; ++i )
    {
      for ( uint_t j = 0; j < nRad; ++j )
        {

          // add new central vertex
          midPhi = ( (real_t)i + 0.5 ) * deltaPhi;
          midRho = rhoMin + ( (real_t)j + 0.5 ) * deltaRho;
          node[0] = midRho * std::cos( midPhi );
          node[1] = midRho * std::sin( midPhi );
          meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), 0 );

          // add four sub-triangles of cell
          meshInfo.addFace( Face( { getIDX( ( i )     ,  j  ), id, getIDX( (i+1)%nTan,  j  ) }, 0 ) );
          meshInfo.addFace( Face( { getIDX( (i+1)%nTan,  j  ), id, getIDX( (i+1)%nTan, j+1 ) }, 0 ) );
          meshInfo.addFace( Face( { getIDX( (i+1)%nTan, j+1 ), id, getIDX( ( i )     , j+1 ) }, 0 ) );
          meshInfo.addFace( Face( { getIDX( ( i )     , j+1 ), id, getIDX( ( i )     ,  j  ) }, 0 ) );
          ++id;
        }
    }

  // generate edges from faces
  meshInfo.deriveEdgesForFullAnnulus( rhoMin + 0.1 * deltaRho, rhoMax - 0.1 * deltaRho );

  // done
  return meshInfo;
}


void MeshInfo:: deriveEdgesForFullAnnulus( real_t minTol, real_t maxTol )
{

  MeshInfo::FaceContainer faces = this->getFaces();
  MeshInfo::VertexContainer verts = this->getVertices();

  uint_t edgeBoundaryFlag = 0;
  Point3D node1, node2;
  real_t radius1, radius2;

  for ( const auto & it : faces )
    {
      // extract the three nodes of the face
      std::vector<IDType> fNode = it.second.getVertices();

      // determine their position w.r.t. the boundary
      std::vector< uint_t > meshBoundaryFlags( 3 );
      meshBoundaryFlags[0] = verts.find( fNode[0] )->second.getBoundaryFlag();
      meshBoundaryFlags[1] = verts.find( fNode[1] )->second.getBoundaryFlag();
      meshBoundaryFlags[2] = verts.find( fNode[2] )->second.getBoundaryFlag();

      // set the three edges of triangle, edge is on boundary, if both
      // its vertices are
      edgeBoundaryFlag = 0;
      if( meshBoundaryFlags[0] == 1 && meshBoundaryFlags[1] == 1 )
        {
          radius1 = verts.find( fNode[0] )->second.getCoordinates().norm();
          radius2 = verts.find( fNode[1] )->second.getCoordinates().norm();
          if( (radius1 < minTol && radius2 < minTol) ||
              (radius1 > maxTol && radius2 < maxTol) )
            {
              edgeBoundaryFlag = 1;
            }
        }
      this->addEdge( Edge( { fNode[0], fNode[1] }, edgeBoundaryFlag ) );

      edgeBoundaryFlag = 0;
      if( meshBoundaryFlags[0] == 1 && meshBoundaryFlags[2] == 1 )
        {
          radius1 = verts.find( fNode[0] )->second.getCoordinates().norm();
          radius2 = verts.find( fNode[2] )->second.getCoordinates().norm();
          if( (radius1 < minTol && radius2 < minTol) ||
              (radius1 > maxTol && radius2 < maxTol) )
            {
              edgeBoundaryFlag = 1;
            }
        }
      this->addEdge( Edge( { fNode[0], fNode[2] }, edgeBoundaryFlag ) );

      edgeBoundaryFlag = 0;
      if( meshBoundaryFlags[1] == 1 && meshBoundaryFlags[2] == 1 )
        {
          radius1 = verts.find( fNode[1] )->second.getCoordinates().norm();
          radius2 = verts.find( fNode[2] )->second.getCoordinates().norm();
          if( (radius1 < minTol && radius2 < minTol) ||
              (radius1 > maxTol && radius2 < maxTol) )
            {
              edgeBoundaryFlag = 1;
            }
        }
      this->addEdge( Edge( { fNode[1], fNode[2] }, edgeBoundaryFlag ) );
    }
}

} // namespace hyteg
