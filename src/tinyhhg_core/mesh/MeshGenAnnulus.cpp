
#include "tinyhhg_core/mesh/MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

#include <array>
#include <vector>

using walberla::math::PI;

namespace hhg {

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
  DoFType dofType;
  real_t rho, phi;
  for ( size_t id = 0; id < meshInfo.vertices_.size(); ++id )
    {
      node = meshInfo.vertices_[id].getCoordinates();
      dofType = meshInfo.vertices_[id].getDoFType();
      rho = node[0];
      phi = node[1];
      node[0] = rho*std::cos(phi);
      node[1] = rho*std::sin(phi);
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), dofType );
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

  real_t deltaPhi = 2.0 * PI / (real_t)nTan;
  real_t deltaRho = (rhoMax - rhoMin) / (real_t)nRad;

  for ( uint_t i = 0; i < nTan; ++i )
    {
      // inner boundary
      node[0] = rhoMin * std::cos( (real_t)i * deltaPhi );
      node[1] = rhoMin * std::sin( (real_t)i * deltaPhi );
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), DirichletBoundary );
      indices2id.insert( { {i,0}, id } );
      ++id;

      // interior nodes
      for ( uint_t j = 1; j < nRad; ++j )
        {
          node[0] = (rhoMin + (real_t)j * deltaRho) * std::cos( (real_t)i * deltaPhi );
          node[1] = (rhoMin + (real_t)j * deltaRho) * std::sin( (real_t)i * deltaPhi );
          meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), Inner );
          indices2id.insert( { {i,j}, id } );
          ++id;
        }

      // outer boundary
      node[0] = rhoMax * std::cos( (real_t)i * deltaPhi );
      node[1] = rhoMax * std::sin( (real_t)i * deltaPhi );
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), DirichletBoundary );
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
          meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), Inner );

          // add four sub-triangles of cell
          meshInfo.addFace( Face( { getIDX( ( i )     ,  j  ), id, getIDX( (i+1)%nTan,  j  ) }, Inner ) );
          meshInfo.addFace( Face( { getIDX( (i+1)%nTan,  j  ), id, getIDX( (i+1)%nTan, j+1 ) }, Inner ) );
          meshInfo.addFace( Face( { getIDX( (i+1)%nTan, j+1 ), id, getIDX( ( i )     , j+1 ) }, Inner ) );
          meshInfo.addFace( Face( { getIDX( ( i )     , j+1 ), id, getIDX( ( i )     ,  j  ) }, Inner ) );
          ++id;
        }
    }

  // generate edges from faces
  meshInfo.deriveEdges();

  // done
  return meshInfo;
}

} // namespace hhg
