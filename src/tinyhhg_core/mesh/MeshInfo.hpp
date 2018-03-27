
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"

#include <array>
#include <map>
#include <set>
#include <vector>

namespace hhg {

using walberla::uint_t;
using walberla::real_t;

/// \brief Contains information about a mesh
/// \author Daniel Drzisga (drzisga@ma.tum.de) \n
///         Nils Kohl (nils.kohl@fau.de)
///         Marcus Mohr (marcus.mohr@lmu.de)
///
/// The \ref MeshInfo class is used to store information about meshes that are
/// either constructed by reading mesh files or generated internally for some
/// standard geometries. It provides static methods to construct meshes by
/// - importing a file in Gmsh format (see method \ref fromGmshFile)
/// - internal mesh generation for rectangles
/// - internal mesh generation for full or partial annuli
///
/// \ref MeshInfo instance can then be used to create domains via \ref SetupPrimitiveStorage.
/// <b>Details of inline mesh generators for rectangles:</b>
///
/// The inline mesh generator expects to input arguments nx and ny for the discretisation.
/// It supports generation of four different flavours of meshes.
/// - For CRISS, CROSS and CRISSCROSS the rectangle is split into a regular mesh of nx by ny
///   cells, which are then subdivived into two (CRISS and CROSS) or four (CRISSCROSS)
///   triangles.
/// - For DIAMOND we use the same number and position of vertices as for CRISSCROSS, but
///   these are connected in a layered fashion from the outside in.
///
/// | flavour    |       #vertices       | #triangles |
/// |:-----------|:---------------------:|:----------:|
/// | CRISS      | (nx+1)*(ny+1)         | 2*nx*ny    |
/// | CROSS      | (nx+1)*(ny+1)         | 2*nx*ny    |
/// | CRISSCROSS | (nx+1)*(ny+1) + nx*ny | 4*nx*ny    |
/// | DIAMOND    | (nx+1)*(ny+1) + nx*ny | 4*nx*ny    |
///
/*! \htmlonly
  <center>
  <table>
  <tr>
  <td colspan="4" align="center">Sample mesh generated for a rectangle using (nx=3, ny=2)</td>
  </tr>
  <tr>
  <td><img src="../ExtraPics/Mesh_RectangleCriss.png" width="100%"/></td>
  <td><img src="../ExtraPics/Mesh_RectangleCross.png" width="100%"/></td>
  <td><img src="../ExtraPics/Mesh_RectangleCrissCross.png" width="100%"/></td>
  <td><img src="../ExtraPics/Mesh_RectangleDiamond.png" width="100%"/></td>
  </tr>
  <tr>
  <td align="center">CRISS</td>
  <td align="center">CROSS</td>
  <td align="center">CRISSCROSS</td>
  <td align="center">DIAMOND</td>
  </tr>
  </table>
  </center>
  \endhtmlonly
*/
/// <b>Details of inline mesh generators for annuli:</b>
///
/// Meshing of a partial annulus is (conceptually) handled by meshing the correspondig
/// rectangle in cartesian coordinates. In case of a partial annulus this is given by
/// lower left vertex (rhoMin, phiMin) and upper right vertex (rhoMax, phiMax).
/// - For a partial annulus the same four flavours as for rectangles can be specified.
/// - A full annulus is meshed using a CRISSCROSS pattern resulting in
///
/// |              |     #vertices    | #triangles  |
/// |:-------------|:----------------:|:-----------:|
/// | full annulus |  nTan*(2*nRad+1) | 4*nTan*nRad |
///
/*! \htmlonly
  <center>
  <table>
  <tr>
  <td align="center"><img src="../ExtraPics/Mesh_AnnulusPartial.png" width="50%"/></td>
  <td align="center"><img src="../ExtraPics/Mesh_AnnulusFull.png" width="50%"/></td>
  </tr>
  <tr>
  <td align="center">partial annulus (nTan=4, nRad=2)</td>
  <td align="center">full annulus (nTan=15, nRad=2)</td>
  </tr>
  </table>
  </center>
  \endhtmlonly
*/
///
/// \note The inline mesh generators currently set all vertex and edge primitives on the domain boundary to
///       **DoFType DirichletBoundary** and those inside the domain, and of course all face primitives, to
///       **DoFType Inner**.
class MeshInfo
{
public:

  typedef uint_t IDType;

  typedef enum { CRISS, CROSS, CRISSCROSS, DIAMOND } meshFlavour;

  class Vertex
  {
  public:
    Vertex() : id_( 0 ), coordinates_( Point3D() ), dofType_( Inner ) {};
    Vertex( const IDType & id, const Point3D & coordinates, const DoFType & dofType ) :
      id_( id ), coordinates_( coordinates ), dofType_( dofType )
    {}

    const IDType  getID()          const { return id_; }
    const Point3D getCoordinates() const { return coordinates_; }
    const DoFType getDoFType()     const { return dofType_; }

  private:
    IDType  id_;
    Point3D coordinates_;
    DoFType dofType_;
  };

  class Edge
  {
  public:
    Edge() : dofType_( Inner ) {};
    Edge( const std::array< IDType, 2 > & vertices, const DoFType & dofType ) :
      vertices_( vertices ), dofType_( dofType )
    {}

    const std::array< IDType, 2 > getVertices() const { return vertices_; }
    const DoFType                 getDoFType()  const { return dofType_; }

  private:
    std::array< IDType, 2 > vertices_;
    DoFType                 dofType_;
  };

  class Face
  {
  public:
    Face() : dofType_( Inner ) {};
    Face( const std::vector< IDType > & vertices, const DoFType & dofType ) :
      vertices_( vertices ), dofType_( dofType )
    {}

    const std::vector< IDType > getVertices() const { return vertices_; }
    const DoFType               getDoFType()  const { return dofType_; }

  private:
    std::vector< IDType > vertices_;
    DoFType               dofType_;
  };

  class Cell
  {
  public:
    Cell() : dofType_( Inner ) {};
    Cell( const std::vector< IDType > & vertices, const DoFType & dofType ) :
      vertices_( vertices ), dofType_( dofType )
    {}

    const std::vector< IDType > getVertices() const { return vertices_; }
    const DoFType               getDoFType()  const { return dofType_; }

  private:
    std::vector< IDType > vertices_;
    DoFType               dofType_;
  };


  typedef std::map< IDType,                  Vertex > VertexContainer;
  typedef std::map< std::array< IDType, 2 >, Edge   > EdgeContainer;
  typedef std::map< std::vector< IDType >,   Face   > FaceContainer;
  typedef std::map< std::vector< IDType >,   Cell   > CellContainer;

  /// Construct empty MeshInfo (for testing)
  static MeshInfo emptyMeshInfo() { return MeshInfo(); }

  /// Construct a MeshInfo object from a file in Gmsh format
  static MeshInfo fromGmshFile( const std::string & meshFileName );

  /// Construct a MeshInfo object for a rectangular domain

  /// \param lowerLeft    coordinates of lower left corner of rectangle
  /// \param upperRigth   coordinates of upper right corner of rectangle
  /// \param flavour      meshing strategy (CRISS, CROSS, CRISSCROSS or DIAMOND)
  /// \param nx           (nx+1) gives the number of vertices along the top and bottom edge of the rectangle
  /// \param ny           (ny+1) gives the number of vertices along the left and right edge of the rectangle
  static MeshInfo meshRectangle( const Point2D lowerLeft, const Point2D upperRight, const meshFlavour flavour,
                                 uint_t nx, uint_t ny );

  /// Construct a MeshInfo object for a partial annulus

  /// \param rhoMin       radius of inner circle of partial annulus
  /// \param rhoMax       radius of outer circle of partial annulus
  /// \param phiMin       smaller angle of radial boundary in polar coordinates
  /// \param phiMax       larger angle of radial boundary in polar coordinates
  /// \param flavour      meshing strategy (CRISS, CROSS, CRISSCROSS or DIAMOND)
  /// \param nTan         number of tangential subdivisions (along inner and outer circle)
  /// \param nRad         number of radial subdivisions (along left and right radial boundary)
  static MeshInfo meshAnnulus( const real_t rhoMin, const real_t rhoMax, const real_t phiLeft,
                               const real_t phiRight, const meshFlavour flavour, uint_t nTan, uint_t nRad );

  /// Construct a MeshInfo object for a full annulus
  static MeshInfo meshAnnulus( const real_t rmin, const real_t rmax, uint_t nTan, uint_t nRad );

  /// Returns vertices of the mesh
  const VertexContainer & getVertices() const { return vertices_; };

  /// Returns edges of the mesh
  const EdgeContainer & getEdges() const { return edges_; };

  /// Returns faces of the mesh
  const FaceContainer & getFaces() const { return faces_; }

  /// Returns cells of the mesh
  const CellContainer & getCells() const { return cells_; }

private:

  MeshInfo() {};

  /// Adds edge in ascending index order and performs checks
  void addEdge( const Edge & edge );

  /// Adds face in ascending index order and performs checks
  void addFace( const Face & face );

  /// Construct a MeshInfo for a rectangular domain using diamond approach
  static MeshInfo meshRectangleDiamond( const Point2D lowerLeft, const Point2D upperRight, uint_t nx, uint_t ny );

  /// Derive information on edges from vertices and faces (for rectangles)

  /// This method is used in the 2D inline mesh generators for rectangles. The latter
  /// provide only information on vertices and faces. The information on the edges is
  /// derived from the faces, as each face has three edges.
  /// \param lowerLeft   lower left vertex of rectangle
  /// \param upperRight  upper right vertex of rectangle
  /// \param tol         parameter used to determine, whether an edge is part of the boundary, or not
  void deriveEdgesForRectangles( const Point2D lowerLeft, const Point2D upperRight, real_t tol );

  /// Derive information on edges from vertices and faces (for full annulus)

  /// This method is used in the 2D inline mesh generator for the full annulus. The latter
  /// provides only information on vertices and faces. The information on the edges is then
  /// derived from the faces, as each face has three edges.
  /// \param minTol  parameter used to determine, whether an edge is part of the inner boundary, or not
  /// \param maxTol  parameter used to determine, whether an edge is part of the outer boundary, or not
  void deriveEdgesForFullAnnulus( real_t minTol, real_t maxTol );

  VertexContainer vertices_;
  EdgeContainer   edges_;
  FaceContainer   faces_;
  CellContainer   cells_;

};

}
