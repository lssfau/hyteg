/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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

#pragma once

#include <array>
#include <map>
#include <set>
#include <vector>

#include "core/Abort.h"
#include "core/debug/Debug.h"

#include "hyteg/types/pointnd.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

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
/// - internal mesh generation for spherical shells
///
/// \ref MeshInfo instance can then be used to create domains via \ref SetupPrimitiveStorage.
///
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
/// | flavour    |      \#vertices       |\#triangles |
/// |:-----------|:---------------------:|:----------:|
/// | CRISS      | (nx+1)*(ny+1)         | 2*nx*ny    |
/// | CROSS      | (nx+1)*(ny+1)         | 2*nx*ny    |
/// | CRISSCROSS | (nx+1)*(ny+1) + nx*ny | 4*nx*ny    |
/// | DIAMOND    | (nx+1)*(ny+1) + nx*ny | 4*nx*ny    |
///
/*! \htmlonly
  <center>
  <table border="1">
  <tr>
  <td colspan="4" align="center">Sample mesh generated for a rectangle using (nx=3, ny=2)</td>
  </tr>
  <tr>
  <td><img src="Mesh_RectangleCriss.png"      width="600" height="600"/></td>
  <td><img src="Mesh_RectangleCross.png"      width="600" height="600"/></td>
  <td><img src="Mesh_RectangleCrissCross.png" width="600" height="600"/></td>
  <td><img src="Mesh_RectangleDiamond.png"    width="600" height="600"/></td>
  </tr>
  <tr>
  <td align="center">CRISS</td>
  <td align="center">CROSS</td>
  <td align="center">CRISSCROSS</td>
  <td align="center">DIAMOND</td>
  </tr>
  </table>
  </center></br>
  \endhtmlonly
*/
/// <b>Details of inline mesh generators for annuli:</b>
///
/// Meshing of a partial annulus is (conceptually) handled by meshing the corresponding
/// rectangle in cartesian coordinates. In case of a partial annulus this is given by
/// lower left vertex (rhoMin, phiMin) and upper right vertex (rhoMax, phiMax).
/// In the case of a full annulus the rectangle is then "glued" together.
///
/// For both a full and a partial annulus the same four  meshing flavours as for rectangles
/// can be specified. Note, however, that blending only works together with CRISS or CROSS,
/// but not for CRISSCROSS or DIAMOND.
///
/// Number of mesh entities for different flavour and a full annulus
/// |    \#vertices    |\#triangles  |   flavour  |
/// |:----------------:|:-----------:|:----------:|
/// |  nTan*(2*nRad+1) | 4*nTan*nRad | CRISSCROSS |
/// |  nTan*(nRad+1)   | 2*nTan*nRad | CRISS      |
/// |  nTan*(nRad+1)   | 2*nTan*nRad | CROSS      |
/// |  nTan*(2*nRad+1) | 4*nTan*nRad | DIAMOND    |
///
/// The boundaryFlags for the full annulus are set to the appropriate values of #hollowFlag, i.e. #flagOuterBoundary,
/// #flagInnerBoundary or #flagInterior.

/*! \htmlonly
  <center>
  <table>
  <tr>
  <td align="center"><img src="Mesh_AnnulusPartial.png" width="50%"/></td>
  <td align="center"><img src="Mesh_AnnulusFull.png" width="50%"/></td>
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
///
/// <b>Details of inline mesh generators for spherical shells:</b>
///
/// This class generates a tetrahedral mesh for a thick spherical shell. It
/// is based on the icosahedral grid approach, i.e. one starts by mapping
/// an icosahedron onto the unit sphere, which results in 20 spherical
/// triangles. Two of these triangles are then conceptually combined resulting
/// in 10 spherical diamonds. These are then further subdivided and the
/// resulting scaled appropriately to obtain radial layers. This leads to
/// a mesh consisting of prismatoidal elements with spherical triangles on
/// top and bottom.
///
/// This class uses the vertices of this mesh and from each prismatoid
/// generates 3 tetrahedrons. Given the following input parameters for the
/// constructor
///
/// | parameter | meaning                                  |
/// |:----------|:-----------------------------------------|
/// | ntan      | number of nodes along a diamond edge     |
/// | nrad      | number of radial layers                  |
/// | rmin      | radius of interior boundary of the shell |
/// | rmax      | radius of exterior boundary of the shell |
///
/// this results in
///
/// - no. of vertices: (10 * (ntan-1) * (ntan-1) + 2) * nrad
/// - no. of tetrahedrons: 60 * (ntan-1) * (ntan-1) * (nrad-1)
///
/*! \htmlonly
     <center>
       <img src="ShellMesh-ntan4-nrad3.png"
       width="50%"/><br/>
     Mesh resulting from parameter choice (ntan,nrad,rmin,rmax)=(5,3,1.0,2.0)
     </center>
     \endhtmlonly
*/
///
///
/// We start by generating a mesh for the unit sphere. This mesh is generated
/// by first splitting the four outer arcs of each diamond into (ntan-1)
/// sections of equal arc length and then splitting the west-to-east arcs
/// connecting these new points into sections of equal arc length such that we
/// end up with ntan x ntan nodes on each diamond.
///
/// The thick spherical shell is then meshed by scaling the spherical mesh to
/// obtain nrad radial layers whose radii split the interval [rmin,rmax] into
/// equidistant sub-intervals. Alternatively the radii of the layers can be
/// computed externally and provided to the class via the corresponding
/// constructor.
///
/*! \htmlonly
     <center>
       <img src="Diamonds-and-SphericalIndices.png"
       width="50%"/><br/>
     Numbering of diamonds and direction of tangential indices depending on
     hemisphere
     </center>
     \endhtmlonly
*/
///
/*! \htmlonly
     <center>
       <img src="Tets-per-local-Cell.png" width="50%"/><br/>
       Splitting a local cell into six tetrahedra
     </center>
     \endhtmlonly
*/
///
/// Every tetrahedron of the mesh can uniquely be addressed by a five-tuple
/// of indices \f$\Big(i_t,i_s^{(1)},i_s^{(2)},i_d,i_r\Big)\f$. These
/// indices have the following meaning an (by our convention) the following
/// ranges:
///
/// |      index       |                           range                 |             indicates              |
/// |:----------------:|:-----------------------------------------------:|:-----------------------------------|
/// | \f$i_t      \f$  | \f$\Big\{ 0, 1, \ldots, 5             \Big\}\f$ | index of tetrahedron in grid cell  |
/// | \f$i_s^{(1)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}-1\Big\}\f$ | first spherical node index         |
/// | \f$i_s^{(2)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}-1\Big\}\f$ | second spherical node index        |
/// | \f$i_d      \f$  | \f$\Big\{ 0, 1, \ldots, 9             \Big\}\f$ | index of diamond                   |
/// | \f$i_r      \f$  | \f$\Big\{ 0, 1, \ldots, n_\text{rad}-1\Big\}\f$ | index of radial layer (inside-out) |
///
/// We denote by
///
/// \f[ \mathcal{M}_\text{elem} : \Big(i_t,i_s^{(1)},i_s^{(2)},i_d,i_r\Big) \mapsto j_\text{elem} \f]
///
/// the mapping that assigns to each such five-tuple a unique one-dimensional
/// index. We construct the mapping by counting the tetrahedra starting at
/// zero and going through the indices in the given order, with \f$i_t\f$
/// being the fastest running and so on.
/// In a similar fashion each vertex of the grid can be addressed by
/// a four-tuple \f$\Big(i_s^{(1)},i_s^{(2)},i_d,i_r\Big)\f$. Here the indices
/// have an identical meaning as before, but slightly different ranges:
///
/// |      index       |                           range               |             indicates              |
/// |:----------------:|:---------------------------------------------:|:-----------------------------------|
/// | \f$i_s^{(1)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}\Big\}\f$ | first spherical node index         |
/// | \f$i_s^{(2)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}\Big\}\f$ | second spherical node index        |
/// | \f$i_d      \f$  | \f$\Big\{ 0, 1, \ldots, 9           \Big\}\f$ | index of diamond                   |
/// | \f$i_r      \f$  | \f$\Big\{ 0, 1, \ldots, n_\text{rad}\Big\}\f$ | index of radial layer (inside-out) |
///
/// Additionally the representation of a vertex by such a four-tuple is
/// non-unique. Considering the spherical grid, we see that the north
/// and south pole belong to five diamonds each, while the remaining
/// ten points of the base icosahedron (pentagonal nodes) belong to
/// three diamonds each, and the remaining nodes along a diamond edge
/// always belong to two diamonds. To construct a unique mapping
///
/// \f[ \mathcal{M}_\text{vert} : \Big(i_s^{(1)},i_s^{(2)},i_d,i_r\Big) \mapsto j_\text{vert} \f]
///
/// we introduce the following conventions:
///
/// - The spherical indices of north and south pole are given by (0,0)
/// - A diamond with index \f$i_d\f$ owns all nodes, which satisfy
///   \f{equation*}{\begin{split}
///   i_s^{(1)} &\in \Big\{ 0, 1, \ldots, n_\text{tan}-1\Big\}\\
///   i_s^{(2)} &\in \Big\{ 1, \ldots, n_\text{tan}\Big\}
///   \end{split}\f}
///   on the northern hemisphere this excludes the upper and lower
///   left edge of the diamond.
/// - The above assignment leads to the poles not belonging to any
///   diamond. We assign the north pole to diamond #1 and the south
///   pole to diamond #6.
///
/// The mapping \f$\mathcal{M}_\text{vert}\f$ is then constructed as
/// follows:
///
/// - We first index the north pole on each radial layer, going from
///   \f$i_r=0\f$ up to \f$i_r=n_\text{rad}\f$.
/// - Then we do the same for all south poles.
/// - Then we go through all vertices assigning them indices
///   \f$j_\text{vert}\f$ starting from \f$2\cdot n_\text{rad}\f$.
///   The ordering again follows the given index ordering with,
///   this time, \f$i_s^{(1)}\f$ being the fastest running index.
/*! \htmlonly
     <center>
       <img src="tuple2VertIndex.png" width="50%"/><br/>
     Visualisation of vertex ownership convention; vertices marked yellow and
     green are owned by neighbouring diamonds in the indicated direction
     </center>
     \endhtmlonly
*/
///
///
/// \note The inline mesh generators currently set all vertex and edge primitives' boundary flags on the domain boundary to
///       1 and those inside the domain, and of course all face primitives, to 0
class MeshInfo
{
 public:
   typedef uint_t IDType;

   typedef enum
   {
      CRISS,
      CROSS,
      CRISSCROSS,
      DIAMOND
   } meshFlavour;

   // Precise type of meshing approach
   typedef enum
   {
      SHELLMESH_ON_THE_FLY, //!< meshing is done on-the-fly
      SHELLMESH_CLASSIC     //!< meshing by midpoint refinement
   } shellMeshType;

   /// Possible boundary flags for a hollow body
   typedef enum
   {
      flagInterior      = 0,
      flagInnerBoundary = 1,
      flagOuterBoundary = 2
   } hollowFlag;

   class Vertex
   {
    public:
      Vertex()
      : id_( 0 )
      , coordinates_( Point3D() )
      , boundaryFlag_( 0 ){};

      Vertex( const IDType& id, const Point3D& coordinates, const uint_t& boundaryFlag )
      : id_( id )
      , coordinates_( coordinates )
      , boundaryFlag_( boundaryFlag )
      {}

      IDType  getID() const { return id_; }
      Point3D getCoordinates() const { return coordinates_; }
      uint_t  getBoundaryFlag() const { return boundaryFlag_; }

      void setCoordinates( const Point3D& coords ) { coordinates_ = coords; }
      void setBoundaryFlag( uint_t boundaryFlag ) { boundaryFlag_ = boundaryFlag; }

    private:
      IDType  id_;
      Point3D coordinates_;
      uint_t  boundaryFlag_;
   };

   class Edge
   {
    public:
      Edge()
      : boundaryFlag_( 0 ){};

      Edge( const std::array< IDType, 2 >& vertices, const uint_t& boundaryFlag )
      : vertices_( vertices )
      , boundaryFlag_( boundaryFlag )
      {}

      std::array< IDType, 2 > getVertices() const { return vertices_; }
      uint_t                  getBoundaryFlag() const { return boundaryFlag_; }

      void setBoundaryFlag( uint_t boundaryFlag ) { boundaryFlag_ = boundaryFlag; }

    private:
      std::array< IDType, 2 > vertices_;
      uint_t                  boundaryFlag_;
   };

   class Face
   {
    public:
      Face()
      : boundaryFlag_( 0 ){};

      Face( const std::vector< IDType >& vertices, const uint_t& boundaryFlag )
      : vertices_( vertices )
      , boundaryFlag_( boundaryFlag )
      {}

      std::vector< IDType > getVertices() const { return vertices_; }
      uint_t                getBoundaryFlag() const { return boundaryFlag_; }

      void setBoundaryFlag( uint_t boundaryFlag ) { boundaryFlag_ = boundaryFlag; }

    private:
      std::vector< IDType > vertices_;
      uint_t                boundaryFlag_;
   };

   class Cell
   {
    public:
      Cell()
      : boundaryFlag_( 0 ){};

      Cell( const std::vector< IDType >& vertices, const uint_t& boundaryFlag )
      : vertices_( vertices )
      , boundaryFlag_( boundaryFlag )
      {}

      std::vector< IDType > getVertices() const { return vertices_; }
      uint_t                getBoundaryFlag() const { return boundaryFlag_; }

      void setBoundaryFlag( uint_t boundaryFlag ) { boundaryFlag_ = boundaryFlag; }

    private:
      std::vector< IDType > vertices_;
      uint_t                boundaryFlag_;
   };

   /// \brief Applies a custom function to all vertices of the mesh, transforming their coordinates.
   void applyCoordinateMap( const std::function< Point3D( const Point3D& ) >& map )
   {
      for ( auto& v : vertices_ )
      {
         auto newCoords = map( v.second.getCoordinates() );
         v.second.setCoordinates( newCoords );
      }
   }

   typedef std::map< IDType, Vertex >                VertexContainer;
   typedef std::map< std::array< IDType, 2 >, Edge > EdgeContainer;
   typedef std::map< std::vector< IDType >, Face >   FaceContainer;
   typedef std::map< std::vector< IDType >, Cell >   CellContainer;

   void clear()
   {
      vertices_.clear();
      edges_.clear();
      faces_.clear();
      cells_.clear();
   }

   /// Construct empty MeshInfo (for testing)
   static MeshInfo emptyMeshInfo() { return MeshInfo(); }

   /// Construct a MeshInfo object from a file in Gmsh format
   static MeshInfo fromGmshFile( const std::string& meshFileName );

   /// Construct a MeshInfo object for a rectangular domain

   /// \param lowerLeft    coordinates of lower left corner of rectangle
   /// \param upperRight   coordinates of upper right corner of rectangle
   /// \param flavour      meshing strategy (CRISS, CROSS, CRISSCROSS or DIAMOND)
   /// \param nx           (nx+1) gives the number of vertices along the top and bottom edge of the rectangle
   /// \param ny           (ny+1) gives the number of vertices along the left and right edge of the rectangle
   static MeshInfo
       meshRectangle( const Point2D lowerLeft, const Point2D upperRight, const meshFlavour flavour, uint_t nx, uint_t ny );

   /// Construct a MeshInfo object for a partial annulus

   /// \param rhoMin       radius of inner circle of partial annulus
   /// \param rhoMax       radius of outer circle of partial annulus
   /// \param phiLeft      smaller angle of radial boundary in polar coordinates
   /// \param phiRight     larger angle of radial boundary in polar coordinates
   /// \param flavour      meshing strategy (CRISS, CROSS, CRISSCROSS or DIAMOND)
   /// \param nTan         number of tangential subdivisions (along inner and outer circle)
   /// \param nRad         number of radial subdivisions (along left and right radial boundary)
   static MeshInfo meshAnnulus( const real_t      rhoMin,
                                const real_t      rhoMax,
                                const real_t      phiLeft,
                                const real_t      phiRight,
                                const meshFlavour flavour,
                                uint_t            nTan,
                                uint_t            nRad );

   /// Construct a MeshInfo object for a full annulus
   static MeshInfo meshAnnulus( const real_t rmin, const real_t rmax, uint_t nTan, uint_t nRad );

   /// Construct a MeshInfo object for a full annulus
   static MeshInfo meshAnnulus( const real_t rmin, const real_t rmax, const meshFlavour flavour, uint_t nTan, uint_t nRad );

   /// Constuct a MeshInfo describing a unit cube discretized by 2 * 4^{level} macro-faces
   static MeshInfo meshUnitSquare( uint_t level );

   /// Constructs a MeshInfo object for a spherical shell (equidistant radial layers)
   ///
   /// \param ntan    number of nodes along spherical diamond edge
   /// \param nrad    number of radial layers
   /// \param rmin    radius of innermost shell (core-mantle-boundary)
   /// \param rmax    radius of outermost shell
   static MeshInfo meshSphericalShell( uint_t        ntan,
                                       uint_t        nrad,
                                       double        rmin,
                                       double        rmax,
                                       shellMeshType meshType = shellMeshType::SHELLMESH_CLASSIC );

   /// Constructs a MeshInfo object for a spherical shell (externally computed radial layers)
   ///
   /// \param ntan    number of nodes along spherical diamond edge
   /// \param layers  vector that gives the radii of all layers, sorted from the
   ///                CMB outwards
   static MeshInfo meshSphericalShell( uint_t                       ntan,
                                       const std::vector< double >& layers,
                                       shellMeshType                meshType = shellMeshType::SHELLMESH_CLASSIC );

   /// Constructs a MeshInfo object for a chain of triangles.
   ///
   /// Starting from the left side, numFaces faces are connected to each other in an alternating fashion
   /// to build a channel-like domain.
   /// There are no inner vertices in this mesh - each face is connected to its two neighbors only.
   /// If numFaces is not even, the channel is not rectangular but a trapezoid.
   ///
   /// \param numFaces number of faces in the mesh
   /// \param width width of the domain
   /// \param height height of the domain
   static MeshInfo meshFaceChain( uint_t numFaces, real_t width, real_t height );

   /// Constructs a MeshInfo object for a chain of triangles with regular size.
   static MeshInfo meshFaceChain( uint_t numFaces );

   /// Construct a MeshInfo object for a rectangular cuboid
   ///
   /// \param lowerLeftFront   coordinates of lower left front corner of cuboid
   /// \param upperRightBack   coordinates of upper right back corner of cuboid
   /// \param nx               (nx+1) gives the number of vertices along cuboid edges in x-direction
   /// \param ny               (ny+1) gives the number of vertices along cuboid edges in y-direction
   /// \param nz               (nz+1) gives the number of vertices along cuboid edges in z-direction
   static MeshInfo meshCuboid( const Point3D lowerLeftFront, const Point3D upperRightBack, uint_t nx, uint_t ny, uint_t nz );

   /// Construct a MeshInfo object for a symmetric, rectangular cuboid.
   /// This version of the cuboid is made up of numCubesX * numCubesY * numCubesZ 24-element-cubes,
   /// and has the following properties:
   /// - point symmetric at the center,
   /// - no edge at the corner points towards its interior.
   ///
   /// \param lowerLeftFront   coordinates of lower left front corner of cuboid
   /// \param upperRightBack   coordinates of upper right back corner of cuboid
   /// \param numCubesX        gives the number of 24-element cubes in x-direction
   /// \param numCubesY        gives the number of 24-element cubes in y-direction
   /// \param numCubesZ        gives the number of 24-element cubes in z-direction
   static MeshInfo meshSymmetricCuboid( const Point3D lowerLeftFront,
                                        const Point3D upperRightBack,
                                        uint_t        numCubesX,
                                        uint_t        numCubesY,
                                        uint_t        numCubesZ );

   /// Constructs a cubed domain from a list of cube "coordinates".
   ///
   /// \param cubeCoordinates contains the x, y, and z coordinates of all cubes that shall be meshed
   /// \param cubeType 0: "smallest" possible cubes with 6 tetrahedrons, 1: 24 tets per cube, has inner edges at all vertices
   static MeshInfo meshCubedDomain( const std::set< std::array< int, 3 > >& cubeCoordinates, int cubeType = 0 );

   /// \brief Meshes a torus around the z-axis.
   ///
   /// The GeometryMap TokamakMap provides a blending function that maps this mesh onto a tokamak geometry.
   /// The very same map can be parameterized so that is reduces to the special case of a torus.
   ///
   /// \param setupStorage the SetupPrimitiveStorage instance
   /// \param numToroidalSlices number of prisms in toroidal direction (along the ring)
   /// \param numPoloidalSlices number of vertices on the boundary of a slice through the tube
   /// \param radiusOriginToCenterOfTube distance from origin to the center of the tube
   /// \param tubeLayerRadii list of radii of layers of the sliced tube - the last element defines the actual radius of the tube
   /// \param toroidalStartAngle angle (in radians) by which the domain shall be rotated about the z-axis
   /// \param poloidalStartAngle angle (in radians) by which the domain shall be rotated about the ring through the center of the tube
   ///
   static MeshInfo meshTorus( uint_t                numToroidalSlices,
                              uint_t                numPoloidalSlices,
                              real_t                radiusOriginToCenterOfTube,
                              std::vector< real_t > tubeLayerRadii,
                              real_t                toroidalStartAngle = 0,
                              real_t                poloidalStartAngle = 0 );

   /// \brief Create a mesh composed of a single triangle
   ///
   /// For testing and performance checks this convenience method allows to generate
   /// a mesh composed of a single triangle by providing its vertices.
   static MeshInfo singleTriangle( const Point2D& v1, const Point2D& v2, const Point2D& v3 );

   /// \brief Create a mesh composed of a single tetrahedron
   ///
   /// For testing and performance checks this convenience method allows to generate
   /// a mesh composed of a single tetrahedron by providing its vertices.
   static MeshInfo singleTetrahedron( const std::array< Point3D, 4 >& vertices );

   /// \brief Creates a finer coarse mesh from a given mesh
   ///
   /// Takes a given MeshInfo and refines it with the default refinment algorithm.
   /// Afterwards a new MeshInfo based on that refinement is created.
   /// This allows to refine the coarse mesh
   ///
   /// \param oldMesh Original MeshIfno
   /// \param refinmentSteps number of refinements
   static MeshInfo refinedCoarseMesh( const MeshInfo& oldMesh, uint_t refinementSteps );

   /// Returns vertices of the mesh
   const VertexContainer& getVertices() const { return vertices_; };

   /// Returns edges of the mesh
   const EdgeContainer& getEdges() const { return edges_; };

   /// Returns faces of the mesh
   const FaceContainer& getFaces() const { return faces_; }

   /// Returns cells of the mesh
   const CellContainer& getCells() const { return cells_; }

 private:
   MeshInfo(){};

   /// Adds vertex to the mesh
   void addVertex( const Vertex& vertex );

   /// Adds edge in ascending index order and performs checks
   void addEdge( const Edge& edge );

   /// Adds face in ascending index order and performs checks
   void addFace( const Face& face );

   /// Adds cell and all edges and faces
   void addCellAndAllEdgesAndFaces( const Cell& cell );

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

   /// Set the boundaryFlag of all edges automatically from that of their vertices

   /// Calling the method results in the boundaryFlag of all edges in the MeshInfo being (re)set.
   /// The flag of a single edge is set to that of its two vertices, if the latter agree. Otherwise
   /// it will be set to the value of the undecided argument.
   void deduceEdgeFlagsFromVertices( uint_t flagInconsistent = 0 );

   /// Set the boundaryFlag of all faces automatically from that of their vertices

   /// Calling the method results in the boundaryFlag of all faces in the MeshInfo being (re)set.
   /// The flag of a single face is set to that of its three vertices, if the latter agree. Otherwise
   /// it will be set to the value of the undecided argument.
   void deduceFaceFlagsFromVertices( uint_t flagInconsistent = 0 );

   VertexContainer vertices_;
   EdgeContainer   edges_;
   FaceContainer   faces_;
   CellContainer   cells_;
};

} // namespace hyteg
