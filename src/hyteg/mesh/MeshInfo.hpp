/*
 * Copyright (c) 2017-2025 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann, Andreas Burkhart.
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

#include "hyteg/mesh/HyTeGMeshDir.hpp"
#include "hyteg/types/PointND.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

// forward declare friend classes
class GmshReaderForMSH22;
class GmshReaderForMSH41;

/// \brief Contains information about a mesh
///
/// The \ref MeshInfo class is used to store information about meshes that are
/// either constructed by reading mesh files or generated internally for some
/// standard geometries. It provides static methods to construct meshes by
/// - importing a file in Gmsh format (see method \ref fromGmshFile)
/// - internal mesh generation for rectangles
/// - internal mesh generation for full or partial annuli
/// - internal mesh generation for spherical shells
/// and more.
///
/// \ref A MeshInfo instance can then be used to create domains via \ref SetupPrimitiveStorage.
///
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
      , coordinates_( Point3D::Zero() )
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
   ///
   /// HyTeG supports reading Gmesh files in MSH2.2 and MSH4.1 format.
   ///
   /// \note If the mesh is in MSH4.1 format it is possible to have the reader set the primitives'
   /// Primitive::meshBoundaryFlag_ from the Gmsh's physicalTags. However, for this to work, the
   /// file needs to contain an Entities section giving physical tags for all entities, i.e. Nodes,
   /// Curves, Surfaces and in 3D also Volumes, which are used in the mesh itself.
   /// See Gmsh's documentation on `Physical Curve` etc. If there is not 'Entities' section in the
   /// file or not every entity has at least one physical tag associated with it, the reader
   /// will abort.
   ///
   /// For importPhysicalTags = false, the reader will flag all primitives as 0, i.e. inner, and
   /// marking primitives and degrees of freedom to be e.g. DirichletBoundary needs to be handled
   /// inside of the application.
   ///
   /// \note The MSH2.2 reader always assumes that there are precisely two flags for each element
   /// in the Element section, of which we take the first one to be the boundary flag, if the
   /// element is an edge or a face. Vertices and tetrahedra, are always flagged as 0, i.e. inner.
   /// Yeah, this seems not to be very consistent.
   static MeshInfo fromGmshFile( const std::string& meshFileName, bool importPhysicalTags = true );

   /// @name Friend classes
   /// The readers for MSH format need access to the internals of MeshInfo
   /// @{
   friend class GmshReaderForMSH22;
   friend class GmshReaderForMSH41;
   /// @}

   /// Construct a MeshInfo object for a rectangular domain
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
   \endhtmlonly */
   ///
   /// \param lowerLeft    coordinates of lower left corner of rectangle
   /// \param upperRight   coordinates of upper right corner of rectangle
   /// \param flavour      meshing strategy (CRISS, CROSS, CRISSCROSS or DIAMOND)
   /// \param nx           (nx+1) gives the number of vertices along the top and bottom edge of the rectangle
   /// \param ny           (ny+1) gives the number of vertices along the left and right edge of the rectangle
   static MeshInfo
       meshRectangle( const Point2D lowerLeft, const Point2D upperRight, const meshFlavour flavour, uint_t nx, uint_t ny );

   /// Construct a MeshInfo object for a partial annulus
   ///
   /// See documentation of meshAnnulus( const real_t, const real_t, const meshFlavour, uint_t, uint_t )
   ///
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
   ///
   /// Meshing of a partial annulus
   /// (see
   ///   meshAnnulus( const real_t,
   ///                const real_t,
   ///                const real_t,
   ///                const real_t,
   ///                const meshFlavour,
   ///                uint_t,
   ///                uint_t ) )
   /// is (conceptually) handled by meshing the corresponding
   /// rectangle in cartesian coordinates. In case of a partial annulus this is given by
   /// lower left vertex (rhoMin, phiMin) and upper right vertex (rhoMax, phiMax).
   ///
   /// In the case of a full annulus (covered by this function) the rectangle is then "glued" together.
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
   <td align="center">partial annulus (nTan=4, nRad=2, flavour=CRISSCROSS)</td>
   <td align="center">full annulus (nTan=15, nRad=2, flavour=CRISSCROSS)</td>
   </tr>
   </table>
   </center>
   \endhtmlonly */
   ///
   static MeshInfo meshAnnulus( const real_t rmin, const real_t rmax, const meshFlavour flavour, uint_t nTan, uint_t nRad );

   /// Constuct a MeshInfo describing a unit cube discretized by 2 * 4^{level} macro-faces
   static MeshInfo meshUnitSquare( uint_t level );

   /** Constructs a MeshInfo object for a spherical shell (equidistant radial layers)
    *
    *  <b>Overview</b>
    *
    *  This function generates a tetrahedral mesh for a thick spherical shell. It
    *  is based on the icosahedral grid approach, i.e. one starts by mapping
    *  an icosahedron onto the unit sphere, which results in 20 spherical
    *  triangles. Pairs of these triangles are then conceptually combined resulting
    *  in 10 spherical diamonds. These are then further subdivided and
    *  scaled appropriately to obtain radial layers. This leads to
    *  a mesh consisting of prismatoidal elements with spherical triangles on
    *  top and bottom.
    *
    *  This function uses the vertices of this mesh and from each prismatoid
    *  generates 3 tetrahedrons. Given the following input parameters for the
    *  constructor
    *
    *  | parameter | constraints                                           | meaning                                  |
    *  |:----------|:------------------------------------------------------|:-----------------------------------------|
    *  | ntan      | (ntan-1) must be a power of 2 (for SHELLMESH_CLASSIC) | number of nodes along a diamond edge     |
    *  | nrad      | nrad >= 2                                             | number of radial layers (a radial layer here is a 2D-manifold in 3D space such that each tet is enclosed by two layers)                 |
    *  | rmin      | rmin < rmax                                           | radius of interior boundary of the shell (innermost layer) |
    *  | rmax      | rmin < rmax                                           | radius of exterior boundary of the shell (outermost layer) |
    *
    *  this results in
    *
    *  - no. of vertices: (10 * (ntan-1) * (ntan-1) + 2) * nrad
    *  - no. of tetrahedrons: 60 * (ntan-1) * (ntan-1) * (nrad-1)
    *
    *  <br/>
    * \htmlonly
   <center>
   <img src="SphericalShell_coarse_ntan_5_diamond.png" width="50%"/><br/>
   Surface of the spherical shell mesh for ntan = 5, highlighting the spherical diamond (split into two spherical triangles)
   and the ntan = 5 vertices along one of the edges of the diamond.
   </center>
   \endhtmlonly
   *
   * \htmlonly
   <center>
   <img src="SphericalShell_coarse_ntan_5_nrad_2_clip.png" width="30%"/>
   <img src="SphericalShell_coarse_ntan_5_nrad_3_clip.png" width="30%"/>
   <img src="SphericalShell_coarse_ntan_5_nrad_4_clip.png" width="30%"/>
   <br/>
   Cut through spherical shell meshes with ntan = 5, rmin = 0.5, rmax = 1.0, and nrad = {2, 3, 4} (from left to right increasing).
   </center>
   \endhtmlonly
    *  <br/>
    *
    *  <b>Boundary flags</b>
    *
    *  The boundaryFlags for the shell mesh are set to the appropriate values of #hollowFlag, i.e. #flagOuterBoundary,
    *  #flagInnerBoundary or #flagInterior.
    *
    *  <b>Details</b>
    *
    *  We start by generating a mesh for the unit sphere. This mesh is generated
    *  by first splitting the four outer arcs of each diamond into (ntan-1)
    *  sections of equal arc length and then splitting the west-to-east arcs
    *  connecting these new points into sections of equal arc length such that we
    *  end up with ntan x ntan nodes on each diamond.
    *
    *  The thick spherical shell is then meshed by scaling the spherical mesh to
    *  obtain nrad radial layers whose radii split the interval [rmin,rmax] into
    *  equidistant sub-intervals. Alternatively the radii of the layers can be
    *  computed externally and provided to the class via the corresponding
    *  constructor.
    *
    * \htmlonly
   <center>
   <img src="Diamonds-and-SphericalIndices.png" width="50%"/><br/>
   Numbering of diamonds and direction of tangential indices depending on the hemisphere
   </center>
   \endhtmlonly
    *
    * \htmlonly
   <center>
   <img src="Tets-per-local-Cell.png" width="50%"/><br/>
   Splitting a local cell into six tetrahedra
   </center>
   \endhtmlonly
    *  <br/>
    *
    *  Every tetrahedron of the mesh can uniquely be addressed by a five-tuple
    *  of indices \f$\Big(i_t,i_s^{(1)},i_s^{(2)},i_d,i_r\Big)\f$. These
    *  indices have the following meaning an (by our convention) the following
    *  ranges:
    *
    *  |      index       |                           range                 |             indicates              |
    *  |:----------------:|:-----------------------------------------------:|:-----------------------------------|
    *  | \f$i_t      \f$  | \f$\Big\{ 0, 1, \ldots, 5             \Big\}\f$ | index of tetrahedron in grid cell  |
    *  | \f$i_s^{(1)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}-1\Big\}\f$ | first spherical node index         |
    *  | \f$i_s^{(2)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}-1\Big\}\f$ | second spherical node index        |
    *  | \f$i_d      \f$  | \f$\Big\{ 0, 1, \ldots, 9             \Big\}\f$ | index of diamond                   |
    *  | \f$i_r      \f$  | \f$\Big\{ 0, 1, \ldots, n_\text{rad}-1\Big\}\f$ | index of radial layer (inside-out) |
    *
    *  We denote by
    *
    *  \f[ \mathcal{M}_\text{elem} : \Big(i_t,i_s^{(1)},i_s^{(2)},i_d,i_r\Big) \mapsto j_\text{elem} \f]
    *
    *  the mapping that assigns to each such five-tuple a unique one-dimensional
    *  index. We construct the mapping by counting the tetrahedra starting at
    *  zero and going through the indices in the given order, with \f$i_t\f$
    *  being the fastest running and so on.
    *  In a similar fashion each vertex of the grid can be addressed by
    *  a four-tuple \f$\Big(i_s^{(1)},i_s^{(2)},i_d,i_r\Big)\f$. Here the indices
    *  have an identical meaning as before, but slightly different ranges:
    *
    *  |      index       |                           range               |             indicates              |
    *  |:----------------:|:---------------------------------------------:|:-----------------------------------|
    *  | \f$i_s^{(1)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}\Big\}\f$ | first spherical node index         |
    *  | \f$i_s^{(2)}\f$  | \f$\Big\{ 0, 1, \ldots, n_\text{tan}\Big\}\f$ | second spherical node index        |
    *  | \f$i_d      \f$  | \f$\Big\{ 0, 1, \ldots, 9           \Big\}\f$ | index of diamond                   |
    *  | \f$i_r      \f$  | \f$\Big\{ 0, 1, \ldots, n_\text{rad}\Big\}\f$ | index of radial layer (inside-out) |
    *
    *  Additionally the representation of a vertex by such a four-tuple is
    *  non-unique. Considering the spherical grid, we see that the north
    *  and south pole belong to five diamonds each, while the remaining
    *  ten points of the base icosahedron (pentagonal nodes) belong to
    *  three diamonds each, and the remaining nodes along a diamond edge
    *  always belong to two diamonds. To construct a unique mapping
    *
    *  \f[ \mathcal{M}_\text{vert} : \Big(i_s^{(1)},i_s^{(2)},i_d,i_r\Big) \mapsto j_\text{vert} \f]
    *
    *  we introduce the following conventions:
    *
    *  - The spherical indices of north and south pole are given by (0,0)
    *  - A diamond with index \f$i_d\f$ owns all nodes, which satisfy
    *    \f{equation*}{\begin{split}
    *    i_s^{(1)} &\in \Big\{ 0, 1, \ldots, n_\text{tan}-1\Big\}\\
    *    i_s^{(2)} &\in \Big\{ 1, \ldots, n_\text{tan}\Big\}
    *    \end{split}\f}
    *    on the northern hemisphere this excludes the upper and lower
    *    left edge of the diamond.
    *  - The above assignment leads to the poles not belonging to any
    *    diamond. We assign the north pole to diamond #1 and the south
    *    pole to diamond #6.
    *
    *  The mapping \f$\mathcal{M}_\text{vert}\f$ is then constructed as
    *  follows:
    *
    *  - We first index the north pole on each radial layer, going from
    *    \f$i_r=0\f$ up to \f$i_r=n_\text{rad}\f$.
    *  - Then we do the same for all south poles.
    *  - Then we go through all vertices assigning them indices
    *    \f$j_\text{vert}\f$ starting from \f$2\cdot n_\text{rad}\f$.
    *    The ordering again follows the given index ordering with,
    *    this time, \f$i_s^{(1)}\f$ being the fastest running index.
    * \htmlonly
   <center>
   <img src="tuple2VertIndex.png" width="50%"/><br/>
   Visualisation of vertex ownership convention; vertices marked yellow and
   green are owned by neighbouring diamonds in the indicated direction
   </center>
   \endhtmlonly
    *
    *
    *  \param ntan      number of nodes along spherical diamond edge
    *  \param nrad      number of radial layers (>= 2)
    *                   (a radial layer here is a 2D-manifold in 3D space such that each tet is enclosed by two layers)
    *  \param rmin      radius of innermost shell (core-mantle-boundary)
    *  \param rmax      radius of outermost shell
    *  \param meshType  allows selecting meshing strategy, defaults to SHELLMESH_CLASSIC (only change this, if you know what you
    *                  are doing)
   **/
   static MeshInfo meshSphericalShell( uint_t        ntan,
                                       uint_t        nrad,
                                       real_t        rmin,
                                       real_t        rmax,
                                       shellMeshType meshType = shellMeshType::SHELLMESH_CLASSIC );

   /// Constructs a MeshInfo object for a spherical shell (externally computed positions of radial layers)
   ///
   /// See documentation of meshSphericalShell( uint_t, uint_t, real_t, real_t, shellMeshType )
   ///
   /// \param ntan      number of nodes along spherical diamond edge
   /// \param layers    vector that gives the radii of all layers, sorted from the
   ///                  CMB outwards
   /// \param meshType  allows selecting meshing strategy, defaults to SHELLMESH_CLASSIC (only change this, if you know what you
   ///                  are doing)
   static MeshInfo meshSphericalShell( uint_t                       ntan,
                                       const std::vector< real_t >& layers,
                                       shellMeshType                meshType = shellMeshType::SHELLMESH_CLASSIC );

   /// Constructs a MeshInfo object for a thin spherical shell
   ///
   /// The method creates an icosahedral mesh for a thin shell, i.e. a 2D manifold, of given radius. It uses the
   /// SHELLMESH_CLASSIC approach for this. All primitives are marked as being hollowFlag::flagInterior.
   ///
   /// \param ntan    number of nodes along spherical diamond edge
   /// \param radius  the shell's radius
   static MeshInfo meshThinSphericalShell( uint_t ntan, real_t radius );

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
   /// \param toroidalResolution number of prisms in toroidal direction (along the ring) in a complete (360 degree) ring
   /// \param poloidalResolution number of vertices on the boundary of a slice through the tube
   /// \param radiusOriginToCenterOfTube distance from origin to the center of the tube
   /// \param tubeLayerRadii list of radii of layers of the sliced tube - the last element defines the actual radius of the tube
   /// \param toroidalStartAngle angle (in radians) by which the domain shall be rotated about the z-axis
   /// \param poloidalStartAngle angle (in radians) by which the domain shall be rotated about the ring through the center of the tube
   /// \param numToroidalSlices number of prisms in toroidal direction in the final mesh (the default value 0 produces a complete ring)
   ///
   static MeshInfo meshTorus( uint_t                toroidalResolution,
                              uint_t                poloidalResolution,
                              real_t                radiusOriginToCenterOfTube,
                              std::vector< real_t > tubeLayerRadii,
                              real_t                toroidalStartAngle = 0,
                              real_t                poloidalStartAngle = 0,
                              uint_t                numToroidalSlices  = 0 );

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
   /// Takes a given MeshInfo and refines it with the default refinement algorithm.
   /// Afterwards a new MeshInfo based on that refinement is created.
   /// This allows to refine the coarse mesh
   ///
   /// \param oldMesh Original MeshIfno
   /// \param refinementSteps number of refinements
   static MeshInfo refinedCoarseMesh( const MeshInfo& oldMesh, uint_t refinementSteps );

   static MeshInfo refinedCoarseMesh2D( const MeshInfo& oldMesh, uint_t refinementSteps );

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

   /// Auxilliary function for meshSphericalShell
   void computeSphericalShellVertices( uint_t ntan, const std::vector< real_t >& layers, MeshInfo::shellMeshType meshType );

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

   VertexContainer vertices_{};
   EdgeContainer   edges_{};
   FaceContainer   faces_{};
   CellContainer   cells_{};

   /// Function for second pass over primitives in the MSH-readers
   ///
   /// The function inserts the primitive information found in the file's Elements
   /// section into MeshInfo object
   void processPrimitivesFromGmshFile( const EdgeContainer& parsedEdges,
                                       const FaceContainer& parsedFaces,
                                       const CellContainer& parsedCells,
                                       bool                 inheritParentBoundaryFlag = false );
};

} // namespace hyteg
