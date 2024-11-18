/*
* Copyright (c) 2017-2024 Nils Kohl, Marcus Mohr.
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

#include <variant>

#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

namespace hyteg::micromesh {

/// \brief Container for FE functions that define the geometric positioning of micro elements.
///
/// This data structure can be used to work with super-, sub-, and isoparametric finite elements.
/// Currently implements conforming (C^0, piecewise polynomial) linear and quadratic mappings using Lagrangian basis functions.
///
/// There are currently essentially three approaches to approximate domains in HyTeG:
///
/// 1. No further approximation (sometimes referred to 'affine' - although somewhat misleading).
///    In that case the domain has to be approximated only through the coarse mesh. This enables extremely fast kernels but only
///    limited approximations to complex geometries if the coarse mesh is actually kept coarse. Clearly the coarse mesh can be
///    of arbitrarily high resolution - but this diminishes the advantages of a block-structured grid.
///
/// 2. "Blending"
///    An analytical mapping (referred to as below "blending map") can be defined that maps from the coarse mesh to the actual
///    domain. That mapping has to be differentiable and its Jacobian needs to be available to the compute kernels. If the domain
///    can be described analytically, this approach is arbitrarily exact (the error is technically only limited by the quadrature
///    error) and requires no additional memory. On the downside, it is hard to construct mappings for arbitrary domains, and the
///    evaluation of the Jacobian can be expensive for complicated maps.
///
/// 3. Parametric mappings (is this what we want to call it?)
///    The geometry is approximated by piecewise polynomials that typically match the degree of the finite element approximation
///    used for the solution (then called isoparametric). Essentially, the node positions are themselves stored as a finite
///    element function of corresponding degree. Note that the degree of the mesh approximation and the finite element solution
///    need not be the same (then referred to as sub- or superparametric).
///    The advantages and disadvantages are somehow reversed compared to the blending approach. The advantage of this approach is
///    that arbitrary geometries can be approximated and no analytical expression needs to be available. The computation of the
///    Jacobians only depends on the degree of the mesh approximation, and especially for linear/low order mappings is comparably
///    cheap. On the downside, three scalars have to be stored per node, and the accuracy is limited by the chosen polynomial
///    degree.
///
/// The MicroMesh class is a container that stores the mesh to implement the third approach outlined above.
///
/// After construction of a MicroMesh instance, all vertices are mapped to the positions they have in the standard, refined coarse
/// mesh. In other words, if no map is applied to the MicroMesh after construction, then the result is what is described in
/// approach 1, but more expensive.
///
/// \note You need to use this feature with appropriate operators. See tests and examples.
///
/// \note Combining step 2 and 3 is not implemented. If a MicroMesh is added to the storage, blending functions are ignored.
///       However, be careful and just use only approach 2 OR approach 3 in the first place to avoid unexpected issues.
///
/// \note Moving meshes (aka interpolating new values into the mesh repeatedly) are technically available through this
///       implementation but need to be handled carefully since the _reference_ frame changes with the mesh.
///       Concretely, if a map f(x) is interpolated, initially x is the coordinate induced by the refined coarse mesh.
///       During the second call to interpolate, you get f(f(x)), i.e., the input coordinate of the function you pass is the
///       result of the previously interpolated value.
///       To avoid these issues, you can _first_ interpolate the original refined coarse mesh and _then_ the apply a map to the
///       resulting, 'original' mesh. This is currently not happening automatically, especially, since the behavior described
///       above might be desired.
///
class MicroMesh
{
 public:
   /// \brief Allocates data for a MicroMesh and interpolates the node locations of the refined coarse grid.
   ///
   /// Note that if a GeometryMap is present, the node location after application of the GeometryMap is applied.
   /// If this is not desired, i.e., if the "affine" coordinates of the original refined mesh shall be applied, either the
   /// GeometryMaps have to be removed first (or need to be overwritten with the IdentityMap) or the mesh nodes have to be
   /// overwritten by a subsequent call to micromesh::interpolateRefinedCoarseMesh() with a disabled blending parameter.
   ///
   /// \param storage the underlying PrimitiveStorage
   /// \param minLevel minimum refinement level
   /// \param maxLevel maximum refinement level
   /// \param polynomialDegree polynomial degree of the piecewise mesh approximation - currently only 1 and 2 are supported
   /// \param dimension dimension of the space the mesh is embedded in - this is independent of the PrimitiveStorage, e.g., useful
   ///                  for 2D manifolds in 3D space
   MicroMesh( const std::shared_ptr< PrimitiveStorage >& storage,
              uint_t                                     minLevel,
              uint_t                                     maxLevel,
              uint_t                                     polynomialDegree,
              uint_t                                     dimension );

   /// Construct a MicroMesh directly from some already allocated vector function.
   explicit MicroMesh( const std::shared_ptr< P1VectorFunction< real_t > >& mesh );

   /// Construct a MicroMesh directly from some already allocated vector function.
   explicit MicroMesh( const std::shared_ptr< P2VectorFunction< real_t > >& mesh );

   /// Returns the polynomial degree of the mesh approximation.
   [[nodiscard]] uint_t polynomialDegree() const;

   /// Returns the dimension of the space the mesh is embedded in.
   [[nodiscard]] uint_t dimension() const;

   /// Returns the P1 mesh. Could be a nullptr if not allocated.
   [[nodiscard]] std::shared_ptr< P1VectorFunction< real_t > > p1Mesh() const;

   /// Returns the P2 mesh. Could be a nullptr if not allocated.
   [[nodiscard]] std::shared_ptr< P2VectorFunction< real_t > > p2Mesh() const;

   /// Returns the mesh. Could be a nullptr if not allocated.
   [[nodiscard]] std::variant< std::shared_ptr< P1VectorFunction< real_t > >, std::shared_ptr< P2VectorFunction< real_t > > >
       mesh() const;

 private:
   std::shared_ptr< P1VectorFunction< real_t > > p1_;
   std::shared_ptr< P2VectorFunction< real_t > > p2_;
};

/// \brief Returns the position of any micro-vertex of the MicroMesh.
///
/// If no MicroMesh was allocated and added to the PrimitiveStorage, it defaults to returning the position with respect to the
/// refined and potentially blended coarse mesh. Thus, this function can (and should) be called safely whenever the position of
/// a micro-vertex is requested.
///
Point3D microVertexPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                             PrimitiveID                                primitiveId,
                             uint_t                                     level,
                             const indexing::Index&                     microVertexIndex );

/// \brief Returns the position of the center of a micro-edge.
///
/// If no MicroMesh was allocated and added to the PrimitiveStorage, it defaults to returning the position with respect to the
/// refined and potentially blended coarse mesh. Thus, this function can (and should) be called safely whenever the position of
/// the center of a micro-edge is requested.
///
Point3D microEdgeCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                primitiveId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microVertexIndexA,
                                 const indexing::Index&                     microVertexIndexB );

/// \brief Returns the position of the center of a micro-edge.
///
/// If no MicroMesh was allocated and added to the PrimitiveStorage, it defaults to returning the position with respect to the
/// refined and potentially blended coarse mesh. Thus, this function can (and should) be called safely whenever the position of
/// the center of a micro-edge is requested.
///
Point3D microEdgeCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                primitiveId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microEdgeIndex,
                                 const edgedof::EdgeDoFOrientation&         microEdgeOrientation );

/// \brief Returns the position of the "center" of a micro-face.
///
/// The method determines and returns the position of the "center" of the micro-face specified by the given microFaceIndex.
/// In this context we use as "center" of the micro-face the position to which the barycenter of the reference triangle is
/// mapped by which ever mapping (affine, blending, isoparametric) is in use. Thus, if no MicroMesh was allocated and added
/// to the PrimitiveStorage either an affine mapping is employed, giving the barycenter of the triangular element on the
/// physical domain, or the result of applying the blending map to the barycenter of the triangular element on the
/// computational domain.
Point3D microFaceCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                faceId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microFaceIndex,
                                 facedof::FaceType                          faceType );

/// Communicates the MicroMesh such that all ghost-layers are updated.
/// Relevant, e.g., for VTK output.
void communicate( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level );

/// Communicates the MicroMesh such that all ghost-layers are updated.
/// Relevant, e.g., for VTK output.
void communicate( MicroMesh& microMesh, uint_t level );

/// Interpolates a map onto the mesh function. The input point is the current position of the respective mesh node.
/// Thus, repeated evaluation of the same map might result in different meshes each time.
/// To interpolate the mesh node coordinates of the refined coarse mesh, use interpolateRefinedCoarseMesh().
void interpolate( MicroMesh& microMesh, const std::vector< std::function< real_t( const Point3D& ) > >& map, uint_t level );

/// Interpolates a map onto the mesh function. The input point is the current position of the respective mesh node.
/// Thus, repeated evaluation of the same map might result in different meshes each time.
/// To interpolate the mesh node coordinates of the refined coarse mesh, use interpolateRefinedCoarseMesh().
void interpolate( const std::shared_ptr< PrimitiveStorage >&                      storage,
                  const std::vector< std::function< real_t( const Point3D& ) > >& map,
                  uint_t                                                          level );

/// Combines interpolate() and communicate() for convenience.
void interpolateAndCommunicate( MicroMesh&                                                      microMesh,
                                const std::vector< std::function< real_t( const Point3D& ) > >& map,
                                uint_t                                                          level );

/// Combines interpolate() and communicate() for convenience.
void interpolateAndCommunicate( const std::shared_ptr< PrimitiveStorage >&                      storage,
                                const std::vector< std::function< real_t( const Point3D& ) > >& map,
                                uint_t                                                          level );

/// Interpolates the node locations of the refined coarse mesh with or without blending.
/// Can be interpreted as "resetting" the mesh if withBlending is set to false (which ignores present GeometryMaps).
/// Can be used to interpolate the node positions of a "blended" mesh to the micro mesh.
void interpolateRefinedCoarseMesh( MicroMesh& microMesh, uint_t level, bool withBlending );

/// Interpolates the node locations of the refined coarse mesh with or without blending.
/// Can be interpreted as "resetting" the mesh if withBlending is set to false (which ignores present GeometryMaps).
/// Can be used to interpolate the node positions of a "blended" mesh to the micro mesh.
void interpolateRefinedCoarseMesh( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level, bool withBlending );

/// Convenience function that wraps a GeometryMap's evalF() method into a vector of std::functions.
/// Helpful to (re-)use an existing GeometryMap via the parametric approach. Might be computationally cheaper for complicated
/// GeometryMaps since the GeometryMap is only approximated and possibly expensive Jacobian evaluations are replaced with the
/// Jacobian evaluations of the P1 and P2 parametric mappings respectively.
std::vector< std::function< real_t( const Point3D& ) > > microMeshMapFromGeometryMap( const GeometryMap& geometryMap );

} // namespace hyteg::micromesh
