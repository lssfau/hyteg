/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

namespace hyteg {

class VTKOutput;

class VTKMeshWriter
{
 public:
   /// \brief Writes the point coordinates for all micro vertices.
   ///
   /// \param write2D       flag for toggling 2D and 3D
   /// \param dstStream     an object that behaves like an output stream, i.e. we can write to it using the << operator
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level to write
   /// \param discontinuous if true, vertices are written for _each_ element - that means eventually each vertex is written
   ///                      multiple times so that discontinuous elements can be output
   template < typename dstStream_t >
   static void writePointsForMicroVertices( bool                                       write2D,
                                            dstStream_t&                               dstStream,
                                            const std::shared_ptr< PrimitiveStorage >& storage,
                                            uint_t                                     level,
                                            bool                                       discontinuous = false );

   /// \brief Writes the point coordinates for all micro edges, i.e. midpoint coordinates, of a certain edge type
   ///
   /// \param write2D       flag for toggling 2D and 3D
   /// \param dstStream     an object that behaves like an output stream, i.e. we can write to it using the << operator
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level to write
   /// \param dofType       selects the type of edge, e.g. horizontal, for which midpoints get written
   template < typename dstStream_t >
   static void writePointsForMicroEdges( bool                                       write2D,
                                         dstStream_t&                               dstStream,
                                         const std::shared_ptr< PrimitiveStorage >& storage,
                                         uint_t                                     level,
                                         const vtk::DoFType&                        dofType );

   /// \brief Writes the point coordinates for all micro face centers.
   ///
   /// \param dstStream     an object that behaves like an output stream, i.e. we can write to it using the << operator
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level to write
   template < typename dstStream_t >
   static void writePointsForMicroFaceCenters( dstStream_t&                               dstStream,
                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                               uint_t                                     level );

   /// \brief Writes the 2D cells.
   ///
   /// \param vtkDataFormat format in which data is stored in the VTK file
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param faceWidth     number of micro-vertices per macro-edge
   /// \param discontinuous should match what is specified in VTKMeshWriter::writePointsForMicroVertices()
   static void writeCells2D( vtk::DataFormat                            vtkDataFormat,
                             std::ostream&                              output,
                             const std::shared_ptr< PrimitiveStorage >& storage,
                             uint_t                                     faceWidth,
                             bool                                       discontinuous = false );

   /// \brief Writes connectivity information for 2D cells of type VTK_QUADRATIC_TRIANGLE
   ///
   /// \param vtkDataFormat format in which data is stored in the VTK file
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   /// \param discontinuous only false is currently supported
   static void writeConnectivityP2Triangles( vtk::DataFormat                            vtkDataFormat,
                                             std::ostream&                              output,
                                             const std::shared_ptr< PrimitiveStorage >& storage,
                                             uint_t                                     level,
                                             bool                                       discontinuous = false );

   /// \brief Writes connectivity information for exporting P2PlusBubbleFunctions on triangles
   ///
   /// For visualisation purposes we split the each triangular microface into six subfaces by connecting
   /// the vertices and edge midpoints with the center of the face. The resulting subtriangles are stored
   /// as 2D cells of type VTK_TRIANGLE.
   ///
   /// \param vtkDataFormat format in which data is stored in the VTK file
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   static void writeConnectivityP2TrianglesWithBubble( vtk::DataFormat                            vtkDataFormat,
                                                       std::ostream&                              output,
                                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                                       uint_t                                     level );

   /// \brief Writes the 3D cells.
   ///
   /// \param vtkDataFormat format in which data is stored in the VTK file
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param faceWidth     number of micro-vertices per macro-edge
   /// \param discontinuous should match what is specified in VTKMeshWriter::writePointsForMicroVertices()
   static void writeCells3D( vtk::DataFormat                            vtkDataFormat,
                             std::ostream&                              output,
                             const std::shared_ptr< PrimitiveStorage >& storage,
                             uint_t                                     width,
                             bool                                       discontinuous = false );

   /// \brief Writes connectivity information for 3D cells of type VTK_QUADRATIC_TETRA
   ///
   /// \param vtkDataFormat format in which data is stored in the VTK file
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   /// \param discontinuous only false is currently supported
   static void writeConnectivityP2Tetrahedrons( vtk::DataFormat                            vtkDataFormat,
                                                std::ostream&                              output,
                                                const std::shared_ptr< PrimitiveStorage >& storage,
                                                uint_t                                     level,
                                                bool                                       discontinuous = false );

   /// \brief Writes element <-> node association information for 2D cells of type VTK_QUADRATIC_TRIANGLE
   ///
   /// \param dstStream     stream-like object to write data to, must support operator<<
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   template < typename dstStream_t >
   static void writeElementNodeAssociationP2Triangles( dstStream_t&                               dstStream,
                                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                                       uint_t                                     level );

   /// \brief Writes element <-> node association information for 2D cells with a bubble
   ///
   /// \param dstStream     stream-like object to write data to, must support operator<<
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   template < typename dstStream_t >
   static void writeElementNodeAssociationP2TrianglesWithBubble( dstStream_t&                               dstStream,
                                                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                                                 uint_t                                     level );

   /// \brief Writes element <-> node association information for 3D cells of type VTK_QUADRATIC_TETRA
   ///
   /// \param dstStream     stream-like object to write data to, must support operator<<
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   template < typename dstStream_t >
   static void writeElementNodeAssociationP2Tetrahedrons( dstStream_t&                               dstStream,
                                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     level );

   /// \brief Writes element <-> node association information for 2D cells of type VTK_TRIANGLE
   ///
   /// \param dstStream     stream-like object to write data to, must support operator<<
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   /// \param discontinuous if true, vertices are written for _each_ element - that means eventually each vertex is written
   ///                      multiple times so that discontinuous elements can be output
   template < typename dstStream_t >
   static void writeElementNodeAssociationP1Triangles( dstStream_t&                               dstStream,
                                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                                       uint_t                                     faceWidth,
                                                       bool                                       discontinuous = false );

   /// \brief Writes element <-> node association information for 3D cells of type VTK_TETRAHEDRON
   ///
   /// \param dstStream     stream-like object to write data to, must support operator<<
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   /// \param discontinuous if true, vertices are written for _each_ element - that means eventually each vertex is written
   ///                      multiple times so that discontinuous elements can be output
   template < typename dstStream_t >
   static void writeElementNodeAssociationP1Tetrahedrons( dstStream_t&                               dstStream,
                                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     width,
                                                          bool                                       discontinuous = false );
};

} // namespace hyteg
