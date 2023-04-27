/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput.hpp"

namespace hyteg {

class VTKOutput;

class VTKMeshWriter
{
 public:
   /// \brief Writes the point coordinates for all micro vertices.
   ///
   /// \param mgr           the corresponding VTKOutput instance
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level to write
   /// \param discontinuous if true, vertices are written for _each_ element - that means eventually each vertex is written
   ///                      multiple times so that discontinuous elements can be output
   static void writePointsForMicroVertices( const VTKOutput&                           mgr,
                                            std::ostream&                              output,
                                            const std::shared_ptr< PrimitiveStorage >& storage,
                                            uint_t                                     level,
                                            bool                                       discontinuous = false );

   static void writePointsForMicroEdges( const VTKOutput&                           mgr,
                                         std::ostream&                              output,
                                         const std::shared_ptr< PrimitiveStorage >& storage,
                                         uint_t                                     level,
                                         const vtk::DoFType&                        dofType );

   /// \brief Writes the 2D cells.
   ///
   /// \param mgr           the corresponding VTKOutput instance
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param faceWidth     number of micro-vertices per macro-edge
   /// \param discontinuous should match what is specified in VTKMeshWriter::writePointsForMicroVertices()
   static void writeCells2D( const VTKOutput&                           mgr,
                             std::ostream&                              output,
                             const std::shared_ptr< PrimitiveStorage >& storage,
                             uint_t                                     faceWidth,
                             bool                                       discontinuous = false );

   /// \brief Writes connectivity information for 2D cells of type VTK_QUADRATIC_TRIANGLE
   ///
   /// \param mgr           the corresponding VTKOutput instance
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param level         refinement level for output
   /// \param discontinuous should match what is specified in VTKMeshWriter::writePointsForMicroVertices()
   static void writeConnectivityP2Triangles( const VTKOutput&                           mgr,
                                             std::ostream&                              output,
                                             const std::shared_ptr< PrimitiveStorage >& storage,
                                             uint_t                                     level,
                                             bool                                       discontinuous = false );

   /// \brief Writes the 3D cells.
   ///
   /// \param mgr           the corresponding VTKOutput instance
   /// \param output        an output stream to write to
   /// \param storage       the associated PrimitiveStorage
   /// \param faceWidth     number of micro-vertices per macro-edge
   /// \param discontinuous should match what is specified in VTKMeshWriter::writePointsForMicroVertices()
   static void writeCells3D( const VTKOutput&                           mgr,
                             std::ostream&                              output,
                             const std::shared_ptr< PrimitiveStorage >& storage,
                             uint_t                                     width,
                             bool                                       discontinuous = false );
};

} // namespace hyteg
