/*
 * Copyright (c) 2017-2026 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput/VTKP3Writer.hpp"

#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"

#include "vtk/UtilityFunctions.h"

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#endif

namespace hyteg {

using walberla::vtk::typeToString;

void VTKP3Writer::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::P3 ) == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   uint_t numberOfPoints      = 0;
   uint_t numberOfCells       = 0;
   uint_t numberOfDoFsPerUnit = 0;
   if ( mgr.write2D_ )
   {
      numberOfDoFsPerUnit = levelinfo::num_microvertices_per_face( level ) + levelinfo::num_microfaces_per_face( level ) +
                            2u * levelinfo::num_microedges_per_face( level );
      numberOfPoints = storage->getNumberOfLocalFaces() * numberOfDoFsPerUnit;
      numberOfCells  = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level );
   }
   else
   {
      // NOTE: There is no levelinfo::num_microfaces_per_cell() at the moment!
      //
      // numberOfDoFsPerUnit = levelinfo::num_microvertices_per_cell( level ) + levelinfo::num_microfaces_per_cell( level ) +
      //                       2u * levelinfo::num_microedges_per_cell( level );
      numberOfDoFsPerUnit = 0;
      numberOfPoints      = storage->getNumberOfLocalCells() * numberOfDoFsPerUnit;
      numberOfCells       = storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level );
   }

   vtk::writePieceHeader( output, numberOfPoints, numberOfCells );

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   {
      VTKStreamWriter< real_t > streamWriter( mgr.vtkDataFormat_ );

      // write nodes for Vertex DoFs
      VTKMeshWriter::writePointsForMicroVertices( mgr.write2D_, streamWriter, storage, level );

      if ( mgr.write2D_ )
      {
         // write nodes for Blue Edge DoFs
         VTKMeshWriter::writePointsOnMicroEdges(
             mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_X, real_c( 1.0 / 3.0 ) );
         VTKMeshWriter::writePointsOnMicroEdges(
             mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XY, real_c( 1.0 / 3.0 ) );
         VTKMeshWriter::writePointsOnMicroEdges(
             mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Y, real_c( 1.0 / 3.0 ) );

         // write nodes for Red Edge DoFs
         VTKMeshWriter::writePointsOnMicroEdges(
             mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_X, real_c( 2.0 / 3.0 ) );
         VTKMeshWriter::writePointsOnMicroEdges(
             mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XY, real_c( 2.0 / 3.0 ) );
         VTKMeshWriter::writePointsOnMicroEdges(
             mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Y, real_c( 2.0 / 3.0 ) );
      }

      else
      {
         WALBERLA_ABORT( "VTKP3Writer::write() missing 3D implementation!" );
      }

      // write nodes for Face DoFs
      VTKMeshWriter::writePointsForMicroFaceCenters( streamWriter, storage, level );

      streamWriter.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeConnectivityP3Triangles( mgr.vtkDataFormat_, output, storage, level );
   }
   else
   {
      // VTKMeshWriter::writeConnectivityP3Tetrahedrons( mgr.vtkDataFormat_, output, storage, level, false );
   }

   output << "<PointData>\n";

   // write all scalar P3Functions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP3Functions().getFunctions< double >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP3Functions().getFunctions< float >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP3Functions().getFunctions< int32_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP3Functions().getFunctions< int64_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   // write all P3VectorFunctions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP3VectorFunctions().getFunctions< double >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP3VectorFunctions().getFunctions< float >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP3VectorFunctions().getFunctions< int32_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP3VectorFunctions().getFunctions< int64_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "</PointData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKP3Writer::writeScalarFunction( std::ostream&                              output,
                                       const P3Function< value_t >&               function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level,
                                       bool                                       write2D,
                                       vtk::DataFormat                            vtkDataFormat )

{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), 1, vtkDataFormat );

   VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   writeP3FunctionData( write2D, streamWriter, function, storage, level );
   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

template < typename value_t >
void VTKP3Writer::writeVectorFunction( std::ostream&                              output,
                                       const P3VectorFunction< value_t >&         function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level,
                                       bool                                       write2D,
                                       vtk::DataFormat                            vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   uint_t dim = function.getDimension();
   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), dim, vtkDataFormat );

   VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   writeP3VectorFunctionData( write2D, streamWriter, function, storage, level );
   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

template < typename dstStream_t, typename value_t >
void VTKP3Writer::writeP3FunctionData( bool                                       write2D,
                                       dstStream_t&                               dstStream,
                                       const P3Function< value_t >&               function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level )
{
   if ( write2D )
   {
      // first write data for vertex dofs
      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : vertexdof::macroface::Iterator( level ) )
         {
            dstStream
                << face.getData( function.getVertexDoFFunction().getFaceDataID() )
                       ->getPointer(
                           level )[vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C )];
         }
      }

      const auto orientations = {
          edgedof::EdgeDoFOrientation::X, edgedof::EdgeDoFOrientation::XY, edgedof::EdgeDoFOrientation::Y };

      // next come the data for the blue edge dofs
      for ( const auto orientation : orientations )
      {
         for ( const auto& itFaces : storage->getFaces() )
         {
            const Face& face = *itFaces.second;

            for ( const auto& it : edgedof::macroface::Iterator( level ) )
            {
               dstStream << face.getData( function.getEdgeDoFFunctionBlue().getFaceDataID() )
                                ->getPointer( level )[edgedof::macroface::index( level, it.x(), it.y(), orientation )];
            }
         }
      }

      // followed by data for the red edge dofs
      for ( const auto orientation : orientations )
      {
         for ( const auto& itFaces : storage->getFaces() )
         {
            const Face& face = *itFaces.second;

            for ( const auto& it : edgedof::macroface::Iterator( level ) )
            {
               dstStream << face.getData( function.getEdgeDoFFunctionRed().getFaceDataID() )
                                ->getPointer( level )[edgedof::macroface::index( level, it.x(), it.y(), orientation )];
            }
         }
      }

      // and finally the data of face dofs
      for ( const auto& itFaces : storage->getFaces() )
      {
         auto       faceDoFData = function.getFaceDoFFunction().dofMemory( itFaces.first, level );
         const auto memLayout   = function.getFaceDoFFunction().memoryLayout();

         for ( const auto& faceType : facedof::allFaceTypes )
         {
            for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
            {
               dstStream << faceDoFData[volumedofspace::indexing::index(
                   microFace.x(), microFace.y(), faceType, 0, 1, level, memLayout )];
            }
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "VTKP3Writer::writeP3FunctionData() lacks support for 3D!" );
   }
}

template < typename dstStream_t, typename value_t >
void VTKP3Writer::writeP3VectorFunctionData( bool                                       write2D,
                                             dstStream_t&                               dstStream,
                                             const P3VectorFunction< value_t >&         function,
                                             const std::shared_ptr< PrimitiveStorage >& storage,
                                             const uint_t&                              level )

{
   WALBERLA_ABORT( "VTKP3Writer::writeP3VectorFunctionData() not implemented!" );
}

// ========================
//  explicit instantiation
// ========================
#ifdef HYTEG_BUILD_WITH_ADIOS2
template void VTKP3Writer::writeP3FunctionData( bool                                       write2D,
                                                AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                const P3Function< real_t >&                function,
                                                const std::shared_ptr< PrimitiveStorage >& storage,
                                                const uint_t&                              level );

template void VTKP3Writer::writeP3FunctionData( bool                                         write2D,
                                                AdiosWriter::StreamAccessBuffer< uint32_t >& dstStream,
                                                const P3Function< uint32_t >&                function,
                                                const std::shared_ptr< PrimitiveStorage >&   storage,
                                                const uint_t&                                level );

template void VTKP3Writer::writeP3FunctionData( bool                                         write2D,
                                                AdiosWriter::StreamAccessBuffer< uint64_t >& dstStream,
                                                const P3Function< uint64_t >&                function,
                                                const std::shared_ptr< PrimitiveStorage >&   storage,
                                                const uint_t&                                level );

template void VTKP3Writer::writeP3VectorFunctionData( bool                                       write2D,
                                                      AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                      const P3VectorFunction< real_t >&          function,
                                                      const std::shared_ptr< PrimitiveStorage >& storage,
                                                      const uint_t&                              level );

template void VTKP3Writer::writeP3VectorFunctionData( bool                                         write2D,
                                                      AdiosWriter::StreamAccessBuffer< uint32_t >& dstStream,
                                                      const P3VectorFunction< uint32_t >&          function,
                                                      const std::shared_ptr< PrimitiveStorage >&   storage,
                                                      const uint_t&                                level );

template void VTKP3Writer::writeP3VectorFunctionData( bool                                         write2D,
                                                      AdiosWriter::StreamAccessBuffer< uint64_t >& dstStream,
                                                      const P3VectorFunction< uint64_t >&          function,
                                                      const std::shared_ptr< PrimitiveStorage >&   storage,
                                                      const uint_t&                                level );

#endif

} // namespace hyteg
