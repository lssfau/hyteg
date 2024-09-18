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

#include "hyteg/dataexport/VTKOutput/VTKP2Writer.hpp"

#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"

#include "vtk/UtilityFunctions.h"

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#endif

namespace hyteg {

using walberla::vtk::typeToString;

void VTKP2Writer::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::P2 ) == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   const uint_t numberOfPoints = mgr.write2D_ ?
                                     storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level + 1 ) :
                                     storage->getNumberOfLocalCells() * levelinfo::num_microvertices_per_cell( level + 1 );

   const uint_t faceCountLevel = level;
   const uint_t cellCountLevel = level;

   const uint_t numberOfCells = mgr.write2D_ ?
                                    storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( faceCountLevel ) :
                                    storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( cellCountLevel );

   vtk::writePieceHeader( output, numberOfPoints, numberOfCells );

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   {
      VTKStreamWriter< real_t > streamWriter( mgr.vtkDataFormat_ );

      VTKMeshWriter::writePointsForMicroVertices( mgr.write2D_, streamWriter, storage, level );
      if ( mgr.write2D_ )
      {
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_X );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XY );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Y );
      }

      else
      {
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_X );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Y );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Z );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XY );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XZ );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_YZ );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XYZ );
      }

      streamWriter.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeConnectivityP2Triangles( mgr.vtkDataFormat_, output, storage, level, false );
   }
   else
   {
      VTKMeshWriter::writeConnectivityP2Tetrahedrons( mgr.vtkDataFormat_, output, storage, level, false );
   }

   output << "<PointData>\n";

   // write all scalar P2Functions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP2Functions().getFunctions< double >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2Functions().getFunctions< float >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2Functions().getFunctions< int32_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2Functions().getFunctions< int64_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   // write all P2VectorFunctions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP2VectorFunctions().getFunctions< double >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2VectorFunctions().getFunctions< float >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2VectorFunctions().getFunctions< int32_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2VectorFunctions().getFunctions< int64_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "</PointData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKP2Writer::writeScalarFunction( std::ostream&                              output,
                                       const P2Function< value_t >&               function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level,
                                       bool                                       write2D,
                                       vtk::DataFormat                            vtkDataFormat )

{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), 1, vtkDataFormat );

   VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   writeP2FunctionData( write2D, streamWriter, function, storage, level );
   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

template < typename value_t >
void VTKP2Writer::writeVectorFunction( std::ostream&                              output,
                                       const P2VectorFunction< value_t >&         function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level,
                                       bool                                       write2D,
                                       vtk::DataFormat                            vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   uint_t dim = function.getDimension();
   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), dim, vtkDataFormat );

   VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   writeP2VectorFunctionData( write2D, streamWriter, function, storage, level );
   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

template < typename dstStream_t, typename value_t >
void VTKP2Writer::writeP2FunctionData( bool                                       write2D,
                                       dstStream_t&                               dstStream,
                                       const P2Function< value_t >&               function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level )
{
   if ( write2D )
   {
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
      for ( const auto orientation : orientations )
      {
         for ( const auto& itFaces : storage->getFaces() )
         {
            const Face& face = *itFaces.second;

            for ( const auto& it : edgedof::macroface::Iterator( level ) )
            {
               dstStream << face.getData( function.getEdgeDoFFunction().getFaceDataID() )
                                ->getPointer( level )[edgedof::macroface::index( level, it.x(), it.y(), orientation )];
            }
         }
      }
   }
   else
   {
      for ( const auto& [pid, cell] : storage->getCells() )
      {
         auto vertexData = cell->getData( function.getVertexDoFFunction().getCellDataID() )->getPointer( level );

         for ( const auto& it : vertexdof::macrocell::Iterator( level ) )
         {
            dstStream
                << vertexData[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];
         }
      }

      for ( auto orientation : edgedof::allEdgeDoFOrientations )
      {
         for ( const auto& [pid, cell] : storage->getCells() )
         {
            auto edgeData = cell->getData( function.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

            if ( orientation == edgedof::EdgeDoFOrientation::XYZ )
            {
               for ( const auto& it : edgedof::macrocell::IteratorXYZ( level ) )
               {
                  dstStream << edgeData[edgedof::macrocell::index( level, it.x(), it.y(), it.z(), orientation )];
               }
            }
            else
            {
               for ( const auto& it : edgedof::macrocell::Iterator( level ) )
               {
                  dstStream << edgeData[edgedof::macrocell::index( level, it.x(), it.y(), it.z(), orientation )];
               }
            }
         }
      }
   }
}

template < typename dstStream_t, typename value_t >
void VTKP2Writer::writeP2VectorFunctionData( bool                                       write2D,
                                             dstStream_t&                               dstStream,
                                             const P2VectorFunction< value_t >&         function,
                                             const std::shared_ptr< PrimitiveStorage >& storage,
                                             const uint_t&                              level )

{
   uint_t dim = function.getDimension();

   if ( write2D )
   {
      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : vertexdof::macroface::Iterator( level, 0 ) )
         {
            for ( uint_t idx = 0; idx < dim; ++idx )
            {
               dstStream
                   << face.getData( function.component( idx ).getVertexDoFFunction().getFaceDataID() )
                          ->getPointer(
                              level )[vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C )];
            }
         }
      }

      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
         {
            for ( uint_t idx = 0; idx < dim; ++idx )
            {
               dstStream << face.getData( function.component( idx ).getEdgeDoFFunction().getFaceDataID() )
                                ->getPointer(
                                    level )[edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::X )];
            }
         }
      }

      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
         {
            for ( uint_t idx = 0; idx < dim; ++idx )
            {
               dstStream << face.getData( function.component( idx ).getEdgeDoFFunction().getFaceDataID() )
                                ->getPointer(
                                    level )[edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::XY )];
            }
         }
      }

      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
         {
            for ( uint_t idx = 0; idx < dim; ++idx )
            {
               dstStream << face.getData( function.component( idx ).getEdgeDoFFunction().getFaceDataID() )
                                ->getPointer(
                                    level )[edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::Y )];
            }
         }
      }
   }
   else
   {
      for ( const auto& [pid, cell] : storage->getCells() )
      {
         for ( const auto& it : vertexdof::macrocell::Iterator( level ) )
         {
            for ( uint_t idx = 0; idx < dim; ++idx )
            {
               auto vertexData =
                   cell->getData( function.component( idx ).getVertexDoFFunction().getCellDataID() )->getPointer( level );
               dstStream << vertexData[vertexdof::macrocell::indexFromVertex(
                   level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];
            }
         }
      }

      for ( auto orientation : edgedof::allEdgeDoFOrientations )
      {
         for ( const auto& [pid, cell] : storage->getCells() )
         {
            if ( orientation == edgedof::EdgeDoFOrientation::XYZ )
            {
               for ( const auto& it : edgedof::macrocell::IteratorXYZ( level ) )
               {
                  for ( uint_t idx = 0; idx < dim; ++idx )
                  {
                     auto edgeData =
                         cell->getData( function.component( idx ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
                     dstStream << edgeData[edgedof::macrocell::index( level, it.x(), it.y(), it.z(), orientation )];
                  }
               }
            }
            else
            {
               for ( const auto& it : edgedof::macrocell::Iterator( level ) )
               {
                  for ( uint_t idx = 0; idx < dim; ++idx )
                  {
                     auto edgeData =
                         cell->getData( function.component( idx ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
                     dstStream << edgeData[edgedof::macrocell::index( level, it.x(), it.y(), it.z(), orientation )];
                  }
               }
            }
         }
      }
   }
}

// ========================
//  explicit instantiation
// ========================
#ifdef HYTEG_BUILD_WITH_ADIOS2
template void VTKP2Writer::writeP2FunctionData( bool                                       write2D,
                                                AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                const P2Function< real_t >&                function,
                                                const std::shared_ptr< PrimitiveStorage >& storage,
                                                const uint_t&                              level );

template void VTKP2Writer::writeP2FunctionData( bool                                         write2D,
                                                AdiosWriter::StreamAccessBuffer< uint32_t >& dstStream,
                                                const P2Function< uint32_t >&                function,
                                                const std::shared_ptr< PrimitiveStorage >&   storage,
                                                const uint_t&                                level );

template void VTKP2Writer::writeP2FunctionData( bool                                         write2D,
                                                AdiosWriter::StreamAccessBuffer< uint64_t >& dstStream,
                                                const P2Function< uint64_t >&                function,
                                                const std::shared_ptr< PrimitiveStorage >&   storage,
                                                const uint_t&                                level );

template void VTKP2Writer::writeP2VectorFunctionData( bool                                       write2D,
                                                      AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                      const P2VectorFunction< real_t >&          function,
                                                      const std::shared_ptr< PrimitiveStorage >& storage,
                                                      const uint_t&                              level );

template void VTKP2Writer::writeP2VectorFunctionData( bool                                         write2D,
                                                      AdiosWriter::StreamAccessBuffer< uint32_t >& dstStream,
                                                      const P2VectorFunction< uint32_t >&          function,
                                                      const std::shared_ptr< PrimitiveStorage >&   storage,
                                                      const uint_t&                                level );

template void VTKP2Writer::writeP2VectorFunctionData( bool                                         write2D,
                                                      AdiosWriter::StreamAccessBuffer< uint64_t >& dstStream,
                                                      const P2VectorFunction< uint64_t >&          function,
                                                      const std::shared_ptr< PrimitiveStorage >&   storage,
                                                      const uint_t&                                level );

#endif

} // namespace hyteg
