/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {

template < typename ValueType >
void P1toP2Conversion( const P1Function< ValueType >& src,
                       const P2Function< ValueType >& dst,
                       const uint_t&                  P2Level,
                       const DoFType&                 flag = All )
{
   const auto        P1Level = P2Level + 1;
   auto              storage = src.getStorage();
   BoundaryCondition bc      = src.getBoundaryCondition();

   storage->getTimingTree()->start( "P1toP2Conversion" );

   for ( const auto& it : storage->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( bc.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data = vertex.getData( src.getVertexDataID() )->getPointer( P1Level );
         auto P2Data = vertex.getData( dst.getVertexDoFFunction().getVertexDataID() )->getPointer( P2Level );
         P2Data[0]   = P1Data[0];
      }
   }

   for ( const auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( bc.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data   = edge.getData( src.getEdgeDataID() )->getPointer( P1Level );
         auto P2Data_v = edge.getData( dst.getVertexDoFFunction().getEdgeDataID() )->getPointer( P2Level );
         auto P2Data_e = edge.getData( dst.getEdgeDoFFunction().getEdgeDataID() )->getPointer( P2Level );

         for ( auto itIdx : vertexdof::macroedge::Iterator( P2Level ) )
         {
            P2Data_v[vertexdof::macroedge::index( P2Level, itIdx.x() )] =
                P1Data[vertexdof::macroedge::index( P1Level, itIdx.x() * 2 )];
         }
         for ( auto itIdx : edgedof::macroedge::Iterator( P2Level ) )
         {
            P2Data_e[edgedof::macroedge::index( P2Level, itIdx.x() )] =
                P1Data[vertexdof::macroedge::index( P1Level, itIdx.x() * 2 + 1 )];
         }
      }
   }

   for ( const auto& it : storage->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( bc.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data   = face.getData( src.getFaceDataID() )->getPointer( P1Level );
         auto P2Data_v = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( P2Level );
         auto P2Data_e = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( P2Level );

         for ( auto itIdx : vertexdof::macroface::Iterator( P2Level ) )
         {
            P2Data_v[vertexdof::macroface::index( P2Level, itIdx.x(), itIdx.y() )] =
                P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 )];
         }
         for ( auto itIdx : edgedof::macroface::Iterator( P2Level ) )
         {
            P2Data_e[edgedof::macroface::index( P2Level, itIdx.x(), itIdx.y(), edgedof::EdgeDoFOrientation::X )] =
                P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 )];

            P2Data_e[edgedof::macroface::index( P2Level, itIdx.x(), itIdx.y(), edgedof::EdgeDoFOrientation::XY )] =
                P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 + 1 )];

            P2Data_e[edgedof::macroface::index( P2Level, itIdx.x(), itIdx.y(), edgedof::EdgeDoFOrientation::Y )] =
                P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 + 1 )];
         }
      }
   }

   for ( const auto& it : storage->getCells() )
   {
      Cell& cell = *it.second;

      if ( testFlag( bc.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data   = cell.getData( src.getCellDataID() )->getPointer( P1Level );
         auto P2Data_v = cell.getData( dst.getVertexDoFFunction().getCellDataID() )->getPointer( P2Level );
         auto P2Data_e = cell.getData( dst.getEdgeDoFFunction().getCellDataID() )->getPointer( P2Level );

         for ( auto itIdx : vertexdof::macrocell::Iterator( P2Level ) )
         {
            P2Data_v[vertexdof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z() )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2, itIdx.z() * 2 )];
         }
         for ( auto itIdx : edgedof::macrocell::Iterator( P2Level ) )
         {
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::X )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2, itIdx.z() * 2 )];

            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::Y )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 + 1, itIdx.z() * 2 )];

            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::Z )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2, itIdx.z() * 2 + 1 )];

            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::XY )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 + 1, itIdx.z() * 2 )];

            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::XZ )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2, itIdx.z() * 2 + 1 )];

            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::YZ )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 + 1, itIdx.z() * 2 + 1 )];
         }

         for ( auto itIdx : edgedof::macrocell::IteratorXYZ( P2Level ) )
         {
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::XYZ )] =
                P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 + 1, itIdx.z() * 2 + 1 )];
         }
      }
   }

   storage->getTimingTree()->stop( "P1toP2Conversion" );
}

template < typename ValueType >
void P1toP2Conversion( const P1VectorFunction< ValueType >& src,
                       const P2VectorFunction< ValueType >& dst,
                       const uint_t&                        P2Level,
                       const DoFType&                       flag )
{
   for ( uint_t k = 0; k < src.getDimension(); k++ )
   {
      P1toP2Conversion( src[k], dst[k], P2Level, flag );
   }
}

// -------------------------
//  Explicit instantiations
// -------------------------
template void P1toP2Conversion< real_t >( const P1Function< real_t >& src,
                                          const P2Function< real_t >& dst,
                                          const uint_t&               P2Level,
                                          const DoFType&              flag );

template void P1toP2Conversion( const P1VectorFunction< real_t >& src,
                                const P2VectorFunction< real_t >& dst,
                                const uint_t&                     P2Level,
                                const DoFType&                    flag );

} // namespace hyteg
