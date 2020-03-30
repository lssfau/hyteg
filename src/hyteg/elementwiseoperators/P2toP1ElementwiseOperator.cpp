/*
 * Copyright (c) 2017-2020 Marcus Mohr, Nils Kohl.
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

#include "P2toP1ElementwiseOperator.hpp"

namespace hyteg {

template < class P2toP1Form >
P2toP1ElementwiseOperator< P2toP1Form >::P2toP1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    size_t                                     minLevel,
                                                                    size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
{}

template < class P2toP1Form >
void P2toP1ElementwiseOperator< P2toP1Form >::apply( const P2Function< real_t >& src,
                                                     const P1Function< real_t >& dst,
                                                     size_t                      level,
                                                     DoFType                     flag,
                                                     UpdateType                  updateType ) const
{
   this->startTiming( "apply" );

   // Make sure that halos are up-to-date (can we improve communication here?)
   communication::syncP2FunctionBetweenPrimitives( src, level );

   if ( updateType == Replace )
   {
      // We need to zero the destination array (including halos)
      dst.setToZero( level );
   }

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< real_t >, Cell > dstVertexDoFIdx = dst.getCellDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getVertexDoFFunction().getCellDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcEdgeDoFIdx = src.getEdgeDoFFunction().getCellDataID();

         real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = cell.getData( srcEdgeDoFIdx )->getPointer( level );

         if ( updateType == Add )
         {
            // Zero out dst halos only - then during additive comm we skip zeroing
            // the data on the lower-dim primitives.

            for ( const auto& idx : vertexdof::macrocell::Iterator( level ) )
            {
               if ( !vertexdof::macrocell::isOnCellFace( idx, level ).empty() )
               {
                  auto arrayIdx           = vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() );
                  dstVertexData[arrayIdx] = real_c( 0 );
               }
            }
         }

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixVectorMultiply3D( cell, level, micro, cType, srcVertexData, srcEdgeData, dstVertexData );
            }
         }
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
   }

   else
   {
      // we only perform computations on face primitives
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         Point3D x0( face.coords[0] );
         Point3D x1( face.coords[1] );
         Point3D x2( face.coords[2] );

         uint_t                   rowsize       = levelinfo::num_microvertices_per_edge( level );
         uint_t                   inner_rowsize = rowsize;
         uint_t                   xIdx, yIdx;
         Point3D                  v0, v1, v2;
         indexing::Index          nodeIdx;
         indexing::IndexIncrement offset;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = face.getData( srcEdgeDoFIdx )->getPointer( level );

         if ( updateType == Add )
         {
            // Zero out dst halos only - then during additive comm we skip zeroing
            // the data on the lower-dim primitives.

            for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
            {
               if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
               {
                  auto arrayIdx           = vertexdof::macroface::index( level, idx.x(), idx.y() );
                  dstVertexData[arrayIdx] = real_c( 0 );
               }
            }
         }

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixVectorMultiply2D(
                   face, level, xIdx, yIdx, P2Elements::P2Face::elementN, srcVertexData, srcEdgeData, dstVertexData );
               localMatrixVectorMultiply2D(
                   face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexData, srcEdgeData, dstVertexData );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixVectorMultiply2D(
                face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexData, srcEdgeData, dstVertexData );
         }

         // top north-west micro-element not treated, yet
         localMatrixVectorMultiply2D(
             face, level, 1, yIdx, P2Elements::P2Face::elementNW, srcVertexData, srcEdgeData, dstVertexData );
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
   }

   this->stopTiming( "apply" );
}

template < class P2toP1Form >
void P2toP1ElementwiseOperator< P2toP1Form >::localMatrixVectorMultiply2D( const Face&                  face,
                                                                           const uint_t                 level,
                                                                           const uint_t                 xIdx,
                                                                           const uint_t                 yIdx,
                                                                           const P2Elements::P2Element& element,
                                                                           const real_t* const          srcVertexData,
                                                                           const real_t* const          srcEdgeData,
                                                                           real_t* const                dstVertexData ) const
{
   WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

   Matrix< real_t, 3, 6 >   elMat;
   Point6D                  elVecOld;
   Point3D                  elVecNew;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P2toP1Form               form;

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( {v0, v1, v2}, elMat );

   // assemble local element vector (note the tweaked ordering to go along with FEniCS indexing)
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
   dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
   dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

   elVecOld[0] = srcVertexData[dofDataIdx[0]];
   elVecOld[1] = srcVertexData[dofDataIdx[1]];
   elVecOld[2] = srcVertexData[dofDataIdx[2]];

   elVecOld[3] = srcEdgeData[dofDataIdx[3]];
   elVecOld[4] = srcEdgeData[dofDataIdx[4]];
   elVecOld[5] = srcEdgeData[dofDataIdx[5]];

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   dstVertexData[dofDataIdx[0]] += elVecNew[0];
   dstVertexData[dofDataIdx[1]] += elVecNew[1];
   dstVertexData[dofDataIdx[2]] += elVecNew[2];
}

template < class P2toP1Form >
void P2toP1ElementwiseOperator< P2toP1Form >::localMatrixVectorMultiply3D( const Cell&             cell,
                                                                           const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           const real_t* const     srcVertexData,
                                                                           const real_t* const     srcEdgeData,
                                                                           real_t* const           dstVertexData ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix< real_t, 4, 10 > elMat;
   P2toP1Form              form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector
   Point10D elVecOld;
   Point4D  elVecNew;
   for ( uint_t k = 0; k < 4; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[k]];
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[k - 4]];
   }

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   for ( uint_t k = 0; k < 4; ++k )
   {
      dstVertexData[vertexDoFIndices[k]] += elVecNew[k];
   }
}

#ifdef HYTEG_BUILD_WITH_PETSC

// Assemble operator as sparse matrix for PETSc
template < class P2toP1Form >
void P2toP1ElementwiseOperator< P2toP1Form >::assembleLocalMatrix( Mat&                          mat,
                                                                   const P2Function< PetscInt >& src,
                                                                   const P1Function< PetscInt >& dst,
                                                                   uint_t                        level,
                                                                   DoFType                       flag ) const
{
   // We currently ignore the flag provided!
   WALBERLA_UNUSED( flag );

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two indexing functions
         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > dstVertexDoFIdx = dst.getCellDataID();
         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > srcVertexDoFIdx = src.getVertexDoFFunction().getCellDataID();

         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > srcEdgeDoFIdx = src.getEdgeDoFFunction().getCellDataID();

         PetscInt* srcVertexIndices = cell.getData( srcVertexDoFIdx )->getPointer( level );
         PetscInt* dstVertexIndices = cell.getData( dstVertexDoFIdx )->getPointer( level );

         PetscInt* srcEdgeIndices = cell.getData( srcEdgeDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixAssembly3D( mat, cell, level, micro, cType, srcVertexIndices, srcEdgeIndices, dstVertexIndices );
            }
         }
      }
   }

   else
   {
      // we only perform computations on face primitives
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         Point3D x0( face.coords[0] );
         Point3D x1( face.coords[1] );
         Point3D x2( face.coords[2] );

         uint_t                   rowsize       = levelinfo::num_microvertices_per_edge( level );
         uint_t                   inner_rowsize = rowsize;
         uint_t                   xIdx, yIdx;
         Point3D                  v0, v1, v2;
         indexing::Index          nodeIdx;
         indexing::IndexIncrement offset;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< PetscInt >, Face > dstVertexDoFIdx = dst.getFaceDataID();
         PrimitiveDataID< FunctionMemory< PetscInt >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< PetscInt >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         PetscInt* srcVertexIndices = face.getData( srcVertexDoFIdx )->getPointer( level );
         PetscInt* dstVertexIndices = face.getData( dstVertexDoFIdx )->getPointer( level );

         PetscInt* srcEdgeIndices = face.getData( srcEdgeDoFIdx )->getPointer( level );

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixAssembly2D( mat,
                                      face,
                                      level,
                                      xIdx,
                                      yIdx,
                                      P2Elements::P2Face::elementN,
                                      srcVertexIndices,
                                      srcEdgeIndices,
                                      dstVertexIndices );
               localMatrixAssembly2D( mat,
                                      face,
                                      level,
                                      xIdx,
                                      yIdx,
                                      P2Elements::P2Face::elementNW,
                                      srcVertexIndices,
                                      srcEdgeIndices,
                                      dstVertexIndices );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixAssembly2D(
                mat, face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexIndices, srcEdgeIndices, dstVertexIndices );
         }

         // top north-west micro-element not treated, yet
         localMatrixAssembly2D(
             mat, face, level, 1, yIdx, P2Elements::P2Face::elementNW, srcVertexIndices, srcEdgeIndices, dstVertexIndices );
      }
   }
}

template < class P2toP1Form >
void P2toP1ElementwiseOperator< P2toP1Form >::localMatrixAssembly2D( Mat&                         mat,
                                                                     const Face&                  face,
                                                                     const uint_t                 level,
                                                                     const uint_t                 xIdx,
                                                                     const uint_t                 yIdx,
                                                                     const P2Elements::P2Element& element,
                                                                     const PetscInt* const        srcVertexIdx,
                                                                     const PetscInt* const        srcEdgeIdx,
                                                                     const PetscInt* const        dstVertexIdx ) const

{
   Matrix< real_t, 3, 6 >   elMat;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P2toP1Form               form;

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( {v0, v1, v2}, elMat );

   // determine global indices of our local DoFs (note the tweaked ordering to go along with FEniCS indexing)
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
   dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
   dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

   PetscInt rowIdx[3];
   rowIdx[0] = dstVertexIdx[dofDataIdx[0]];
   rowIdx[1] = dstVertexIdx[dofDataIdx[1]];
   rowIdx[2] = dstVertexIdx[dofDataIdx[2]];

   PetscInt colIdx[6];
   colIdx[0] = srcVertexIdx[dofDataIdx[0]];
   colIdx[1] = srcVertexIdx[dofDataIdx[1]];
   colIdx[2] = srcVertexIdx[dofDataIdx[2]];

   colIdx[3] = srcEdgeIdx[dofDataIdx[3]];
   colIdx[4] = srcEdgeIdx[dofDataIdx[4]];
   colIdx[5] = srcEdgeIdx[dofDataIdx[5]];

   // add local matrix into global matrix
   PetscErrorCode ierr = MatSetValues( mat, 3, rowIdx, 6, colIdx, elMat.data(), ADD_VALUES );
   WALBERLA_ASSERT_EQUAL( ierr, 0 )
   WALBERLA_UNUSED( ierr );
}

template < class P2toP1Form >
void P2toP1ElementwiseOperator< P2toP1Form >::localMatrixAssembly3D( Mat&                    mat,
                                                                     const Cell&             cell,
                                                                     const uint_t            level,
                                                                     const indexing::Index&  microCell,
                                                                     const celldof::CellType cType,
                                                                     const PetscInt* const   srcVertexIdx,
                                                                     const PetscInt* const   srcEdgeIdx,
                                                                     const PetscInt* const   dstVertexIdx ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix< real_t, 4, 10 > elMat;
   P2toP1Form              form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   PetscInt rowIdx[4];
   PetscInt colIdx[10];

   for ( uint_t k = 0; k < 4; ++k )
   {
      rowIdx[k] = dstVertexIdx[vertexDoFIndices[k]];
      colIdx[k] = srcVertexIdx[vertexDoFIndices[k]];
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      colIdx[k] = srcEdgeIdx[edgeDoFIndices[k - 4]];
   }

   // add local matrix into global matrix
   PetscErrorCode ierr = MatSetValues( mat, 4, rowIdx, 10, colIdx, elMat.data(), ADD_VALUES );
   WALBERLA_ASSERT_EQUAL( ierr, 0 )
   WALBERLA_UNUSED( ierr );
}

#endif

template class P2toP1ElementwiseOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;

template class P2toP1ElementwiseOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;

template class P2toP1ElementwiseOperator<
    P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

} // namespace hyteg
