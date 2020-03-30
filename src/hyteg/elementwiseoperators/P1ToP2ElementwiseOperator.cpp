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

#include "P1ToP2ElementwiseOperator.hpp"

namespace hyteg {

template < class P1toP2Form >
P1ToP2ElementwiseOperator< P1toP2Form >::P1ToP2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    size_t                                     minLevel,
                                                                    size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
{}

template < class P1toP2Form >
void P1ToP2ElementwiseOperator< P1toP2Form >::apply( const P1Function< real_t >& src,
                                                     const P2Function< real_t >& dst,
                                                     size_t                      level,
                                                     DoFType                     flag,
                                                     UpdateType                  updateType ) const
{
   this->startTiming( "apply" );

   // Make sure that halos are up-to-date (can we improve communication here?)
   communication::syncFunctionBetweenPrimitives( src, level );

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
         PrimitiveDataID< FunctionMemory< real_t >, Cell > dstVertexDoFIdx = dst.getVertexDoFFunction().getCellDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getCellDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Cell > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getCellDataID();

         real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* dstEdgeData = cell.getData( dstEdgeDoFIdx )->getPointer( level );

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

            for ( const auto& idx : edgedof::macrocell::Iterator( level ) )
            {
               for ( const auto& orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
               {
                  if ( !edgedof::macrocell::isInnerEdgeDoF( level, idx, orientation ) )
                  {
                     auto arrayIdx         = edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), orientation );
                     dstEdgeData[arrayIdx] = real_c( 0 );
                  }
               }
            }
         }

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixVectorMultiply3D( cell, level, micro, cType, srcVertexData, dstVertexData, dstEdgeData );
            }
         }
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Cell, Face >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getVertexDoFFunction().communicateAdditively< Cell, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getVertexDoFFunction().communicateAdditively< Cell, Vertex >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getEdgeDoFFunction().communicateAdditively< Cell, Face >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getEdgeDoFFunction().communicateAdditively< Cell, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* dstEdgeData = face.getData( dstEdgeDoFIdx )->getPointer( level );

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

            for ( const auto& idx : edgedof::macroface::Iterator( level ) )
            {
               for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
               {
                  if ( !edgedof::macroface::isInnerEdgeDoF( level, idx, orientation ) )
                  {
                     auto arrayIdx         = edgedof::macroface::index( level, idx.x(), idx.y(), orientation );
                     dstEdgeData[arrayIdx] = real_c( 0 );
                  }
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
                   face, level, xIdx, yIdx, P2Elements::P2Face::elementN, srcVertexData, dstVertexData, dstEdgeData );
               localMatrixVectorMultiply2D(
                   face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexData, dstVertexData, dstEdgeData );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixVectorMultiply2D(
                face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexData, dstVertexData, dstEdgeData );
         }

         // top north-west micro-element not treated, yet
         localMatrixVectorMultiply2D(
             face, level, 1, yIdx, P2Elements::P2Face::elementNW, srcVertexData, dstVertexData, dstEdgeData );
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getVertexDoFFunction().communicateAdditively< Face, Vertex >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getEdgeDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
   }

   this->stopTiming( "apply" );
}

template < class P1toP2Form >
void P1ToP2ElementwiseOperator< P1toP2Form >::localMatrixVectorMultiply2D( const Face&                  face,
                                                                           const uint_t                 level,
                                                                           const uint_t                 xIdx,
                                                                           const uint_t                 yIdx,
                                                                           const P2Elements::P2Element& element,
                                                                           const real_t* const          srcVertexData,
                                                                           real_t* const                dstVertexData,
                                                                           real_t* const                dstEdgeData ) const
{
   WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

   Matrixr< 6, 3 >          elMat;
   Point3D                  elVecOld;
   Point6D                  elVecNew;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P1toP2Form               form;

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

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   dstVertexData[dofDataIdx[0]] += elVecNew[0];
   dstVertexData[dofDataIdx[1]] += elVecNew[1];
   dstVertexData[dofDataIdx[2]] += elVecNew[2];

   dstEdgeData[dofDataIdx[3]] += elVecNew[3];
   dstEdgeData[dofDataIdx[4]] += elVecNew[4];
   dstEdgeData[dofDataIdx[5]] += elVecNew[5];
}

template < class P1toP2Form >
void P1ToP2ElementwiseOperator< P1toP2Form >::localMatrixVectorMultiply3D( const Cell&             cell,
                                                                           const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           const real_t* const     srcVertexData,
                                                                           real_t* const           dstVertexData,
                                                                           real_t* const           dstEdgeData ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrixr< 10, 4 > elMat;
   P1toP2Form       form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector
   Point4D  elVecOld;
   Point10D elVecNew;
   for ( uint_t k = 0; k < 4; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[k]];
   }

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   for ( uint_t k = 0; k < 4; ++k )
   {
      dstVertexData[vertexDoFIndices[k]] += elVecNew[k];
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      dstEdgeData[edgeDoFIndices[k - 4]] += elVecNew[k];
   }
}

#ifdef HYTEG_BUILD_WITH_PETSC

// Assemble operator as sparse matrix for PETSc
template < class P1toP2Form >
void P1ToP2ElementwiseOperator< P1toP2Form >::assembleLocalMatrix( Mat&                          mat,
                                                                   const P1Function< PetscInt >& src,
                                                                   const P2Function< PetscInt >& dst,
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
         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > dstVertexDoFIdx = dst.getVertexDoFFunction().getCellDataID();
         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > srcVertexDoFIdx = src.getCellDataID();

         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getCellDataID();

         PetscInt* srcVertexIndices = cell.getData( srcVertexDoFIdx )->getPointer( level );
         PetscInt* dstVertexIndices = cell.getData( dstVertexDoFIdx )->getPointer( level );

         PetscInt* dstEdgeIndices = cell.getData( dstEdgeDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixAssembly3D( mat, cell, level, micro, cType, srcVertexIndices, dstVertexIndices, dstEdgeIndices );
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
         PrimitiveDataID< FunctionMemory< PetscInt >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< PetscInt >, Face > srcVertexDoFIdx = src.getFaceDataID();

         PrimitiveDataID< FunctionMemory< PetscInt >, Face > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getFaceDataID();

         PetscInt* srcVertexIndices = face.getData( srcVertexDoFIdx )->getPointer( level );
         PetscInt* dstVertexIndices = face.getData( dstVertexDoFIdx )->getPointer( level );

         PetscInt* dstEdgeIndices = face.getData( dstEdgeDoFIdx )->getPointer( level );

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
                                      dstVertexIndices,
                                      dstEdgeIndices );
               localMatrixAssembly2D( mat,
                                      face,
                                      level,
                                      xIdx,
                                      yIdx,
                                      P2Elements::P2Face::elementNW,
                                      srcVertexIndices,
                                      dstVertexIndices,
                                      dstEdgeIndices );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixAssembly2D(
                mat, face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexIndices, dstVertexIndices, dstEdgeIndices );
         }

         // top north-west micro-element not treated, yet
         localMatrixAssembly2D(
             mat, face, level, 1, yIdx, P2Elements::P2Face::elementNW, srcVertexIndices, dstVertexIndices, dstEdgeIndices );
      }
   }
}

template < class P1toP2Form >
void P1ToP2ElementwiseOperator< P1toP2Form >::localMatrixAssembly2D( Mat&                         mat,
                                                                     const Face&                  face,
                                                                     const uint_t                 level,
                                                                     const uint_t                 xIdx,
                                                                     const uint_t                 yIdx,
                                                                     const P2Elements::P2Element& element,
                                                                     const PetscInt* const        srcVertexIdx,
                                                                     const PetscInt* const        dstVertexIdx,
                                                                     const PetscInt* const        dstEdgeIdx ) const

{
   Matrixr< 6, 3 >          elMat;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P1toP2Form               form;

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

   PetscInt rowIdx[6];
   rowIdx[0] = dstVertexIdx[dofDataIdx[0]];
   rowIdx[1] = dstVertexIdx[dofDataIdx[1]];
   rowIdx[2] = dstVertexIdx[dofDataIdx[2]];

   rowIdx[3] = dstEdgeIdx[dofDataIdx[3]];
   rowIdx[4] = dstEdgeIdx[dofDataIdx[4]];
   rowIdx[5] = dstEdgeIdx[dofDataIdx[5]];

   PetscInt colIdx[3];
   colIdx[0] = srcVertexIdx[dofDataIdx[0]];
   colIdx[1] = srcVertexIdx[dofDataIdx[1]];
   colIdx[2] = srcVertexIdx[dofDataIdx[2]];

   // add local matrix into global matrix
   PetscErrorCode ierr = MatSetValues( mat, 6, rowIdx, 3, colIdx, elMat.data(), ADD_VALUES );
   WALBERLA_ASSERT_EQUAL( ierr, 0 )
   WALBERLA_UNUSED( ierr );
}

template < class P1toP2Form >
void P1ToP2ElementwiseOperator< P1toP2Form >::localMatrixAssembly3D( Mat&                    mat,
                                                                     const Cell&             cell,
                                                                     const uint_t            level,
                                                                     const indexing::Index&  microCell,
                                                                     const celldof::CellType cType,
                                                                     const PetscInt* const   srcVertexIdx,
                                                                     const PetscInt* const   dstVertexIdx,
                                                                     const PetscInt* const   dstEdgeIdx ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrixr< 10, 4 > elMat;
   P1toP2Form       form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   PetscInt rowIdx[10];
   PetscInt colIdx[4];

   for ( uint_t k = 0; k < 4; ++k )
   {
      rowIdx[k] = dstVertexIdx[vertexDoFIndices[k]];
      colIdx[k] = srcVertexIdx[vertexDoFIndices[k]];
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      rowIdx[k] = dstEdgeIdx[edgeDoFIndices[k - 4]];
   }

   // add local matrix into global matrix
   PetscErrorCode ierr = MatSetValues( mat, 10, rowIdx, 4, colIdx, elMat.data(), ADD_VALUES );
   WALBERLA_ASSERT_EQUAL( ierr, 0 )
   WALBERLA_UNUSED( ierr );
}

#endif

template class P1ToP2ElementwiseOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >;

template class P1ToP2ElementwiseOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >;

template class P1ToP2ElementwiseOperator<
    P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;

} // namespace hyteg
