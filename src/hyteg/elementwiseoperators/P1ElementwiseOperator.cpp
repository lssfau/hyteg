/*
 * Copyright (c) 2017-2019 Marcus Mohr.
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

#include "P1ElementwiseOperator.hpp"

namespace hyteg {

template < class P1Form >
P1ElementwiseOperator< P1Form >::P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        size_t                                     minLevel,
                                                        size_t                                     maxLevel,
                                                        bool                                       needsDiagEntries )
: Operator( storage, minLevel, maxLevel )
{
   if ( needsDiagEntries )
   {
      diagonalValues_ =
          std::unique_ptr< P1Function< real_t > >( new P1Function< real_t >( "diagonal entries", storage, minLevel, maxLevel ) );
      computeDiagonalOperatorValues( maxLevel, true );
   }
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::apply( const P1Function< real_t >& src,
                                             const P1Function< real_t >& dst,
                                             size_t                      level,
                                             DoFType                     flag,
                                             UpdateType                  updateType ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "apply" );

   // Make sure that halos are up-to-date (can we improve communication here?)
   communication::syncFunctionBetweenPrimitives( src, level );

   if ( updateType == Add )
   {
      WALBERLA_ABORT( "P1ElementwiseOperator::apply does not support additive update!" );
      // Note: By zeroing only the dofs in the cells' or faces' halos of dst we could
      //       actually support an additive update type. That's incorrect, as the additive
      //       communication zeros before doing the accumulation.
      //
      //       "Add" could be done, if we perform the computations on all primitives, which
      //       means that the additive communication is saved.
   }
   else
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
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getCellDataID();

         real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixVectorMultiply3D( cell, level, micro, cType, srcVertexData, dstVertexData );
            }
         }
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_ );
      dst.communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_ );
      dst.communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage_ );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         // the explicit uint_c cast prevents a segfault in intel compiler 2018.4
         // now loop over micro-faces of macro-face
         for ( yIdx = uint_c( 0 ); yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = uint_c( 1 ); xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixVectorMultiply2D(
                   face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementN, srcVertexData, dstVertexData );
               localMatrixVectorMultiply2D(
                   face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, srcVertexData, dstVertexData );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixVectorMultiply2D(
                face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, srcVertexData, dstVertexData );
         }

         // top north-west micro-element not treated, yet
         localMatrixVectorMultiply2D( face, level, 1, yIdx, P1Elements::P1Elements2D::elementNW, srcVertexData, dstVertexData );
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_ );
      dst.communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_ );
   }

   this->stopTiming( "apply" );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::smooth_jac( const P1Function< real_t >& dst,
                                                  const P1Function< real_t >& rhs,
                                                  const P1Function< real_t >& src,
                                                  real_t                      omega,
                                                  size_t                      level,
                                                  DoFType                     flag ) const
{
   this->startTiming( "smooth_jac" );

   // compute the current residual
   this->apply( src, dst, level, flag );
   dst.assign( {real_c( 1 ), real_c( -1 )}, {rhs, dst}, level, flag );

   // perform Jacobi update step
   dst.multElementwise( {( *diagonalValues_ ), dst}, level, flag );
   dst.assign( {1.0, omega}, {src, dst}, level, flag );

   this->stopTiming( "smooth_jac" );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::localMatrixVectorMultiply2D( const Face&                                face,
                                                                   const uint_t                               level,
                                                                   const uint_t                               xIdx,
                                                                   const uint_t                               yIdx,
                                                                   const P1Elements::P1Elements2D::P1Element& element,
                                                                   const real_t* const                        srcVertexData,
                                                                   real_t* const dstVertexData ) const
{
   WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

   Matrix3r                 elMat;
   Point3D                  elVecOld, elVecNew;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 3 >  dofDataIdx;
   P1Form                   form;

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

   // assemble local element vector
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   elVecOld[0] = srcVertexData[dofDataIdx[0]];
   elVecOld[1] = srcVertexData[dofDataIdx[1]];
   elVecOld[2] = srcVertexData[dofDataIdx[2]];

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   dstVertexData[dofDataIdx[0]] += elVecNew[0];
   dstVertexData[dofDataIdx[1]] += elVecNew[1];
   dstVertexData[dofDataIdx[2]] += elVecNew[2];
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::localMatrixVectorMultiply3D( const Cell&             cell,
                                                                   const uint_t            level,
                                                                   const indexing::Index&  microCell,
                                                                   const celldof::CellType cType,
                                                                   const real_t* const     srcVertexData,
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
   Matrix4r elMat;
   P1Form   form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   // assemble local element vector
   Point4D elVecOld, elVecNew;
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
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeDiagonalOperatorValues( uint_t level, bool invert )
{
   WALBERLA_ASSERT_GREATER_EQUAL( level, minLevel_ );
   WALBERLA_ASSERT_LESS_EQUAL( level, maxLevel_ );
   WALBERLA_ASSERT_NOT_NULLPTR( diagonalValues_.get() );

   // Make sure that halos are up-to-date (can we improve communication here?)
   communication::syncFunctionBetweenPrimitives( *diagonalValues_, level );

   // Zero destination before performing additive computation
   diagonalValues_->setToZero( level );

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data
         PrimitiveDataID< FunctionMemory< real_t >, Cell > diagVertexDoFIdx = diagonalValues_->getCellDataID();
         real_t*                                           diagVertexData = cell.getData( diagVertexDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               computeLocalDiagonalContributions3D( cell, level, micro, cType, diagVertexData );
            }
         }
      }

      // Push result to lower-dimensional primitives
      diagonalValues_->communicateAdditively< Cell, Face >( level );
      diagonalValues_->communicateAdditively< Cell, Edge >( level );
      diagonalValues_->communicateAdditively< Cell, Vertex >( level );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > vertexDoFIdx = diagonalValues_->getFaceDataID();
         real_t*                                           vertexData   = face.getData( vertexDoFIdx )->getPointer( level );

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               computeLocalDiagonalContributions2D( face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementN, vertexData );
               computeLocalDiagonalContributions2D( face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, vertexData );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            computeLocalDiagonalContributions2D( face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, vertexData );
         }

         // top north-west micro-element not treated, yet
         computeLocalDiagonalContributions2D( face, level, 1, yIdx, P1Elements::P1Elements2D::elementNW, vertexData );
      }

      // Push result to lower-dimensional primitives
      diagonalValues_->communicateAdditively< Face, Edge >( level );
      diagonalValues_->communicateAdditively< Face, Vertex >( level );

      // Retrieve assembled data values
      diagonalValues_->communicate< Vertex, Edge >( level );
      diagonalValues_->communicate< Edge, Face >( level );

      // Invert values if desired (note: using false below means we only invert in the interior of the primitives,
      // the values in the halos are untouched; should be okay for using diagonalValue_ in smoothers)
      if ( invert )
      {
         diagonalValues_->invertElementwise( level, All, false );
      }
   }
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeLocalDiagonalContributions2D( const Face&                                face,
                                                                           const uint_t                               level,
                                                                           const uint_t                               xIdx,
                                                                           const uint_t                               yIdx,
                                                                           const P1Elements::P1Elements2D::P1Element& element,
                                                                           real_t* const dstVertexData )
{
   Matrix3r                 elMat;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P1Form                   form;

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

   // get global indices for local dofs
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   // add local contributions to diagonal entries
   dstVertexData[dofDataIdx[0]] += elMat( 0, 0 );
   dstVertexData[dofDataIdx[1]] += elMat( 1, 1 );
   dstVertexData[dofDataIdx[2]] += elMat( 2, 2 );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeLocalDiagonalContributions3D( const Cell&             cell,
                                                                           const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           real_t* const           vertexData )
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix4r elMat;
   P1Form   form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   // add contributions for central stencil weights
   for ( uint_t k = 0; k < 4; ++k )
   {
      vertexData[vertexDoFIndices[k]] += elMat( k, k );
   }
}

#ifdef HYTEG_BUILD_WITH_PETSC

// Assemble operator as sparse matrix for PETSc
template < class P1Form >
void P1ElementwiseOperator< P1Form >::assembleLocalMatrix( Mat&                          mat,
                                                           const P1Function< PetscInt >& src,
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
         PrimitiveDataID< FunctionMemory< PetscInt >, Cell > srcVertexDoFIdx = src.getCellDataID();

         PetscInt* srcIdx = cell.getData( srcVertexDoFIdx )->getPointer( level );
         PetscInt* dstIdx = cell.getData( dstVertexDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixAssembly3D( mat, cell, level, micro, cType, srcIdx, dstIdx );
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
         PrimitiveDataID< FunctionMemory< PetscInt >, Face > srcVertexDoFIdx = src.getFaceDataID();

         PetscInt* srcIndices = face.getData( srcVertexDoFIdx )->getPointer( level );
         PetscInt* dstIndices = face.getData( dstVertexDoFIdx )->getPointer( level );

         // the explicit uint_c cast prevents a segfault in intel compiler 2018.4
         // now loop over micro-faces of macro-face
         for ( yIdx = uint_c( 0 ); yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = uint_c( 1 ); xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixAssembly2D( mat, face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementN, srcIndices, dstIndices );
               localMatrixAssembly2D( mat, face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, srcIndices, dstIndices );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixAssembly2D( mat, face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, srcIndices, dstIndices );
         }

         // top north-west micro-element not treated, yet
         localMatrixAssembly2D( mat, face, level, 1, yIdx, P1Elements::P1Elements2D::elementNW, srcIndices, dstIndices );
      }
   }
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::localMatrixAssembly2D( Mat&                                       mat,
                                                             const Face&                                face,
                                                             const uint_t                               level,
                                                             const uint_t                               xIdx,
                                                             const uint_t                               yIdx,
                                                             const P1Elements::P1Elements2D::P1Element& element,
                                                             const PetscInt* const                      srcIdx,
                                                             const PetscInt* const                      dstIdx ) const
{
   Matrix3r                 elMat;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 3 >  dofDataIdx;
   P1Form                   form;

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

   // determine global indices of our local DoFs
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   PetscInt rowIdx[3];
   rowIdx[0] = dstIdx[dofDataIdx[0]];
   rowIdx[1] = dstIdx[dofDataIdx[1]];
   rowIdx[2] = dstIdx[dofDataIdx[2]];

   PetscInt colIdx[3];
   colIdx[0] = srcIdx[dofDataIdx[0]];
   colIdx[1] = srcIdx[dofDataIdx[1]];
   colIdx[2] = srcIdx[dofDataIdx[2]];

   // add local matrix into global matrix
   PetscErrorCode ierr = MatSetValues( mat, 3, rowIdx, 3, colIdx, elMat.data(), ADD_VALUES );
   WALBERLA_ASSERT_EQUAL( ierr, 0 )
   WALBERLA_UNUSED( ierr );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::localMatrixAssembly3D( Mat&                    mat,
                                                             const Cell&             cell,
                                                             const uint_t            level,
                                                             const indexing::Index&  microCell,
                                                             const celldof::CellType cType,
                                                             const PetscInt* const   srcIdx,
                                                             const PetscInt* const   dstIdx ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix4r elMat;
   P1Form   form;
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFDataIdx;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFDataIdx );

   PetscInt rowIdx[4];
   PetscInt colIdx[4];
   for ( uint_t k = 0; k < 4; ++k )
   {
      rowIdx[k] = dstIdx[vertexDoFDataIdx[k]];
      colIdx[k] = srcIdx[vertexDoFDataIdx[k]];
   }

   // add local matrix into global matrix
   PetscErrorCode ierr = MatSetValues( mat, 4, rowIdx, 4, colIdx, elMat.data(), ADD_VALUES );
   WALBERLA_ASSERT_EQUAL( ierr, 0 )
   WALBERLA_UNUSED( ierr );
}

#endif

// P1ElementwiseLaplaceOperator
template class P1ElementwiseOperator<
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >;

// P1ElementwisePolarLaplaceOperator
template class P1ElementwiseOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > >;

// P1ElementwiseMassOperator
template class P1ElementwiseOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >;

// P1ElementwiseBlendingMassOperator
template class P1ElementwiseOperator< P1Form_mass >;

// P1ElementwiseBlendingMassOperator3D
template class P1ElementwiseOperator< P1Form_mass3D >;

// P1ElementwiseBlendingLaplaceOperator
template class P1ElementwiseOperator< P1Form_laplace >;

} // namespace hyteg
