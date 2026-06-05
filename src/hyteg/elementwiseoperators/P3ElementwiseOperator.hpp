/*
 * Copyright (c) 2017-2026 Marcus Mohr, Andreas Burkhart.
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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_hyteg_base/P3FormHyTeG.hpp"
#include "hyteg/forms/form_hyteg_generated/p3/p3_diffusion_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p3/p3_mass_affine_qe.hpp"
#include "hyteg/forms/form_hyteg_generated/p3/p3_mass_blending_q6.hpp"
#include "hyteg/operators/Operator.hpp"

// Just for testing compilation, remove this during implementation!!!!
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p3functionspace/P3Function.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::real_t;

template < class P3FormHyTeG >
class P3ElementwiseOperator : public Operator< P3Function< real_t >, P3Function< real_t > >,
                              public WeightedJacobiSmoothable< P3Function< real_t > >,
                              public OperatorWithInverseDiagonal< P3Function< real_t > >
{
 public:
   P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          const P3FormHyTeG&                         form );

   P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          bool                                       needsInverseDiagEntries );

   P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          const P3FormHyTeG&                         form,
                          bool                                       needsInverseDiagEntries );

   void gemv( const real_t&               alpha,
              const P3Function< real_t >& src,
              const real_t&               beta,
              const P3Function< real_t >& dst,
              uint_t                      level,
              DoFType                     flag ) const override final;

   void applyScaled( const real_t&               alpha,
                     const P3Function< real_t >& src,
                     const P3Function< real_t >& dst,
                     uint_t                      level,
                     DoFType                     flag,
                     UpdateType                  updateType = Replace ) const override final;

   void apply( const P3Function< real_t >& src,
               const P3Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false, false, real_c( 1 ) ); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() override final { computeDiagonalOperatorValues( true, false, real_c( 1 ) ); }

   void computeLumpedDiagonalOperatorValues() { computeDiagonalOperatorValues( false, true, real_c( 1 ) ); }
   void computeLumpedInverseDiagonalOperatorValues() { computeDiagonalOperatorValues( true, true, real_c( 1 ) ); }

   /// Trigger (re)computation of scaled inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValuesScaled( const real_t& alpha ) override final
   {
      computeDiagonalOperatorValues( true, false, alpha );
   }

   void computeLumpedInverseDiagonalOperatorValuesScaled( const real_t& alpha )
   {
      computeDiagonalOperatorValues( true, true, alpha );
   }

   /// Trigger (re)computation of scaled diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValuesScaled( const real_t& alpha ) { computeDiagonalOperatorValues( false, false, alpha ); }

   void computeLumpedDiagonalOperatorValuesScaled( const real_t& alpha ) { computeDiagonalOperatorValues( false, true, alpha ); }

   std::shared_ptr< P3Function< real_t > > getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          diagonalValues_,
          "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return diagonalValues_;
   };

   std::shared_ptr< P3Function< real_t > > getInverseDiagonalValues() const override
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          inverseDiagonalValues_,
          "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return inverseDiagonalValues_;
   };

   std::shared_ptr< P3Function< real_t > > getLumpedDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          lumpedDiagonalValues_,
          "Lumped diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return lumpedDiagonalValues_;
   };

   std::shared_ptr< P3Function< real_t > > getLumpedInverseDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          lumpedInverseDiagonalValues_,
          "Lumped inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return lumpedInverseDiagonalValues_;
   };

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   void computeAndStoreLocalElementMatrices();

   void smooth_jac_scaled( const real_t&               alpha,
                           const P3Function< real_t >& dst,
                           const P3Function< real_t >& rhs,
                           const P3Function< real_t >& src,
                           real_t                      omega,
                           uint_t                      level,
                           DoFType                     flag ) const override;

   void smooth_jac( const P3Function< real_t >& dst,
                    const P3Function< real_t >& rhs,
                    const P3Function< real_t >& src,
                    real_t                      omega,
                    size_t                      level,
                    DoFType                     flag ) const override;

   void smooth_gs( const P3Function< real_t >&, const P3Function< real_t >&, size_t, DoFType ) const
   {
      WALBERLA_ABORT( "Gauss-Seidel not implemented for P3ElementwiseOperator." )
   }

   void smooth_sor( const P3Function< real_t >&, const P3Function< real_t >&, real_t, size_t, DoFType ) const
   {
      WALBERLA_ABORT( "SOR not implemented for P3ElementwiseOperator." )
   }

   /// Assemble operator as sparse matrix with scaling
   ///
   /// \param alpha constant scaling of the matrix
   /// \param mat   a sparse matrix proxy
   /// \param src   P3Function for determining column indices
   /// \param dst   P3Function for determining row indices
   /// \param level le2el in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrixScaled( const real_t&                               alpha,
                        const std::shared_ptr< SparseMatrixProxy >& mat,
                        const P3Function< idx_t >&                  src,
                        const P3Function< idx_t >&                  dst,
                        uint_t                                      level,
                        DoFType                                     flag ) const override;

   /// Assemble operator as sparse matrix.
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   P3Function for determining column indices
   /// \param dst   P3Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P3Function< idx_t >&                  src,
                  const P3Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override;

   P3FormHyTeG getForm() const;

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microFace      index associated with the current element = micro-face
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   /// \param elMat          the 10x10 element matrix to be multiplied
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply2D( uint_t                 level,
                                     const indexing::Index& microFace,
                                     facedof::FaceType      fType,
                                     const real_t* const    srcVertexData,
                                     const real_t* const    srcEdgeDataBlue,
                                     const real_t* const    srcEdgeDataRed,
                                     const real_t* const    srcFaceData,
                                     real_t* const          dstVertexData,
                                     real_t* const          dstEdgeDataBlue,
                                     real_t* const          dstEdgeDataRed,
                                     real_t* const          dstFaceData,
                                     const Matrix10r&       elMat,
                                     const real_t&          alpha ) const;

   /// Compute contributions to operator diagonal for given micro-face
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   void computeLocalDiagonalContributions2D( const Face&             face,
                                             const uint_t            level,
                                             const indexing::Index&  microFace,
                                             const facedof::FaceType faceType,
                                             real_t* const           dstVertexData,
                                             real_t* const           dstEdgeDataBlue,
                                             real_t* const           dstEdgeDataRed,
                                             real_t* const           dstFaceData,
                                             const real_t&           alpha,
                                             bool                    lumped );

   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const indexing::Index&                      microFace,
                               const facedof::FaceType                     faceType,
                               std::array< const idx_t*, 8 >               indexDataPointer,
                               const real_t&                               alpha ) const;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert, bool lumped, const real_t& alpha );

   std::shared_ptr< P3Function< real_t > > diagonalValues_;
   std::shared_ptr< P3Function< real_t > > inverseDiagonalValues_;
   std::shared_ptr< P3Function< real_t > > lumpedDiagonalValues_;
   std::shared_ptr< P3Function< real_t > > lumpedInverseDiagonalValues_;

   P3FormHyTeG form_;

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrix10r& localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType )
   {
      WALBERLA_ASSERT( !storage_->hasGlobalCells(), "Retrieving local element matrix for 2D in 3D run. Why?" )
      const auto idx = facedof::macroface::index( level, microFace.x(), microFace.y(), fType );
      WALBERLA_ASSERT( localElementMatrices2D_.count( face.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices2D_.at( face.getID() ).at( level ).empty() )
      return localElementMatrices2D_[face.getID()][level][idx];
   }

   /// \brief Returns a const reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   const Matrix10r&
       localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType ) const
   {
      WALBERLA_ASSERT( !storage_->hasGlobalCells(), "Retrieving local element matrix for 2D in 3D run. Why?" )
      const auto idx = facedof::macroface::index( level, microFace.x(), microFace.y(), fType );
      WALBERLA_ASSERT( localElementMatrices2D_.count( face.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices2D_.at( face.getID() ).at( level ).empty() )
      return localElementMatrices2D_.at( face.getID() ).at( level ).at( idx );
   }

   bool localElementMatricesPrecomputed_;

   /// Pre-computed local element matrices.
   /// localElementMatrices2D_[macroCellID][level][cellIdx] = mat10x10
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrix10r, Eigen::aligned_allocator< Matrix10r > > > >
       localElementMatrices2D_;
};

template < class P3FormHyTeG >
void assembleLocalElementMatrix2D( const Face&            face,
                                   uint_t                 level,
                                   const indexing::Index& microFace,
                                   facedof::FaceType      fType,
                                   P3FormHyTeG            form,
                                   Matrix10r&             elMat )
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 3 > verts = facedof::macroface::getMicroVerticesFromMicroFace< true >( microFace, fType );
   std::array< Point3D, 3 >         coords;
   for ( uint_t k = 0; k < 3; ++k )
   {
      coords[k] = vertexdof::macroface::coordinateFromIndex( level, face, verts[k] );
   }

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( coords, elMat );
}

typedef P3ElementwiseOperator< forms::p3_mass_affine_qe > P3ElementwiseMassOperator;
typedef P3ElementwiseOperator< forms::p3_mass_blending_q6 > P3ElementwiseBlendingMassOperator;
typedef P3ElementwiseOperator< forms::p3_diffusion_affine_q4 > P3ElementwiseDiffusionOperator;

} // namespace hyteg
