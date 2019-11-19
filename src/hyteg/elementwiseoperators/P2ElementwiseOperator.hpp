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
#pragma once

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

using walberla::real_t;

template < class P2Form >
class P2ElementwiseOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          bool                                       needsDiagEntries = true );

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

   std::unique_ptr< P2Function< real_t > > diagonalValues_;

   void smooth_jac( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const P2Function< real_t >& src,
                    real_t                      omega,
                    size_t                      level,
                    DoFType                     flag ) const;

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   ///
   /// \note The src and dst data array must not be identical.
   void localMatrixVectorMultiply2D( const Face&                  face,
                                     const uint_t                 level,
                                     const uint_t                 xIdx,
                                     const uint_t                 yIdx,
                                     const P2Elements::P2Element& element,
                                     const real_t* const          srcVertexData,
                                     const real_t* const          srcEdgeData,
                                     real_t* const                dstVertexData,
                                     real_t* const                dstEdgeData ) const;

   void computeDiagonalOperatorValues( uint_t level );

   /// Compute contributions to operator diagonal for given micro-face
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   void   computeLocalDiagonalContributions2D( Face&                        face,
                                               const uint_t                 level,
                                               const uint_t                 xIdx,
                                               const uint_t                 yIdx,
                                               const P2Elements::P2Element& element,
                                               real_t* const                dstVertexData,
                                               real_t* const                dstEdgeData );
   P2Form form_;

   // std::unique_ptr< P2Function< real_t > > diagonalValues_;

   bool recomputeDiagonal_ = true;
};

typedef P2ElementwiseOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >
    P2ElementwiseLaplaceOperator;
typedef P2ElementwiseOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
    P2ElementwiseMassOperator;

} // namespace hyteg
