/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFRotationOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1RotationOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

using walberla::real_t;
/**
 * @brief Class that performs rotations on the finite element vector field with respect to the normal function
 *        The rotation matrix is calculated as given in, Engelman 1982.
 */
class P2RotationOperator : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >
{
 public:
   /**
    * @brief Constructs a P2RotationOperator object.
    * @param storage, minLevel, maxLevel Must be trivial by now!
    * @param normalFunction Lambda function evaluating the normal (vector normalized) on the boundary.
    */
   P2RotationOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                       size_t                                                   minLevel,
                       size_t                                                   maxLevel,
                       const std::function< void( const Point3D&, Point3D& ) >& normalFunction );

   ~P2RotationOperator() override = default;

   /**
    * @brief Rotates the vector function and transform into poalr/spherical coordinates.
    *        In 2D: (ux, uy) --> (uTangential, uRadial)
    *        In 3D: (ux, uy, uz) --> (uTheta1, uTheta2, uRadial)
    *               
    *        For 3D, the two theta components may not be exactly what is intended, 
    *        this is because we could choose infinitely many pairs of orthogonal 
    *        directions on a plane (here the tangential plane to the surface at a point)
    *        Here we choose however is done in Engelman 1982
    */
   void rotate( const P2Function< real_t >& dst_u,
                const P2Function< real_t >& dst_v,
                const P2Function< real_t >& dst_w,
                size_t                      level,
                DoFType                     flag,
                bool                        transpose = false ) const;

   void rotate( const P2VectorFunction< real_t >& dst, size_t level, DoFType flag, bool transpose = false ) const
   {
      rotate( dst[0], dst[1], dst[2], level, flag, transpose );
   }

   void rotate( const P2P1TaylorHoodFunction< real_t >& dst, size_t level, DoFType flag, bool transpose = false ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  numU,
                  const P2Function< idx_t >&                  numV,
                  const P2Function< idx_t >&                  numW,
                  uint_t                                      level,
                  DoFType                                     flag,
                  bool                                        transpose = false ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            num,
                  uint_t                                      level,
                  DoFType                                     flag,
                  bool                                        transpose = false ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P2VectorFunction< idx_t >&            dst,
                  uint_t                                      level,
                  DoFType                                     flag,
                  bool                                        transpose ) const
   {
      WALBERLA_UNUSED( dst );
      toMatrix( mat, src, level, flag, transpose );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&,
                  const P2VectorFunction< idx_t >&,
                  const P2VectorFunction< idx_t >&,
                  uint_t,
                  DoFType ) const override
   {
      WALBERLA_ABORT( "Do not use this variant of toMatrix for RotationOperator, please use the one with transpose flag" );
   }

 private:
   P1RotationOperator      p1Operator;
   EdgeDoFRotationOperator edgeDoFOperator;
};

} // namespace hyteg
