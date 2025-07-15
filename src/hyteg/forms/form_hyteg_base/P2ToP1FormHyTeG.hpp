/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "core/Abort.h"

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/Form.hpp"
#include "hyteg/forms/P2Form.hpp"
#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

class P2ToP1FormHyTeG : public Form
{
 public:
   virtual ~P2ToP1FormHyTeG() {}

   virtual void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 3, 6 >& elMat ) const = 0;

   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 4, 10 >& elMat ) const = 0;

   /// Transitional routines to allow HyTeG forms inplace of FEniCS forms until we clean up the interfaces
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix< real_t, 3, 6 > localStiffnessMatrix{ Matrix< real_t, 3, 6 >::Zero() };
      integrateAll( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix< real_t, 3, 6 > localStiffnessMatrix{ Matrix< real_t, 3, 6 >::Zero() };
      integrateAll( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 3 );
      out[1] = localStiffnessMatrix( 0, 4 );
      out[2] = localStiffnessMatrix( 0, 5 );
   }

   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      Matrix< real_t, 4, 10 > localStiffnessMatrix{ Matrix< real_t, 4, 10 >::Zero() };
      integrateAll( coords, localStiffnessMatrix );
      int rowIdx = fenics::P2DoFMap[0][0];
      out[0]     = localStiffnessMatrix( rowIdx, fenics::P2DoFMap[0][0] );
      out[1]     = localStiffnessMatrix( rowIdx, fenics::P2DoFMap[1][1] );
      out[2]     = localStiffnessMatrix( rowIdx, fenics::P2DoFMap[2][2] );
      out[3]     = localStiffnessMatrix( rowIdx, fenics::P2DoFMap[3][3] );
   }

   real_t integrate( const std::array< Point3D, 4 >&     coords,
                     const P2Form::dofPosByVertexPair3D& cntrPos,
                     const P2Form::dofPosByVertexPair3D& leafPos ) const
   {
      Matrix< real_t, 4, 10 > elMat{ Matrix< real_t, 4, 10 >::Zero() };
      integrateAll( coords, elMat );
      int rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
      int colIdx = fenics::P2DoFMap[leafPos[0]][leafPos[1]];

      return walberla::real_c( elMat( rowIdx, colIdx ) );
   }

   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                    const P2Form::dofPosByVertexPair3D&                cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      WALBERLA_ABORT( "Missing implementation" );
   }

 private:
   virtual void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const
   {
      WALBERLA_ABORT( "not implemented" );
   };

   virtual void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const
   {
      WALBERLA_ABORT( "not implemented" );
   };
};

} // namespace hyteg
