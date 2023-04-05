/*
 * Copyright (c) 2017-2022 Marcus Mohr.
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

// #include "hyteg/forms/Form.hpp"
#include "hyteg/forms/P1Form.hpp"

namespace hyteg {

class P1FormHyTeG : public P1Form
{
 public:
   virtual ~P1FormHyTeG() {}

   // implemented here to allow using the forms in form_hyteg_generated with the P1ElementwiseOperator
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const override
   {
      WALBERLA_ABORT( "integrateAll() for 2D not implemented by current HyTeG form." );
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const override
   {
      WALBERLA_ABORT( "integrateAll() for 3D not implemented by current HyTeG form." );
   }

   /// Transitional routine to allow 2D HyTeG forms inplace of FEniCS forms until we clean up the interfaces
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      Matrix3r elMat;
      this->integrateAll( coords, elMat );
      out[0] = elMat( 0, 0 );
      out[1] = elMat( 0, 1 );
      out[2] = elMat( 0, 2 );
   }

   /// Transitional routine to allow 3D HyTeG forms inplace of FEniCS forms until we clean up the interfaces
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override
   {
      Matrix4r elMat;
      this->integrateAll( coords, elMat );
      out[0] = elMat( 0, 0 );
      out[1] = elMat( 0, 1 );
      out[2] = elMat( 0, 2 );
      out[3] = elMat( 0, 3 );
   }
};

} // namespace hyteg
