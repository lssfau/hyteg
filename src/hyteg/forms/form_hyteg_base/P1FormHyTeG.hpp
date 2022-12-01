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
      Point3D                  matrixRow;
      std::array< Point3D, 3 > vertexCoords( coords );

      // integrate first row
      this->integrate( vertexCoords, matrixRow );
      elMat( 0, 0 ) = matrixRow[0];
      elMat( 0, 1 ) = matrixRow[1];
      elMat( 0, 2 ) = matrixRow[2];

      // integrate second row
      vertexCoords[0] = coords[1];
      vertexCoords[1] = coords[2];
      vertexCoords[2] = coords[0];
      this->integrate( vertexCoords, matrixRow );
      elMat( 1, 0 ) = matrixRow[2];
      elMat( 1, 1 ) = matrixRow[0];
      elMat( 1, 2 ) = matrixRow[1];

      // integrate third row
      vertexCoords[0] = coords[2];
      vertexCoords[1] = coords[0];
      vertexCoords[2] = coords[1];
      this->integrate( vertexCoords, matrixRow );
      elMat( 2, 0 ) = matrixRow[1];
      elMat( 2, 1 ) = matrixRow[2];
      elMat( 2, 2 ) = matrixRow[0];
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const override
   {
      WALBERLA_ABORT( "integrateAll() for 3D not implemented by current HyTeG form." );
   }

   /// Transitional routine to allow 2D HyTeG forms inplace of FEniCS forms until we clean up the interfaces
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override {
     Matrix3r elMat;
     this->integrateAll( coords, elMat );
     out[0] = elMat( 0, 0 );
     out[1] = elMat( 0, 1 );
     out[2] = elMat( 0, 2 );
   }

   /// Transitional routine to allow 3D HyTeG forms inplace of FEniCS forms until we clean up the interfaces
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override {
     Matrix4r elMat;
     this->integrateAll( coords, elMat );
     out[0] = elMat( 0, 0 );
     out[1] = elMat( 0, 1 );
     out[2] = elMat( 0, 2 );
     out[3] = elMat( 0, 3 );
   }

   // We'd need to implement that in each child as we partially have separate 2D and 3D forms for P1 elements
   // at the moment; although the P1ElementwiseOperator does not make use of these anyway
   bool assemble2D() const override
   {
      WALBERLA_ABORT( "Don't call assemble2D on a P1FormHyteG child" );
      return false;
   };
   bool assemble3D() const override
   {
      WALBERLA_ABORT( "Don't call assemble3D on a P1FormHyteG child" );
      return false;
   };
   bool assembly2DDefined() const override
   {
      WALBERLA_ABORT( "Don't call assembly2DDefined on a P1FormHyteG child" );
      return false;
   };
   bool assembly3DDefined() const override
   {
      WALBERLA_ABORT( "Don't call assembly3DDefined on a P1FormHyteG child" );
      return false;
   };

};

} // namespace hyteg
