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
// #include "hyteg/forms/P1Form.hpp"

namespace hyteg {

class P1FormHyTeG
{
 public:
   virtual ~P1FormHyTeG() {}

   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const
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

   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const { WALBERLA_ABORT( "integrateAll() for 3D not implemented by current HyTeG form." ); }

   void setGeometryMap( const std::shared_ptr< GeometryMap > map ) { this->geometryMap_ = map; }

 protected:
   std::shared_ptr< GeometryMap > geometryMap_;
};

} // namespace hyteg
