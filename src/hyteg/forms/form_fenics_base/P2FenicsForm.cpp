/*
* Copyright (c) 2017-2022 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "P2FenicsForm.hpp"

namespace hyteg {

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
{
   Matrix6r localStiffnessMatrix;
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 0, 0 );
   out[1] = localStiffnessMatrix( 0, 1 );
   out[2] = localStiffnessMatrix( 0, 2 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateEdgeToVertex( const std::array< Point3D, 3 >& coords,
                                                                          Point3D&                        out ) const
{
   Matrix6r localStiffnessMatrix;
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 0, 3 );
   out[1] = localStiffnessMatrix( 0, 4 );
   out[2] = localStiffnessMatrix( 0, 5 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateVertexToEdge( const std::array< Point3D, 3 >& coords,
                                                                          Point3D&                        out ) const
{
   Matrix6r localStiffnessMatrix;
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 5, 0 );
   out[1] = localStiffnessMatrix( 5, 1 );
   out[2] = localStiffnessMatrix( 5, 2 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateEdgeToEdge( const std::array< Point3D, 3 >& coords,
                                                                        Point3D&                        out ) const
{
   Matrix6r localStiffnessMatrix;
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 5, 3 );
   out[1] = localStiffnessMatrix( 5, 4 );
   out[2] = localStiffnessMatrix( 5, 5 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const
{
   computeLocalStiffnessMatrix( coords, elMat );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
{
   Matrix10r localStiffnessMatrix;
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 0, 0 );
   out[1] = localStiffnessMatrix( 0, 1 );
   out[2] = localStiffnessMatrix( 0, 2 );
   out[3] = localStiffnessMatrix( 0, 3 );
}

template < class UFCOperator2D, class UFCOperator3D >
real_t P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrate( const std::array< Point3D, 4 >&     coords,
                                                                const P2Form::dofPosByVertexPair3D& cntrPos,
                                                                const P2Form::dofPosByVertexPair3D& leafPos ) const
{
   Matrix10r elMat;
   computeLocalStiffnessMatrix( coords, elMat );
   uint_t rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
   uint_t colIdx = fenics::P2DoFMap[leafPos[0]][leafPos[1]];

   return walberla::real_c( elMat( rowIdx, colIdx ) );
}

template < class UFCOperator2D, class UFCOperator3D >
std::vector< real_t >
    P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrate( const std::array< Point3D, 4 >&                    coords,
                                                             const P2Form::dofPosByVertexPair3D&                cntrPos,
                                                             const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
{
   Matrix10r elMat;
   computeLocalStiffnessMatrix( coords, elMat );
   std::vector< real_t > matRow( leafPos.size(), 0 );

   uint_t rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
   uint_t colIdx = 0;

   for ( uint_t k = 0; k < leafPos.size(); ++k )
   {
      colIdx    = fenics::P2DoFMap[leafPos[k][0]][leafPos[k][1]];
      matRow[k] = walberla::real_c( elMat( rowIdx, colIdx ) );
   }

   return matRow;
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const
{
   computeLocalStiffnessMatrix( coords, elMat );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords,
                                                                                Matrix6r& localStiffnessMatrix ) const
{
   real_t fenicsCoords[6];
   fenicsCoords[0] = coords[0][0];
   fenicsCoords[1] = coords[0][1];
   fenicsCoords[2] = coords[1][0];
   fenicsCoords[3] = coords[1][1];
   fenicsCoords[4] = coords[2][0];
   fenicsCoords[5] = coords[2][1];
   UFCOperator2D gen;
   gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords,
                                                                                Matrix10r& localStiffnessMatrix ) const
{
   real_t fenicsCoords[12];
   for ( int node = 0; node < 4; ++node )
   {
      for ( int dim = 0; dim < 3; ++dim )
      {
         fenicsCoords[node * 3 + dim] = coords[node][dim];
      }
   }
   UFCOperator3D gen;
   gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
}

} // namespace hyteg