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

#include "hyteg/forms/form_fenics_generated/p2_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_div.h"
#include "hyteg/forms/form_fenics_generated/p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p2_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_pspg.h"
#include "hyteg/forms/form_fenics_generated/p2_stokes_epsilon.h"
#include "hyteg/forms/form_fenics_generated/p2_stokes_full.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_div_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_divt_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_pspg_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_stokes_epsilon_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_stokes_full_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_polar_laplacian.h"

namespace hyteg {

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
{
   Matrix6r localStiffnessMatrix{ Matrix6r::Zero() };
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 0, 0 );
   out[1] = localStiffnessMatrix( 0, 1 );
   out[2] = localStiffnessMatrix( 0, 2 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateEdgeToVertex( const std::array< Point3D, 3 >& coords,
                                                                          Point3D&                        out ) const
{
   Matrix6r localStiffnessMatrix{ Matrix6r::Zero() };
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 0, 3 );
   out[1] = localStiffnessMatrix( 0, 4 );
   out[2] = localStiffnessMatrix( 0, 5 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateVertexToEdge( const std::array< Point3D, 3 >& coords,
                                                                          Point3D&                        out ) const
{
   Matrix6r localStiffnessMatrix{ Matrix6r::Zero() };
   computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
   out[0] = localStiffnessMatrix( 5, 0 );
   out[1] = localStiffnessMatrix( 5, 1 );
   out[2] = localStiffnessMatrix( 5, 2 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrateEdgeToEdge( const std::array< Point3D, 3 >& coords,
                                                                        Point3D&                        out ) const
{
   Matrix6r localStiffnessMatrix{ Matrix6r::Zero() };
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
   Matrix10r localStiffnessMatrix{ Matrix10r::Zero() };
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
   Matrix10r elMat{ Matrix10r::Zero() };
   computeLocalStiffnessMatrix( coords, elMat );
   int rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
   int colIdx = fenics::P2DoFMap[leafPos[0]][leafPos[1]];

   return walberla::real_c( elMat( rowIdx, colIdx ) );
}

template < class UFCOperator2D, class UFCOperator3D >
std::vector< real_t >
    P2FenicsForm< UFCOperator2D, UFCOperator3D >::integrate( const std::array< Point3D, 4 >&                    coords,
                                                             const P2Form::dofPosByVertexPair3D&                cntrPos,
                                                             const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
{
   Matrix10r elMat{ Matrix10r::Zero() };
   computeLocalStiffnessMatrix( coords, elMat );
   std::vector< real_t > matRow( leafPos.size(), 0 );

   int rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
   int colIdx = 0;

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
   double fenicsCoords[6];
   fenicsCoords[0] = coords[0][0];
   fenicsCoords[1] = coords[0][1];
   fenicsCoords[2] = coords[1][0];
   fenicsCoords[3] = coords[1][1];
   fenicsCoords[4] = coords[2][0];
   fenicsCoords[5] = coords[2][1];
   UFCOperator2D gen;
   Eigen::Matrix< double, 6, 6, Eigen::RowMajor > tmp = localStiffnessMatrix.cast<double>();
   gen.tabulate_tensor( tmp.data(), nullptr, fenicsCoords, 0 );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2FenicsForm< UFCOperator2D, UFCOperator3D >::computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords,
                                                                                Matrix10r& localStiffnessMatrix ) const
{
   double fenicsCoords[12];
   for ( int node = 0; node < 4; ++node )
   {
      for ( int dim = 0; dim < 3; ++dim )
      {
         fenicsCoords[node * 3 + dim] = coords[walberla::uint_c( node )][walberla::uint_c( dim )];
      }
   }
   UFCOperator3D gen;
   Eigen::Matrix< double, 10, 10, Eigen::RowMajor > tmp = localStiffnessMatrix.cast<double>();
   gen.tabulate_tensor( tmp.data(), nullptr, fenicsCoords, 0 );
}

template class P2FenicsForm<hyteg::fenics::NoAssemble, hyteg::fenics::NoAssemble>;

template class P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >;
template class P2FenicsForm< p2_mass_cell_integral_0_otherwise,      p2_tet_mass_cell_integral_0_otherwise >;
template class P2FenicsForm< p2_divt_cell_integral_0_otherwise,      p2_tet_divt_tet_cell_integral_0_otherwise >;
template class P2FenicsForm< p2_divt_cell_integral_1_otherwise,      p2_tet_divt_tet_cell_integral_1_otherwise >;
template class P2FenicsForm< fenics::NoAssemble,                     p2_tet_divt_tet_cell_integral_2_otherwise >;
template class P2FenicsForm< p2_div_cell_integral_0_otherwise,       p2_tet_div_tet_cell_integral_0_otherwise >;
template class P2FenicsForm< p2_div_cell_integral_1_otherwise,       p2_tet_div_tet_cell_integral_1_otherwise >;
template class P2FenicsForm< fenics::NoAssemble,                     p2_tet_div_tet_cell_integral_2_otherwise >;
template class P2FenicsForm< p2_pspg_cell_integral_0_otherwise,      p2_tet_pspg_tet_cell_integral_0_otherwise >;

template class P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise >;
template class P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise >;
template class P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise >;
template class P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise >;

template class P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise >;
template class P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise >;
template class P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise >;
template class P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise >;
template class P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise >;

template class P2FenicsForm< p2_polar_laplacian_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >;
template class P2FenicsForm<p2_mass_cell_integral_0_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_diffusion_cell_integral_0_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise >;
template class P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise >;
template class P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise >;
template class P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise >;
template class P2FenicsForm<p2_stokes_full_cell_integral_0_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_full_cell_integral_1_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_full_cell_integral_2_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_full_cell_integral_3_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_epsilon_cell_integral_0_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_epsilon_cell_integral_1_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_epsilon_cell_integral_2_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<p2_stokes_epsilon_cell_integral_3_otherwise, hyteg::fenics::NoAssemble>;
template class P2FenicsForm<hyteg::fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise>;
template class P2FenicsForm<hyteg::fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_1_otherwise>;
template class P2FenicsForm<hyteg::fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_3_otherwise>;
template class P2FenicsForm<hyteg::fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_4_otherwise>;
template class P2FenicsForm<hyteg::fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise>;
template class P2FenicsForm<hyteg::fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_0_otherwise>;

} // namespace hyteg