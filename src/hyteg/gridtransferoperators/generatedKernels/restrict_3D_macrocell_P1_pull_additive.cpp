/*
 * Copyright (c) 2019 Nils Kohl, Dominik Thoennes.
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

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "restrict_3D_macrocell_P1_pull_additive.hpp"

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

template < typename ValueType >
static void restrict_3D_macrocell_P1_pull_additive_level_any( ValueType* RESTRICT             _data_vertexCoarseDst,
                                                              ValueType const* RESTRICT const _data_vertexFineSrc,
                                                              int                             coarse_level,
                                                              ValueType                       num_neighbor_cells_edge0,
                                                              ValueType                       num_neighbor_cells_edge1,
                                                              ValueType                       num_neighbor_cells_edge2,
                                                              ValueType                       num_neighbor_cells_edge3,
                                                              ValueType                       num_neighbor_cells_edge4,
                                                              ValueType                       num_neighbor_cells_edge5,
                                                              ValueType                       num_neighbor_cells_face0,
                                                              ValueType                       num_neighbor_cells_face1,
                                                              ValueType                       num_neighbor_cells_face2,
                                                              ValueType                       num_neighbor_cells_face3,
                                                              ValueType                       num_neighbor_cells_vertex0,
                                                              ValueType                       num_neighbor_cells_vertex1,
                                                              ValueType                       num_neighbor_cells_vertex2,
                                                              ValueType                       num_neighbor_cells_vertex3 )
{
   const ValueType xi_21  = 1 / ( num_neighbor_cells_edge0 );
   const ValueType xi_22  = 1 / ( num_neighbor_cells_edge1 );
   const ValueType xi_23  = 1 / ( num_neighbor_cells_edge3 );
   const ValueType xi_24  = 1 / ( num_neighbor_cells_vertex0 );
   const ValueType xi_70  = 1 / ( num_neighbor_cells_face0 );
   const ValueType xi_71  = 1 / ( num_neighbor_cells_face1 );
   const ValueType xi_34  = 1 / ( num_neighbor_cells_edge2 );
   const ValueType xi_35  = 1 / ( num_neighbor_cells_edge4 );
   const ValueType xi_36  = 1 / ( num_neighbor_cells_vertex1 );
   const ValueType xi_85  = 1 / ( num_neighbor_cells_face2 );
   const ValueType xi_100 = 1 / ( num_neighbor_cells_face3 );
   const ValueType xi_47  = 1 / ( num_neighbor_cells_edge5 );
   const ValueType xi_48  = 1 / ( num_neighbor_cells_vertex2 );
   const ValueType xi_60  = 1 / ( num_neighbor_cells_vertex3 );
   {
      for ( int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1 )
      {
         for ( int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1 )
         {
            // vertex 0
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_26 =
                   xi_21 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_27 =
                   xi_22 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_28 =
                   xi_23 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_29 =
                   static_cast< ValueType >( 1.0 ) * xi_24 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_26 + xi_27 + xi_28 + xi_29;
            }
            // edge 0
            for ( int ctr_1 = 1; ctr_1 < ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_73 =
                   xi_21 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_74 =
                   xi_21 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_79 =
                   static_cast< ValueType >( 1.0 ) * xi_21 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_75 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_76 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_77 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_78 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79;
            }
            // vertex 1
            for ( int ctr_1 = ( 1 << ( coarse_level ) ); ctr_1 < ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
            {
               const ValueType xi_38 =
                   xi_21 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_39 =
                   xi_34 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_40 =
                   xi_35 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_41 =
                   static_cast< ValueType >( 1.0 ) * xi_36 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_38 + xi_39 + xi_40 + xi_41;
            }
         }
         for ( int ctr_2 = 1; ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 += 1 )
         {
            // edge 1
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_87 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_88 =
                   xi_22 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_89 =
                   xi_22 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_94 =
                   static_cast< ValueType >( 1.0 ) * xi_22 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_90 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_91 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_92 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_93 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94;
            }
            // face 0
            for ( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_157 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_158 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_159 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_160 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_161 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_162 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_163 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_164 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_165 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_166 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_167 =
                   static_cast< ValueType >( 1.0 ) * xi_70 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] =
                   xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162 + xi_163 + xi_164 + xi_165 + xi_166 + xi_167;
            }
            // edge 2
            for ( int ctr_1 = -ctr_2 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_2 + ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
            {
               const ValueType xi_104 =
                   xi_34 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_105 =
                   xi_34 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_108 =
                   static_cast< ValueType >( 1.0 ) * xi_34 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_106 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_107 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_102 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_103 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_102 + xi_103 + xi_104 + xi_105 + xi_106 + xi_107 + xi_108;
            }
         }
         for ( int ctr_2 = -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ) + 1; ctr_2 += 1 )
         {
            // vertex 2
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_50 =
                   xi_22 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_51 =
                   xi_34 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_52 =
                   xi_47 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_53 =
                   static_cast< ValueType >( 1.0 ) * xi_48 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_50 + xi_51 + xi_52 + xi_53;
            }
         }
      }
      for ( int ctr_3 = 1; ctr_3 < ( 1 << ( coarse_level ) ); ctr_3 += 1 )
      {
         for ( int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1 )
         {
            // edge 3
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_116 =
                   xi_23 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_117 =
                   xi_23 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_122 =
                   static_cast< ValueType >( 1.0 ) * xi_23 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_118 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_119 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_120 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_121 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
            }
            // face 1
            for ( int ctr_1 = 1; ctr_1 < -ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_173 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_174 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_175 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_176 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_177 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_178 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_179 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_180 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_181 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_182 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_183 =
                   static_cast< ValueType >( 1.0 ) * xi_71 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] =
                   xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183;
            }
            // edge 4
            for ( int ctr_1 = -ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_3 + ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
            {
               const ValueType xi_130 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_131 =
                   xi_35 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_132 =
                   xi_35 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_137 =
                   static_cast< ValueType >( 1.0 ) * xi_35 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_133 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_134 =
                   xi_71 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_135 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_136 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_130 + xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137;
            }
         }
         for ( int ctr_2 = 1; ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 += 1 )
         {
            // face 2
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_189 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_190 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_191 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_192 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_193 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_194 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_195 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_196 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_197 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_198 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_199 =
                   static_cast< ValueType >( 1.0 ) * xi_85 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] =
                   xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199;
            }
            // cell (inner)
            for ( int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_3 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_4 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_5 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_6 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_7 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_8 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_9 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_10 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_11 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_12 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_13 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_14 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_15 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_16 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_17 =
                   static_cast< ValueType >( 1.0 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] =
                   xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_17 + xi_3 + xi_4 + xi_5 + xi_6 + xi_7 + xi_8 + xi_9;
            }
            // face 3
            for ( int ctr_1 = -ctr_2 - ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_2 - ctr_3 + ( 1 << ( coarse_level ) ) + 1;
                  ctr_1 += 1 )
            {
               const ValueType xi_205 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_206 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_207 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_208 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_209 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_210 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_211 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_212 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_213 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_214 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_215 =
                   static_cast< ValueType >( 1.0 ) * xi_100 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] =
                   xi_205 + xi_206 + xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
            }
         }
         for ( int ctr_2 = -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ) + 1; ctr_2 += 1 )
         {
            // edge 5
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_145 =
                   xi_47 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_146 =
                   xi_47 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_151 =
                   static_cast< ValueType >( 1.0 ) * xi_47 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_147 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_148 =
                   xi_85 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_149 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_150 =
                   xi_100 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_145 + xi_146 + xi_147 + xi_148 + xi_149 + xi_150 + xi_151;
            }
         }
      }
      for ( int ctr_3 = ( 1 << ( coarse_level ) ); ctr_3 < ( 1 << ( coarse_level ) ) + 1; ctr_3 += 1 )
      {
         for ( int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1 )
         {
            // vertex 3
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_62 =
                   xi_23 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_63 =
                   xi_35 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_64 =
                   xi_47 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_65 =
                   static_cast< ValueType >( 1.0 ) * xi_60 *
                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexCoarseDst[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                     ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                     ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) ) -
                                     ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                         ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                       ( 6 ) )] = xi_62 + xi_63 + xi_64 + xi_65;
            }
         }
      }
   }
}

template < typename ValueType >
void restrict_3D_macrocell_P1_pull_additive( ValueType* RESTRICT             _data_vertexCoarseDst,
                                             ValueType const* RESTRICT const _data_vertexFineSrc,
                                             int                             coarse_level,
                                             ValueType                       num_neighbor_cells_edge0,
                                             ValueType                       num_neighbor_cells_edge1,
                                             ValueType                       num_neighbor_cells_edge2,
                                             ValueType                       num_neighbor_cells_edge3,
                                             ValueType                       num_neighbor_cells_edge4,
                                             ValueType                       num_neighbor_cells_edge5,
                                             ValueType                       num_neighbor_cells_face0,
                                             ValueType                       num_neighbor_cells_face1,
                                             ValueType                       num_neighbor_cells_face2,
                                             ValueType                       num_neighbor_cells_face3,
                                             ValueType                       num_neighbor_cells_vertex0,
                                             ValueType                       num_neighbor_cells_vertex1,
                                             ValueType                       num_neighbor_cells_vertex2,
                                             ValueType                       num_neighbor_cells_vertex3 )
{
   switch ( coarse_level )
   {
   default:
      restrict_3D_macrocell_P1_pull_additive_level_any( _data_vertexCoarseDst,
                                                        _data_vertexFineSrc,
                                                        coarse_level,
                                                        num_neighbor_cells_edge0,
                                                        num_neighbor_cells_edge1,
                                                        num_neighbor_cells_edge2,
                                                        num_neighbor_cells_edge3,
                                                        num_neighbor_cells_edge4,
                                                        num_neighbor_cells_edge5,
                                                        num_neighbor_cells_face0,
                                                        num_neighbor_cells_face1,
                                                        num_neighbor_cells_face2,
                                                        num_neighbor_cells_face3,
                                                        num_neighbor_cells_vertex0,
                                                        num_neighbor_cells_vertex1,
                                                        num_neighbor_cells_vertex2,
                                                        num_neighbor_cells_vertex3 );
      break;
   }
}

// ========================
//  explicit instantiation
// ========================
template void
    restrict_3D_macrocell_P1_pull_additive< walberla::float64 >( walberla::float64* RESTRICT             _data_vertexCoarseDst,
                                                                 walberla::float64 const* RESTRICT const _data_vertexFineSrc,
                                                                 int                                     coarse_level,
                                                                 walberla::float64                       num_neighbor_cells_edge0,
                                                                 walberla::float64                       num_neighbor_cells_edge1,
                                                                 walberla::float64                       num_neighbor_cells_edge2,
                                                                 walberla::float64                       num_neighbor_cells_edge3,
                                                                 walberla::float64                       num_neighbor_cells_edge4,
                                                                 walberla::float64                       num_neighbor_cells_edge5,
                                                                 walberla::float64                       num_neighbor_cells_face0,
                                                                 walberla::float64                       num_neighbor_cells_face1,
                                                                 walberla::float64                       num_neighbor_cells_face2,
                                                                 walberla::float64                       num_neighbor_cells_face3,
                                                                 walberla::float64 num_neighbor_cells_vertex0,
                                                                 walberla::float64 num_neighbor_cells_vertex1,
                                                                 walberla::float64 num_neighbor_cells_vertex2,
                                                                 walberla::float64 num_neighbor_cells_vertex3 );
template void
    restrict_3D_macrocell_P1_pull_additive< walberla::float32 >( walberla::float32* RESTRICT             _data_vertexCoarseDst,
                                                                 walberla::float32 const* RESTRICT const _data_vertexFineSrc,
                                                                 int                                     coarse_level,
                                                                 walberla::float32                       num_neighbor_cells_edge0,
                                                                 walberla::float32                       num_neighbor_cells_edge1,
                                                                 walberla::float32                       num_neighbor_cells_edge2,
                                                                 walberla::float32                       num_neighbor_cells_edge3,
                                                                 walberla::float32                       num_neighbor_cells_edge4,
                                                                 walberla::float32                       num_neighbor_cells_edge5,
                                                                 walberla::float32                       num_neighbor_cells_face0,
                                                                 walberla::float32                       num_neighbor_cells_face1,
                                                                 walberla::float32                       num_neighbor_cells_face2,
                                                                 walberla::float32                       num_neighbor_cells_face3,
                                                                 walberla::float32 num_neighbor_cells_vertex0,
                                                                 walberla::float32 num_neighbor_cells_vertex1,
                                                                 walberla::float32 num_neighbor_cells_vertex2,
                                                                 walberla::float32 num_neighbor_cells_vertex3 );
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void
    restrict_3D_macrocell_P1_pull_additive< walberla::float16 >( walberla::float16* RESTRICT             _data_vertexCoarseDst,
                                                                 walberla::float16 const* RESTRICT const _data_vertexFineSrc,
                                                                 int                                     coarse_level,
                                                                 walberla::float16                       num_neighbor_cells_edge0,
                                                                 walberla::float16                       num_neighbor_cells_edge1,
                                                                 walberla::float16                       num_neighbor_cells_edge2,
                                                                 walberla::float16                       num_neighbor_cells_edge3,
                                                                 walberla::float16                       num_neighbor_cells_edge4,
                                                                 walberla::float16                       num_neighbor_cells_edge5,
                                                                 walberla::float16                       num_neighbor_cells_face0,
                                                                 walberla::float16                       num_neighbor_cells_face1,
                                                                 walberla::float16                       num_neighbor_cells_face2,
                                                                 walberla::float16                       num_neighbor_cells_face3,
                                                                 walberla::float16 num_neighbor_cells_vertex0,
                                                                 walberla::float16 num_neighbor_cells_vertex1,
                                                                 walberla::float16 num_neighbor_cells_vertex2,
                                                                 walberla::float16 num_neighbor_cells_vertex3 );
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg