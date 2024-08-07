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

#include "prolongate_3D_macrocell_P1_push_additive.hpp"

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

template < typename ValueType >
static void prolongate_3D_macrocell_P1_push_additive_level_any( ValueType const* RESTRICT const _data_vertexCoarseSrc,
                                                                ValueType* RESTRICT             _data_vertexFineDst,
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
   const ValueType xi_36  = 1 / ( num_neighbor_cells_edge0 );
   const ValueType xi_37  = 1 / ( num_neighbor_cells_edge1 );
   const ValueType xi_38  = 1 / ( num_neighbor_cells_edge3 );
   const ValueType xi_39  = 1 / ( num_neighbor_cells_vertex0 );
   const ValueType xi_101 = 1 / ( num_neighbor_cells_face0 );
   const ValueType xi_102 = 1 / ( num_neighbor_cells_face1 );
   const ValueType xi_53  = 1 / ( num_neighbor_cells_edge2 );
   const ValueType xi_54  = 1 / ( num_neighbor_cells_edge4 );
   const ValueType xi_55  = 1 / ( num_neighbor_cells_vertex1 );
   const ValueType xi_123 = 1 / ( num_neighbor_cells_face2 );
   const ValueType xi_146 = 1 / ( num_neighbor_cells_face3 );
   const ValueType xi_70  = 1 / ( num_neighbor_cells_edge5 );
   const ValueType xi_71  = 1 / ( num_neighbor_cells_vertex2 );
   const ValueType xi_87  = 1 / ( num_neighbor_cells_vertex3 );
   {
      for ( int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1 )
      {
         for ( int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1 )
         {
            // vertex 0
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_42 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_44 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_46 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_48 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_43 =
                   xi_36 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_41 =
                   xi_37 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_45 =
                   xi_38 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_47 =
                   static_cast< ValueType >( 1.0 ) * xi_39 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_41 + xi_42;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_43 + xi_44;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_45 + xi_46;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_47 + xi_48;
            }
            // edge 0
            for ( int ctr_1 = 1; ctr_1 < ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_105 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
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
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_109 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_111 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_113 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_115 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_117 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_112 =
                   xi_36 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_116 =
                   static_cast< ValueType >( 1.0 ) * xi_36 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_110 =
                   xi_101 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_114 =
                   xi_102 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_105 + xi_112;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_107 + xi_110;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_109 + xi_114;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_110 + xi_111;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_112 + xi_113;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_114 + xi_115;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_116 + xi_117;
            }
            // vertex 1
            for ( int ctr_1 = ( 1 << ( coarse_level ) ); ctr_1 < ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
            {
               const ValueType xi_58 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_60 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_62 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_64 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_57 =
                   xi_36 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_59 =
                   xi_53 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_61 =
                   xi_54 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_63 =
                   static_cast< ValueType >( 1.0 ) * xi_55 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_57 + xi_58;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_59 + xi_60;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_61 + xi_62;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_63 + xi_64;
            }
         }
         for ( int ctr_2 = 1; ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 += 1 )
         {
            // edge 1
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_126 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_128 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_130 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_132 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_134 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_136 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_137 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_138 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_140 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_131 =
                   xi_37 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_139 =
                   static_cast< ValueType >( 1.0 ) * xi_37 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_133 =
                   xi_101 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_135 =
                   xi_123 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_126 + xi_131;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_128 + xi_133;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_130 + xi_135;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_131 + xi_132;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_133 + xi_134;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_135 + xi_136;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_137 + xi_138;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_139 + xi_140;
            }
            // face 0
            for ( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_233 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_235 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_250 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_237 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_239 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_241 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_243 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_245 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_247 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_249 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_251 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_253 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_246 =
                   xi_101 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_252 =
                   static_cast< ValueType >( 1.0 ) * xi_101 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_233 + xi_246;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_235 + xi_246;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_237 + xi_250;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_239 + xi_246;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_241 + xi_246;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_243 + xi_250;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_245 + xi_246;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_246 + xi_247;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_249 + xi_250;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_250 + xi_251;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_252 + xi_253;
            }
            // edge 2
            for ( int ctr_1 = -ctr_2 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_2 + ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
            {
               const ValueType xi_149 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_151 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_153 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_155 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_157 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_159 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_161 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_156 =
                   xi_53 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_160 =
                   static_cast< ValueType >( 1.0 ) * xi_53 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_154 =
                   xi_101 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_158 =
                   xi_146 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_149 + xi_154;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_151 + xi_156;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_153 + xi_158;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_154 + xi_155;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_156 + xi_157;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_158 + xi_159;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_160 + xi_161;
            }
         }
         for ( int ctr_2 = -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ) + 1; ctr_2 += 1 )
         {
            // vertex 2
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_74 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_76 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_78 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_80 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_73 =
                   xi_37 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_75 =
                   xi_53 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_77 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_79 =
                   static_cast< ValueType >( 1.0 ) * xi_71 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_73 + xi_74;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_75 + xi_76;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_77 + xi_78;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_79 + xi_80;
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
               const ValueType xi_170 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_172 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_174 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_176 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_178 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_180 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_182 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_179 =
                   xi_38 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_181 =
                   static_cast< ValueType >( 1.0 ) * xi_38 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_177 =
                   xi_102 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_175 =
                   xi_123 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_170 + xi_179;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_172 + xi_177;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_174 + xi_175;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_175 + xi_176;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_177 + xi_178;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_179 + xi_180;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_181 + xi_182;
            }
            // face 1
            for ( int ctr_1 = 1; ctr_1 < -ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_260 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_262 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_277 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_264 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_266 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_268 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_270 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_272 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_274 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_276 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_278 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_280 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_275 =
                   xi_102 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_279 =
                   static_cast< ValueType >( 1.0 ) * xi_102 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_260 + xi_275;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_262 + xi_275;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_264 + xi_277;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_266 + xi_275;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_268 + xi_277;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_270 + xi_275;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_272 + xi_277;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_274 + xi_275;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_275 + xi_276;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_277 + xi_278;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_279 + xi_280;
            }
            // edge 4
            for ( int ctr_1 = -ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_3 + ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
            {
               const ValueType xi_191 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_193 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_195 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_197 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_199 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_201 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_202 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_203 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_205 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_200 =
                   xi_54 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_204 =
                   static_cast< ValueType >( 1.0 ) * xi_54 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_196 =
                   xi_102 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_198 =
                   xi_146 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_191 + xi_196;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_193 + xi_200;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_195 + xi_198;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_196 + xi_197;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_198 + xi_199;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_200 + xi_201;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_202 + xi_203;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_204 + xi_205;
            }
         }
         for ( int ctr_2 = 1; ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 += 1 )
         {
            // face 2
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_287 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_304 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_289 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_291 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_293 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_295 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_297 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_299 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_301 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_303 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_305 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_307 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_302 =
                   xi_123 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_306 =
                   static_cast< ValueType >( 1.0 ) * xi_123 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_287 + xi_302;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_289 + xi_304;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_291 + xi_302;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_293 + xi_302;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_295 + xi_304;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_297 + xi_302;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_299 + xi_302;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_301 + xi_304;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_302 + xi_303;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_304 + xi_305;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_306 + xi_307;
            }
            // cell (inner)
            for ( int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
            {
               const ValueType xi_29 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_4 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_6 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_8 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_10 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_12 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
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
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_16 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_18 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_20 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_22 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_24 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_26 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_28 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_30 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_31 =
                   static_cast< ValueType >( 1.0 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_32 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_29 + xi_4;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_29 + xi_6;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_29 + xi_8;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_10 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_12 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_14 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_16 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_18 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_20 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_22 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_24 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_26 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_28 + xi_29;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_29 + xi_30;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_31 + xi_32;
            }
            // face 3
            for ( int ctr_1 = -ctr_2 - ctr_3 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_2 - ctr_3 + ( 1 << ( coarse_level ) ) + 1;
                  ctr_1 += 1 )
            {
               const ValueType xi_331 =
                   static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_314 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_316 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_318 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_320 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_322 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_324 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_326 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_328 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_330 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_332 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) -
                                       1];
               const ValueType xi_334 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_329 =
                   xi_146 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_333 =
                   static_cast< ValueType >( 1.0 ) * xi_146 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_314 + xi_331;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_316 + xi_329;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_318 + xi_329;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_320 + xi_331;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_322 + xi_329;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_324 + xi_329;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_326 + xi_331;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_328 + xi_329;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_329 + xi_330;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) -
                                   1]         = xi_331 + xi_332;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_333 + xi_334;
            }
         }
         for ( int ctr_2 = -ctr_3 + ( 1 << ( coarse_level ) ); ctr_2 < -ctr_3 + ( 1 << ( coarse_level ) ) + 1; ctr_2 += 1 )
         {
            // edge 5
            for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
            {
               const ValueType xi_214 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_216 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_218 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_220 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_222 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_224 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                         ( 6 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_226 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_223 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_225 =
                   static_cast< ValueType >( 1.0 ) * xi_70 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_219 =
                   xi_123 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_221 =
                   xi_146 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_214 + xi_219;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_216 + xi_221;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_218 + xi_223;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_219 + xi_220;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_221 + xi_222;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) ) /
                                     ( 6 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_223 + xi_224;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_225 + xi_226;
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
               const ValueType xi_90 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_92 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) ) +
                                       1];
               const ValueType xi_94 =
                   _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                         ( 6 ) )];
               const ValueType xi_96 =
                   _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                       ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) ) -
                                       ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                           ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                         ( 6 ) )];
               const ValueType xi_89 =
                   xi_38 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_91 =
                   xi_54 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_93 =
                   xi_70 * static_cast< ValueType >( 0.5 ) *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               const ValueType xi_95 =
                   static_cast< ValueType >( 1.0 ) * xi_87 *
                   _data_vertexCoarseSrc[ctr_1 + ctr_2 * ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) -
                                         ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) ) +
                                         ( ( ( ( 1 << ( coarse_level ) ) + 1 ) * ( ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) ) -
                                         ( ( ( -ctr_3 + ( 1 << ( coarse_level ) ) + 1 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 2 ) *
                                             ( -ctr_3 + ( 1 << ( coarse_level ) ) + 3 ) ) /
                                           ( 6 ) )];
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_89 + xi_90;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) ) +
                                   1]         = xi_91 + xi_92;
               _data_vertexFineDst[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) -
                                   ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 4 ) ) /
                                     ( 6 ) )] = xi_93 + xi_94;
               _data_vertexFineDst[2 * ctr_1 + 2 * ctr_2 * ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                   ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) +
                                   ( ( ( ( 1 << ( coarse_level + 1 ) ) + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) ) -
                                   ( ( ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 1 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 2 ) *
                                       ( -2 * ctr_3 + ( 1 << ( coarse_level + 1 ) ) + 3 ) ) /
                                     ( 6 ) )] = xi_95 + xi_96;
            }
         }
      }
   }
}

template < typename ValueType >
void prolongate_3D_macrocell_P1_push_additive( ValueType const* RESTRICT const _data_vertexCoarseSrc,
                                               ValueType* RESTRICT             _data_vertexFineDst,
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
      prolongate_3D_macrocell_P1_push_additive_level_any( _data_vertexCoarseSrc,
                                                          _data_vertexFineDst,
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
    prolongate_3D_macrocell_P1_push_additive< walberla::float64 >( walberla::float64 const* RESTRICT const _data_vertexCoarseSrc,
                                                                   walberla::float64* RESTRICT             _data_vertexFineDst,
                                                                   int                                     coarse_level,
                                                                   walberla::float64 num_neighbor_cells_edge0,
                                                                   walberla::float64 num_neighbor_cells_edge1,
                                                                   walberla::float64 num_neighbor_cells_edge2,
                                                                   walberla::float64 num_neighbor_cells_edge3,
                                                                   walberla::float64 num_neighbor_cells_edge4,
                                                                   walberla::float64 num_neighbor_cells_edge5,
                                                                   walberla::float64 num_neighbor_cells_face0,
                                                                   walberla::float64 num_neighbor_cells_face1,
                                                                   walberla::float64 num_neighbor_cells_face2,
                                                                   walberla::float64 num_neighbor_cells_face3,
                                                                   walberla::float64 num_neighbor_cells_vertex0,
                                                                   walberla::float64 num_neighbor_cells_vertex1,
                                                                   walberla::float64 num_neighbor_cells_vertex2,
                                                                   walberla::float64 num_neighbor_cells_vertex3 );
template void
    prolongate_3D_macrocell_P1_push_additive< walberla::float32 >( walberla::float32 const* RESTRICT const _data_vertexCoarseSrc,
                                                                   walberla::float32* RESTRICT             _data_vertexFineDst,
                                                                   int                                     coarse_level,
                                                                   walberla::float32 num_neighbor_cells_edge0,
                                                                   walberla::float32 num_neighbor_cells_edge1,
                                                                   walberla::float32 num_neighbor_cells_edge2,
                                                                   walberla::float32 num_neighbor_cells_edge3,
                                                                   walberla::float32 num_neighbor_cells_edge4,
                                                                   walberla::float32 num_neighbor_cells_edge5,
                                                                   walberla::float32 num_neighbor_cells_face0,
                                                                   walberla::float32 num_neighbor_cells_face1,
                                                                   walberla::float32 num_neighbor_cells_face2,
                                                                   walberla::float32 num_neighbor_cells_face3,
                                                                   walberla::float32 num_neighbor_cells_vertex0,
                                                                   walberla::float32 num_neighbor_cells_vertex1,
                                                                   walberla::float32 num_neighbor_cells_vertex2,
                                                                   walberla::float32 num_neighbor_cells_vertex3 );
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void
    prolongate_3D_macrocell_P1_push_additive< walberla::float16 >( walberla::float16 const* RESTRICT const _data_vertexCoarseSrc,
                                                                   walberla::float16* RESTRICT             _data_vertexFineDst,
                                                                   int                                     coarse_level,
                                                                   walberla::float16 num_neighbor_cells_edge0,
                                                                   walberla::float16 num_neighbor_cells_edge1,
                                                                   walberla::float16 num_neighbor_cells_edge2,
                                                                   walberla::float16 num_neighbor_cells_edge3,
                                                                   walberla::float16 num_neighbor_cells_edge4,
                                                                   walberla::float16 num_neighbor_cells_edge5,
                                                                   walberla::float16 num_neighbor_cells_face0,
                                                                   walberla::float16 num_neighbor_cells_face1,
                                                                   walberla::float16 num_neighbor_cells_face2,
                                                                   walberla::float16 num_neighbor_cells_face3,
                                                                   walberla::float16 num_neighbor_cells_vertex0,
                                                                   walberla::float16 num_neighbor_cells_vertex1,
                                                                   walberla::float16 num_neighbor_cells_vertex2,
                                                                   walberla::float16 num_neighbor_cells_vertex3 );
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg