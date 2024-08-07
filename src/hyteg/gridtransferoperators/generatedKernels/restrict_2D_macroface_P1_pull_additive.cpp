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

#include "restrict_2D_macroface_P1_pull_additive.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

template < typename ValueType >
static void restrict_2D_macroface_P1_pull_additive_level_any( ValueType* RESTRICT             _data_vertexCoarseDst,
                                                              ValueType const* RESTRICT const _data_vertexFineSrc,
                                                              int                             coarse_level,
                                                              ValueType                       num_neighbor_faces_edge0,
                                                              ValueType                       num_neighbor_faces_edge1,
                                                              ValueType                       num_neighbor_faces_edge2,
                                                              ValueType                       num_neighbor_faces_vertex0,
                                                              ValueType                       num_neighbor_faces_vertex1,
                                                              ValueType                       num_neighbor_faces_vertex2 )
{
   const ValueType xi_13 = 1 / ( num_neighbor_faces_edge0 );
   const ValueType xi_14 = 1 / ( num_neighbor_faces_edge1 );
   const ValueType xi_15 = 1 / ( num_neighbor_faces_vertex0 );
   const ValueType xi_24 = 1 / ( num_neighbor_faces_edge2 );
   const ValueType xi_25 = 1 / ( num_neighbor_faces_vertex1 );
   const ValueType xi_35 = 1 / ( num_neighbor_faces_vertex2 );
   {
      for ( int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1 )
      {
         // bottom left vertex
         for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
         {
            const ValueType xi_17 = xi_13 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_18 = xi_14 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) )];
            const ValueType xi_19 = static_cast< ValueType >( 1.0 ) * xi_15 *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_17 + xi_18 + xi_19;
         }
         // bottom edge
         for ( int ctr_1 = 1; ctr_1 < ( 1 << ( coarse_level ) ); ctr_1 += 1 )
         {
            const ValueType xi_49 = static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) - 1];
            const ValueType xi_50 = static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) )];
            const ValueType xi_51 = xi_13 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) - 1];
            const ValueType xi_52 = xi_13 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_53 = static_cast< ValueType >( 1.0 ) * xi_13 *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_49 + xi_50 + xi_51 + xi_52 + xi_53;
         }
         // bottom right vertex
         for ( int ctr_1 = ( 1 << ( coarse_level ) ); ctr_1 < ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
         {
            const ValueType xi_27 = xi_13 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) - 1];
            const ValueType xi_28 = xi_24 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) - 1];
            const ValueType xi_29 = static_cast< ValueType >( 1.0 ) * xi_25 *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_27 + xi_28 + xi_29;
         }
      }
      for ( int ctr_2 = 1; ctr_2 < ( 1 << ( coarse_level ) ); ctr_2 += 1 )
      {
         // left edge
         for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
         {
            const ValueType xi_59 = static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_60 = static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_61 = xi_14 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) )];
            const ValueType xi_62 = xi_14 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) )];
            const ValueType xi_63 = static_cast< ValueType >( 1.0 ) * xi_14 *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_59 + xi_60 + xi_61 + xi_62 + xi_63;
         }
         // inner triangle
         for ( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << ( coarse_level ) ); ctr_1 += 1 )
         {
            const ValueType xi_3 = static_cast< ValueType >( 0.5 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) - 1];
            const ValueType xi_4 = static_cast< ValueType >( 0.5 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) - 1];
            const ValueType xi_5 = static_cast< ValueType >( 0.5 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) )];
            const ValueType xi_6 = static_cast< ValueType >( 0.5 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) )];
            const ValueType xi_7 = static_cast< ValueType >( 0.5 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_8 = static_cast< ValueType >( 0.5 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_9 = static_cast< ValueType >( 1.0 ) *
                                   _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                       ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_3 + xi_4 + xi_5 + xi_6 + xi_7 + xi_8 + xi_9;
         }
         // diagonal edge
         for ( int ctr_1 = -ctr_2 + ( 1 << ( coarse_level ) ); ctr_1 < -ctr_2 + ( 1 << ( coarse_level ) ) + 1; ctr_1 += 1 )
         {
            const ValueType xi_69 = static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) ) - 1];
            const ValueType xi_70 = static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) )];
            const ValueType xi_71 = xi_24 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 + 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( ( 2 * ctr_2 + 1 ) * ( 2 * ctr_2 + 2 ) ) / ( 2 ) ) - 1];
            const ValueType xi_72 = xi_24 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_73 = static_cast< ValueType >( 1.0 ) * xi_24 *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_69 + xi_70 + xi_71 + xi_72 + xi_73;
         }
      }
      for ( int ctr_2 = ( 1 << ( coarse_level ) ); ctr_2 < ( 1 << ( coarse_level ) ) + 1; ctr_2 += 1 )
      {
         // top vertex
         for ( int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1 )
         {
            const ValueType xi_37 = xi_14 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) )];
            const ValueType xi_38 = xi_24 * static_cast< ValueType >( 0.5 ) *
                                    _data_vertexFineSrc[2 * ctr_1 + ( 2 * ctr_2 - 1 ) * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 - 1 ) ) / ( 2 ) ) + 1];
            const ValueType xi_39 = static_cast< ValueType >( 1.0 ) * xi_35 *
                                    _data_vertexFineSrc[2 * ctr_1 + 2 * ctr_2 * ( ( 1 << ( coarse_level + 1 ) ) + 2 ) -
                                                        ( ( 2 * ctr_2 * ( 2 * ctr_2 + 1 ) ) / ( 2 ) )];
            _data_vertexCoarseDst[ctr_1 + ctr_2 * ( ( 1 << ( coarse_level ) ) + 2 ) - ( ( ctr_2 * ( ctr_2 + 1 ) ) / ( 2 ) )] =
                xi_37 + xi_38 + xi_39;
         }
      }
   }
}

template < typename ValueType >
void restrict_2D_macroface_P1_pull_additive( ValueType* RESTRICT             _data_vertexCoarseDst,
                                             ValueType const* RESTRICT const _data_vertexFineSrc,
                                             int                             coarse_level,
                                             ValueType                       num_neighbor_faces_edge0,
                                             ValueType                       num_neighbor_faces_edge1,
                                             ValueType                       num_neighbor_faces_edge2,
                                             ValueType                       num_neighbor_faces_vertex0,
                                             ValueType                       num_neighbor_faces_vertex1,
                                             ValueType                       num_neighbor_faces_vertex2 )
{
   switch ( coarse_level )
   {
   default:
      restrict_2D_macroface_P1_pull_additive_level_any( _data_vertexCoarseDst,
                                                        _data_vertexFineSrc,
                                                        coarse_level,
                                                        num_neighbor_faces_edge0,
                                                        num_neighbor_faces_edge1,
                                                        num_neighbor_faces_edge2,
                                                        num_neighbor_faces_vertex0,
                                                        num_neighbor_faces_vertex1,
                                                        num_neighbor_faces_vertex2 );
      break;
   }
}

// ========================
//  explicit instantiation
// ========================
template void
    restrict_2D_macroface_P1_pull_additive< walberla::float64 >( walberla::float64* RESTRICT             _data_vertexCoarseDst,
                                                                 walberla::float64 const* RESTRICT const _data_vertexFineSrc,
                                                                 int                                     coarse_level,
                                                                 walberla::float64                       num_neighbor_faces_edge0,
                                                                 walberla::float64                       num_neighbor_faces_edge1,
                                                                 walberla::float64                       num_neighbor_faces_edge2,
                                                                 walberla::float64 num_neighbor_faces_vertex0,
                                                                 walberla::float64 num_neighbor_faces_vertex1,
                                                                 walberla::float64 num_neighbor_faces_vertex2 );
template void
    restrict_2D_macroface_P1_pull_additive< walberla::float32 >( walberla::float32* RESTRICT             _data_vertexCoarseDst,
                                                                 walberla::float32 const* RESTRICT const _data_vertexFineSrc,
                                                                 int                                     coarse_level,
                                                                 walberla::float32                       num_neighbor_faces_edge0,
                                                                 walberla::float32                       num_neighbor_faces_edge1,
                                                                 walberla::float32                       num_neighbor_faces_edge2,
                                                                 walberla::float32 num_neighbor_faces_vertex0,
                                                                 walberla::float32 num_neighbor_faces_vertex1,
                                                                 walberla::float32 num_neighbor_faces_vertex2 );
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void
    restrict_2D_macroface_P1_pull_additive< walberla::float16 >( walberla::float16* RESTRICT             _data_vertexCoarseDst,
                                                                 walberla::float16 const* RESTRICT const _data_vertexFineSrc,
                                                                 int                                     coarse_level,
                                                                 walberla::float16                       num_neighbor_faces_edge0,
                                                                 walberla::float16                       num_neighbor_faces_edge1,
                                                                 walberla::float16                       num_neighbor_faces_edge2,
                                                                 walberla::float16 num_neighbor_faces_vertex0,
                                                                 walberla::float16 num_neighbor_faces_vertex1,
                                                                 walberla::float16 num_neighbor_faces_vertex2 );
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg