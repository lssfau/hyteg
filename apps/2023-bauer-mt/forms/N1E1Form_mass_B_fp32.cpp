/*
 * Copyright (c) 2023 Daniel Bauer.
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

#include "N1E1Form_mass_B_fp32.hpp"

namespace hyteg {
namespace n1e1 {

void N1E1Form_mass_B_fp32::integrateAll( const std::array< PointND< float, 3 >, 4 >& coords,
                                         const std::array< walberla::int16_t, 6 >&   edgeDirections,
                                         Matrix< float, 6, 6 >&                      elMat ) const
{
   float p_affine_0_0 = coords[0][0];
   float p_affine_0_1 = coords[0][1];
   float p_affine_0_2 = coords[0][2];
   float p_affine_1_0 = coords[1][0];
   float p_affine_1_1 = coords[1][1];
   float p_affine_1_2 = coords[1][2];
   float p_affine_2_0 = coords[2][0];
   float p_affine_2_1 = coords[2][1];
   float p_affine_2_2 = coords[2][2];
   float p_affine_3_0 = coords[3][0];
   float p_affine_3_1 = coords[3][1];
   float p_affine_3_2 = coords[3][2];
   float tmp_0        = -p_affine_0_0 + p_affine_1_0;
   float tmp_1        = -p_affine_0_1 + p_affine_2_1;
   float tmp_2        = tmp_0 * tmp_1;
   float tmp_3        = -p_affine_0_0 + p_affine_2_0;
   float tmp_4        = -p_affine_0_1 + p_affine_1_1;
   float tmp_5        = tmp_3 * tmp_4;
   float tmp_6        = tmp_2 - tmp_5;
   float tmp_7        = ( tmp_6 * tmp_6 );
   float tmp_8        = -p_affine_0_2 + p_affine_3_2;
   float tmp_9        = -p_affine_0_1 + p_affine_3_1;
   float tmp_10       = -p_affine_0_2 + p_affine_1_2;
   float tmp_11       = -p_affine_0_0 + p_affine_3_0;
   float tmp_12       = -p_affine_0_2 + p_affine_2_2;
   float tmp_13       = tmp_12 * tmp_4;
   float tmp_14       = tmp_0 * tmp_9;
   float tmp_15       = tmp_10 * tmp_11;
   float tmp_16 =
       1.0f /
       ( ( -tmp_1 * tmp_15 + tmp_10 * tmp_3 * tmp_9 + tmp_11 * tmp_13 - tmp_12 * tmp_14 + tmp_2 * tmp_8 - tmp_5 * tmp_8 ) *
         ( -tmp_1 * tmp_15 + tmp_10 * tmp_3 * tmp_9 + tmp_11 * tmp_13 - tmp_12 * tmp_14 + tmp_2 * tmp_8 - tmp_5 * tmp_8 ) );
   float tmp_17 = p_affine_0_0 * p_affine_1_1;
   float tmp_18 = p_affine_0_0 * p_affine_1_2;
   float tmp_19 = p_affine_2_1 * p_affine_3_2;
   float tmp_20 = p_affine_0_1 * p_affine_1_0;
   float tmp_21 = p_affine_0_1 * p_affine_1_2;
   float tmp_22 = p_affine_2_2 * p_affine_3_0;
   float tmp_23 = p_affine_0_2 * p_affine_1_0;
   float tmp_24 = p_affine_0_2 * p_affine_1_1;
   float tmp_25 = p_affine_2_0 * p_affine_3_1;
   float tmp_26 = p_affine_2_2 * p_affine_3_1;
   float tmp_27 = p_affine_2_0 * p_affine_3_2;
   float tmp_28 = p_affine_2_1 * p_affine_3_0;
   float tmp_29 = std::abs( p_affine_0_0 * tmp_19 - p_affine_0_0 * tmp_26 + p_affine_0_1 * tmp_22 - p_affine_0_1 * tmp_27 +
                            p_affine_0_2 * tmp_25 - p_affine_0_2 * tmp_28 - p_affine_1_0 * tmp_19 + p_affine_1_0 * tmp_26 -
                            p_affine_1_1 * tmp_22 + p_affine_1_1 * tmp_27 - p_affine_1_2 * tmp_25 + p_affine_1_2 * tmp_28 +
                            p_affine_2_0 * tmp_21 - p_affine_2_0 * tmp_24 - p_affine_2_1 * tmp_18 + p_affine_2_1 * tmp_23 +
                            p_affine_2_2 * tmp_17 - p_affine_2_2 * tmp_20 - p_affine_3_0 * tmp_21 + p_affine_3_0 * tmp_24 +
                            p_affine_3_1 * tmp_18 - p_affine_3_1 * tmp_23 - p_affine_3_2 * tmp_17 + p_affine_3_2 * tmp_20 );
   float tmp_30 = tmp_16 * tmp_29;
   float tmp_31 = ( 1.0f / 60.0f ) * tmp_30;
   float tmp_32 = static_cast< float >( edgeDirections[0] * edgeDirections[0] ) * tmp_31;
   float tmp_33 = tmp_11 * tmp_4 - tmp_14;
   float tmp_34 = ( tmp_33 * tmp_33 );
   float tmp_35 = -tmp_0 * tmp_12 + tmp_10 * tmp_3;
   float tmp_36 = ( tmp_35 * tmp_35 );
   float tmp_37 = tmp_0 * tmp_8 - tmp_15;
   float tmp_38 = ( tmp_37 * tmp_37 );
   float tmp_39 = -tmp_1 * tmp_10 + tmp_13;
   float tmp_40 = ( tmp_39 * tmp_39 );
   float tmp_41 = tmp_10 * tmp_9 - tmp_4 * tmp_8;
   float tmp_42 = ( tmp_41 * tmp_41 );
   float tmp_43 = tmp_33 * tmp_6;
   float tmp_44 = tmp_35 * tmp_37;
   float tmp_45 = tmp_39 * tmp_41;
   float tmp_46 = static_cast< float >( edgeDirections[0] ) * tmp_7;
   float tmp_47 = ( 1.0f / 120.0f ) * tmp_30;
   float tmp_48 = static_cast< float >( edgeDirections[1] ) * tmp_47;
   float tmp_49 = static_cast< float >( edgeDirections[0] ) * tmp_48;
   float tmp_50 = -tmp_1 * tmp_11 + tmp_3 * tmp_9;
   float tmp_51 = tmp_50 * tmp_6;
   float tmp_52 = tmp_33 * tmp_50;
   float tmp_53 = static_cast< float >( edgeDirections[0] ) * tmp_31;
   float tmp_54 = static_cast< float >( edgeDirections[1] ) * tmp_53;
   float tmp_55 = tmp_11 * tmp_12 - tmp_3 * tmp_8;
   float tmp_56 = tmp_35 * tmp_55;
   float tmp_57 = tmp_37 * tmp_55;
   float tmp_58 = static_cast< float >( edgeDirections[1] ) * tmp_57;
   float tmp_59 = tmp_1 * tmp_8 - tmp_12 * tmp_9;
   float tmp_60 = tmp_39 * tmp_59;
   float tmp_61 = tmp_41 * tmp_59;
   float tmp_62 = tmp_36 * tmp_49 + tmp_40 * tmp_49 - tmp_43 * tmp_49 - tmp_44 * tmp_49 - tmp_45 * tmp_49 + tmp_46 * tmp_48 -
                  tmp_49 * tmp_51 - tmp_49 * tmp_56 - tmp_49 * tmp_60 + tmp_52 * tmp_54 + tmp_53 * tmp_58 + tmp_54 * tmp_61;
   float tmp_63 = static_cast< float >( edgeDirections[2] ) * tmp_47;
   float tmp_64 = static_cast< float >( edgeDirections[0] ) * tmp_63;
   float tmp_65 = static_cast< float >( edgeDirections[0] ) * tmp_38;
   float tmp_66 = static_cast< float >( edgeDirections[2] ) * tmp_53;
   float tmp_67 = -tmp_34 * tmp_64 - tmp_42 * tmp_64 + tmp_43 * tmp_64 + tmp_44 * tmp_64 + tmp_45 * tmp_64 - tmp_51 * tmp_66 +
                  tmp_52 * tmp_64 - tmp_56 * tmp_66 + tmp_57 * tmp_64 - tmp_60 * tmp_66 + tmp_61 * tmp_64 - tmp_63 * tmp_65;
   float tmp_68 = static_cast< float >( edgeDirections[3] ) * tmp_31;
   float tmp_69 = static_cast< float >( edgeDirections[0] ) * tmp_68;
   float tmp_70 =
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[0] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_35 * tmp_55 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[0] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_36 +
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[0] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_39 * tmp_59 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[0] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_40 +
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[0] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_50 * tmp_6 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[0] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_7 -
       tmp_34 * tmp_69 - tmp_42 * tmp_69 - tmp_43 * tmp_69 - tmp_44 * tmp_69 - tmp_45 * tmp_69 - tmp_52 * tmp_69 -
       tmp_57 * tmp_69 - tmp_61 * tmp_69 - tmp_65 * tmp_68;
   float tmp_71 = static_cast< float >( edgeDirections[4] ) * tmp_31;
   float tmp_72 = static_cast< float >( edgeDirections[0] ) * tmp_71;
   float tmp_73 = static_cast< float >( edgeDirections[4] ) * tmp_47;
   float tmp_74 = static_cast< float >( edgeDirections[0] ) * tmp_73;
   float tmp_75 = static_cast< float >( edgeDirections[0] ) * tmp_57;
   float tmp_76 = -tmp_34 * tmp_72 + tmp_36 * tmp_72 + tmp_40 * tmp_72 - tmp_42 * tmp_72 + tmp_43 * tmp_72 + tmp_44 * tmp_72 +
                  tmp_45 * tmp_72 + tmp_46 * tmp_71 + tmp_51 * tmp_72 - tmp_52 * tmp_74 + tmp_56 * tmp_72 + tmp_60 * tmp_72 -
                  tmp_61 * tmp_74 - tmp_65 * tmp_71 - tmp_73 * tmp_75;
   float tmp_77 = static_cast< float >( edgeDirections[5] ) * tmp_47;
   float tmp_78 = static_cast< float >( edgeDirections[0] ) * tmp_77;
   float tmp_79 = static_cast< float >( edgeDirections[5] ) * tmp_31;
   float tmp_80 = static_cast< float >( edgeDirections[0] ) * tmp_79;
   float tmp_81 = -tmp_34 * tmp_78 + tmp_36 * tmp_78 + tmp_40 * tmp_78 - tmp_42 * tmp_78 + tmp_46 * tmp_77 + tmp_51 * tmp_80 -
                  tmp_52 * tmp_80 + tmp_56 * tmp_80 + tmp_60 * tmp_80 - tmp_61 * tmp_80 - tmp_65 * tmp_77 - tmp_75 * tmp_79;
   float tmp_82 = static_cast< float >( edgeDirections[1] * edgeDirections[1] ) * tmp_31;
   float tmp_83 = ( tmp_50 * tmp_50 );
   float tmp_84 = ( tmp_55 * tmp_55 );
   float tmp_85 = ( tmp_59 * tmp_59 );
   float tmp_86 = static_cast< float >( edgeDirections[2] ) * tmp_48;
   float tmp_87 = static_cast< float >( edgeDirections[1] ) * tmp_43;
   float tmp_88 = static_cast< float >( edgeDirections[2] ) * tmp_31;
   float tmp_89 = static_cast< float >( edgeDirections[1] ) * tmp_88;
   float tmp_90 = tmp_44 * tmp_89 + tmp_45 * tmp_89 - tmp_51 * tmp_86 - tmp_52 * tmp_86 - tmp_56 * tmp_86 - tmp_57 * tmp_86 -
                  tmp_60 * tmp_86 - tmp_61 * tmp_86 + tmp_83 * tmp_86 + tmp_84 * tmp_86 + tmp_85 * tmp_86 + tmp_87 * tmp_88;
   float tmp_91 = static_cast< float >( edgeDirections[1] ) * tmp_68;
   float tmp_92 =
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[1] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_33 * tmp_6 +
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[1] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_35 * tmp_37 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[1] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_36 +
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[1] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_39 * tmp_41 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[1] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_40 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[1] * edgeDirections[3] ) * tmp_16 * tmp_29 * tmp_7 -
       tmp_51 * tmp_91 - tmp_52 * tmp_91 - tmp_56 * tmp_91 - tmp_57 * tmp_91 - tmp_60 * tmp_91 - tmp_61 * tmp_91 -
       tmp_83 * tmp_91 - tmp_84 * tmp_91 - tmp_85 * tmp_91;
   float tmp_93 = static_cast< float >( edgeDirections[4] ) * tmp_48;
   float tmp_94 = static_cast< float >( edgeDirections[1] ) * tmp_71;
   float tmp_95 = tmp_36 * tmp_93 + tmp_40 * tmp_93 + tmp_44 * tmp_94 + tmp_45 * tmp_94 - tmp_52 * tmp_94 - tmp_58 * tmp_71 -
                  tmp_61 * tmp_94 + tmp_7 * tmp_93 + tmp_71 * tmp_87 - tmp_83 * tmp_93 - tmp_84 * tmp_93 - tmp_85 * tmp_93;
   float tmp_96 = static_cast< float >( edgeDirections[1] ) * tmp_79;
   float tmp_97 = static_cast< float >( edgeDirections[5] ) * tmp_48;
   float tmp_98 = tmp_36 * tmp_96 + tmp_40 * tmp_96 + tmp_43 * tmp_96 + tmp_44 * tmp_96 + tmp_45 * tmp_96 + tmp_51 * tmp_96 -
                  tmp_52 * tmp_97 + tmp_56 * tmp_96 - tmp_57 * tmp_97 + tmp_60 * tmp_96 - tmp_61 * tmp_97 + tmp_7 * tmp_96 -
                  tmp_83 * tmp_96 - tmp_84 * tmp_96 - tmp_85 * tmp_96;
   float tmp_99  = static_cast< float >( edgeDirections[2] * edgeDirections[2] ) * tmp_31;
   float tmp_100 = static_cast< float >( edgeDirections[3] ) * tmp_63;
   float tmp_101 = static_cast< float >( edgeDirections[2] ) * tmp_68;
   float tmp_102 = tmp_100 * tmp_34 + tmp_100 * tmp_38 + tmp_100 * tmp_42 - tmp_100 * tmp_83 - tmp_100 * tmp_84 -
                   tmp_100 * tmp_85 + tmp_101 * tmp_43 + tmp_101 * tmp_44 + tmp_101 * tmp_45 - tmp_101 * tmp_51 -
                   tmp_101 * tmp_56 - tmp_101 * tmp_60;
   float tmp_103 = static_cast< float >( edgeDirections[2] ) * tmp_71;
   float tmp_104 =
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[2] * edgeDirections[4] ) * tmp_16 * tmp_29 * tmp_33 * tmp_6 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[2] * edgeDirections[4] ) * tmp_16 * tmp_29 * tmp_34 +
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[2] * edgeDirections[4] ) * tmp_16 * tmp_29 * tmp_35 * tmp_37 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[2] * edgeDirections[4] ) * tmp_16 * tmp_29 * tmp_38 +
       ( 1.0f / 120.0f ) * static_cast< float >( edgeDirections[2] * edgeDirections[4] ) * tmp_16 * tmp_29 * tmp_39 * tmp_41 +
       ( 1.0f / 60.0f ) * static_cast< float >( edgeDirections[2] * edgeDirections[4] ) * tmp_16 * tmp_29 * tmp_42 -
       tmp_103 * tmp_51 - tmp_103 * tmp_52 - tmp_103 * tmp_56 - tmp_103 * tmp_57 - tmp_103 * tmp_60 - tmp_103 * tmp_61 -
       tmp_103 * tmp_83 - tmp_103 * tmp_84 - tmp_103 * tmp_85;
   float tmp_105 = static_cast< float >( edgeDirections[2] ) * tmp_79;
   float tmp_106 = static_cast< float >( edgeDirections[5] ) * tmp_63;
   float tmp_107 = tmp_105 * tmp_34 + tmp_105 * tmp_38 + tmp_105 * tmp_42 + tmp_105 * tmp_43 + tmp_105 * tmp_44 +
                   tmp_105 * tmp_45 + tmp_105 * tmp_52 + tmp_105 * tmp_57 + tmp_105 * tmp_61 - tmp_105 * tmp_83 -
                   tmp_105 * tmp_84 - tmp_105 * tmp_85 - tmp_106 * tmp_51 - tmp_106 * tmp_56 - tmp_106 * tmp_60;
   float tmp_108 = static_cast< float >( edgeDirections[3] * edgeDirections[3] );
   float tmp_109 = ( 1.0f / 20.0f ) * tmp_30;
   float tmp_110 = tmp_108 * tmp_109;
   float tmp_111 = tmp_108 * tmp_31;
   float tmp_112 = ( 1.0f / 30.0f ) * tmp_30;
   float tmp_113 = tmp_108 * tmp_112;
   float tmp_114 = static_cast< float >( edgeDirections[4] ) * tmp_68;
   float tmp_115 = static_cast< float >( edgeDirections[3] * edgeDirections[4] ) * tmp_47;
   float tmp_116 = static_cast< float >( edgeDirections[3] * edgeDirections[4] );
   float tmp_117 = tmp_109 * tmp_116;
   float tmp_118 = ( 1.0f / 40.0f ) * tmp_30;
   float tmp_119 = tmp_116 * tmp_118;
   float tmp_120 = tmp_114 * tmp_34 + tmp_114 * tmp_36 + tmp_114 * tmp_38 + tmp_114 * tmp_40 + tmp_114 * tmp_42 +
                   tmp_114 * tmp_7 + tmp_115 * tmp_83 + tmp_115 * tmp_84 + tmp_115 * tmp_85 + tmp_117 * tmp_43 +
                   tmp_117 * tmp_44 + tmp_117 * tmp_45 + tmp_119 * tmp_51 + tmp_119 * tmp_52 + tmp_119 * tmp_56 +
                   tmp_119 * tmp_57 + tmp_119 * tmp_60 + tmp_119 * tmp_61;
   float tmp_121 = static_cast< float >( edgeDirections[5] ) * tmp_68;
   float tmp_122 = static_cast< float >( edgeDirections[3] ) * tmp_77;
   float tmp_123 = static_cast< float >( edgeDirections[3] * edgeDirections[5] );
   float tmp_124 = tmp_118 * tmp_123;
   float tmp_125 = tmp_109 * tmp_123;
   float tmp_126 = tmp_121 * tmp_36 + tmp_121 * tmp_40 + tmp_121 * tmp_7 + tmp_121 * tmp_83 + tmp_121 * tmp_84 +
                   tmp_121 * tmp_85 + tmp_122 * tmp_34 + tmp_122 * tmp_38 + tmp_122 * tmp_42 + tmp_124 * tmp_43 +
                   tmp_124 * tmp_44 + tmp_124 * tmp_45 + tmp_124 * tmp_52 + tmp_124 * tmp_57 + tmp_124 * tmp_61 +
                   tmp_125 * tmp_51 + tmp_125 * tmp_56 + tmp_125 * tmp_60;
   float tmp_127 = static_cast< float >( edgeDirections[4] * edgeDirections[4] );
   float tmp_128 = tmp_127 * tmp_31;
   float tmp_129 = tmp_109 * tmp_127;
   float tmp_130 = tmp_112 * tmp_127;
   float tmp_131 = static_cast< float >( edgeDirections[4] ) * tmp_77;
   float tmp_132 = static_cast< float >( edgeDirections[5] ) * tmp_71;
   float tmp_133 = static_cast< float >( edgeDirections[4] * edgeDirections[5] );
   float tmp_134 = tmp_118 * tmp_133;
   float tmp_135 = tmp_109 * tmp_133;
   float tmp_136 = tmp_131 * tmp_36 + tmp_131 * tmp_40 + tmp_131 * tmp_7 + tmp_132 * tmp_34 + tmp_132 * tmp_38 +
                   tmp_132 * tmp_42 + tmp_132 * tmp_83 + tmp_132 * tmp_84 + tmp_132 * tmp_85 + tmp_134 * tmp_43 +
                   tmp_134 * tmp_44 + tmp_134 * tmp_45 + tmp_134 * tmp_51 + tmp_134 * tmp_56 + tmp_134 * tmp_60 +
                   tmp_135 * tmp_52 + tmp_135 * tmp_57 + tmp_135 * tmp_61;
   float tmp_137 = static_cast< float >( edgeDirections[5] * edgeDirections[5] );
   float tmp_138 = tmp_137 * tmp_31;
   float tmp_139 = tmp_109 * tmp_137;
   float tmp_140 = tmp_112 * tmp_137;
   float a_0_0   = tmp_32 * tmp_34 + tmp_32 * tmp_36 + tmp_32 * tmp_38 + tmp_32 * tmp_40 + tmp_32 * tmp_42 - tmp_32 * tmp_43 -
                 tmp_32 * tmp_44 - tmp_32 * tmp_45 + tmp_32 * tmp_7;
   float a_0_1 = tmp_62;
   float a_0_2 = tmp_67;
   float a_0_3 = tmp_70;
   float a_0_4 = tmp_76;
   float a_0_5 = tmp_81;
   float a_1_0 = tmp_62;
   float a_1_1 = tmp_36 * tmp_82 + tmp_40 * tmp_82 - tmp_51 * tmp_82 - tmp_56 * tmp_82 - tmp_60 * tmp_82 + tmp_7 * tmp_82 +
                 tmp_82 * tmp_83 + tmp_82 * tmp_84 + tmp_82 * tmp_85;
   float a_1_2 = tmp_90;
   float a_1_3 = tmp_92;
   float a_1_4 = tmp_95;
   float a_1_5 = tmp_98;
   float a_2_0 = tmp_67;
   float a_2_1 = tmp_90;
   float a_2_2 = tmp_34 * tmp_99 + tmp_38 * tmp_99 + tmp_42 * tmp_99 - tmp_52 * tmp_99 - tmp_57 * tmp_99 - tmp_61 * tmp_99 +
                 tmp_83 * tmp_99 + tmp_84 * tmp_99 + tmp_85 * tmp_99;
   float a_2_3 = tmp_102;
   float a_2_4 = tmp_104;
   float a_2_5 = tmp_107;
   float a_3_0 = tmp_70;
   float a_3_1 = tmp_92;
   float a_3_2 = tmp_102;
   float a_3_3 = tmp_110 * tmp_36 + tmp_110 * tmp_40 + tmp_110 * tmp_43 + tmp_110 * tmp_44 + tmp_110 * tmp_45 + tmp_110 * tmp_51 +
                 tmp_110 * tmp_56 + tmp_110 * tmp_60 + tmp_110 * tmp_7 + tmp_111 * tmp_34 + tmp_111 * tmp_38 + tmp_111 * tmp_42 +
                 tmp_111 * tmp_83 + tmp_111 * tmp_84 + tmp_111 * tmp_85 + tmp_113 * tmp_52 + tmp_113 * tmp_57 + tmp_113 * tmp_61;
   float a_3_4 = tmp_120;
   float a_3_5 = tmp_126;
   float a_4_0 = tmp_76;
   float a_4_1 = tmp_95;
   float a_4_2 = tmp_104;
   float a_4_3 = tmp_120;
   float a_4_4 = tmp_128 * tmp_36 + tmp_128 * tmp_40 + tmp_128 * tmp_7 + tmp_128 * tmp_83 + tmp_128 * tmp_84 + tmp_128 * tmp_85 +
                 tmp_129 * tmp_34 + tmp_129 * tmp_38 + tmp_129 * tmp_42 + tmp_129 * tmp_43 + tmp_129 * tmp_44 + tmp_129 * tmp_45 +
                 tmp_129 * tmp_52 + tmp_129 * tmp_57 + tmp_129 * tmp_61 + tmp_130 * tmp_51 + tmp_130 * tmp_56 + tmp_130 * tmp_60;
   float a_4_5 = tmp_136;
   float a_5_0 = tmp_81;
   float a_5_1 = tmp_98;
   float a_5_2 = tmp_107;
   float a_5_3 = tmp_126;
   float a_5_4 = tmp_136;
   float a_5_5 = tmp_138 * tmp_34 + tmp_138 * tmp_36 + tmp_138 * tmp_38 + tmp_138 * tmp_40 + tmp_138 * tmp_42 + tmp_138 * tmp_7 +
                 tmp_139 * tmp_51 + tmp_139 * tmp_52 + tmp_139 * tmp_56 + tmp_139 * tmp_57 + tmp_139 * tmp_60 + tmp_139 * tmp_61 +
                 tmp_139 * tmp_83 + tmp_139 * tmp_84 + tmp_139 * tmp_85 + tmp_140 * tmp_43 + tmp_140 * tmp_44 + tmp_140 * tmp_45;
   elMat( 0, 0 ) = a_0_0;
   elMat( 0, 1 ) = a_0_1;
   elMat( 0, 2 ) = a_0_2;
   elMat( 0, 3 ) = a_0_3;
   elMat( 0, 4 ) = a_0_4;
   elMat( 0, 5 ) = a_0_5;
   elMat( 1, 0 ) = a_1_0;
   elMat( 1, 1 ) = a_1_1;
   elMat( 1, 2 ) = a_1_2;
   elMat( 1, 3 ) = a_1_3;
   elMat( 1, 4 ) = a_1_4;
   elMat( 1, 5 ) = a_1_5;
   elMat( 2, 0 ) = a_2_0;
   elMat( 2, 1 ) = a_2_1;
   elMat( 2, 2 ) = a_2_2;
   elMat( 2, 3 ) = a_2_3;
   elMat( 2, 4 ) = a_2_4;
   elMat( 2, 5 ) = a_2_5;
   elMat( 3, 0 ) = a_3_0;
   elMat( 3, 1 ) = a_3_1;
   elMat( 3, 2 ) = a_3_2;
   elMat( 3, 3 ) = a_3_3;
   elMat( 3, 4 ) = a_3_4;
   elMat( 3, 5 ) = a_3_5;
   elMat( 4, 0 ) = a_4_0;
   elMat( 4, 1 ) = a_4_1;
   elMat( 4, 2 ) = a_4_2;
   elMat( 4, 3 ) = a_4_3;
   elMat( 4, 4 ) = a_4_4;
   elMat( 4, 5 ) = a_4_5;
   elMat( 5, 0 ) = a_5_0;
   elMat( 5, 1 ) = a_5_1;
   elMat( 5, 2 ) = a_5_2;
   elMat( 5, 3 ) = a_5_3;
   elMat( 5, 4 ) = a_5_4;
   elMat( 5, 5 ) = a_5_5;
}

} // namespace n1e1
} // namespace hyteg
