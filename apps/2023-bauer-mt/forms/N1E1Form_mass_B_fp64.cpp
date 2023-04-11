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

#include "N1E1Form_mass_B_fp64.hpp"

namespace hyteg {
namespace n1e1 {

void N1E1Form_mass_B_fp64::integrateAll( const std::array< PointND< double, 3 >, 4 >& coords,
                                         const std::array< walberla::int16_t, 6 >&    edgeDirections,
                                         Matrix< double, 6, 6 >&                      elMat ) const
{
   double p_affine_0_0 = coords[0][0];
   double p_affine_0_1 = coords[0][1];
   double p_affine_0_2 = coords[0][2];
   double p_affine_1_0 = coords[1][0];
   double p_affine_1_1 = coords[1][1];
   double p_affine_1_2 = coords[1][2];
   double p_affine_2_0 = coords[2][0];
   double p_affine_2_1 = coords[2][1];
   double p_affine_2_2 = coords[2][2];
   double p_affine_3_0 = coords[3][0];
   double p_affine_3_1 = coords[3][1];
   double p_affine_3_2 = coords[3][2];
   double tmp_0        = -p_affine_0_0 + p_affine_1_0;
   double tmp_1        = -p_affine_0_1 + p_affine_2_1;
   double tmp_2        = tmp_0 * tmp_1;
   double tmp_3        = -p_affine_0_0 + p_affine_2_0;
   double tmp_4        = -p_affine_0_1 + p_affine_1_1;
   double tmp_5        = tmp_3 * tmp_4;
   double tmp_6        = tmp_2 - tmp_5;
   double tmp_7        = ( tmp_6 * tmp_6 );
   double tmp_8        = -p_affine_0_2 + p_affine_3_2;
   double tmp_9        = -p_affine_0_1 + p_affine_3_1;
   double tmp_10       = -p_affine_0_2 + p_affine_1_2;
   double tmp_11       = -p_affine_0_0 + p_affine_3_0;
   double tmp_12       = -p_affine_0_2 + p_affine_2_2;
   double tmp_13       = tmp_12 * tmp_4;
   double tmp_14       = tmp_0 * tmp_9;
   double tmp_15       = tmp_10 * tmp_11;
   double tmp_16 =
       1.0 / ( ( -tmp_1 * tmp_15 + tmp_10 * tmp_3 * tmp_9 + tmp_11 * tmp_13 - tmp_12 * tmp_14 + tmp_2 * tmp_8 - tmp_5 * tmp_8 ) *
               ( -tmp_1 * tmp_15 + tmp_10 * tmp_3 * tmp_9 + tmp_11 * tmp_13 - tmp_12 * tmp_14 + tmp_2 * tmp_8 - tmp_5 * tmp_8 ) );
   double tmp_17 = p_affine_0_0 * p_affine_1_1;
   double tmp_18 = p_affine_0_0 * p_affine_1_2;
   double tmp_19 = p_affine_2_1 * p_affine_3_2;
   double tmp_20 = p_affine_0_1 * p_affine_1_0;
   double tmp_21 = p_affine_0_1 * p_affine_1_2;
   double tmp_22 = p_affine_2_2 * p_affine_3_0;
   double tmp_23 = p_affine_0_2 * p_affine_1_0;
   double tmp_24 = p_affine_0_2 * p_affine_1_1;
   double tmp_25 = p_affine_2_0 * p_affine_3_1;
   double tmp_26 = p_affine_2_2 * p_affine_3_1;
   double tmp_27 = p_affine_2_0 * p_affine_3_2;
   double tmp_28 = p_affine_2_1 * p_affine_3_0;
   double tmp_29 = std::abs( p_affine_0_0 * tmp_19 - p_affine_0_0 * tmp_26 + p_affine_0_1 * tmp_22 - p_affine_0_1 * tmp_27 +
                             p_affine_0_2 * tmp_25 - p_affine_0_2 * tmp_28 - p_affine_1_0 * tmp_19 + p_affine_1_0 * tmp_26 -
                             p_affine_1_1 * tmp_22 + p_affine_1_1 * tmp_27 - p_affine_1_2 * tmp_25 + p_affine_1_2 * tmp_28 +
                             p_affine_2_0 * tmp_21 - p_affine_2_0 * tmp_24 - p_affine_2_1 * tmp_18 + p_affine_2_1 * tmp_23 +
                             p_affine_2_2 * tmp_17 - p_affine_2_2 * tmp_20 - p_affine_3_0 * tmp_21 + p_affine_3_0 * tmp_24 +
                             p_affine_3_1 * tmp_18 - p_affine_3_1 * tmp_23 - p_affine_3_2 * tmp_17 + p_affine_3_2 * tmp_20 );
   double tmp_30 = tmp_16 * tmp_29;
   double tmp_31 = ( 1.0 / 60.0 ) * tmp_30;
   double tmp_32 = ( edgeDirections[0] * edgeDirections[0] ) * tmp_31;
   double tmp_33 = tmp_11 * tmp_4 - tmp_14;
   double tmp_34 = ( tmp_33 * tmp_33 );
   double tmp_35 = -tmp_0 * tmp_12 + tmp_10 * tmp_3;
   double tmp_36 = ( tmp_35 * tmp_35 );
   double tmp_37 = tmp_0 * tmp_8 - tmp_15;
   double tmp_38 = ( tmp_37 * tmp_37 );
   double tmp_39 = -tmp_1 * tmp_10 + tmp_13;
   double tmp_40 = ( tmp_39 * tmp_39 );
   double tmp_41 = tmp_10 * tmp_9 - tmp_4 * tmp_8;
   double tmp_42 = ( tmp_41 * tmp_41 );
   double tmp_43 = tmp_33 * tmp_6;
   double tmp_44 = tmp_35 * tmp_37;
   double tmp_45 = tmp_39 * tmp_41;
   double tmp_46 = edgeDirections[0] * tmp_7;
   double tmp_47 = ( 1.0 / 120.0 ) * tmp_30;
   double tmp_48 = edgeDirections[1] * tmp_47;
   double tmp_49 = edgeDirections[0] * tmp_48;
   double tmp_50 = -tmp_1 * tmp_11 + tmp_3 * tmp_9;
   double tmp_51 = tmp_50 * tmp_6;
   double tmp_52 = tmp_33 * tmp_50;
   double tmp_53 = edgeDirections[0] * tmp_31;
   double tmp_54 = edgeDirections[1] * tmp_53;
   double tmp_55 = tmp_11 * tmp_12 - tmp_3 * tmp_8;
   double tmp_56 = tmp_35 * tmp_55;
   double tmp_57 = tmp_37 * tmp_55;
   double tmp_58 = edgeDirections[1] * tmp_57;
   double tmp_59 = tmp_1 * tmp_8 - tmp_12 * tmp_9;
   double tmp_60 = tmp_39 * tmp_59;
   double tmp_61 = tmp_41 * tmp_59;
   double tmp_62 = tmp_36 * tmp_49 + tmp_40 * tmp_49 - tmp_43 * tmp_49 - tmp_44 * tmp_49 - tmp_45 * tmp_49 + tmp_46 * tmp_48 -
                   tmp_49 * tmp_51 - tmp_49 * tmp_56 - tmp_49 * tmp_60 + tmp_52 * tmp_54 + tmp_53 * tmp_58 + tmp_54 * tmp_61;
   double tmp_63 = edgeDirections[2] * tmp_47;
   double tmp_64 = edgeDirections[0] * tmp_63;
   double tmp_65 = edgeDirections[0] * tmp_38;
   double tmp_66 = edgeDirections[2] * tmp_53;
   double tmp_67 = -tmp_34 * tmp_64 - tmp_42 * tmp_64 + tmp_43 * tmp_64 + tmp_44 * tmp_64 + tmp_45 * tmp_64 - tmp_51 * tmp_66 +
                   tmp_52 * tmp_64 - tmp_56 * tmp_66 + tmp_57 * tmp_64 - tmp_60 * tmp_66 + tmp_61 * tmp_64 - tmp_63 * tmp_65;
   double tmp_68 = edgeDirections[3] * tmp_31;
   double tmp_69 = edgeDirections[0] * tmp_68;
   double tmp_70 = ( 1.0 / 120.0 ) * edgeDirections[0] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_35 * tmp_55 +
                   ( 1.0 / 60.0 ) * edgeDirections[0] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_36 +
                   ( 1.0 / 120.0 ) * edgeDirections[0] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_39 * tmp_59 +
                   ( 1.0 / 60.0 ) * edgeDirections[0] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_40 +
                   ( 1.0 / 120.0 ) * edgeDirections[0] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_50 * tmp_6 +
                   ( 1.0 / 60.0 ) * edgeDirections[0] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_7 - tmp_34 * tmp_69 -
                   tmp_42 * tmp_69 - tmp_43 * tmp_69 - tmp_44 * tmp_69 - tmp_45 * tmp_69 - tmp_52 * tmp_69 - tmp_57 * tmp_69 -
                   tmp_61 * tmp_69 - tmp_65 * tmp_68;
   double tmp_71 = edgeDirections[4] * tmp_31;
   double tmp_72 = edgeDirections[0] * tmp_71;
   double tmp_73 = edgeDirections[4] * tmp_47;
   double tmp_74 = edgeDirections[0] * tmp_73;
   double tmp_75 = edgeDirections[0] * tmp_57;
   double tmp_76 = -tmp_34 * tmp_72 + tmp_36 * tmp_72 + tmp_40 * tmp_72 - tmp_42 * tmp_72 + tmp_43 * tmp_72 + tmp_44 * tmp_72 +
                   tmp_45 * tmp_72 + tmp_46 * tmp_71 + tmp_51 * tmp_72 - tmp_52 * tmp_74 + tmp_56 * tmp_72 + tmp_60 * tmp_72 -
                   tmp_61 * tmp_74 - tmp_65 * tmp_71 - tmp_73 * tmp_75;
   double tmp_77 = edgeDirections[5] * tmp_47;
   double tmp_78 = edgeDirections[0] * tmp_77;
   double tmp_79 = edgeDirections[5] * tmp_31;
   double tmp_80 = edgeDirections[0] * tmp_79;
   double tmp_81 = -tmp_34 * tmp_78 + tmp_36 * tmp_78 + tmp_40 * tmp_78 - tmp_42 * tmp_78 + tmp_46 * tmp_77 + tmp_51 * tmp_80 -
                   tmp_52 * tmp_80 + tmp_56 * tmp_80 + tmp_60 * tmp_80 - tmp_61 * tmp_80 - tmp_65 * tmp_77 - tmp_75 * tmp_79;
   double tmp_82 = ( edgeDirections[1] * edgeDirections[1] ) * tmp_31;
   double tmp_83 = ( tmp_50 * tmp_50 );
   double tmp_84 = ( tmp_55 * tmp_55 );
   double tmp_85 = ( tmp_59 * tmp_59 );
   double tmp_86 = edgeDirections[2] * tmp_48;
   double tmp_87 = edgeDirections[1] * tmp_43;
   double tmp_88 = edgeDirections[2] * tmp_31;
   double tmp_89 = edgeDirections[1] * tmp_88;
   double tmp_90 = tmp_44 * tmp_89 + tmp_45 * tmp_89 - tmp_51 * tmp_86 - tmp_52 * tmp_86 - tmp_56 * tmp_86 - tmp_57 * tmp_86 -
                   tmp_60 * tmp_86 - tmp_61 * tmp_86 + tmp_83 * tmp_86 + tmp_84 * tmp_86 + tmp_85 * tmp_86 + tmp_87 * tmp_88;
   double tmp_91 = edgeDirections[1] * tmp_68;
   double tmp_92 = ( 1.0 / 120.0 ) * edgeDirections[1] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_33 * tmp_6 +
                   ( 1.0 / 120.0 ) * edgeDirections[1] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_35 * tmp_37 +
                   ( 1.0 / 60.0 ) * edgeDirections[1] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_36 +
                   ( 1.0 / 120.0 ) * edgeDirections[1] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_39 * tmp_41 +
                   ( 1.0 / 60.0 ) * edgeDirections[1] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_40 +
                   ( 1.0 / 60.0 ) * edgeDirections[1] * edgeDirections[3] * tmp_16 * tmp_29 * tmp_7 - tmp_51 * tmp_91 -
                   tmp_52 * tmp_91 - tmp_56 * tmp_91 - tmp_57 * tmp_91 - tmp_60 * tmp_91 - tmp_61 * tmp_91 - tmp_83 * tmp_91 -
                   tmp_84 * tmp_91 - tmp_85 * tmp_91;
   double tmp_93 = edgeDirections[4] * tmp_48;
   double tmp_94 = edgeDirections[1] * tmp_71;
   double tmp_95 = tmp_36 * tmp_93 + tmp_40 * tmp_93 + tmp_44 * tmp_94 + tmp_45 * tmp_94 - tmp_52 * tmp_94 - tmp_58 * tmp_71 -
                   tmp_61 * tmp_94 + tmp_7 * tmp_93 + tmp_71 * tmp_87 - tmp_83 * tmp_93 - tmp_84 * tmp_93 - tmp_85 * tmp_93;
   double tmp_96 = edgeDirections[1] * tmp_79;
   double tmp_97 = edgeDirections[5] * tmp_48;
   double tmp_98 = tmp_36 * tmp_96 + tmp_40 * tmp_96 + tmp_43 * tmp_96 + tmp_44 * tmp_96 + tmp_45 * tmp_96 + tmp_51 * tmp_96 -
                   tmp_52 * tmp_97 + tmp_56 * tmp_96 - tmp_57 * tmp_97 + tmp_60 * tmp_96 - tmp_61 * tmp_97 + tmp_7 * tmp_96 -
                   tmp_83 * tmp_96 - tmp_84 * tmp_96 - tmp_85 * tmp_96;
   double tmp_99  = ( edgeDirections[2] * edgeDirections[2] ) * tmp_31;
   double tmp_100 = edgeDirections[3] * tmp_63;
   double tmp_101 = edgeDirections[2] * tmp_68;
   double tmp_102 = tmp_100 * tmp_34 + tmp_100 * tmp_38 + tmp_100 * tmp_42 - tmp_100 * tmp_83 - tmp_100 * tmp_84 -
                    tmp_100 * tmp_85 + tmp_101 * tmp_43 + tmp_101 * tmp_44 + tmp_101 * tmp_45 - tmp_101 * tmp_51 -
                    tmp_101 * tmp_56 - tmp_101 * tmp_60;
   double tmp_103 = edgeDirections[2] * tmp_71;
   double tmp_104 = ( 1.0 / 120.0 ) * edgeDirections[2] * edgeDirections[4] * tmp_16 * tmp_29 * tmp_33 * tmp_6 +
                    ( 1.0 / 60.0 ) * edgeDirections[2] * edgeDirections[4] * tmp_16 * tmp_29 * tmp_34 +
                    ( 1.0 / 120.0 ) * edgeDirections[2] * edgeDirections[4] * tmp_16 * tmp_29 * tmp_35 * tmp_37 +
                    ( 1.0 / 60.0 ) * edgeDirections[2] * edgeDirections[4] * tmp_16 * tmp_29 * tmp_38 +
                    ( 1.0 / 120.0 ) * edgeDirections[2] * edgeDirections[4] * tmp_16 * tmp_29 * tmp_39 * tmp_41 +
                    ( 1.0 / 60.0 ) * edgeDirections[2] * edgeDirections[4] * tmp_16 * tmp_29 * tmp_42 - tmp_103 * tmp_51 -
                    tmp_103 * tmp_52 - tmp_103 * tmp_56 - tmp_103 * tmp_57 - tmp_103 * tmp_60 - tmp_103 * tmp_61 -
                    tmp_103 * tmp_83 - tmp_103 * tmp_84 - tmp_103 * tmp_85;
   double tmp_105 = edgeDirections[2] * tmp_79;
   double tmp_106 = edgeDirections[5] * tmp_63;
   double tmp_107 = tmp_105 * tmp_34 + tmp_105 * tmp_38 + tmp_105 * tmp_42 + tmp_105 * tmp_43 + tmp_105 * tmp_44 +
                    tmp_105 * tmp_45 + tmp_105 * tmp_52 + tmp_105 * tmp_57 + tmp_105 * tmp_61 - tmp_105 * tmp_83 -
                    tmp_105 * tmp_84 - tmp_105 * tmp_85 - tmp_106 * tmp_51 - tmp_106 * tmp_56 - tmp_106 * tmp_60;
   double tmp_108 = ( edgeDirections[3] * edgeDirections[3] );
   double tmp_109 = ( 1.0 / 20.0 ) * tmp_30;
   double tmp_110 = tmp_108 * tmp_109;
   double tmp_111 = tmp_108 * tmp_31;
   double tmp_112 = ( 1.0 / 30.0 ) * tmp_30;
   double tmp_113 = tmp_108 * tmp_112;
   double tmp_114 = edgeDirections[4] * tmp_68;
   double tmp_115 = edgeDirections[3] * edgeDirections[4] * tmp_47;
   double tmp_116 = edgeDirections[3] * edgeDirections[4];
   double tmp_117 = tmp_109 * tmp_116;
   double tmp_118 = ( 1.0 / 40.0 ) * tmp_30;
   double tmp_119 = tmp_116 * tmp_118;
   double tmp_120 = tmp_114 * tmp_34 + tmp_114 * tmp_36 + tmp_114 * tmp_38 + tmp_114 * tmp_40 + tmp_114 * tmp_42 +
                    tmp_114 * tmp_7 + tmp_115 * tmp_83 + tmp_115 * tmp_84 + tmp_115 * tmp_85 + tmp_117 * tmp_43 +
                    tmp_117 * tmp_44 + tmp_117 * tmp_45 + tmp_119 * tmp_51 + tmp_119 * tmp_52 + tmp_119 * tmp_56 +
                    tmp_119 * tmp_57 + tmp_119 * tmp_60 + tmp_119 * tmp_61;
   double tmp_121 = edgeDirections[5] * tmp_68;
   double tmp_122 = edgeDirections[3] * tmp_77;
   double tmp_123 = edgeDirections[3] * edgeDirections[5];
   double tmp_124 = tmp_118 * tmp_123;
   double tmp_125 = tmp_109 * tmp_123;
   double tmp_126 = tmp_121 * tmp_36 + tmp_121 * tmp_40 + tmp_121 * tmp_7 + tmp_121 * tmp_83 + tmp_121 * tmp_84 +
                    tmp_121 * tmp_85 + tmp_122 * tmp_34 + tmp_122 * tmp_38 + tmp_122 * tmp_42 + tmp_124 * tmp_43 +
                    tmp_124 * tmp_44 + tmp_124 * tmp_45 + tmp_124 * tmp_52 + tmp_124 * tmp_57 + tmp_124 * tmp_61 +
                    tmp_125 * tmp_51 + tmp_125 * tmp_56 + tmp_125 * tmp_60;
   double tmp_127 = ( edgeDirections[4] * edgeDirections[4] );
   double tmp_128 = tmp_127 * tmp_31;
   double tmp_129 = tmp_109 * tmp_127;
   double tmp_130 = tmp_112 * tmp_127;
   double tmp_131 = edgeDirections[4] * tmp_77;
   double tmp_132 = edgeDirections[5] * tmp_71;
   double tmp_133 = edgeDirections[4] * edgeDirections[5];
   double tmp_134 = tmp_118 * tmp_133;
   double tmp_135 = tmp_109 * tmp_133;
   double tmp_136 = tmp_131 * tmp_36 + tmp_131 * tmp_40 + tmp_131 * tmp_7 + tmp_132 * tmp_34 + tmp_132 * tmp_38 +
                    tmp_132 * tmp_42 + tmp_132 * tmp_83 + tmp_132 * tmp_84 + tmp_132 * tmp_85 + tmp_134 * tmp_43 +
                    tmp_134 * tmp_44 + tmp_134 * tmp_45 + tmp_134 * tmp_51 + tmp_134 * tmp_56 + tmp_134 * tmp_60 +
                    tmp_135 * tmp_52 + tmp_135 * tmp_57 + tmp_135 * tmp_61;
   double tmp_137 = ( edgeDirections[5] * edgeDirections[5] );
   double tmp_138 = tmp_137 * tmp_31;
   double tmp_139 = tmp_109 * tmp_137;
   double tmp_140 = tmp_112 * tmp_137;
   double a_0_0   = tmp_32 * tmp_34 + tmp_32 * tmp_36 + tmp_32 * tmp_38 + tmp_32 * tmp_40 + tmp_32 * tmp_42 - tmp_32 * tmp_43 -
                  tmp_32 * tmp_44 - tmp_32 * tmp_45 + tmp_32 * tmp_7;
   double a_0_1 = tmp_62;
   double a_0_2 = tmp_67;
   double a_0_3 = tmp_70;
   double a_0_4 = tmp_76;
   double a_0_5 = tmp_81;
   double a_1_0 = tmp_62;
   double a_1_1 = tmp_36 * tmp_82 + tmp_40 * tmp_82 - tmp_51 * tmp_82 - tmp_56 * tmp_82 - tmp_60 * tmp_82 + tmp_7 * tmp_82 +
                  tmp_82 * tmp_83 + tmp_82 * tmp_84 + tmp_82 * tmp_85;
   double a_1_2 = tmp_90;
   double a_1_3 = tmp_92;
   double a_1_4 = tmp_95;
   double a_1_5 = tmp_98;
   double a_2_0 = tmp_67;
   double a_2_1 = tmp_90;
   double a_2_2 = tmp_34 * tmp_99 + tmp_38 * tmp_99 + tmp_42 * tmp_99 - tmp_52 * tmp_99 - tmp_57 * tmp_99 - tmp_61 * tmp_99 +
                  tmp_83 * tmp_99 + tmp_84 * tmp_99 + tmp_85 * tmp_99;
   double a_2_3 = tmp_102;
   double a_2_4 = tmp_104;
   double a_2_5 = tmp_107;
   double a_3_0 = tmp_70;
   double a_3_1 = tmp_92;
   double a_3_2 = tmp_102;
   double a_3_3 = tmp_110 * tmp_36 + tmp_110 * tmp_40 + tmp_110 * tmp_43 + tmp_110 * tmp_44 + tmp_110 * tmp_45 +
                  tmp_110 * tmp_51 + tmp_110 * tmp_56 + tmp_110 * tmp_60 + tmp_110 * tmp_7 + tmp_111 * tmp_34 + tmp_111 * tmp_38 +
                  tmp_111 * tmp_42 + tmp_111 * tmp_83 + tmp_111 * tmp_84 + tmp_111 * tmp_85 + tmp_113 * tmp_52 +
                  tmp_113 * tmp_57 + tmp_113 * tmp_61;
   double a_3_4 = tmp_120;
   double a_3_5 = tmp_126;
   double a_4_0 = tmp_76;
   double a_4_1 = tmp_95;
   double a_4_2 = tmp_104;
   double a_4_3 = tmp_120;
   double a_4_4 = tmp_128 * tmp_36 + tmp_128 * tmp_40 + tmp_128 * tmp_7 + tmp_128 * tmp_83 + tmp_128 * tmp_84 + tmp_128 * tmp_85 +
                  tmp_129 * tmp_34 + tmp_129 * tmp_38 + tmp_129 * tmp_42 + tmp_129 * tmp_43 + tmp_129 * tmp_44 +
                  tmp_129 * tmp_45 + tmp_129 * tmp_52 + tmp_129 * tmp_57 + tmp_129 * tmp_61 + tmp_130 * tmp_51 +
                  tmp_130 * tmp_56 + tmp_130 * tmp_60;
   double a_4_5 = tmp_136;
   double a_5_0 = tmp_81;
   double a_5_1 = tmp_98;
   double a_5_2 = tmp_107;
   double a_5_3 = tmp_126;
   double a_5_4 = tmp_136;
   double a_5_5 = tmp_138 * tmp_34 + tmp_138 * tmp_36 + tmp_138 * tmp_38 + tmp_138 * tmp_40 + tmp_138 * tmp_42 + tmp_138 * tmp_7 +
                  tmp_139 * tmp_51 + tmp_139 * tmp_52 + tmp_139 * tmp_56 + tmp_139 * tmp_57 + tmp_139 * tmp_60 +
                  tmp_139 * tmp_61 + tmp_139 * tmp_83 + tmp_139 * tmp_84 + tmp_139 * tmp_85 + tmp_140 * tmp_43 +
                  tmp_140 * tmp_44 + tmp_140 * tmp_45;
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
