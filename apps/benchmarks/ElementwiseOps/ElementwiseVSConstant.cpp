/*
 * Copyright (c) 2017-2020 Dominik Thoennes.
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
#include "core/Environment.h"
#include "core/Format.hpp"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/generatedKernels/apply_3D_macrocell_edgedof_to_edgedof_replace.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_3D_macrocell_edgedof_to_vertexdof_add.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_3D_macrocell_vertexdof_to_edgedof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/types/pointnd.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

static void evalQuadraturePoint3D(const Point3D& x_hat, const Point3D& x_tilde, const std::array<Point3D,4>& coords, const Matrix3r& DF, real_t w, Matrix10r& elMat)
{
   real_t tmp0 = DF(0,0)*DF(1,1);
   real_t tmp1 = DF(0,1)*DF(1,2);
   real_t tmp2 = DF(0,2)*DF(1,0);
   real_t tmp3 = DF(0,0)*DF(1,2);
   real_t tmp4 = DF(0,1)*DF(1,0);
   real_t tmp5 = DF(0,2)*DF(1,1);
   real_t tmp6 = DF(2,0)*tmp1 - DF(2,0)*tmp5 + DF(2,1)*tmp2 - DF(2,1)*tmp3 + DF(2,2)*tmp0 - DF(2,2)*tmp4;
   real_t tmp7 = coords[0][0]*coords[1][1];
   real_t tmp8 = coords[0][0]*coords[3][1];
   real_t tmp9 = coords[0][0]*coords[3][2];
   real_t tmp10 = coords[0][1]*coords[1][0];
   real_t tmp11 = coords[0][1]*coords[1][2];
   real_t tmp12 = coords[0][1]*coords[3][0];
   real_t tmp13 = coords[0][2]*coords[1][0];
   real_t tmp14 = coords[1][1]*coords[3][0];
   real_t tmp15 = coords[0][2]*coords[3][1];
   real_t tmp16 = coords[1][0]*coords[3][1];
   real_t tmp17 = coords[1][1]*coords[3][2];
   real_t tmp18 = coords[1][2]*coords[3][0];
   real_t tmp19 = coords[0][0]*coords[1][2];
   real_t tmp20 = coords[0][1]*coords[3][2];
   real_t tmp21 = coords[0][2]*coords[1][1];
   real_t tmp22 = coords[0][2]*coords[3][0];
   real_t tmp23 = coords[1][0]*coords[3][2];
   real_t tmp24 = coords[1][2]*coords[3][1];
   real_t tmp25 = coords[0][2]*tmp14 - coords[0][2]*tmp16 - coords[1][2]*tmp12 + coords[1][2]*tmp8 + coords[2][0]*tmp11 + coords[2][0]*tmp15 + coords[2][0]*tmp17 - coords[2][0]*tmp20 - coords[2][0]*tmp21 - coords[2][0]*tmp24 + coords[2][1]*tmp13 + coords[2][1]*tmp18 - coords[2][1]*tmp19 - coords[2][1]*tmp22 - coords[2][1]*tmp23 + coords[2][1]*tmp9 - coords[2][2]*tmp10 + coords[2][2]*tmp12 - coords[2][2]*tmp14 + coords[2][2]*tmp16 + coords[2][2]*tmp7 - coords[2][2]*tmp8 + coords[3][2]*tmp10 - coords[3][2]*tmp7;
   real_t tmp26 = std::fabs(tmp25*tmp6);
   real_t tmp27 = DF(0,0)*DF(2,1) - DF(0,1)*DF(2,0);
   real_t tmp28 = -tmp10 + tmp7;
   real_t tmp29 = tmp12 - tmp8;
   real_t tmp30 = -tmp14 + tmp16 + tmp28 + tmp29;
   real_t tmp31 = tmp27*tmp30;
   real_t tmp32 = DF(0,0)*DF(2,2) - DF(0,2)*DF(2,0);
   real_t tmp33 = -tmp13 + tmp19;
   real_t tmp34 = tmp22 - tmp9;
   real_t tmp35 = -tmp18 + tmp23 + tmp33 + tmp34;
   real_t tmp36 = tmp32*tmp35;
   real_t tmp37 = DF(0,1)*DF(2,2) - DF(0,2)*DF(2,1);
   real_t tmp38 = tmp11 - tmp21;
   real_t tmp39 = tmp15 - tmp20;
   real_t tmp40 = tmp17 - tmp24 + tmp38 + tmp39;
   real_t tmp41 = tmp37*tmp40;
   real_t tmp42 = coords[0][0]*coords[2][1];
   real_t tmp43 = coords[0][1]*coords[2][0];
   real_t tmp44 = coords[2][0]*coords[3][1] - coords[2][1]*coords[3][0] + tmp29 + tmp42 - tmp43;
   real_t tmp45 = coords[0][0]*coords[2][2];
   real_t tmp46 = coords[0][2]*coords[2][0];
   real_t tmp47 = coords[2][0]*coords[3][2] - coords[2][2]*coords[3][0] + tmp34 + tmp45 - tmp46;
   real_t tmp48 = coords[0][1]*coords[2][2];
   real_t tmp49 = coords[0][2]*coords[2][1];
   real_t tmp50 = coords[2][1]*coords[3][2] - coords[2][2]*coords[3][1] + tmp39 + tmp48 - tmp49;
   real_t tmp51 = tmp27*tmp44 + tmp32*tmp47 + tmp37*tmp50;
   real_t tmp52 = coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0] + tmp28 - tmp42 + tmp43;
   real_t tmp53 = coords[1][0]*coords[2][2] - coords[1][2]*coords[2][0] + tmp33 - tmp45 + tmp46;
   real_t tmp54 = coords[1][1]*coords[2][2] - coords[1][2]*coords[2][1] + tmp38 - tmp48 + tmp49;
   real_t tmp55 = tmp27*tmp52 + tmp32*tmp53 + tmp37*tmp54;
   real_t tmp56 = -tmp31 - tmp36 - tmp41 + tmp51 + tmp55;
   real_t tmp57 = 1.0 / (std::pow(tmp25, 2)*std::pow(tmp6, 2));
   real_t tmp58 = 16.0*tmp57;
   real_t tmp59 = tmp58*std::pow(x_hat[0] + x_hat[1] + x_hat[2] - 0.75, 2);
   real_t tmp60 = tmp0 - tmp4;
   real_t tmp61 = tmp52*tmp60;
   real_t tmp62 = tmp44*tmp60;
   real_t tmp63 = -tmp2 + tmp3;
   real_t tmp64 = tmp53*tmp63;
   real_t tmp65 = tmp47*tmp63;
   real_t tmp66 = tmp1 - tmp5;
   real_t tmp67 = tmp54*tmp66;
   real_t tmp68 = tmp50*tmp66;
   real_t tmp69 = tmp30*tmp60 + tmp35*tmp63 + tmp40*tmp66;
   real_t tmp70 = -tmp61 - tmp62 - tmp64 - tmp65 - tmp67 - tmp68 + tmp69;
   real_t tmp71 = DF(1,0)*DF(2,1) - DF(1,1)*DF(2,0);
   real_t tmp72 = tmp52*tmp71;
   real_t tmp73 = tmp44*tmp71;
   real_t tmp74 = DF(1,0)*DF(2,2) - DF(1,2)*DF(2,0);
   real_t tmp75 = tmp53*tmp74;
   real_t tmp76 = tmp47*tmp74;
   real_t tmp77 = DF(1,1)*DF(2,2) - DF(1,2)*DF(2,1);
   real_t tmp78 = tmp54*tmp77;
   real_t tmp79 = tmp50*tmp77;
   real_t tmp80 = tmp30*tmp71 + tmp35*tmp74 + tmp40*tmp77;
   real_t tmp81 = -tmp72 - tmp73 - tmp75 - tmp76 - tmp78 - tmp79 + tmp80;
   real_t tmp82 = 4.0*x_hat[0];
   real_t tmp83 = tmp82 - 1.0;
   real_t tmp84 = tmp51*tmp83;
   real_t tmp85 = 4.0*x_hat[2];
   real_t tmp86 = 4.0*x_hat[1];
   real_t tmp87 = tmp82 + tmp86;
   real_t tmp88 = tmp85 + tmp87 - 3.0;
   real_t tmp89 = tmp56*tmp88;
   real_t tmp90 = tmp57*tmp89;
   real_t tmp91 = tmp62 + tmp65 + tmp68;
   real_t tmp92 = tmp83*tmp91;
   real_t tmp93 = tmp70*tmp88;
   real_t tmp94 = tmp57*tmp93;
   real_t tmp95 = tmp73 + tmp76 + tmp79;
   real_t tmp96 = tmp83*tmp95;
   real_t tmp97 = tmp81*tmp88;
   real_t tmp98 = tmp57*tmp97;
   real_t tmp99 = tmp26*(tmp84*tmp90 - tmp92*tmp94 - tmp96*tmp98);
   real_t tmp100 = tmp31 + tmp36 + tmp41;
   real_t tmp101 = tmp86 - 1.0;
   real_t tmp102 = tmp100*tmp101;
   real_t tmp103 = tmp101*tmp69;
   real_t tmp104 = tmp101*tmp80;
   real_t tmp105 = tmp26*(-tmp102*tmp90 + tmp103*tmp94 + tmp104*tmp98);
   real_t tmp106 = tmp85 - 1.0;
   real_t tmp107 = tmp106*tmp55;
   real_t tmp108 = tmp61 + tmp64 + tmp67;
   real_t tmp109 = tmp106*tmp108;
   real_t tmp110 = tmp72 + tmp75 + tmp78;
   real_t tmp111 = tmp106*tmp110;
   real_t tmp112 = tmp26*(tmp107*tmp90 - tmp109*tmp94 - tmp111*tmp98);
   real_t tmp113 = tmp55*x_hat[1];
   real_t tmp114 = tmp100*x_hat[2];
   real_t tmp115 = tmp113 - tmp114;
   real_t tmp116 = 4.0*tmp115;
   real_t tmp117 = tmp69*x_hat[2];
   real_t tmp118 = tmp108*x_hat[1];
   real_t tmp119 = tmp117 - tmp118;
   real_t tmp120 = 4.0*tmp119;
   real_t tmp121 = tmp80*x_hat[2];
   real_t tmp122 = tmp110*x_hat[1];
   real_t tmp123 = tmp121 - tmp122;
   real_t tmp124 = 4.0*tmp123;
   real_t tmp125 = tmp26*(tmp116*tmp90 + tmp120*tmp94 + tmp124*tmp98);
   real_t tmp126 = tmp55*x_hat[0];
   real_t tmp127 = tmp51*x_hat[2];
   real_t tmp128 = tmp126 + tmp127;
   real_t tmp129 = 4.0*tmp128;
   real_t tmp130 = tmp108*tmp82;
   real_t tmp131 = tmp85*tmp91;
   real_t tmp132 = tmp57*(tmp130 + tmp131);
   real_t tmp133 = tmp110*tmp82;
   real_t tmp134 = tmp85*tmp95;
   real_t tmp135 = tmp57*(tmp133 + tmp134);
   real_t tmp136 = tmp26*(tmp129*tmp90 - tmp132*tmp93 - tmp135*tmp97);
   real_t tmp137 = tmp51*x_hat[1];
   real_t tmp138 = tmp100*x_hat[0];
   real_t tmp139 = tmp137 - tmp138;
   real_t tmp140 = 4.0*tmp139;
   real_t tmp141 = tmp69*x_hat[0];
   real_t tmp142 = tmp91*x_hat[1];
   real_t tmp143 = tmp141 - tmp142;
   real_t tmp144 = 4.0*tmp143;
   real_t tmp145 = tmp80*x_hat[0];
   real_t tmp146 = tmp95*x_hat[1];
   real_t tmp147 = tmp145 - tmp146;
   real_t tmp148 = 4.0*tmp147;
   real_t tmp149 = tmp26*(tmp140*tmp90 + tmp144*tmp94 + tmp148*tmp98);
   real_t tmp150 = tmp87 + 8.0*x_hat[2] - 4.0;
   real_t tmp151 = tmp150*tmp55;
   real_t tmp152 = tmp57*(tmp100*tmp85 - tmp151 - tmp51*tmp85);
   real_t tmp153 = tmp108*tmp150;
   real_t tmp154 = tmp131 + tmp153 - tmp69*tmp85;
   real_t tmp155 = tmp154*tmp57;
   real_t tmp156 = tmp110*tmp150;
   real_t tmp157 = tmp134 + tmp156 - tmp80*tmp85;
   real_t tmp158 = tmp157*tmp57;
   real_t tmp159 = tmp26*(tmp152*tmp89 + tmp155*tmp93 + tmp158*tmp97);
   real_t tmp160 = tmp85 - 4.0;
   real_t tmp161 = tmp160 + tmp82 + 8.0*x_hat[1];
   real_t tmp162 = tmp100*tmp161;
   real_t tmp163 = tmp162 - tmp51*tmp86 - tmp55*tmp86;
   real_t tmp164 = tmp163*tmp57;
   real_t tmp165 = tmp161*tmp69;
   real_t tmp166 = tmp108*tmp86 - tmp165 + tmp86*tmp91;
   real_t tmp167 = tmp166*tmp57;
   real_t tmp168 = tmp161*tmp80;
   real_t tmp169 = tmp110*tmp86 - tmp168 + tmp86*tmp95;
   real_t tmp170 = tmp169*tmp57;
   real_t tmp171 = tmp26*(tmp164*tmp89 + tmp167*tmp93 + tmp170*tmp97);
   real_t tmp172 = tmp160 + tmp86 + 8.0*x_hat[0];
   real_t tmp173 = tmp172*tmp51;
   real_t tmp174 = tmp100*tmp82 - tmp173 - tmp55*tmp82;
   real_t tmp175 = tmp174*tmp57;
   real_t tmp176 = tmp172*tmp91;
   real_t tmp177 = tmp130 + tmp176 - tmp69*tmp82;
   real_t tmp178 = tmp177*tmp57;
   real_t tmp179 = tmp172*tmp95;
   real_t tmp180 = tmp133 + tmp179 - tmp80*tmp82;
   real_t tmp181 = tmp180*tmp57;
   real_t tmp182 = tmp26*(tmp175*tmp89 + tmp178*tmp93 + tmp181*tmp97);
   real_t tmp183 = tmp58*std::pow(x_hat[0] - 0.25, 2);
   real_t tmp184 = tmp57*tmp92;
   real_t tmp185 = tmp57*tmp84;
   real_t tmp186 = tmp57*tmp96;
   real_t tmp187 = tmp26*(-tmp102*tmp185 - tmp103*tmp184 - tmp104*tmp186);
   real_t tmp188 = tmp26*(tmp107*tmp185 + tmp109*tmp184 + tmp111*tmp186);
   real_t tmp189 = tmp26*(tmp116*tmp185 - tmp120*tmp184 - tmp124*tmp186);
   real_t tmp190 = tmp26*(tmp129*tmp185 + tmp132*tmp92 + tmp135*tmp96);
   real_t tmp191 = tmp26*(tmp140*tmp185 - tmp144*tmp184 - tmp148*tmp186);
   real_t tmp192 = tmp26*(tmp152*tmp84 - tmp155*tmp92 - tmp158*tmp96);
   real_t tmp193 = tmp26*(tmp164*tmp84 - tmp167*tmp92 - tmp170*tmp96);
   real_t tmp194 = tmp26*(tmp175*tmp84 - tmp178*tmp92 - tmp181*tmp96);
   real_t tmp195 = tmp58*std::pow(x_hat[1] - 0.25, 2);
   real_t tmp196 = tmp103*tmp57;
   real_t tmp197 = tmp102*tmp57;
   real_t tmp198 = tmp104*tmp57;
   real_t tmp199 = tmp26*(-tmp107*tmp197 - tmp109*tmp196 - tmp111*tmp198);
   real_t tmp200 = tmp26*(-tmp116*tmp197 + tmp120*tmp196 + tmp124*tmp198);
   real_t tmp201 = tmp26*(-tmp103*tmp132 - tmp104*tmp135 - tmp129*tmp197);
   real_t tmp202 = tmp26*(-tmp140*tmp197 + tmp144*tmp196 + tmp148*tmp198);
   real_t tmp203 = tmp26*(-tmp102*tmp152 + tmp103*tmp155 + tmp104*tmp158);
   real_t tmp204 = tmp26*(-tmp102*tmp164 + tmp103*tmp167 + tmp104*tmp170);
   real_t tmp205 = tmp26*(-tmp102*tmp175 + tmp103*tmp178 + tmp104*tmp181);
   real_t tmp206 = tmp58*std::pow(x_hat[2] - 0.25, 2);
   real_t tmp207 = tmp109*tmp57;
   real_t tmp208 = tmp107*tmp57;
   real_t tmp209 = tmp111*tmp57;
   real_t tmp210 = tmp26*(tmp116*tmp208 - tmp120*tmp207 - tmp124*tmp209);
   real_t tmp211 = tmp26*(tmp109*tmp132 + tmp111*tmp135 + tmp129*tmp208);
   real_t tmp212 = tmp26*(tmp140*tmp208 - tmp144*tmp207 - tmp148*tmp209);
   real_t tmp213 = tmp26*(tmp107*tmp152 - tmp109*tmp155 - tmp111*tmp158);
   real_t tmp214 = tmp26*(tmp107*tmp164 - tmp109*tmp167 - tmp111*tmp170);
   real_t tmp215 = tmp26*(tmp107*tmp175 - tmp109*tmp178 - tmp111*tmp181);
   real_t tmp216 = tmp115*tmp58;
   real_t tmp217 = tmp26*(-tmp120*tmp132 - tmp124*tmp135 + tmp128*tmp216);
   real_t tmp218 = tmp26*(tmp119*tmp143*tmp58 + tmp123*tmp147*tmp58 + tmp139*tmp216);
   real_t tmp219 = tmp26*(tmp116*tmp152 + tmp120*tmp155 + tmp124*tmp158);
   real_t tmp220 = tmp26*(tmp116*tmp164 + tmp120*tmp167 + tmp124*tmp170);
   real_t tmp221 = tmp26*(tmp116*tmp175 + tmp120*tmp178 + tmp124*tmp181);
   real_t tmp222 = tmp108*x_hat[0];
   real_t tmp223 = tmp91*x_hat[2];
   real_t tmp224 = tmp110*x_hat[0];
   real_t tmp225 = tmp95*x_hat[2];
   real_t tmp226 = tmp26*(tmp128*tmp139*tmp58 - tmp132*tmp144 - tmp135*tmp148);
   real_t tmp227 = tmp26*(tmp129*tmp152 - tmp132*tmp154 - tmp135*tmp157);
   real_t tmp228 = tmp26*(tmp129*tmp164 - tmp132*tmp166 - tmp135*tmp169);
   real_t tmp229 = tmp26*(tmp129*tmp175 - tmp132*tmp177 - tmp135*tmp180);
   real_t tmp230 = tmp26*(tmp140*tmp152 + tmp144*tmp155 + tmp148*tmp158);
   real_t tmp231 = tmp26*(tmp140*tmp164 + tmp144*tmp167 + tmp148*tmp170);
   real_t tmp232 = tmp26*(tmp140*tmp175 + tmp144*tmp178 + tmp148*tmp181);
   real_t tmp233 = tmp26*(tmp152*tmp163 + tmp155*tmp166 + tmp158*tmp169);
   real_t tmp234 = tmp26*(tmp152*tmp174 + tmp155*tmp177 + tmp158*tmp180);
   real_t tmp235 = tmp26*(tmp164*tmp174 + tmp167*tmp177 + tmp170*tmp180);
   elMat(0,0) += w * tmp26*(std::pow(tmp56, 2)*tmp59 + tmp59*std::pow(tmp70, 2) + tmp59*std::pow(tmp81, 2));
   elMat(0,1) += w * tmp99;
   elMat(0,2) += w * tmp105;
   elMat(0,3) += w * tmp112;
   elMat(0,4) += w * tmp125;
   elMat(0,5) += w * tmp136;
   elMat(0,6) += w * tmp149;
   elMat(0,7) += w * tmp159;
   elMat(0,8) += w * tmp171;
   elMat(0,9) += w * tmp182;
   elMat(1,0) += w * tmp99;
   elMat(1,1) += w * tmp26*(tmp183*std::pow(tmp51, 2) + tmp183*std::pow(tmp91, 2) + tmp183*std::pow(tmp95, 2));
   elMat(1,2) += w * tmp187;
   elMat(1,3) += w * tmp188;
   elMat(1,4) += w * tmp189;
   elMat(1,5) += w * tmp190;
   elMat(1,6) += w * tmp191;
   elMat(1,7) += w * tmp192;
   elMat(1,8) += w * tmp193;
   elMat(1,9) += w * tmp194;
   elMat(2,0) += w * tmp105;
   elMat(2,1) += w * tmp187;
   elMat(2,2) += w * tmp26*(std::pow(tmp100, 2)*tmp195 + tmp195*std::pow(tmp69, 2) + tmp195*std::pow(tmp80, 2));
   elMat(2,3) += w * tmp199;
   elMat(2,4) += w * tmp200;
   elMat(2,5) += w * tmp201;
   elMat(2,6) += w * tmp202;
   elMat(2,7) += w * tmp203;
   elMat(2,8) += w * tmp204;
   elMat(2,9) += w * tmp205;
   elMat(3,0) += w * tmp112;
   elMat(3,1) += w * tmp188;
   elMat(3,2) += w * tmp199;
   elMat(3,3) += w * tmp26*(std::pow(tmp108, 2)*tmp206 + std::pow(tmp110, 2)*tmp206 + tmp206*std::pow(tmp55, 2));
   elMat(3,4) += w * tmp210;
   elMat(3,5) += w * tmp211;
   elMat(3,6) += w * tmp212;
   elMat(3,7) += w * tmp213;
   elMat(3,8) += w * tmp214;
   elMat(3,9) += w * tmp215;
   elMat(4,0) += w * tmp125;
   elMat(4,1) += w * tmp189;
   elMat(4,2) += w * tmp200;
   elMat(4,3) += w * tmp210;
   elMat(4,4) += w * tmp26*(std::pow(tmp115, 2)*tmp58 + std::pow(tmp119, 2)*tmp58 + std::pow(tmp123, 2)*tmp58);
   elMat(4,5) += w * tmp217;
   elMat(4,6) += w * tmp218;
   elMat(4,7) += w * tmp219;
   elMat(4,8) += w * tmp220;
   elMat(4,9) += w * tmp221;
   elMat(5,0) += w * tmp136;
   elMat(5,1) += w * tmp190;
   elMat(5,2) += w * tmp201;
   elMat(5,3) += w * tmp211;
   elMat(5,4) += w * tmp217;
   elMat(5,5) += w * tmp26*(std::pow(tmp128, 2)*tmp58 + tmp58*std::pow(tmp222 + tmp223, 2) + tmp58*std::pow(tmp224 + tmp225, 2));
   elMat(5,6) += w * tmp226;
   elMat(5,7) += w * tmp227;
   elMat(5,8) += w * tmp228;
   elMat(5,9) += w * tmp229;
   elMat(6,0) += w * tmp149;
   elMat(6,1) += w * tmp191;
   elMat(6,2) += w * tmp202;
   elMat(6,3) += w * tmp212;
   elMat(6,4) += w * tmp218;
   elMat(6,5) += w * tmp226;
   elMat(6,6) += w * tmp26*(std::pow(tmp139, 2)*tmp58 + std::pow(tmp143, 2)*tmp58 + std::pow(tmp147, 2)*tmp58);
   elMat(6,7) += w * tmp230;
   elMat(6,8) += w * tmp231;
   elMat(6,9) += w * tmp232;
   elMat(7,0) += w * tmp159;
   elMat(7,1) += w * tmp192;
   elMat(7,2) += w * tmp203;
   elMat(7,3) += w * tmp213;
   elMat(7,4) += w * tmp219;
   elMat(7,5) += w * tmp227;
   elMat(7,6) += w * tmp230;
   elMat(7,7) += w * tmp26*(tmp58*std::pow(tmp114 - tmp127 - 0.25*tmp151, 2) + tmp58*std::pow(-tmp117 + 0.25*tmp153 + tmp223, 2) + tmp58*std::pow(-tmp121 + 0.25*tmp156 + tmp225, 2));
   elMat(7,8) += w * tmp233;
   elMat(7,9) += w * tmp234;
   elMat(8,0) += w * tmp171;
   elMat(8,1) += w * tmp193;
   elMat(8,2) += w * tmp204;
   elMat(8,3) += w * tmp214;
   elMat(8,4) += w * tmp220;
   elMat(8,5) += w * tmp228;
   elMat(8,6) += w * tmp231;
   elMat(8,7) += w * tmp233;
   elMat(8,8) += w * tmp26*(tmp58*std::pow(-tmp113 - tmp137 + 0.25*tmp162, 2) + tmp58*std::pow(tmp118 + tmp142 - 0.25*tmp165, 2) + tmp58*std::pow(tmp122 + tmp146 - 0.25*tmp168, 2));
   elMat(8,9) += w * tmp235;
   elMat(9,0) += w * tmp182;
   elMat(9,1) += w * tmp194;
   elMat(9,2) += w * tmp205;
   elMat(9,3) += w * tmp215;
   elMat(9,4) += w * tmp221;
   elMat(9,5) += w * tmp229;
   elMat(9,6) += w * tmp232;
   elMat(9,7) += w * tmp234;
   elMat(9,8) += w * tmp235;
   elMat(9,9) += w * tmp26*(tmp58*std::pow(-tmp126 + tmp138 - 0.25*tmp173, 2) + tmp58*std::pow(-tmp141 + 0.25*tmp176 + tmp222, 2) + tmp58*std::pow(-tmp145 + 0.25*tmp179 + tmp224, 2));
}
template < class P2ElementwiseOperator >
void performBenchmark( hyteg::P2Function< double >&      src,
                       hyteg::P2Function< double >&      dstConst,
                       hyteg::P2Function< double >&      dstElem,
                       hyteg::P2ConstantLaplaceOperator& constantOperator,
                       P2ElementwiseOperator&            elementwiseOperator,
                       hyteg::Cell&                      cell,
                       walberla::WcTimingTree&           timingTree,
                       uint_t                            level,
                       uint_t                            startiterations )
{
   const std::string benchInfoString = "level" + ( level < 10 ? "0" + std::to_string( level ) : std::to_string( level ) ) +
                                       "-numProcs" + std::to_string( walberla::mpi::MPIManager::instance()->numProcesses() );

   const std::string cOString   = "constOperator-" + benchInfoString;
   const std::string eOString   = "elementOperator-" + benchInfoString;
   uint_t            iterations = startiterations;

   auto dstConstVertexPtr  = cell.getData( dstConst.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto dstConstEdgePtr    = cell.getData( dstConst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto dstElemVertexPtr   = cell.getData( dstElem.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto dstElemEdgePtr     = cell.getData( dstElem.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto srcVertexPtr       = cell.getData( src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto srcEdgePtr         = cell.getData( src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto const_v2v_opr_data = cell.getData( constantOperator.getVertexToVertexOpr().getCellStencilID() )->getData( level );
   auto const_v2e_opr_data = cell.getData( constantOperator.getVertexToEdgeOpr().getCellStencilID() )->getData( level );
   auto const_e2v_opr_data = cell.getData( constantOperator.getEdgeToVertexOpr().getCellStencilID() )->getData( level );
   auto const_e2e_opr_data = cell.getData( constantOperator.getEdgeToEdgeOpr().getCellStencilID() )->getData( level );

   typedef hyteg::edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >                     firstEdgeIdx;
   for ( auto e : hyteg::edgedof::allEdgeDoFOrientations )
      firstEdgeIdx[e] = hyteg::edgedof::macrocell::index( level, 0, 0, 0, e );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s|%10s|%15s|%15s", "Operator", "Time (s)", "iteration", "timer/iter" ) )

   do
   {
      timingTree.start( cOString );
      ///only works with likwid 4.3.3 and higher
      LIKWID_MARKER_RESET( cOString.c_str() );
      LIKWID_MARKER_START( cOString.c_str() );

      for ( uint_t iter = 0; iter < iterations; ++iter )
      {
         hyteg::vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
             dstConstVertexPtr, srcVertexPtr, static_cast< int32_t >( level ), const_v2v_opr_data );

         hyteg::edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_replace(
             &dstConstEdgePtr[firstEdgeIdx[eo::X]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XY]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XYZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Y]],
             &dstConstEdgePtr[firstEdgeIdx[eo::YZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Z]],
             &srcEdgePtr[firstEdgeIdx[eo::X]],
             &srcEdgePtr[firstEdgeIdx[eo::XY]],
             &srcEdgePtr[firstEdgeIdx[eo::XYZ]],
             &srcEdgePtr[firstEdgeIdx[eo::XZ]],
             &srcEdgePtr[firstEdgeIdx[eo::Y]],
             &srcEdgePtr[firstEdgeIdx[eo::YZ]],
             &srcEdgePtr[firstEdgeIdx[eo::Z]],
             const_e2e_opr_data,
             static_cast< int32_t >( level ) );

         hyteg::VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_add(
             &dstConstEdgePtr[firstEdgeIdx[eo::X]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XY]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XYZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Y]],
             &dstConstEdgePtr[firstEdgeIdx[eo::YZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Z]],
             srcVertexPtr,
             static_cast< int32_t >( level ),
             const_v2e_opr_data );

         hyteg::EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_add( &srcEdgePtr[firstEdgeIdx[eo::X]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::XY]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::XYZ]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::XZ]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::Y]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::YZ]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::Z]],
                                                                                            dstConstVertexPtr,
                                                                                            const_e2v_opr_data,
                                                                                            static_cast< int32_t >( level ) );
      }
      WALBERLA_MPI_BARRIER()
      LIKWID_MARKER_STOP( cOString.c_str() );
      timingTree.stop( cOString );
      iterations *= 2;
   } while ( timingTree[cOString].last() < 0.5 );
   iterations /= 2;
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s|%10.5f|%15u|%15.5f",
                                                "P2Constant",
                                                timingTree[cOString].last(),
                                                iterations,
                                                timingTree[cOString].last() / real_c( iterations ) ) )

   iterations = startiterations;

   do
   {
      timingTree.start( eOString );
      ///only works with likwid 4.3.3 and higher
      LIKWID_MARKER_RESET( eOString.c_str() );
      LIKWID_MARKER_START( eOString.c_str() );
      for ( uint_t iter = 0; iter < iterations; ++iter )
      {
         dstElem.interpolate( []( const hyteg::Point3D& ) { return 0; }, level, hyteg::All );
         for ( const auto& idx : hyteg::vertexdof::macrocell::Iterator( level ) )
         {
            if ( !hyteg::vertexdof::macrocell::isOnCellFace( idx, level ).empty() )
            {
               auto arrayIdx              = hyteg::vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() );
               dstElemVertexPtr[arrayIdx] = real_c( 0 );
            }
         }

         for ( const auto& idx : hyteg::edgedof::macrocell::Iterator( level ) )
         {
            for ( const auto& orientation : hyteg::edgedof::allEdgeDoFOrientationsWithoutXYZ )
            {
               if ( !hyteg::edgedof::macrocell::isInnerEdgeDoF( level, idx, orientation ) )
               {
                  auto arrayIdx            = hyteg::edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), orientation );
                  dstElemEdgePtr[arrayIdx] = real_c( 0 );
               }
            }
         }

         // loop over micro-cells
         int end = 0;
         auto form = elementwiseOperator.getForm();
         auto geometryMap_ = cell.getGeometryMap();
         for ( const auto& cType : hyteg::celldof::allCellTypes )
         {
            if ( cType == hyteg::celldof::CellType::WHITE_UP )
            {
               end = 0;
            }
            else if ( cType == hyteg::celldof::CellType::WHITE_DOWN )
            {
               end = 2;
            }
            else
            {
               end = 1;
            }
            for ( int ctr_3 = 0; ctr_3 < ( 1 << ( level ) ) - end; ctr_3 += 1 )
            {
               for ( int ctr_2 = 0; ctr_2 < -ctr_3 + ( 1 << ( level ) ) - end; ctr_2 += 1 )
               {
                  // cell (inner)
                  for ( int ctr_1 = 0; ctr_1 < -ctr_2 - ctr_3 + ( 1 << ( level ) ) - end; ctr_1 += 1 )
                  {
                     hyteg::indexing::Index microCell( {ctr_1, ctr_2, ctr_3} );
                     // determine coordinates of vertices of micro-element
                     std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
                     std::array< Point3D, 4 >         coords;
                     for ( uint_t k = 0; k < 4; ++k )
                     {
                        coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
                     }

                     // assemble local element matrix
                     Matrix10r elMat;

                     Point3D x_hat;
                     Point3D x_tilde;
                     Matrix3r DF;
                     elMat(0,0) = 0;
                     elMat(0,1) = 0;
                     elMat(0,2) = 0;
                     elMat(0,3) = 0;
                     elMat(0,4) = 0;
                     elMat(0,5) = 0;
                     elMat(0,6) = 0;
                     elMat(0,7) = 0;
                     elMat(0,8) = 0;
                     elMat(0,9) = 0;
                     elMat(1,0) = 0;
                     elMat(1,1) = 0;
                     elMat(1,2) = 0;
                     elMat(1,3) = 0;
                     elMat(1,4) = 0;
                     elMat(1,5) = 0;
                     elMat(1,6) = 0;
                     elMat(1,7) = 0;
                     elMat(1,8) = 0;
                     elMat(1,9) = 0;
                     elMat(2,0) = 0;
                     elMat(2,1) = 0;
                     elMat(2,2) = 0;
                     elMat(2,3) = 0;
                     elMat(2,4) = 0;
                     elMat(2,5) = 0;
                     elMat(2,6) = 0;
                     elMat(2,7) = 0;
                     elMat(2,8) = 0;
                     elMat(2,9) = 0;
                     elMat(3,0) = 0;
                     elMat(3,1) = 0;
                     elMat(3,2) = 0;
                     elMat(3,3) = 0;
                     elMat(3,4) = 0;
                     elMat(3,5) = 0;
                     elMat(3,6) = 0;
                     elMat(3,7) = 0;
                     elMat(3,8) = 0;
                     elMat(3,9) = 0;
                     elMat(4,0) = 0;
                     elMat(4,1) = 0;
                     elMat(4,2) = 0;
                     elMat(4,3) = 0;
                     elMat(4,4) = 0;
                     elMat(4,5) = 0;
                     elMat(4,6) = 0;
                     elMat(4,7) = 0;
                     elMat(4,8) = 0;
                     elMat(4,9) = 0;
                     elMat(5,0) = 0;
                     elMat(5,1) = 0;
                     elMat(5,2) = 0;
                     elMat(5,3) = 0;
                     elMat(5,4) = 0;
                     elMat(5,5) = 0;
                     elMat(5,6) = 0;
                     elMat(5,7) = 0;
                     elMat(5,8) = 0;
                     elMat(5,9) = 0;
                     elMat(6,0) = 0;
                     elMat(6,1) = 0;
                     elMat(6,2) = 0;
                     elMat(6,3) = 0;
                     elMat(6,4) = 0;
                     elMat(6,5) = 0;
                     elMat(6,6) = 0;
                     elMat(6,7) = 0;
                     elMat(6,8) = 0;
                     elMat(6,9) = 0;
                     elMat(7,0) = 0;
                     elMat(7,1) = 0;
                     elMat(7,2) = 0;
                     elMat(7,3) = 0;
                     elMat(7,4) = 0;
                     elMat(7,5) = 0;
                     elMat(7,6) = 0;
                     elMat(7,7) = 0;
                     elMat(7,8) = 0;
                     elMat(7,9) = 0;
                     elMat(8,0) = 0;
                     elMat(8,1) = 0;
                     elMat(8,2) = 0;
                     elMat(8,3) = 0;
                     elMat(8,4) = 0;
                     elMat(8,5) = 0;
                     elMat(8,6) = 0;
                     elMat(8,7) = 0;
                     elMat(8,8) = 0;
                     elMat(8,9) = 0;
                     elMat(9,0) = 0;
                     elMat(9,1) = 0;
                     elMat(9,2) = 0;
                     elMat(9,3) = 0;
                     elMat(9,4) = 0;
                     elMat(9,5) = 0;
                     elMat(9,6) = 0;
                     elMat(9,7) = 0;
                     elMat(9,8) = 0;
                     elMat(9,9) = 0;
                     x_hat[0] = 0.1381966011250105;
                     x_hat[1] = 0.1381966011250105;
                     x_hat[2] = 0.5854101966249684;
                     x_tilde[0] = 0.138196601125011*coords[0][0] + 0.138196601125011*coords[1][0] + 0.138196601125011*coords[2][0] + 0.585410196624968*coords[3][0];
                     x_tilde[1] = 0.138196601125011*coords[0][1] + 0.138196601125011*coords[1][1] + 0.138196601125011*coords[2][1] + 0.585410196624968*coords[3][1];
                     x_tilde[2] = 0.138196601125011*coords[0][2] + 0.138196601125011*coords[1][2] + 0.138196601125011*coords[2][2] + 0.585410196624968*coords[3][2];
                     geometryMap_->evalDF(x_tilde, DF);
                     evalQuadraturePoint3D(x_hat, x_tilde, coords, DF, 0.04166666666666666, elMat);
                     x_hat[0] = 0.1381966011250105;
                     x_hat[1] = 0.5854101966249684;
                     x_hat[2] = 0.1381966011250105;
                     x_tilde[0] = 0.138196601125011*coords[0][0] + 0.138196601125011*coords[1][0] + 0.585410196624968*coords[2][0] + 0.138196601125011*coords[3][0];
                     x_tilde[1] = 0.138196601125011*coords[0][1] + 0.138196601125011*coords[1][1] + 0.585410196624968*coords[2][1] + 0.138196601125011*coords[3][1];
                     x_tilde[2] = 0.138196601125011*coords[0][2] + 0.138196601125011*coords[1][2] + 0.585410196624968*coords[2][2] + 0.138196601125011*coords[3][2];
                     geometryMap_->evalDF(x_tilde, DF);
                     evalQuadraturePoint3D(x_hat, x_tilde, coords, DF, 0.04166666666666666, elMat);
                     x_hat[0] = 0.5854101966249684;
                     x_hat[1] = 0.1381966011250105;
                     x_hat[2] = 0.1381966011250105;
                     x_tilde[0] = 0.138196601125011*coords[0][0] + 0.585410196624968*coords[1][0] + 0.138196601125011*coords[2][0] + 0.138196601125011*coords[3][0];
                     x_tilde[1] = 0.138196601125011*coords[0][1] + 0.585410196624968*coords[1][1] + 0.138196601125011*coords[2][1] + 0.138196601125011*coords[3][1];
                     x_tilde[2] = 0.138196601125011*coords[0][2] + 0.585410196624968*coords[1][2] + 0.138196601125011*coords[2][2] + 0.138196601125011*coords[3][2];
                     geometryMap_->evalDF(x_tilde, DF);
                     evalQuadraturePoint3D(x_hat, x_tilde, coords, DF, 0.04166666666666666, elMat);
                     x_hat[0] = 0.1381966011250105;
                     x_hat[1] = 0.1381966011250105;
                     x_hat[2] = 0.1381966011250105;
                     x_tilde[0] = 0.585410196624968*coords[0][0] + 0.138196601125011*coords[1][0] + 0.138196601125011*coords[2][0] + 0.138196601125011*coords[3][0];
                     x_tilde[1] = 0.585410196624968*coords[0][1] + 0.138196601125011*coords[1][1] + 0.138196601125011*coords[2][1] + 0.138196601125011*coords[3][1];
                     x_tilde[2] = 0.585410196624968*coords[0][2] + 0.138196601125011*coords[1][2] + 0.138196601125011*coords[2][2] + 0.138196601125011*coords[3][2];
                     geometryMap_->evalDF(x_tilde, DF);
                     evalQuadraturePoint3D(x_hat, x_tilde, coords, DF, 0.04166666666666666, elMat);

                     // obtain data indices of dofs associated with micro-cell
                     std::array< uint_t, 4 > vertexDoFIndices;
                     vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

                     std::array< uint_t, 6 > edgeDoFIndices;
                     edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

                     // assemble local element vector
                     Point10D elVecOld, elVecNew;
                     for ( uint_t k = 0; k < 4; ++k )
                     {
                        elVecOld[k] = srcVertexPtr[vertexDoFIndices[k]];
                     }
                     for ( uint_t k = 4; k < 10; ++k )
                     {
                        elVecOld[k] = srcEdgePtr[edgeDoFIndices[k - 4]];
                     }

                     // apply matrix (operator locally)
                     elVecNew = elMat.mul( elVecOld );

                     // redistribute result from "local" to "global vector"
                     for ( uint_t k = 0; k < 4; ++k )
                     {
                        dstElemVertexPtr[vertexDoFIndices[k]] += elVecNew[k];
                     }
                     for ( uint_t k = 4; k < 10; ++k )
                     {
                        dstElemEdgePtr[edgeDoFIndices[k - 4]] += elVecNew[k];
                     }
                  }
               }
            }
         }
      }
      WALBERLA_MPI_BARRIER()
      LIKWID_MARKER_STOP( eOString.c_str() );
      timingTree.stop( eOString );
      iterations *= 2;
   } while ( timingTree[eOString].last() < 0.5 );
   iterations /= 2;
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s|%10.5f|%15u|%15.5f",
                                                "P2Elementwise",
                                                timingTree[eOString].last(),
                                                iterations,
                                                timingTree[eOString].last() / real_c( iterations ) ) )
}
int main( int argc, char** argv )
{
   uint_t minLevel   = 4;
   uint_t maxLevel   = 4;
   uint_t benchLevel = 4;

   LIKWID_MARKER_INIT;
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   walberla::WcTimingTree timingTree;

   auto meshInfo = hyteg::MeshInfo::meshCuboid( hyteg::Point3D( {0, 0, 0} ), hyteg::Point3D( {1, 1, 1} ), 1, 1, 1 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   //setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return std::sin( walberla::math::pi * xx[0] ) + std::cos( walberla::math::pi * xx[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > zeros = []( const hyteg::Point3D& ) { return 0; };

   hyteg::P2Function< double > src( "src", storage, minLevel, maxLevel );
   hyteg::P2Function< double > dstConst( "dstConst", storage, minLevel, maxLevel );
   hyteg::P2Function< double > dstElem( "dstElem", storage, minLevel, maxLevel );
   hyteg::P2Function< double > diff( "diff", storage, minLevel, maxLevel );

   hyteg::P2ConstantLaplaceOperator constantOperator( storage, minLevel, maxLevel );
   //hyteg::P2ElementwiseLaplaceOperator elementWiseOperator( storage, minLevel, maxLevel );
   //hyteg::P2ElementwiseBlendingLaplaceOperator elementWiseOperator( storage, minLevel, maxLevel );
   hyteg::P2ElementwiseDivKGradBlendingOperator elementWiseOperator( storage, minLevel, maxLevel );

   src.interpolate( exact, benchLevel, hyteg::All );
   dstConst.interpolate( zeros, benchLevel, hyteg::All );
   dstElem.interpolate( zeros, benchLevel, hyteg::All );
   diff.interpolate( zeros, benchLevel, hyteg::All );
   hyteg::communication::syncP2FunctionBetweenPrimitives( src, benchLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( dstConst, benchLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( dstElem, benchLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( diff, benchLevel );

   //each process will get its first tet here
   std::vector< hyteg::PrimitiveID > macroCells = storage->getCellIDs();
   WALBERLA_CHECK_GREATER_EQUAL( macroCells.size(), 1 )
   hyteg::Cell* cell = storage->getCell( macroCells.front() );

   performBenchmark( src, dstConst, dstElem, constantOperator, elementWiseOperator, *cell, timingTree, benchLevel, 1 );

   //compar results
   diff.assign( {1.0, -1.0}, {dstConst, dstElem}, benchLevel, hyteg::All );
   WALBERLA_LOG_INFO_ON_ROOT( "Diff Magnitude: " << diff.getMaxMagnitude( benchLevel ) )

   LIKWID_MARKER_CLOSE;
   // hyteg::VTKOutput vtkOutput( ".", "P2ConstantVSP2Elementwise", storage );
   // vtkOutput.add( src );
   // vtkOutput.add( diff );
   // vtkOutput.add( dstConst );
   // vtkOutput.add( dstElem );
   // vtkOutput.write( benchLevel, 0 );
}
