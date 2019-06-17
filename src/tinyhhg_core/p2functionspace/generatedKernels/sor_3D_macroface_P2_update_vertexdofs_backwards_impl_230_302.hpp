
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

void sor_3D_macroface_P2_update_vertexdofs_backwards_impl_230_302(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double const * RESTRICT const _data_edgeFaceDst_gl1_X, double const * RESTRICT const _data_edgeFaceDst_gl1_XY, double const * RESTRICT const _data_edgeFaceDst_gl1_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl1_XZ, double const * RESTRICT const _data_edgeFaceDst_gl1_Y, double const * RESTRICT const _data_edgeFaceDst_gl1_YZ, double const * RESTRICT const _data_edgeFaceDst_gl1_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double * RESTRICT _data_vertexFaceDst_gl1, double const * RESTRICT const _data_vertexFaceRhs, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil_fused_face_0, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil_fused_face_1, int32_t level, double relax, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_1);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg