
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

void sor_3D_macroface_P2_update_edgedofs_impl_013_320(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceDst_gl0_X, double * RESTRICT _data_edgeFaceDst_gl0_XYZ, double * RESTRICT _data_edgeFaceDst_gl0_XZ, double * RESTRICT _data_edgeFaceDst_gl0_Y, double * RESTRICT _data_edgeFaceDst_gl0_YZ, double * RESTRICT _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_edgeFaceDst_gl1_X, double * RESTRICT _data_edgeFaceDst_gl1_XYZ, double * RESTRICT _data_edgeFaceDst_gl1_XZ, double * RESTRICT _data_edgeFaceDst_gl1_Y, double * RESTRICT _data_edgeFaceDst_gl1_YZ, double * RESTRICT _data_edgeFaceDst_gl1_Z, double const * RESTRICT const _data_edgeFaceRhs_X, double const * RESTRICT const _data_edgeFaceRhs_XY, double const * RESTRICT const _data_edgeFaceRhs_Y, double const * RESTRICT const _data_vertexFaceDst, double const * RESTRICT const _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceDst_gl1, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2e_cell_stencil_fused_face_0, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2e_cell_stencil_fused_face_1, int32_t level, double relax, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2e_cell_stencil_fused_face_0, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2e_cell_stencil_fused_face_1);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg