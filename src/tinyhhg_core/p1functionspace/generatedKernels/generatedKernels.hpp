#pragma once

#include "core/DataTypes.h"

void faceApplyReplace(double * fd_p1FaceDst, const double * fd_p1FaceSrc, const double * fd_p1FaceStencil, const walberla::uint_t level);
void faceApplyAdd(double * fd_p1FaceDst, const double * fd_p1FaceSrc, const double * fd_p1FaceStencil, const walberla::uint_t level);