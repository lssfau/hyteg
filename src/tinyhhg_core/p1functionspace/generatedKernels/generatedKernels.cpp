#include "generatedKernels.hpp"

void faceApplyReplace( double*                fd_p1FaceDst,
                       const double*          fd_p1FaceSrc,
                       const double*          fd_p1FaceStencil,
                       const walberla::uint_t level )
{
   for( int ctr_2 = 1; ctr_2 < ( 1 << level ) - 1; ctr_2 += 1 )
      for( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << level ); ctr_1 += 1 )
      {
         fd_p1FaceDst[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] =
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 ) - 1] *
                 fd_p1FaceStencil[5] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 )] *
                 fd_p1FaceStencil[6] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 ) + 1] *
                 fd_p1FaceStencil[1] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 )] * fd_p1FaceStencil[0] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) + 1] * fd_p1FaceStencil[4] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) - 1] * fd_p1FaceStencil[2] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] * fd_p1FaceStencil[3];
      }
}

void faceApplyAdd( double*                fd_p1FaceDst,
                   const double*          fd_p1FaceSrc,
                   const double*          fd_p1FaceStencil,
                   const walberla::uint_t level )
{
   for( int ctr_2 = 1; ctr_2 < ( 1 << level ) - 1; ctr_2 += 1 )
      for( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << level ); ctr_1 += 1 )
      {
         fd_p1FaceDst[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] +=
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 ) - 1] *
                 fd_p1FaceStencil[5] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 )] *
                 fd_p1FaceStencil[6] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 ) + 1] *
                 fd_p1FaceStencil[1] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 )] * fd_p1FaceStencil[0] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) + 1] * fd_p1FaceStencil[4] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) - 1] * fd_p1FaceStencil[2] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] * fd_p1FaceStencil[3];
      }
}