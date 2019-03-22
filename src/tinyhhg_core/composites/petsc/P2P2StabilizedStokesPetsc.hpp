#pragma once

#include "tinyhhg_core/composites/P2P2StokesFunction.hpp"
#include "tinyhhg_core/p2functionspace/P2Petsc.hpp"

namespace hhg {
namespace petsc {

inline void createVectorFromFunction( const P2P2StokesFunction< PetscScalar >& function,
                                      const P2P2StokesFunction< PetscInt >&    numerator,
                                      Vec&                                     vec,
                                      uint_t                                   level,
                                      DoFType                                  flag )
{
   createVectorFromFunction( function.u, numerator.u, vec, level, flag );
   createVectorFromFunction( function.v, numerator.v, vec, level, flag );
   if ( function.u.getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented for 3D" );
      // createVectorFromFunction( function.w, numerator.w, vec, level, flag );
   }
   createVectorFromFunction( function.p, numerator.p, vec, level, flag );
}

inline void createFunctionFromVector( const P2P2StokesFunction< PetscScalar >& function,
                                      const P2P2StokesFunction< PetscInt >&    numerator,
                                      Vec&                                     vec,
                                      uint_t                                   level,
                                      DoFType                                  flag )
{
   createFunctionFromVector( function.u, numerator.u, vec, level, flag );
   createFunctionFromVector( function.v, numerator.v, vec, level, flag );
   if ( function.u.getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented for 3D" );
      // createFunctionFromVector( function.w, numerator.w, vec, level, flag );
   }
   createFunctionFromVector( function.p, numerator.p, vec, level, flag );
}

inline void applyDirichletBC( const P2P2StokesFunction< PetscInt >& numerator, std::vector< PetscInt >& mat, uint_t level )
{
   applyDirichletBC( numerator.u, mat, level );
   applyDirichletBC( numerator.v, mat, level );
   if ( numerator.u.getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented for 3D" );
      // applyDirichletBC( numerator.w, mat, level );
   }
   //  applyDirichletBC(numerator.p, mat, level);
}

template < class OperatorType >
inline void createMatrix( const OperatorType&                   opr,
                          const P2P2StokesFunction< PetscInt >& src,
                          const P2P2StokesFunction< PetscInt >& dst,
                          Mat&                                  mat,
                          size_t                                level,
                          DoFType                               flag )
{
   createMatrix( opr.A, src.u, dst.u, mat, level, flag );
   createMatrix( opr.divT_x, src.p, dst.u, mat, level, flag );

   createMatrix( opr.A, src.v, dst.v, mat, level, flag );
   createMatrix( opr.divT_y, src.p, dst.v, mat, level, flag );

   if ( src.u.getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented for 3D" );
      // createMatrix( opr.A, src.w, dst.w, mat, level, flag );
      // createMatrix( opr.divT_z, src.p, dst.w, mat, level, flag );
   }

   createMatrix( opr.div_x, src.u, dst.p, mat, level, flag | DirichletBoundary );
   createMatrix( opr.div_y, src.v, dst.p, mat, level, flag | DirichletBoundary );
   if ( src.u.getStorage()->hasGlobalCells() )
   {
      // createMatrix( opr.div_z, src.w, dst.p, mat, level, flag | DirichletBoundary );
   }
   createMatrix( opr.pspg, src.p, dst.p, mat, level, flag | DirichletBoundary );
}

} // namespace petsc
} // namespace hhg