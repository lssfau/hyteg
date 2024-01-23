/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constantStencilOperator/P1ConstantOperator.cpp"
#include "mixedOperator/EGOperators.hpp"
#include "mixedOperator/EGOperatorsNitscheBC.hpp"
#include "mixedOperator/P2P1TaylorHoodStokesOperator.hpp"

using walberla::real_t;
typedef std::function< real_t( const hyteg::PointND< real_t, 3 >& p ) > ScalarLambda;
namespace hyteg {

void EGApplyTest( ScalarLambda       srcLambda,
                  const std::string& testName,
                  uint_t             level,
                  const MeshInfo&    meshInfo,
                  real_t             eps,
                  bool               writeVTK = false )
{
   using namespace dg::eg;

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGFunction< real_t > src( "src", storage, level, level );
   EGFunction< real_t > tmp( "tmp", storage, level, level );
   EGFunction< real_t > hytegDst( "hytegDst", storage, level, level );
   EGFunction< real_t > petscDst( "petscDst", storage, level, level );
   EGFunction< real_t > err( "error", storage, level, level );
   EGFunction< idx_t >  numerator( "numerator", storage, level, level );

   EGSIPGLaplaceOperator L( storage, level, level );

   // PETSc apply
   PETScVector< real_t, EGFunction >          srcPetscVec;
   PETScVector< real_t, EGFunction >          dstPetscVec;
   PETScSparseMatrix< EGSIPGLaplaceOperator > L_Matrix;
   PETScSparseMatrix< EGMassOperator >        M_Matrix;

   std::function< real_t( const hyteg::Point3D& ) > srcFunction = srcLambda;
   src.interpolate( { srcFunction, srcFunction, srcFunction }, level, All );
   src.getDiscontinuousPart()->interpolate( srcFunction, level, All );

   numerator.copyBoundaryConditionFromFunction( src );
   numerator.enumerate( level );

   srcPetscVec.createVectorFromFunction( src, numerator, level );
   dstPetscVec.createVectorFromFunction( petscDst, numerator, level );
   L_Matrix.createMatrixFromOperator( L, level, numerator );
   // L_Matrix.print( "EGApplyTest_L.m", false, PETSC_VIEWER_ASCII_MATLAB );
   //  M_Matrix.print( "EGApplyTest_M.m", false, PETSC_VIEWER_ASCII_MATLAB );
   L.apply( src, hytegDst, level, All, Replace );

   // WALBERLA_CHECK( L_Matrix.isSymmetric() );

   MatMult( L_Matrix.get(), srcPetscVec.get(), dstPetscVec.get() );

   dstPetscVec.createFunctionFromVector( petscDst, numerator, level );

   // compare
   err.assign( { 1.0, -1.0 }, { hytegDst, petscDst }, level, All );
   WALBERLA_LOG_INFO_ON_ROOT( "||e_disc|| = "
                              << sqrt( err.getDiscontinuousPart()->dotGlobal( *err.getDiscontinuousPart(), level, All ) /
                                       real_c( numberOfGlobalDoFs( *err.getDiscontinuousPart(), level ) ) )
                              << ", ||e_conf|| = "
                              << sqrt( err.getConformingPart()->dotGlobal( *err.getConformingPart(), level, All ) /
                                       real_c( numberOfGlobalDoFs( *err.getConformingPart(), level ) ) ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output", testName, storage );

      vtk.add( src );
      vtk.add( *src.getConformingPart() );
      vtk.add( *src.getDiscontinuousPart() );
      vtk.add( hytegDst );
      vtk.add( *hytegDst.getConformingPart() );
      vtk.add( *hytegDst.getDiscontinuousPart() );
      vtk.add( petscDst );
      vtk.add( *petscDst.getConformingPart() );
      vtk.add( *petscDst.getDiscontinuousPart() );
      vtk.add( err );
      vtk.add( *err.getConformingPart() );
      vtk.add( *err.getDiscontinuousPart() );
      vtk.write( level, 0 );
   }

   auto maxMag = err.getMaxMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||e||_max = " << maxMag );

   WALBERLA_CHECK_LESS( maxMag, eps );
}

template < typename EGOperatorType >
void EGApplyNitscheBCTest( ScalarLambda       srcLambda,
                           const std::string& testName,
                           uint_t             level,
                           const MeshInfo&    meshInfo,
                           real_t             eps,
                           bool               writeVTK = false )
{
   using namespace dg::eg;

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGFunction< real_t > src( "src", storage, level, level );
   EGFunction< real_t > tmp( "tmp", storage, level, level );
   EGFunction< real_t > hytegDst( "hytegDst", storage, level, level );
   EGFunction< real_t > petscDst( "petscDst", storage, level, level );
   EGFunction< real_t > err( "error", storage, level, level );
   EGFunction< idx_t >  numerator( "numerator", storage, level, level );

   // PETSc apply
   PETScVector< real_t, EGFunction >   srcPetscVec;
   PETScVector< real_t, EGFunction >   dstPetscVec;
   PETScSparseMatrix< EGOperatorType > L_Matrix;
   PETScSparseMatrix< EGMassOperator > M_Matrix;
   std::shared_ptr< EGOperatorType >   L;

   std::function< real_t( const hyteg::Point3D& ) > srcFunction = srcLambda;
   src.interpolate( { srcFunction, srcFunction, srcFunction }, level, All );
   src.getDiscontinuousPart()->interpolate( srcFunction, level, All );

   numerator.copyBoundaryConditionFromFunction( src );
   numerator.enumerate( level );

   srcPetscVec.createVectorFromFunction( src, numerator, level );
   dstPetscVec.createVectorFromFunction( petscDst, numerator, level );

   if constexpr ( std::is_same< EGOperatorType, EGEpsilonOperatorNitscheBC >::value ||
                  std::is_same< EGOperatorType, EGEpsilonOperator >::value )
   {
      L = std::make_shared< EGOperatorType >( storage, level, level, []( const Point3D& ) -> real_t { return 1; } );
   }
   else
   {
      L = std::make_shared< EGOperatorType >( storage, level, level );
   }
   L_Matrix.createMatrixFromOperator( *L, level, numerator );
   L->apply( src, hytegDst, level, All, Replace );
   MatMult( L_Matrix.get(), srcPetscVec.get(), dstPetscVec.get() );
   // L_Matrix.print( "EGApplyTest_L.m", false, PETSC_VIEWER_ASCII_MATLAB );
   //  M_Matrix.print( "EGApplyTest_M.m", false, PETSC_VIEWER_ASCII_MATLAB );

   // WALBERLA_CHECK( L_Matrix.isSymmetric() );

   dstPetscVec.createFunctionFromVector( petscDst, numerator, level );

   // compare
   err.assign( { 1.0, -1.0 }, { hytegDst, petscDst }, level, All );
   WALBERLA_LOG_INFO_ON_ROOT( "||e_disc|| = "
                              << sqrt( err.getDiscontinuousPart()->dotGlobal( *err.getDiscontinuousPart(), level, All ) /
                                       real_c( numberOfGlobalDoFs( *err.getDiscontinuousPart(), level ) ) )
                              << ", ||e_conf|| = "
                              << sqrt( err.getConformingPart()->dotGlobal( *err.getConformingPart(), level, All ) /
                                       real_c( numberOfGlobalDoFs( *err.getConformingPart(), level ) ) )
                              << ", ||e_conf_inner|| = "
                              << sqrt( err.getConformingPart()->dotGlobal( *err.getConformingPart(), level, Inner ) /
                                       real_c( numberOfGlobalDoFs( *err.getConformingPart(), level ) ) ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output", testName, storage );

      vtk.add( src );
      vtk.add( *src.getConformingPart() );
      vtk.add( *src.getDiscontinuousPart() );
      vtk.add( hytegDst );
      vtk.add( *hytegDst.getConformingPart() );
      vtk.add( *hytegDst.getDiscontinuousPart() );
      vtk.add( petscDst );
      vtk.add( *petscDst.getConformingPart() );
      vtk.add( *petscDst.getDiscontinuousPart() );
      vtk.add( err );
      vtk.add( *err.getConformingPart() );
      vtk.add( *err.getDiscontinuousPart() );
      vtk.write( level, 0 );
   }

   auto maxMag = err.getMaxMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||e||_max = " << maxMag );

   WALBERLA_CHECK_LESS( maxMag, eps );
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   using namespace hyteg::dg::eg;
   using hyteg::Point3D;
   using walberla::math::pi;
   const bool writeVTK = false;
   real_t     eps      = 1e-13;

   ScalarLambda srcLambda1 = []( const hyteg::Point3D& x ) { return std::sin( 3 * pi * x[0] ) * std::sin( 3 * pi * x[1] ); };
   ScalarLambda srcLambda2 = []( const hyteg::Point3D& ) { return 1; };
   ScalarLambda srcLambda3 = []( const hyteg::Point3D& x ) { return x[0] * x[0] * x[0] * std::sin( 3 * pi * x[1] ); };
   ScalarLambda srcLambda4 = []( const hyteg::Point3D& x ) {
      return std::sin( 3 * pi * x[0] ) * std::sin( 3 * pi * x[1] ) * std::sin( 3 * pi * x[2] );
   };
   ScalarLambda srcLambda5 = []( const hyteg::Point3D& x ) {
      return x[0] * x[0] * x[0] * std::sin( 3 * pi * x[1] ) * std::sin( 3 * pi * x[2] );
   };

   hyteg::EGApplyNitscheBCTest< EGLaplaceOperatorNitscheBC >( srcLambda1,
                                                              "EGLaplaceOperatorNitscheBC_tri_1el_4_src1",
                                                              4,
                                                              hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ),
                                                              eps,
                                                              writeVTK );

   hyteg::EGApplyNitscheBCTest< EGLaplaceOperatorNitscheBC >( srcLambda3,
                                                              "EGLaplaceOperatorNitscheBC_tri_1el_4_src3",
                                                              4,
                                                              hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ),
                                                              eps,
                                                              writeVTK );

   hyteg::EGApplyNitscheBCTest< EGLaplaceOperatorNitscheBC >( srcLambda1,
                                                              "EGLaplaceOperatorNitscheBC_quad_4el_4_src1",
                                                              4,
                                                              hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ),
                                                              eps,
                                                              writeVTK );

   hyteg::EGApplyNitscheBCTest< EGLaplaceOperatorNitscheBC >( srcLambda3,
                                                              "EGLaplaceOperatorNitscheBC_quad_4el_4_src3",
                                                              4,
                                                              hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ),
                                                              eps,
                                                              writeVTK );

   hyteg::EGApplyNitscheBCTest< EGLaplaceOperatorNitscheBC >(
       srcLambda4,
       "EGLaplaceOperatorNitscheBC_cube_6el_3_src4",
       3,
       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ),
       eps,
       writeVTK );

   hyteg::EGApplyNitscheBCTest< EGLaplaceOperatorNitscheBC >(
       srcLambda5,
       "EGLaplaceOperatorNitscheBC_cube_6el_3_src5",
       3,
       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ),
       eps,
       writeVTK );

   hyteg::EGApplyNitscheBCTest< EGEpsilonOperatorNitscheBC >( srcLambda1,
                                                              "EGEpsilonOperatorNitscheBC_quad_4el_3_src1",
                                                              4,
                                                              hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ),
                                                              eps,
                                                              writeVTK );

   hyteg::EGApplyNitscheBCTest< EGEpsilonOperatorNitscheBC >( srcLambda3,
                                                              "EGEpsilonOperatorNitscheBC_quad_4el_3_src3",
                                                              4,
                                                              hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ),
                                                              eps,
                                                              writeVTK );

   hyteg::EGApplyNitscheBCTest< EGEpsilonOperatorNitscheBC >(
       srcLambda4,
       "EGEpsilonOperatorNitscheBC_cube_6el_3_src4",
       3,
       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ),
       eps,
       writeVTK );

   hyteg::EGApplyNitscheBCTest< EGEpsilonOperatorNitscheBC >(
       srcLambda5,
       "EGEpsilonOperatorNitscheBC_cube_6el_3_src5",
       3,
       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ),
       eps,
       writeVTK );

   hyteg::EGApplyTest(
       srcLambda1, "tri_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "tri_2el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "quad_4el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), eps, writeVTK );

   hyteg::EGApplyTest(
       srcLambda2, "tri_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tri_2el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "quad_4el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), eps, writeVTK );

   hyteg::EGApplyTest(
       srcLambda3, "tri_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tri_2el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "quad_4el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), eps, writeVTK );

   hyteg::EGApplyTest(
       srcLambda1, "tet_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tet_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), eps, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tet_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), eps, writeVTK );

   hyteg::EGApplyTest( srcLambda1,
                       "pyramid_2el_3_src1",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       2.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda2,
                       "pyramid_2el_3_src2",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       2.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda3,
                       "pyramid_2el_3_src3",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       2.0e-15,
                       writeVTK );

   hyteg::EGApplyTest( srcLambda1,
                       "pyramid_4el_3_src1",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       2.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda2,
                       "pyramid_4el_3_src2",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       2.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda3,
                       "pyramid_4el_3_src3",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       2.0e-15,
                       writeVTK );

   /*
      * // Test divT operator
     EGApplyTestDivt(
         "EGApplyDivT3D", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), true );

    EGApplyTestDivt( "EGApplyDivT2D", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), true );


    // test div operator
    EGApplyTestDiv("EGApplyDiv_2D_Ones", 3, hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh"), 0, true);
    EGApplyTestDiv("EGApplyDiv_2D_AllOnes", 3, hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh"), 1, true);
    EGApplyTestDiv(
        "EGApplyDiv_2D_Sinusoidal", 3,
        hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh"), 2, true);

    EGApplyTestDiv(
        "EGApplyDiv_3D_Ones", 3,
        hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_center_at_origin_24el.msh"), 0, true);
    EGApplyTestDiv(
        "EGApplyDiv_3D_AllOnes", 3,
        hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_center_at_origin_24el.msh"), 1, true);
    EGApplyTestDiv(
        "EGApplyDiv_3D_Sinusoidal", 3,
        hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_center_at_origin_24el.msh"), 2, true);
    */

   return EXIT_SUCCESS;
}
