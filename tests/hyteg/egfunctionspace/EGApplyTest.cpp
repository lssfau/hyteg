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
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/dg1functionspace/DG1Operator.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

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

   EGLaplaceOperator L( storage, level, level );

   // PETSc apply
   PETScVector< real_t, EGFunction >      srcPetscVec;
   PETScVector< real_t, EGFunction >      dstPetscVec;
   PETScSparseMatrix< EGLaplaceOperator > L_Matrix;
   PETScSparseMatrix< EGMassOperator >    M_Matrix;

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

void EGApplyTestDiv( const std::string& testName, uint_t level, const MeshInfo& meshInfo, uint_t source, bool writeVTK = false )
{
   using namespace dg::eg;

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGFunction< real_t > src( "src", storage, level, level, BoundaryCondition::create0123BC() );
   P0Function< real_t > hytegDst( "hytegDst", storage, level, level );
   P0Function< real_t > petscDst( "petscDst", storage, level, level );
   P0Function< real_t > err( "error", storage, level, level );
   EGFunction< idx_t >  numeratorP1Vec( "numeratorP1Vec", storage, level, level, BoundaryCondition::create0123BC() );
   P0Function< idx_t >  numeratorP0( "numeratorP0", storage, level, level );

   if ( source == 0 )
   {
      src.interpolate( 1, level, Inner );

      src.getDiscontinuousPart()->interpolate( 1, level, All );
   }
   else
   {
      if ( storage->hasGlobalCells() )
      {
         src.interpolate( { []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                            []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                            []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); } },
                          level,
                          All );
      }
      else
      {
         src.interpolate( { []( const hyteg::Point3D& xx ) { return real_c( 20 ) * xx[0] * std::pow( xx[1], 3.0 ); },
                            []( const hyteg::Point3D& xx ) {
                               return real_c( 5 ) * std::pow( xx[0], 4.0 ) - real_c( 5 ) * std::pow( xx[1], 4.0 );
                            } },
                          level,
                          All );
      }
   }
   EGToP0DivOperator Div( storage, level, level );

   // PETSc apply
   PETScVector< real_t, EGFunction >      srcPetscVec;
   PETScVector< real_t, P0Function >      dstPetscVec;
   PETScSparseMatrix< EGToP0DivOperator > DivMatrix;

   numeratorP1Vec.enumerate( level );
   numeratorP0.enumerate( level );

   srcPetscVec.createVectorFromFunction( src, numeratorP1Vec, level );
   dstPetscVec.createVectorFromFunction( petscDst, numeratorP0, level );

   DivMatrix.createMatrixFromOperator( Div, level, numeratorP1Vec, numeratorP0 );
   Div.apply( src, hytegDst, level, All, Replace );

   MatMult( DivMatrix.get(), srcPetscVec.get(), dstPetscVec.get() );

   dstPetscVec.createFunctionFromVector( petscDst, numeratorP0, level );

   // compare
   err.assign( { 1.0, -1.0 }, { hytegDst, petscDst }, level, All );

   if ( writeVTK )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing vtk..." );
      VTKOutput vtk( "/mnt/c/Users/Fabia/OneDrive/Desktop/hyteg_premerge/hyteg-build/output", testName, storage );
      vtk.add( src );
      vtk.add( *src.getConformingPart() );
      vtk.add( *src.getDiscontinuousPart() );

      vtk.add( hytegDst );
      vtk.add( petscDst );
      vtk.add( err );
      vtk.write( level, 0 );
   }
   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||e|| ="
                                       << sqrt( err.dotGlobal( err, level ) ) / real_c( numberOfGlobalDoFs( err, level ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||hytegDst|| ="
                                       << sqrt( hytegDst.dotGlobal( hytegDst, level ) ) /
                                              real_c( numberOfGlobalDoFs( hytegDst, level ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||petscDst|| ="
                                       << sqrt( petscDst.dotGlobal( petscDst, level ) ) /
                                              real_c( numberOfGlobalDoFs( petscDst, level ) ) );
}

void EGApplyTestDivt(  const std::string& testName, uint_t level, const MeshInfo& meshInfo, bool writeVTK = false  )
{
   using namespace dg::eg;

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   P0Function< real_t >       src( "src", storage, level, level );
   EGFunction< real_t > hytegDst( "hytegDst", storage, level, level, BoundaryCondition::create0123BC() );
   EGFunction< real_t > petscDst( "petscDst", storage, level, level, BoundaryCondition::create0123BC() );
   EGFunction< real_t > err( "error", storage, level, level, BoundaryCondition::create0123BC() );
   EGFunction< idx_t >  numeratorP1Vec( "numeratorP1Vec", storage, level, level, BoundaryCondition::create0123BC() );
   P0Function< idx_t >        numeratorP0( "numeratorP0", storage, level, level );

   P0ToEGDivTOperator Divt( storage, level, level );

   // PETSc apply
   PETScVector< real_t, P0Function >                  srcPetscVec;
   PETScVector< real_t, EGFunction >            dstPetscVec;
   PETScSparseMatrix< P0ToEGDivTOperator > Divt_Matrix;

   src.interpolate( 1, level, All );

   numeratorP1Vec.enumerate( level );
   numeratorP0.enumerate( level );

   srcPetscVec.createVectorFromFunction( src, numeratorP0, level );
   dstPetscVec.createVectorFromFunction( petscDst, numeratorP1Vec, level );

   Divt_Matrix.createMatrixFromOperator( Divt, level, numeratorP0, numeratorP1Vec );
   Divt.apply( src, hytegDst, level, All, Replace );

   // WALBERLA_CHECK( L_Matrix.isSymmetric() );

   MatMult( Divt_Matrix.get(), srcPetscVec.get(), dstPetscVec.get() );

   dstPetscVec.createFunctionFromVector( petscDst, numeratorP1Vec, level );

   // compare
   err.assign( { 1.0, -1.0 }, { hytegDst, petscDst }, level, All );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output", testName, storage );

      vtk.add( src );
      vtk.add( hytegDst );
      vtk.add( petscDst );
      vtk.add( err );
      vtk.write( level, 0 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||e||=" << sqrt( err.dotGlobal( err, level ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||hytegDst|| ="
                                       << sqrt( hytegDst.dotGlobal( hytegDst, level ) ) /
                                              real_c( numberOfGlobalDoFs( hytegDst, level ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( testName << ": ||petscDst|| ="
                                       << sqrt( petscDst.dotGlobal( petscDst, level ) ) /
                                              real_c( numberOfGlobalDoFs( petscDst, level ) ) );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   using hyteg::Point3D;
   using walberla::math::pi;
   const bool writeVTK = true;

   ScalarLambda srcLambda1 = []( const hyteg::Point3D& x ) { return std::sin( 3 * pi * x[0] ) * std::sin( 3 * pi * x[1] ); };
   ScalarLambda srcLambda2 = []( const hyteg::Point3D& ) { return 1; };
   ScalarLambda srcLambda3 = []( const hyteg::Point3D& x ) { return x[0] * x[0] * x[0] * std::sin( 3 * pi * x[1] ); };

   EGApplyTestDivt(
       "EGApplyDivT3D", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), true );

   EGApplyTestDivt( "EGApplyDivT2D", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), true );

   EGApplyTestDiv(
       "EGApplyDiv3D_src0", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ),  0,true );

   EGApplyTestDiv( "EGApplyDiv2D_src0", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 0,true );

   EGApplyTestDiv(
       "EGApplyDiv3D_src1", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ),  1,true );

   EGApplyTestDiv( "EGApplyDiv2D_src1", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1,true );
   /*
   hyteg::EGApplyTest(
       srcLambda1, "tri_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "tri_2el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "quad_4el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   hyteg::EGApplyTest(
       srcLambda2, "tri_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tri_2el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "quad_4el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   hyteg::EGApplyTest(
       srcLambda3, "tri_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tri_2el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "quad_4el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   hyteg::EGApplyTest(
       srcLambda1, "tet_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tet_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tet_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   
   hyteg::EGApplyTest( srcLambda1,
                       "pyramid_2el_3_src1",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda2,
                       "pyramid_2el_3_src2",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda3,
                       "pyramid_2el_3_src3",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       1.0e-15,
                       writeVTK );
   
   hyteg::EGApplyTest( srcLambda1,
                       "pyramid_4el_3_src1",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda2,
                       "pyramid_4el_3_src2",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda3,
                       "pyramid_4el_3_src3",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       1.0e-15,
                       writeVTK );
*/
   return EXIT_SUCCESS;
}
