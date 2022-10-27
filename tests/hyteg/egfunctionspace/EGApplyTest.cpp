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
   numerator.enumerate( level );

   EGLaplaceOperator L( storage, level, level );
   EGMassOperator    M( storage, level, level );

   // PETSc apply
   PETScVector< real_t, EGFunction >      srcPetscVec;
   PETScVector< real_t, EGFunction >      dstPetscVec;
   PETScSparseMatrix< EGLaplaceOperator > L_Matrix;
   PETScSparseMatrix< EGMassOperator >    M_Matrix;

   std::function< real_t( const hyteg::Point3D& ) > srcFunction = srcLambda;
   src.interpolate( srcFunction, level, All );

   srcPetscVec.createVectorFromFunction( src, numerator, level );
   dstPetscVec.createVectorFromFunction( petscDst, numerator, level );
   L_Matrix.createMatrixFromOperator( L, level, numerator );
   M_Matrix.createMatrixFromOperator( M, level, numerator );

   // L_Matrix.print( "EGApplyTest_L.m", false, PETSC_VIEWER_ASCII_MATLAB );
   //  M_Matrix.print( "EGApplyTest_M.m", false, PETSC_VIEWER_ASCII_MATLAB );

   L.apply( src, hytegDst, level, All, Replace );

   // WALBERLA_CHECK( petscMatrix.isSymmetric() );

   MatMult( L_Matrix.get(), srcPetscVec.get(), dstPetscVec.get() );

   dstPetscVec.createFunctionFromVector( petscDst, numerator, level );

   // compare
   err.assign( { 1.0, -1.0 }, { hytegDst, petscDst }, level, Inner );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output", testName, storage );
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

   using walberla::math::pi;
   const bool   writeVTK   = true;
   ScalarLambda srcLambda1 = []( const hyteg::Point3D& x ) { return std::sin( 3 * pi * x[0] ) * std::sin( 3 * pi * x[1] ); };
   hyteg::EGApplyTest(
       srcLambda1, "tet_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "tri_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "tri_2el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "quad_4el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   ScalarLambda srcLambda2 = []( const hyteg::Point3D& x ) { return 1; };
   hyteg::EGApplyTest(
       srcLambda2, "tet_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tri_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tri_2el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "quad_4el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   ScalarLambda srcLambda3 = []( const hyteg::Point3D& x ) { return x[0] * x[0] * x[0] * std::sin( 3 * pi * x[1] ); };
   hyteg::EGApplyTest(
       srcLambda3, "tet_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tri_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tri_2el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "quad_4el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   /* 
   hyteg::EGApplyTest( "pyramid_2el_2", 2, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ), 1.0e-15 );
   hyteg::EGApplyTest( "pyramid_2el_4", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ), 1.0e-15 );
   hyteg::EGApplyTest( "pyramid_4el_2", 2, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ), 1.0e-15 );
   hyteg::EGApplyTest( "pyramid_4el_4", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ), 1.0e-15 );
*/
   return EXIT_SUCCESS;
}
