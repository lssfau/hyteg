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

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_affine_q2.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

void testDirichletBC()
{
   using namespace dg;

   auto f = []( const Point3D& x ) { return 1 - x[0] - x[1]; };

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >( f );

   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsElement( 3 );
   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsFacet0( 2 );
   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsFacet1( 2 );
   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsFacet2( 2 );
   Eigen::Matrix< real_t, 3, 1 >                normal0, normal1, normal2;

   coordsElement[0]( 0 ) = 0;
   coordsElement[0]( 1 ) = 0;
   coordsElement[0]( 2 ) = 0;

   coordsElement[1]( 0 ) = 2;
   coordsElement[1]( 1 ) = 0;
   coordsElement[1]( 2 ) = 0;

   coordsElement[2]( 0 ) = 0;
   coordsElement[2]( 1 ) = 2;
   coordsElement[2]( 2 ) = 0;

   coordsFacet0[0] = coordsElement[0];
   coordsFacet0[1] = coordsElement[1];

   coordsFacet1[0] = coordsElement[0];
   coordsFacet1[1] = coordsElement[2];

   coordsFacet2[0] = coordsElement[1];
   coordsFacet2[1] = coordsElement[2];

   normal0( 0 ) = 0;
   normal0( 1 ) = -1;
   normal0( 2 ) = 0;

   normal1( 0 ) = -1;
   normal1( 1 ) = 0;
   normal1( 2 ) = 0;

   normal2( 0 ) = 1;
   normal2( 1 ) = 1;
   normal2( 2 ) = 0;
   normal2.normalize();

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > mat;
   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > tmpMat;
   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > rhs;
   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > tmpRhs;
   mat.resize( 3, 3 );
   tmpMat.resize( 3, 3 );
   rhs.resize( 3, 1 );
   tmpRhs.resize( 3, 1 );

   mat.setZero();
   tmpMat.setZero();
   rhs.setZero();
   tmpRhs.setZero();

   laplaceForm->integrateVolume( 2, coordsElement, *basis, *basis, 1, 1, tmpMat );
   mat += tmpMat;
   tmpMat.setZero();

   laplaceForm->integrateFacetDirichletBoundary(
       2, coordsElement, coordsFacet0, coordsElement[2], normal0, *basis, *basis, 1, 1, tmpMat );
   mat += tmpMat;
   tmpMat.setZero();
   laplaceForm->integrateFacetDirichletBoundary(
       2, coordsElement, coordsFacet1, coordsElement[1], normal1, *basis, *basis, 1, 1, tmpMat );
   mat += tmpMat;
   tmpMat.setZero();
   laplaceForm->integrateFacetDirichletBoundary(
       2, coordsElement, coordsFacet2, coordsElement[0], normal2, *basis, *basis, 1, 1, tmpMat );
   mat += tmpMat;

   laplaceForm->integrateRHSDirichletBoundary( 2, coordsElement, coordsFacet0, coordsElement[2], normal0, *basis, 1, tmpRhs );
   rhs += tmpRhs;
   tmpRhs.setZero();
   laplaceForm->integrateRHSDirichletBoundary( 2, coordsElement, coordsFacet1, coordsElement[1], normal1, *basis, 1, tmpRhs );
   rhs += tmpRhs;
   tmpRhs.setZero();
   laplaceForm->integrateRHSDirichletBoundary( 2, coordsElement, coordsFacet2, coordsElement[0], normal2, *basis, 1, tmpRhs );
   rhs += tmpRhs;

   WALBERLA_LOG_INFO_ON_ROOT( mat );
   WALBERLA_LOG_INFO_ON_ROOT( rhs );

   WALBERLA_LOG_INFO_ON_ROOT( mat.inverse() * rhs );

} // namespace hyteg

void testDiffForm()
{
   using namespace dg;

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >();

   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsElement;
   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsNeighborElement;
   std::vector< Eigen::Matrix< real_t, 3, 1 > > coordsFacet;
   Eigen::Matrix< real_t, 3, 1 >                normal;

   coordsElement[0]( 0 ) = 0;
   coordsElement[0]( 1 ) = 0;

   coordsElement[1]( 0 ) = 1;
   coordsElement[1]( 1 ) = 0;

   coordsElement[2]( 0 ) = 0;
   coordsElement[2]( 1 ) = 1;

   coordsNeighborElement[0]( 0 ) = 1;
   coordsNeighborElement[0]( 1 ) = -1;

   coordsNeighborElement[1]( 0 ) = 1;
   coordsNeighborElement[1]( 1 ) = 0;

   coordsNeighborElement[2]( 0 ) = 0;
   coordsNeighborElement[2]( 1 ) = 0;

   coordsFacet[0] = coordsElement[0];
   coordsFacet[1] = coordsElement[1];

   normal( 0 ) = 0;
   normal( 1 ) = -1;

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > elMat;
   elMat.resize( 3, 3 );

   laplaceForm->integrateFacetCoupling( 2,
                                        coordsElement,
                                        coordsNeighborElement,
                                        coordsFacet,
                                        coordsElement[2],
                                        coordsNeighborElement[0],
                                        normal,
                                        *basis,
                                        *basis,
                                        1,
                                        1,
                                        elMat );

   WALBERLA_LOG_INFO_ON_ROOT( elMat );

   elMat.setZero();
   laplaceForm->integrateFacetCoupling( 2,
                                        coordsNeighborElement,
                                        coordsElement,
                                        coordsFacet,
                                        coordsNeighborElement[0],
                                        coordsElement[2],
                                        -normal,
                                        *basis,
                                        *basis,
                                        1,
                                        1,
                                        elMat );

   WALBERLA_LOG_INFO_ON_ROOT( elMat );

   elMat.setZero();
   laplaceForm->integrateFacetInner( 2,
                                     coordsElement,
                                     { coordsElement[1], coordsElement[2] },
                                     coordsElement[0],
                                     Eigen::Matrix< real_t, 3, 1 >( 1, 1 ).normalized(),
                                     *basis,
                                     *basis,
                                     1,
                                     1,
                                     elMat );

   WALBERLA_LOG_INFO_ON_ROOT( elMat );
}

void testDiffOp( uint_t maxLevel )
{
   using namespace dg;

   // MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );
   MeshInfo meshInfo = MeshInfo::meshFaceChain( 2 );

//   auto ff = []( const Point3D& p ) {
//      Matrix3r mat;
//      mat( 0, 0 ) = .34;
//      mat( 0, 1 ) = .678;
//      mat( 0, 2 ) = 0;
//
//      mat( 1, 0 ) = .282;
//      mat( 1, 1 ) = .35353;
//      mat( 1, 2 ) = 0;
//
//      mat( 2, 0 ) = 0;
//      mat( 2, 1 ) = 0;
//      mat( 2, 2 ) = 0;
//
//      return mat.mul( p );
//   };

   // meshInfo.applyCoordinateMap( ff );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   const uint_t minLevel = 2;

   std::function< real_t( const Point3D& ) > linearFunc0  = []( const Point3D& ) { return 0; };
   std::function< real_t( const Point3D& ) > linearFuncX  = []( const Point3D& x ) { return x[0]; };
   std::function< real_t( const Point3D& ) > linearFuncY  = []( const Point3D& x ) { return x[1]; };
   std::function< real_t( const Point3D& ) > linearFuncXY = []( const Point3D& x ) { return 1 - x[0] - x[1]; };
   std::function< real_t( const Point3D& ) > linearFuncC  = []( const Point3D& ) { return 1; };

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >( linearFunc0 );
   auto massForm    = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > u( "u", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > f( "f", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > sol( "sol", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, basis, 1 );

   DGFunction< real_t > err_x( "err_x", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > err_y( "err_y", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > err_xy( "err_xy", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > err_c( "err_c", storage, minLevel, maxLevel, basis, 1 );

   DGFunction< idx_t > numerator( "numerator", storage, minLevel, maxLevel, basis, 1 );
   numerator.enumerate( maxLevel );

   DGOperator A( storage, minLevel, maxLevel, laplaceForm );
   DGOperator M( storage, minLevel, maxLevel, massForm );

   std::string                     fileName = "../../output/diffusion.m";
   PETScSparseMatrix< DGOperator > mat;
   mat.createMatrixFromOperator( A, maxLevel, numerator, numerator );
   mat.print( fileName, false, PETSC_VIEWER_ASCII_MATLAB );

   tmp.evaluateLinearFunctional( linearFuncX, maxLevel );
   PETScCGSolver< DGOperator > solverM( storage, maxLevel, numerator );
   solverM.solve( M, sol, tmp, maxLevel );
   A.apply( sol, err_x, maxLevel, All, Replace );
   auto discrL2 = sqrt( err_x.dotGlobal( err_x, maxLevel ) / real_c( numberOfGlobalDoFs( u, maxLevel ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "x: level " << maxLevel << ": L2 error " << discrL2 );

   tmp.evaluateLinearFunctional( linearFuncY, maxLevel );
   solverM.solve( M, sol, tmp, maxLevel );
   A.apply( sol, err_y, maxLevel, All, Replace );
   discrL2 = sqrt( err_y.dotGlobal( err_y, maxLevel ) / real_c( numberOfGlobalDoFs( u, maxLevel ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "y: level " << maxLevel << ": L2 error " << discrL2 );

   tmp.evaluateLinearFunctional( linearFuncXY, maxLevel );
   solverM.solve( M, sol, tmp, maxLevel );
   A.apply( sol, err_xy, maxLevel, All, Replace );
   discrL2 = sqrt( err_xy.dotGlobal( err_xy, maxLevel ) / real_c( numberOfGlobalDoFs( u, maxLevel ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "xy: level " << maxLevel << ": L2 error " << discrL2 );

   tmp.evaluateLinearFunctional( linearFuncC, maxLevel );
   solverM.solve( M, sol, tmp, maxLevel );
   A.apply( sol, err_c, maxLevel, All, Replace );
   discrL2 = sqrt( err_c.dotGlobal( err_c, maxLevel ) / real_c( numberOfGlobalDoFs( u, maxLevel ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "c: level " << maxLevel << ": L2 error " << discrL2 );

   writeDomainPartitioningVTK( *storage, "../../output", "DGSmokeTestDiffOp_Domain" );

   VTKOutput vtkOutput( "../../output", "DGSmokeTestDiffOp", storage );

   vtkOutput.add( u );
   vtkOutput.add( err_x );
   vtkOutput.add( err_y );
   vtkOutput.add( err_xy );
   vtkOutput.add( err_c );
   vtkOutput.add( sol );

   vtkOutput.write( maxLevel );
}

void testMass()
{
   using namespace dg;

   // MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );
   MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );
   // MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 5;

   // std::function< real_t( const Point3D& ) > f = []( const Point3D& x ) { return sqrt( x[0] * x[0] + x[1] * x[1] ) < 0.3 ? 1 : 0; };
   // std::function< real_t( const Point3D& ) > f = []( const Point3D& x ) { return tanh( 8 * x[0] ); };
   std::function< real_t( const Point3D& ) > f = []( const Point3D& x ) { return sin( 2 * pi * x[0] ); };

   auto basis = std::make_shared< DGBasisLinearLagrange_Example >();

   DGFunction< real_t > u( "u", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, basis, 1 );

   tmp.evaluateLinearFunctional( f, maxLevel );

   auto       mass = std::make_shared< DGMassForm_Example >();
   DGOperator M( storage, minLevel, maxLevel, mass );

   std::string                     fileName = "/tmp/mass.m";
   PETScSparseMatrix< DGOperator > mat;
   DGFunction< idx_t >             numeratorSrc( "numeratorSrc", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< idx_t >             numeratorDst( "numeratorDst", storage, minLevel, maxLevel, basis, 1 );
   numeratorSrc.enumerate( maxLevel );
   numeratorDst.enumerate( maxLevel );
   mat.createMatrixFromOperator( M, maxLevel, numeratorSrc, numeratorDst );
   mat.print( fileName, false, PETSC_VIEWER_ASCII_MATLAB );

   PETScCGSolver< DGOperator > solver( storage, maxLevel, numeratorSrc );
   solver.solve( M, u, tmp, maxLevel );

   VTKOutput vtkOutput( "../../output", "DGSmokeTestMass", storage );

   vtkOutput.add( u );
   vtkOutput.add( tmp );

   vtkOutput.write( maxLevel );
}

void testDiffusion( uint_t maxLevel )
{
   using namespace dg;

   MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

//   auto ff = []( const Point3D& p ) {
//      Matrix3r mat;
//      mat( 0, 0 ) = .34;
//      mat( 0, 1 ) = .678;
//      mat( 0, 2 ) = 0;
//
//      mat( 1, 0 ) = .282;
//      mat( 1, 1 ) = .35353;
//      mat( 1, 2 ) = 0;
//
//      mat( 2, 0 ) = 0;
//      mat( 2, 1 ) = 0;
//      mat( 2, 2 ) = 0;
//
//      return mat.mul( p );
//   };

   // meshInfo.applyCoordinateMap( ff );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 2, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;

   std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
      return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
   };

   std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
      return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
   };

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >( solFunc );
   auto massForm    = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > u( "u", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > f( "f", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > sol( "sol", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > err( "err", storage, minLevel, maxLevel, basis, 1 );

   DGFunction< idx_t > numerator( "numerator", storage, minLevel, maxLevel, basis, 1 );
   numerator.enumerate( maxLevel );

   DGOperator A( storage, minLevel, maxLevel, laplaceForm );
   DGOperator M( storage, minLevel, maxLevel, massForm );

   std::string                     fileName = "/tmp/diffusion.m";
   PETScSparseMatrix< DGOperator > mat;
   mat.createMatrixFromOperator( A, maxLevel, numerator, numerator );
   mat.print( fileName, false, PETSC_VIEWER_ASCII_MATLAB );

   f.evaluateLinearFunctional( rhsFunc, maxLevel );
   tmp.evaluateLinearFunctional( solFunc, maxLevel );

   PETScCGSolver< DGOperator > solverM( storage, maxLevel, numerator );
   solverM.solve( M, sol, tmp, maxLevel );
   PETScCGSolver< DGOperator > solverA( storage, maxLevel, numerator, 1e-12, 1e-12, 10000 );
   solverA.solve( A, u, f, maxLevel );

   err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
   auto discrL2 = sqrt( err.dotGlobal( err, maxLevel ) / real_c( numberOfGlobalDoFs( u, maxLevel ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( "level " << maxLevel << ": L2 error " << discrL2 );

   VTKOutput vtkOutput( "../../output", "DGSmokeTestLaplace", storage );

   vtkOutput.add( u );
   vtkOutput.add( err );
   vtkOutput.add( sol );

   vtkOutput.write( maxLevel );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   // hyteg::testDiffForm();

   // hyteg::testMass();

   // hyteg::testDiffusion( 3 );
   // hyteg::testDiffusion( 4 );
   // hyteg::testDiffusion( 5 );
   // hyteg::testDiffusion( 6 );
#if 0
   hyteg::testDiffOp( 3 );
   // hyteg::testDiffOp( 4 );
   // hyteg::testDiffOp( 5 );
#endif

   hyteg::testDirichletBC();
   return EXIT_SUCCESS;
}
