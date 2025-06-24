/*
* Copyright (c) 2017-2025 Andreas Wagner, Marcus Mohr.
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

#include "hyteg/gridtransferoperators/DGProlongation.hpp"

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/forms/form_hyteg_dg/DGMassForm_Example.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/gridtransferoperators/DGRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;
using namespace dg;

void testDG1Prolongation( const std::string& meshFile )
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = minLevel + 1;

   const uint_t degree = 1;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto dgBasis  = std::make_shared< DGBasisLinearLagrange_Example >();
   auto massForm = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > test_function( "test_function", storage, minLevel, maxLevel, dgBasis, degree );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, dgBasis, degree );
   DGFunction< real_t > solution_function( "solution_function", storage, maxLevel, maxLevel, dgBasis, degree );

   DGFunction< idx_t > numerator( "numerator", storage, minLevel, maxLevel, dgBasis, degree );
   numerator.enumerate( minLevel );
   numerator.enumerate( maxLevel );
   DGOperator M( storage, minLevel, maxLevel, massForm );
   tmp.evaluateLinearFunctional( []( const Point3D& p ) { return 0.5 * p[0] + p[1] - 0.25 * p[2]; }, minLevel );
   tmp.evaluateLinearFunctional( []( const Point3D& p ) { return 0.5 * p[0] + p[1] - 0.25 * p[2]; }, maxLevel );
   // interpolate test_function on lower level
   {
      PETScCGSolver< DGOperator > solverM( storage, minLevel, numerator );
      solverM.solve( M, test_function, tmp, minLevel );
   }
   // interpolate solution function on upper level
   {
      PETScCGSolver< DGOperator > solverM( storage, maxLevel, numerator );
      solverM.solve( M, solution_function, tmp, maxLevel );
   }

   DG1Prolongation prolongation;

   prolongation.prolongate( test_function, minLevel, All );

   tmp.assign( { +1., -1. }, { test_function, solution_function }, maxLevel, All );
   const real_t error = tmp.dotGlobal( tmp, maxLevel, All );

   WALBERLA_CHECK_FLOAT_EQUAL( error, 0. );

   // VTKOutput vtk( "../../output/", "DGProlongationTest", storage );
   // vtk.add( test_function );
   // vtk.write( minLevel );
   // vtk.write( maxLevel );
}

void testDG0Prolongation()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = minLevel + 1;

   MeshInfo              mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto bc = BoundaryCondition::createAllInnerBC();

   P0Function< real_t > test_function( "test_function", storage, minLevel, maxLevel, bc );
   P0Function< real_t > solution_function( "solution_function", storage, maxLevel, maxLevel, bc );
   P0Function< real_t > tmp( "tmp", storage, maxLevel, maxLevel, bc );

   auto fun = []( const Point3D& p ) {
      const real_t value1 = p[0] < 0.5 ? 0.5 : 0.2;
      const real_t value2 = p[1] < 0.5 ? -1. : -2.;
      return value1 + value2;
   };

   test_function.interpolate( fun, minLevel );
   solution_function.interpolate( fun, maxLevel );

   DG0Prolongation prolongation;

   prolongation.prolongate( *test_function.getDGFunction(), minLevel, All );

   tmp.assign( { +1., -1. }, { test_function, solution_function }, maxLevel, All );
   const real_t error = tmp.dotGlobal( tmp, maxLevel, All );

   WALBERLA_CHECK_FLOAT_EQUAL( error, 0. );

   // VTKOutput vtk( "../../output/", "DGProlongationTest", storage );
   // vtk.add( test_function );
   // vtk.write( minLevel );
   // vtk.write( maxLevel );
}

void compareProlongationWithRestrictionDG1( const std::string& meshFile )
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = minLevel + 1;

   const uint_t degree = 1;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto dgBasis  = std::make_shared< DGBasisLinearLagrange_Example >();
   auto massForm = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > fine_function( "test_function", storage, minLevel, maxLevel, dgBasis, degree );
   DGFunction< real_t > coarse_function( "test_function", storage, minLevel, maxLevel, dgBasis, degree );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, dgBasis, degree );

   DGFunction< idx_t > numerator( "numerator", storage, minLevel, maxLevel, dgBasis, degree );
   numerator.enumerate( minLevel );
   numerator.enumerate( maxLevel );
   DGOperator M( storage, minLevel, maxLevel, massForm );
   tmp.evaluateLinearFunctional( []( const Point3D& p ) { return 0.5 * p[0] + p[1]; }, minLevel );
   tmp.evaluateLinearFunctional( []( const Point3D& p ) { return 1.2 * p[0] - p[1]; }, maxLevel );
   // interpolate test_function on lower level
   {
      PETScCGSolver< DGOperator > solverM( storage, minLevel, numerator );
      solverM.solve( M, coarse_function, tmp, minLevel );
   }
   // interpolate solution function on upper level
   {
      PETScCGSolver< DGOperator > solverM( storage, maxLevel, numerator );
      solverM.solve( M, fine_function, tmp, maxLevel );
   }

   DG1Restriction  restriction;
   DG1Prolongation prolongation;

   // value_restrict = coarse_function * ( restriction * fine_function )
   restriction.restrict( fine_function, maxLevel, All );
   const real_t value_restrict = coarse_function.dotGlobal( fine_function, minLevel, All );

   // value_prolongate = fine_function * ( prolongation * coarse_function )
   prolongation.prolongate( coarse_function, minLevel, All );
   const real_t value_prolongate = fine_function.dotGlobal( coarse_function, maxLevel, All );

   const real_t error = value_restrict - value_prolongate;

   WALBERLA_CHECK_FLOAT_EQUAL( error, 0. );
}

void compareProlongationWithRestrictionDG0( const std::string& meshFile )
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = minLevel + 1;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto bc = BoundaryCondition::createAllInnerBC();

   P0Function< real_t > fine_function( "test_function", storage, minLevel, maxLevel, bc );
   P0Function< real_t > coarse_function( "test_function", storage, minLevel, maxLevel, bc );

   coarse_function.interpolate( []( const Point3D& p ) { return 0.5 * p[0] + p[1]; }, minLevel );
   fine_function.interpolate( []( const Point3D& p ) { return 1.2 * p[0] - p[1]; }, maxLevel );

   DG0Restriction  restriction;
   DG0Prolongation prolongation;

   // value_restrict = coarse_function * ( restriction * fine_function )
   restriction.restrict( *fine_function.getDGFunction(), maxLevel, All );
   const real_t value_restrict = coarse_function.dotGlobal( fine_function, minLevel, All );

   // value_prolongate = fine_function * ( prolongation * coarse_function )
   prolongation.prolongate( *coarse_function.getDGFunction(), minLevel, All );
   const real_t value_prolongate = fine_function.dotGlobal( coarse_function, maxLevel, All );

   const real_t error = value_restrict - value_prolongate;

   WALBERLA_CHECK_FLOAT_EQUAL( error, 0. );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   testDG1Prolongation( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   testDG1Prolongation( hyteg::prependHyTeGMeshDir( "2D/circle.msh" ) );

   testDG1Prolongation( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ) );
   testDG1Prolongation( hyteg::prependHyTeGMeshDir( "3D/cube_6el.msh" ) );

   testDG0Prolongation();

   testDG1Prolongation( hyteg::prependHyTeGMeshDir( "2D/circle.msh" ) );

   // 2d comparisons
   compareProlongationWithRestrictionDG1( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   compareProlongationWithRestrictionDG1( hyteg::prependHyTeGMeshDir( "2D/circle.msh" ) );

   compareProlongationWithRestrictionDG0( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   compareProlongationWithRestrictionDG0( hyteg::prependHyTeGMeshDir( "2D/circle.msh" ) );

   // 3d comparisons
   compareProlongationWithRestrictionDG1( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ) );
   compareProlongationWithRestrictionDG1( hyteg::prependHyTeGMeshDir( "3D/cube_6el.msh" ) );

   // uncommented, when DG0 interpolation in 3D is implemented
   // compareProlongationWithRestrictionDG0( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ));
   // compareProlongationWithRestrictionDG0( hyteg::prependHyTeGMeshDir( "3D/cube_6el.msh" ));

   return 0;
}
