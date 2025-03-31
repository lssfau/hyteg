/*
 * Copyright (c) 2017-2025 Marcus Mohr.
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
#include "hyteg/experimental/P2PlusBubbleFunction.hpp"

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/eigen/EigenSparseDirectSolver.hpp"
#include "hyteg/experimental/P2PlusBubbleFunction.hpp"
#include "hyteg/experimental/P2PlusBubbleOperators/P2PlusBubbleElementwiseDiffusion_float64.hpp"
#include "hyteg/experimental/P2PlusBubbleOperators/P2PlusBubbleElementwiseMass_AffineMap2D_float64.hpp"
#include "hyteg/experimental/P2PlusBubbleOperators/P2PlusBubbleElementwiseMass_AnnulusMap_float64.hpp"
#include "hyteg/experimental/P2PlusBubbleOperators/P2PlusBubbleElementwiseMass_float64.hpp"
#include "hyteg/experimental/P2PlusBubbleVectorFunction.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg_operators/operators/diffusion/P2ElementwiseDiffusion.hpp"

namespace hyteg {

enum class MeshType
{
   SIMPLE_MESH,
   ANNULUS,
   UNIT_SQUARE,
   OTHER_UNIT_SQUARE_AFFINE_BLENDING,
   QUADRATIC_MAPPED_SQUARE
};

// Define the matrix for the OTHER_UNIT_SQUARE_AFFINE_BLENDING meshing type
Matrix2r getMatrix()
{
   Matrix2r mat;
   real_t   phi = real_c( 2.0 / 9.0 ) * walberla::math::pi;
   mat( 0, 0 )  = +std::cos( phi );
   mat( 0, 1 )  = -std::sin( phi );
   mat( 1, 0 )  = +std::sin( phi ) * real_c( 2.25 );
   mat( 1, 1 )  = +std::cos( phi ) * real_c( 2.25 );

   return mat;
}

std::shared_ptr< PrimitiveStorage > generateStorage( MeshType meshType )
{
   std::shared_ptr< PrimitiveStorage > storage;
   uint_t                              additionalHaloDepth{ 0 };

   switch ( meshType )
   {
   case MeshType::SIMPLE_MESH: {
      // simple mesh
      MeshInfo              mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      storage = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
      break;
   }

   case MeshType::ANNULUS: {
      // annulus
      MeshInfo              mesh = MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), MeshInfo::CROSS, 4, 2 );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      AnnulusMap::setMap( setupStorage );
      storage = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
      break;
   }

   case MeshType::UNIT_SQUARE: {
      // unit square
      MeshInfo mesh = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      storage = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
      break;
   }

   case MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING: {
      // unit square with affine blending
      MeshInfo mesh = MeshInfo::meshRectangle(
          Point2D( real_c( -1 ), real_c( -1 ) ), Point2D( real_c( +1 ), real_c( +1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // define our affine map (rotation + scaling + shift)
      Matrix2r mat = getMatrix();
      Point2D  vec( real_c( -7.0 ), real_c( 3.0 ) );
      AffineMap2D::setMap( setupStorage, mat, vec );
      storage = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
      break;
   }

   case MeshType::QUADRATIC_MAPPED_SQUARE: {
      // quadratic mapping
      MeshInfo mesh = MeshInfo::meshRectangle(
          Point2D( real_c( -1 ), real_c( -1 ) ), Point2D( real_c( +1 ), real_c( +1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::vector< std::function< real_t( const hyteg::Point3D& ) > > mapper = {
          []( const Point3D& x ) { return x[0]; },
          []( const Point3D& x ) {
             return x[1] + ( real_c( 1 ) - x[0] * x[0] ) * x[1] * ( std::sqrt( real_c( 2 ) ) - real_c( 1 ) );
          } };
      storage                     = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
      uint_t     level            = 4;
      const auto microMeshDegree2 = std::make_shared< micromesh::MicroMesh >( storage, level, level, 2, 2 );
      micromesh::interpolateAndCommunicate( *microMeshDegree2, mapper, level );
      storage->setMicroMesh( microMeshDegree2 );
      break;
   }
   }

   return storage;
}

void testInterpolation( std::shared_ptr< PrimitiveStorage > storage, std::array< real_t, 3 > tolerance, bool doVTKOutput = false )
{
   const uint_t minLevel = 4;
   const uint_t maxLevel = minLevel;

   auto function0 = P2PlusBubbleFunction< real_t >( "poly0", storage, minLevel, maxLevel );
   auto function1 = P2PlusBubbleFunction< real_t >( "poly1", storage, minLevel, maxLevel );
   auto function2 = P2PlusBubbleFunction< real_t >( "poly2", storage, minLevel, maxLevel );
   auto function3 = P2PlusBubbleFunction< real_t >( "poly3", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > polyDegree0 = []( const Point3D& ) -> real_t { return real_c( 1 ); };

   std::function< real_t( const hyteg::Point3D& ) > polyDegree1 = []( const Point3D& x ) -> real_t {
      return real_c( 2 ) * x[1] - x[0];
   };

   std::function< real_t( const hyteg::Point3D& ) > polyDegree2 = []( const Point3D& x ) -> real_t {
      return x[0] * x[0] + real_c( 2 ) * x[0] * x[1];
   };

   std::function< real_t( const hyteg::Point3D& ) > polyDegree3 = []( const Point3D& x ) -> real_t {
      return x[0] * x[1] * ( real_c( 1 ) - x[0] - x[1] ) * real_c( 27 );
   };

   function0.interpolate( polyDegree0, maxLevel, All );
   function1.interpolate( polyDegree1, maxLevel, All );
   function2.interpolate( polyDegree2, maxLevel, All );
   function3.interpolate( polyDegree3, maxLevel, All );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "p2_plus_bubble", storage );
      vtkOutput.add( function0 );
      vtkOutput.add( function1 );
      vtkOutput.add( function2 );
      vtkOutput.add( function3 );
      vtkOutput.write( maxLevel );
   }

   // for the quadratic polynomial all bubble dofs need to be zero
   const auto& bubbleFunc0 = function0.getVolumeDoFFunction();
   const auto& bubbleFunc1 = function1.getVolumeDoFFunction();
   const auto& bubbleFunc2 = function2.getVolumeDoFFunction();
   const auto& bubbleFunc3 = function3.getVolumeDoFFunction();

   real_t checkValue0 = bubbleFunc0.getMaxDoFMagnitude( maxLevel );
   real_t checkValue1 = bubbleFunc1.getMaxDoFMagnitude( maxLevel );
   real_t checkValue2 = bubbleFunc2.getMaxDoFMagnitude( maxLevel );
   real_t checkValue3 = bubbleFunc3.getMaxDoFMagnitude( maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 0: max magnitude of bubble dofs = " << checkValue0 )
   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 1: max magnitude of bubble dofs = " << checkValue1 )
   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 2: max magnitude of bubble dofs = " << checkValue2 )
   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 3: max magnitude of bubble dofs = " << checkValue3 )

   WALBERLA_CHECK_FLOAT_EQUAL( checkValue0, real_c( 0 ) );
   WALBERLA_CHECK_LESS( checkValue1, tolerance[0] );
   WALBERLA_CHECK_LESS( checkValue2, tolerance[1] );
   WALBERLA_CHECK_LESS( checkValue3, tolerance[2] );
}

void testInterpolationForBCs( bool doVTKOutput = false )
{
   std::shared_ptr< PrimitiveStorage > storage = generateStorage( MeshType::UNIT_SQUARE );

   const uint_t minLevel = 4;
   const uint_t maxLevel = minLevel;

   auto function = P2PlusBubbleFunction< real_t >( "p2+", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > polyDegree2 = []( const Point3D& x ) -> real_t {
      return x[0] * x[0] + real_c( 2 ) * x[0] * x[1];
   };

   function.interpolate( polyDegree2, maxLevel, DirichletBoundary );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "p2_plus_bubble_boundary_interpolation", storage );
      vtkOutput.add( function );
      vtkOutput.write( maxLevel );
   }

   // check that interior values are all zero
   real_t checkVertexDoFs = function.getVertexDoFFunction().getMaxDoFMagnitude( maxLevel, Inner );
   real_t checkEdgeDoFs   = function.getEdgeDoFFunction().getMaxDoFMagnitude( maxLevel, Inner );
   real_t checkBubbleDoFs = function.getVolumeDoFFunction().getMaxDoFMagnitude( maxLevel, Inner );

   WALBERLA_LOG_INFO_ON_ROOT( " * max magnitude of vertex DoFs in the interior = " << checkVertexDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( " * max magnitude of edge   DoFs in the interior = " << checkEdgeDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( " * max magnitude of bubble DoFs in the interior = " << checkBubbleDoFs );

   WALBERLA_CHECK_FLOAT_EQUAL( checkVertexDoFs, real_c( 0 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( checkEdgeDoFs, real_c( 0 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( checkBubbleDoFs, real_c( 0 ) );

   // check values on the boundary values
   checkVertexDoFs = function.getVertexDoFFunction().getMaxDoFMagnitude( maxLevel, DirichletBoundary );
   checkEdgeDoFs   = function.getEdgeDoFFunction().getMaxDoFMagnitude( maxLevel, DirichletBoundary );

   WALBERLA_LOG_INFO_ON_ROOT( " * max magnitude of vertex DoFs on boundary = " << checkVertexDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( " * max magnitude of edge   DoFs on boundary = " << std::scientific << checkEdgeDoFs );

   WALBERLA_CHECK_FLOAT_EQUAL( checkVertexDoFs, real_c( 3 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( checkEdgeDoFs, real_c( 2.937500e+00 ) );
}

void runInterpolationTests()
{
   WALBERLA_LOG_INFO_ON_ROOT( "==========================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " P2PlusBubbleFunction: Interpolation Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "==========================================" );

   std::shared_ptr< PrimitiveStorage > storage1 = generateStorage( MeshType::SIMPLE_MESH );
   std::shared_ptr< PrimitiveStorage > storage2 = generateStorage( MeshType::ANNULUS );
   std::shared_ptr< PrimitiveStorage > storage3 = generateStorage( MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING );
   std::shared_ptr< PrimitiveStorage > storage4 = generateStorage( MeshType::QUADRATIC_MAPPED_SQUARE );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on unmodified mesh" );
   testInterpolation( storage1, { real_c( 3e-16 ), real_c( 3e-16 ), real_c( 2.5e-4 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on mesh with affine blending" );
   testInterpolation( storage3, { real_c( 8e-15 ), real_c( 3e-14 ), real_c( 2.1e-3 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation with quadratic micromesh" );
   testInterpolation( storage4, { real_c( 9e-16 ), real_c( 4.5e-5 ), real_c( 3e-3 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on mesh with annulus blending" );
   testInterpolation( storage2, { real_c( 1.8e-4 ), real_c( 5e-4 ), real_c( 0.03 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on boundary of an unmodified mesh" );
   testInterpolationForBCs();
}

template < template < typename > typename func_t, typename operator_t >
real_t solveDiffusionProblem( bool doVTKOutput = false )
{
   std::string feFamily = FunctionTrait< func_t< real_t > >::getTypeName();
   WALBERLA_LOG_INFO_ON_ROOT( " - solving diffusion problem with " << feFamily );

   // setup domain as unit square (0,1)^2
   MeshInfo mesh = MeshInfo::meshRectangle(
       Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( +1 ), real_c( +1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // pick an analytical solution for the homogeneous diffusion equation
   real_t freqNum{ real_c( 1 ) };

   std::function< real_t( const hyteg::Point3D& ) > expr_analytic = [&freqNum]( const Point3D& x ) {
      using walberla::math::pi;
      return std::sinh( freqNum * pi * x[1] ) / std::sinh( freqNum * pi ) * std::sin( freqNum * pi * x[0] );
      // return real_c(1.0);
      // return real_c(3)*x[0];
   };

   // setup some P2+ functions
   uint_t minLevel = 0;
   uint_t maxLevel = 5;

   func_t< real_t > u_analytic( "true solution", storage, minLevel, maxLevel );
   u_analytic.interpolate( expr_analytic, maxLevel, All );

   func_t< real_t > u_discrete( "discrete solution", storage, minLevel, maxLevel );
   u_discrete.interpolate( expr_analytic, maxLevel, DirichletBoundary );

   func_t< real_t > weak_diffusion( "weak diffusion", storage, minLevel, maxLevel );

   // analytical solution is an eigenfunction for zero eigenvalue, so rhs = 0
   func_t< real_t > rhs( "rhs", storage, minLevel, maxLevel );

   operator_t laplacian( storage, maxLevel, maxLevel );
   laplacian.apply( u_analytic, weak_diffusion, maxLevel, All, Replace );

   // we want to solve this small 2D problem directly
   EigenSparseDirectSolver< operator_t > eigenLU( storage, maxLevel );
   eigenLU.solve( laplacian, u_discrete, rhs, maxLevel );

   // check error (only regression testing, so use solve level for comparison)
   uint_t           errLevel = maxLevel;
   func_t< real_t > error( "error", storage, errLevel, errLevel );
   error.assign( { real_c( +1 ), real_c( -1 ) }, { u_analytic, u_discrete }, errLevel, All );

   real_t errorNorm = std::sqrt( error.dotGlobal( error, errLevel, All ) );
   errorNorm /= real_c( error.getNumberOfGlobalDoFs( errLevel ) );

   // output stuff for visualisation
   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", feFamily + "_diffusion", storage );
      vtkOutput.add( u_analytic );
      vtkOutput.add( u_discrete );
      vtkOutput.add( weak_diffusion );
      vtkOutput.add( rhs );
      vtkOutput.add( error );
      vtkOutput.write( maxLevel );
   }

   WALBERLA_LOG_INFO_ON_ROOT( " - error is approximately = " << errorNorm );

   return errorNorm;
}

void runSingleMassTest( MeshType meshType )
{
   std::shared_ptr< PrimitiveStorage > storage;
   real_t                              volumeCtrl{ real_c( 0 ) };

   if ( meshType == MeshType::UNIT_SQUARE )
   {
      storage    = generateStorage( meshType );
      volumeCtrl = real_c( 1 );
   }
   else if ( meshType == MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING )
   {
      storage         = generateStorage( meshType );
      volumeCtrl      = real_c( 4 );
      Matrix2r mat    = getMatrix();
      real_t   jacDet = std::abs( mat.determinant() );
      volumeCtrl *= jacDet;
   }
   else if ( meshType == MeshType::ANNULUS )
   {
      storage    = generateStorage( meshType );
      volumeCtrl = real_c( 3 ) * walberla::math::pi;
   }
   else
   {
      WALBERLA_ABORT( "Request for unsupported MeshType!" );
   }

   uint_t level = 3;

   P2PlusBubbleFunction< real_t > aux( "aux", storage, level, level );
   P2PlusBubbleFunction< real_t > vecOfOnes( "vecOfOnes", storage, level, level );
   P2PlusBubbleFunction< real_t > integrand( "integrand", storage, level, level );

   // test 1: volume of the domain
   vecOfOnes.interpolate( real_c( 1 ), level, All );

   if ( meshType == MeshType::UNIT_SQUARE )
   {
      using MassOperator = operatorgeneration::P2PlusBubbleElementwiseMass_float64;
      MassOperator massOp( storage, level, level );
      massOp.apply( vecOfOnes, aux, level, All );
   }
   else if ( meshType == MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING )
   {
      using MassOperator = operatorgeneration::P2PlusBubbleElementwiseMass_AffineMap2D_float64;
      MassOperator massOp( storage, level, level );
      massOp.apply( vecOfOnes, aux, level, All );
   }
   else if ( meshType == MeshType::ANNULUS )
   {
      using MassOperator = operatorgeneration::P2PlusBubbleElementwiseMass_AnnulusMap_float64;
      MassOperator massOp( storage, level, level );
      massOp.apply( vecOfOnes, aux, level, All );
   }

   real_t measure = vecOfOnes.dotGlobal( aux, level );

   WALBERLA_LOG_INFO_ON_ROOT( " - computed volume of physical domain = " << std::scientific << measure );
   WALBERLA_LOG_INFO_ON_ROOT( " - analytic volume of physical domain = " << std::scientific << volumeCtrl );

   if ( meshType == MeshType::UNIT_SQUARE || meshType == MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( measure, volumeCtrl );
   }
   else if ( meshType == MeshType::ANNULUS )
   {
      WALBERLA_CHECK_LESS_EQUAL( std::abs( measure - volumeCtrl ), real_c( 1e-7 ) );
      return;
   }

   // test 2: "integrate" polynomial
   std::function< real_t( const hyteg::Point3D& ) > cubicPoly = []( const hyteg::Point3D& coords ) {
      real_t x = coords( 0 );
      real_t y = coords( 1 );
      return x * x * y + real_c( 2 ) * x * x + real_c( 3 ) * y * y + real_c( 4 ) * x + real_c( 5 ) * x * y - y + real_c( 2 );
   };
   vecOfOnes.interpolate( real_c( 1 ), level, All );
   integrand.interpolate( cubicPoly, level, All );

   if ( meshType == MeshType::UNIT_SQUARE )
   {
      using MassOperator = operatorgeneration::P2PlusBubbleElementwiseMass_float64;
      MassOperator massOp( storage, level, level );
      massOp.apply( integrand, aux, level, All );
   }
   else if ( meshType == MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING )
   {
      using MassOperator = operatorgeneration::P2PlusBubbleElementwiseMass_AffineMap2D_float64;
      MassOperator massOp( storage, level, level );
      massOp.apply( integrand, aux, level, All );
   }

   measure = vecOfOnes.dotGlobal( aux, level );

   const auto& bubbleFunc = integrand.getVolumeDoFFunction();
   real_t      controlMax = bubbleFunc.getMaxDoFValue( level );
   real_t      controlMin = bubbleFunc.getMinDoFValue( level );

   WALBERLA_LOG_INFO_ON_ROOT( " - integral over cubic polynomial = " << measure );
   WALBERLA_LOG_INFO_ON_ROOT( " - control bubble dof (max/min) = (" << controlMax << ", " << controlMin << ")" );

   if ( meshType == MeshType::UNIT_SQUARE )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " - check value = " << real_c( 79.0 / 12.0 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, real_c( 79.0 / 12.0 ) );
   }
   else if ( meshType == MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING )
   {
      Matrix2r mat  = getMatrix();
      real_t   ctrl = real_c( 1447289930 ) / real_c( 2499997 ) * std::abs( mat.determinant() );
      WALBERLA_LOG_INFO_ON_ROOT( " - check value = " << ctrl );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, ctrl );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

void runMassTests()
{
   WALBERLA_LOG_INFO_ON_ROOT( "=========================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " P2PlusBubbleFunction: Mass Problem Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "=========================================" );

   WALBERLA_LOG_INFO_ON_ROOT( " Executing test with unit square" );
   runSingleMassTest( MeshType::UNIT_SQUARE );

   WALBERLA_LOG_INFO_ON_ROOT( " Executing test with unit square + affine blending" );
   runSingleMassTest( MeshType::OTHER_UNIT_SQUARE_AFFINE_BLENDING );

   WALBERLA_LOG_INFO_ON_ROOT( " Executing test with annulus" );
   runSingleMassTest( MeshType::ANNULUS );
}

void runDiffusionTest( bool doVTKOutput = true )
{
   WALBERLA_LOG_INFO_ON_ROOT( "==============================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " P2PlusBubbleFunction: Diffusion Problem Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "==============================================" );

   real_t value1 =
       solveDiffusionProblem< P2PlusBubbleFunction, operatorgeneration::P2PlusBubbleElementwiseDiffusion_float64 >( doVTKOutput );
   real_t value2 = solveDiffusionProblem< P2Function, operatorgeneration::P2ElementwiseDiffusion >( doVTKOutput );

   WALBERLA_CHECK_LESS( value1, 2e-9 );
   WALBERLA_CHECK_LESS( value2, 0.5e-9 );
}

void runVectorFunctionTest( bool doVTKOutput = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "========================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " P2PlusBubbleVectorFunction: Basic Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "========================================" );

   WALBERLA_LOG_INFO_ON_ROOT( "* function object creation" );
   std::shared_ptr< PrimitiveStorage > storage = generateStorage( MeshType::SIMPLE_MESH );

   const uint_t level = 3;

   auto vecFunc = P2PlusBubbleVectorFunction< real_t >( "P2+ Vector-Function", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exprX = []( const Point3D& coords ) {
      real_t x = coords( 0 );
      real_t y = coords( 1 );
      return x * x * y + real_c( 2 ) * x * x + real_c( 3 ) * y * y + real_c( 4 ) * x + real_c( 5 ) * x * y - y + real_c( 2 );
   };

   std::function< real_t( const hyteg::Point3D& ) > exprY = []( const Point3D& x ) { return real_c( 1.5 ) * x[0] + x[1]; };

   WALBERLA_LOG_INFO_ON_ROOT( "* function interpolation" );
   vecFunc.interpolate( { exprX, exprY }, level, All );

   if ( doVTKOutput )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "* VTK export" );
      VTKOutput vtkOutput( ".", "p2_plus_bubble_vector", storage );
      vtkOutput.add( vecFunc );
      vtkOutput.write( level );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   using namespace hyteg;

   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   runInterpolationTests();
   runMassTests();
   runDiffusionTest();
   runVectorFunctionTest();

   return EXIT_SUCCESS;
}
