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
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

/// Adaptively refines a domain and interpolates and evaluates a function in DG space.
///
/// The storage (i.e. all allocated memory) is freed after each run. So this is a static, _not_ a dynamic refinement test.
///
/// \param minLevel min micro-refinement level
/// \param maxLevel max micro-refinement level
/// \param degree   element poly degree
void test( uint_t                                           dim,
           uint_t                                           level,
           uint_t                                           coarseRefinements,
           uint_t                                           degree,
           const std::function< real_t( const Point3D& ) >& f,
           real_t                                           maxPointwiseError )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running interpolate + evaluate test with AMR." );
   WALBERLA_LOG_INFO_ON_ROOT( " + dim:                     " << dim );
   WALBERLA_LOG_INFO_ON_ROOT( " + level (micro):           " << level );
   WALBERLA_LOG_INFO_ON_ROOT( " + num refinements (macro): " << coarseRefinements );
   WALBERLA_LOG_INFO_ON_ROOT( " + degree:                  " << degree );

   using namespace dg;

   const bool writeVTK = true;

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/quad_16el.msh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Refine near (0, 0).
   for ( uint_t r = 0; r < coarseRefinements; r++ )
   {
      std::vector< PrimitiveID > faceIDsToRefine;

      for ( const auto& faceID : storage->getFaceIDs() )
      {
         auto coordinates = storage->getFace( faceID )->getCoordinates();
         auto refine      = std::all_of( coordinates.begin(), coordinates.end(), [coarseRefinements, r]( const Point3D& x ) {
            return x.norm() < ( 0.6 / real_c( coarseRefinements ) ) * real_c( coarseRefinements - r );
         } );

         if ( refine )
         {
            faceIDsToRefine.push_back( faceID );
         }
      }

      storage->refinementAndCoarseningHanging( faceIDsToRefine, {} );
   }

   writeDomainPartitioningVTK(
       storage, "../../output", "DGAdaptiveRefinementInterpolateEvaluateTest_Domain_" + std::to_string( coarseRefinements ) );

   walberla::math::seedRandomGenerator( 12345678 );
   const uint_t numRandomEvaluations = 200;

   auto basis    = std::make_shared< DGBasisLinearLagrange_Example >();
   auto massForm = std::make_shared< DGMassForm_Example >();

   P1Function< real_t > interpolatedP1( "u_P1", storage, level, level );
   DGFunction< real_t > u( "u", storage, level, level, basis, degree );
   DGFunction< real_t > tmp( "tmp", storage, level, level, basis, degree );

   DGFunction< idx_t > numerator( "numerator", storage, level, level, basis, degree );
   numerator.enumerate( level );

   DGOperator M( storage, level, level, massForm );

   // Interpolate solution into u
   tmp.evaluateLinearFunctional( f, level );
   PETScCGSolver< DGOperator > solverM( storage, level, numerator );
   solverM.solve( M, u, tmp, level );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "DGAdaptiveRefinementInterpolateEvaluateTest", storage );
      vtk.add( u );
      vtk.add( tmp );
      vtk.write( level );
   }

   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      Point3D coordinates( Point3D::Zero() );
      coordinates[0] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
      coordinates[1] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
      coordinates[2] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

      real_t     value;
      const bool success = u.evaluate( coordinates, level, value, 1e-14 );
      WALBERLA_CHECK( success, "Could not evaluate successfully." );

      const real_t err = std::abs( value - f( coordinates ) );

      WALBERLA_CHECK_LESS(
          err, maxPointwiseError, "Failed at " << coordinates << ", evaluated value: " << value << ", error: " << err );
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   using hyteg::MeshInfo;
   using hyteg::Point2D;
   using hyteg::Point3D;
   using walberla::real_t;

   const uint_t coarseRefinements = 4;

   for ( uint_t r = 0; r < coarseRefinements; r++ )
   {
      hyteg::test(
          2, 3, r, 1, []( const hyteg::Point3D& ) { return 1; }, 1e-9 );
      hyteg::test(
          2, 3, r, 1, []( const hyteg::Point3D& x ) { return x[0] - 2 * x[1]; }, 2e-9 );
   }

#if 0
   hyteg::test(
       3, 3, 1, []( const hyteg::Point3D& ) { return 1; }, 1e-11 );
   hyteg::test(
       3, 3, 1, []( const hyteg::Point3D& x ) { return x[0] - 2 * x[1] + 3 * x[2]; }, 1e-11 );
#endif

   return EXIT_SUCCESS;
}
