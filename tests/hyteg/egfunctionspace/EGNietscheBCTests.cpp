/*
* Copyright (c) 2022 Andreas Wagner.
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

#include <hyteg/dgfunctionspace/P1WithDGFormOperator.hpp>

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/dg1functionspace/DG1Operator.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/egfunctionspace/EGOperatorsNew.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

namespace hyteg {
using walberla::real_t;
using walberla::math::pi;

void applyDirichletBC( const std::shared_ptr< dg::DGForm >& form, const P1Function< real_t >& f, const uint_t level )
{
   if ( f.getStorage()->hasGlobalCells() )
      WALBERLA_ABORT( "3D not implemented yet" );

   using ValueType = real_t;

   const int dim = 2;

   auto basis = std::make_shared< DGBasisLinearLagrange_Example >();

   auto boundaryCondition = f.getBoundaryCondition();
   auto storage           = f.getStorage();

   for ( const auto& faceIt : f.getStorage()->getFaces() )
   {
      const auto faceId = faceIt.first;
      const Face& face   = *faceIt.second;

      // ?
      // if ( !face.hasData( f.getFaceDataID() ) )
      //    continue;

      const auto polyDegree = 1;
      const auto numDofs    = 3;
      const auto dofMemory = face.getData( f.getFaceDataID() )->getPointer( level );

      for ( const auto& [n, _] : face.getIndirectNeighborFaceIDsOverEdges() )
      {
         auto data     = face.getData( f.getFaceGLDataID( n ) );
         auto data_ptr = data->getPointer( level );

         for ( uint_t i = 0; i < data->getSize( level ); i += 1 )
            data_ptr[i] = 0.;
      }

      // zero out halo
      for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
      {
         if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
         {
            auto arrayIdx       = vertexdof::macroface::index( level, idx.x(), idx.y() );
            dofMemory[arrayIdx] = real_c( 0 );
            // WALBERLA_LOG_INFO_ON_ROOT(idx.x() << " " << idx.y());
         }
      }

      for ( auto faceType : facedof::allFaceTypes )
      {
         for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
         {
            volumedofspace::indexing::ElementNeighborInfo neighborInfo(
                elementIdx, faceType, level, boundaryCondition, faceId, storage );

            auto vertexDoFIndicesArray = facedof::macroface::getMicroVerticesFromMicroFace( elementIdx, faceType );

            for ( uint_t n = 0; n < 3; n++ )
            {
               if ( neighborInfo.atMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == DirichletBoundary )
               {
                  Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > localMat;
                  localMat.resize( Eigen::Index( numDofs ), 1 );
                  localMat.setZero();
                  form->integrateRHSDirichletBoundary( dim,
                                                       neighborInfo.elementVertexCoords(),
                                                       neighborInfo.interfaceVertexCoords( n ),
                                                       neighborInfo.oppositeVertexCoords( n ),
                                                       neighborInfo.outwardNormal( n ),
                                                       *basis,
                                                       polyDegree,
                                                       localMat );

                  for ( uint_t dofIdx = 0; dofIdx < numDofs; dofIdx++ )
                  {
                     const uint_t p1Index = vertexdof::macroface::index(
                         level, vertexDoFIndicesArray[dofIdx].x(), vertexDoFIndicesArray[dofIdx].y() );
                     dofMemory[p1Index] += ValueType( localMat( Eigen::Index( dofIdx ), 0 ) );
                  }
               }
            }
         }
      }
   }

   if ( dim == 2 )
   {
      f.template communicateAdditively< Face, Edge >( level, All ^ DirichletBoundary, *storage,  false );
      f.template communicateAdditively< Face, Vertex >( level, All ^ DirichletBoundary, *storage,  false );

      // f.template communicateAdditively< Face, Edge >( level, All, *storage, false );
      // f.template communicateAdditively< Face, Vertex >( level, All, *storage, false );
   }
   else
   {
      f.template communicateAdditively< Cell, Face >( level, All ^ DirichletBoundary, *storage, false );
      f.template communicateAdditively< Cell, Edge >( level, All ^ DirichletBoundary, *storage, false );
      f.template communicateAdditively< Cell, Vertex >( level, All ^ DirichletBoundary, *storage, false );
   }
}

void applyDirichletBC( const std::shared_ptr< dg::DGForm >& p1Form1,
                       const std::shared_ptr< dg::DGForm >& p1Form2,
                       const std::shared_ptr< dg::DGForm >& dgForm,
                       const EGFunction< real_t >&          f,
                       const uint_t                         level )
{
   applyDirichletBC( p1Form1, f.getConformingPart()->component( 0 ), level );
   applyDirichletBC( p1Form2, f.getConformingPart()->component( 1 ), level );

   // correct ?
   // does not work, wrong basis :)
   f.getDiscontinuousPart()->getDGFunction()->applyDirichletBoundaryConditions( dgForm, level );
}

/// Returns the scaled L2 error.
real_t testP1( uint_t                                    level,
               MeshInfo                                  meshInfo,
               std::function< real_t( const Point3D& ) > solFunc,
               std::function< real_t( const Point3D& ) > rhsFunc,
               bool                                      writeVTK = true )
{
   using namespace dg;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   real_t beta_0 = storage->hasGlobalCells() ? 0.5 : 1.0;

   auto basis                                                = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm0                                         = std::make_shared< EGVectorLaplaceFormP1P1_new_00 >();
   auto laplaceForm1                                         = std::make_shared< EGVectorLaplaceFormP1P1_new_11 >();
   auto laplaceFormDG                                        = std::make_shared< EGVectorLaplaceFormEDGEDG_new >();
   laplaceForm1->callback_Scalar_Variable_Coefficient_2D_g1  = solFunc;
   laplaceForm0->callback_Scalar_Variable_Coefficient_2D_g0  = solFunc;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_2D_g0 = solFunc;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_2D_g1 = solFunc;
   auto massForm                                             = std::make_shared< DGMassForm_Example >();

   EGFunction< real_t > u( "u", storage, level, level );
   EGFunction< real_t > f( "f", storage, level, level );
   EGFunction< real_t > ff( "ff", storage, level, level );
   EGFunction< real_t > sol( "sol", storage, level, level );
   EGFunction< real_t > tmp( "tmp", storage, level, level );
   EGFunction< real_t > err( "err", storage, level, level );
   EGFunction< real_t > Merr( "Merr", storage, level, level );
   EGFunction< real_t > rhs_int( "rhs_int", storage, level, level );

   eg::EGLaplaceOperatorNew A( storage, level, level );
   eg::EGMassOperator       M( storage, level, level );

   DGToP1Operator< P1ToDG1InterpolationForm, real_t > opDGToP1Real(
       storage, level, level, std::make_shared< P1ToDG1InterpolationForm >() );

   // Assemble RHS.
   DG1Function< real_t > fDG( "fDG", storage, level, level );
   fDG.evaluateLinearFunctional( rhsFunc, level );
   fDG.applyDirichletBoundaryConditions( laplaceForm0, level ); // TODO
   f.interpolate( 0, level, All );
   // opDGToP1Real.apply(*fDG.getDGFunction(), f.getConformingPart()->component(0), level, All, Replace);
   // opDGToP1Real.apply(*fDG.getDGFunction(), f.getConformingPart()->component(1), level, All, Replace);
   opDGToP1Real.apply( *fDG.getDGFunction(), ff.getConformingPart()->component( 0 ), level, All, Replace );
   opDGToP1Real.apply( *fDG.getDGFunction(), ff.getConformingPart()->component( 1 ), level, All, Replace );
   f.evaluateLinearFunctional( rhsFunc, rhsFunc, level );
   {
      VTKOutput out( "../../output", "test1", storage );
      out.add( f );
      out.write( level );
   }
   applyDirichletBC( laplaceForm0, laplaceForm1, laplaceFormDG, f, level );
   {
      VTKOutput out( "../../output", "test2", storage );
      out.add( f );
      out.write( level );
   }

   {
      writeDomainPartitioningVTK(*storage, "../../output", "partitioning");
   }

   // rhs_int.interpolate(0, level, All);
   // rhs_int.getConformingPart()->interpolate({rhsFunc, rhsFunc}, level, All);
   // M.apply(rhs_int, f, level, All, Replace);


   // Interpolate solution
   sol.interpolate( { solFunc, solFunc }, level, All );

   // Solve system.
   PETScCGSolver< eg::EGLaplaceOperatorNew > solverA( storage, level, 1e-6, 1e-6, 10000 );
   solverA.disableApplicationBC( true );
   solverA.solve( A, u, f, level );
   // u.getDiscontinuousPart()->interpolate(0, level, All); // TODO

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "DGPoisson2DConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.add( *f.getConformingPart() );
      vtk.add( *ff.getConformingPart() );
      vtk.add( ff );
      vtk.add( ff );
      vtk.add( fDG );
      vtk.write( level );
   }

   return discrL2;
}

void runTest( uint_t                                           minLevel,
              uint_t                                           maxLevel,
              const MeshInfo&                                  meshInfo,
              const std::function< real_t( const Point3D& ) >& solFunc,
              const std::function< real_t( const Point3D& ) >& rhsFunc )
{
   auto l2ConvRate  = std::pow( 2, -( int( 1 ) + 1 ) );
   auto convRateEps = l2ConvRate * 0.1;
   auto err         = hyteg::testP1( minLevel, meshInfo, solFunc, rhsFunc );
   WALBERLA_LOG_INFO_ON_ROOT( " expected L2 rate: " << l2ConvRate << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
   for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
   {
      auto errFiner     = hyteg::testP1( l, meshInfo, solFunc, rhsFunc );
      auto computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

      WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
      err = errFiner;
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
   using walberla::math::pi;

   /*
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, inhom. BC, rhs = 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
        return 2 * x[0]  - x[1];
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
        return 0.;
      };

      hyteg::runTest( 3, 5, meshInfo, solFunc, rhsFunc );
   }
    */

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
      };

      hyteg::runTest( 3, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 8 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) );
      };

      hyteg::runTest( 3, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, inhom. BC, rhs = 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

      hyteg::runTest( 2, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, inhom. BC, rhs != 0 ###" );

      MeshInfo meshInfo =
          hyteg::MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), hyteg::MeshInfo::CRISS, 2, 2 );

      std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
         return std::exp( -x[0] - ( x[1] * x[1] ) );
      };
      std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
         return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
      };

      hyteg::runTest( 3, 4, meshInfo, solFunc, rhsFunc );
   }

   return EXIT_SUCCESS;
}
