/*
 * Copyright (c) 2023 Ponsuganth Ilangovan P
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

#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"

#include "mixed_operator/VectorMassOperator.hpp"
// #include "hyteg/p2functionspace/P2SurfaceDeltaOperator.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
// #include "hyteg/elementwiseoperators/P2ToP1ElementwiseBlendingGradientOperator.hpp"
#include "hyteg/eigen/EigenSparseDirectSolver.hpp"
#include "hyteg/operatorgeneration/generated/EvaluateViscosityViscoplastic/P2EvaluateViscosityViscoplastic.hpp"
#include "hyteg/operatorgeneration/generated/FullStokesViscoplastic/P2VectorElementwiseFullStokesViscoplastic.hpp"
#include "hyteg/operatorgeneration/generated/Ones/P2Ones.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

namespace hyteg {
namespace operatorgeneration {

using P2P1StokesViscoplastic =
    detail::P2P1StokesNonlinViscOperatorTemplate< operatorgeneration::P2VectorElementwiseFullStokesViscoplastic,
                                                  operatorgeneration::P1ToP2GradientOperator,
                                                  operatorgeneration::P2ToP1DivergenceOperator >;
} // namespace operatorgeneration
} // namespace hyteg

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

typedef operatorgeneration::P2P1StokesFullP1ViscosityOperator StokesOperator;
typedef P2P1ElementwiseBlendingStokesOperator      StokesOperatorLinear;
typedef operatorgeneration::P2P1StokesViscoplastic StokesOperatorNonlinear;

typedef StrongFreeSlipWrapper< StokesOperator, P2ProjectNormalOperator, true > StokesOperatorFS;

class P2TransportTimesteppingOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2TransportTimesteppingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                    const uint_t                               minLevel,
                                    const uint_t                               maxLevel,
                                    real_t                                     k )
   : Operator( storage, minLevel, maxLevel )
   , diffusionOperator( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   , k_( k )
   {}

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      diffusionOperator.apply( src, dst, level, flag, updateType );
      dst.assign( { k_ * dt }, { dst }, level, flag );
      massOperator.apply( src, dst, level, flag, Add );
   }

   void setDt( real_t dt_ ) { dt = dt_; }

 private:
   P2ElementwiseBlendingLaplaceOperator diffusionOperator;
   P2ElementwiseBlendingMassOperator    massOperator;

   real_t dt = 0.01, k_ = 1.0;
};

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./TosiNonlinearUnitSquare.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   const uint_t nx = mainConf.getParameter< uint_t >( "nx" );
   const uint_t ny = mainConf.getParameter< uint_t >( "ny" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISSCROSS, nx, ny );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   enum BoundaryMarkers
   {
      Bottom = 23,
      Right,
      Left,
      Top,
      Corners
   };

   const real_t threshold = 1e-6;

   std::function< bool( const Point3D& ) > bottomMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[1] ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > rightMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[0] - 1.0 ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > leftMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[0] ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > topMarker = [threshold]( const Point3D& x ) {
      if ( std::abs( x[1] - 1.0 ) < threshold )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > cornerMarker = [&]( const Point3D& x ) {
      if ( ( topMarker( x ) && rightMarker( x ) ) || ( topMarker( x ) && leftMarker( x ) ) ||
           ( bottomMarker( x ) && leftMarker( x ) ) || ( bottomMarker( x ) && rightMarker( x ) ) )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Bottom, bottomMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Left, leftMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Right, rightMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Top, topMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Corners, cornerMarker );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   BoundaryCondition bcTemp, bcVelocity, bcVelocityX, bcVelocityY;

   bcTemp.createAllInnerBC();
   bcTemp.createDirichletBC( "DirichletTopAndBottom",
                             { BoundaryMarkers::Bottom, BoundaryMarkers::Top, BoundaryMarkers::Corners } );
   bcTemp.createNeumannBC( "NeumannSides", { BoundaryMarkers::Left, BoundaryMarkers::Right } );

   bcVelocity.createAllInnerBC();
   bcVelocity.createDirichletBC( "CornersDirichlet", { BoundaryMarkers::Corners } );
   bcVelocity.createFreeslipBC(
       "FreeslipAll", { BoundaryMarkers::Bottom, BoundaryMarkers::Top, BoundaryMarkers::Left, BoundaryMarkers::Right } );

   bcVelocityX.createAllInnerBC();
   bcVelocityX.createDirichletBC( "DirichletLeftAndRight",
                                  { BoundaryMarkers::Left, BoundaryMarkers::Right, BoundaryMarkers::Corners } );
   bcVelocityX.createNeumannBC( "NeumannTopAndBottom", { BoundaryMarkers::Bottom, BoundaryMarkers::Top } );

   bcVelocityY.createAllInnerBC();
   bcVelocityY.createDirichletBC( "DirichletTopAndBottom",
                                  { BoundaryMarkers::Bottom, BoundaryMarkers::Top, BoundaryMarkers::Corners } );
   bcVelocityY.createNeumannBC( "NeumannLeftAndRight", { BoundaryMarkers::Left, BoundaryMarkers::Right } );

   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > uOp( "uOp", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > uPrev( "uPrev", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel, bcVelocity );

   P2P1TaylorHoodFunction< real_t > uDirect( "uDirect", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > fDirect( "fDirect", storage, minLevel, maxLevel );

   uDirect.uvw().component( 0 ).setBoundaryCondition( bcVelocityX );
   uDirect.uvw().component( 1 ).setBoundaryCondition( bcVelocityY );

   fDirect.uvw().component( 0 ).setBoundaryCondition( bcVelocityX );
   fDirect.uvw().component( 1 ).setBoundaryCondition( bcVelocityY );

   P2P1TaylorHoodFunction< real_t > fTemp( "fTemp", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > fStrong( "fStrong", storage, minLevel, maxLevel, bcVelocity );

   P2Function< real_t >  T( "T", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > TPrev( "TPrev", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > fT( "fT", storage, minLevel, maxLevel, bcTemp );

   P1VectorFunction< real_t > uxDer( "uxDer", storage, minLevel, maxLevel, bcVelocity );
   P1VectorFunction< real_t > uyDer( "uyDer", storage, minLevel, maxLevel, bcVelocity );

   real_t diffusivity = mainConf.getParameter< real_t >( "diffusivity" );

   P1Function< real_t > onesP1( "onesP1", storage, minLevel, maxLevel );
   P1Function< real_t > viscosityP1( "viscosityP1", storage, minLevel, maxLevel );
   P1Function< real_t > viscosityP1Nonlinear( "viscosityP1Nonlinear", storage, minLevel, maxLevel );
   P1Function< real_t > viscosityP1Out( "viscosityP1Out", storage, minLevel, maxLevel );

   onesP1.interpolate( 1.0, maxLevel, All );
   viscosityP1.interpolate( 1.0, maxLevel, All );
   viscosityP1Nonlinear.interpolate( 1.0, maxLevel, All );
   viscosityP1Out.interpolate( 1.0, maxLevel, All );

   real_t etaStar = real_c( 0.001 );
   real_t sigmaY  = real_c( 1.0 );

   auto stokesOperator = std::make_shared< StokesOperator >(storage, minLevel, maxLevel, viscosityP1Nonlinear);
   auto stokesOperatorLinear = std::make_shared< StokesOperatorLinear >( storage, minLevel, maxLevel );
   auto stokesOperatorNonlinear = std::make_shared< StokesOperatorNonlinear >( storage,
                                                             minLevel,
                                                             maxLevel,
                                                             viscosityP1,
                                                             uOp.uvw().component( 0u ).getVertexDoFFunction(),
                                                             uOp.uvw().component( 1u ).getVertexDoFFunction(),
                                                             etaStar,
                                                             sigmaY );

   P2TransportTimesteppingOperator         transportOperator( storage, minLevel, maxLevel, diffusivity );
   P2ElementwiseBlendingMassOperator       massOperator( storage, minLevel, maxLevel );
   P2ElementwiseBlendingVectorMassOperator vectorMassOperator( storage, minLevel, maxLevel );

   operatorgeneration::P2Ones onesOp( storage, minLevel, maxLevel );
   onesOp.computeInverseDiagonalOperatorValues();

   operatorgeneration::P2EvaluateViscosityViscoplastic evaluateVisc( storage,
                                                                     minLevel,
                                                                     maxLevel,
                                                                     viscosityP1,
                                                                     uOp.uvw().component( 0u ).getVertexDoFFunction(),
                                                                     uOp.uvw().component( 1u ).getVertexDoFFunction(),
                                                                     etaStar,
                                                                     sigmaY );

   // P2ToP1ElementwiseBlendingGradientOperator gradientOperator(storage, minLevel, maxLevel);

   std::function< void( const Point3D&, Point3D& ) > normalsFS = [&]( const Point3D& x, Point3D& normal ) {
      if ( rightMarker( x ) )
      {
         normal[0] = 1.0;
         normal[1] = 0.0;
      }
      else if ( leftMarker( x ) )
      {
         normal[0] = -1.0;
         normal[1] = 0.0;
      }
      else if ( topMarker( x ) )
      {
         normal[0] = 0.0;
         normal[1] = 1.0;
      }
      else if ( bottomMarker( x ) )
      {
         normal[0] = 0.0;
         normal[1] = -1.0;
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Probably shoudln't be here!" );
      }
   };

   auto projectionNormal = std::make_shared< P2ProjectNormalOperator >( storage, minLevel, maxLevel, normalsFS );

   auto mmocTransport = MMOCTransport< P2Function< real_t > >( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   StokesOperatorFS stokesOperatorFS( stokesOperator, projectionNormal, FreeslipBoundary );

   uint_t stokesIter   = mainConf.getParameter< uint_t >( "stokesIter" );
   real_t stokesRelTol = mainConf.getParameter< real_t >( "stokesRelTol" );

   uint_t transportIter   = mainConf.getParameter< uint_t >( "transportIter" );
   real_t transportRelTol = mainConf.getParameter< real_t >( "transportRelTol" );

   MinResSolver< StokesOperatorFS >            minresSolver( storage, minLevel, maxLevel, stokesIter, stokesRelTol );
   CGSolver< P2TransportTimesteppingOperator > transportSolver( storage, minLevel, maxLevel, transportIter, transportRelTol );

   // EigenSparseDirectSolver<StokesOperator> eigenDirect(storage, maxLevel);
   // eigenDirect.setReassembleMatrix(true);

   PETScLUSolver< StokesOperator > petscDirect( storage, maxLevel );
   PETScLUSolver< StokesOperatorLinear > petscDirectLinear( storage, maxLevel );
   PETScLUSolver< StokesOperatorNonlinear > petscDirectNonlinear( storage, maxLevel );
   petscDirect.setReassembleMatrix( true );

   bool verbose = mainConf.getParameter< bool >( "verbose" );
   minresSolver.setPrintInfo( verbose );

   const real_t Ra     = mainConf.getParameter< real_t >( "RayleighNumber" );
   const real_t deltaT = mainConf.getParameter< real_t >( "deltaT" );
   const real_t deltaZ = mainConf.getParameter< real_t >( "deltaZ" );

   const real_t AiniPerturb = 0.05;

   std::string outputPath     = mainConf.getParameter< std::string >( "outputPath" );
   std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );

   bool useAdios2 = mainConf.getParameter< bool >( "useAdios2" );

   // auto clockTime  = std::time( nullptr );
   // auto clockTimeM = *std::localtime( &clockTime );

   std::ostringstream ossVtkName;
   ossVtkName << outputFilename; // << "_" << std::put_time( &clockTimeM, "%d-%m-%Y_%H-%M-%S" );

   VTKOutput vtkOutput( outputPath, ossVtkName.str(), storage );

#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::string                    adiosXmlConfig = mainConf.getParameter< std::string >( "adiosXmlConfig" );
   std::shared_ptr< AdiosWriter > adios2Output =
       std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );
#endif

   std::function< real_t( const Point3D& ) > TIni = [&]( const Point3D& x ) {
      return ( 1 - x[1] ) + AiniPerturb * std::cos( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] );
   };

   T.interpolate( TIni, maxLevel, All );

   real_t dt = 0.0001;

   real_t cflMax = mainConf.getParameter< real_t >( "cflMax" );

   real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );

   real_t stokesResidualInitialPicard = 0.0, stokesResidualInitial = 0.0, stokesResidualFinal = 0.0,
          stokesResidualBeforeViscRecalc = 0.0;

   std::function< real_t() > calculateResidual = [&]() {
      stokesOperator->apply( u, fTemp, maxLevel, Inner );
      fTemp.assign( { 1.0, -1.0 }, { fTemp, f }, maxLevel, Inner );

      return fTemp.uvw().dotGlobal( fTemp.uvw(), maxLevel, Inner );
   };

   real_t gammaT = std::log( deltaT );
   real_t gammaZ = std::log( deltaZ );

   //  real_t etaStar = real_c(0.001);
   //  real_t sigmaY = real_c(1.0);

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > viscFunc = [&]( const Point3D&               x,
                                                                                           const std::vector< real_t >& vals ) {
      real_t Tval = vals[0];

      real_t etaLin = std::exp( -gammaT * Tval + gammaZ * ( 1 - x[1] ) );

      return etaLin;
   };

   std::function< void() > calculateNewViscosity = [&]() { 
      viscosityP1.interpolate( viscFunc, { T.getVertexDoFFunction() }, maxLevel, All );
      viscosityP1Nonlinear.assign({1.0}, {viscosityP1}, maxLevel, All);
   };

   std::function< void() > calculateNewViscosityNonlinear = [&]() { 
      viscosityP1.interpolate( viscFunc, { T.getVertexDoFFunction() }, maxLevel, All );
      evaluateVisc.apply(onesP1, viscosityP1Nonlinear, maxLevel, All);
      viscosityP1Nonlinear.multElementwise({viscosityP1Nonlinear, *(onesOp.getInverseDiagonalValues())}, maxLevel, All);
   };

   std::function< void() > solveU = [&]() {
      u.interpolate( 0.0, maxLevel, DirichletBoundary );
      fStrong.uvw().component( 1U ).interpolate( Ra, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
      fStrong.uvw().component( 1U ).multElementwise(
          { fStrong.uvw().component( 1U ), T }, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );

      vectorMassOperator.apply( fStrong.uvw(), f.uvw(), maxLevel, All );

      calculateNewViscosity();

      stokesResidualInitial = calculateResidual();

      fDirect.assign( { 1.0 }, { f }, maxLevel, All );
      // eigenDirect.solve(*stokesOperator, uDirect, fDirect, maxLevel);
      petscDirect.solve( *stokesOperator, uDirect, fDirect, maxLevel );
      // petscDirectLinear.solve(*stokesOperatorLinear, uDirect, fDirect, maxLevel);
      u.assign( { 1.0 }, { uDirect }, maxLevel, All );

      // projectionNormal->project( f, maxLevel, FreeslipBoundary );
      // minresSolver.solve( stokesOperatorFS, u, f, maxLevel );

      stokesResidualBeforeViscRecalc = calculateResidual();

      calculateNewViscosity();

      stokesResidualFinal = calculateResidual();
   };

   uint_t nPicard = mainConf.getParameter< uint_t >( "nPicard" );

   real_t reltolPicard = mainConf.getParameter< real_t >( "reltolPicard" );
   real_t abstolPicard = mainConf.getParameter< real_t >( "abstolPicard" );

   std::function< void() > solveUNonlinear = [&]() {
      u.interpolate( 0.0, maxLevel, DirichletBoundary );
      fStrong.uvw().component( 1U ).interpolate( Ra, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
      fStrong.uvw().component( 1U ).multElementwise(
          { fStrong.uvw().component( 1U ), T }, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );

      vectorMassOperator.apply( fStrong.uvw(), f.uvw(), maxLevel, All );

      for ( uint_t iPicard = 0U; iPicard < nPicard; iPicard++ )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Picard, iter %d\n\n", iPicard ) );

         uOp.assign( { 1.0 }, { u }, maxLevel, All );

         calculateNewViscosityNonlinear();

         stokesResidualInitial = calculateResidual();

         if ( iPicard == 0U )
         {
            stokesResidualInitialPicard = stokesResidualInitial;
         }

         fDirect.assign( { 1.0 }, { f }, maxLevel, All );
         // eigenDirect.solve(*stokesOperator, uDirect, fDirect, maxLevel);
         petscDirect.solve( *stokesOperator, uDirect, fDirect, maxLevel );
         u.assign( { 1.0 }, { uDirect }, maxLevel, All );

         uOp.assign( { 1.0 }, { u }, maxLevel, All );

         // projectionNormal->project( f, maxLevel, FreeslipBoundary );
         // minresSolver.solve( stokesOperatorFS, u, f, maxLevel );

         stokesResidualBeforeViscRecalc = calculateResidual();

         calculateNewViscosityNonlinear();

         stokesResidualFinal = calculateResidual();

         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( "\n\nInitial residual = %4.7e\nMid-life residual = %4.7e\nFinal residual = %4.7e",
                               stokesResidualInitial,
                               stokesResidualBeforeViscRecalc,
                               stokesResidualFinal ) );

         if ( stokesResidualFinal / stokesResidualInitialPicard < reltolPicard )
         {
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Picard relative residual reached, exiting!" ) );
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nEnding Picard, iter %d", iPicard ) );
            break;
         }

         if ( stokesResidualFinal < abstolPicard )
         {
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Picard absolute residual reached, exiting!" ) );
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nEnding Picard, iter %d", iPicard ) );
            break;
         }

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nEnding Picard, iter %d", iPicard ) );
      }
   };

   std::function< void() > solveT = [&]() {
      real_t vMax = u.uvw().getMaxComponentMagnitude( maxLevel, All );

      dt = cflMax * hMin / vMax;

      T.interpolate( TIni, maxLevel, DirichletBoundary );
      mmocTransport.step( T, u.uvw(), uPrev.uvw(), maxLevel, Inner | NeumannBoundary | FreeslipBoundary, dt, 1 );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "dt = %4.4e", dt ) );

      transportOperator.setDt( dt );

      T.interpolate( TIni, maxLevel, DirichletBoundary );
      massOperator.apply( T, fT, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
      transportSolver.solve( transportOperator, T, fT, maxLevel );
   };

   //    std::function< void() > solvePicard = [&]() {
   //       solveU();
   //    };

   uint_t NTimesteps = mainConf.getParameter< uint_t >( "NTimesteps" );

   if ( useAdios2 )
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      adios2Output->add( u );
      adios2Output->add( T );
      adios2Output->add( viscosityP1 );
      adios2Output->add( viscosityP1Out );
      adios2Output->add( viscosityP1Nonlinear );
      adios2Output->add( *( onesOp.getInverseDiagonalValues() ) );
#else
      WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
   }
   else
   {
      vtkOutput.add( u );
      vtkOutput.add( T );
      vtkOutput.add( viscosityP1 );
      vtkOutput.add( viscosityP1Out );
      vtkOutput.add( viscosityP1Nonlinear );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Starting initial Stokes solve\n\n" );

   solveU();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Initial Stokes solve done" ) );

   // WALBERLA_LOG_INFO_ON_ROOT(walberla::format( "Initial residual = %4.7e, Final residual = %4.7e", stokesResidualInitial, stokesResidualFinal ));

   uPrev.assign( { 1.0 }, { u }, maxLevel, All );
   uOp.assign( { 1.0 }, { u }, maxLevel, All );

   if ( useAdios2 )
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      adios2Output->write( maxLevel, 0U );
#else
      WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
   }
   else
   {
      vtkOutput.write( maxLevel, 0U );
   }

   real_t simulationTime = 0.0;

   real_t endTime = mainConf.getParameter< real_t >( "endTime" );

   uint_t iTimestep = 1U;

   uint_t logNsFreq = mainConf.getParameter< uint_t >( "logNsFreq" );

   uint_t vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );

   for ( iTimestep = 1U; iTimestep <= NTimesteps; iTimestep++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Started step %d at time %4.7e", iTimestep, simulationTime ) );

      solveT();

      simulationTime += dt;

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Transport done", iTimestep ) );

      uPrev.assign( { 1.0 }, { u }, maxLevel, All );
      uOp.assign( { 1.0 }, { u }, maxLevel, All );

      evaluateVisc.apply(onesP1, viscosityP1Out, maxLevel, All);
      viscosityP1Out.multElementwise({viscosityP1Out, *(onesOp.getInverseDiagonalValues())}, maxLevel, All);

      WALBERLA_LOG_INFO_ON_ROOT( "Starting Stokes\n\n" );

      solveUNonlinear();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Stokes done", iTimestep ) );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Ended step %d", iTimestep ) );

      if ( simulationTime > endTime )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Specified endTime reached!" );
         break;
      }

      if ( iTimestep % vtkWriteFrequency == 0 )
      {
         if ( useAdios2 )
         {
#ifdef HYTEG_BUILD_WITH_ADIOS2
            adios2Output->write( maxLevel, iTimestep );
#else
            WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
         }
         else
         {
            vtkOutput.write( maxLevel, iTimestep );
         }
      }

      if ( iTimestep % logNsFreq == 0 )
      {
         uint_t numDoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Number of global DoFs = %d", numDoFs ) );

         real_t nusseltNumber = nusseltcalc::calculateNusseltNumber2D( T, maxLevel, 0.001, 1e-6, 101 );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "nusseltNumber = %4.7e", nusseltNumber ) );
      }
   }

   uint_t numDoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Number of global DoFs = %d", numDoFs ) );

   real_t nusseltNumber = nusseltcalc::calculateNusseltNumber2D( T, maxLevel, 0.001, 1e-6, 101 );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "nusseltNumber = %4.7e", nusseltNumber ) );

   return 0;
}