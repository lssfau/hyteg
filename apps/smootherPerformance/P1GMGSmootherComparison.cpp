/*
 * Copyright (c) 2020 Andreas Wagner
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
#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/numerictools/SpectrumEstimation.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/SORSmoother.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

namespace hyteg {

class P1JacobiSmoother : public Solver< P1ConstantLaplaceOperator >
{
 public:
   using FunctionType = typename P1ConstantLaplaceOperator::srcType;

   P1JacobiSmoother( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : flag_( hyteg::Inner | hyteg::NeumannBoundary )
   , tmp_( "jac_tmp", storage, minLevel, maxLevel )
   , weight_( 1. )
   {}

   void setWeight( const real_t& weight ) { weight_ = weight; }

   void solve( const P1ConstantLaplaceOperator&                   A,
               const typename P1ConstantLaplaceOperator::srcType& x,
               const typename P1ConstantLaplaceOperator::dstType& b,
               const walberla::uint_t                                level ) override
   {
      tmp_.assign( {1.0}, {x}, level, hyteg::All );
      A.smooth_jac( x, b, tmp_, weight_, level, flag_ );
   }

 private:
   DoFType      flag_;
   FunctionType tmp_;
   real_t       weight_;
};

class P1RichardsonSmoother : public Solver< P1ConstantLaplaceOperator >
{
 public:
   using FunctionType = typename P1ConstantLaplaceOperator::srcType;

   P1RichardsonSmoother( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : flag_( hyteg::Inner | hyteg::NeumannBoundary )
   , tmp_( "ric_tmp", storage, minLevel, maxLevel )
   , weight_( 1. )
   {}

   void setWeight( const real_t& weight ) { weight_ = weight; }

   void solve( const P1ConstantLaplaceOperator&                   A,
               const typename P1ConstantLaplaceOperator::srcType& x,
               const typename P1ConstantLaplaceOperator::dstType& b,
               const walberla::uint_t                                level ) override
   {
      tmp_.assign( {1.0}, {x}, level, hyteg::All );

      // update halos
      x.communicate< Face, Edge >( level );
      x.communicate< Edge, Face >( level );

      // dst = A x
      A.apply( tmp_, x, level, flag_, Replace );

      // dst = b - A x
      // (assigning to dst from dst should work here, because the operation is completely local)
      x.assign( {1., -1.}, {b, x}, level, flag_ );

      // x = x + weight (b - A x)
      // (assigning to x from dst should work here, because the operation is completely local)
      x.assign( {1., weight_}, {tmp_, x}, level, flag_ );
   }

 private:
   DoFType      flag_;
   FunctionType tmp_;
   real_t       weight_;
};

std::shared_ptr< Solver< P1ConstantLaplaceOperator > > setupSmoother( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                         const uint_t&                              minLevel,
                                                                         const uint_t&                              maxLevel,
                                                                         const std::string&                         smootherType,
                                                                         const uint_t&                              order,
                                                                         P1ConstantLaplaceOperator& laplaceOperator,
                                                                         P1Function< real_t >          function,
                                                                         P1Function< real_t >          tmpFunction )
{
   if ( smootherType == "richardson" )
   {
      auto       smoother       = std::make_shared< P1RichardsonSmoother >( storage, minLevel, maxLevel );
      const auto spectralRadius = estimateSpectralRadiusWithPowerIteration( laplaceOperator, function, tmpFunction, 40, storage, minLevel );
      smoother->setWeight( 1 / spectralRadius );
      return smoother;
   }
   else if ( smootherType == "jacobi" )
   {
      laplaceOperator.computeInverseDiagonalOperatorValues();
      auto smoother = std::make_shared< P1JacobiSmoother >( storage, minLevel, maxLevel );
      smoother->setWeight( real_c( 2. / 3. ) );
      return smoother;
   }
   else if ( smootherType == "chebyshev" )
   {
      laplaceOperator.computeInverseDiagonalOperatorValues();
      auto       smoother = std::make_shared< ChebyshevSmoother< P1ConstantLaplaceOperator > >( storage, minLevel, maxLevel );
      const auto spectralRadius = chebyshev::estimateRadius( laplaceOperator, minLevel, 100, storage, function, tmpFunction );
      smoother->setupCoefficients( order, spectralRadius );
      return smoother;
   }
   else if ( smootherType == "sor" )
   {
      auto smoother = std::make_shared< SORSmoother< P1ConstantLaplaceOperator > >( 1.0 );
      return smoother;
   }
   else
   {
      WALBERLA_ABORT( "unknown smoother type " << smootherType );
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./P1GMGSmootherComparison.prm" );
   }
   else
   {
      cfg = env.config();
   }

   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   parameters.listParameters();

   const uint_t minLevel         = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel         = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t max_outer_iter   = parameters.getParameter< uint_t >( "maxOuterIter" );
   const uint_t max_coarse_iter  = parameters.getParameter< uint_t >( "maxCoarseIter" );
   const real_t coarse_tolerance = parameters.getParameter< real_t >( "coarseTolerance" );
   const uint_t smoothingSteps   = parameters.getParameter< uint_t >( "smoothingSteps" );
   const uint_t order            = parameters.getParameter< uint_t >( "order" );

   const std::string smootherType = parameters.getParameter< std::string >( "smootherType" );

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( -1, -1 ), Point2D( 1., 1. ), MeshInfo::CRISSCROSS, 4, 4 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1Function< real_t > residual( "residual", storage, minLevel, maxLevel );
   P1Function< real_t > rightHandSide( "rightHandSide", storage, minLevel, maxLevel );
   P1Function< real_t > function( "function", storage, minLevel, maxLevel );
   P1Function< real_t > laplaceTimesFunction( "laplaceTimesFunction", storage, minLevel, maxLevel );

   P1ConstantLaplaceOperator laplaceOperator( storage, minLevel, maxLevel );

   const auto calculateResiduum = [&]() -> real_t {
      laplaceOperator.apply( function, laplaceTimesFunction, maxLevel, Inner );
      residual.assign( {1.0, -1.0}, {rightHandSide, laplaceTimesFunction}, maxLevel, Inner );
      return std::sqrt( residual.dotGlobal( residual, maxLevel, Inner ) );
   };

   const auto analytic_solution = []( const hyteg::Point3D& p ) -> real_t {
      return std::sin( 2 * pi * p[0] ) * std::cos( 3 * pi * p[1] );
   };

   const auto rhs = []( const hyteg::Point3D& p ) {
      return ( 4 + 9 ) * pi * pi * std::sin( 2 * pi * p[0] ) * std::cos( 3 * pi * p[1] );
   };

   function.interpolate( analytic_solution, maxLevel, DirichletBoundary );

   // assemble rhs:
   // we emulate this by interpolating the rhs pointwise and then multiply with the mass matrix
   P1ConstantMassOperator massOperator( storage, minLevel, maxLevel );

   residual.interpolate( rhs, maxLevel, DoFType::All );

   massOperator.apply( residual, rightHandSide, maxLevel, DoFType::All, UpdateType::Replace );

   // create smoothers
   P1Function< real_t > eigenvector( "eigenvector", storage, minLevel, maxLevel );
   P1Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   eigenvector.interpolate( analytic_solution, minLevel, DirichletBoundary );
   auto smoother = setupSmoother( storage, minLevel, maxLevel, smootherType, order, laplaceOperator, eigenvector, tmp );

   auto coarseGridSolver = std::make_shared< CGSolver< P1ConstantLaplaceOperator > >(
       storage, minLevel, minLevel, max_coarse_iter, coarse_tolerance );
   auto restrictionOperator  = std::make_shared< P1toP1LinearRestriction<> >();
   auto prolongationOperator = std::make_shared< P1toP1LinearProlongation<> >();

   auto multiGridSolver = GeometricMultigridSolver< P1ConstantLaplaceOperator >( storage,
                                                                                    smoother,
                                                                                    coarseGridSolver,
                                                                                    restrictionOperator,
                                                                                    prolongationOperator,
                                                                                    minLevel,
                                                                                    maxLevel,
                                                                                    smoothingSteps,
                                                                                    smoothingSteps );

   std::vector< real_t > resNormList;

   for ( uint_t i = 0; i < max_outer_iter; ++i )
   {
      multiGridSolver.solve( laplaceOperator, function, rightHandSide, maxLevel );
      const real_t resNorm = calculateResiduum();
      WALBERLA_LOG_INFO_ON_ROOT( "iter " << i << ": residual: " << resNorm );
      resNormList.push_back( resNorm );
   }

   if ( parameters.getParameter< bool >( "writeCSV" ) && walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      // construct a filename
      std::stringstream filename;
      filename << smootherType << "-";
      if ( smootherType == "chebyshev" )
         filename << order << "order-";
      filename << smoothingSteps << "steps-";
      filename << maxLevel << "level";
      filename << ".csv";

      // write the results to csv
      std::ofstream ofs( "csv/" + filename.str() );
      ofs << "iter,resNorm\n";
      for ( uint_t i = 0; i < resNormList.size(); i += 1 )
      {
         ofs << i << "," << resNormList[i] << "\n";
      }
   }
}
