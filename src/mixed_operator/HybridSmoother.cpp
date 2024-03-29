/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "mixed_operator/HybridSmoother.hpp"

#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/gridtransferoperators/N1E1toP1Lifting.hpp"
#include "hyteg/gridtransferoperators/P1toN1E1Gradient.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

namespace hyteg {
namespace n1e1 {

template < class N1E1OperatorType, class P1LaplaceOperatorType >
HybridSmoother< N1E1OperatorType, P1LaplaceOperatorType >::HybridSmoother(
    const std::shared_ptr< PrimitiveStorage >&         storage,
    const n1e1::N1E1VectorFunction< real_t >&          tmpFunction,
    std::shared_ptr< P1LaplaceOperatorType >           p1LaplaceOperator,
    std::shared_ptr< Solver< N1E1OperatorType > >      n1e1Smoother,
    std::shared_ptr< Solver< P1LaplaceOperatorType > > p1Smoother,
    const uint_t                                       minLevel,
    const uint_t                                       maxLevel,
    const uint_t                                       n1e1SmoothSteps,
    const uint_t                                       p1SmoothSteps )
: n1e1SmoothSteps_( n1e1SmoothSteps )
, p1SmoothSteps_( p1SmoothSteps )
, flag_( hyteg::Inner | hyteg::NeumannBoundary )
, p1LaplaceOperator_( p1LaplaceOperator )
, n1e1Smoother_( n1e1Smoother )
, p1Smoother_( p1Smoother )
, vectorResidual_( tmpFunction )
, scalarResidual_( "hybridSmootherScalarResidual", storage, minLevel, maxLevel )
, scalarPotential_( "hybridSmootherScalarPotential", storage, minLevel, maxLevel )
, timingTree_( storage->getTimingTree() )
{}

template < class N1E1OperatorType, class P1LaplaceOperatorType >
HybridSmoother< N1E1OperatorType, P1LaplaceOperatorType >::HybridSmoother(
    const std::shared_ptr< PrimitiveStorage >&         storage,
    std::shared_ptr< P1LaplaceOperatorType >           p1LaplaceOperator,
    std::shared_ptr< Solver< N1E1OperatorType > >      n1e1Smoother,
    std::shared_ptr< Solver< P1LaplaceOperatorType > > p1Smoother,
    const uint_t                                       minLevel,
    const uint_t                                       maxLevel,
    const uint_t                                       n1e1SmoothSteps,
    const uint_t                                       p1SmoothSteps )
: HybridSmoother( storage,
                  N1E1VectorFunction< real_t >( "hybridSmootherVectorResidual", storage, minLevel, maxLevel ),
                  p1LaplaceOperator,
                  n1e1Smoother,
                  p1Smoother,
                  minLevel,
                  maxLevel,
                  n1e1SmoothSteps,
                  p1SmoothSteps )
{}

template < class N1E1OperatorType, class P1LaplaceOperatorType >
void HybridSmoother< N1E1OperatorType, P1LaplaceOperatorType >::solve( const N1E1OperatorType&             A,
                                                                       const N1E1VectorFunction< real_t >& x,
                                                                       const N1E1VectorFunction< real_t >& b,
                                                                       const uint_t                        level )
{
   timingTree_->start( "Hybrid Smoother" );

   vectorResidual_.copyBoundaryConditionFromFunction( x );
   scalarResidual_.setBoundaryCondition( x.getBoundaryCondition() );
   // By lifting to potential space Dirichlet BCs become Neumann BCs.
   // Since there is currently no easy way to modify the BCs of `x` accordingly, we handle all DoFs like Inner DoFs.
   // This is Ok, since Neumann and Inner DoFs are treated equally.
   scalarPotential_.setBoundaryCondition( BoundaryCondition::createAllInnerBC() );

   timingTree_->start( "Smoother in N(curl)^⊥" );
   for ( uint_t i = 0; i < n1e1SmoothSteps_; ++i )
   {
      n1e1Smoother_->solve( A, x, b, level );
   }
   timingTree_->stop( "Smoother in N(curl)^⊥" );

   A.apply( x, vectorResidual_, level, flag_ );
   vectorResidual_.assign( { 1.0, -1.0 }, { b, vectorResidual_ }, level, flag_ );

   timingTree_->start( "Lifting" );
   N1E1toP1Lifting( vectorResidual_, scalarResidual_, level, flag_ );
   timingTree_->stop( "Lifting" );

   scalarPotential_.setToZero( level );

   timingTree_->start( "Smoother in N(curl)" );
   for ( uint_t i = 0; i < p1SmoothSteps_; ++i )
   {
      p1Smoother_->solve( *p1LaplaceOperator_, scalarPotential_, scalarResidual_, level );
   }
   timingTree_->stop( "Smoother in N(curl)" );

   timingTree_->start( "Gradient" );
   P1toN1E1Gradient( scalarPotential_, vectorResidual_, level, flag_ );
   timingTree_->stop( "Gradient" );

   x.add( { 1.0 }, { vectorResidual_ }, level, flag_ );

   timingTree_->stop( "Hybrid Smoother" );
}

template class HybridSmoother< N1E1ElementwiseLinearCombinationOperator, P1ConstantLaplaceOperator >;
template class HybridSmoother< N1E1ElementwiseLinearCombinationOperator, P1ConstantLinearCombinationOperator >;
template class HybridSmoother< N1E1ElementwiseLinearCombinationOperator, P1ElementwiseBlendingLaplaceOperator >;
template class HybridSmoother< N1E1ElementwiseCurlCurlPlusMassOperatorQ2, P1ConstantLaplaceOperator >;
template class HybridSmoother< N1E1ElementwiseBlendingCurlCurlPlusMassOperatorQ2, P1ElementwiseBlendingLaplaceOperatorQ2 >;
template class HybridSmoother< N1E1ElementwiseBlendingCurlCurlPlusMassOperatorQ2, P1ElementwiseBlendingLaplaceOperator >;

} // namespace n1e1
} // namespace hyteg
