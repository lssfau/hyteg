/*
 * Copyright (c) 2021 Andreas Wagner.
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
#include <utility>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

template < typename value_t >
class P1CahnHilliardFunction : public BlockFunction< value_t >
{
 public:
   P1CahnHilliardFunction( const std::string&                         name,
                           const std::shared_ptr< PrimitiveStorage >& storage,
                           size_t                                     minLevel,
                           size_t                                     maxLevel )
   : BlockFunction< value_t >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_mu", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_phi", storage, minLevel, maxLevel ) );
   };

   [[nodiscard]] const P1Function< value_t >& getMu() const
   {
      return this->getSubFunction( 0 ).template unwrap< P1Function< value_t > >();
   }

   [[nodiscard]] const P1Function< value_t >& getPhi() const
   {
      return this->getSubFunction( 1 ).template unwrap< P1Function< value_t > >();
   }

   // TODO: Move this into BlockFunction
   template < typename OtherFunctionValueType >
   void copyBoundaryConditionFromFunction( const P1CahnHilliardFunction< OtherFunctionValueType >& other )
   {
      for ( uint_t k = 0; k < BlockFunction< value_t >::getNumberOfBlocks(); ++k )
      {
         BlockFunction< value_t >::getSubFunction( k ).setBoundaryCondition( other.getSubFunction( k ).getBoundaryCondition() );
      }
   }
};

class P1CahnHilliardFunctionTag
{};

template < typename VType >
struct FunctionTrait< P1CahnHilliardFunction< VType > >
{
   typedef VType                     ValueType;
   typedef P1CahnHilliardFunctionTag Tag;

   static std::string getTypeName() { return "CahnHilliardFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

using chType = P1CahnHilliardFunction< real_t >;

using MassFormType    = hyteg::P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >;
using LaplaceFormType = hyteg::P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >;

std::shared_ptr< BlockOperator< chType, chType > > create_lhs_operator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                        uint_t                                     minLevel,
                                                                        uint_t                                     maxLevel,
                                                                        real_t                                     dt,
                                                                        real_t                                     epsilon )
{
   P1LinearCombinationForm form00;
   form00.addOwnedForm< LaplaceFormType >( dt );

   P1LinearCombinationForm form11;
   form11.addOwnedForm< LaplaceFormType >( -std::pow( epsilon, 2 ) );
   form11.addOwnedForm< MassFormType >( -2 );

   auto op = std::make_shared< BlockOperator< chType, chType > >( storage, minLevel, maxLevel, 2, 2 );
   op->createSubOperator< P1ConstantLinearCombinationOperator >( 0, 0, storage, minLevel, maxLevel, form00 );
   op->createSubOperator< P1ConstantLinearCombinationOperator >( 1, 1, storage, minLevel, maxLevel, form11 );
   op->createSubOperator< P1ConstantMassOperator >( 0, 1, storage, minLevel, maxLevel );
   op->createSubOperator< P1ConstantMassOperator >( 1, 0, storage, minLevel, maxLevel );

   return op;
}

class RHSAssembler
{
 public:
   RHSAssembler( std::shared_ptr< PrimitiveStorage > storage, uint_t minLevel, uint_t maxLevel, real_t epsilon )
   : storage_( std::move( storage ) )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , epsilon_( epsilon )
   {}

   void assemble( const chType& u_prev, const chType& rhs, uint_t level, DoFType flag )
   {
      auto kappa = [&]( const Point3D& p ) -> real_t {
         real_t phi_prev_value;
         WALBERLA_CHECK( u_prev.getPhi().evaluate( p, level, phi_prev_value ) );

         return -3. + std::pow( phi_prev_value, 2 );
      };

      forms::p1_k_mass_affine_q4 form( kappa, kappa );

      P1ElementwiseKMassOperator op_phi( storage_, minLevel_, maxLevel_, form );

      op_phi.apply( u_prev.getPhi(), rhs.getPhi(), level, flag );

      P1ElementwiseMassOperator op_mu( storage_, minLevel_, maxLevel_ );

      op_mu.apply( u_prev.getPhi(), rhs.getMu(), level, flag );
   }

 protected:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   real_t epsilon_;
};

std::shared_ptr< hyteg::PrimitiveStorage > create_storage()
{
   hyteg::Point2D lowerRight( { -1, -1 } );
   hyteg::Point2D upperLeft( { +1, +1 } );

   auto meshInfo     = hyteg::MeshInfo::meshRectangle( lowerRight, upperLeft, hyteg::MeshInfo::CROSS, 2, 2 );
   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // Neumann BC everywhere
   setupStorage->setMeshBoundaryFlagsOnBoundary( 2, 0, true );

   return std::make_shared< hyteg::PrimitiveStorage >( *setupStorage );
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "./CahnHilliard.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   parameters.listParameters();

   const uint_t minLevel = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = parameters.getParameter< uint_t >( "maxLevel" );

   const real_t tau     = 5e-3;
   const real_t epsilon = 4e-2;

   const real_t t_end = 1;

   const uint_t max_time_steps = 300;

   const auto storage = create_storage();

   chType u_now( "u_now", storage, minLevel, maxLevel );
   chType u_prev( "u_prev", storage, minLevel, maxLevel );
   chType rhs( "rhs", storage, minLevel, maxLevel );

   auto solver = std::make_shared< MinResSolver< BlockOperator< chType, chType > > >( storage, minLevel, maxLevel );

   auto lhsOp = create_lhs_operator( storage, minLevel, maxLevel, tau, epsilon );

   RHSAssembler rhsAssembler( storage, minLevel, maxLevel, epsilon );

   // one species has a slight majority
   u_prev.getPhi().interpolate(
       []( const Point3D& ) { return walberla::math::realRandom( real_t( -0.01 ), real_t( 0.01 ) ); }, maxLevel, All );
   u_prev.getMu().interpolate( 0, maxLevel, All );

   real_t t_now = 0;

   // output
   VTKOutput vtkOutput( ".", "CahnHilliard", storage );

   for ( uint_t k = 0; k < max_time_steps; ++k )
   {
      t_now += tau;

      if ( t_now > t_end )
         break;

      WALBERLA_LOG_INFO( "iter = " << k << " t = " << t_now );

      rhsAssembler.assemble( u_prev, rhs, maxLevel, All );

      solver->solve( *lhsOp, u_now, rhs, maxLevel );

      u_prev.assign( { 1. }, { u_now }, maxLevel, All );

      // vtk output
      vtkOutput.add( u_now );
      vtkOutput.write( maxLevel, k );
   }
}
