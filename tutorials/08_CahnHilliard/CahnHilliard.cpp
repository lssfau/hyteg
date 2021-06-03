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

/**
 * \page 08_CahnHilliard A full app for the Cahn-Hilliard equation.
 *
 * \dontinclude tutorials/08_CahnHilliard/CahnHilliard.cpp
 *
 * \brief In this tutorial we will use the LinearCombinationForm to define specialized bilinear forms for the Cahn-Hilliard (CH) equations.
 *        In combination with block-functions and block-operators we will discretize the CH-Equations.
 *        We will speed up our solvers by implementing a simple diagonal preconditioner.
 *
 * We already showed how to solve the Poisson- and the Stokes-Equation, which already are well supported by HyTeG.
 * In this tutorial, we will implement a solver for a completely unrelated problem, which is not a core application for HyTeG (yet ;)).
 * In addition, we will learn about block-functions, block-operators and the composition of finite-element forms.
 *
 * \section Problem
 *
 * As an example we will try to solve the Cahn-Hilliard equation,
 * which simulates the phase separation of a mixture of two fluids.
 * For convenience, lets call these fluids A and B.
 * They are represented by an indicator function \f$\phi \in [-1, +1]\f$,
 * such that \f$\phi(x) = +1\f$ if only fluid A is present and \f$\phi(x) = -1\f$  if only fluid B is present.
 * All intermediate values correspond to mixtures of the fluid, where either A (\f$\phi(x) < 1\f$) or B dominates (\f$\phi(x) > 1\f$).
 *
 * We want to assign an energy functional to these equations by given by
 * \f{align*}{
 *     F(\phi) = \int \boldsymbol \psi(\phi) dx + \frac{\epsilon^2}{2}\int |\nabla \phi |^2 dx
 * \f}
 * where \f$\psi(\phi) = \frac{1}{4}(1-\phi^2)^2\f$ is a double well potential, which penalizes mixed states,
 * while the second term models a surface energy and encourages small interfaces.
 * A gradient-flow for this functional is given by the Cahn-Hilliard equations
 * \f{align*}{
 *     \partial_t \phi &= \Delta \mu \\
 *     \mu &= \boldsymbol \psi'(\phi) - \epsilon^2 \Delta \phi
 *     ,
 * \f}
 * which minimizes \f$F(\phi(t))\f$ as time passes, i.e. \f$F(\phi(t)) \leq F(\phi(s))\f$ for \f$t \leq s\f$.
 * The equations are a mixed system of the fields \f$\phi\f$ and \f$\mu\f$.
 *
 * \section Time-discretization
 *
 * The energy minimization of $F$ is usually something we want to enforce by our time discretization.
 *
 * For this we split up the potential $\boldsymbol \psi$ into a convex and a concave part, such that
 * \f{align*}{
 *     \boldsymbol \psi = \boldsymbol \psi_c - \boldsymbol \psi_e
 * \f}
 * where \f$\boldsymbol \psi_e\f$, \f$\boldsymbol \psi_c\f$ are convex.
 * It can be shown that if we treat \f$\psi_e\f$ implicitly and \f$\psi_c\f$ explicitly the resulting scheme will be unconditionally gradient stable.
 * That means that \f$F(\phi^{n+1}) \leq F(\phi^{n})\f$ for our time-discrete solutions \f$\phi^n\f$ independent of the time-step size \f$\tau\f$.
 * The resulting scheme is
 *
 * \f{align}{
 *     \phi^{(n+1)}-\phi^{(n)} &= \tau \Delta \mu \\
 *     \mu &= \boldsymbol \psi_c'(\phi^{(n+1)}) - \boldsymbol \psi_e'(\phi^{(n)}) - \epsilon^2 \Delta \phi^{(n+1)}
 *     .
 * \f}
 *
 * The decomposition of the potential in a convex and concave function is typically not unique.
 * A possible choice is
 * \f{align}{
 *     \boldsymbol \psi_c(\phi) = \frac{1}{4}
 *     \left( 4 \phi^2 + 1 \right)
 *     \quad \text{ and } \quad
 *     \boldsymbol \psi_e(\phi) = \frac{1}{4}
 *     \left( 6 \phi^2 - \phi^4 \right)
 * \f}
 * yielding the derivatives
 * \f{align}{
 *     \boldsymbol \psi_c'(\phi) =
 *     2 \phi
 *     \quad \text{ and } \quad
 *     \boldsymbol \psi_e'(\phi) =
 *     3 \phi - \phi^3
 *     .
 * \f}
 *
 * \section Space-discretization
 *
 * Applying a FEM discretization on the equations yields the following system of equations
 * \f{align}{
 *     \begin{pmatrix}
 *         \tau K & M\\
 *         M & - \epsilon^2 K - 2 M
 *     \end{pmatrix}
 *     \begin{pmatrix}
 *         \boldsymbol \mu \\
 *         \boldsymbol{\phi}^{(n+1)}
 *     \end{pmatrix}
 *     =
 *     \begin{pmatrix}
 *         \boldsymbol f^{(n)} \\ \boldsymbol g^{(n)}
 *     \end{pmatrix}
 *     ,
 * \f}
 * with the stiffness matrix \f$K\f$, the mass matrix \f$M\f$ and the right-hand side given by
 * \f{align}{
 *     \boldsymbol{f}^{(n)}_i = \left\langle\phi^n_h, \nu_{h,i}\right\rangle_{L_2}
 * \f}
 * and
 * \f{align}{
 *     \boldsymbol g^{(n)}_i = \left\langle-3 \phi^n_h + (\phi^n_h)^3, \psi_{h,i} \right\rangle_{L_2}
 * \f}
 * with test functions \f$\nu_{h,i}\f$ and \f$\psi_{h,i}\f$.
 * We consider a \f$\mathcal P^1\f$ (continuous, piecewise linear polynomials) discretization of FE-shape functions.
 *
 *
 * \section Preconditioner
 *
 * To solve our problem efficiently we will need a preconditioner.
 * Here, we will follow the approach in [BRENNER2018] and get simply by rescaling the equations
 * \f{align}{
 *     \begin{pmatrix}
 *         \tau K & M\\
 *         M & - \epsilon^2 K - 2 M
 *     \end{pmatrix}^{-1}
 *     =
 *     \begin{pmatrix}
 *         \tau^{-1/4} & \\
 *         & \tau^{1/4}
 *     \end{pmatrix}
 *     \begin{pmatrix}
 *         \sqrt{\tau} K & M\\
 *         M & - \sqrt \tau \epsilon^2 K - 2 \sqrt \tau M
 *     \end{pmatrix}^{-1}
 *     \begin{pmatrix}
 *         \tau^{-1/4} & \\
 *         & \tau^{1/4}
 *     \end{pmatrix}
 *     .
 * \f}
 * We define
 * \f{align}{
 *     B:=
 *     \begin{pmatrix}
 *         \sqrt{\tau} K & M\\
 *         M & - \sqrt \tau \epsilon^2 K - 2 \sqrt \tau M
 *     \end{pmatrix}
 *     .
 * \f}
 * As shown in [Brenner2018] \f$P^{-1} B\f$ can be solved efficiently with MINRES, where
 * \f{align}{
 *     P^{-1}
 *     :=
 *     \begin{pmatrix}
 *         (\sqrt{\tau} K + M)^{-1} & \\
 *         &
 *         (\sqrt{\tau} \epsilon^2 K + M)^{-1}
 *     \end{pmatrix}
 *     .
 * \f}
 *
 * To get an optimal solver, we can use a geometric multigrid to invert the diagonal blocks.
 *
 * \section Implementation
 *
 * We start with defining a block function representing the state of our system at a given time step by inheriting from the hyteg::BlockFunction.
 *
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction definition
 *
 * The two components for the \f$\phi \f$ and \f$\mu\f$ are defined inside the constructor by pushing a hyteg::P1Function onto
 * `subFunc_` vector.
 *
 * For our own convenience we define two helper functions `getMu` and `getPhi`, which return the subcomponents of our vector.
 *
 * TODO: Write why we need hyteg::FunctionWrapper
 *
 * The compiler needs to be able to infer some static information about each function.
 * For instance its name, its value type and it needs a tag for certain meta programming features.
 * This can be achieved by extending the given trait system:
 *
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction meta-information
 *
 * If now one of HyTeG's internal algorithms needs to now the value type of a unknown function type `UnknownType`, it can do this by calling
 * `hyteg::FunctionTrait< UnknownType >::ValueType`.
 * This is implemented for hyteg::P1Function, hyteg::P2Function, hyteg::P2P1TaylorHoodFunction, ... and now also for our
 * `CahnHilliardFunction`.
 *
 * A more in-depth tutorial about traits can be found at [TRAITS].
 *
 *
 *
 * TODO: Write the rest
 *
 * The full code is:
 *
 * \include tutorials/08_CahnHilliard/CahnHilliard.cpp
 *
 * \section References
 *
 * [Brenner2018] BRENNER, Susanne C.; DIEGEL, Amanda E.; SUNG, Li-Yeng.
 *               A robust solver for a mixed finite element method for the Cahn–Hilliard equation.
 *               Journal of Scientific Computing, 2018, 77. Jg., Nr. 2, S. 1234-1249.
 *
 * [TRAITS] https://www.boost.org/doc/libs/1_53_0/libs/geometry/doc/html/geometry/design.html
 *
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

/// [CahnHilliardFunction definition]
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
/// [CahnHilliardFunction definition]

/// [CahnHilliardFunction meta-information]
class P1CahnHilliardFunctionTag
{};

namespace hyteg {

template < typename VType >
struct FunctionTrait< P1CahnHilliardFunction< VType > >
{
   typedef VType                     ValueType;
   typedef P1CahnHilliardFunctionTag Tag;

   static std::string getTypeName() { return "CahnHilliardFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

} // namespace hyteg
/// [CahnHilliardFunction meta-information]

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

std::shared_ptr< MinResSolver< BlockOperator< chType, chType > > >
    create_solver( const std::shared_ptr< PrimitiveStorage >& storage,
                   uint_t                                     minLevel,
                   uint_t                                     maxLevel,
                   real_t                                     tau,
                   real_t                                     epsilon )
{
   WALBERLA_UNUSED( tau );
   WALBERLA_UNUSED( epsilon );
   return std::make_shared< MinResSolver< BlockOperator< chType, chType > > >( storage, minLevel, maxLevel );
}

std::shared_ptr< BlockOperator< chType, chType > > create_lhs_operator_v2( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           uint_t                                     minLevel,
                                                                           uint_t                                     maxLevel,
                                                                           real_t                                     tau,
                                                                           real_t                                     epsilon )
{
   P1LinearCombinationForm form00;
   form00.addOwnedForm< LaplaceFormType >( std::sqrt( tau ) );

   P1LinearCombinationForm form11;
   form11.addOwnedForm< LaplaceFormType >( -std::pow( epsilon, 2 ) * std::sqrt( tau ) );
   form11.addOwnedForm< MassFormType >( -2 * std::sqrt( tau ) );

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

class CahnHilliardEvolutionOperator
{
 public:
   CahnHilliardEvolutionOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     minLevel,
                                  uint_t                                     maxLevel,
                                  real_t                                     tau,
                                  real_t                                     epsilon )
   : epsilon_( epsilon )
   , tau_( tau )
   , lhsOp_( create_lhs_operator_v2( storage, minLevel, maxLevel, tau, epsilon ) )
   , solver_( create_solver( storage, minLevel, maxLevel, tau, epsilon ) )
   , rhsAssembler_( storage, minLevel, maxLevel, epsilon )
   , rhs_( "rhs", storage, minLevel, maxLevel )
   {}

   void apply( const chType& u_old, const chType& u_now, uint_t level )
   {
      rhsAssembler_.assemble( u_old, rhs_, level, All );

      scale( rhs_, rhs_, level, All );

      // solve the system
      u_now.assign( { 1 }, { u_old }, level, All );
      solver_->solve( *lhsOp_, u_now, rhs_, level );

      scale( u_now, u_now, level, All );
   }

 protected:
   void scale( const chType& src, const chType& dst, uint_t level, DoFType flag ) const
   {
      const real_t factor = std::pow( tau_, 0.25 );
      dst.getMu().assign( { 1. / factor }, { src.getMu() }, level, flag );
      dst.getPhi().assign( { factor }, { src.getPhi() }, level, flag );
   }

 protected:
   real_t epsilon_;
   real_t tau_;

   std::shared_ptr< BlockOperator< chType, chType > > lhsOp_;

   std::shared_ptr< MinResSolver< BlockOperator< chType, chType > > > solver_;

   RHSAssembler rhsAssembler_;

   chType rhs_;
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

   CahnHilliardEvolutionOperator chOperator( storage, minLevel, maxLevel, tau, epsilon );

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
      chOperator.apply( u_prev, u_now, maxLevel );

      u_prev.assign( { 1. }, { u_now }, maxLevel, All );

      // vtk output
      vtkOutput.add( u_now );
      vtkOutput.write( maxLevel, k );
   }
}
