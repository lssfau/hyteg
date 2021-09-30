/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#pragma once

#include <core/math/MatrixMxN.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

enum class DiffusionTimeIntegrator
{
   ImplicitEuler,
   CrankNicolson
};

/// \brief Unsteady diffusion operator for a P2 finite element discretization of the solution.
///
/// This is a composite operator used to solve the parabolic PDE that describes unsteady diffusion:
///
/// \f$ \frac{\partial u}{\partial t} = D  \Delta u + f\f$
///
/// where \f$ D \f$ is a constant diffusion parameter.
///
/// This operator can be used for an implicit Euler or Crank-Nicolson time-integrator.
/// In particular it equals
///
/// \f$ M + \theta\ dt\ D\ \mathcal{L} \f$
///
/// where \f$ D \f$ is again the diffusivity constant, \f$ M \f$ the finite element mass matrix, \f$ dt \f$ a constant time step,
/// \f$ \mathcal{L} = - \Delta \f$ the negative Laplacian, and \f$ \theta \in \{0.5, 1\} \f$ for Crank-Nicolson or implicit Euler.
///
/// To solve the unsteady diffusion equation, see UnsteadyDiffusion.
template < typename FunctionType,
           template < class >
           class Operator_T,
           typename LaplaceForm_T,
           typename MassForm_T,
           typename LinearCombinationForm_T >
class UnsteadyDiffusionOperator : public Operator< FunctionType, FunctionType >,
                                  public SORSmoothable< FunctionType >,
                                  public SORBackwardsSmoothable< FunctionType >,
                                  public GSSmoothable< FunctionType >,
                                  public GSBackwardsSmoothable< FunctionType >,
                                  public WeightedJacobiSmoothable< FunctionType >
{
 public:
   UnsteadyDiffusionOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                              const uint_t&                              minLevel,
                              const uint_t&                              maxLevel,
                              const real_t&                              dt,
                              const real_t&                              diffusionCoefficient,
                              const DiffusionTimeIntegrator&             timeIntegrator )
   : Operator< FunctionType, FunctionType >( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , diffusionCoefficient_( diffusionCoefficient )
   , laplaceForm_( new LaplaceForm_T() )
   , massForm_( new MassForm_T() )
   , timeIntegrator_( timeIntegrator )
   , dt_( dt )
   {
      setDt( dt );
   }

   void apply( const FunctionType& src,
               const FunctionType& dst,
               const uint_t&       level,
               const DoFType&      flag,
               const UpdateType&   updateType = Replace ) const
   {
      unsteadyDiffusionOperator_->apply( src, dst, level, flag, updateType );
   }

   void smooth_sor( const FunctionType& dst, const FunctionType& rhs, real_t relax, size_t level, DoFType flag ) const override
   {
      smooth_sor( dst, rhs, relax, level, flag, false );
   }

   void smooth_sor_backwards( const FunctionType& dst,
                              const FunctionType& rhs,
                              real_t              relax,
                              size_t              level,
                              DoFType             flag ) const override
   {
      smooth_sor( dst, rhs, relax, level, flag, true );
   }

   void smooth_sor( const FunctionType& dst,
                    const FunctionType& rhs,
                    real_t              relax,
                    size_t              level,
                    DoFType             flag,
                    const bool&         backwards ) const
   {
      if ( !backwards )
      {
         if ( const auto* op_sor = dynamic_cast< const SORSmoothable< FunctionType >* >( unsteadyDiffusionOperator_.get() ) )
         {
            op_sor->smooth_sor( dst, rhs, relax, level, flag );
         }
         else
         {
            throw std::runtime_error( "The UnsteadyDiffusionOperator requires the SORSmoothable interface." );
         }
      }
      else
      {
         if ( const auto* op_sor =
                  dynamic_cast< const SORBackwardsSmoothable< FunctionType >* >( unsteadyDiffusionOperator_.get() ) )
         {
            op_sor->smooth_sor_backwards( dst, rhs, relax, level, flag );
         }
         else
         {
            throw std::runtime_error( "The UnsteadyDiffusionOperator requires the SORBackwardsSmoothable interface." );
         }
      }
   }

   void smooth_gs( const FunctionType& dst, const FunctionType& rhs, size_t level, DoFType flag, const bool& backwards ) const
   {
      smooth_sor( dst, rhs, 1.0, level, flag, backwards );
   }

   void smooth_gs( const FunctionType& dst, const FunctionType& rhs, size_t level, DoFType flag ) const override
   {
      smooth_gs( dst, rhs, level, flag, false );
   }

   void smooth_gs_backwards( const FunctionType& dst, const FunctionType& rhs, size_t level, DoFType flag ) const override
   {
      smooth_gs( dst, rhs, level, flag, true );
   }

   void smooth_jac( const FunctionType& dst,
                    const FunctionType& rhs,
                    const FunctionType& src,
                    real_t              relax,
                    size_t              level,
                    DoFType             flag ) const override
   {
      if ( const auto* op_jac =
               dynamic_cast< const WeightedJacobiSmoothable< FunctionType >* >( unsteadyDiffusionOperator_.get() ) )
      {
         op_jac->smooth_jac( dst, rhs, src, relax, level, flag );
      }
      else
      {
         throw std::runtime_error( "The UnsteadyDiffusionOperator requires the WeightedJacobiSmoothable interface." );
      }
   }

   real_t dt() const { return dt_; }

   void setDt( const real_t& dt )
   {
      dt_                        = dt;
      unsteadyDiffusionOperator_ = std::make_shared< Operator_T< LinearCombinationForm_T > >(
          storage_,
          minLevel_,
          maxLevel_,
          LinearCombinationForm_T( {1.0, dtScaling() * dt * diffusionCoefficient_}, {massForm_.get(), laplaceForm_.get()} ) );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                                                    mat,
                  const typename Operator_T< LinearCombinationForm_T >::srcType::template FunctionType< idx_t >& src,
                  const typename Operator_T< LinearCombinationForm_T >::srcType::template FunctionType< idx_t >& dst,
                  size_t                                                                                         level,
                  DoFType                                                                                        flag ) const
   {
      unsteadyDiffusionOperator_->toMatrix( mat, src, dst, level, flag );
   }

   DiffusionTimeIntegrator getTimeIntegrator() const { return timeIntegrator_; }

   const Operator_T< LinearCombinationForm_T >& getOperator() const { return *unsteadyDiffusionOperator_; }
   Operator_T< LinearCombinationForm_T >&       getOperator() { return *unsteadyDiffusionOperator_; }

 private:
   real_t dtScaling()
   {
      switch ( timeIntegrator_ )
      {
      case DiffusionTimeIntegrator::ImplicitEuler:
         return 1.0;
      case DiffusionTimeIntegrator::CrankNicolson:
         return 0.5;
      default:
         WALBERLA_ABORT( "Invalid time integrator" )
      }
   }

   std::shared_ptr< PrimitiveStorage >                      storage_;
   uint_t                                                   minLevel_;
   uint_t                                                   maxLevel_;
   real_t                                                   diffusionCoefficient_;
   std::shared_ptr< LaplaceForm_T >                         laplaceForm_;
   std::shared_ptr< MassForm_T >                            massForm_;
   DiffusionTimeIntegrator                                  timeIntegrator_;
   std::shared_ptr< Operator_T< LinearCombinationForm_T > > unsteadyDiffusionOperator_;

   real_t dt_;
};

template < typename P1Form >
using P1ConstantOperatorSingleTemplateParamter = P1ConstantOperator< P1Form, false, false, false >;

typedef UnsteadyDiffusionOperator<
    P1Function< real_t >,
    P1ConstantOperatorSingleTemplateParamter,
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >,
    P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
    P1LinearCombinationForm >
    P1ConstantUnsteadyDiffusionOperator;

typedef UnsteadyDiffusionOperator<
    P2Function< real_t >,
    P2ConstantOperator,
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >,
    P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise >,
    P2LinearCombinationForm >
    P2ConstantUnsteadyDiffusionOperator;

/// This unsteady diffusion operator supports blending.
typedef UnsteadyDiffusionOperator< P2Function< real_t >,
                                   P2ElementwiseOperator,
                                   P2Form_laplace,
                                   P2Form_mass,
                                   P2LinearCombinationForm >
    P2ElementwiseUnsteadyDiffusionOperator;

/// \brief Wrapper class to solve the unsteady diffusion equation in time.
///
/// \tparam FunctionType function discretization type (e.g. P2Function< real_t >)
/// \tparam UnsteadyDiffusionOperatorType type of the unsteady diffusion operator, must match the function type (e.g. P2UnsteadyDiffusionOperator)
/// \tparam LaplaceOperatorType Laplace operator matrix type matching the function (e.g. P2ConstantLaplaceOperator)
/// \tparam MassOperatorType mass matrix type matching the function (e.g. P2ConstantMassOperator)
///
/// Performs the implicit Euler or Crank-Nicholson method to advance the solution of the unsteady diffusion equation.
/// Therefore in each time step, the passed UnsteadyDiffusionOperatorType instance must be inverted.
///
/// Note that the selection of the time-integrator is done during the construction of the diffusion operator.
template < typename FunctionType,
           typename UnsteadyDiffusionOperatorType,
           typename LaplaceOperatorType,
           typename MassOperatorType >
class UnsteadyDiffusion
{
 public:
   UnsteadyDiffusion( const std::shared_ptr< PrimitiveStorage >&                        storage,
                      const uint_t&                                                     minLevel,
                      const uint_t&                                                     maxLevel,
                      const std::shared_ptr< Solver< UnsteadyDiffusionOperatorType > >& diffusionSolver )
   : storage_( storage )
   , uOld_( "uOld", storage, minLevel, maxLevel )
   , fWeak_( "fWeak", storage, minLevel, maxLevel )
   , solver_( diffusionSolver )
   {}

   void setSolver( const std::shared_ptr< Solver< UnsteadyDiffusionOperatorType > >& diffusionSolver )
   {
      solver_ = diffusionSolver;
   }

   /// \brief Performs one implicit Euler step to advance the solution of the PDE from time step n to n+1.
   ///
   /// \param A the unsteady diffusion operator instance (dt, diffusivity and time integrator are defined by this operator)
   /// \param L the Laplacian operator
   /// \param M the finite element mass operator
   /// \param u the solution function of the next time step n+1
   ///          (must have corresponding BCs of the next time step interpolated on the boundary)
   /// \param uOld the solution of the previous time step n
   ///             (the interpolated Dirichlet BCs must be equal to those of the previous time step)
   /// \param f the right hand side of time step n+1
   ///          (in "strong" form - it is multiplied by the mass in this method)
   /// \param fOld the right hand side of time step n
   ///             (in "strong" form - it is multiplied by the mass in this method)
   /// \param level the refinement level
   /// \param flag where to solve
   ///
   /// This method solves for \f$ u^{n+1}\f$:
   ///
   /// \f$ (M + \theta\ dt\ D\ \mathcal{L}) u^{n+1} := (M - (1 - \theta)\ dt\ D\ \mathcal{L}) u^n + \theta\ dt\ M f^{n+1} + (1-\theta)\ dt\ M f^{n} \f$
   ///
   /// where \f$ A = M + \theta\ dt\ D\ \mathcal{L} \f$ the unsteady diffusion operator (see e.g. P2UnsteadyDiffusionOperator).
   ///
   void step( const UnsteadyDiffusionOperatorType& A,
              const LaplaceOperatorType&           L,
              const MassOperatorType&              M,
              const FunctionType&                  u,
              const FunctionType&                  uOld,
              const FunctionType&                  f,
              const FunctionType&                  fOld,
              const uint_t&                        level,
              const DoFType&                       flag )
   {
      uOld_.copyBoundaryConditionFromFunction( u );
      fWeak_.copyBoundaryConditionFromFunction( u );

      if ( A.getTimeIntegrator() == DiffusionTimeIntegrator::ImplicitEuler )
      {
         // implicit Euler
         M.apply( f, fWeak_, level, flag );
         M.apply( uOld, uOld_, level, flag );
         uOld_.assign( {1.0, A.dt()}, {uOld_, fWeak_}, level, flag );
         solver_->solve( A, u, uOld_, level );
      }
      else if ( A.getTimeIntegrator() == DiffusionTimeIntegrator::CrankNicolson )
      {
         // Crank-Nicholson
         M.apply( f, fWeak_, level, flag );
         M.apply( fOld, fWeak_, level, flag, Add );
         M.apply( uOld, uOld_, level, flag );
         uOld_.assign( {1.0, 0.5 * A.dt()}, {uOld_, fWeak_}, level, flag );
         L.apply( uOld, fWeak_, level, flag );
         uOld_.assign( {1.0, -0.5 * A.dt()}, {uOld_, fWeak_}, level, flag );
         solver_->solve( A, u, uOld_, level );
      }
   }

   /// \brief Same step implementation but for rhs == 0.
   void step( const UnsteadyDiffusionOperatorType& A,
              const LaplaceOperatorType&           L,
              const MassOperatorType&              M,
              const FunctionType&                  u,
              const FunctionType&                  uOld,
              const uint_t&                        level,
              const DoFType&                       flag )
   {
      uOld_.copyBoundaryConditionFromFunction( u );
      fWeak_.copyBoundaryConditionFromFunction( u );

      if ( A.getTimeIntegrator() == DiffusionTimeIntegrator::ImplicitEuler )
      {
         M.apply( uOld, uOld_, level, flag );
         uOld_.assign( {1.0}, {uOld_}, level, flag );
         solver_->solve( A, u, uOld_, level );
      }
      else if ( A.getTimeIntegrator() == DiffusionTimeIntegrator::CrankNicolson )
      {
         // Crank-Nicholson
         M.apply( uOld, uOld_, level, flag );
         uOld_.assign( {1.0}, {uOld_}, level, flag );
         L.apply( uOld, fWeak_, level, flag );
         uOld_.assign( {1.0, -0.5 * A.dt()}, {uOld_, fWeak_}, level, flag );
         solver_->solve( A, u, uOld_, level );
      }
   }

   /// \brief Calculates the residual of the computed solution u^n+1 and writes it to r.
   void calculateResidual( const UnsteadyDiffusionOperatorType& A,
                           const LaplaceOperatorType&           L,
                           const MassOperatorType&              M,
                           const FunctionType&                  u,
                           const FunctionType&                  uOld,
                           const FunctionType&                  f,
                           const FunctionType&                  fOld,
                           FunctionType&                        r,
                           const uint_t&                        level,
                           const DoFType&                       flag )
   {
      uOld_.copyBoundaryConditionFromFunction( u );
      fWeak_.copyBoundaryConditionFromFunction( u );

      if ( A.getTimeIntegrator() == DiffusionTimeIntegrator::ImplicitEuler )
      {
         M.apply( f, fWeak_, level, flag );
         M.apply( uOld, uOld_, level, flag );
         uOld_.assign( {1.0, A.dt()}, {uOld_, fWeak_}, level, flag );
      }
      else if ( A.getTimeIntegrator() == DiffusionTimeIntegrator::CrankNicolson )
      {
         M.apply( f, fWeak_, level, flag );
         M.apply( fOld, fWeak_, level, flag, Add );
         M.apply( uOld, uOld_, level, flag );
         uOld_.assign( {1.0, 0.5 * A.dt()}, {uOld_, fWeak_}, level, flag );
         L.apply( uOld, fWeak_, level, flag );
         uOld_.assign( {1.0, -0.5 * A.dt()}, {uOld_, fWeak_}, level, flag );
      }
      A.apply( u, fWeak_, level, flag );
      r.assign( {1.0, -1.0}, {uOld_, fWeak_}, level, flag );
   }

 private:
   const std::shared_ptr< PrimitiveStorage >                  storage_;
   FunctionType                                               uOld_;
   FunctionType                                               fWeak_;
   std::shared_ptr< Solver< UnsteadyDiffusionOperatorType > > solver_;
};

} // namespace hyteg
