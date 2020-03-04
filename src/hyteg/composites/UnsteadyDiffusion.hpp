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
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

/// \brief Unsteady diffusion operator for a P2 finite element discretization of the solution.
///
/// This is a composite operator used to solve the parabolic PDE that describes unsteady diffusion:
///
/// \f$ \frac{\partial u}{\partial t} = D  \Delta u + f\f$
///
/// where \f$ D \f$ is a constant diffusion parameter.
///
/// This operator can be used for an implicit Euler discretization in time.
/// In particular it equals
///
/// \f$ M + dt\ D\ \mathcal{L} \f$
///
/// where \f$ D \f$ is again the diffusivity constant, \f$ M \f$ the finite element mass matrix, \f$ dt \f$ a constant time step,
/// and \f$ \mathcal{L} = - \Delta \f$ the negative Laplacian.
///
/// To solve the unsteady diffusion equation, see UnsteadyDiffusion.
class P2UnsteadyDiffusionOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2UnsteadyDiffusionOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                const uint_t&                              minLevel,
                                const uint_t&                              maxLevel,
                                const real_t&                              dt,
                                const real_t&                              diffusionCoefficient )
   : Operator( storage, minLevel, maxLevel )
   , laplaceForm_( new LaplaceForm_T )
   , massForm_( new MassForm_T )
   , unsteadyDiffusionOperator_(
         storage,
         minLevel,
         maxLevel,
         P2LinearCombinationForm( {1.0, dt * diffusionCoefficient}, {massForm_.get(), laplaceForm_.get()} ) )
   , dt_( dt )
   {}

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               const uint_t&               level,
               const DoFType&              flag,
               const UpdateType&           updateType = Replace ) const
   {
      unsteadyDiffusionOperator_.apply( src, dst, level, flag, updateType );
   }

   void smooth_sor( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const real_t&               relax,
                    size_t                      level,
                    DoFType                     flag,
                    const bool&                 backwards = false ) const
   {
      unsteadyDiffusionOperator_.smooth_sor( dst, rhs, relax, level, flag, backwards );
   }

   void smooth_gs( const P2Function< real_t >& dst,
                   const P2Function< real_t >& rhs,
                   size_t                      level,
                   DoFType                     flag,
                   const bool&                 backwards = false ) const
   {
      smooth_sor( dst, rhs, 1.0, level, flag, backwards );
   }

   real_t dt() const { return dt_; }

 private:
   typedef P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > LaplaceForm_T;
   typedef P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise >           MassForm_T;

   std::shared_ptr< LaplaceForm_T >              laplaceForm_;
   std::shared_ptr< MassForm_T >                 massForm_;
   P2ConstantOperator< P2LinearCombinationForm > unsteadyDiffusionOperator_;

   real_t dt_;
};

/// \brief Wrapper class to solve the unsteady diffusion equation in time.
///
/// \tparam FunctionType function discretization type (e.g. P2Function< real_t >)
/// \tparam UnsteadyDiffusionOperatorType type of the unsteady diffusion operator, must match the function type (e.g. P2UnsteadyDiffusionOperator)
/// \tparam MassOperatorType mass matrix type matching the function (e.g. P2ConstantMassOperator)
///
/// Performs the implicit Euler method to advance the solution of the unsteady diffusion equation.
/// Therefore in each time step, the passed UnsteadyDiffusionOperatorType instance must be inverted.
///
template < typename FunctionType, typename UnsteadyDiffusionOperatorType, typename MassOperatorType >
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

   /// \brief Performs one implicit Euler step to advance the solution of the PDE.
   ///
   /// \param A the unsteady diffusion operator instance (dt and the diffusivity are defined by this operator)
   /// \param M the finite element mass matrix
   /// \param u the solution function
   /// \param f the right hand side (in "strong" form - it is multiplied by the mass in this method)
   /// \param level the refinement level
   /// \param flag where to solve
   ///
   /// This method performs
   ///
   /// \f$ u^{n+1} := (M + dt\ D\ \mathcal{L})^{-1} M (dt\ f + u^n) \f$
   ///
   /// where \f$ A = M + dt\ D\ \mathcal{L} \f$ the unsteady diffusion operator (see e.g. P2UnsteadyDiffusionOperator).
   ///
   void step( const UnsteadyDiffusionOperatorType& A,
              const MassOperatorType&              M,
              const FunctionType&                  u,
              const FunctionType&                  f,
              const uint_t&                        level,
              const DoFType&                       flag )
   {
      M.apply( f, fWeak_, level, flag );
      M.apply( u, uOld_, level, flag );
      uOld_.assign( {1.0, A.dt()}, {uOld_, fWeak_}, level, flag );
      solver_->solve( A, u, uOld_, level );
   }

 private:
   const std::shared_ptr< PrimitiveStorage >                  storage_;
   FunctionType                                               uOld_;
   FunctionType                                               fWeak_;
   std::shared_ptr< Solver< UnsteadyDiffusionOperatorType > > solver_;
};

} // namespace hyteg
