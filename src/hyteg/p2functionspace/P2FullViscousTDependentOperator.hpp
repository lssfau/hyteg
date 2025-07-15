/*
 * Copyright (c) 2023-2025 Andreas Burkhart.
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

#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "elementwise_dof_value_operator/generated/p1_full_stokes_T_p1_dep_eta_centroid_blending_q3_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_full_stokes_T_p2_dep_eta_blending_q4_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_full_stokes_T_p2_dep_eta_centroid_blending_q4_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_blending_q4_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_full_stokes_eta_p1_blending_q4_ElementwiseOperator.hpp"
#include "mixed_operator/VectorToVectorOperator.hpp"

namespace hyteg {

using walberla::real_t;

class P2ElementwiseBlendingFullViscousTDependentOperator
: public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
  public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >,
  public WeightedJacobiSmoothable< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseBlendingFullViscousTDependentOperator( const std::shared_ptr< PrimitiveStorage >&        storage,
                                                       size_t                                            minLevel,
                                                       size_t                                            maxLevel,
                                                       P2Function< real_t >&                             T,
                                                       std::function< real_t( const Point3D&, real_t ) > viscosity,
                                                       bool computeInverseDiagEntries                        = true,
                                                       std::shared_ptr< P2ProjectNormalOperator > projection = nullptr )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   , projection_( projection )
   {
      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_0_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_0_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][2] = std::make_shared< p2_full_stokes_T_p2_dep_eta_0_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[1][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_1_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_1_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][2] = std::make_shared< p2_full_stokes_T_p2_dep_eta_1_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[2][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_2_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[2][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_2_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[2][2] = std::make_shared< p2_full_stokes_T_p2_dep_eta_2_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_0_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_0_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[1][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_1_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_1_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
      }
   }

   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const override
   {
      VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::apply( src, dst, level, flag, updateType );
      if ( projection_ != nullptr )
      {
         projection_->project( dst, level, FreeslipBoundary );
      }
   };

   void smooth_jac( const P2VectorFunction< real_t >& dst,
                    const P2VectorFunction< real_t >& rhs,
                    const P2VectorFunction< real_t >& src,
                    real_t                            relax,
                    size_t                            level,
                    DoFType                           flag ) const override
   {
      for ( uint_t k = 0; k < this->dim_; ++k )
      {
         if ( const auto subOp =
                  std::dynamic_pointer_cast< WeightedJacobiSmoothable< P2Function< real_t > > >( this->subOper_[k][k] ) )
         {
            subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
         }
         else
         {
            throw std::runtime_error(
                "Jacobi smoothing of a P2ElementwiseBlendingFullViscousTDependentOperator requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
         }
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::computeInverseDiagonalOperatorValues();
   }

   void computeAndStoreLocalElementMatrices()
   {
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_0_0_blending_q4_ElementwiseOperator* >( this->subOper_[0][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_0_1_blending_q4_ElementwiseOperator* >( this->subOper_[0][1].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_1_0_blending_q4_ElementwiseOperator* >( this->subOper_[1][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_1_1_blending_q4_ElementwiseOperator* >( this->subOper_[1][1].get() ) )
          ->computeAndStoreLocalElementMatrices();

      if ( this->dim_ == 3 )
      {
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_0_2_blending_q4_ElementwiseOperator* >( this->subOper_[0][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_1_2_blending_q4_ElementwiseOperator* >( this->subOper_[1][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_2_0_blending_q4_ElementwiseOperator* >( this->subOper_[2][0].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_2_1_blending_q4_ElementwiseOperator* >( this->subOper_[2][1].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_2_2_blending_q4_ElementwiseOperator* >( this->subOper_[2][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
      }
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

class P2ElementwiseBlendingFullViscousTDependentOperator_Centroid
: public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
  public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >,
  public WeightedJacobiSmoothable< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseBlendingFullViscousTDependentOperator_Centroid(
       const std::shared_ptr< PrimitiveStorage >&                                        storage,
       size_t                                                                            minLevel,
       size_t                                                                            maxLevel,
       P2Function< real_t >&                                                             T,
       std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > viscosity,
       bool                                                                              computeInverseDiagEntries = true,
       std::shared_ptr< P2ProjectNormalOperator >                                        projection                = nullptr )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   , projection_( projection )
   {
      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_0_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_0_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][2] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_0_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[1][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_1_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_1_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][2] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_1_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[2][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_2_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[2][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_2_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[2][2] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_2_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_0_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_0_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[1][0] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_1_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][1] = std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_1_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
      }
   }

   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const override
   {
      VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::apply( src, dst, level, flag, updateType );
      if ( projection_ != nullptr )
      {
         projection_->project( dst, level, FreeslipBoundary );
      }
   };

   void smooth_jac( const P2VectorFunction< real_t >& dst,
                    const P2VectorFunction< real_t >& rhs,
                    const P2VectorFunction< real_t >& src,
                    real_t                            relax,
                    size_t                            level,
                    DoFType                           flag ) const override
   {
      for ( uint_t k = 0; k < this->dim_; ++k )
      {
         if ( const auto subOp =
                  std::dynamic_pointer_cast< WeightedJacobiSmoothable< P2Function< real_t > > >( this->subOper_[k][k] ) )
         {
            subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
         }
         else
         {
            throw std::runtime_error(
                "Jacobi smoothing of a P2ElementwiseBlendingFullViscousTDependentOperator_Centroid requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
         }
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::computeInverseDiagonalOperatorValues();
   }

   void computeAndStoreLocalElementMatrices()
   {
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_0_0_blending_q4_ElementwiseOperator* >( this->subOper_[0][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_0_1_blending_q4_ElementwiseOperator* >( this->subOper_[0][1].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_1_0_blending_q4_ElementwiseOperator* >( this->subOper_[1][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_1_1_blending_q4_ElementwiseOperator* >( this->subOper_[1][1].get() ) )
          ->computeAndStoreLocalElementMatrices();

      if ( this->dim_ == 3 )
      {
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_0_2_blending_q4_ElementwiseOperator* >(
               this->subOper_[0][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_1_2_blending_q4_ElementwiseOperator* >(
               this->subOper_[1][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_2_0_blending_q4_ElementwiseOperator* >(
               this->subOper_[2][0].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_2_1_blending_q4_ElementwiseOperator* >(
               this->subOper_[2][1].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_2_2_blending_q4_ElementwiseOperator* >(
               this->subOper_[2][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
      }
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

class P2ElementwiseBlendingFullViscousTDependentOperator_Centroid_ScaledDivDiv
: public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
  public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >,
  public WeightedJacobiSmoothable< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseBlendingFullViscousTDependentOperator_Centroid_ScaledDivDiv(
       const std::shared_ptr< PrimitiveStorage >&                                        storage,
       size_t                                                                            minLevel,
       size_t                                                                            maxLevel,
       P2Function< real_t >&                                                             T,
       std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > viscosity,
       std::function< real_t( const Point3D& ) >                                         divdivScaling,
       bool                                                                              computeInverseDiagEntries = true,
       std::shared_ptr< P2ProjectNormalOperator >                                        projection                = nullptr )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   , projection_( projection )
   {
      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_0_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[0][1] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_1_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[0][2] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_2_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );

         this->subOper_[1][0] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_0_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[1][1] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_1_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[1][2] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_2_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );

         this->subOper_[2][0] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_2_0_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[2][1] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_2_1_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[2][2] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_2_2_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
      }
      else
      {
         this->subOper_[0][0] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_0_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[0][1] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_1_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );

         this->subOper_[1][0] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_0_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
         this->subOper_[1][1] =
             std::make_shared< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_1_blending_q4_ElementwiseOperator >(
                 storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity, divdivScaling );
      }
   }

   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const override
   {
      if ( projection_ != nullptr )
      {
         projection_->project( src, level, FreeslipBoundary );
      }
      VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::apply( src, dst, level, flag, updateType );
      if ( projection_ != nullptr )
      {
         projection_->project( dst, level, FreeslipBoundary );
      }
   };

   void smooth_jac( const P2VectorFunction< real_t >& dst,
                    const P2VectorFunction< real_t >& rhs,
                    const P2VectorFunction< real_t >& src,
                    real_t                            relax,
                    size_t                            level,
                    DoFType                           flag ) const override
   {
      for ( uint_t k = 0; k < this->dim_; ++k )
      {
         if ( const auto subOp =
                  std::dynamic_pointer_cast< WeightedJacobiSmoothable< P2Function< real_t > > >( this->subOper_[k][k] ) )
         {
            subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
         }
         else
         {
            throw std::runtime_error(
                "Jacobi smoothing of a P2ElementwiseBlendingFullViscousTDependentOperator_centroid_scaled_divdiv requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
         }
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::computeInverseDiagonalOperatorValues();
   }

   void computeAndStoreLocalElementMatrices()
   {
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_0_blending_q4_ElementwiseOperator* >(
            this->subOper_[0][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_1_blending_q4_ElementwiseOperator* >(
            this->subOper_[0][1].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_0_blending_q4_ElementwiseOperator* >(
            this->subOper_[1][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_1_blending_q4_ElementwiseOperator* >(
            this->subOper_[1][1].get() ) )
          ->computeAndStoreLocalElementMatrices();

      if ( this->dim_ == 3 )
      {
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_0_2_blending_q4_ElementwiseOperator* >(
               this->subOper_[0][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_1_2_blending_q4_ElementwiseOperator* >(
               this->subOper_[1][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_2_0_blending_q4_ElementwiseOperator* >(
               this->subOper_[2][0].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_2_1_blending_q4_ElementwiseOperator* >(
               this->subOper_[2][1].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_T_p2_dep_eta_centroid_scaled_divdiv_2_2_blending_q4_ElementwiseOperator* >(
               this->subOper_[2][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
      }
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

class P2ElementwiseBlendingFullViscousFEMOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
                                                    public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >,
                                                    public WeightedJacobiSmoothable< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseBlendingFullViscousFEMOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                size_t                                     minLevel,
                                                size_t                                     maxLevel,
                                                P1Function< real_t >&                      eta,
                                                bool                                       computeInverseDiagEntries = true,
                                                std::shared_ptr< P2ProjectNormalOperator > projection                = nullptr )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   , projection_( projection )
   {
      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = std::make_shared< p2_full_stokes_eta_p1_0_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[0][1] = std::make_shared< p2_full_stokes_eta_p1_0_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[0][2] = std::make_shared< p2_full_stokes_eta_p1_0_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );

         this->subOper_[1][0] = std::make_shared< p2_full_stokes_eta_p1_1_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[1][1] = std::make_shared< p2_full_stokes_eta_p1_1_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[1][2] = std::make_shared< p2_full_stokes_eta_p1_1_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );

         this->subOper_[2][0] = std::make_shared< p2_full_stokes_eta_p1_2_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[2][1] = std::make_shared< p2_full_stokes_eta_p1_2_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[2][2] = std::make_shared< p2_full_stokes_eta_p1_2_2_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< p2_full_stokes_eta_p1_0_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[0][1] = std::make_shared< p2_full_stokes_eta_p1_0_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );

         this->subOper_[1][0] = std::make_shared< p2_full_stokes_eta_p1_1_0_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
         this->subOper_[1][1] = std::make_shared< p2_full_stokes_eta_p1_1_1_blending_q4_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, eta );
      }
   }

   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const override
   {
      if ( projection_ != nullptr )
      {
         projection_->project( src, level, FreeslipBoundary );
      }
      VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::apply( src, dst, level, flag, updateType );
      if ( projection_ != nullptr )
      {
         projection_->project( dst, level, FreeslipBoundary );
      }
   };

   void smooth_jac( const P2VectorFunction< real_t >& dst,
                    const P2VectorFunction< real_t >& rhs,
                    const P2VectorFunction< real_t >& src,
                    real_t                            relax,
                    size_t                            level,
                    DoFType                           flag ) const override
   {
      for ( uint_t k = 0; k < this->dim_; ++k )
      {
         if ( const auto subOp =
                  std::dynamic_pointer_cast< WeightedJacobiSmoothable< P2Function< real_t > > >( this->subOper_[k][k] ) )
         {
            subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
         }
         else
         {
            throw std::runtime_error(
                "Jacobi smoothing of a P2ElementwiseBlendingFullViscousTDependentOperator_centroid_scaled_divdiv requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
         }
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::computeInverseDiagonalOperatorValues();
   }

   void computeAndStoreLocalElementMatrices()
   {
      ( dynamic_cast< p2_full_stokes_eta_p1_0_0_blending_q4_ElementwiseOperator* >( this->subOper_[0][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_eta_p1_0_1_blending_q4_ElementwiseOperator* >( this->subOper_[0][1].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_eta_p1_1_0_blending_q4_ElementwiseOperator* >( this->subOper_[1][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p2_full_stokes_eta_p1_1_1_blending_q4_ElementwiseOperator* >( this->subOper_[1][1].get() ) )
          ->computeAndStoreLocalElementMatrices();

      if ( this->dim_ == 3 )
      {
         ( dynamic_cast< p2_full_stokes_eta_p1_0_2_blending_q4_ElementwiseOperator* >( this->subOper_[0][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_eta_p1_1_2_blending_q4_ElementwiseOperator* >( this->subOper_[1][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_eta_p1_2_0_blending_q4_ElementwiseOperator* >( this->subOper_[2][0].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_eta_p1_2_1_blending_q4_ElementwiseOperator* >( this->subOper_[2][1].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p2_full_stokes_eta_p1_2_2_blending_q4_ElementwiseOperator* >( this->subOper_[2][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
      }
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

class P1ElementwiseBlendingFullViscousTDependentOperator_Centroid
: public VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >,
  public OperatorWithInverseDiagonal< P1VectorFunction< real_t > >,
  public WeightedJacobiSmoothable< P1VectorFunction< real_t > >
{
 public:
   P1ElementwiseBlendingFullViscousTDependentOperator_Centroid(
       const std::shared_ptr< PrimitiveStorage >&                                        storage,
       size_t                                                                            minLevel,
       size_t                                                                            maxLevel,
       P1Function< real_t >&                                                             T,
       std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > viscosity,
       bool                                                                              computeInverseDiagEntries = true,
       std::shared_ptr< P1ProjectNormalOperator >                                        projection                = nullptr )
   : VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >( storage, minLevel, maxLevel )
   , projection_( projection )
   {
      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_0_0_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][1] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_0_1_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][2] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_0_2_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[1][0] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_1_0_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][1] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_1_1_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][2] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_1_2_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[2][0] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_2_0_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[2][1] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_2_1_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[2][2] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_2_2_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_0_0_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[0][1] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_0_1_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );

         this->subOper_[1][0] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_1_0_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
         this->subOper_[1][1] = std::make_shared< p1_full_stokes_T_p1_dep_eta_centroid_1_1_blending_q3_ElementwiseOperator >(
             storage, minLevel, maxLevel, computeInverseDiagEntries, T, viscosity );
      }
   }

   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const override
   {
      VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >::apply( src, dst, level, flag, updateType );
      if ( projection_ != nullptr )
      {
         projection_->project( dst, level, FreeslipBoundary );
      }
   };

   void smooth_jac( const P1VectorFunction< real_t >& dst,
                    const P1VectorFunction< real_t >& rhs,
                    const P1VectorFunction< real_t >& src,
                    real_t                            relax,
                    size_t                            level,
                    DoFType                           flag ) const override
   {
      for ( uint_t k = 0; k < this->dim_; ++k )
      {
         if ( const auto subOp =
                  std::dynamic_pointer_cast< WeightedJacobiSmoothable< P1Function< real_t > > >( this->subOper_[k][k] ) )
         {
            subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
         }
         else
         {
            throw std::runtime_error(
                "Jacobi smoothing of a P2ElementwiseBlendingFullViscousTDependentOperator_Centroid requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
         }
      }
   }

   std::shared_ptr< P1VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >::computeInverseDiagonalOperatorValues();
   }

   void computeAndStoreLocalElementMatrices()
   {
      ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_0_0_blending_q3_ElementwiseOperator* >( this->subOper_[0][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_0_1_blending_q3_ElementwiseOperator* >( this->subOper_[0][1].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_1_0_blending_q3_ElementwiseOperator* >( this->subOper_[1][0].get() ) )
          ->computeAndStoreLocalElementMatrices();
      ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_1_1_blending_q3_ElementwiseOperator* >( this->subOper_[1][1].get() ) )
          ->computeAndStoreLocalElementMatrices();

      if ( this->dim_ == 3 )
      {
         ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_0_2_blending_q3_ElementwiseOperator* >(
               this->subOper_[0][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_1_2_blending_q3_ElementwiseOperator* >(
               this->subOper_[1][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_2_0_blending_q3_ElementwiseOperator* >(
               this->subOper_[2][0].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_2_1_blending_q3_ElementwiseOperator* >(
               this->subOper_[2][1].get() ) )
             ->computeAndStoreLocalElementMatrices();
         ( dynamic_cast< p1_full_stokes_T_p1_dep_eta_centroid_2_2_blending_q3_ElementwiseOperator* >(
               this->subOper_[2][2].get() ) )
             ->computeAndStoreLocalElementMatrices();
      }
   }

 private:
   std::shared_ptr< P1ProjectNormalOperator > projection_;
};

} // namespace hyteg
