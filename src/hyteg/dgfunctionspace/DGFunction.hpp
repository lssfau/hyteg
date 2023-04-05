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

#pragma once

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGForm.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"
namespace hyteg {
namespace dg {

using volumedofspace::VolumeDoFFunction;
using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType >
class DGFunction final : public Function< DGFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = DGFunction< VType >;

   /// \brief Instantiates a DGFunction.
   ///
   /// \param name              name of the function for output / debugging purposes
   /// \param storage           PrimitiveStorage that this lives on
   /// \param minLevel          min refinement level to allocate
   /// \param maxLevel          max refinement level to allocate
   /// \param basis             polynomial basis information through corresponding DGBasisInfo object
   /// \param initialPolyDegree polynomial degree to initialize with globally
   DGFunction( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               const std::shared_ptr< DGBasisInfo >&      basis,
               uint_t                                     initialPolyDegree,
               BoundaryCondition                          boundaryCondition = BoundaryCondition::create0123BC() );


   void multElementwise( const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                         uint_t                                                                        level,
                         DoFType                                                                       flag = All ) const
   {
      WALBERLA_UNUSED( functions );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "Not implemented." );
   }

   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const
   {
      WALBERLA_UNUSED( constant );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "Not implemented." );
   };

   void interpolate( const std::function< ValueType( const hyteg::Point3D& ) >& expr, uint_t level, DoFType flag = All ) const
   {
      WALBERLA_UNUSED( expr );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "Not implemented." );
   };

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_UNUSED( expressions );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "Not implemented." );
   };

   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const { volumeDoFFunction_->add( scalar, level, flag ); };

   void add( const std::vector< ValueType >                                                scalars,
             const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
             uint_t                                                                        level,
             DoFType                                                                       flag = All ) const
   {
      WALBERLA_UNUSED( flag );
      std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > > vFunctions;
      for ( const auto& f : functions )
      {
         vFunctions.push_back( *( f.get().volumeDoFFunction_ ) );
      }
      volumeDoFFunction_->add( scalars, vFunctions, level );
   };

   /// \brief Assigns a linear combination of multiple VolumeDoFFunctions to this.
   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All )
   {
      WALBERLA_UNUSED( flag );
      std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > > vFunctions;
      for ( const auto& f : functions )
      {
         vFunctions.push_back( *( f.get().volumeDoFFunction_ ) );
      }

      volumeDoFFunction_->assign( scalars, vFunctions, level );
   }

   /// \brief Evaluates the dot product on all local DoFs. No communication is involved and the results may be different on each
   /// process.
   ValueType dotLocal( const DGFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const
   {
      WALBERLA_UNUSED( flag );
      return volumeDoFFunction_->dotLocal( *rhs.volumeDoFFunction_, level );
   }

   /// \brief Evaluates the (global) dot product. Involves communication and has to be called collectively.
   ValueType dotGlobal( const DGFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const
   {
      WALBERLA_UNUSED( flag );
      return volumeDoFFunction_->dotGlobal( *rhs.volumeDoFFunction_, level );
   }

   /// \brief Evaluates the sum on all DoFs
   ValueType sumGlobal( uint_t level, DoFType flag = All ) const
   {
      WALBERLA_UNUSED( flag );
      return volumeDoFFunction_->sumGlobal( level );
   }

   /// \brief Evaluate finite element function at a specific coordinates.
   ///
   /// In a parallel setting, the specified coordinate might not lie in the local subdomain.
   ///
   /// Evaluation is performed in two steps:
   ///
   ///   1. For all volume primitives of the local subdomain:
   ///      If a point-tet (point-triangle in 2D) inclusion test succeeds, the function returns true
   ///      and the finite-element function is evaluated.
   ///
   ///   2. Skipped, if radius is negative.
   ///      For all volume primitives of the local subdomain:
   ///      A sphere-tet (circle-triangle) intersection
   ///      test is performed, if successful returns true, and the finite-element function is extrapolated
   ///      to the specified coordinate and evaluated, depending on the radius, this might introduce (large) errors.
   ///
   /// If both tests fail, this function returns false, and no evaluation is performed (i.e. the returned, evaluated
   /// value is not set to anything meaningful).
   ///
   /// Note that two parallel processes that return true, may return _different_ values.
   ///
   /// No communication is performed in this function.
   /// -> Does not need to be called collectively.
   /// -> Different values are returned on each process.
   ///
   /// \param physicalCoords coordinates in physical domain where the function is to be evaluated
   /// \param level refinement level
   /// \param value function value at the coordinate if search was successful
   /// \param searchToleranceRadius radius of the sphere (circle) for the second search phase, skipped if negative
   /// \return true if the function was evaluated successfully, false otherwise
   ///
   bool evaluate( const Point3D& physicalCoords, uint_t level, ValueType& value, real_t searchToleranceRadius = 1e-05 ) const;

   /// \brief Evaluate finite element function on a specific micro-face.
   ///
   /// The standard evaluate() function does not take the discontinuity of the function into account.
   /// If the function is evaluated at a discontinuity with evaluate(), it is kind of "random" which of the neighboring elements
   /// is used to evaluate the polynom.
   ///
   /// This is for example an issue for accurate visualization of discontinuities.
   ///
   /// The present method solves this issue by taking a specific (micro-)element as input argument. This way the function can be
   /// evaluated on _different_ elements at the same coordinate. Note that regardless of whether the specified coordinate is
   /// located on the element or not - the local polynomial is extrapolated.
   ///
   /// \param coordinates  where the function shall be evaluated
   /// \param level        refinement level
   /// \param faceID       the macro-face where the (micro-)element is located on
   /// \param elementIndex the logical index of the micro-element
   /// \param faceType     the type of the local micro-element
   /// \param value        the evaluation
   void evaluateOnMicroElement( const Point3D&         coordinates,
                                uint_t                 level,
                                const PrimitiveID&     faceID,
                                hyteg::indexing::Index elementIndex,
                                facedof::FaceType      faceType,
                                ValueType&             value ) const;

   /// \brief Evaluate finite element function on a specific micro-cell.
   ///
   /// The standard evaluate() function does not take the discontinuity of the function into account.
   /// If the function is evaluated at a discontinuity with evaluate(), it is kind of "random" which of the neighboring elements
   /// is used to evaluate the polynom.
   ///
   /// This is for example an issue for accurate visualization of discontinuities.
   ///
   /// The present method solves this issue by taking a specific (micro-)element as input argument. This way the function can be
   /// evaluated on _different_ elements at the same coordinate. Note that regardless of whether the specified coordinate is
   /// located on the element or not - the local polynomial is extrapolated.
   ///
   /// \param coordinates  where the function shall be evaluated
   /// \param level        refinement level
   /// \param cellID       the macro-cell where the (micro-)element is located on
   /// \param elementIndex the logical index of the micro-element
   /// \param cellType     the type of the local micro-element
   /// \param value        the evaluation
   void evaluateOnMicroElement( const Point3D&         coordinates,
                                uint_t                 level,
                                const PrimitiveID&     cellID,
                                hyteg::indexing::Index elementIndex,
                                celldof::CellType      cellType,
                                ValueType&             value ) const;

   /// Evaluates the linear functional
   ///
   ///   l( v ) = \int_\Omega f * v
   ///
   /// by integration over the local basis functions and writes the result into the vector, i.e.
   ///
   ///   u_i <- \int_T f * \phi_i
   ///
   /// where \phi_i is the basis function associated with the DoF u_i and f a given analytical function.
   void evaluateLinearFunctional( const std::function< real_t( const Point3D& ) >& f, uint_t level );

   /// \brief Applies the Dirichlet boundary conditions to this function, treating it as the right-hand side of the linear system
   ///        that corresponds to the form object.
   void applyDirichletBoundaryConditions( const std::shared_ptr< DGForm >& form, uint_t level );

   /// \brief Returns the internally stored VolumeDoFFunction.
   [[nodiscard]] std::shared_ptr< volumedofspace::VolumeDoFFunction< ValueType > > volumeDoFFunction() const
   {
      return volumeDoFFunction_;
   }

   /// \brief Returns the associated DGBasisInfo object.
   [[nodiscard]] std::shared_ptr< DGBasisInfo > basis() const { return basis_; }

   /// \brief Returns the polynomial degree of the function on the passed primitive.
   [[nodiscard]] int polynomialDegree( const PrimitiveID& primitiveID ) const
   {
      WALBERLA_ASSERT( storage_->primitiveExistsLocally( primitiveID ),
                       "No information on the polynomial degree on this primitive (primitive does not exist locally)." );
      WALBERLA_ASSERT( polyDegreesPerPrimitive_.count( primitiveID ) > 0,
                       "No information on the polynomial degree on this primitive (no information stored in the map)." );
      return polyDegreesPerPrimitive_.at( primitiveID );
   }

   /// \brief Assigns unique values to all data points. Increments by 1 per DoF.
   ///
   /// Mainly used for sparse matrix assembly to get global integer identifiers for all DoFs.
   ///
   /// \param level refinement level to enumerate
   void enumerate( uint_t level ) const;

   /// \brief Like enumerate() but starting with a value given by the offset parameter.
   ///
   /// \param level  refinement level to enumerate
   /// \param offset value to start from, the next "free" value is returned
   void enumerate( uint_t level, ValueType& offset ) const;

   template < typename OtherValueType >
   void copyBoundaryConditionFromFunction( const DGFunction< OtherValueType >& other )
   {
      boundaryCondition_ = other.getBoundaryCondition();
   }

   /// \brief Returns the number of DoFs that are allocated on this process.
   uint_t getNumberOfLocalDoFs( uint_t level ) const;

   /// \brief Returns the number of DoFs. Performs global reduction, must be called collectively.
   ///
   /// \param level        refinement level
   /// \param communicator if required, a custom communicator can be passed
   /// \param onRootOnly   if true, the result is only returned on the root process
   /// \return
   uint_t getNumberOfGlobalDoFs( uint_t          level,
                                 const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                 const bool&     onRootOnly   = false ) const;

   uint_t getDimension() const { return storage_->hasGlobalCells() ? 3 : 2; };
   /// \brief Updates ghost-layers.
   void communicate( uint_t level ) const { volumeDoFFunction_->communicate( level ); }

   /// \brief Returns the boundary conditions of this function.
   BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

   /// \brief Sets the boundary conditions of this function.
   void setBoundaryCondition( BoundaryCondition bc ) { boundaryCondition_ = bc; }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const DGFunction< idx_t >&            numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const;

   void fromVector( const DGFunction< idx_t >&            numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const;
   /// @}

   void copyFrom( const DGFunction< ValueType >&                 other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const
   {
      WALBERLA_UNUSED( other );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( localPrimitiveIDsToRank );
      WALBERLA_UNUSED( otherPrimitiveIDsToRank );
      WALBERLA_ABORT( "DGFunction::copyFrom() not implemented." );
   }

   void swap( const DGFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      WALBERLA_UNUSED( flag );
      volumeDoFFunction_->swap(*(other.volumeDoFFunction()), level);
   }

   /// \brief Returns the max absolute DoF.
   ///
   /// \param level     refinement level
   /// \param mpiReduce if true, reduces over all processes (global max magnitude), if false returns the process local value
   ValueType getMaxMagnitude( uint_t level, bool mpiReduce = true ) const;

 private:
   using Function< DGFunction< ValueType > >::communicators_;
   using Function< DGFunction< ValueType > >::additiveCommunicators_;

   std::string                         name_;
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;
   std::shared_ptr< DGBasisInfo >      basis_;

   std::shared_ptr< volumedofspace::VolumeDoFFunction< ValueType > > volumeDoFFunction_;

   std::map< PrimitiveID, uint_t > polyDegreesPerPrimitive_;

   BoundaryCondition boundaryCondition_;
};

} // namespace dg

void createVectorFromFunction( const dg::DGFunction< real_t >&       function,
                               const dg::DGFunction< idx_t >&        numerator,
                               const std::shared_ptr< VectorProxy >& vec,
                               uint_t                                level,
                               DoFType                               flag );

void createFunctionFromVector( const dg::DGFunction< real_t >&       function,
                               const dg::DGFunction< idx_t >&        numerator,
                               const std::shared_ptr< VectorProxy >& vec,
                               uint_t                                level,
                               DoFType                               flag );

void applyDirichletBC( const dg::DGFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level );

} // namespace hyteg
