/*
* Copyright (c) 2017-2024 Nils Kohl, Marcus Mohr.
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

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using namespace dg;

template < typename ValueType >
class DG1Function : public Function< DG1Function< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = DG1Function< VType >;

   DG1Function( const std::string&                         name,
                const std::shared_ptr< PrimitiveStorage >& storage,
                uint_t                                     minLevel,
                uint_t                                     maxLevel,
                BoundaryCondition                          bc )
   : Function< DG1Function< ValueType > >( name, storage, minLevel, maxLevel )
   {
      auto basis  = std::make_shared< DGBasisLinearLagrange_Example >();
      dgFunction_ = std::make_shared< DGFunction< ValueType > >( name, storage, minLevel, maxLevel, basis, 1, bc );
   }

   DG1Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : DG1Function( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   const std::shared_ptr< DGFunction< ValueType > > getDGFunction() const { return dgFunction_; }

   void setBoundaryCondition( BoundaryCondition bc ) { dgFunction_->setBoundaryCondition( bc ); }

   BoundaryCondition getBoundaryCondition() const { return dgFunction_->getBoundaryCondition(); }

   uint_t getDimension() const override final { return dgFunction_->getDimension(); };

   void communicate( const uint_t& level ) const { dgFunction_->communicate( level ); }

   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const { dgFunction_->add( scalar, level, flag ); };

   void add( const std::vector< ValueType >                                                 scalars,
             const std::vector< std::reference_wrapper< const DG1Function< ValueType > > >& functions,
             uint_t                                                                         level,
             DoFType                                                                        flag = All ) const
   {
      std::vector< ValueType > new_scalars( scalars );
      new_scalars.push_back( 1.0 );
      std::vector< std::reference_wrapper< const DG1Function< ValueType > > > new_functions( functions );
      new_functions.push_back( *this );
      assign( new_scalars, new_functions, level, flag );
   };

   void multElementwise( const std::vector< std::reference_wrapper< const DG1Function< ValueType > > >& functions,
                         uint_t                                                                         level,
                         DoFType                                                                        flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void interpolate( ValueType constant, uint_t level, DoFType dofType = All ) const { WALBERLA_ABORT( "Not implemented." ); }

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType dofType = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const override final { dgFunction_->setToZero( level ); };

   void swap( const DG1Function< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   void copyFrom( const DG1Function< ValueType >&        other,
                  const uint_t&                          level,
                  const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   bool evaluate( const Point3D& coordinates, uint_t level, ValueType& value, real_t searchToleranceRadius = 1e-05 ) const
   {
      return dgFunction_->evaluate( coordinates, level, value, searchToleranceRadius );
   }

   void assign( const std::vector< ValueType >&                                                scalars,
                const std::vector< std::reference_wrapper< const DG1Function< ValueType > > >& functions,
                uint_t                                                                         level,
                DoFType                                                                        flag = All ) const
   {
      std::vector< std::reference_wrapper< const DGFunction< ValueType > > > dgFunctions;
      for ( auto f : functions )
      {
         dgFunctions.push_back( *f.get().getDGFunction() );
      }
      dgFunction_->assign( scalars, dgFunctions, level, flag );
   }

   ValueType dotGlobal( const DG1Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      return dgFunction_->dotGlobal( *rhs.getDGFunction(), level );
   }

   ValueType dotLocal( const DG1Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      return dgFunction_->dotLocal( *rhs.getDGFunction(), level );
   }

   ValueType sumGlobal( uint_t level, const DoFType& flag = All ) const { return dgFunction_->sumGlobal( level ); }

   void enumerate( uint_t level ) const { dgFunction_->enumerate( level ); }

   void enumerate( uint_t level, ValueType& offset ) const { dgFunction_->enumerate( level, offset ); }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const DG1Function< idx_t >&           numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      dgFunction_->toVector( *numerator.getDGFunction(), vec, level, flag );
   }

   void fromVector( const DG1Function< idx_t >&           numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      dgFunction_->fromVector( *numerator.getDGFunction(), vec, level, flag );
   }
   /// @}

   uint_t getNumberOfLocalDoFs( uint_t level ) const { return dgFunction_->getNumberOfLocalDoFs( level ); }

   uint_t getNumberOfGlobalDoFs( uint_t          level,
                                 const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                 const bool&     onRootOnly   = false ) const
   {
      return dgFunction_->getNumberOfGlobalDoFs( level, communicator, onRootOnly );
   }

   ValueType getMaxDoFMagnitude( uint_t level, bool mpiReduce = true ) const
   {
      return dgFunction_->getMaxDoFMagnitude( level, mpiReduce );
   }

   ValueType getMaxDoFValue( uint_t level, bool mpiReduce = true ) const
   {
      return dgFunction_->getMaxDoFValue( level, mpiReduce );
   }

   ValueType getMinDoFValue( uint_t level, bool mpiReduce = true ) const
   {
      return dgFunction_->getMinDoFValue( level, mpiReduce );
   }

   template < typename OtherValueType >
   void copyBoundaryConditionFromFunction( const DG1Function< OtherValueType >& other )
   {
      dgFunction_->copyBoundaryConditionFromFunction( *other.getDGFunction() );
   }

   /// Evaluates the linear functional
   /// \f[
   ///   l( v ) = \int_\Omega f \cdot v
   /// \f]
   /// by integration over the local basis functions and writes the result into the vector, i.e.
   /// \f[
   ///   u_i \leftarrow \int_T f \cdot \phi_i
   /// \f]
   /// where \f$\phi_i\f$ is the basis function associated with the DoF \f$u_i\f$ and \f$f\f$ a given analytical function.
   void evaluateLinearFunctional( const std::function< real_t( const Point3D& ) >& f, uint_t level )
   {
      dgFunction_->evaluateLinearFunctional( f, level );
   }

   /// \brief Applies the Dirichlet boundary conditions to this function, treating it as the right-hand side of the linear system
   ///        that corresponds to the form object.
   void applyDirichletBoundaryConditions( const std::shared_ptr< DGForm >& form, uint_t level )
   {
      dgFunction_->applyDirichletBoundaryConditions( form, level );
   }

 private:
   std::shared_ptr< DGFunction< ValueType > > dgFunction_;
};

void applyDirichletBC( const DG1Function< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   WALBERLA_UNUSED( numerator );
   WALBERLA_UNUSED( mat );
   WALBERLA_UNUSED( level );
}
} // namespace hyteg
