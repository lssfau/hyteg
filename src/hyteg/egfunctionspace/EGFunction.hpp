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

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/egfunctionspace/EGBasis.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

template < typename ValueType >
class EGFunction final : public Function< EGFunction< ValueType > >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = EGFunction< VType >;

   EGFunction( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          boundaryCondition = BoundaryCondition::create0123BC() );

   std::shared_ptr< P1VectorFunction< ValueType > > getConformingPart() const { return u_conforming_; }

   std::shared_ptr< P0Function< ValueType > > getDiscontinuousPart() const { return u_discontinuous_; }

   void setBoundaryCondition( BoundaryCondition bc )
   {
      u_conforming_->setBoundaryCondition( bc );
      u_discontinuous_->setBoundaryCondition( bc );
   }

   [[nodiscard]] BoundaryCondition getBoundaryCondition() const { return u_conforming_->getBoundaryCondition(); }

   template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const
   {
      u_conforming_->add( scalar, level, flag );
      u_discontinuous_->add( scalar, level, flag );
   };

   void add( const std::vector< ValueType >                                                scalars,
             const std::vector< std::reference_wrapper< const EGFunction< ValueType > > >& functions,
             uint_t                                                                        level,
             DoFType                                                                       flag = All ) const
   {
      u_conforming_->add( scalars, filter_conforming( functions ), level, flag );
      u_discontinuous_->add( scalars, filter_discontinuous( functions ), level, flag );
   };

   void multElementwise( const std::vector< std::reference_wrapper< const EGFunction< ValueType > > >& functions,
                         uint_t                                                                        level,
                         DoFType                                                                       flag = All ) const
   {
      u_conforming_->multElementwise( filter_conforming( functions ), level, flag );
      u_discontinuous_->multElementwise( filter_discontinuous( functions ), level, flag );
   }

   void interpolate( ValueType constant, uint_t level, DoFType dofType = All ) const
   {
      u_conforming_->interpolate( constant, level, dofType );
      u_discontinuous_->interpolate( 0, level, dofType );
   }

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType dofType = All ) const
   {
      u_conforming_->interpolate( expr, level, dofType );
      u_discontinuous_->interpolate( 0, level, dofType );
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      u_conforming_->interpolate( expressions, level, flag );
      u_discontinuous_->interpolate( 0, level, flag );
   };

   void swap( const EGFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      u_conforming_->swap( *other.getConformingPart(), level, flag );
      u_discontinuous_->swap( *other.getDiscontinuousPart(), level, flag );
   };

   void copyFrom( const EGFunction< ValueType >&                 other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      u_conforming_->copyFrom( *other.getConformingPart(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      u_discontinuous_->copyFrom( *other.getDiscontinuousPart(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   };

   bool evaluate( const Point3D& coordinates, uint_t level, ValueType& value, real_t searchToleranceRadius = 1e-05 ) const
   {
      WALBERLA_UNUSED( coordinates );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( value );
      WALBERLA_UNUSED( searchToleranceRadius );
      WALBERLA_ABORT( "Not implemented." );
   }

   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const EGFunction< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All ) const
   {
      u_conforming_->assign( scalars, filter_conforming( functions ), level, flag );
      u_discontinuous_->assign( scalars, filter_discontinuous( functions ), level, flag );
   }

   ValueType dotGlobal( const EGFunction< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      real_t result = 0;
      result += u_conforming_->dotGlobal( *rhs.getConformingPart(), level, flag );
      result += u_discontinuous_->dotGlobal( *rhs.getDiscontinuousPart(), level, flag );
      return result;
   }

   ValueType dotLocal( const EGFunction< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      real_t result = 0;
      result += u_conforming_->dotLocal( *rhs.getConformingPart(), level, flag );
      result += u_discontinuous_->dotLocal( *rhs.getDiscontinuousPart(), level, flag );
      return result;
   }

   void enumerate( uint_t level ) const
   {
      ValueType start = 0;
      enumerate( level, start );
      WALBERLA_UNUSED( start );
   }

   void enumerate( uint_t level, ValueType& offset ) const
   {
      u_conforming_->enumerate( level );
      offset += u_conforming_->getNumberOfGlobalDoFs( level );
      u_discontinuous_->enumerate( level, offset );
   }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const EGFunction< idx_t >&            numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      u_conforming_->toVector( *numerator.getConformingPart(), vec, level, flag );
      u_discontinuous_->toVector( *numerator.getDiscontinuousPart(), vec, level, flag );
   }

   void fromVector( const EGFunction< idx_t >&            numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      u_conforming_->fromVector( *numerator.getConformingPart(), vec, level, flag );
      u_discontinuous_->fromVector( *numerator.getDiscontinuousPart(), vec, level, flag );
   }
   /// @}

   [[nodiscard]] uint_t getNumberOfLocalDoFs( uint_t level ) const
   {
      return u_conforming_->getNumberOfLocalDoFs( level ) + u_discontinuous_->getNumberOfLocalDoFs( level );
   }

   [[nodiscard]] uint_t getNumberOfGlobalDoFs( uint_t          level,
                                               const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                               const bool&     onRootOnly   = false ) const
   {
      return u_conforming_->getNumberOfGlobalDoFs( level ) + u_discontinuous_->getNumberOfGlobalDoFs( level );
   }

   /// \brief Returns the max absolute DoF.
   ValueType getMaxMagnitude( uint_t level, bool mpiReduce = true ) const
   {
      return std::max( u_discontinuous_->getMaxMagnitude( level, mpiReduce ),
                       u_conforming_->getMaxComponentMagnitude( level, All, mpiReduce ) );
   }

   void evaluateLinearFunctional( const std::function< real_t( const Point3D& ) >& f0,
                                  const std::function< real_t( const Point3D& ) >& f1,
                                  uint_t                                           level )
   {
      if ( u_conforming_->getStorage()->hasGlobalCells() )
         WALBERLA_ABORT( "evaluation of linear functional not supported in 3D" );

      std::vector< std::function< real_t( const Point3D& ) > > f = { f0, f1 };

      for ( auto& it : this->getStorage()->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         // P1
         {
            const auto degree  = 1;
            const auto numDofs = 3;

            std::vector< uint_t > vertexDoFIndices( numDofs );
            std::vector< real_t > dofValues( numDofs );

            for ( int d = 0; d < 2; d += 1 )
            {
               auto dofs = u_conforming_->getStorage()
                               ->getFace( faceID )
                               ->template getData( u_conforming_->component( d ).getFaceDataID() )
                               ->getPointer( level );

               for ( auto faceType : facedof::allFaceTypes )
               {
                  for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
                  {
                     const std::array< indexing::Index, 3 > vertexIndices =
                         facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
                     std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > elementVertices;
                     for ( uint_t i = 0; i < 3; i++ )
                     {
                        const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                        elementVertices[i]( 0 )  = elementVertex[0];
                        elementVertices[i]( 1 )  = elementVertex[1];
                     }

                     vertexdof::getVertexDoFDataIndicesFromMicroFace( idxIt, faceType, level, vertexDoFIndices );
                     basis_conforming_.integrateBasisFunction( degree, elementVertices, f[d], dofValues );
                     for ( uint_t i = 0; i < numDofs; i++ )
                     {
                        dofs[vertexDoFIndices[i]] = ValueType( dofValues[i] );
                     }
                  }
               }
            }
         }

         // DGE
         {
            const auto degree  = 0;
            const auto numDofs = 1;

            std::vector< uint_t > vertexDoFIndices( numDofs );
            std::vector< real_t > dofValues( numDofs );

            auto       dofs      = u_discontinuous_->getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );
            const auto memLayout = u_discontinuous_->getDGFunction()->volumeDoFFunction()->memoryLayout();

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  const std::array< indexing::Index, 3 > vertexIndices =
                      facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
                  std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > elementVertices;
                  for ( uint_t i = 0; i < 3; i++ )
                  {
                     const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                  }

                  basis_discontinuous_.integrateBasisFunction( degree, elementVertices, f[0], f[1], dofValues );
                  for ( uint_t i = 0; i < numDofs; i++ )
                  {
                     dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )] =
                         ValueType( dofValues[i] );
                  }
               }
            }
         }
      }
   };

   uint_t getDimension() const override
   {
      if ( u_conforming_->getStorage()->hasGlobalCells() )
         return 3;
      else
         return 2;
   };

   void evaluateOnMicroElement( const Point3D&         coordinates,
                                uint_t                 level,
                                const PrimitiveID&     faceID,
                                hyteg::indexing::Index elementIndex,
                                facedof::FaceType      faceType,
                                uint_t                 componentIdx,
                                ValueType&             value ) const
   {
      // 2D

      auto storage = u_conforming_->getStorage();

      WALBERLA_ASSERT( !storage->hasGlobalCells() );

      Point2D coordinates2D( { coordinates[0], coordinates[1] } );

      WALBERLA_ASSERT( storage->faceExistsLocally( faceID ) );
      const Face& face = *storage->getFace( faceID );

      const uint_t polyDegree = 1;
      const uint_t ndofsP1    = 3;

      Eigen::Matrix< real_t, 2, 1 > affineCoordinates( { coordinates[0], coordinates[1] } );

      std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > affineElementVertices;
      auto vertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( elementIndex, faceType );
      for ( uint_t i = 0; i < 3; i++ )
      {
         const auto coord              = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
         affineElementVertices[i]( 0 ) = coord[0];
         affineElementVertices[i]( 1 ) = coord[1];
      }

      // trafo from affine to reference space
      Eigen::Matrix< real_t, 2, 2 > A;
      A( 0, 0 )       = ( affineElementVertices[1] - affineElementVertices[0] )( 0 );
      A( 0, 1 )       = ( affineElementVertices[2] - affineElementVertices[0] )( 0 );
      A( 1, 0 )       = ( affineElementVertices[1] - affineElementVertices[0] )( 1 );
      A( 1, 1 )       = ( affineElementVertices[2] - affineElementVertices[0] )( 1 );
      const auto Ainv = A.inverse();

      const Eigen::Matrix< real_t, 2, 1 > affineCoordsTranslated = affineCoordinates - affineElementVertices[0];

      const Eigen::Matrix< real_t, 2, 1 > refPos = Ainv * affineCoordsTranslated;

      Point2D midpoint( { 0., 0. } );
      for ( uint_t i = 0; i < 3; i++ )
         for ( uint_t d = 0; d < 2; d++ )
            midpoint[d] += affineElementVertices[i][d] / 3.;

      // evaluate P1 function
      std::array< uint_t, ndofsP1 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroFace( elementIndex, faceType, level, vertexDoFIndices );

      ValueType* p1Data = face.getData( u_conforming_->component( componentIdx ).getFaceDataID() )->getPointer( level );

      std::vector< real_t > dofs( ndofsP1 );
      for ( uint_t i = 0; i < ndofsP1; i++ )
      {
         dofs[i] = real_t( p1Data[vertexDoFIndices[i]] );
      }

      real_t value_conforming;
      basis_conforming_.evaluate( polyDegree, refPos, dofs, value_conforming );

      // evaluate DG function
      auto dof_discontinuous =
          u_discontinuous_->getDGFunction()->volumeDoFFunction()->dof( faceID, elementIndex, 0, faceType, level );
      real_t value_discontinuous = ValueType( dof_discontinuous ) * ( coordinates[componentIdx] - midpoint[componentIdx] );

      value = ValueType( value_conforming + value_discontinuous );
   }

   template < typename OtherValueType >
   void copyBoundaryConditionFromFunction( const EGFunction< OtherValueType >& other )
   {
      u_conforming_->copyBoundaryConditionFromFunction( *( other.getConformingPart() ) );
      u_discontinuous_->copyBoundaryConditionFromFunction( *( other.getDiscontinuousPart() ) );
   }

 protected:
   static std::vector< std::reference_wrapper< const P0Function< ValueType > > >
       filter_discontinuous( const std::vector< std::reference_wrapper< const EGFunction< ValueType > > >& functions )
   {
      std::vector< std::reference_wrapper< const P0Function< ValueType > > > dg_list;
      for ( auto& f : functions )
      {
         dg_list.push_back( *( f.get().getDiscontinuousPart() ) );
      }
      return dg_list;
   }

   static std::vector< std::reference_wrapper< const P1VectorFunction< ValueType > > >
       filter_conforming( const std::vector< std::reference_wrapper< const EGFunction< ValueType > > >& functions )
   {
      std::vector< std::reference_wrapper< const P1VectorFunction< ValueType > > > conforming_list;
      for ( auto& f : functions )
      {
         conforming_list.push_back( *( f.get().getConformingPart() ) );
      }
      return conforming_list;
   }

 protected:
   std::shared_ptr< P1VectorFunction< ValueType > > u_conforming_;
   std::shared_ptr< P0Function< ValueType > >       u_discontinuous_;

   dg::DGBasisLinearLagrange_Example basis_conforming_;
   EGBasis                           basis_discontinuous_;
};

template < typename ValueType >
EGFunction< ValueType >::EGFunction( const std::string&                                name,
                                     const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                     uint_t                                            minLevel,
                                     uint_t                                            maxLevel,
                                     BoundaryCondition                                 bc )
: Function< EGFunction< ValueType > >( name, storage, minLevel, maxLevel )
, u_conforming_{ std::make_shared< P1VectorFunction< ValueType > >( name + "_conforming",
                                                                    storage,
                                                                    minLevel,
                                                                    maxLevel,
                                                                    bc,
                                                                    true ) }
, u_discontinuous_{ std::make_shared< P0Function< ValueType > >( name + "_discontinuous", storage, minLevel, maxLevel, bc ) }
, basis_conforming_()
, basis_discontinuous_()
{}

void applyDirichletBC( const EGFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level );

} // namespace hyteg
