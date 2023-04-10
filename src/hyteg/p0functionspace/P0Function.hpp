/*
* Copyright (c) 2017-2023 Ponsuganth Ilangovan, Nils Kohl, Marcus Mohr.
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
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using namespace dg;

template < typename ValueType >
class P0Function : public Function< P0Function< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = P0Function< VType >;

   P0Function( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          bc )
   : Function< P0Function< ValueType > >( name, storage, minLevel, maxLevel )
   {
      auto basis  = std::make_shared< DGBasisLinearLagrange_Example >();
      dgFunction_ = std::make_shared< DGFunction< ValueType > >( name, storage, minLevel, maxLevel, basis, 0, bc );
   }

   P0Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : P0Function( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}


   const std::shared_ptr< DGFunction< ValueType > > getDGFunction() const { return dgFunction_; }

   void setBoundaryCondition( BoundaryCondition bc ) { dgFunction_->setBoundaryCondition( bc ); }

   BoundaryCondition getBoundaryCondition() const { return dgFunction_->getBoundaryCondition(); }

   uint_t getDimension() const { return dgFunction_->getDimension(); };
   void setDoNotWarnOnInterpolateFlag() { doNotWarnOnInterpolate_ = true; }

   // template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const { dgFunction_->communicate( level ); }


   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const { dgFunction_->add( scalar, level, flag ); };

   void add( const std::vector< ValueType >                                                scalars,
             const std::vector< std::reference_wrapper< const P0Function< ValueType > > >& functions,
             uint_t                                                                        level,
             DoFType                                                                       flag = All ) const
   {
      std::vector< ValueType > new_scalars( scalars );
      new_scalars.push_back( 1.0 );
      std::vector< std::reference_wrapper< const P0Function< ValueType > > > new_functions( functions );
      new_functions.push_back( *this );
      assign( new_scalars, new_functions, level, flag );
   };

   void multElementwise( const std::vector< std::reference_wrapper< const P0Function< ValueType > > >& functions,
                         uint_t                                                                        level,
                         DoFType                                                                       flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void interpolate( ValueType constant, uint_t level, DoFType dofType = All ) const
   {
      if( !doNotWarnOnInterpolate_ ) WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::interpolate() 'interpolates' values at the centroid." );
      if ( this->storage_->hasGlobalCells() )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            const auto  cellID = it.first;
            const auto& cell   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( cellID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level );

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
                      constant;
               }
            }
         }
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            const auto  faceID = it.first;
            const auto& face   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( faceID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] = constant;
               }
            }
         }
      }
   }

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType dofType = All ) const
   {
      WALBERLA_UNUSED( dofType );
      if( !doNotWarnOnInterpolate_ ) WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::interpolate() 'interpolates' values at the centroid." );

      if ( this->storage_->hasGlobalCells() )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            const auto  cellID = it.first;
            const auto& cell   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( cellID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level );

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  const std::array< indexing::Index, 4 > vertexIndices =
                      celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );
                  std::array< Eigen::Matrix< real_t, 3, 1 >, 4 > elementVertices;
                  for ( uint_t i = 0; i < 4; i++ )
                  {
                     const auto elementVertex = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                     elementVertices[i]( 2 )  = elementVertex[2];
                  }

                  const Eigen::Matrix< real_t, 3, 1 > centroid =
                      ( elementVertices[0] + elementVertices[1] + elementVertices[2] + elementVertices[3] ) / real_c( 4 );

                  const auto val = expr( Point3D( centroid( 0 ), centroid( 1 ), centroid( 2 ) ) );

                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
                      ValueType( val );
               }
            }
         }
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            const auto  faceID = it.first;
            const auto& face   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( faceID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

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

                  const Eigen::Matrix< real_t, 2, 1 > centroid =
                      ( elementVertices[0] + elementVertices[1] + elementVertices[2] ) / real_c( 3 );

                  const auto val = expr( Point3D( centroid( 0 ), centroid( 1 ), 0 ) );

                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] =
                      ValueType( val );
               }
            }
         }
      }
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expression,
                     const std::vector< std::reference_wrapper< const P0Function< ValueType > > >&        srcFunctions,
                     uint_t                                                                               level,
                     DoFType                                                                              flag = All ) const
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::interpolate() 'interpolates' values at the centroid." );
      if ( this->storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "Not implemented" );
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            const auto  faceID = it.first;
            const auto& face   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( faceID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

            std::vector< ValueType* > dofs;
            dofs.reserve( srcFunctions.size() + 1 );

            std::vector< volumedofspace::indexing::VolumeDoFMemoryLayout > memLayouts;
            memLayouts.reserve( srcFunctions.size() + 1 );

            std::vector< ValueType > srcValues;
            srcValues.reserve( srcFunctions.size() );

            dofs.push_back( getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level ) );
            memLayouts.push_back( getDGFunction()->volumeDoFFunction()->memoryLayout() );

            for ( const auto src : srcFunctions )
            {
               dofs.push_back( src.get().getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level ) );
               memLayouts.push_back( src.get().getDGFunction()->volumeDoFFunction()->memoryLayout() );
            }

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

                  const Eigen::Matrix< real_t, 2, 1 > centroid =
                      ( elementVertices[0] + elementVertices[1] + elementVertices[2] ) / 3.;

                  for ( size_t k = 0; k < srcFunctions.size(); ++k )
                  {
                     srcValues[k] =
                         dofs[k + 1]
                             [volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayouts[k + 1] )];
                  }

                  const auto newValue = expression( Point3D( centroid( 0 ), centroid( 1 ), 0 ), srcValues );

                  dofs[0][volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayouts[0] )] =
                      ValueType( newValue );
               }
            }
         }
      }
   };

   void swap( const P0Function< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
       WALBERLA_UNUSED(flag);
       dgFunction_->swap(*other.getDGFunction(),level,flag);
   };

   void copyFrom( const P0Function< ValueType >&         other,
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

   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const P0Function< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All ) const
   {
      std::vector< std::reference_wrapper< const DGFunction< ValueType > > > dgFunctions;
      for ( auto f : functions )
      {
         dgFunctions.push_back( *f.get().getDGFunction() );
      }
      dgFunction_->assign( scalars, dgFunctions, level, flag );
   }

   ValueType dotGlobal( const P0Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      return dgFunction_->dotGlobal( *rhs.getDGFunction(), level );
   }

   ValueType sumGlobal( uint_t level, const DoFType& flag = All ) const { return dgFunction_->sumGlobal( level ); }

   ValueType dotLocal( const P0Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      return dgFunction_->dotLocal( *rhs.getDGFunction(), level );
   }

   void enumerate( uint_t level ) const { dgFunction_->enumerate( level ); }

   void enumerate( uint_t level, ValueType& offset ) const { dgFunction_->enumerate( level, offset ); }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const P0Function< idx_t >&            numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      dgFunction_->toVector( *numerator.getDGFunction(), vec, level, flag );
   }

   void fromVector( const P0Function< idx_t >&            numerator,
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

    template < typename OtherValueType >
    void copyBoundaryConditionFromFunction( const P0Function< OtherValueType >& other )
    {
        dgFunction_->copyBoundaryConditionFromFunction( *other.getDGFunction() );
    }

   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const
   {
      if ( flag != All && flag != Inner )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::getMaxMagnitude -> DoFType flag will be ignored!" );
      }
      return dgFunction_->getMaxMagnitude( level, mpiReduce );
   }

   ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const
   {
      if ( flag != All && flag != Inner )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::getMaxValue -> DoFType flag will be ignored!" );
      }
      return dgFunction_->getMax( level, mpiReduce );
   }

   ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const
   {
      if ( flag != All && flag != Inner )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::getMinValue -> DoFType flag will be ignored!" );
      }
      return dgFunction_->getMin( level, mpiReduce );
   }

 private:
   std::shared_ptr< DGFunction< ValueType > > dgFunction_;
   bool doNotWarnOnInterpolate_ = false;

};

namespace dg {
/// \brief removes the mean value of func to project it onto the space of functions with mean value 0
///
/// \param func             [in] function to project
/// \param level            [in] refinement level
inline void projectMean( P0Function< real_t >& func, const uint_t& level )
{
   if ( func.isDummy() )
   {
      return;
   }
   const uint_t numGlobalVertices = func.getNumberOfGlobalDoFs( level );
   const real_t sum               = func.sumGlobal( level, Inner );
   func.add( -sum / ( real_c( numGlobalVertices ) ), level, Inner );
}
} // namespace dg

} // namespace hyteg
