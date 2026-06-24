/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/types/Concepts.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_c;

// forward declare VolumeDoFFunction class to avoid cyclic dependencies with MicroMesh
namespace volumedofspace {
template < typename >
class VolumeDoFFunction;
}

/// Class representing an FE function of continuous Lagrangre type of order 3
///
/// The P3Function clas simplements support for the continuous Lagrange element of order 3
/// in a nodal basis representation. It is currently limited to triangles.
///
/// The P3 Lagrange element for triangles has 10 local degrees of freedom. In a nodal basis
/// representation these are associated with the three vertices, two equidistantly spaced
/// nodes on each edge and the barycenter of the triangle.
/// 
/// In order to distinguish between the two degrees of freedom and have consistency between
/// neighbouring micro-cells, we use the following approach, see also Scroggs et al.,
/// "Construction of Arbitrary Order Finite Element Degree-of-Freedom Maps on Polygonal and
/// Polyhedral Cell Meshes" (2022).
/// 
/// We index the three vertices of a micro-cell in such a way that for gray triangles their
/// ordering corresponds to the ordering of the vertices of the macro-face. For blue
/// triangles we assign an index of 0 to the lower vertex and then go anti-clockwise.
/// 
/// We make the definition that an edge is always oriented from the micro-vertex with the
/// smaller index to that with the larger index. Hence, the indexing described above will
/// give us a consistent edge orientation for the micro-edgesbetween neighbouring micro-faces.
/// 
/// Along an edge we call the degree freedom closer to the vertex with the smaller index as
/// blue and the one farther away as red. We index these degrees of freedom in the same fashion
/// used by FEniCS for the P2 element, going first over the blue dofs and then over the red ones.
/// The face degree of freedom finally gets index 9.
/// <center>
/// <img src="P3_Local_Indexing.png" width="600">
/// </center>
template < typename ValueType >
class P3Function final : public Function< P3Function< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = P3Function< VType >;

   P3Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel );

   P3Function( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          boundaryCondition );

   uint_t getDimension() const override final { return 1u; }

   P3Function< ValueType >& operator[]( uint_t idx )
   {
      WALBERLA_ASSERT( idx == 0 );
      WALBERLA_UNUSED( idx );
      return *this;
   }

   const P3Function< ValueType >& operator[]( uint_t idx ) const
   {
      WALBERLA_ASSERT( idx == 0 );
      WALBERLA_UNUSED( idx );
      return *this;
   }

   inline const vertexdof::VertexDoFFunction< ValueType >&      getVertexDoFFunction() const { return vertexDoFFunction_; }
   inline const EdgeDoFFunction< ValueType >&                   getEdgeDoFFunctionBlue() const { return edgeDoFFunctionBlue_; }
   inline const EdgeDoFFunction< ValueType >&                   getEdgeDoFFunctionRed() const { return edgeDoFFunctionRed_; }
   inline const volumedofspace::VolumeDoFFunction< ValueType >& getFaceDoFFunction() const { return faceDoFFunction_; }

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const override final;

   /// @name Member functions for interpolation using DoFType flags
   ///@{
   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::vector< std::function< ValueType( const Point3D& ) > >& expr,
                     uint_t                                                             level,
                     DoFType                                                            flag = All ) const
   {
      WALBERLA_ASSERT_EQUAL( expr.size(), 1 );
      this->interpolate( expr[0], level, flag );
   };
   ///@}

   /// @name Member functions for (MPI or local) communication
   ///@{
   template < typename SenderType, typename ReceiverType >
   inline void startCommunication( const uint_t& level ) const
   {
      WALBERLA_CHECK_EQUAL( communicators_.count( level ),
                            1,
                            "No communicator found for level = " << level << ".\nDoes function '" << this->functionName_
                                                                 << "' exist on this level?" );
      communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void endCommunication( const uint_t& level ) const
   {
      WALBERLA_CHECK_EQUAL( communicators_.count( level ),
                            1,
                            "No communicator found for level = " << level << ".\nDoes function '" << this->functionName_
                                                                 << "' exist on this level?" );
      communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const
   {
      startCommunication< SenderType, ReceiverType >( level );
      endCommunication< SenderType, ReceiverType >( level );
   }
   ///@}

   /// Compute the product of several functions in an elementwise fashion
   ///
   /// The method takes as input a collection of functions. These are multiplied together in an elementwise fashion.
   /// The latter is to be understood not in a FE context, but in the sense of element-wise operators in matrix/array
   /// oriented languages, i.e. the product is a function of the same type as the inputs and its DoFs are formed as
   /// product of the corresponding DoFs of the input functions. The result is stored in the function object on which
   /// the method is invoked, overwritting its contents. It is safe, if the destination function is part of the product.
   ///
   /// \param functions  the functions forming the product
   /// \param level      level on which the multiplication should be computed
   /// \param flag       marks those primitives which are partaking in the computation of the product
   void multElementwise( const std::vector< std::reference_wrapper< const P3Function< ValueType > > >& functions,
                         uint_t                                                                        level,
                         DoFType                                                                       flag = All ) const;

   /// Replace values of the function by their inverses in an elementwise fashion
   void invertElementwise( uint_t level, DoFType flag = All, bool workOnHalos = false ) const;

   void add( ValueType scalar, uint_t level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                               scalars,
             const std::vector< std::reference_wrapper< const P3Function< ValueType > > >& functions,
             uint_t                                                                        level,
             DoFType                                                                       flag = All ) const;

   void swap( const P3Function< ValueType >& other, const uint_t& level, const DoFType& dofType = All ) const;

   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const P3Function< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All ) const;

   ValueType dotLocal( const P3Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const;

   ValueType dotGlobal( const P3Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const;

   /// @name Member functions for accessing/manipulating boundary conditions
   ///@{
   BoundaryCondition getBoundaryCondition() const;

   void setBoundaryCondition( BoundaryCondition bc );

   template < typename OtherFunctionValueType >
   void copyBoundaryConditionFromFunction( const P3Function< OtherFunctionValueType >& other )
   {
      setBoundaryCondition( other.getBoundaryCondition() );
   }
   ///@}

   /// \brief Copies all values function data from other to this.
   ///
   /// This method can be used safely if the other function is located on a different PrimitiveStorage.
   /// This method also works, if the storages are distributed differently.
   ///
   /// \param other another function
   /// \param level the refinement level
   /// \param localPrimitiveIDsToRank Map that contains as keys all primitive IDs of all primitives that are local regarding the
   ///                                storage of this function, and as values the MPI ranks of the processes that own these
   ///                                primitives regarding the storage of the other function
   /// \param otherPrimitiveIDsToRank Map that contains as keys all primitive IDs of all primitives that are local regarding the
   ///                                storage of the other function, and as values the MPI ranks of the processes that own these
   ///                                primitives regarding the storage this function lives on.
   ///
   void copyFrom( const P3Function< ValueType >&         other,
                  const uint_t&                          level,
                  const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const;

   /// @name Member functions for conversion to/from linear algebra representation
   ///@{
   void toVector( const P3Function< idx_t >&            numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const;

   void fromVector( const P3Function< idx_t >&            numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const;

   void enumerate( uint_t level ) const;

   void enumerate( uint_t level, ValueType& offset ) const;
   ///@}

   /// \brief Returns the number of DoFs. Performs global reduction, must be called collectively.
   ///
   /// \param level        refinement level
   /// \param communicator if required, a custom communicator can be passed
   /// \param onRootOnly   if true, the result is only returned on the root process
   /// \return
   uint_t getNumberOfGlobalDoFs( uint_t          level,
                                 const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                 const bool&     onRootOnly   = false ) const
   {
      return numberOfGlobalDoFs< P3FunctionTag >( *this->storage_, level, communicator, onRootOnly );
   }

   /// @Name Functions for getting meta-information on DoF values
   ///@{
   ValueType getMaxDoFValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   ValueType getMaxDoFMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   ValueType getMinDoFValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ///@}

 private:
   using Function< P3Function< ValueType > >::communicators_;
   using Function< P3Function< ValueType > >::additiveCommunicators_;

   vertexdof::VertexDoFFunction< ValueType > vertexDoFFunction_;
   EdgeDoFFunction< ValueType >              edgeDoFFunctionBlue_;
   EdgeDoFFunction< ValueType >              edgeDoFFunctionRed_;

   // TODO: In 2D this is fine, but in 3D we need a true FaceDoFFunction!!!!
   volumedofspace::VolumeDoFFunction< ValueType > faceDoFFunction_;

   /// Auxilliary class with static method for functionality not available in VolumeDoFFunction
   struct faceDoFHelpers
   {
      static void interpolate( const volumedofspace::VolumeDoFFunction< ValueType >& function,
                               ValueType                                             constant,
                               uint_t                                                level,
                               DoFType                                               flag );

      static void interpolate( const volumedofspace::VolumeDoFFunction< ValueType >& function,
                               const std::function< ValueType( const Point3D& ) >&   expr,
                               uint_t                                                level,
                               DoFType                                               flag );

      static void invertElementwise( const volumedofspace::VolumeDoFFunction< ValueType >& function, uint_t level, DoFType flag );

      static void multElementwise(
          const volumedofspace::VolumeDoFFunction< ValueType >&                                                function,
          const std::vector< std::reference_wrapper< const volumedofspace::VolumeDoFFunction< ValueType > > >& srcFunctions,
          uint_t                                                                                               level,
          DoFType                                                                                              flag );
   };
};

extern template class P3Function< double >;
extern template class P3Function< int >;
extern template class P3Function< idx_t >;

} //namespace hyteg
