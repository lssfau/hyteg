/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGVectorFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/experimental/P2PlusBubbleFunction.hpp"
#include "hyteg/experimental/P2PlusBubbleVectorFunction.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionMultiStore.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

/// Class to hold different Finite Element functions
///
/// This class allows registering Finite Element functions of different kind
/// and different value types (real_t, int64_t, ...). This functionality is
/// for example used by the VTKOutput class. For every supported kind of
/// FE function there exists a corresponding getter method which returns a
/// reference to the associated FunctionMultiStore onject.
class FEFunctionRegistry
{
 public:
   /// Add an FE Function to the registry
   template < template < typename > class func_t, typename value_t >
   inline void add( const func_t< value_t >& function )
   {
      // -------------
      //  CGFunctions
      // -------------

      // P1Functions
      if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > )
      {
         p1Functions_.add( function );
      }

      // P1VectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
      {
         p1VecFunctions_.add( function );
      }

      // P2Functions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > )
      {
         p2Functions_.add( function );
      }

      // P2VectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
      {
         p2VecFunctions_.add( function );
      }

      // -------------
      //  DGFunctions
      // -------------

      // P0Functions
      else if constexpr ( std::is_same_v< func_t< value_t >, P0Function< value_t > > )
      {
         dgFunctions_.add( *function.getDGFunction() );
      }

      // DG1Functions
      else if constexpr ( std::is_same_v< func_t< value_t >, DG1Function< value_t > > )
      {
         dgFunctions_.add( *function.getDGFunction() );
      }

      // DGFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, dg::DGFunction< value_t > > )
      {
         dgFunctions_.add( function );
      }

      // DGVectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, dg::DGVectorFunction< value_t > > )
      {
         dgVecFunctions_.add( function );
      }

      // ---------------------------------
      //  Special and Composite Functions
      // ---------------------------------

      // P2PlusBubbleFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2PlusBubbleFunction< value_t > > )
      {
         p2PlusBubbleFunctions_.add( function );
      }

      // P2PlusBubbleVectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2PlusBubbleVectorFunction< value_t > > )
      {
         p2PlusBubbleVecFunctions_.add( function );
      }

      // EdgeDoFFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, EdgeDoFFunction< value_t > > )
      {
         edgeDoFFunctions_.add( function );
      }

      // N1E1VectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, n1e1::N1E1VectorFunction< value_t > > )
      {
         n1e1Functions_.add( function );
      }

      // EGFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, EGFunction< value_t > > )
      {
         p1dgeVecFunctions_.add( function );
      }

      // P1StokesFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, P1StokesFunction< value_t > > )
      {
         this->add( function.uvw() );
         this->add( function.p() );
      }

      // P2P1TaylorHoodFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
      {
         this->add( function.uvw() );
         this->add( function.p() );
      }

      // EGP0StokesFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, EGP0StokesFunction< value_t > > )
      {
         this->add( function.uvw() );
         this->add( function.p() );
      }

      // -----------------------
      //  "Technical" Functions
      // -----------------------

      // BlockFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, BlockFunction< value_t > > )
      {
         for ( uint_t k = 0; k < function.getNumberOfBlocks(); k++ )
         {
            this->add( function[k] );
         }
      }

      else if constexpr ( std::is_base_of_v< BlockFunction< value_t >, func_t< value_t > > )
      {
         for ( uint_t k = 0; k < function.getNumberOfBlocks(); k++ )
         {
            this->add( function[k] );
         }
      }

      // GenericFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, GenericFunction< value_t > > )
      {
         addGenericFunction( function );
      }

      // NO MATCH !!!
      else
      {
         WALBERLA_ABORT( "Could not add function of type '" << FunctionTrait< func_t< value_t > >::getTypeName()
                                                            << "' to FEFunctionRegistry!" );
      }
   }

   /// Remove an FE Function from the registry
   template < template < typename > class func_t, typename value_t >
   inline void remove( const func_t< value_t >& function )
   {
      // -------------
      //  CGFunctions
      // -------------

      // P1Functions
      if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > )
      {
         p1Functions_.remove( function );
      }

      // P1VectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
      {
         p1VecFunctions_.remove( function );
      }

      // P2Functions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > )
      {
         p2Functions_.remove( function );
      }

      // P2VectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
      {
         p2VecFunctions_.remove( function );
      }

      // -------------
      //  DGFunctions
      // -------------

      // P0Functions
      else if constexpr ( std::is_same_v< func_t< value_t >, P0Function< value_t > > )
      {
         dgFunctions_.remove( *function.getDGFunction() );
      }

      // DG1Functions
      else if constexpr ( std::is_same_v< func_t< value_t >, DG1Function< value_t > > )
      {
         dgFunctions_.remove( *function.getDGFunction() );
      }

      // DGFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, dg::DGFunction< value_t > > )
      {
         dgFunctions_.remove( function );
      }

      // DGVectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, dg::DGVectorFunction< value_t > > )
      {
         dgVecFunctions_.remove( function );
      }

      // ---------------------------------
      //  Special and Composite Functions
      // ---------------------------------

      // P2PlusBubbleFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2PlusBubbleFunction< value_t > > )
      {
         p2PlusBubbleFunctions_.remove( function );
      }

      // P2PlusBubbleVectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, P2PlusBubbleVectorFunction< value_t > > )
      {
         p2PlusBubbleVecFunctions_.remove( function );
      }

      // EdgeDoFFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, EdgeDoFFunction< value_t > > )
      {
         edgeDoFFunctions_.remove( function );
      }

      // N1E1VectorFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, n1e1::N1E1VectorFunction< value_t > > )
      {
         n1e1Functions_.remove( function );
      }

      // EGFunctions
      else if constexpr ( std::is_same_v< func_t< value_t >, EGFunction< value_t > > )
      {
         p1dgeVecFunctions_.remove( function );
      }

      // P1StokesFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, P1StokesFunction< value_t > > )
      {
         this->remove( function.uvw() );
         this->remove( function.p() );
      }

      // P2P1TaylorHoodFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
      {
         this->remove( function.uvw() );
         this->remove( function.p() );
      }

      // EGP0StokesFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, EGP0StokesFunction< value_t > > )
      {
         this->remove( function.uvw() );
         this->remove( function.p() );
      }

      // -----------------------
      //  "Technical" Functions
      // -----------------------

      // BlockFunction
      else if constexpr ( std::is_same_v< func_t< value_t >, BlockFunction< value_t > > )
      {
         for ( uint_t k = 0; k < function.getNumberOfBlocks(); k++ )
         {
            this->remove( function[k] );
         }
      }

      else if constexpr ( std::is_base_of_v< BlockFunction< value_t >, func_t< value_t > > )
      {
         for ( uint_t k = 0; k < function.getNumberOfBlocks(); k++ )
         {
            this->remove( function[k] );
         }
      }

      // NO MATCH !!!
      else
      {
         WALBERLA_ABORT( "Could not remove function of type '" << FunctionTrait< func_t< value_t > >::getTypeName()
                                                               << "' from FEFunctionRegistry!" );
      }
   }

   // clang-format off
   const FunctionMultiStore< P1Function >&                 getP1Functions()                 const { return p1Functions_;              }
   const FunctionMultiStore< P2Function >&                 getP2Functions()                 const { return p2Functions_;              }
   const FunctionMultiStore< P1VectorFunction >&           getP1VectorFunctions()           const { return p1VecFunctions_;           }
   const FunctionMultiStore< P2VectorFunction >&           getP2VectorFunctions()           const { return p2VecFunctions_;           }
   const FunctionMultiStore< EdgeDoFFunction >&            getEdgeDoFFunctions()            const { return edgeDoFFunctions_;         }
   const FunctionMultiStore< dg::DGFunction >&             getDGFunctions()                 const { return dgFunctions_;              }
   const FunctionMultiStore< dg::DGVectorFunction >&       getDGVectorFunctions()           const { return dgVecFunctions_;           }
   const FunctionMultiStore< n1e1::N1E1VectorFunction >&   getN1E1VectorFunctions()         const { return n1e1Functions_;            }
   const FunctionMultiStore< EGFunction >&                 getEGFunctions()                 const { return p1dgeVecFunctions_;        }
   const FunctionMultiStore< P2PlusBubbleFunction >&       getP2PlusBubbleFunctions()       const { return p2PlusBubbleFunctions_;    }
   const FunctionMultiStore< P2PlusBubbleVectorFunction >& getP2PlusBubbleVectorFunctions() const { return p2PlusBubbleVecFunctions_; }
   // clang-format on

   /// return the total number of registered functions
   ///
   /// \note The returned number is the number of functions stored inside the object.
   /// This number can be different than the sum of functions which were registered,
   /// since e.g. composite functions are disassembled into components when placed into
   /// the object's internal FunctionMultiStore objects.
   uint_t getNumberOfFunctionsInRegistry() const
   {
      uint_t num{ 0u };
      num += p1Functions_.size();
      num += p2Functions_.size();
      num += p2PlusBubbleFunctions_.size();
      num += p1VecFunctions_.size();
      num += p2VecFunctions_.size();
      num += p2PlusBubbleVecFunctions_.size();
      num += edgeDoFFunctions_.size();
      num += dgFunctions_.size();
      num += dgVecFunctions_.size();
      num += n1e1Functions_.size();
      num += p1dgeVecFunctions_.size();
      return num;
   }

   template < template < class > class func_t >
   const FunctionMultiStore< func_t >& getFunctions()
   {
      if constexpr ( std::is_same_v< func_t< real_t >, P1Function< real_t > > )
      {
         return p1Functions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, P2Function< real_t > > )
      {
         return p2Functions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, P2PlusBubbleFunction< real_t > > )
      {
         return p2PlusBubbleFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, P1VectorFunction< real_t > > )
      {
         return p1VecFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, P2VectorFunction< real_t > > )
      {
         return p2VecFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, P2PlusBubbleVectorFunction< real_t > > )
      {
         return p2PlusBubbleVecFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, EdgeDoFFunction< real_t > > )
      {
         return edgeDoFFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, dg::DGFunction< real_t > > )
      {
         return dgFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, dg::DGVectorFunction< real_t > > )
      {
         return dgVecFunctions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, n1e1::N1E1VectorFunction< real_t > > )
      {
         return n1e1Functions_;
      }
      else if constexpr ( std::is_same_v< func_t< real_t >, EGFunction< real_t > > )
      {
         return p1dgeVecFunctions_;
      }
      else
      {
         WALBERLA_ABORT( "Unsupported Function Kind in FEFunctionRegistry::getFunctions()!" );
      }
   }

   /// Append names of all registered functions of a certain kind to the provided vector.
   void extractFunctionNames( std::vector< std::string >& names, functionTraits::FunctionKind funcKind ) const
   {
      std::vector< std::string > namesFound;

      switch ( funcKind )
      {
      case functionTraits::P1_FUNCTION:
         namesFound = p1Functions_.getFunctionNames();
         break;
      case functionTraits::P1_VECTOR_FUNCTION:
         namesFound = p1VecFunctions_.getFunctionNames();
         break;
      case functionTraits::P2_FUNCTION:
         namesFound = p2Functions_.getFunctionNames();
         break;
      case functionTraits::P2_PLUS_BUBBLE_FUNCTION:
         namesFound = p2PlusBubbleFunctions_.getFunctionNames();
         break;
      case functionTraits::P2_VECTOR_FUNCTION:
         namesFound = p2VecFunctions_.getFunctionNames();
         break;
      case functionTraits::P2_PLUS_BUBBLE_VECTOR_FUNCTION:
         namesFound = p2PlusBubbleVecFunctions_.getFunctionNames();
         break;
      case functionTraits::EDGE_DOF_FUNCTION:
         namesFound = edgeDoFFunctions_.getFunctionNames();
         break;
      case functionTraits::DG_FUNCTION:
         namesFound = dgFunctions_.getFunctionNames();
         break;
      case functionTraits::DG_VECTOR_FUNCTION:
         namesFound = dgVecFunctions_.getFunctionNames();
         break;
      case functionTraits::N1E1_VECTOR_FUNCTION:
         namesFound = n1e1Functions_.getFunctionNames();
         break;
      case functionTraits::EG_FUNCTION:
         namesFound = p1dgeVecFunctions_.getFunctionNames();
         break;
      default:
         WALBERLA_ABORT( "Unimplemented case found in FEFunctionRegistry::extractFunctionNames()!" );
      }

      names.insert( names.end(), namesFound.begin(), namesFound.end() );
   }

 private:
   FunctionMultiStore< P1Function >                 p1Functions_;
   FunctionMultiStore< P2Function >                 p2Functions_;
   FunctionMultiStore< P1VectorFunction >           p1VecFunctions_;
   FunctionMultiStore< P2VectorFunction >           p2VecFunctions_;
   FunctionMultiStore< EdgeDoFFunction >            edgeDoFFunctions_;
   FunctionMultiStore< dg::DGFunction >             dgFunctions_;
   FunctionMultiStore< dg::DGVectorFunction >       dgVecFunctions_;
   FunctionMultiStore< n1e1::N1E1VectorFunction >   n1e1Functions_;
   FunctionMultiStore< EGFunction >                 p1dgeVecFunctions_;
   FunctionMultiStore< P2PlusBubbleFunction >       p2PlusBubbleFunctions_;
   FunctionMultiStore< P2PlusBubbleVectorFunction > p2PlusBubbleVecFunctions_;

   template < typename value_t >
   void addGenericFunction( const GenericFunction< value_t >& function )
   {
      bool matchFound = false;
      switch ( function.getFunctionKind() )
      {
      case functionTraits::P1_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< P1Function< value_t > > >( function );
         break;

      case functionTraits::P2_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< P2Function< value_t > > >( function );
         break;

      case functionTraits::P2_PLUS_BUBBLE_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< P2PlusBubbleFunction< value_t > > >( function );
         break;

      case functionTraits::P1_VECTOR_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< P1VectorFunction< value_t > > >( function );
         break;

      case functionTraits::P2_VECTOR_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< P2VectorFunction< value_t > > >( function );
         break;

      case functionTraits::P2_PLUS_BUBBLE_VECTOR_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< P2PlusBubbleVectorFunction< value_t > > >( function );
         break;

      case functionTraits::DG_VECTOR_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< dg::DGVectorFunction< value_t > > >( function );
         break;

      case functionTraits::EDGE_DOF_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< EdgeDoFFunction< value_t > > >( function );
         break;

      case functionTraits::DG_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< dg::DGFunction< value_t > > >( function );
         break;

      case functionTraits::N1E1_VECTOR_FUNCTION:
         matchFound = tryUnwrapAndAdd< FunctionWrapper< n1e1::N1E1VectorFunction< value_t > > >( function );
         break;

      default:
         matchFound = false;
      }

      if ( !matchFound )
      {
         WALBERLA_ABORT( "FEFunctionRegistry: Failed to add GenericFunction object!" );
      }
   }

   template < typename WrapperFunc, typename value_t >
   bool tryUnwrapAndAdd( const GenericFunction< value_t >& function )
   {
      bool               success = false;
      const WrapperFunc* aux     = dynamic_cast< const WrapperFunc* >( &function );
      if ( aux != nullptr )
      {
         add( aux->unwrap() );
         success = true;
      }
      return success;
   }
};

} // namespace hyteg
