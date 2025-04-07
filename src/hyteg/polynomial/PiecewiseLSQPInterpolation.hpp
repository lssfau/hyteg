/*
* Copyright (c) 2025 Nils Kohl.
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

#include <iostream>
#include <vector>

#include "core/extern/json.hpp"
#include "core/math/extern/exprtk.h"

#include "hyteg/dataexport/VTKOutput/VTKHexahedraOutput.hpp"

#include "LSQPInterpolator.hpp"

/// @file PiecewiseLSQPInterpolation.hpp
///
/// Piecewise least squares polynomial function approximation.
///
/// This file contains data structures to produce a piecewise polynomial function approximation.
/// A single high-order polynomial might be insufficient and the coefficient computation can be slow.
/// A kd-tree is used for space partitioning (subdividing an axis-aligned bounding box (AABB) that captures
/// the targeted domain). Each subdomain then approximates the function locally with a polynomial.
/// To improve accuracy, the polynomial interpolation can be configured to also include points that are located
/// outside the bounds of the local subdomain by temporary extension of the subdomain for construction
/// of the polynomial. During evaluation at some point P, the unique enclosing subdomain is found with a tree-search
/// and then the local polynomial is evaluated. No continuity between subdomains is enforced, though.
///
/// One application is the memory efficient approximate transfer of data between different meshes.
/// The KDTree equipped with the polynomials can be serialized and deserialized from and to file.

namespace hyteg {

/// @class AABB
/// @brief Represents an axis-aligned bounding box (AABB) in n-dimensional space.
template < int Dim >
class AABB
{
 public:
   static constexpr int D = Dim;

   /// @brief Constructs an AABB with specified minimum and maximum coordinates.
   ///
   /// @param aabbMin The minimum coordinates of the bounding box.
   /// @param aabbMax The maximum coordinates of the bounding box.
   AABB( const PointND< real_t, Dim >& aabbMin, const PointND< real_t, Dim >& aabbMax )
   : min_( aabbMin )
   , max_( aabbMax )
   {
      computeLargestAxis();
   }

   PointND< real_t, Dim > min() const { return min_; }
   PointND< real_t, Dim > max() const { return max_; }
   int                    largestAxis() const { return largestAxis_; }

   /// @brief Checks if a given point is contained within the AABB.
   ///
   /// @param point The point to check.
   /// @return True if the point is within the AABB, false otherwise.
   bool contains( const PointND< real_t, Dim >& point ) const
   {
      for ( size_t i = 0; i < min_.size(); ++i )
      {
         if ( point[i] < min_[i] || point[i] > max_[i] )
         {
            return false;
         }
      }
      return true;
   }

   /// @brief Extends the AABB by a given factor. The center of the AABB does not change.
   ///
   /// @param factor The factor by which to extend the AABB.
   /// @return A new AABB that is extended by the specified factor.
   AABB extend( real_t factor ) const
   {
      PointND< real_t, Dim > newMin;
      PointND< real_t, Dim > newMax;
      for ( size_t i = 0; i < min_.size(); ++i )
      {
         real_t center     = ( min_[i] + max_[i] ) / 2;
         real_t halfExtent = ( max_[i] - min_[i] ) / 2 * factor;
         newMin[i]         = center - halfExtent;
         newMax[i]         = center + halfExtent;
      }
      return AABB( newMin, newMax );
   }

   /// @brief Returns the vertices of an AABB.
   ///
   /// 2D: in counter-clockwise order
   /// 3D: first the bottom vertices in counter-clockwise order, then the top vertices in counter-clockwise order
   ///
   /// @return A vector of vectors, where each inner vector represents a vertex of the AABB.
   std::vector< PointND< real_t, Dim > > getVertices() const
   {
      if constexpr ( Dim == 2 )
      {
         return { { min_[0], min_[1] }, { min_[0], max_[1] }, { max_[0], max_[1] }, { max_[0], min_[1] } };
      }
      else if constexpr ( Dim == 3 )
      {
         return { { min_[0], min_[1], min_[2] },
                  { min_[0], max_[1], min_[2] },
                  { max_[0], max_[1], min_[2] },
                  { max_[0], min_[1], min_[2] },
                  { min_[0], min_[1], max_[2] },
                  { min_[0], max_[1], max_[2] },
                  { max_[0], max_[1], max_[2] },
                  { max_[0], min_[1], max_[2] } };
      }

      WALBERLA_ABORT( "Invalid AABB dim." )
   }

   void serialize( SendBuffer& sendBuffer ) const
   {
      sendBuffer << min_;
      sendBuffer << max_;
      sendBuffer << largestAxis_;
   }

   void deserialize( RecvBuffer& recvBuffer )
   {
      recvBuffer >> min_;
      recvBuffer >> max_;
      recvBuffer >> largestAxis_;
   }

 private:
   void computeLargestAxis()
   {
      largestAxis_         = 0;
      real_t largestExtent = max_[0] - min_[0];
      for ( int i = 1; i < Dim; ++i )
      {
         real_t extent = max_[i] - min_[i];
         if ( extent > largestExtent )
         {
            largestExtent = extent;
            largestAxis_  = i;
         }
      }
   }

   PointND< real_t, Dim > min_; ///< Minimum coordinates of the bounding box.
   PointND< real_t, Dim > max_; ///< Maximum coordinates of the bounding box.
   int                    largestAxis_;
};

/// @class KDTreeNode
/// @brief Represents a node in a KD-Tree, which is a space-partitioning data structure for organizing points in k-dimensional space.
template < typename T, int AABBDim >
class KDTreeNode
{
 public:
   using ValueType = T;
   using AABBType  = AABB< AABBDim >;

   std::unique_ptr< KDTreeNode< T, AABBDim > > left;   ///< Left child node.
   std::unique_ptr< KDTreeNode< T, AABBDim > > right;  ///< Right child node.
   AABB< AABBDim >                             bounds; ///< Bounding box of the node.
   bool                                        isLeaf; ///< Indicates if the node is a leaf.
   int                                         depth;  ///< Depth of the node in the tree.
   std::optional< T >                          data;   ///< Optional data of type T associated with the node.

   /// @brief Constructs a KDTreeNode with specified bounds and depth.
   ///
   /// @param nodeBounds The bounding box of the node.
   /// @param nodeDepth The depth of the node in the tree.
   KDTreeNode( const AABB< AABBDim >& nodeBounds, int nodeDepth )
   : bounds( nodeBounds )
   , isLeaf( true )
   , depth( nodeDepth )
   , data( std::nullopt )
   {}

   void serialize( SendBuffer& sendBuffer ) const
   {
      bounds.serialize( sendBuffer );
      sendBuffer << isLeaf;
      sendBuffer << depth;

      sendBuffer << data.has_value();
      if ( data.has_value() )
      {
         data.value().serialize( sendBuffer );
      }

      sendBuffer << !!left;
      if ( left )
      {
         left->serialize( sendBuffer );
      }

      sendBuffer << !!right;
      if ( right )
      {
         right->serialize( sendBuffer );
      }
   }

   void deserialize( RecvBuffer& recvBuffer )
   {
      bounds.deserialize( recvBuffer );
      recvBuffer >> isLeaf;
      recvBuffer >> depth;

      bool hasData;
      recvBuffer >> hasData;
      if ( hasData )
      {
         // T has to implement a ctor
         //
         //     explicit T( RecvBuffer & recvBuffer )
         //
         // Too lazy to write a C++ concept.
         data = std::optional< T >( T( recvBuffer ) );
      }

      bool hasLeft;
      recvBuffer >> hasLeft;
      if ( hasLeft )
      {
         // Dummy arguments
         left = std::make_unique< KDTreeNode >( bounds, depth );
         left->deserialize( recvBuffer );
      }

      bool hasRight;
      recvBuffer >> hasRight;
      if ( hasRight )
      {
         // Dummy arguments
         right = std::make_unique< KDTreeNode >( bounds, depth );
         right->deserialize( recvBuffer );
      }
   }
};

/// @class KDTree
/// @brief Represents a KD-Tree, a space-partitioning data structure for organizing points in k-dimensional space.
template < typename T, int AABBDim >
class KDTree
{
 public:
   using NodeType = KDTreeNode< T, AABBDim >;

   std::unique_ptr< KDTreeNode< T, AABBDim > > root;       ///< Root node of the KD-Tree.
   int                                         maxDepth;   ///< Maximum depth of the tree.
   int                                         dimensions; ///< Number of dimensions.

   /// @brief Constructs a KDTree with specified bounds and maximum depth.
   ///
   /// @param bounds The bounding box of the root node.
   /// @param treeMaxDepth The maximum depth of the tree.
   KDTree( const AABB< AABBDim >& bounds, int treeMaxDepth )
   : maxDepth( treeMaxDepth )
   , dimensions( bounds.min().size() )
   {
      root = std::make_unique< KDTreeNode< T, AABBDim > >( bounds, 0 );
      buildTree( root, 0 );
   }

   /// @brief Finds the leaf node containing a given point.
   ///
   /// @param point The point to search for.
   /// @return A pointer to the leaf node containing the point, or nullptr if not found.
   KDTreeNode< T, AABBDim >* findLeafNode( const PointND< real_t, AABBDim >& point ) const
   {
      return findLeafNodeRecursive( root.get(), point );
   }

   /// @brief Finds all leaf nodes whose extended bounds contain a given point.
   ///
   /// @param point The point to search for.
   /// @param factor The factor by which to extend the bounds.
   /// @return A vector of pointers to leaf nodes containing the point.
   std::vector< KDTreeNode< T, AABBDim >* > findLeafNodesWithExtendedBounds( const PointND< real_t, AABBDim >& point,
                                                                             real_t                            factor ) const
   {
      std::vector< KDTreeNode< T, AABBDim >* > result;
      findLeafNodesWithExtendedBoundsRecursive( root.get(), point, factor, result );
      return result;
   }

   /// @brief Iterates over all nodes and applies a callback function.
   ///
   /// @param callback The callback function to apply to each node.
   void iterate( const std::function< void( KDTreeNode< T, AABBDim >* ) >& callback ) const
   {
      iterateRecursive( root.get(), callback );
   }

   /// @brief Prints the structure of the KD-Tree.
   void printTree() const { printTreeRecursive( root.get(), "" ); }

   /// @brief Gets the number of leaf nodes in the KD-Tree.
   ///
   /// @return The number of leaf nodes.
   int getNumberOfLeaves() const { return 1 << maxDepth; }

   /// @brief Serializes the KDTree into a JSON object.
   ///
   /// @return A JSON object representing the KDTree.
   void serialize( SendBuffer& sendBuffer ) const
   {
      sendBuffer << maxDepth;
      sendBuffer << dimensions;
      root->serialize( sendBuffer );
   }

   /// @brief Deserializes a KDTree from a JSON object.
   ///
   /// @param j The JSON object representing the KDTree.
   /// @return A deserialized KDTree object.
   void deserialize( RecvBuffer& recvBuffer )
   {
      recvBuffer >> maxDepth;
      recvBuffer >> dimensions;
      root->deserialize( recvBuffer );
   }

 private:
   /// @brief Recursively builds the KD-Tree.
   ///
   /// @param node The current node being built.
   /// @param depth The current depth in the tree.
   void buildTree( std::unique_ptr< KDTreeNode< T, AABBDim > >& node, int depth )
   {
      if ( depth >= maxDepth )
      {
         return;
      }

      node->isLeaf = false;
      int    axis  = node->bounds.largestAxis();
      real_t mid   = ( node->bounds.min()[axis] + node->bounds.max()[axis] ) / 2;

      auto leftBoundsMin  = node->bounds.min();
      auto leftBoundsMax  = node->bounds.max();
      leftBoundsMax[axis] = mid;
      AABB leftBounds( leftBoundsMin, leftBoundsMax );
      node->left = std::make_unique< KDTreeNode< T, AABBDim > >( leftBounds, depth + 1 );
      buildTree( node->left, depth + 1 );

      auto rightBoundsMin  = node->bounds.min();
      auto rightBoundsMax  = node->bounds.max();
      rightBoundsMin[axis] = mid;
      AABB rightBounds( rightBoundsMin, rightBoundsMax );
      node->right = std::make_unique< KDTreeNode< T, AABBDim > >( rightBounds, depth + 1 );
      buildTree( node->right, depth + 1 );
   }

   /// @brief Recursively finds the leaf node containing a given point.
   ///
   /// @param node The current node.
   /// @param point The point to search for.
   /// @return A pointer to the leaf node containing the point, or nullptr if not found.
   KDTreeNode< T, AABBDim >* findLeafNodeRecursive( KDTreeNode< T, AABBDim >*         node,
                                                    const PointND< real_t, AABBDim >& point ) const
   {
      if ( node == nullptr || !node->bounds.contains( point ) )
      {
         return nullptr;
      }

      if ( node->isLeaf )
      {
         return node;
      }

      int axis = node->bounds.largestAxis();
      if ( point[axis] <= ( node->bounds.min()[axis] + node->bounds.max()[axis] ) / 2 )
      {
         return findLeafNodeRecursive( node->left.get(), point );
      }
      else
      {
         return findLeafNodeRecursive( node->right.get(), point );
      }
   }

   /// @brief Recursively finds all leaf nodes whose extended bounds contain a given point.
   ///
   /// @param node The current node.
   /// @param point The point to search for.
   /// @param factor The factor by which to extend the bounds.
   /// @param result A vector to store the resulting leaf nodes.
   void findLeafNodesWithExtendedBoundsRecursive( KDTreeNode< T, AABBDim >*                 node,
                                                  const PointND< real_t, AABBDim >&         point,
                                                  real_t                                    factor,
                                                  std::vector< KDTreeNode< T, AABBDim >* >& result ) const
   {
      if ( node == nullptr )
         return;

      AABB extendedBounds = node->bounds.extend( factor );
      if ( extendedBounds.contains( point ) )
      {
         if ( node->isLeaf )
         {
            result.push_back( node );
         }
         else
         {
            findLeafNodesWithExtendedBoundsRecursive( node->left.get(), point, factor, result );
            findLeafNodesWithExtendedBoundsRecursive( node->right.get(), point, factor, result );
         }
      }
   }

   /// @brief Recursively iterates over all nodes and applies a callback function.
   ///
   /// @param node The current node.
   /// @param callback The callback function to apply to each node.
   void iterateRecursive( KDTreeNode< T, AABBDim >*                                 node,
                          const std::function< void( KDTreeNode< T, AABBDim >* ) >& callback ) const
   {
      if ( node == nullptr )
         return;

      callback( node );

      iterateRecursive( node->left.get(), callback );
      iterateRecursive( node->right.get(), callback );
   }

   /// @brief Recursively prints the structure of the KD-Tree.
   ///
   /// @param node The current node.
   /// @param prefix The prefix string for formatting the output.
   void printTreeRecursive( KDTreeNode< T, AABBDim >* node, const std::string& prefix ) const
   {
      if ( node == nullptr )
         return;

      std::cout << prefix;

      std::string nodeType = node->isLeaf ? "Leaf" : "Internal";
      std::cout << "|-- " << nodeType << " Node at depth " << node->depth << "\n";

      std::cout << prefix << "    Bounds: ";
      for ( size_t i = 0; i < node->bounds.min.size(); ++i )
      {
         std::cout << "[" << node->bounds.min[i] << ", " << node->bounds.max[i] << "]";
         if ( i < node->bounds.min.size() - 1 )
         {
            std::cout << " x ";
         }
      }
      std::cout << "\n";

      if ( node->isLeaf && node->data.has_value() )
      {
         // std::cout << prefix << "    Data: " << node->data.value() << "\n";
         std::cout << prefix << "    has data\n";
      }

      printTreeRecursive( node->left.get(), prefix + "    " );
      printTreeRecursive( node->right.get(), prefix + "    " );
   }
};

/// @brief Simple struct carrying an LSQP object and a polynomial.
/// To be used as data item in the KDTree for piecewise interpolation.
struct LSQPPolyPair2D
{
   using LSQP_T = LSQPInterpolator< MonomialBasis2D, LSQPType::VERTEX >;
   LSQP_T                          lsqp;
   Polynomial2D< MonomialBasis2D > poly;

   explicit LSQPPolyPair2D( uint_t degree )
   : lsqp( degree )
   , poly( degree )
   {}

   explicit LSQPPolyPair2D( RecvBuffer& recvBuffer )
   : lsqp( 0 )
   , poly( 0 )
   {
      deserialize( recvBuffer );
   }

   void serialize( SendBuffer& sendBuffer ) const { sendBuffer << poly; }

   void deserialize( RecvBuffer& recvBuffer ) { recvBuffer >> poly; }
};

/// @brief Simple struct carrying an LSQP object and a polynomial.
/// To be used as data item in the KDTree for piecewise interpolation.
struct LSQPPolyPair3D
{
   using LSQP_T = LSQPInterpolator3D< MonomialBasis3D, LSQPType::VERTEX >;
   LSQP_T                          lsqp;
   Polynomial3D< MonomialBasis3D > poly;

   explicit LSQPPolyPair3D( uint_t degree )
   : lsqp( degree )
   , poly( degree )
   {}

   explicit LSQPPolyPair3D( RecvBuffer& recvBuffer )
   : lsqp( 0 )
   , poly( 0 )
   {
      deserialize( recvBuffer );
   }

   void serialize( SendBuffer& sendBuffer ) const { sendBuffer << poly; }

   void deserialize( RecvBuffer& recvBuffer ) { recvBuffer >> poly; }
};

/// @brief Type alias for a KDTree equipped with necessary data for piecewise interpolation (2D polynomials).
using PiecewiseLSQPPolyKDTree2D = KDTree< LSQPPolyPair2D, 2 >;
/// @brief Type alias for a KDTree equipped with necessary data for piecewise interpolation (3D polynomials).
using PiecewiseLSQPPolyKDTree3D = KDTree< LSQPPolyPair3D, 3 >;

/// @brief Helper function to initialize the polynomial data on all leaves.
///
/// You probably want to use setupPiecewiseLSQPPolyKDTree().
template < typename PiecewiseLSQPPolyKDTree >
void addPolynomialsToLeaves( PiecewiseLSQPPolyKDTree& kdTree, uint_t degree )
{
   auto addPolyToLeaf = [&degree]( typename PiecewiseLSQPPolyKDTree::NodeType* node ) {
      if ( node->isLeaf )
      {
         node->data = typename PiecewiseLSQPPolyKDTree::NodeType::ValueType( degree );
      }
   };
   kdTree.iterate( addPolyToLeaf );
}

/// @brief Returns a function to be passed to interpolate that adds all points to the LSQP of the corresponding KDTree leave
/// (i.e., subdomain).
///
/// You probably want to use setupPiecewiseLSQPPolyKDTree().
///
/// Example usage:
/// @code
///    auto lsqpPointAdder = interpolationPointAdder( kdTree, aabbExtensionFactor );
///    u.interpolate( lsqpPointAdder, { u }, level );
/// @endcode
template < typename PiecewiseLSQPPolyKDTree >
std::function< real_t( const Point3D& p, const std::vector< real_t >& feFuncData ) >
    interpolationPointAdder( PiecewiseLSQPPolyKDTree& kdTree, real_t aabbExtensionFactor )
{
   auto func = [&kdTree, aabbExtensionFactor]( const Point3D& x, const std::vector< real_t >& us ) {
      if constexpr ( std::is_same_v< PiecewiseLSQPPolyKDTree, PiecewiseLSQPPolyKDTree2D > )
      {
         auto containingLeaves = kdTree.findLeafNodesWithExtendedBounds( { x( 0 ), x( 1 ) }, aabbExtensionFactor );
         for ( auto containingLeaf : containingLeaves )
         {
            containingLeaf->data.value().lsqp.addInterpolationPoint( Point2D( x( 0 ), x( 1 ) ), us[0] );
         }
      }
      else
      {
         auto containingLeaves = kdTree.findLeafNodesWithExtendedBounds( { x( 0 ), x( 1 ), x( 2 ) }, aabbExtensionFactor );
         for ( auto containingLeaf : containingLeaves )
         {
            containingLeaf->data.value().lsqp.addInterpolationPoint( x, us[0] );
         }
      }
      return us[0];
   };

   return func;
}

/// @brief Solves all LSQP problems in the leaves.
///
/// You probably want to use setupPiecewiseLSQPPolyKDTree().
template < typename PiecewiseLSQPPolyKDTree >
void solveLSQPs( PiecewiseLSQPPolyKDTree& kdTree )
{
   auto constructPolys = []( typename PiecewiseLSQPPolyKDTree::NodeType* node ) {
      if ( !node->isLeaf )
      {
         return;
      }
      if ( node->data.value().lsqp.numInterpolationPoints() == 0 )
      {
         return;
      }
      node->data.value().lsqp.interpolate( node->data.value().poly );
   };
   kdTree.iterate( constructPolys );
}

/// @brief Sets up the piecewise LSQP kd-tree.
///
/// Essentially just a helper that calls
/// - addPolynomialsToLeaves( ... )         // initializes polynomials
/// - interpolationPointAdder( ... )        // adds interpolation points
/// - solveLSQPs( ... )                     // solves the LSQPs on all leaves of the kd-tree
/// to set up the tree for evaluation.
///
/// @tparam PiecewiseLSQPPolyKDTree The type of the piecewise LSQP polynomial kd-tree.
/// @tparam FEFunctionType The type of the finite element function used for interpolation.
///
/// @param kdTree The kd-tree to set up for interpolation.
/// @param degree The polynomial degree to use for the piecewise polynomials.
/// @param aabbExtensionFactor The factor by which to extend the axis-aligned bounding box (AABB) around each leaf.
/// @param level The FE function refinement level.
/// @param feFunction The finite element function that shall be approximated.
template < typename PiecewiseLSQPPolyKDTree, typename FEFunctionType >
void setupPiecewiseLSQPPolyKDTree( PiecewiseLSQPPolyKDTree& kdTree,
                                   uint_t                   degree,
                                   real_t                   aabbExtensionFactor,
                                   uint_t                   level,
                                   FEFunctionType&          feFunction )
{
   addPolynomialsToLeaves( kdTree, degree );
   auto lsqpPointAdder = interpolationPointAdder( kdTree, aabbExtensionFactor );
   feFunction.interpolate( lsqpPointAdder, { feFunction }, level );
   solveLSQPs( kdTree );
}

/// @brief Returns a function to evaluate a piecewise KDTree polynomial.
///
/// Example usage:
/// @code
///    u.interpolate( piecewiseKDTreePolyEvaluator( kdTree ), level );
/// @endcode
template < typename PiecewiseLSQPPolyKDTree >
std::function< real_t( const Point3D& ) > piecewiseKDTreePolyEvaluator( const PiecewiseLSQPPolyKDTree& kdTree )
{
   auto evalKDPoly = [&kdTree]( const Point3D& x ) {
      if constexpr ( std::is_same_v< PiecewiseLSQPPolyKDTree, PiecewiseLSQPPolyKDTree2D > )
      {
         auto leaf = kdTree.findLeafNode( { x( 0 ), x( 1 ) } );
         WALBERLA_ASSERT( leaf->data.has_value() );
         return leaf->data.value().poly.eval( Point2D( x( 0 ), x( 1 ) ) );
      }
      else
      {
         auto leaf = kdTree.findLeafNode( { x( 0 ), x( 1 ), x( 2 ) } );
         WALBERLA_ASSERT( leaf->data.has_value() );
         return leaf->data.value().poly.eval( x );
      }
   };
   return evalKDPoly;
}

/// @brief Writes a VTK file in parallel with the AABBs of the KDTree with the specified depth.
///
/// If you only want to write one tree - only call this on root.
template < typename T, int AABBDim >
void writeKDTreeAABBsToVTK( const KDTree< T, AABBDim >& kdTree,
                            uint_t                      depth,
                            const std::string&          dir,
                            const std::string&          filename )
{
   VTKHexahedraOutput vtk( dir, filename );
   auto               addNodeAABB = [&vtk, depth]( typename KDTree< T, AABBDim >::NodeType* node ) {
      if ( node->depth == depth )
      {
         if ( AABBDim == 2 )
         {
            const auto               vertices = node->bounds.getVertices();
            std::array< Point3D, 4 > rectPoints;
            for ( int i = 0; i < 4; ++i )
            {
               for ( int d = 0; d < 2; ++d )
               {
                  rectPoints[i][d] = vertices[i][d];
               }
            }
            vtk.addRectangle( rectPoints );
         }
         else if ( AABBDim == 3 )
         {
            const auto               vertices = node->bounds.getVertices();
            std::array< Point3D, 8 > hexaPoints;
            for ( int i = 0; i < 8; ++i )
            {
               for ( int d = 0; d < 3; ++d )
               {
                  hexaPoints[i][d] = vertices[i][d];
               }
            }
            vtk.addHexahedron( hexaPoints );
         }
      }
   };
   kdTree.iterate( addNodeAABB );
   vtk.write();
}

} // namespace hyteg