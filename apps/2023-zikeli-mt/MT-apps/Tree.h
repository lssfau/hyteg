/*
 * Copyright (c) 2024 Michael Zikeli
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

#include <iostream>
#include <vector>

#include "core/extern/json.hpp"

using json = nlohmann::json;

/**
 * @brief Represents a node in a tree data structure.
 *
 * @tparam T Type of data stored in the node.
 */
template < typename T >
class TreeNode
{
 public:
   T                        data;       /**< Data associated with the node */
   std::string              label;      /**< Label of the node */
   std::string              groupLabel; /**< Label of the group this node belongs to */
   std::vector< TreeNode* > children;   /**< Children nodes */

   /**
     * @brief Constructor for TreeNode class.
     *
     * @param value Data value to be stored in the node.
     * @param lbl Label of the node. Default is "root".
     * @param lblGroup Label of the group this node belongs to.
     */
   TreeNode( T value, const std::string& lbl = "root", const std::string& lblGroup = "root" )
   : data( value )
   , label( lbl )
   , groupLabel( lblGroup )
   {}
};

/**
 * @brief Represents a tree data structure.
 *
 * @tparam T Type of data stored in the tree nodes.
 */
template < typename T >
class Tree
{
 private:
   TreeNode< T >* root; /**< Pointer to the root node of the tree */

   /**
     * @brief Recursively destroys the tree nodes starting from the given node.
     *
     * @param node Pointer to the starting node.
     */
   void destroy( TreeNode< T >* node )
   {
      if ( node )
      {
         for ( TreeNode< T >* child : node->children )
         {
            destroy( child );
         }
         delete node;
      }
   }

 public:
   /**
     * @brief Constructor for Tree class.
     */
   Tree()
   : root( nullptr )
   {}

   /**
     * @brief Destructor for Tree class.
     *
     * Calls the destroy function to deallocate memory of all nodes in the tree.
     */
   ~Tree() { destroy( root ); }

   /**
     * @brief Inserts a new node with the given value and label as a child of the specified parent node.
     *
     * @param value Data value to be stored in the new node.
     * @param label Label of the new node.
     * @param groupLabel Label of the group this node belongs to.
     * @param parent Pointer to the parent node.
     * @return Pointer to the newly inserted node.
     */
   TreeNode< T >* insert( const T& value, const std::string& label, const std::string& groupLabel, TreeNode< T >* parent )
   {
      TreeNode< T >* newNode = new TreeNode< T >( value, label, groupLabel );
      if ( !root )
      {
         root = newNode;
      }
      else
      {
         parent->children.push_back( newNode );
      }
      return newNode;
   }

   /**
     * @brief Traverses the tree starting from the given node and prints the label and data of each node.
     *
     * @param node Pointer to the starting node.
     * @param depth Depth of the node in the tree.
     */
   void traverse( TreeNode< T >* node, int depth = 0 ) const
   {
      if ( node )
      {
         std::cout << std::string( depth * 2, ' ' ) << node->label << ":\t" << node->data << std::endl;
         for ( TreeNode< T >* child : node->children )
         {
            traverse( child, depth + 1 );
         }
      }
   }

   /**
     * @brief Gets the pointer to the root node of the tree.
     *
     * @return Pointer to the root node.
     */
   TreeNode< T >* getRoot() const { return root; }

   /**
     * @brief Converts the tree starting from the given node into JSON format.
     *
     * @param node Pointer to the starting node.
     * @return JSON representation of the tree.
     */
   json toJSON( TreeNode< T >* node ) const
   {
      json j;
      j["data"]       = node->data;
      j["label"]      = node->label;
      j["groupLabel"] = node->groupLabel;
      for ( TreeNode< T >* child : node->children )
      {
         j["children"].push_back( toJSON( child ) );
      }
      return j;
   }

   /**
     * @brief Converts the entire tree into JSON format.
     *
     * @return JSON representation of the tree.
     */
   json toJSON() const
   {
      if ( !root )
         return json();
      return toJSON( root );
   }
};
