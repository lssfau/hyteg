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

#include "terraneo/helpers/ConvectionToolbox.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

namespace convectionToolbox {

///////////////////////////////////////////////////////
/////////////////////// General ///////////////////////
///////////////////////////////////////////////////////

std::ostream& operator<<( std::ostream& out, const ConvectionBC& cbc )
{
   for ( const auto& [i, v] : cbc.DoFTypes_ )
   {
      out << '[' << i << "] = " << ( (size_t) v ) << "; ";
   }
   return out;
}

// 0 = X, 1 = Y, 2 = Z, 3 = T
std::vector< BoundaryUID > ConvectionBC::getUIDsWithName( std::string name, std::set< uint_t > types )
{
   std::transform( name.begin(), name.end(), name.begin(), []( unsigned char c ) { return std::tolower( c ); } );

   std::vector< BoundaryUID > ids;

   for ( const auto& [i, v] : BoundaryUIDs_ )
   {
      if ( ( types.find( v.first ) != types.end() ) && ( i.find( name ) != std::string::npos ) )
      {
         ids.push_back( v.second );
      }
   }

   return ids;
}

std::pair< P2P1TaylorHoodFunction< real_t >, P2Function< real_t > >
    ConvectionBC::createVisualisation( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   P2P1TaylorHoodFunction< real_t > up = createP2P1TaylorHoodFunction( "BC_up", storage, level, level, *this );
   P2Function< real_t >             T  = createP2TemperatureFunction( "BC_T", storage, level, level, *this );

   // up.interpolate(0,level,Inner);
   // up.interpolate(2,level,DirichletBoundary);
   // up.interpolate(4,level,NeumannBoundary);
   // up.interpolate(8,level,FreeslipBoundary);

   // T.interpolate(0,level,Inner);
   // T.interpolate(2,level,DirichletBoundary);
   // T.interpolate(4,level,NeumannBoundary);
   // T.interpolate(8,level,FreeslipBoundary);

   for ( const auto& [i, v] : BoundaryUIDs_ )
   {
      if ( v.first <= 2 )
      {
         up.uvw()[v.first].interpolate( (real_t) DoFTypes_.at( i ), level, v.second );
      }
      T.interpolate( (real_t) DoFTypes_.at( i ), level, v.second );
   }

   return std::pair< P2P1TaylorHoodFunction< real_t >, P2Function< real_t > >( up, T );
}

ConvectionBC::ConvectionBC( ConvectionMeshType meshType )
: meshType_( meshType )
{}

ConvectionBC::ConvectionBC( ConvectionMeshType meshType,
                            BoundaryCondition  bcVelocityX,
                            BoundaryCondition  bcVelocityY,
                            BoundaryCondition  bcVelocityZ,
                            BoundaryCondition  bcTemperature )
: meshType_( meshType )
, bcVelocityX_( bcVelocityX )
, bcVelocityY_( bcVelocityY )
, bcVelocityZ_( bcVelocityZ )
, bcTemperature_( bcTemperature )
{}

void ConvectionBC::setBC( BoundaryCondition bcVelocityX,
                          BoundaryCondition bcVelocityY,
                          BoundaryCondition bcVelocityZ,
                          BoundaryCondition bcTemperature )
{
   bcVelocityX_   = bcVelocityX;
   bcVelocityY_   = bcVelocityY;
   bcVelocityZ_   = bcVelocityZ;
   bcTemperature_ = bcTemperature;
}

real_t projectPressureZeroMean( const P1Function< real_t >& p,
                                const P1Function< real_t >& One,
                                uint_t                      level,
                                real_t                      VolumeOmega,
                                DoFType                     flag )
{
   real_t integral = One.dotGlobal( p, level, flag );

   return integral / VolumeOmega;
}

void createVelocityBoundaryConditionAndUID( std::string                             name,
                                            DoFType                                 t,
                                            uint_t                                  flag,
                                            std::pair< uint_t, BoundaryCondition >& bc,
                                            ConvectionBC&                           cbc )
{
   BoundaryUID id;

   switch ( t )
   {
   case DoFType::DirichletBoundary:
      id = bc.second.createDirichletBC( name.c_str(), flag );
      break;
   case DoFType::NeumannBoundary:
      id = bc.second.createNeumannBC( name.c_str(), flag );
      break;
   case DoFType::FreeslipBoundary:
      id = bc.second.createFreeslipBC( name.c_str(), flag );
      break;
   default:
      id = bc.second.createInnerBC( name.c_str(), flag );
      break;
   }

   cbc.BoundaryUIDs_.insert( { name, { bc.first, id } } );
   cbc.DoFTypes_.insert( { name, t } );
}

// t1 is the type of the edge
void create2DEdgeVelocityBoundaryConditionAndUID( std::string                             name,
                                                  DoFType                                 t,
                                                  uint_t                                  flag,
                                                  std::pair< uint_t, BoundaryCondition >& bc1,
                                                  std::pair< uint_t, BoundaryCondition >& bc2,
                                                  ConvectionBC&                           cbc,
                                                  bool                                    invert )
{
   std::string name1 = name + "X";
   std::string name2 = name + "Y";
   if ( t == DoFType::FreeslipBoundary )
   {
      if ( invert )
      {
         createVelocityBoundaryConditionAndUID( name2.c_str(), DoFType::DirichletBoundary, flag, bc2, cbc );
         createVelocityBoundaryConditionAndUID( name1.c_str(), DoFType::Inner, flag, bc1, cbc );
      }
      else
      {
         createVelocityBoundaryConditionAndUID( name1.c_str(), DoFType::DirichletBoundary, flag, bc1, cbc );
         createVelocityBoundaryConditionAndUID( name2.c_str(), DoFType::Inner, flag, bc2, cbc );
      }
   }
   else
   {
      createVelocityBoundaryConditionAndUID( name1.c_str(), t, flag, bc1, cbc );
      createVelocityBoundaryConditionAndUID( name2.c_str(), t, flag, bc2, cbc );
   }
}

// t1 and t2 are the types of the edges joined at the corner
void create2DCornerVelocityBoundaryConditionAndUID( std::string                             name,
                                                    DoFType                                 t1,
                                                    DoFType                                 t2,
                                                    uint_t                                  flag,
                                                    std::pair< uint_t, BoundaryCondition >& bc1,
                                                    std::pair< uint_t, BoundaryCondition >& bc2,
                                                    ConvectionBC&                           cbc )
{
   std::string name1 = name + "X";
   std::string name2 = name + "Y";
   if ( t1 == DoFType::FreeslipBoundary && t2 == DoFType::NeumannBoundary )
   {
      createVelocityBoundaryConditionAndUID( name1.c_str(), DoFType::DirichletBoundary, flag, bc1, cbc );
      createVelocityBoundaryConditionAndUID( name2.c_str(), t2, flag, bc2, cbc );
   }
   else if ( t1 == DoFType::NeumannBoundary && t2 == DoFType::FreeslipBoundary )
   {
      createVelocityBoundaryConditionAndUID( name1.c_str(), t1, flag, bc1, cbc );
      createVelocityBoundaryConditionAndUID( name2.c_str(), DoFType::DirichletBoundary, flag, bc2, cbc );
   }
   else
   {
      // if one of the two edges is Dirichlet or both are Freeslip then the edge is Dirichlet
      createVelocityBoundaryConditionAndUID( name1.c_str(), DoFType::DirichletBoundary, flag, bc1, cbc );
      createVelocityBoundaryConditionAndUID( name2.c_str(), DoFType::DirichletBoundary, flag, bc2, cbc );
   }
}

// t1 is the type of the face, axis 1 is the normal direction of the face
void create3DFaceVelocityBoundaryConditionAndUID( std::string                             name,
                                                  DoFType                                 t,
                                                  uint_t                                  flag,
                                                  std::pair< uint_t, BoundaryCondition >& bc1,
                                                  std::pair< uint_t, BoundaryCondition >& bc2,
                                                  std::pair< uint_t, BoundaryCondition >& bc3,
                                                  ConvectionBC&                           cbc,
                                                  uint_t                                  axis1 = 1,
                                                  uint_t                                  axis2 = 2,
                                                  uint_t                                  axis3 = 3 )
{
   std::string name1 = name + "X";
   std::string name2 = name + "Y";
   std::string name3 = name + "Z";

   if ( t == DoFType::FreeslipBoundary )
   {
      std::map< uint_t, std::pair< uint_t, BoundaryCondition >& > permBC   = { { 1, bc1 }, { 2, bc2 }, { 3, bc3 } };
      std::map< uint_t, std::string& >                            permName = { { 1, name1 }, { 2, name2 }, { 3, name3 } };

      createVelocityBoundaryConditionAndUID(
          permName.at( axis1 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis1 ), cbc );
      createVelocityBoundaryConditionAndUID( permName.at( axis2 ).c_str(), DoFType::Inner, flag, permBC.at( axis2 ), cbc );
      createVelocityBoundaryConditionAndUID( permName.at( axis3 ).c_str(), DoFType::Inner, flag, permBC.at( axis3 ), cbc );
   }
   else
   {
      createVelocityBoundaryConditionAndUID( name1.c_str(), t, flag, bc1, cbc );
      createVelocityBoundaryConditionAndUID( name2.c_str(), t, flag, bc2, cbc );
      createVelocityBoundaryConditionAndUID( name3.c_str(), t, flag, bc3, cbc );
   }
}

// t1 and t2 are the types of the faces joined by the edge
// axis 1 is the normal direction of face1, axis 2 is the normal direction of face2, axis3 is the remaining one
void create3DEdgeVelocityBoundaryConditionAndUID( std::string                             name,
                                                  DoFType                                 t1,
                                                  DoFType                                 t2,
                                                  uint_t                                  flag,
                                                  std::pair< uint_t, BoundaryCondition >& bc1,
                                                  std::pair< uint_t, BoundaryCondition >& bc2,
                                                  std::pair< uint_t, BoundaryCondition >& bc3,
                                                  ConvectionBC&                           cbc,
                                                  uint_t                                  axis1,
                                                  uint_t                                  axis2,
                                                  uint_t                                  axis3 )
{
   std::string                                                 name1    = name + "X";
   std::string                                                 name2    = name + "Y";
   std::string                                                 name3    = name + "Z";
   std::map< uint_t, std::pair< uint_t, BoundaryCondition >& > permBC   = { { 1, bc1 }, { 2, bc2 }, { 3, bc3 } };
   std::map< uint_t, std::string& >                            permName = { { 1, name1 }, { 2, name2 }, { 3, name3 } };

   if ( t1 == DoFType::FreeslipBoundary && t2 == DoFType::FreeslipBoundary )
   {
      createVelocityBoundaryConditionAndUID(
          permName.at( axis1 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis1 ), cbc );
      createVelocityBoundaryConditionAndUID(
          permName.at( axis2 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis2 ), cbc );
      createVelocityBoundaryConditionAndUID( permName.at( axis3 ).c_str(), DoFType::Inner, flag, permBC.at( axis3 ), cbc );
   }
   else
   {
      if ( t1 == DoFType::FreeslipBoundary && t2 == DoFType::NeumannBoundary )
      {
         createVelocityBoundaryConditionAndUID(
             permName.at( axis1 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis1 ), cbc );
         createVelocityBoundaryConditionAndUID(
             permName.at( axis2 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis2 ), cbc );
         createVelocityBoundaryConditionAndUID(
             permName.at( axis3 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis3 ), cbc );
      }
      else if ( t1 == DoFType::NeumannBoundary && t2 == DoFType::FreeslipBoundary )
      {
         createVelocityBoundaryConditionAndUID(
             permName.at( axis1 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis1 ), cbc );
         createVelocityBoundaryConditionAndUID(
             permName.at( axis2 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis2 ), cbc );
         createVelocityBoundaryConditionAndUID(
             permName.at( axis3 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis3 ), cbc );
      }
      else
      {
         // if one of the two faces is Dirichlet then the edge is Dirichlet
         createVelocityBoundaryConditionAndUID(
             permName.at( axis1 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis1 ), cbc );
         createVelocityBoundaryConditionAndUID(
             permName.at( axis2 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis2 ), cbc );
         createVelocityBoundaryConditionAndUID(
             permName.at( axis3 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis3 ), cbc );
      }
   }
}

// t1, t2 and t3 are the types of the faces joined at the corner
// axis 1 is the normal direction of face1, axis 2 is the normal direction of face2, axis3 is the remaining one
void create3DCornerVelocityBoundaryConditionAndUID( std::string                             name,
                                                    DoFType                                 t1,
                                                    DoFType                                 t2,
                                                    DoFType                                 t3,
                                                    uint_t                                  flag,
                                                    std::pair< uint_t, BoundaryCondition >& bc1,
                                                    std::pair< uint_t, BoundaryCondition >& bc2,
                                                    std::pair< uint_t, BoundaryCondition >& bc3,
                                                    ConvectionBC&                           cbc,
                                                    uint_t                                  axis1,
                                                    uint_t                                  axis2,
                                                    uint_t                                  axis3 )
{
   std::string                                                 name1    = name + "X";
   std::string                                                 name2    = name + "Y";
   std::string                                                 name3    = name + "Z";
   std::map< uint_t, std::pair< uint_t, BoundaryCondition >& > permBC   = { { 1, bc1 }, { 2, bc2 }, { 3, bc3 } };
   std::map< uint_t, std::string& >                            permName = { { 1, name1 }, { 2, name2 }, { 3, name3 } };

   if ( ( t1 == DoFType::DirichletBoundary || t2 == DoFType::DirichletBoundary || t3 == DoFType::DirichletBoundary ) ||
        ( t1 == DoFType::FreeslipBoundary && t2 == DoFType::FreeslipBoundary && t3 == DoFType::FreeslipBoundary ) )
   {
      // if one of the three faces is Dirichlet then the corner is Dirichlet
      // if all three faces are freeslip than the corner is Dirichlet
      createVelocityBoundaryConditionAndUID(
          permName.at( axis1 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis1 ), cbc );
      createVelocityBoundaryConditionAndUID(
          permName.at( axis2 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis2 ), cbc );
      createVelocityBoundaryConditionAndUID(
          permName.at( axis3 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis3 ), cbc );
   }
   else
   {
      switch ( t1 )
      {
      case DoFType::FreeslipBoundary:
         createVelocityBoundaryConditionAndUID(
             permName.at( axis1 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis1 ), cbc );
         break;
      case DoFType::NeumannBoundary:
         createVelocityBoundaryConditionAndUID(
             permName.at( axis1 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis1 ), cbc );
         break;
      default:
         break;
      }

      switch ( t2 )
      {
      case DoFType::FreeslipBoundary:
         createVelocityBoundaryConditionAndUID(
             permName.at( axis2 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis2 ), cbc );
         break;
      case DoFType::NeumannBoundary:
         createVelocityBoundaryConditionAndUID(
             permName.at( axis2 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis2 ), cbc );
         break;
      default:
         break;
      }

      switch ( t3 )
      {
      case DoFType::FreeslipBoundary:
         createVelocityBoundaryConditionAndUID(
             permName.at( axis3 ).c_str(), DoFType::DirichletBoundary, flag, permBC.at( axis3 ), cbc );
         break;
      case DoFType::NeumannBoundary:
         createVelocityBoundaryConditionAndUID(
             permName.at( axis3 ).c_str(), DoFType::NeumannBoundary, flag, permBC.at( axis3 ), cbc );
         break;
      default:
         break;
      }
   }
}

void createTemperatureBoundaryConditionAndUID( std::string        name,
                                               DoFType            t,
                                               uint_t             flag,
                                               BoundaryCondition& bc,
                                               ConvectionBC&      cbc )
{
   BoundaryUID id;
   std::string nameT = name + "T";

   switch ( t )
   {
   case DoFType::DirichletBoundary:
      id = bc.createDirichletBC( nameT.c_str(), flag );
      break;
   case DoFType::NeumannBoundary:
      id = bc.createNeumannBC( nameT.c_str(), flag );
      break;
   case DoFType::FreeslipBoundary:
      id = bc.createFreeslipBC( nameT.c_str(), flag );
      break;
   default:
      id = BoundaryUID( nameT );
      break;
   }

   cbc.BoundaryUIDs_.insert( { nameT, { 3, id } } );
   cbc.DoFTypes_.insert( { nameT, t } );
}

P2P1TaylorHoodFunction< real_t > createP2P1TaylorHoodFunction( const std::string&                         name,
                                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                                               size_t                                     minLevel,
                                                               size_t                                     maxLevel,
                                                               ConvectionBC&                              cbc )
{
   P2P1TaylorHoodFunction< real_t > f( name, storage, minLevel, maxLevel );

   f.uvw().component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f.uvw().component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f.uvw().component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

P1StokesFunction< real_t > createP1StokesFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   size_t                                     minLevel,
                                                   size_t                                     maxLevel,
                                                   ConvectionBC&                              cbc )
{
   P1StokesFunction< real_t > f( name, storage, minLevel, maxLevel );

   f.uvw().component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f.uvw().component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f.uvw().component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

P1VectorFunction< real_t > createP1VectorVelocityFunction( const std::string&                         name,
                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                           size_t                                     minLevel,
                                                           size_t                                     maxLevel,
                                                           ConvectionBC&                              cbc )
{
   P1VectorFunction< real_t > f( name, storage, minLevel, maxLevel );

   f.component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f.component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f.component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

P2VectorFunction< real_t > createP2VectorVelocityFunction( const std::string&                         name,
                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                           size_t                                     minLevel,
                                                           size_t                                     maxLevel,
                                                           ConvectionBC&                              cbc )
{
   P2VectorFunction< real_t > f( name, storage, minLevel, maxLevel );

   f.component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f.component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f.component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

P1Function< real_t > createP1TemperatureFunction( const std::string&                         name,
                                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                                  size_t                                     minLevel,
                                                  size_t                                     maxLevel,
                                                  ConvectionBC&                              cbc )
{
   P1Function< real_t > f( name, storage, minLevel, maxLevel );

   f.setBoundaryCondition( cbc.bcTemperature_ );

   return f;
}

P2Function< real_t > createP2TemperatureFunction( const std::string&                         name,
                                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                                  size_t                                     minLevel,
                                                  size_t                                     maxLevel,
                                                  ConvectionBC&                              cbc )
{
   P2Function< real_t > f( name, storage, minLevel, maxLevel );

   f.setBoundaryCondition( cbc.bcTemperature_ );

   return f;
}

////////////////

std::shared_ptr< P2P1TaylorHoodFunction< real_t > >
    createSharedP2P1TaylorHoodFunction( const std::string&                         name,
                                        const std::shared_ptr< PrimitiveStorage >& storage,
                                        size_t                                     minLevel,
                                        size_t                                     maxLevel,
                                        ConvectionBC&                              cbc )
{
   auto f = std::make_shared< P2P1TaylorHoodFunction< real_t > >( name, storage, minLevel, maxLevel );

   f->uvw().component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f->uvw().component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f->uvw().component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

std::shared_ptr< P1VectorFunction< real_t > >
    createSharedP1VectorVelocityFunction( const std::string&                         name,
                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                          size_t                                     minLevel,
                                          size_t                                     maxLevel,
                                          ConvectionBC&                              cbc )
{
   auto f = std::make_shared< P1VectorFunction< real_t > >( name, storage, minLevel, maxLevel );

   f->component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f->component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f->component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

std::shared_ptr< P2VectorFunction< real_t > >
    createSharedP2VectorVelocityFunction( const std::string&                         name,
                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                          size_t                                     minLevel,
                                          size_t                                     maxLevel,
                                          ConvectionBC&                              cbc )
{
   auto f = std::make_shared< P2VectorFunction< real_t > >( name, storage, minLevel, maxLevel );

   f->component( 0 ).setBoundaryCondition( cbc.bcVelocityX_ );
   f->component( 1 ).setBoundaryCondition( cbc.bcVelocityY_ );
   if ( storage->hasGlobalCells() )
   {
      f->component( 2 ).setBoundaryCondition( cbc.bcVelocityZ_ );
   }

   return f;
}

std::shared_ptr< P1Function< real_t > > createSharedP1TemperatureFunction( const std::string&                         name,
                                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           size_t                                     minLevel,
                                                                           size_t                                     maxLevel,
                                                                           ConvectionBC&                              cbc )
{
   auto f = std::make_shared< P1Function< real_t > >( name, storage, minLevel, maxLevel );

   f->setBoundaryCondition( cbc.bcTemperature_ );

   return f;
}

std::shared_ptr< P2Function< real_t > > createSharedP2TemperatureFunction( const std::string&                         name,
                                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           size_t                                     minLevel,
                                                                           size_t                                     maxLevel,
                                                                           ConvectionBC&                              cbc )
{
   auto f = std::make_shared< P2Function< real_t > >( name, storage, minLevel, maxLevel );

   f->setBoundaryCondition( cbc.bcTemperature_ );

   return f;
}

std::function< bool( const Point3D& x ) > createAllVerticesFct()
{
   return []( const Point3D& x ) {
      WALBERLA_UNUSED( x );
      return true;
   };
}

std::function< real_t( const Point3D& x, real_t temp ) >
    createSpaceDependentViscosityProfileWithJumps( real_t rCMB, real_t rSurface, real_t eta0 )
{
   return [=]( const Point3D& x, real_t temp ) {
      WALBERLA_UNUSED( temp );

      real_t pos = ( x.norm() - rCMB ) / ( rSurface - rCMB );

      return ( ( pos <= 3.0 / 4.0 || pos < 0 ) ?
                   ( 1.25e+23 ) :
                   ( ( pos <= 447.0 / 580.0 ) ?
                         ( 8.1653213196657478e+22 * std::pow( 8.214452314287204e-42, pos - 0.75450192648181957 ) +
                           0.78234362044754846 ) :
                         ( ( pos <= 247.0 / 290.0 ) ?
                               ( 1.5e+22 ) :
                               ( ( pos <= 101.0 / 116.0 ) ?
                                     ( 7.3660295214714627e+21 * std::pow( 2.4378729963174497e-38, pos - 0.85993667827582509 ) +
                                       0.13217190828289649 ) :
                                     ( ( pos <= 136.0 / 145.0 ) ?
                                           ( 2.5e+21 ) :
                                           ( 1.1063585614663961e+54 *
                                                 std::pow( 1.1773916059277528e+37, pos - 1.8185642274387086 ) -
                                             1.4612240520530394e+17 ) ) ) ) ) ) /
             eta0;
   };

   //    if ( pos < 0.5 )
   //    {
   //       return 1.0;
   //    }
   //    else
   //    {
   //       return 1e3;
   //    }
   // };
}

///////////////////////////////////////////////////////
/////////////////// Spherical Shell ///////////////////
///////////////////////////////////////////////////////
bool SphericalShellCMBBoundary( const Point3D& x, real_t rCMB, real_t boundaryTolerance )
{
   return x.norm() < rCMB + boundaryTolerance;
}
bool SphericalShellSurfaceBoundary( const Point3D& x, real_t rSurface, real_t boundaryTolerance )
{
   return x.norm() > rSurface - boundaryTolerance;
}

void SphericalShellBoundaryNormal( const Point3D& x, Point3D& normal, real_t rSurface, real_t rCMB )
{
   real_t radius = x.norm();
   if ( std::abs( radius - rSurface ) < std::abs( radius - rCMB ) )
   {
      normal = Point3D( { x[0] / radius, x[1] / radius, x[2] / radius } );
   }
   else
   {
      normal = Point3D( { -x[0] / radius, -x[1] / radius, -x[2] / radius } );
   }
}

std::function< bool( const Point3D& x ) > createSphericalShellSurfaceBoundaryFct( real_t rSurface, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return SphericalShellSurfaceBoundary( x, rSurface, boundaryTolerance ); };
}

std::function< bool( const Point3D& x ) > createSphericalShellCMBBoundaryFct( real_t rCMB, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return SphericalShellCMBBoundary( x, rCMB, boundaryTolerance ); };
}

std::function< void( const Point3D& in, Point3D& out ) > createSphericalShellBoundaryNormalFct( real_t rSurface, real_t rCMB )
{
   return [=]( const Point3D& x, Point3D& normal ) { SphericalShellBoundaryNormal( x, normal, rSurface, rCMB ); };
}

std::shared_ptr< PrimitiveStorage > createSphericalShellStorage( uint_t nTan,
                                                                 uint_t nRad,
                                                                 real_t rCMB,
                                                                 real_t rSurface,
                                                                 bool   blending,
                                                                 real_t boundaryTolerance,
                                                                 uint_t flagSurface,
                                                                 uint_t flagCMB,
                                                                 uint_t additionalHaloDepth )
{
   // create the spherical shell mesh
   MeshInfo              meshInfo = MeshInfo::meshSphericalShell( nTan, nRad, rCMB, rSurface );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( blending )
   {
      IcosahedralShellMap::setMap( setupStorage );
   }

   // set load balancing
   loadbalancing::roundRobinVolume( setupStorage );

   // set surface and cmb boundary flags by vertex location
   setupStorage.setMeshBoundaryFlagsByVertexLocation( 0, createAllVerticesFct(), true );
   setupStorage.setMeshBoundaryFlagsByVertexLocation(
       flagSurface, createSphericalShellSurfaceBoundaryFct( rSurface, boundaryTolerance ), true );
   setupStorage.setMeshBoundaryFlagsByVertexLocation(
       flagCMB, createSphericalShellCMBBoundaryFct( rCMB, boundaryTolerance ), true );

   // create storage
   return std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
}

ConvectionBC createSphericalShellBoundaryConditions( DoFType typeSurface, DoFType typeCMB, uint_t flagSurface, uint_t flagCMB )
{
   std::pair< uint_t, BoundaryCondition > bcVelocityX( 0, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityY( 1, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityZ( 2, BoundaryCondition() );
   BoundaryCondition                      bcTemperature;

   ConvectionBC cbc( ConvectionMeshType::SphericalShell );

   createVelocityBoundaryConditionAndUID( "surfaceX", typeSurface, flagSurface, bcVelocityX, cbc );
   createVelocityBoundaryConditionAndUID( "surfaceY", typeSurface, flagSurface, bcVelocityY, cbc );
   createVelocityBoundaryConditionAndUID( "surfaceZ", typeSurface, flagSurface, bcVelocityZ, cbc );

   createVelocityBoundaryConditionAndUID( "cmbX", typeCMB, flagCMB, bcVelocityX, cbc );
   createVelocityBoundaryConditionAndUID( "cmbY", typeCMB, flagCMB, bcVelocityY, cbc );
   createVelocityBoundaryConditionAndUID( "cmbZ", typeCMB, flagCMB, bcVelocityZ, cbc );

   // we want to fix the temperature at the surface and cmb
   // -> we always want Dirichlet boundaries for the temperature at the surface and cmb
   createTemperatureBoundaryConditionAndUID( "surface", DoFType::DirichletBoundary, flagSurface, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmb", DoFType::DirichletBoundary, flagCMB, bcTemperature, cbc );

   cbc.setBC( bcVelocityX.second, bcVelocityY.second, bcVelocityZ.second, bcTemperature );

   return cbc;
}

///////////////////////////////////////////////////////
/////////////////////// Annulus ///////////////////////
///////////////////////////////////////////////////////
bool AnnulusCMBBoundary( const Point3D& x, real_t rCMB, real_t boundaryTolerance )
{
   return x.norm() < rCMB + boundaryTolerance;
}
bool AnnulusSurfaceBoundary( const Point3D& x, real_t rSurface, real_t boundaryTolerance )
{
   return x.norm() > rSurface - boundaryTolerance;
}
void AnnulusBoundaryNormal( const Point3D& x, Point3D& normal, real_t rSurface, real_t rCMB )
{
   real_t radius = x.norm();
   if ( std::abs( radius - rSurface ) < std::abs( radius - rCMB ) )
   {
      normal = Point3D( { x[0] / radius, x[1] / radius, x[2] / radius } );
   }
   else
   {
      normal = Point3D( { -x[0] / radius, -x[1] / radius, -x[2] / radius } );
   }
}

std::function< bool( const Point3D& x ) > createAnnulusSurfaceBoundaryFct( real_t rSurface, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return AnnulusSurfaceBoundary( x, rSurface, boundaryTolerance ); };
}

std::function< bool( const Point3D& x ) > createAnnulusCMBBoundaryFct( real_t rCMB, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return AnnulusCMBBoundary( x, rCMB, boundaryTolerance ); };
}

std::function< void( const Point3D& in, Point3D& out ) > createAnnulusBoundaryNormalFct( real_t rSurface, real_t rCMB )
{
   return [=]( const Point3D& x, Point3D& normal ) { AnnulusBoundaryNormal( x, normal, rSurface, rCMB ); };
}

std::shared_ptr< hyteg::PrimitiveStorage > createAnnulusStorage( uint_t nTan,
                                                                 uint_t nRad,
                                                                 real_t rCMB,
                                                                 real_t rSurface,
                                                                 bool   blending,
                                                                 real_t boundaryTolerance,
                                                                 uint_t flagSurface,
                                                                 uint_t flagCMB,
                                                                 uint_t additionalHaloDepth )
{
   // create the annulus mesh
   MeshInfo              meshInfo = MeshInfo::meshAnnulus( rCMB, rSurface, MeshInfo::meshFlavour::CRISS, nTan, nRad );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( blending )
   {
      AnnulusMap::setMap( setupStorage );
   }

   // set load balancing
   loadbalancing::roundRobinVolume( setupStorage );

   // set surface and cmb boundary flags by vertex location
   setupStorage.setMeshBoundaryFlagsByVertexLocation( 0, createAllVerticesFct(), true );
   setupStorage.setMeshBoundaryFlagsByVertexLocation(
       flagSurface, createAnnulusSurfaceBoundaryFct( rSurface, boundaryTolerance ), true );
   setupStorage.setMeshBoundaryFlagsByVertexLocation( flagCMB, createAnnulusCMBBoundaryFct( rCMB, boundaryTolerance ), true );

   // create storage
   return std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
}

ConvectionBC createAnnulusBoundaryConditions( DoFType typeSurface, DoFType typeCMB, uint_t flagSurface, uint_t flagCMB )
{
   std::pair< uint_t, BoundaryCondition > bcVelocityX( 0, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityY( 1, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityZ( 2, BoundaryCondition() );
   BoundaryCondition                      bcTemperature;

   ConvectionBC cbc( ConvectionMeshType::Annulus );

   createVelocityBoundaryConditionAndUID( "surfaceX", typeSurface, flagSurface, bcVelocityX, cbc );
   createVelocityBoundaryConditionAndUID( "surfaceY", typeSurface, flagSurface, bcVelocityY, cbc );

   createVelocityBoundaryConditionAndUID( "cmbX", typeCMB, flagCMB, bcVelocityX, cbc );
   createVelocityBoundaryConditionAndUID( "cmbY", typeCMB, flagCMB, bcVelocityY, cbc );

   // we want to fix the temperature at the surface and cmb
   // -> we always want Dirichlet boundaries for the temperature at the surface and cmb
   createTemperatureBoundaryConditionAndUID( "surface", DoFType::DirichletBoundary, flagSurface, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmb", DoFType::DirichletBoundary, flagCMB, bcTemperature, cbc );

   cbc.setBC( bcVelocityX.second, bcVelocityY.second, bcVelocityZ.second, bcTemperature );

   return cbc;
}

///////////////////////////////////////////////////////
////////////////////// Rectangle //////////////////////
///////////////////////////////////////////////////////
bool RectangleCMBBoundary( const hyteg::Point3D& x, real_t yMin, real_t boundaryTolerance )
{
   return x[1] < yMin + boundaryTolerance;
}

bool RectangleSurfaceBoundary( const hyteg::Point3D& x, real_t yMax, real_t boundaryTolerance )
{
   return x[1] > yMax - boundaryTolerance;
}

bool RectangleLeftBoundary( const hyteg::Point3D& x, real_t xMin, real_t boundaryTolerance )
{
   return x[0] < xMin + boundaryTolerance;
}

bool RectangleRightBoundary( const hyteg::Point3D& x, real_t xMax, real_t boundaryTolerance )
{
   return x[0] > xMax - boundaryTolerance;
}

bool RectangleCornerLeftCMBBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMin, real_t boundaryTolerance )
{
   // (cmb) && (left)
   return RectangleCMBBoundary( x, yMin, boundaryTolerance ) && RectangleLeftBoundary( x, xMin, boundaryTolerance );
}

bool RectangleCornerRightCMBBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMin, real_t boundaryTolerance )
{
   // (cmb) && (right)
   return RectangleCMBBoundary( x, yMin, boundaryTolerance ) && RectangleRightBoundary( x, xMax, boundaryTolerance );
}

bool RectangleCornerRightSurfaceBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMax, real_t boundaryTolerance )
{
   // (surface) && (right)
   return RectangleSurfaceBoundary( x, yMax, boundaryTolerance ) && RectangleRightBoundary( x, xMax, boundaryTolerance );
}

bool RectangleCornerLeftSurfaceBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMax, real_t boundaryTolerance )
{
   // (surface) && (left)
   return RectangleSurfaceBoundary( x, yMax, boundaryTolerance ) && RectangleLeftBoundary( x, xMin, boundaryTolerance );
}

bool RectangleCornerBoundary( const hyteg::Point3D& x,
                              real_t                xMin,
                              real_t                xMax,
                              real_t                yMin,
                              real_t                yMax,
                              real_t                boundaryTolerance )
{
   // (surface OR cmb) && (left OR right)
   return ( RectangleSurfaceBoundary( x, yMax, boundaryTolerance ) || RectangleCMBBoundary( x, yMin, boundaryTolerance ) ) &&
          ( RectangleLeftBoundary( x, xMin, boundaryTolerance ) || RectangleRightBoundary( x, xMax, boundaryTolerance ) );
}

// In case of corners (i.e. normals with more than one component != 0) you have to project out every component of the returned normal individually
// If normal = (0,0,0), x is not on the boundary
//
void RectangleBoundaryNormal( const Point3D& x,
                              Point3D&       normal,
                              real_t         xMin,
                              real_t         xMax,
                              real_t         yMin,
                              real_t         yMax,
                              real_t         boundaryTolerance )
{
   normal = Point3D( { 0, 0, 0 } );

   if ( RectangleSurfaceBoundary( x, yMax, boundaryTolerance ) )
   {
      normal[1] = real_c( 1.0 );
   }
   else if ( RectangleCMBBoundary( x, yMin, boundaryTolerance ) )
   {
      normal[1] = real_c( -1.0 );
   }

   if ( RectangleRightBoundary( x, xMax, boundaryTolerance ) )
   {
      normal[0] = real_c( 1.0 );
   }
   else if ( RectangleLeftBoundary( x, xMin, boundaryTolerance ) )
   {
      normal[0] = real_c( -1.0 );
   }
}

std::function< bool( const hyteg::Point3D& x ) > createRectangleCMBBoundaryFct( real_t yMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleCMBBoundary( x, yMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createRectangleSurfaceBoundaryFct( real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleSurfaceBoundary( x, yMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createRectangleLeftBoundaryFct( real_t xMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleLeftBoundary( x, xMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createRectangleRightBoundaryFct( real_t xMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleRightBoundary( x, xMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerLeftCMBBoundaryFct( real_t xMin, real_t yMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleCornerLeftCMBBoundary( x, xMin, yMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerRightCMBBoundaryFct( real_t xMax, real_t yMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleCornerRightCMBBoundary( x, xMax, yMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerRightSurfaceBoundaryFct( real_t xMax, real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleCornerRightSurfaceBoundary( x, xMax, yMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerLeftSurfaceBoundaryFct( real_t xMin, real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleCornerLeftSurfaceBoundary( x, xMin, yMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerBoundaryFct( real_t xMin, real_t xMax, real_t yMin, real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return RectangleCornerBoundary( x, xMin, xMax, yMin, yMax, boundaryTolerance ); };
}

std::function< void( const Point3D& in, Point3D& out ) >
    createRectangleBoundaryNormalFct( real_t xMin, real_t xMax, real_t yMin, real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x, Point3D& normal ) {
      RectangleBoundaryNormal( x, normal, xMin, xMax, yMin, yMax, boundaryTolerance );
   };
}

std::shared_ptr< hyteg::PrimitiveStorage > createRectangleStorage( uint_t  nX,
                                                                   uint_t  nY,
                                                                   Point2D P1,
                                                                   Point2D P2,
                                                                   real_t  boundaryTolerance,
                                                                   uint_t  flagSurface,
                                                                   uint_t  flagCMB,
                                                                   uint_t  flagLeft,
                                                                   uint_t  flagRight,
                                                                   uint_t  flagCornerLeftCMB,
                                                                   uint_t  flagCornerRightCMB,
                                                                   uint_t  flagCornerRightSurface,
                                                                   uint_t  flagCornerLeftSurface,
                                                                   uint_t  additionalHaloDepth )
{
   // Init mesh and setup storage
   MeshInfo              meshInfo = MeshInfo::meshRectangle( P1, P2, MeshInfo::CRISSCROSS, nX, nY );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // set load balancing
   loadbalancing::roundRobinVolume( setupStorage );

   // set boundary flags by centroid location
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( 0, createAllVerticesFct(), true );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagSurface,
                                                        createRectangleSurfaceBoundaryFct( P2[1], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagCMB, createRectangleCMBBoundaryFct( P1[1], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagLeft, createRectangleLeftBoundaryFct( P1[0], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagRight, createRectangleRightBoundaryFct( P2[0], boundaryTolerance ) );

   setupStorage.setMeshBoundaryFlagsByCentroidLocation(
       flagCornerLeftCMB, createRectangleCornerLeftCMBBoundaryFct( P1[0], P1[1], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(
       flagCornerRightCMB, createRectangleCornerRightCMBBoundaryFct( P2[0], P1[1], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(
       flagCornerRightSurface, createRectangleCornerRightSurfaceBoundaryFct( P2[0], P2[1], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(
       flagCornerLeftSurface, createRectangleCornerLeftSurfaceBoundaryFct( P1[0], P2[1], boundaryTolerance ) );

   // create storage
   return std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
}

ConvectionBC createRectangleBoundaryConditions( DoFType typeSurface,
                                                DoFType typeCMB,
                                                DoFType typeLeft,
                                                DoFType typeRight,
                                                DoFType typeTemperatureSide,
                                                uint_t  flagSurface,
                                                uint_t  flagCMB,
                                                uint_t  flagLeft,
                                                uint_t  flagRight,
                                                uint_t  flagCornerLeftCMB,
                                                uint_t  flagCornerRightCMB,
                                                uint_t  flagCornerRightSurface,
                                                uint_t  flagCornerLeftSurface )
{
   std::pair< uint_t, BoundaryCondition > bcVelocityX( 0, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityY( 1, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityZ( 2, BoundaryCondition() );
   BoundaryCondition                      bcTemperature;

   ConvectionBC cbc( ConvectionMeshType::Rectangle );

   // clang-format off
   create2DEdgeVelocityBoundaryConditionAndUID( "surface", typeSurface, flagSurface, bcVelocityX, bcVelocityY, cbc, true);
   create2DEdgeVelocityBoundaryConditionAndUID( "cmb"    , typeCMB    , flagCMB    , bcVelocityX, bcVelocityY, cbc, true );

   create2DEdgeVelocityBoundaryConditionAndUID( "left" , typeLeft , flagLeft , bcVelocityX, bcVelocityY, cbc );
   create2DEdgeVelocityBoundaryConditionAndUID( "right", typeRight, flagRight, bcVelocityX, bcVelocityY, cbc );

   create2DCornerVelocityBoundaryConditionAndUID( "leftcmbcorner"     , typeLeft , typeCMB    , flagCornerLeftCMB     , bcVelocityX, bcVelocityY, cbc);
   create2DCornerVelocityBoundaryConditionAndUID( "rightcmbcorner"    , typeRight, typeCMB    , flagCornerRightCMB    , bcVelocityX, bcVelocityY, cbc);
   create2DCornerVelocityBoundaryConditionAndUID( "rightsurfacecorner", typeRight, typeSurface, flagCornerRightSurface, bcVelocityX, bcVelocityY, cbc);
   create2DCornerVelocityBoundaryConditionAndUID( "leftsurfacecorner" , typeLeft , typeSurface, flagCornerLeftSurface , bcVelocityX, bcVelocityY, cbc);
   

   // we want to fix the temperature at the surface and cmb (including corners)
   // -> we always want Dirichlet boundaries for the temperature at the surface and cmb
   // we use neumann boundary conditions for the temperature at the sides (reflective symmetry)
   createTemperatureBoundaryConditionAndUID( "surface"           , DoFType::DirichletBoundary, flagSurface           , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmb"               , DoFType::DirichletBoundary, flagCMB               , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "leftcmbcorner"     , DoFType::DirichletBoundary, flagCornerLeftCMB     , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "rightcmbcorner"    , DoFType::DirichletBoundary, flagCornerRightCMB    , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "rightsurfacecorner", DoFType::DirichletBoundary, flagCornerRightSurface, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "leftsurfacecorner" , DoFType::DirichletBoundary, flagCornerLeftSurface , bcTemperature, cbc );

   createTemperatureBoundaryConditionAndUID( "left" , typeTemperatureSide, flagLeft , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "right", typeTemperatureSide, flagRight, bcTemperature, cbc );
   // clang-format on

   cbc.setBC( bcVelocityX.second, bcVelocityY.second, bcVelocityZ.second, bcTemperature );

   return cbc;
}
//////////////////////////////////////////////////////
//////////////////////// Cuboid ////////////////////////
//////////////////////////////////////////////////////
bool CuboidCMBBoundary( const hyteg::Point3D& x, real_t zMin, real_t boundaryTolerance )
{
   return x[2] < zMin + boundaryTolerance;
}

bool CuboidSurfaceBoundary( const hyteg::Point3D& x, real_t zMax, real_t boundaryTolerance )
{
   return x[2] > zMax - boundaryTolerance;
}

bool CuboidLeftBoundary( const hyteg::Point3D& x, real_t xMin, real_t boundaryTolerance )
{
   return x[0] < xMin + boundaryTolerance;
}

bool CuboidRightBoundary( const hyteg::Point3D& x, real_t xMax, real_t boundaryTolerance )
{
   return x[0] > xMax - boundaryTolerance;
}

bool CuboidFrontBoundary( const hyteg::Point3D& x, real_t yMin, real_t boundaryTolerance )
{
   return x[1] < yMin + boundaryTolerance;
}

bool CuboidBackBoundary( const hyteg::Point3D& x, real_t yMax, real_t boundaryTolerance )
{
   return x[1] > yMax - boundaryTolerance;
}

bool CuboidLeftFrontEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMin, real_t boundaryTolerance )
{
   // (( left ) && ( front ))
   return ( ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) && ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) );
}

bool CuboidRightFrontEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMin, real_t boundaryTolerance )
{
   // (( right ) && ( front ))
   return ( ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) && ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) );
}

bool CuboidRightBackEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMax, real_t boundaryTolerance )
{
   // (( right ) && ( back ))
   return ( ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) && ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) );
}

bool CuboidLeftBackEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMax, real_t boundaryTolerance )
{
   // (( left ) && ( back ))
   return ( ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) && ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) );
}

bool CuboidLeftCMBEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t zMin, real_t boundaryTolerance )
{
   // (( cmb ) && ( left ))
   return ( ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) && ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) );
}

bool CuboidRightCMBEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t zMin, real_t boundaryTolerance )
{
   // (( cmb   ) && ( right ))
   return ( ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) && ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) );
}

bool CuboidRightSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t zMax, real_t boundaryTolerance )
{
   // (( surface ) && ( right ))
   return ( ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) && ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) );
}

bool CuboidLeftSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t zMax, real_t boundaryTolerance )
{
   // (( surface  ) && ( left ))
   return ( ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) && ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) );
}

bool CuboidFrontCMBEdgeBoundary( const hyteg::Point3D& x, real_t yMin, real_t zMin, real_t boundaryTolerance )
{
   // (( front  ) && ( cmb ))
   return ( ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) && ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) );
}

bool CuboidBackCMBEdgeBoundary( const hyteg::Point3D& x, real_t yMax, real_t zMin, real_t boundaryTolerance )
{
   // (( back  ) && ( cmb ))
   return ( ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) && ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) );
}

bool CuboidBackSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t yMax, real_t zMax, real_t boundaryTolerance )
{
   // (( back  ) && ( surface ))
   return ( ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) && ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) );
}

bool CuboidFrontSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t yMin, real_t zMax, real_t boundaryTolerance )
{
   // (( front ) && ( surface ))
   return ( ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) && ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) );
}

bool CuboidEdgeBoundary( const hyteg::Point3D& x,
                         real_t                xMin,
                         real_t                xMax,
                         real_t                yMin,
                         real_t                yMax,
                         real_t                zMin,
                         real_t                zMax,
                         real_t                boundaryTolerance )
{
   // (( surface || cmb   ) && ( left    || right )) ||
   // (( left    || right ) && ( front   || back  )) ||
   // (( front   || back  ) && ( surface || cmb   )) ||
   return ( ( ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) || ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) ) &&
            ( ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) || ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) ) ) ||
          ( ( ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) || ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) ) &&
            ( ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) || ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) ) ) ||
          ( ( ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) || ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) ) &&
            ( ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) || ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) ) );
}

bool CuboidSurfaceLeftFrontCornerBoundary( const hyteg::Point3D& x,
                                           real_t                xMin,
                                           real_t                yMin,
                                           real_t                zMax,
                                           real_t                boundaryTolerance )
{
   //(surface) && (left) && (front)
   return ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) && ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) &&
          ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) );
}

bool CuboidSurfaceLeftBackCornerBoundary( const hyteg::Point3D& x,
                                          real_t                xMin,
                                          real_t                yMax,
                                          real_t                zMax,
                                          real_t                boundaryTolerance )
{
   //(surface) && (left) && (back)
   return ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) && ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) &&
          ( CuboidBackBoundary( x, yMax, boundaryTolerance ) );
}

bool CuboidSurfaceRightFrontCornerBoundary( const hyteg::Point3D& x,
                                            real_t                xMax,
                                            real_t                yMin,
                                            real_t                zMax,
                                            real_t                boundaryTolerance )
{
   //(surface) && (right) && (front)
   return ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) && ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) &&
          ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) );
}

bool CuboidSurfaceRightBackCornerBoundary( const hyteg::Point3D& x,
                                           real_t                xMax,
                                           real_t                yMax,
                                           real_t                zMax,
                                           real_t                boundaryTolerance )
{
   //(surface) && (right) && (back)
   return ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) && ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) &&
          ( CuboidBackBoundary( x, yMax, boundaryTolerance ) );
}

bool CuboidCMBLeftFrontCornerBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMin, real_t zMin, real_t boundaryTolerance )
{
   //(cmb) && (left) && (front)
   return ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) && ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) &&
          ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) );
}

bool CuboidCMBLeftBackCornerBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMax, real_t zMin, real_t boundaryTolerance )
{
   //(cmb) && (left) && (back)
   return ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) && ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) &&
          ( CuboidBackBoundary( x, yMax, boundaryTolerance ) );
}

bool CuboidCMBRightFrontCornerBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMin, real_t zMin, real_t boundaryTolerance )
{
   //(cmb) && (right) && (front)
   return ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) && ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) &&
          ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) );
}

bool CuboidCMBRightBackCornerBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMax, real_t zMin, real_t boundaryTolerance )
{
   //(cmb) && (right) && (back)
   return ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) && ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) &&
          ( CuboidBackBoundary( x, yMax, boundaryTolerance ) );
}

bool CuboidCornerBoundary( const hyteg::Point3D& x,
                           real_t                xMin,
                           real_t                xMax,
                           real_t                yMin,
                           real_t                yMax,
                           real_t                zMin,
                           real_t                zMax,
                           real_t                boundaryTolerance )
{
   //(surface || cmb) && (left || right) && (front || back)
   return ( ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) ) || ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) ) ) &&
          ( ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) ) || ( CuboidRightBoundary( x, xMax, boundaryTolerance ) ) ) &&
          ( ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) ) || ( CuboidBackBoundary( x, yMax, boundaryTolerance ) ) );
}

// In case of corners and edges (i.e. normals with more than one component != 0) you have to project out every component of the returned normal individually
// If normal = (0,0,0), x is not on the boundary
void CuboidBoundaryNormal( const Point3D& x,
                           Point3D&       normal,
                           real_t         xMin,
                           real_t         xMax,
                           real_t         yMin,
                           real_t         yMax,
                           real_t         zMin,
                           real_t         zMax,
                           real_t         boundaryTolerance )
{
   normal = Point3D( { 0, 0, 0 } );

   if ( CuboidSurfaceBoundary( x, zMax, boundaryTolerance ) )
   {
      normal[2] = real_c( 1.0 );
   }
   else if ( CuboidCMBBoundary( x, zMin, boundaryTolerance ) )
   {
      normal[2] = real_c( -1.0 );
   }

   if ( CuboidLeftBoundary( x, xMin, boundaryTolerance ) )
   {
      normal[0] = real_c( -1.0 );
   }
   else if ( CuboidRightBoundary( x, xMax, boundaryTolerance ) )
   {
      normal[0] = real_c( 1.0 );
   }

   if ( CuboidFrontBoundary( x, yMin, boundaryTolerance ) )
   {
      normal[1] = real_c( -1.0 );
   }
   else if ( CuboidBackBoundary( x, yMax, boundaryTolerance ) )
   {
      normal[1] = real_c( 1.0 );
   }
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidCMBBoundaryFct( real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidCMBBoundary( x, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidSurfaceBoundaryFct( real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidSurfaceBoundary( x, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidLeftBoundaryFct( real_t xMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidLeftBoundary( x, xMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidRightBoundaryFct( real_t xMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidRightBoundary( x, xMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidFrontBoundaryFct( real_t yMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidFrontBoundary( x, yMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidBackBoundaryFct( real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidBackBoundary( x, yMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftFrontEdgeBoundaryFct( real_t xMin, real_t yMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidLeftFrontEdgeBoundary( x, xMin, yMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightFrontEdgeBoundaryFct( real_t xMax, real_t yMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidRightFrontEdgeBoundary( x, xMax, yMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightBackEdgeBoundaryFct( real_t xMax, real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidRightBackEdgeBoundary( x, xMax, yMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftBackEdgeBoundaryFct( real_t xMin, real_t yMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidLeftBackEdgeBoundary( x, xMin, yMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftCMBEdgeBoundaryFct( real_t xMin, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidLeftCMBEdgeBoundary( x, xMin, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightCMBEdgeBoundaryFct( real_t xMax, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidRightCMBEdgeBoundary( x, xMax, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightSurfaceEdgeBoundaryFct( real_t xMax, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidRightSurfaceEdgeBoundary( x, xMax, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftSurfaceEdgeBoundaryFct( real_t xMin, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidLeftSurfaceEdgeBoundary( x, xMin, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidFrontCMBEdgeBoundaryFct( real_t yMin, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidFrontCMBEdgeBoundary( x, yMin, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidBackCMBEdgeBoundaryFct( real_t yMax, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidBackCMBEdgeBoundary( x, yMax, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidBackSurfaceEdgeBoundaryFct( real_t yMax, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidBackSurfaceEdgeBoundary( x, yMax, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidFrontSurfaceEdgeBoundaryFct( real_t yMin, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidFrontSurfaceEdgeBoundary( x, yMin, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidEdgeBoundaryFct( real_t xMin,
                                                                              real_t xMax,
                                                                              real_t yMin,
                                                                              real_t yMax,
                                                                              real_t zMin,
                                                                              real_t zMax,
                                                                              real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidEdgeBoundary( x, xMin, xMax, yMin, yMax, zMin, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceLeftFrontCornerBoundaryFct( real_t xMin, real_t yMin, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidSurfaceLeftFrontCornerBoundary( x, xMin, yMin, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceLeftBackCornerBoundaryFct( real_t xMin, real_t yMax, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidSurfaceLeftBackCornerBoundary( x, xMin, yMax, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceRightFrontCornerBoundaryFct( real_t xMax, real_t yMin, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidSurfaceRightFrontCornerBoundary( x, xMax, yMin, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceRightBackCornerBoundaryFct( real_t xMax, real_t yMax, real_t zMax, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidSurfaceRightBackCornerBoundary( x, xMax, yMax, zMax, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBLeftFrontCornerBoundaryFct( real_t xMin, real_t yMin, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidCMBLeftFrontCornerBoundary( x, xMin, yMin, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBLeftBackCornerBoundaryFct( real_t xMin, real_t yMax, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidCMBLeftBackCornerBoundary( x, xMin, yMax, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBRightFrontCornerBoundaryFct( real_t xMax, real_t yMin, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidCMBRightFrontCornerBoundary( x, xMax, yMin, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBRightBackCornerBoundaryFct( real_t xMax, real_t yMax, real_t zMin, real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidCMBRightBackCornerBoundary( x, xMax, yMax, zMin, boundaryTolerance ); };
}

std::function< bool( const hyteg::Point3D& x ) > createCuboidCornerBoundaryFct( real_t xMin,
                                                                                real_t xMax,
                                                                                real_t yMin,
                                                                                real_t yMax,
                                                                                real_t zMin,
                                                                                real_t zMax,
                                                                                real_t boundaryTolerance )
{
   return [=]( const Point3D& x ) { return CuboidCornerBoundary( x, xMin, xMax, yMin, yMax, zMin, zMax, boundaryTolerance ); };
}

std::function< void( const Point3D& in, Point3D& out ) > createCuboidBoundaryNormalFct( real_t xMin,
                                                                                        real_t xMax,
                                                                                        real_t yMin,
                                                                                        real_t yMax,
                                                                                        real_t zMin,
                                                                                        real_t zMax,
                                                                                        real_t boundaryTolerance )
{
   return [=]( const Point3D& x, Point3D& normal ) {
      CuboidBoundaryNormal( x, normal, xMin, xMax, yMin, yMax, zMin, zMax, boundaryTolerance );
   };
}

std::shared_ptr< hyteg::PrimitiveStorage > createCuboidStorage( uint_t  nX,
                                                                uint_t  nY,
                                                                uint_t  nZ,
                                                                Point3D P1,
                                                                Point3D P2,
                                                                real_t  boundaryTolerance,
                                                                uint_t  flagSurface,
                                                                uint_t  flagCMB,
                                                                uint_t  flagLeft,
                                                                uint_t  flagRight,
                                                                uint_t  flagFront,
                                                                uint_t  flagBack,

                                                                uint_t flagLeftFrontEdge,
                                                                uint_t flagRightFrontEdge,
                                                                uint_t flagRightBackEdge,
                                                                uint_t flagLeftBackEdge,
                                                                uint_t flagLeftCMBEdge,
                                                                uint_t flagRightCMBEdge,
                                                                uint_t flagRightSurfaceEdge,
                                                                uint_t flagLeftSurfaceEdge,
                                                                uint_t flagFrontCMBEdge,
                                                                uint_t flagBackCMBEdge,
                                                                uint_t flagBackSurfaceEdge,
                                                                uint_t flagFrontSurfaceEdge,

                                                                uint_t flagSurfaceLeftFrontCorner,
                                                                uint_t flagSurfaceLeftBackCorner,
                                                                uint_t flagSurfaceRightFrontCorner,
                                                                uint_t flagSurfaceRightBackCorner,
                                                                uint_t flagCMBLeftFrontCorner,
                                                                uint_t flagCMBLeftBackCorner,
                                                                uint_t flagCMBRightFrontCorner,
                                                                uint_t flagCMBRightBackCorner,
                                                                uint_t additionalHaloDepth )
{
   // Init mesh and setup storage
   // Symmetric cuboid is important, so Taylor-Hood fulfills the inf-sup condition
   MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( P1, P2, nX, nY, nZ );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // set load balancing
   loadbalancing::roundRobinVolume( setupStorage );

   // set boundary flags by centroid location
   // clang-format off
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( 0, createAllVerticesFct(), true );

   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagSurface, createCuboidSurfaceBoundaryFct       ( P2[2], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagCMB,     createCuboidCMBBoundaryFct           ( P1[2], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagLeft,    createCuboidLeftBoundaryFct  ( P1[0], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagRight,   createCuboidRightBoundaryFct ( P2[0], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagFront,   createCuboidFrontBoundaryFct ( P1[1], boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( flagBack,    createCuboidBackBoundaryFct  ( P2[1], boundaryTolerance ) );

   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagLeftFrontEdge,    createCuboidLeftFrontEdgeBoundaryFct    ( P1[0], P1[1] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagRightFrontEdge,   createCuboidRightFrontEdgeBoundaryFct   ( P2[0], P1[1] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagRightBackEdge,    createCuboidRightBackEdgeBoundaryFct    ( P2[0], P2[1] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagLeftBackEdge,     createCuboidLeftBackEdgeBoundaryFct     ( P1[0], P2[1] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagLeftCMBEdge,      createCuboidLeftCMBEdgeBoundaryFct      ( P1[0], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagRightCMBEdge,     createCuboidRightCMBEdgeBoundaryFct     ( P2[0], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagRightSurfaceEdge, createCuboidRightSurfaceEdgeBoundaryFct ( P2[0], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagLeftSurfaceEdge,  createCuboidLeftSurfaceEdgeBoundaryFct  ( P1[0], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagFrontCMBEdge,     createCuboidFrontCMBEdgeBoundaryFct     ( P1[1], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagBackCMBEdge,      createCuboidBackCMBEdgeBoundaryFct      ( P2[1], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagBackSurfaceEdge,  createCuboidBackSurfaceEdgeBoundaryFct  ( P2[1], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagFrontSurfaceEdge, createCuboidFrontSurfaceEdgeBoundaryFct ( P1[1], P2[2] , boundaryTolerance ) );

   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagSurfaceLeftFrontCorner,  createCuboidSurfaceLeftFrontCornerBoundaryFct  (P1[0], P1[1], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagSurfaceLeftBackCorner,   createCuboidSurfaceLeftBackCornerBoundaryFct   (P1[0], P2[1], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagSurfaceRightFrontCorner, createCuboidSurfaceRightFrontCornerBoundaryFct (P2[0], P1[1], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagSurfaceRightBackCorner,  createCuboidSurfaceRightBackCornerBoundaryFct  (P2[0], P2[1], P2[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagCMBLeftFrontCorner,      createCuboidCMBLeftFrontCornerBoundaryFct      (P1[0], P1[1], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagCMBLeftBackCorner,       createCuboidCMBLeftBackCornerBoundaryFct       (P1[0], P2[1], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagCMBRightFrontCorner,     createCuboidCMBRightFrontCornerBoundaryFct     (P2[0], P1[1], P1[2] , boundaryTolerance ) );
   setupStorage.setMeshBoundaryFlagsByCentroidLocation(flagCMBRightBackCorner,      createCuboidCMBRightBackCornerBoundaryFct      (P2[0], P2[1], P1[2] , boundaryTolerance ) );
   // clang-format on

   // create storage
   return std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );
}

ConvectionBC createCuboidBoundaryConditions( DoFType typeSurface,
                                             DoFType typeCMB,
                                             DoFType typeLeft,
                                             DoFType typeRight,
                                             DoFType typeFront,
                                             DoFType typeBack,
                                             DoFType typeTemperatureSide,
                                             uint_t  flagSurface,
                                             uint_t  flagCMB,
                                             uint_t  flagLeft,
                                             uint_t  flagRight,
                                             uint_t  flagFront,
                                             uint_t  flagBack,

                                             uint_t flagLeftFrontEdge,
                                             uint_t flagRightFrontEdge,
                                             uint_t flagRightBackEdge,
                                             uint_t flagLeftBackEdge,
                                             uint_t flagLeftCMBEdge,
                                             uint_t flagRightCMBEdge,
                                             uint_t flagRightSurfaceEdge,
                                             uint_t flagLeftSurfaceEdge,
                                             uint_t flagFrontCMBEdge,
                                             uint_t flagBackCMBEdge,
                                             uint_t flagBackSurfaceEdge,
                                             uint_t flagFrontSurfaceEdge,

                                             uint_t flagSurfaceLeftFrontCorner,
                                             uint_t flagSurfaceLeftBackCorner,
                                             uint_t flagSurfaceRightFrontCorner,
                                             uint_t flagSurfaceRightBackCorner,
                                             uint_t flagCMBLeftFrontCorner,
                                             uint_t flagCMBLeftBackCorner,
                                             uint_t flagCMBRightFrontCorner,
                                             uint_t flagCMBRightBackCorner )
{
   std::pair< uint_t, BoundaryCondition > bcVelocityX( 0, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityY( 1, BoundaryCondition() );
   std::pair< uint_t, BoundaryCondition > bcVelocityZ( 2, BoundaryCondition() );
   BoundaryCondition                      bcTemperature;

   ConvectionBC cbc( ConvectionMeshType::Cuboid );

   // clang-format off

   // faces
   create3DFaceVelocityBoundaryConditionAndUID( "surface", typeSurface, flagSurface, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,2,1);
   create3DFaceVelocityBoundaryConditionAndUID( "cmb"    , typeCMB    , flagCMB    , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,2,1);

   create3DFaceVelocityBoundaryConditionAndUID( "left" , typeLeft , flagLeft , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1,2,3);
   create3DFaceVelocityBoundaryConditionAndUID( "right", typeRight, flagRight, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1,2,3);
   
   create3DFaceVelocityBoundaryConditionAndUID( "front", typeFront, flagFront, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 2,1,3);
   create3DFaceVelocityBoundaryConditionAndUID( "back" , typeBack , flagBack , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 2,1,3);

   // edges
   create3DEdgeVelocityBoundaryConditionAndUID( "leftfrontedge" , typeLeft , typeFront, flagLeftFrontEdge , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 2, 3 );
   create3DEdgeVelocityBoundaryConditionAndUID( "rightfrontedge", typeRight, typeFront, flagRightFrontEdge, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 2, 3 );
   create3DEdgeVelocityBoundaryConditionAndUID( "rightbackedge" , typeRight, typeBack , flagRightBackEdge , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 2, 3 );
   create3DEdgeVelocityBoundaryConditionAndUID( "leftbackedge"  , typeLeft , typeBack , flagLeftBackEdge  , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 2, 3 );
  
   create3DEdgeVelocityBoundaryConditionAndUID( "leftcmbedge"     , typeLeft , typeCMB    , flagLeftCMBEdge     , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 3, 2 );
   create3DEdgeVelocityBoundaryConditionAndUID( "rightcmbedge"    , typeRight, typeCMB    , flagRightCMBEdge    , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 3, 2 );
   create3DEdgeVelocityBoundaryConditionAndUID( "rightsurfaceedge", typeRight, typeSurface, flagRightSurfaceEdge, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 3, 2 );
   create3DEdgeVelocityBoundaryConditionAndUID( "leftsurfaceedge" , typeLeft , typeSurface, flagLeftSurfaceEdge , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 1, 3, 2 );

   create3DEdgeVelocityBoundaryConditionAndUID( "frontcmbedge"    , typeFront, typeCMB    , flagFrontCMBEdge    , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 2, 3, 1 );
   create3DEdgeVelocityBoundaryConditionAndUID( "backcmbedge"     , typeBack , typeCMB    , flagBackCMBEdge     , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 2, 3, 1 );
   create3DEdgeVelocityBoundaryConditionAndUID( "backsurfaceedge" , typeBack , typeSurface, flagBackSurfaceEdge , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 2, 3, 1 );
   create3DEdgeVelocityBoundaryConditionAndUID( "frontsurfaceedge", typeFront, typeSurface, flagFrontSurfaceEdge, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 2, 3, 1 );
   
   // corners
   create3DCornerVelocityBoundaryConditionAndUID( "surfaceleftfrontcorner" , typeSurface, typeLeft , typeFront, flagSurfaceLeftFrontCorner , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);
   create3DCornerVelocityBoundaryConditionAndUID( "surfaceleftbackcorner"  , typeSurface, typeLeft , typeBack , flagSurfaceLeftBackCorner  , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);
   create3DCornerVelocityBoundaryConditionAndUID( "surfacerightfrontcorner", typeSurface, typeRight, typeFront, flagSurfaceRightFrontCorner, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);
   create3DCornerVelocityBoundaryConditionAndUID( "surfacerightbackcorner" , typeSurface, typeRight, typeBack , flagSurfaceRightBackCorner , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);

   create3DCornerVelocityBoundaryConditionAndUID( "cmbleftfrontcorner" , typeCMB, typeLeft , typeFront, flagCMBLeftFrontCorner , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);
   create3DCornerVelocityBoundaryConditionAndUID( "cmbleftbackcorner"  , typeCMB, typeLeft , typeBack , flagCMBLeftBackCorner  , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);
   create3DCornerVelocityBoundaryConditionAndUID( "cmbrightfrontcorner", typeCMB, typeRight, typeFront, flagCMBRightFrontCorner, bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);
   create3DCornerVelocityBoundaryConditionAndUID( "cmbrightbackcorner" , typeCMB, typeRight, typeBack , flagCMBRightBackCorner , bcVelocityX, bcVelocityY, bcVelocityZ, cbc, 3,1,2);

   // we want to fix the temperature at the surface and cmb (including corners and edges)
   // -> we always want Dirichlet boundaries for the temperature at the surface and cmb
   // we use neumann boundary conditions for the temperature at the sides (reflective symmetry)
   createTemperatureBoundaryConditionAndUID("surface", DoFType::DirichletBoundary, flagSurface, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID("cmb",     DoFType::DirichletBoundary, flagCMB, bcTemperature, cbc );

   createTemperatureBoundaryConditionAndUID( "leftcmbedge"     , DoFType::DirichletBoundary, flagLeftCMBEdge     , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "rightcmbedge"    , DoFType::DirichletBoundary, flagRightCMBEdge    , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "rightsurfaceedge", DoFType::DirichletBoundary, flagRightSurfaceEdge, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "leftsurfaceedge" , DoFType::DirichletBoundary, flagLeftSurfaceEdge , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "frontcmbedge"    , DoFType::DirichletBoundary, flagFrontCMBEdge    , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "backcmbedge"     , DoFType::DirichletBoundary, flagBackCMBEdge     , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "backsurfaceedge" , DoFType::DirichletBoundary, flagBackSurfaceEdge , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "frontsurfaceedge", DoFType::DirichletBoundary, flagFrontSurfaceEdge, bcTemperature, cbc );

   createTemperatureBoundaryConditionAndUID( "surfaceleftfrontcorner" , DoFType::DirichletBoundary, flagSurfaceLeftFrontCorner , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "surfaceleftbackcorner"  , DoFType::DirichletBoundary, flagSurfaceLeftBackCorner  , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "surfacerightfrontcorner", DoFType::DirichletBoundary, flagSurfaceRightFrontCorner, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "surfacerightbackcorner" , DoFType::DirichletBoundary, flagSurfaceRightBackCorner , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmbleftfrontcorner"     , DoFType::DirichletBoundary, flagCMBLeftFrontCorner     , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmbleftbackcorner"      , DoFType::DirichletBoundary, flagCMBLeftBackCorner      , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmbrightfrontcorner"    , DoFType::DirichletBoundary, flagCMBRightFrontCorner    , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "cmbrightbackcorner"     , DoFType::DirichletBoundary, flagCMBRightBackCorner     , bcTemperature, cbc );

   createTemperatureBoundaryConditionAndUID( "left",  typeTemperatureSide, flagLeft , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "right", typeTemperatureSide, flagRight, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "front", typeTemperatureSide, flagFront, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "back",  typeTemperatureSide, flagBack , bcTemperature, cbc );

   createTemperatureBoundaryConditionAndUID( "leftfrontedge" , typeTemperatureSide, flagLeftFrontEdge , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "rightfrontedge", typeTemperatureSide, flagRightFrontEdge, bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "rightbackedge" , typeTemperatureSide, flagRightBackEdge , bcTemperature, cbc );
   createTemperatureBoundaryConditionAndUID( "leftbackedge"  , typeTemperatureSide, flagLeftBackEdge  , bcTemperature, cbc );

   // clang-format on

   cbc.setBC( bcVelocityX.second, bcVelocityY.second, bcVelocityZ.second, bcTemperature );

   return cbc;
}

} // namespace convectionToolbox
} // namespace hyteg
