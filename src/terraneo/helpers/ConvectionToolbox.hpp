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

#include <cmath>
#include <iostream>
#include <map>

#include "core/Environment.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {


/// @brief The convection toolbox is a wrapper that allows the easy creation of earth mantle convection type storages on the rectangle, annulus, cuboid, spherical shell.
/// It also provides tools to specifically interpolate only some surfaces (targetable via UID name) and handles boundary conditions on the edges and corners of a rectangle/cuboid correctly.
/// The toolbox uses its own type of ConvectionBC boundary condtion and provides tools to create functions with them.
///
/// This was designed for the MantleConvection app. Please don't change this without checking if it breaks the app.
namespace convectionToolbox {

///////////////////////////////////////////////////////
/////////////////////// General ///////////////////////
///////////////////////////////////////////////////////
enum ConvectionMeshType
{
   SphericalShell,
   Annulus,
   Rectangle,
   Cuboid
};

struct ConvectionBC
{
 public:
   ConvectionBC(){};
   ConvectionBC( ConvectionMeshType meshType );

   ConvectionBC( ConvectionMeshType meshType,
                 BoundaryCondition  bcVelocityX,
                 BoundaryCondition  bcVelocityY,
                 BoundaryCondition  bcVelocityZ,
                 BoundaryCondition  bcTemperature );

   friend std::ostream& operator<<( std::ostream& stream, const ConvectionBC& cbc );

   std::vector< BoundaryUID > getUIDsWithName( std::string name, std::set< uint_t > types );

   // 0 = X, 1 = Y, 2 = Z, 3 = T
   template < class Interpolant_T, class InterpolationFct >
   void interpolateUIDsWithName( InterpolationFct           inter,
                                 uint_t                     level,
                                 Interpolant_T              f,
                                 std::vector< std::string > names,
                                 std::set< uint_t >         types,
                                 DoFType                    flag = hyteg::All )
   {
      for ( std::string name : names )
      {
         std::transform( name.begin(), name.end(), name.begin(), []( unsigned char c ) { return std::tolower( c ); } );

         for ( const auto& [i, v] : BoundaryUIDs_ )
         {
            if ( ( DoFTypes_.find( i ) != DoFTypes_.end() ) && testFlag( DoFTypes_.at( i ), flag ) &&
                 ( types.find( v.first ) != types.end() ) && ( i.find( name ) != std::string::npos ) )
            {
               f.interpolate( { inter }, level, v.second );
            }
         }
      }
   }

   // 0 = X, 1 = Y, 2 = Z, 3 = T
   template < class Interpolant_T, class InterpolationFct >
   void interpolateUIDsWithAllNames( InterpolationFct                          inter,
                                     uint_t                                    level,
                                     Interpolant_T                             f,
                                     std::vector< std::vector< std::string > > names,
                                     std::set< uint_t >                        types,
                                     DoFType                                   flag = hyteg::All )
   {
      for ( const auto& [i, v] : BoundaryUIDs_ )
      {
         for ( std::vector< std::string > sublist : names )
         {
            bool matchName = true;
            for ( std::string name : sublist )
            {
               std::transform( name.begin(), name.end(), name.begin(), []( unsigned char c ) { return std::tolower( c ); } );

               if ( i.find( name ) == std::string::npos )
               {
                  matchName = false;
                  break;
               }
            }

            if ( matchName && ( DoFTypes_.find( i ) != DoFTypes_.end() ) && testFlag( DoFTypes_.at( i ), flag ) &&
                 ( types.find( v.first ) != types.end() ) )
            {
               f.interpolate( { inter }, level, v.second );
            }
         }
      }
   }

   std::pair< P2P1TaylorHoodFunction< real_t >, P2Function< real_t > >
       createVisualisation( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level );

   void setBC( BoundaryCondition bcVelocityX,
               BoundaryCondition bcVelocityY,
               BoundaryCondition bcVelocityZ,
               BoundaryCondition bcTemperature );

 public:
   ConvectionMeshType meshType_;

   BoundaryCondition bcVelocityX_;
   BoundaryCondition bcVelocityY_;
   BoundaryCondition bcVelocityZ_;
   BoundaryCondition bcTemperature_;

   std::map< std::string, DoFType >                          DoFTypes_;
   std::map< std::string, std::pair< uint_t, BoundaryUID > > BoundaryUIDs_;
};

// expects p to already lie in the dual space, e.g. p = M * p_tilde
real_t projectPressureZeroMean( const P1Function< real_t >& p,
                                const P1Function< real_t >& One,
                                uint_t                      level,
                                real_t                      VolumeOmega = 1.0,
                                DoFType                     flag        = All );

void createVelocityBoundaryConditionAndUID( std::string                             name,
                                            DoFType                                 t,
                                            uint_t                                  flag,
                                            std::pair< uint_t, BoundaryCondition >& bc,
                                            ConvectionBC&                           cbc );
void create2DEdgeVelocityBoundaryConditionAndUID( std::string                             name,
                                                  DoFType                                 t,
                                                  uint_t                                  flag,
                                                  std::pair< uint_t, BoundaryCondition >& bc1,
                                                  std::pair< uint_t, BoundaryCondition >& bc2,
                                                  ConvectionBC&                           cbc,
                                                  bool                                    invert = false );
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
                                                  uint_t                                  axis3 );
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
                                                    uint_t                                  axis3 );
void createTemperatureBoundaryConditionAndUID( std::string        name,
                                               DoFType            t,
                                               uint_t             flag,
                                               BoundaryCondition& bc,
                                               ConvectionBC&      cbc );

P2P1TaylorHoodFunction< real_t > createP2P1TaylorHoodFunction( const std::string&                         name,
                                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                                               size_t                                     minLevel,
                                                               size_t                                     maxLevel,
                                                               ConvectionBC&                              cbc );

P1StokesFunction< real_t > createP1StokesFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   size_t                                     minLevel,
                                                   size_t                                     maxLevel,
                                                   ConvectionBC&                              cbc );

P1VectorFunction< real_t > createP1VectorVelocityFunction( const std::string&                         name,
                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                           size_t                                     minLevel,
                                                           size_t                                     maxLevel,
                                                           ConvectionBC&                              cbc );
P2VectorFunction< real_t > createP2VectorVelocityFunction( const std::string&                         name,
                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                           size_t                                     minLevel,
                                                           size_t                                     maxLevel,
                                                           ConvectionBC&                              cbc );
P1Function< real_t >       createP1TemperatureFunction( const std::string&                         name,
                                                        const std::shared_ptr< PrimitiveStorage >& storage,
                                                        size_t                                     minLevel,
                                                        size_t                                     maxLevel,
                                                        ConvectionBC&                              cbc );
P2Function< real_t >       createP2TemperatureFunction( const std::string&                         name,
                                                        const std::shared_ptr< PrimitiveStorage >& storage,
                                                        size_t                                     minLevel,
                                                        size_t                                     maxLevel,
                                                        ConvectionBC&                              cbc );

std::shared_ptr< P2P1TaylorHoodFunction< real_t > >
    createSharedP2P1TaylorHoodFunction( const std::string&                         name,
                                        const std::shared_ptr< PrimitiveStorage >& storage,
                                        size_t                                     minLevel,
                                        size_t                                     maxLevel,
                                        ConvectionBC&                              cbc );
std::shared_ptr< P1VectorFunction< real_t > >
    createSharedP1VectorVelocityFunction( const std::string&                         name,
                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                          size_t                                     minLevel,
                                          size_t                                     maxLevel,
                                          ConvectionBC&                              cbc );
std::shared_ptr< P2VectorFunction< real_t > >
                                        createSharedP2VectorVelocityFunction( const std::string&                         name,
                                                                              const std::shared_ptr< PrimitiveStorage >& storage,
                                                                              size_t                                     minLevel,
                                                                              size_t                                     maxLevel,
                                                                              ConvectionBC&                              cbc );
std::shared_ptr< P1Function< real_t > > createSharedP1TemperatureFunction( const std::string&                         name,
                                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           size_t                                     minLevel,
                                                                           size_t                                     maxLevel,
                                                                           ConvectionBC&                              cbc );
std::shared_ptr< P2Function< real_t > > createSharedP2TemperatureFunction( const std::string&                         name,
                                                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           size_t                                     minLevel,
                                                                           size_t                                     maxLevel,
                                                                           ConvectionBC&                              cbc );

std::function< bool( const Point3D& x ) > createAllVerticesFct();

std::function< real_t( const Point3D& x, real_t temp ) >
    createSpaceDependentViscosityProfileWithJumps( real_t rSurface, real_t rCMB, real_t eta0 );

///////////////////////////////////////////////////////
/////////////////// Spherical Shell ///////////////////
///////////////////////////////////////////////////////
bool SphericalShellCMBBoundary( const Point3D& x, real_t rCMB, real_t boundaryTolerance );
bool SphericalShellSurfaceBoundary( const Point3D& x, real_t rSurface, real_t boundaryTolerance );
void SphericalShellBoundaryNormal( const Point3D& x, Point3D& normal, real_t rSurface, real_t rCMB );

std::function< bool( const hyteg::Point3D& x ) > createSphericalShellCMBBoundaryFct( real_t rCMB, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createSphericalShellSurfaceBoundaryFct( real_t rSurface,
                                                                                         real_t boundaryTolerance );
std::function< void( const Point3D& in, Point3D& out ) > createSphericalShellBoundaryNormalFct( real_t rSurface, real_t rCMB );

std::shared_ptr< hyteg::PrimitiveStorage > createSphericalShellStorage( uint_t nTan,
                                                                        uint_t nRad,
                                                                        real_t rCMB                = 1.0,
                                                                        real_t rSurface            = 2.0,
                                                                        bool   blending            = true,
                                                                        real_t boundaryTolerance   = 1e-10,
                                                                        uint_t flagSurface         = 1,
                                                                        uint_t flagCMB             = 2,
                                                                        uint_t additionalHaloDepth = 1 );

ConvectionBC createSphericalShellBoundaryConditions( DoFType typeSurface = DoFType::DirichletBoundary,
                                                     DoFType typeCMB     = DoFType::FreeslipBoundary,
                                                     uint_t  flagSurface = 1,
                                                     uint_t  flagCMB     = 2 );

///////////////////////////////////////////////////////
/////////////////////// Annulus ///////////////////////
///////////////////////////////////////////////////////
bool AnnulusCMBBoundary( const Point3D& x, real_t rCMB, real_t boundaryTolerance );
bool AnnulusSurfaceBoundary( const Point3D& x, real_t rSurface, real_t boundaryTolerance );
void AnnulusBoundaryNormal( const Point3D& x, Point3D& normal, real_t rSurface, real_t rCMB );

std::function< bool( const hyteg::Point3D& x ) > createAnnulusCMBBoundaryFct( real_t rCMB, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createAnnulusSurfaceBoundaryFct( real_t rSurface, real_t boundaryTolerance );
std::function< void( const Point3D& in, Point3D& out ) > createAnnulusBoundaryNormalFct( real_t rSurface, real_t rCMB );

std::shared_ptr< hyteg::PrimitiveStorage > createAnnulusStorage( uint_t nTan,
                                                                 uint_t nRad,
                                                                 real_t rCMB                = 1.0,
                                                                 real_t rSurface            = 2.0,
                                                                 bool   blending            = true,
                                                                 real_t boundaryTolerance   = 1e-10,
                                                                 uint_t flagSurface         = 1,
                                                                 uint_t flagCMB             = 2,
                                                                 uint_t additionalHaloDepth = 1 );

ConvectionBC createAnnulusBoundaryConditions( DoFType typeSurface = DoFType::DirichletBoundary,
                                              DoFType typeCMB     = DoFType::FreeslipBoundary,
                                              uint_t  flagSurface = 1,
                                              uint_t  flagCMB     = 2 );

///////////////////////////////////////////////////////
////////////////////// Rectangle //////////////////////
///////////////////////////////////////////////////////
bool RectangleCMBBoundary( const hyteg::Point3D& x, real_t yMin, real_t boundaryTolerance );
bool RectangleSurfaceBoundary( const hyteg::Point3D& x, real_t yMax, real_t boundaryTolerance );
bool RectangleLeftBoundary( const hyteg::Point3D& x, real_t xMin, real_t boundaryTolerance );
bool RectangleRightBoundary( const hyteg::Point3D& x, real_t xMax, real_t boundaryTolerance );
bool RectangleCornerLeftCMBBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMin, real_t boundaryTolerance );
bool RectangleCornerRightCMBBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMin, real_t boundaryTolerance );
bool RectangleCornerRightSurfaceBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMax, real_t boundaryTolerance );
bool RectangleCornerLeftSurfaceBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMax, real_t boundaryTolerance );
bool RectangleCornerBoundary( const hyteg::Point3D& x,
                              real_t                xMin,
                              real_t                xMax,
                              real_t                yMin,
                              real_t                yMax,
                              real_t                boundaryTolerance );
void RectangleBoundaryNormal( const Point3D& x,
                              Point3D&       normal,
                              real_t         xMin,
                              real_t         xMax,
                              real_t         yMin,
                              real_t         yMax,
                              real_t         boundaryTolerance );

std::function< bool( const hyteg::Point3D& x ) > createRectangleCMBBoundaryFct( real_t yMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createRectangleSurfaceBoundaryFct( real_t yMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createRectangleLeftBoundaryFct( real_t xMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createRectangleRightBoundaryFct( real_t xMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerLeftCMBBoundaryFct( real_t xMin, real_t yMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerRightCMBBoundaryFct( real_t xMax, real_t yMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerRightSurfaceBoundaryFct( real_t xMax, real_t yMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerLeftSurfaceBoundaryFct( real_t xMin, real_t yMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createRectangleCornerBoundaryFct( real_t xMin, real_t xMax, real_t yMin, real_t yMax, real_t boundaryTolerance );
std::function< void( const Point3D& in, Point3D& out ) >
    createRectangleBoundaryNormalFct( real_t xMin, real_t xMax, real_t yMin, real_t yMax, real_t boundaryTolerance );

std::shared_ptr< hyteg::PrimitiveStorage > createRectangleStorage( uint_t  nX,
                                                                   uint_t  nY,
                                                                   Point2D P1                     = Point2D( { 0, 0 } ),
                                                                   Point2D P2                     = Point2D( { 1, 1 } ),
                                                                   real_t  boundaryTolerance      = 1e-10,
                                                                   uint_t  flagSurface            = 1,
                                                                   uint_t  flagCMB                = 2,
                                                                   uint_t  flagLeft               = 3,
                                                                   uint_t  flagRight              = 4,
                                                                   uint_t  flagCornerLeftCMB      = 5,
                                                                   uint_t  flagCornerRightCMB     = 6,
                                                                   uint_t  flagCornerRightSurface = 7,
                                                                   uint_t  flagCornerLeftSurface  = 8,
                                                                   uint_t  additionalHaloDepth    = 1 );

ConvectionBC createRectangleBoundaryConditions( DoFType typeSurface            = DoFType::DirichletBoundary,
                                                DoFType typeCMB                = DoFType::FreeslipBoundary,
                                                DoFType typeLeft               = DoFType::FreeslipBoundary,
                                                DoFType typeRight              = DoFType::FreeslipBoundary,
                                                DoFType typeTemperatureSide    = DoFType::NeumannBoundary,
                                                uint_t  flagSurface            = 1,
                                                uint_t  flagCMB                = 2,
                                                uint_t  flagLeft               = 3,
                                                uint_t  flagRight              = 4,
                                                uint_t  flagCornerLeftCMB      = 5,
                                                uint_t  flagCornerRightCMB     = 6,
                                                uint_t  flagCornerRightSurface = 7,
                                                uint_t  flagCornerLeftSurface  = 8 );

//////////////////////////////////////////////////////
//////////////////////// Cuboid ////////////////////////
//////////////////////////////////////////////////////
bool CuboidCMBBoundary( const hyteg::Point3D& x, real_t zMin, real_t boundaryTolerance );
bool CuboidSurfaceBoundary( const hyteg::Point3D& x, real_t zMax, real_t boundaryTolerance );
bool CuboidLeftBoundary( const hyteg::Point3D& x, real_t xMin, real_t boundaryTolerance );
bool CuboidRightBoundary( const hyteg::Point3D& x, real_t xMax, real_t boundaryTolerance );
bool CuboidFrontBoundary( const hyteg::Point3D& x, real_t yMin, real_t boundaryTolerance );
bool CuboidBackBoundary( const hyteg::Point3D& x, real_t yMax, real_t boundaryTolerance );

bool CuboidLeftFrontEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMin, real_t boundaryTolerance );
bool CuboidRightFrontEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMin, real_t boundaryTolerance );
bool CuboidRightBackEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMax, real_t boundaryTolerance );
bool CuboidLeftBackEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMax, real_t boundaryTolerance );
bool CuboidLeftCMBEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t zMin, real_t boundaryTolerance );
bool CuboidRightCMBEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t zMin, real_t boundaryTolerance );
bool CuboidRightSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t xMax, real_t zMax, real_t boundaryTolerance );
bool CuboidLeftSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t xMin, real_t zMax, real_t boundaryTolerance );
bool CuboidFrontCMBEdgeBoundary( const hyteg::Point3D& x, real_t yMin, real_t zMin, real_t boundaryTolerance );
bool CuboidBackCMBEdgeBoundary( const hyteg::Point3D& x, real_t yMax, real_t zMin, real_t boundaryTolerance );
bool CuboidBackSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t yMax, real_t zMax, real_t boundaryTolerance );
bool CuboidFrontSurfaceEdgeBoundary( const hyteg::Point3D& x, real_t yMin, real_t zMax, real_t boundaryTolerance );
bool CuboidEdgeBoundary( const hyteg::Point3D& x,
                         real_t                xMin,
                         real_t                xMax,
                         real_t                yMin,
                         real_t                yMax,
                         real_t                zMin,
                         real_t                zMax,
                         real_t                boundaryTolerance );

bool CuboidSurfaceLeftFrontCornerBoundary( const hyteg::Point3D& x,
                                           real_t                xMin,
                                           real_t                yMin,
                                           real_t                zMax,
                                           real_t                boundaryTolerance );
bool CuboidSurfaceLeftBackCornerBoundary( const hyteg::Point3D& x,
                                          real_t                xMin,
                                          real_t                yMax,
                                          real_t                zMax,
                                          real_t                boundaryTolerance );
bool CuboidSurfaceRightFrontCornerBoundary( const hyteg::Point3D& x,
                                            real_t                xMax,
                                            real_t                yMin,
                                            real_t                zMax,
                                            real_t                boundaryTolerance );
bool CuboidSurfaceRightBackCornerBoundary( const hyteg::Point3D& x,
                                           real_t                xMax,
                                           real_t                yMax,
                                           real_t                zMax,
                                           real_t                boundaryTolerance );
bool CuboidCMBLeftFrontCornerBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMin, real_t zMin, real_t boundaryTolerance );
bool CuboidCMBLeftBackCornerBoundary( const hyteg::Point3D& x, real_t xMin, real_t yMax, real_t zMin, real_t boundaryTolerance );
bool CuboidCMBRightFrontCornerBoundary( const hyteg::Point3D& x,
                                        real_t                xMax,
                                        real_t                yMin,
                                        real_t                zMin,
                                        real_t                boundaryTolerance );
bool CuboidCMBRightBackCornerBoundary( const hyteg::Point3D& x, real_t xMax, real_t yMax, real_t zMin, real_t boundaryTolerance );
bool CuboidCornerBoundary( const hyteg::Point3D& x,
                           real_t                xMin,
                           real_t                xMax,
                           real_t                yMin,
                           real_t                yMax,
                           real_t                zMin,
                           real_t                zMax,
                           real_t                boundaryTolerance );
void CuboidBoundaryNormal( const Point3D& x,
                           Point3D&       normal,
                           real_t         xMin,
                           real_t         xMax,
                           real_t         yMin,
                           real_t         yMax,
                           real_t         zMin,
                           real_t         zMax,
                           real_t         boundaryTolerance );

std::function< bool( const hyteg::Point3D& x ) > createCuboidCMBBoundaryFct( real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createCuboidSurfaceBoundaryFct( real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createCuboidLeftBoundaryFct( real_t xMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createCuboidRightBoundaryFct( real_t xMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createCuboidFrontBoundaryFct( real_t yMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createCuboidBackBoundaryFct( real_t yMax, real_t boundaryTolerance );

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftFrontEdgeBoundaryFct( real_t xMin, real_t yMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightFrontEdgeBoundaryFct( real_t xMax, real_t yMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightBackEdgeBoundaryFct( real_t xMax, real_t yMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftBackEdgeBoundaryFct( real_t xMin, real_t yMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftCMBEdgeBoundaryFct( real_t xMin, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightCMBEdgeBoundaryFct( real_t xMax, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidRightSurfaceEdgeBoundaryFct( real_t xMax, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidLeftSurfaceEdgeBoundaryFct( real_t xMin, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidFrontCMBEdgeBoundaryFct( real_t yMin, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidBackCMBEdgeBoundaryFct( real_t yMax, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidBackSurfaceEdgeBoundaryFct( real_t yMax, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidFrontSurfaceEdgeBoundaryFct( real_t yMin, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) > createCuboidEdgeBoundaryFct( real_t xMin,
                                                                              real_t xMax,
                                                                              real_t yMin,
                                                                              real_t yMax,
                                                                              real_t zMin,
                                                                              real_t zMax,
                                                                              real_t boundaryTolerance );

std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceLeftFrontCornerBoundaryFct( real_t xMin, real_t yMin, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceLeftBackCornerBoundaryFct( real_t xMin, real_t yMax, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceRightFrontCornerBoundaryFct( real_t xMax, real_t yMin, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidSurfaceRightBackCornerBoundaryFct( real_t xMax, real_t yMax, real_t zMax, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBLeftFrontCornerBoundaryFct( real_t xMin, real_t yMin, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBLeftBackCornerBoundaryFct( real_t xMin, real_t yMax, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBRightFrontCornerBoundaryFct( real_t xMax, real_t yMin, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >
    createCuboidCMBRightBackCornerBoundaryFct( real_t xMax, real_t yMax, real_t zMin, real_t boundaryTolerance );
std::function< bool( const hyteg::Point3D& x ) >         createCuboidCornerBoundaryFct( real_t xMin,
                                                                                        real_t xMax,
                                                                                        real_t yMin,
                                                                                        real_t yMax,
                                                                                        real_t zMin,
                                                                                        real_t zMax,
                                                                                        real_t boundaryTolerance );
std::function< void( const Point3D& in, Point3D& out ) > createCuboidBoundaryNormalFct( real_t xMin,
                                                                                        real_t xMax,
                                                                                        real_t yMin,
                                                                                        real_t yMax,
                                                                                        real_t zMin,
                                                                                        real_t zMax,
                                                                                        real_t boundaryTolerance );

ConvectionBC createCuboidBoundaryConditions( DoFType typeSurface         = DoFType::DirichletBoundary,
                                             DoFType typeCMB             = DoFType::FreeslipBoundary,
                                             DoFType typeLeft            = DoFType::FreeslipBoundary,
                                             DoFType typeRight           = DoFType::FreeslipBoundary,
                                             DoFType typeFront           = DoFType::FreeslipBoundary,
                                             DoFType typeBack            = DoFType::FreeslipBoundary,
                                             DoFType typeTemperatureSide = DoFType::NeumannBoundary,
                                             uint_t  flagSurface         = 1,
                                             uint_t  flagCMB             = 2,
                                             uint_t  flagLeft            = 3,
                                             uint_t  flagRight           = 4,
                                             uint_t  flagFront           = 5,
                                             uint_t  flagBack            = 6,

                                             uint_t flagLeftFrontEdge    = 7,
                                             uint_t flagRightFrontEdge   = 8,
                                             uint_t flagRightBackEdge    = 9,
                                             uint_t flagLeftBackEdge     = 10,
                                             uint_t flagLeftCMBEdge      = 11,
                                             uint_t flagRightCMBEdge     = 12,
                                             uint_t flagRightSurfaceEdge = 13,
                                             uint_t flagLeftSurfaceEdge  = 14,
                                             uint_t flagFrontCMBEdge     = 15,
                                             uint_t flagBackCMBEdge      = 16,
                                             uint_t flagBackSurfaceEdge  = 17,
                                             uint_t flagFrontSurfaceEdge = 18,

                                             uint_t flagSurfaceLeftFrontCorner  = 19,
                                             uint_t flagSurfaceLeftBackCorner   = 20,
                                             uint_t flagSurfaceRightFrontCorner = 21,
                                             uint_t flagSurfaceRightBackCorner  = 22,
                                             uint_t flagCMBLeftFrontCorner      = 23,
                                             uint_t flagCMBLeftBackCorner       = 24,
                                             uint_t flagCMBRightFrontCorner     = 25,
                                             uint_t flagCMBRightBackCorner      = 26 );

std::shared_ptr< hyteg::PrimitiveStorage > createCuboidStorage( uint_t  nX,
                                                                uint_t  nY,
                                                                uint_t  nZ,
                                                                Point3D P1                = Point3D( { 0, 0, 0 } ),
                                                                Point3D P2                = Point3D( { 1, 1, 1 } ),
                                                                real_t  boundaryTolerance = 1e-10,
                                                                uint_t  flagSurface       = 1,
                                                                uint_t  flagCMB           = 2,
                                                                uint_t  flagLeft          = 3,
                                                                uint_t  flagRight         = 4,
                                                                uint_t  flagFront         = 5,
                                                                uint_t  flagBack          = 6,

                                                                uint_t flagLeftFrontEdge    = 7,
                                                                uint_t flagRightFrontEdge   = 8,
                                                                uint_t flagRightBackEdge    = 9,
                                                                uint_t flagLeftBackEdge     = 10,
                                                                uint_t flagLeftCMBEdge      = 11,
                                                                uint_t flagRightCMBEdge     = 12,
                                                                uint_t flagRightSurfaceEdge = 13,
                                                                uint_t flagLeftSurfaceEdge  = 14,
                                                                uint_t flagFrontCMBEdge     = 15,
                                                                uint_t flagBackCMBEdge      = 16,
                                                                uint_t flagBackSurfaceEdge  = 17,
                                                                uint_t flagFrontSurfaceEdge = 18,

                                                                uint_t flagSurfaceLeftFrontCorner  = 19,
                                                                uint_t flagSurfaceLeftBackCorner   = 20,
                                                                uint_t flagSurfaceRightFrontCorner = 21,
                                                                uint_t flagSurfaceRightBackCorner  = 22,
                                                                uint_t flagCMBLeftFrontCorner      = 23,
                                                                uint_t flagCMBLeftBackCorner       = 24,
                                                                uint_t flagCMBRightFrontCorner     = 25,
                                                                uint_t flagCMBRightBackCorner      = 26,
                                                                uint_t additionalHaloDepth         = 1 );
} // namespace convectionToolbox
} // namespace hyteg
