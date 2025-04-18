/*
 * Copyright (c) 2021-2024 Benjamin Mann
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

#include "mesh.hpp"

#include <algorithm>
#include <assert.h>
#include <core/mpi/Broadcast.h>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <type_traits>

#include "hyteg/Algorithms.hpp"

#include "refine_cell.hpp"
#include "simplexFactory.hpp"

namespace hyteg {
namespace adaptiveRefinement {

template < class K_Simplex >
K_Mesh< K_Simplex >::K_Mesh( const SetupPrimitiveStorage& setupStorage )
: _n_processes( setupStorage.getNumberOfProcesses() )
{
   // copy geometry maps for each primitive in initial setupStorage
   SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
   setupStorage.getSetupPrimitives( setupPrimitives );
   for ( auto& [id, primitive] : setupPrimitives )
   {
      _geometryMap[id] = primitive->getGeometryMap();
   }
   _geometryMap[_invalidID] = nullptr; // used for uninitialized values

   // internal data structures are only required on rank_0
   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      // extract vertices
      _n_vertices = setupStorage.getVertices().size();

      // [0,1,...,n-1]
      std::vector< uint_t > vtxIndices( _n_vertices );
      // convert PrimitiveID::ID to Mesh::vertexID
      std::map< PrimitiveID, uint_t > conversion;

      // initialize vertices
      for ( auto& [id, vtx] : setupStorage.getVertices() )
      {
         auto idx = _V.get_next_idx();
         _coords.append( vtx->getCoordinates() );
         _V.append( VertexData( id, vtx->getMeshBoundaryFlag(), id, { { idx } }, setupStorage.getTargetRank( id ) ) );

         // prepare element setup
         conversion[id]  = idx;
         vtxIndices[idx] = idx;
      }

      // initialize all edges, faces and cells

      SimplexFactory fac( nullptr, vtxIndices );
      // simplex factory does not store cells
      std::set< std::shared_ptr< Simplex3 > > cells;
      std::vector< PrimitiveID >              v;

      for ( auto& [id, edge] : setupStorage.getEdges() )
      {
         edge->getNeighborVertices( v );
         auto myEdge = fac.make_edge( conversion[v[0]], conversion[v[1]] );
         myEdge->setPrimitiveID( id );
         myEdge->setGeometryMap( id );
         myEdge->setBoundaryFlag( edge->getMeshBoundaryFlag() );
         myEdge->setTargetRank( setupStorage.getTargetRank( PrimitiveID( id ) ) );
      }

      for ( auto& [id, face] : setupStorage.getFaces() )
      {
         face->getNeighborVertices( v );
         auto myFace = fac.make_face( conversion[v[0]], conversion[v[1]], conversion[v[2]] );
         myFace->setPrimitiveID( id );
         myFace->setGeometryMap( id );
         myFace->setBoundaryFlag( face->getMeshBoundaryFlag() );
         myFace->setTargetRank( setupStorage.getTargetRank( PrimitiveID( id ) ) );
      }

      for ( auto& [id, cell] : setupStorage.getCells() )
      {
         cell->getNeighborVertices( v );
         auto myCell = fac.make_cell( conversion[v[0]], conversion[v[1]], conversion[v[2]], conversion[v[3]] );
         myCell->setPrimitiveID( id );
         myCell->setGeometryMap( id );
         myCell->setBoundaryFlag( cell->getMeshBoundaryFlag() );
         myCell->setTargetRank( setupStorage.getTargetRank( PrimitiveID( id ) ) );
         cells.insert( myCell );
      }

      // insert volume elements into _T
      if constexpr ( VOL == FACE )
      {
         if ( !cells.empty() )
         {
            WALBERLA_ABORT( "Adaptive 2D mesh requires SetupPrimitiveStorage without any cells!" );
         }
         for ( auto& p : fac.faces() )
         {
            _T.insert( p.second );
         }
      }
      if constexpr ( VOL == CELL )
      {
         if ( cells.empty() )
         {
            WALBERLA_ABORT( "Adaptive 3D mesh requires SetupPrimitiveStorage containing at least one cell!" );
         }

         _T = cells;
      }

      _n_elements = _T.size();
   }

   check_integrity( 0 );

   // broadcast required values to all processes
   walberla::mpi::broadcastObject( _n_elements );
   walberla::mpi::broadcastObject( _n_vertices );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::refine_uniform( uint_t k )
{
   for ( uint_t j = 0; j < k; ++j )
   {
      std::vector< PrimitiveID > all( _T.size() );

      uint_t i = 0;
      for ( auto& el : _T )
      {
         all[i] = el->getPrimitiveID();
         ++i;
      }

      refineRG( all, {} );
   }
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::refineRG( const std::vector< PrimitiveID >& elements_to_refine,
                                    const std::vector< PrimitiveID >& elements_to_coarsen )
{
   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      if ( elements_to_refine.size() + elements_to_coarsen.size() > n_elements() )
      {
         WALBERLA_ABORT( "Can't apply both refinement and coarsening to the same element!" );
      }
      // get elements corresponding to given IDs
      auto [R, Pc] = init_R_Pc( elements_to_refine, elements_to_coarsen );

      // undo last green refinement step to prevent mesh degeneration
      remove_green_edges();

      // for each t in Pc, add t to T and remove the children of t from T
      unrefine( Pc );

      /* iteratively apply red refinement to elements that would otherwise
         be subject to multiple green refinement steps later on
      */
      while ( !R.empty() )
      {
         refine_red( R );
         R = find_elements_for_red_refinement();
      }

      // apply green refinement
      refine_green();

      // update current configuration
      _n_vertices = _coords.size();
      _n_elements = _T.size();
   }

   walberla::mpi::broadcastObject( _n_vertices );
   walberla::mpi::broadcastObject( _n_elements );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::refineRG( const ErrorVector&                                         errors_local,
                                    const std::function< bool( const ErrorVector&, uint_t ) >& criterion_r,
                                    const std::function< bool( const ErrorVector&, uint_t ) >& criterion_c )
{
   ErrorVector err_glob;
   gatherGlobalError( errors_local, err_glob );

   // apply criterion
   std::vector< PrimitiveID > elements_to_refine;
   std::vector< PrimitiveID > elements_to_coarsen;
   for ( uint_t i = 0; i < _n_elements; ++i )
   {
      bool do_r = criterion_r( err_glob, i );
      bool do_c = criterion_c( err_glob, i );
      auto id   = err_glob[i].second;
      if ( do_r && do_c )
      {
         WALBERLA_ABORT( "Can't apply both refinement and coarsening to element " << id << "!" );
      }
      if ( do_r )
      {
         elements_to_refine.push_back( id );
      }
      if ( do_c )
      {
         elements_to_coarsen.push_back( id );
      }
   }

   refineRG( elements_to_refine, elements_to_coarsen );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::refineRG( const ErrorVector& errors_local, Strategy ref_strat, Strategy cors_strat, bool verbose )
{
   ErrorVector err_glob;
   gatherGlobalError( errors_local, err_glob );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Adaptive refinement:" );
      WALBERLA_LOG_INFO_ON_ROOT( " -> min_i e_i = " << err_glob.back().first );
      WALBERLA_LOG_INFO_ON_ROOT( " -> max_i e_i = " << err_glob.front().first );
   }

   auto compute_weighted_mean = [&]( real_t p ) {
      if ( p <= -std::numeric_limits< real_t >::infinity() )
      {
         // p = -∞ => e_mean = e_{n-1} (smallest error), ie. refine all / coarsen none
         return err_glob.back().first;
      }
      if ( p >= std::numeric_limits< real_t >::infinity() )
      {
         // p = ∞ => e_mean = e_{0} (greatest error), ie. refine T_0 / coarsen all except T_0
         return err_glob.front().first;
      }
      // else

      // disable weighting if p is close to zero
      bool weighted = ( std::abs( p ) > real_t( 0.01 ) );

      auto e_mean_w = real_t( 0.0 );
      auto sum_w    = real_t( 0.0 );
      auto n        = err_glob.size();
      for ( uint_t n_i = 1; n_i <= n; ++n_i )
      {
         // iterate in reverse order to add smallest contribution first, reducing round off error
         auto   i   = n - n_i; // n - (n-i) = i
         real_t e_i = err_glob[i].first;
         if ( weighted )
         {
            real_t w_i = std::pow( real_t( n_i ), p );
            e_mean_w += w_i * e_i;
            sum_w += w_i;
         }
         else
         {
            e_mean_w += e_i;
            sum_w = real_t( n );
         }
      }
      return e_mean_w /= sum_w;
   };

   std::function< bool( real_t, uint_t ) > crit_r;
   std::function< bool( real_t, uint_t ) > crit_c;

   if ( ref_strat.t == Strategy::WEIGHTED_MEAN )
   {
      auto e_mean = compute_weighted_mean( ref_strat.p );
      crit_r      = [e_mean]( real_t e_i, uint_t ) -> bool { return e_i >= e_mean; };
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " -> refinemnt strategy: weighted mean with weights w_i = (n-i)^(" << ref_strat.p << ")" );
         WALBERLA_LOG_INFO_ON_ROOT( " -> refining all elements i where e_i >= " << e_mean );
      }
   }
   else if ( ref_strat.t == Strategy::PERCENTILE )
   {
      auto sizeR = uint_t( std::round( real_t( _n_elements ) * ref_strat.p ) );
      crit_r     = [=]( real_t, uint_t i ) -> bool { return i < sizeR; };
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " -> refinemnt strategy: percentile with p = " << ref_strat.p * 100 << "%" );
         WALBERLA_LOG_INFO_ON_ROOT( " -> refining all elements i where e_i >= " << err_glob[sizeR - 1].first );
      }
   }
   else if ( ref_strat.t == Strategy::DOERFLER )
   {
      // compute sum of all errors
      auto e_sum = real_t( 0.0 );
      for ( uint_t n_i = 1; n_i <= _n_elements; ++n_i )
      {
         // iterate in reverse order to add smallest contribution first, reducing round off error
         auto i = _n_elements - n_i; // i = n - (n-i)
         e_sum += err_glob[i].first;
      }
      // find smallest subset R of T where ∑_t∈R e_t >= p * ∑_t∈T e_t
      auto   e_sum_R = real_t( 0.0 );
      uint_t sizeR   = uint_t( 0 );
      while ( sizeR < _n_elements && e_sum_R < ref_strat.p * e_sum )
      {
         e_sum_R += err_glob[sizeR].first;
         ++sizeR;
      }
      // refine all elements in R
      crit_r = [=]( real_t, uint_t i ) -> bool { return i < sizeR; };
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " -> refinemnt strategy: Dörfler marking with p = " << ref_strat.p );
         WALBERLA_LOG_INFO_ON_ROOT( " -> refining all elements i where e_i >= " << err_glob[sizeR - 1].first );
      }
   }
   else
   {
      WALBERLA_ABORT( "Unknown refinement strategy!" );
   }

   std::vector< PrimitiveID > elements_to_refine;
   for ( uint_t i = 0; i < _n_elements; ++i )
   {
      if ( crit_r( err_glob[i].first, i ) )
      {
         elements_to_refine.push_back( err_glob[i].second );
      }
      else
      {
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( " -> " << i << " elements marked for refinement." );
         }
         break;
      }
   }

   if ( cors_strat.t == Strategy::WEIGHTED_MEAN )
   {
      auto e_mean = compute_weighted_mean( cors_strat.p );
      crit_c      = [e_mean]( real_t e_i, uint_t ) -> bool { return e_i < e_mean; };
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " -> coarsening strategy: weighted mean with weights w_i = (n-i)^(" << cors_strat.p << ")" );
         WALBERLA_LOG_INFO_ON_ROOT( " -> unrefine all elements i where e_i < " << e_mean );
      }
   }
   else if ( cors_strat.t == Strategy::PERCENTILE )
   {
      auto sizeC = uint_t( std::round( real_t( _n_elements ) * cors_strat.p ) );
      crit_c     = [=]( real_t, uint_t i ) -> bool { return i >= _n_elements - sizeC; };
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " -> coarsening strategy: percentile with p = " << cors_strat.p * 100 << "%" );
         WALBERLA_LOG_INFO_ON_ROOT( " -> unrefine all elements i where e_i <= " << err_glob[_n_elements - sizeC].first );
      }
   }
   else if ( cors_strat.t == Strategy::MULTIPLE_OF_R )
   {
      auto sizeR = elements_to_refine.size();
      auto sizeC = uint_t( std::round( real_t( sizeR ) * cors_strat.p ) );
      if ( sizeR + sizeC > _n_elements )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "n_C + n_R > n_el for n_C = p*n_R -> using n_C = n_el - n_R instead!" );
         sizeC = _n_elements - sizeR;
      }
      crit_c = [=]( real_t, uint_t i ) -> bool { return i >= _n_elements - sizeC; };
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " -> coarsening strategy: coarsen p*|R| elements where p = " << cors_strat.p );
         WALBERLA_LOG_INFO_ON_ROOT( " -> unrefine all elements i where e_i <= " << err_glob[_n_elements - sizeC].first );
      }
   }
   else
   {
      WALBERLA_ABORT( "Unknown coarsening strategy!" );
   }

   std::vector< PrimitiveID > elements_to_coarsen;
   for ( uint_t n_i = 1; n_i <= _n_elements; ++n_i )
   {
      auto i = _n_elements - n_i; // n - (n-i) = i
      if ( crit_c( err_glob[i].first, i ) )
      {
         elements_to_coarsen.push_back( err_glob[i].second );
      }
      else
      {
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( " -> " << ( n_i - 1 ) << " elements marked for coarsening." );
         }
         break;
      }
   }

   refineRG( elements_to_refine, elements_to_coarsen );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::gatherGlobalError( const ErrorVector& err_loc, ErrorVector& err_glob_sorted ) const
{
   ErrorVector err_other;

   walberla::mpi::SendBuffer send;
   walberla::mpi::RecvBuffer recv;

   send << err_loc;
   walberla::mpi::allGathervBuffer( send, recv );
   for ( uint_t rnk = 0; rnk < _n_processes; ++rnk )
   {
      recv >> err_other;
      err_glob_sorted.insert( err_glob_sorted.end(), err_other.begin(), err_other.end() );
   }

   if ( err_glob_sorted.size() != _n_elements )
   {
      WALBERLA_ABORT( "total number of error values must be equal to number of macro elements (cells/faces)" );
   }

   // sort by errors
   std::sort( err_glob_sorted.begin(), err_glob_sorted.end(), std::greater< std::pair< real_t, PrimitiveID > >() );
}

template < class K_Simplex >
std::shared_ptr< PrimitiveStorage > K_Mesh< K_Simplex >::make_storage()
{
   std::map< PrimitiveID, EdgeData > edges;
   std::map< PrimitiveID, FaceData > faces;
   std::map< PrimitiveID, CellData > cells;

   extract_data( edges, faces, cells );

   walberla::mpi::broadcastObject( edges );
   walberla::mpi::broadcastObject( faces );
   walberla::mpi::broadcastObject( cells );

   return make_localPrimitives( edges, faces, cells );
}

template < class K_Simplex >
MigrationInfo K_Mesh< K_Simplex >::loadbalancing( const Loadbalancing& lb,
                                                  const bool           migrationInfo_required,
                                                  const bool           allow_split_siblings,
                                                  const bool           verbose )
{
   if ( !allow_split_siblings )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "loadbalancing with allow_split_siblings disabled may result in bad load balance!" );
   }

   // neighborhood information
   NeighborhoodMap nbrHood;
   // old and new target rank
   std::map< PrimitiveID, std::pair< uint_t, uint_t > > targetRank;

   // collect neighbor interfaces of volumes and vice versa
   auto addNeighbor = [&]( const PrimitiveID volId, const PrimitiveType& pt, const PrimitiveID& nbrId ) {
      // add nbrId as neighbor to volId and vice versa
      nbrHood[VOL][volId][pt].push_back( nbrId );
      nbrHood[pt][nbrId][VOL].push_back( volId );
      // also add direct volume neighbors, i.e. elements that share a face (edge in 2d)
      if ( pt == ( VOL - 1 ) && nbrHood[pt][nbrId][VOL].size() == 2 )
      {
         auto nbrVolId = nbrHood[pt][nbrId][VOL][0];
         // add nbrVolId as neighbor to volId and vice versa
         nbrHood[VOL][volId][VOL].push_back( nbrVolId );
         nbrHood[VOL][nbrVolId][VOL].push_back( volId );
      }
   };

   // reset target ranks and gather required info
   for ( auto& el : _T )
   {
      auto elId = el->getPrimitiveID();

      targetRank[elId].first = el->getTargetRank();
      el->setTargetRank( _n_processes );

      if constexpr ( VOL == CELL )
      {
         for ( auto& face : el->get_faces() )
         {
            if ( face->getTargetRank() < _n_processes )
            {
               targetRank[face->getPrimitiveID()].first = face->getTargetRank();
               face->setTargetRank( _n_processes );
            }
            addNeighbor( elId, FACE, face->getPrimitiveID() );
         }
      }
      for ( auto& edge : el->get_edges() )
      {
         if ( edge->getTargetRank() < _n_processes )
         {
            targetRank[edge->getPrimitiveID()].first = edge->getTargetRank();
            edge->setTargetRank( _n_processes );
         }
         addNeighbor( elId, EDGE, edge->getPrimitiveID() );
      }
      for ( auto& idx : el->get_vertices() )
      {
         if ( _V[idx].getTargetRank() < _n_processes )
         {
            targetRank[_V[idx].getPrimitiveID()].first = _V[idx].getTargetRank();
            _V[idx].setTargetRank( _n_processes );
         }
         addNeighbor( elId, VTX, _V[idx].getPrimitiveID() );
      }
   }

   std::vector< uint_t > n_vol_on_rnk;

   // assign volume primitives
   if ( lb == ROUND_ROBIN )
   {
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Load balancing strategy: Round Robin" );
      }
      n_vol_on_rnk = loadbalancing_roundRobin( allow_split_siblings );
   }
   else if ( lb == GREEDY )
   {
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Load balancing strategy: Greedy" );
      }
      n_vol_on_rnk = loadbalancing_greedy( nbrHood, allow_split_siblings );
   }
   else
   {
      // todo: implement better loadbalancing
      WALBERLA_LOG_WARNING_ON_ROOT( "loadbalancing scheme not implemented! Using Round Robin instead" );
      n_vol_on_rnk = loadbalancing_roundRobin( allow_split_siblings );
   }

   if ( verbose )
   {
      uint_t max_load = 0;
      uint_t min_load = _n_elements;
      for ( auto& load : n_vol_on_rnk )
      {
         max_load = std::max( load, max_load );
         min_load = std::min( load, min_load );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "min/max/mean load (number of vol. elements): "
                                 << min_load << " / " << max_load << " / " << double( _n_elements ) / double( _n_processes ) );
   }

   // assign interface primitives to processes
   auto computeTargetRank = [&]( uint_t curRank, const std::vector< hyteg::PrimitiveID >& myNbrVolumes ) -> uint_t {
      if ( curRank < _n_processes )
      {
         return curRank; // already assigned
      }

      // loop over all neighboring volume primitives and extract their target ranks
      std::map< uint_t, uint_t > nbrRnks;
      for ( auto& nbrID : myNbrVolumes )
      {
         nbrRnks[targetRank[nbrID].second] += 1;
      }

      // find rank where the most neighboring volume primitives are located
      auto cmp = []( const std::pair< uint_t, uint_t >& a, const std::pair< uint_t, uint_t >& b ) { return a.second < b.second; };
      return std::max_element( nbrRnks.begin(), nbrRnks.end(), cmp )->first;
   };

   for ( auto& el : _T )
   {
      targetRank[el->getPrimitiveID()].second = el->getTargetRank();
   }

   for ( auto& el : _T )
   {
      // neighbor faces (3d only)
      if constexpr ( VOL == CELL )
      {
         for ( auto& face : el->get_faces() )
         {
            auto id      = face->getPrimitiveID();
            auto newRank = computeTargetRank( face->getTargetRank(), nbrHood[FACE].at( id )[VOL] );
            face->setTargetRank( newRank );
            targetRank[id].second = newRank;
         }
      }
      // neighbor edges
      for ( auto& edge : el->get_edges() )
      {
         auto id      = edge->getPrimitiveID();
         auto newRank = computeTargetRank( edge->getTargetRank(), nbrHood[EDGE].at( id )[VOL] );
         edge->setTargetRank( newRank );
         targetRank[id].second = newRank;
      }
      // neighbor vertices
      for ( auto& idx : el->get_vertices() )
      {
         auto id      = _V[idx].getPrimitiveID();
         auto newRank = computeTargetRank( _V[idx].getTargetRank(), nbrHood[VTX].at( id )[VOL] );
         _V[idx].setTargetRank( newRank );
         targetRank[id].second = newRank;
      }
   }

   // gather migration data
   MigrationMap_T migrationMap;
   uint_t         numReceivingPrimitives = 0;
   if ( migrationInfo_required )
   {
      const uint_t rank = uint_t( walberla::mpi::MPIManager::instance()->rank() );
      // broadcast data to all processes
      walberla::mpi::broadcastObject( targetRank );

      for ( auto& [id, rnk] : targetRank )
      {
         if ( rank == rnk.first )
         {
            migrationMap[id] = rnk.second;
         }
         if ( rank == rnk.second )
         {
            ++numReceivingPrimitives;
         }
      }
   }
   return MigrationInfo( migrationMap, numReceivingPrimitives );
}

template < class K_Simplex >
std::vector< uint_t > K_Mesh< K_Simplex >::loadbalancing_roundRobin( const bool allow_split_siblings )
{
   if ( walberla::mpi::MPIManager::instance()->rank() != 0 )
      return std::vector< uint_t >( 0 );

   // apply roundRobin to volume elements
   uint_t                targetRnk = 0;
   std::vector< uint_t > n_vol_on_rnk( _n_processes );
   uint_t                n_vol_max = _n_elements / _n_processes;
   for ( auto el : _T )
   {
      // already assigned
      if ( el->getTargetRank() < _n_processes )
      {
         continue;
      }

      // process saturated -> move to next process
      while ( n_vol_on_rnk[targetRnk] >= n_vol_max )
      {
         ++targetRnk;

         // round robin complete but unassigned elements left
         if ( targetRnk == _n_processes )
         {
            targetRnk = 0;
            ++n_vol_max;
         }
      }

      // assign el to targetRnk
      if ( !allow_split_siblings && el->has_green_edge() )
      {
         /* elements coming from a green step must reside on the same process
            to ensure compatibility with further refinement (interpolation between grids)
         */
         for ( auto s : el->get_siblings() )
         {
            s->setTargetRank( targetRnk );
            ++n_vol_on_rnk[targetRnk];
         }
         // also change rank of parent s.th. next refinement resides on same process
         el->get_parent()->setTargetRank( targetRnk );
      }
      else
      {
         el->setTargetRank( targetRnk );
         ++n_vol_on_rnk[targetRnk];
      }
   }
   return n_vol_on_rnk;
}

template < class K_Simplex >
std::vector< uint_t > K_Mesh< K_Simplex >::loadbalancing_greedy( const NeighborhoodMap& nbrHood, const bool allow_split_siblings )
{
   if ( walberla::mpi::MPIManager::instance()->rank() != 0 )
      return std::vector< uint_t >( 0 );

   std::vector< uint_t > n_vol_on_rnk( _n_processes );
   uint_t                n_vol_max0     = _n_elements / _n_processes;
   uint_t                n_vol_max1     = n_vol_max0 + std::min( uint_t( 1 ), _n_elements % _n_processes );
   uint_t                n_assigned     = 0;
   auto                  random_element = _T.begin();

   std::map< PrimitiveID, K_Simplex* > id_to_el;
   for ( auto& el : _T )
   {
      id_to_el[el->getPrimitiveID()] = el.get();
   }

   for ( bool fill_any = false; n_assigned < _n_elements; fill_any = true )
   {
      // in the first run, all ranks are assigned at least n_vol_max0 elements
      // if there are still unassigned elements left, we distribute them at random (fill_any)

      for ( uint_t targetRnk = 0; targetRnk < _n_processes; ++targetRnk )
      {
         if ( fill_any )
         {
            // when distributing leftover elements, we have to reset these numbers
            n_vol_max0 = _n_elements / _n_processes;
            n_vol_max1 = n_vol_max0 + std::min( uint_t( 1 ), _n_elements % _n_processes );
         }
         else
         {
            // on the initial run, use the number of currently unassigned elements and empty processes
            auto n_unassigned = _n_elements - n_assigned;
            auto n_empty_proc = _n_processes - targetRnk;
            n_vol_max0        = n_unassigned / n_empty_proc;
            n_vol_max1        = n_vol_max0 + std::min( uint_t( 1 ), n_unassigned % n_empty_proc );
         }

         std::queue< std::vector< K_Simplex* > > Q;
         std::unordered_map< K_Simplex*, bool >  visited;
         // add element to queue
         auto add_to_Q = [&]( K_Simplex* el ) -> bool {
            if ( visited[el] || el->getTargetRank() < _n_processes )
            {
               return false;
            }

            std::vector< K_Simplex* > next;
            if ( !allow_split_siblings && el->has_green_edge() )
            {
               /* elements coming from a green step must reside on the same process
                  to ensure compatibility with further refinement (interpolation between grids)
               */
               for ( auto s : el->get_siblings() )
               {
                  next.push_back( s.get() );
                  visited[s.get()] = true;
               }
            }
            else
            {
               next.push_back( el );
               visited[el] = true;
            }

            Q.push( next );

            return true;
         };

         // assign elements to targetRnk until saturated
         while ( n_vol_on_rnk[targetRnk] < n_vol_max1 && n_assigned < _n_elements )
         {
            if ( Q.empty() && n_assigned < _n_elements )
            {
               // put a random element (+siblings) into the queue
               if ( n_vol_on_rnk[targetRnk] < n_vol_max0 || fill_any )
               {
                  // find the next unassigned element
                  while ( ( *random_element )->getTargetRank() < _n_processes )
                  {
                     ++random_element;
                  }

                  add_to_Q( random_element->get() );
               }
               // don't start a new queue for just one element
               else
               {
                  break;
               }
            }

            // take next element (+siblings) from queue
            auto next = Q.front();
            do
            {
               next = Q.front();
               Q.pop();
            } while ( n_vol_on_rnk[targetRnk] + next.size() > n_vol_max1 && !Q.empty() );

            // assign next (+siblings) to targetRnk
            for ( auto el : next )
            {
               if ( el->getTargetRank() < _n_processes )
               {
                  WALBERLA_ABORT( "something went wrong here" );
               }
               el->setTargetRank( targetRnk );
               ++n_vol_on_rnk[targetRnk];
               ++n_assigned;
            }
            // also change rank of parent s.th. future refinement of parent resides on same process
            if ( next[0]->get_parent() )
               next[0]->get_parent()->setTargetRank( targetRnk );

            // add neighbors to queue
            for ( auto el : next )
            {
               for ( auto nbrID : nbrHood[VOL].at( el->getPrimitiveID() )[VOL] )
               {
                  add_to_Q( id_to_el[nbrID] );
               }
            }
         }
      }
   }
   return n_vol_on_rnk;
}

template < class K_Simplex >
std::shared_ptr< PrimitiveStorage > K_Mesh< K_Simplex >::make_localPrimitives( std::map< PrimitiveID, EdgeData >& edges,
                                                                               std::map< PrimitiveID, FaceData >& faces,
                                                                               std::map< PrimitiveID, CellData >& cells )
{
   auto rank = uint_t( walberla::mpi::MPIManager::instance()->rank() );

   // ****** create maps M: vertexIds -> primitiveID ******
   // identify edges with their vertex indices
   std::map< Idx< 2 >, EdgeData* > vertexIDXsToEdge;
   for ( auto& [_, edge] : edges )
   {
      WALBERLA_ASSERT( _ == edge.getPrimitiveID() );
      vertexIDXsToEdge[edge.get_vertices()] = &edge;
   }
   // identify faces with their vertex IDs
   std::map< Idx< 3 >, FaceData* > vertexIDXsToFace;
   for ( auto& [_, face] : faces )
   {
      WALBERLA_ASSERT( _ == face.getPrimitiveID() );
      vertexIDXsToFace[face.get_vertices()] = &face;
   }

   // ****** find primitives required locally or for halos ******

   walberla::mpi::broadcastObject( _V );

   hyteg::MigrationMap_T nbrRanks;

   for ( auto& [_, vtx] : _V )
   {
      vtx.setLocality( ( vtx.getTargetRank() == rank ) ? LOCAL : NONE );
   }

   for ( auto& [_, edge] : edges )
   {
      WALBERLA_ASSERT( _ == edge.getPrimitiveID() );
      edge.setLocality( ( edge.getTargetRank() == rank ) ? LOCAL : NONE );

      auto& v = edge.get_vertices();

      for ( auto& vtxIdx : v )
      {
         auto& vtx = _V[vtxIdx];
         WALBERLA_ASSERT( vtx.get_vertices()[0] == vtxIdx );

         if ( edge.isLocal() && !vtx.isLocal() )
         {
            vtx.setLocality( HALO );
            nbrRanks[vtx.getPrimitiveID()] = vtx.getTargetRank();
         }
         if ( !edge.isLocal() && vtx.isLocal() )
         {
            edge.setLocality( HALO );
            nbrRanks[edge.getPrimitiveID()] = edge.getTargetRank();
         }
      }
   }

   for ( auto& [_, face] : faces )
   {
      WALBERLA_ASSERT( _ == face.getPrimitiveID() );
      face.setLocality( ( face.getTargetRank() == rank ) ? LOCAL : NONE );

      auto& v = face.get_vertices();

      for ( auto& vtxIdx : v )
      {
         auto& vtx = _V[vtxIdx];
         WALBERLA_ASSERT( vtx.get_vertices()[0] == vtxIdx );

         if ( face.isLocal() && !vtx.isLocal() )
         {
            vtx.setLocality( HALO );
            nbrRanks[vtx.getPrimitiveID()] = vtx.getTargetRank();
         }
         if ( !face.isLocal() && vtx.isLocal() )
         {
            face.setLocality( HALO );
            nbrRanks[face.getPrimitiveID()] = face.getTargetRank();
         }
      }

      for ( uint_t i = 0; i < 3; ++i )
      {
         Idx< 2 >  edgeVtxs{ v[i], v[( i + 1 ) % 3] };
         EdgeData* edge = vertexIDXsToEdge[edgeVtxs];

         if ( face.isLocal() && !edge->isLocal() )
         {
            edge->setLocality( HALO );
            nbrRanks[edge->getPrimitiveID()] = edge->getTargetRank();
         }
         if ( !face.isLocal() && edge->isLocal() )
         {
            face.setLocality( HALO );
            nbrRanks[face.getPrimitiveID()] = face.getTargetRank();
         }
      }
   }

   for ( auto& [_, cell] : cells )
   {
      WALBERLA_ASSERT( _ == cell.getPrimitiveID() );
      cell.setLocality( ( cell.getTargetRank() == rank ) ? LOCAL : NONE );

      auto& v = cell.get_vertices();

      for ( auto& vtxIdx : v )
      {
         auto& vtx = _V[vtxIdx];
         WALBERLA_ASSERT( vtx.get_vertices()[0] == vtxIdx );

         if ( cell.isLocal() && !vtx.isLocal() )
         {
            vtx.setLocality( HALO );
            nbrRanks[vtx.getPrimitiveID()] = vtx.getTargetRank();
         }
         if ( !cell.isLocal() && vtx.isLocal() )
         {
            cell.setLocality( HALO );
            nbrRanks[cell.getPrimitiveID()] = cell.getTargetRank();
         }
      }

      for ( uint_t i = 0; i < 4; ++i )
      {
         for ( uint_t j = i + 1; j < 4; ++j )
         {
            Idx< 2 >  edgeVtxs{ v[i], v[j] };
            EdgeData* edge = vertexIDXsToEdge[edgeVtxs];

            if ( cell.isLocal() && !edge->isLocal() )
            {
               edge->setLocality( HALO );
               nbrRanks[edge->getPrimitiveID()] = edge->getTargetRank();
            }
            if ( !cell.isLocal() && edge->isLocal() )
            {
               cell.setLocality( HALO );
               nbrRanks[cell.getPrimitiveID()] = cell.getTargetRank();
            }
         }
      }

      for ( uint_t j = 0; j < 4; ++j )
      {
         Idx< 3 >  faceVtxs{ v[j], v[( j + 1 ) % 4], v[( j + 2 ) % 4] };
         FaceData* face = vertexIDXsToFace[faceVtxs];

         if ( cell.isLocal() && !face->isLocal() )
         {
            face->setLocality( HALO );
            nbrRanks[face->getPrimitiveID()] = face->getTargetRank();
         }
         if ( !cell.isLocal() && face->isLocal() )
         {
            cell.setLocality( HALO );
            nbrRanks[cell.getPrimitiveID()] = cell.getTargetRank();
         }
      }
   }

   // ****** create primitives ******

   // broadcast vertices to all processes
   walberla::mpi::broadcastObject( _coords );

   // create new vertex and add it to map
   auto add_vertex = [&]( PrimitiveStorage::VertexMap& map, const VertexData& vtx ) {
      // vertex coordinate
      auto coord = vtx.get_coordinates( _coords )[0];
      // add new vertex
      auto& id = vtx.getPrimitiveID();
      map[id]  = std::make_shared< Vertex >( id, coord );
      // add properties
      map[id]->meshBoundaryFlag_ = vtx.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( vtx.getGeometryMap() );
   };
   // create new edge and add it to map
   auto add_edge = [&]( PrimitiveStorage::EdgeMap& map, const EdgeData& edge ) {
      constexpr uint_t K = 1;

      // vertex coordinates and IDs
      auto v      = edge.get_vertices();
      auto coords = edge.get_coordinates( _coords );

      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         auto id      = _V[v[i]].getPrimitiveID();
         vertexIDs[i] = id;
      }

      // add new edge
      auto& id = edge.getPrimitiveID();
      map[id]  = std::make_shared< Edge >( id, vertexIDs[0], vertexIDs[1], coords );

      // add properties
      map[id]->meshBoundaryFlag_ = edge.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( edge.getGeometryMap() );
   };
   // create new face and add it to map
   auto add_face = [&]( PrimitiveStorage::FaceMap& map, const FaceData& face ) {
      constexpr uint_t K = 2;

      // vertex coordinates and IDs
      auto v      = face.get_vertices();
      auto coords = face.get_coordinates( _coords );

      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         auto id      = _V[v[i]].getPrimitiveID();
         vertexIDs[i] = id;
      }

      // ordering of edges
      constexpr std::array< std::array< uint_t, 2 >, K + 1 > edgeOrder{ { { 0ul, 1ul }, { 0ul, 2ul }, { 1ul, 2ul } } };

      // edge IDs and orientation
      std::array< PrimitiveID, K + 1 > edgeIDs;
      std::array< int, K + 1 >         edgeOrientation;
      for ( uint_t i = 0; i <= K; ++i )
      {
         EdgeData* edge = vertexIDXsToEdge[{ v[edgeOrder[i][0]], v[edgeOrder[i][1]] }];

         edgeIDs[i] = edge->getPrimitiveID();

         auto edgeVtx0      = edge->get_vertices()[0];
         edgeOrientation[i] = ( edgeVtx0 == v[edgeOrder[i][0]] ) ? 1 : -1;
      }

      // add new face
      auto& id = face.getPrimitiveID();
      map[id]  = std::make_shared< Face >( id, vertexIDs, edgeIDs, edgeOrientation, coords );

      // add properties
      map[id]->meshBoundaryFlag_ = face.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( face.getGeometryMap() );
   };
   // create new cell and add it to map
   auto add_cell = [&]( PrimitiveStorage::CellMap& map, const CellData& cell ) {
      constexpr uint_t K = 3;

      // vertex coordinates and IDs
      auto v      = cell.get_vertices();
      auto coords = cell.get_coordinates( _coords );

      std::vector< PrimitiveID > vertexIDs( K + 1 );
      for ( uint_t i = 0; i <= K; ++i )
      {
         auto id      = _V[v[i]].getPrimitiveID();
         vertexIDs[i] = id;
      }

      // find cell edges

      constexpr std::array< std::array< uint_t, 2 >, 6 > edgeOrder{
          { { 0ul, 1ul }, { 0ul, 2ul }, { 1ul, 2ul }, { 0ul, 3ul }, { 1ul, 3ul }, { 2ul, 3ul } } };

      std::vector< PrimitiveID >                  edgeIDs( 6 );
      std::array< std::map< uint_t, uint_t >, 6 > edgeLocalVertexToCellLocalVertexMaps;

      for ( uint_t cellLocalEdgeID = 0; cellLocalEdgeID < 6; ++cellLocalEdgeID )
      {
         EdgeData* edge           = vertexIDXsToEdge[{ v[edgeOrder[cellLocalEdgeID][0]], v[edgeOrder[cellLocalEdgeID][1]] }];
         edgeIDs[cellLocalEdgeID] = edge->getPrimitiveID();
         auto& edgeVtxs           = edge->get_vertices();

         for ( const auto& cellLocalVtxID : edgeOrder[cellLocalEdgeID] )
         {
            auto edgeLocalVtxID = uint_t( std::find( edgeVtxs.begin(), edgeVtxs.end(), v[cellLocalVtxID] ) - edgeVtxs.begin() );
            WALBERLA_ASSERT_LESS( edgeLocalVtxID, edgeVtxs.size() );
            edgeLocalVertexToCellLocalVertexMaps[cellLocalEdgeID][edgeLocalVtxID] = cellLocalVtxID;
         }
      }

      // find cell faces

      constexpr std::array< std::array< uint_t, 3 >, K + 1 > faceOrder{
          { { 0ul, 1ul, 2ul }, { 0ul, 1ul, 3ul }, { 0ul, 2ul, 3ul }, { 1ul, 2ul, 3ul } } };

      std::vector< PrimitiveID >                  faceIDs( K + 1 );
      std::array< std::map< uint_t, uint_t >, 4 > faceLocalVertexToCellLocalVertexMaps;

      for ( uint_t cellLocalFaceID = 0; cellLocalFaceID < K + 1; ++cellLocalFaceID )
      {
         FaceData* face           = vertexIDXsToFace[{
             v[faceOrder[cellLocalFaceID][0]], v[faceOrder[cellLocalFaceID][1]], v[faceOrder[cellLocalFaceID][2]] }];
         faceIDs[cellLocalFaceID] = face->getPrimitiveID();
         auto& faceVtxs           = face->get_vertices();

         for ( auto& cellLocalVtxID : faceOrder[cellLocalFaceID] )
         {
            auto faceLocalVtxID = uint_t( std::find( faceVtxs.begin(), faceVtxs.end(), v[cellLocalVtxID] ) - faceVtxs.begin() );
            WALBERLA_ASSERT_LESS( faceLocalVtxID, faceVtxs.size() );
            faceLocalVertexToCellLocalVertexMaps[cellLocalFaceID][faceLocalVtxID] = cellLocalVtxID;
         }
      }

      // add new cell
      auto& id = cell.getPrimitiveID();
      map[id]  = std::make_shared< Cell >(
          id, vertexIDs, edgeIDs, faceIDs, coords, edgeLocalVertexToCellLocalVertexMaps, faceLocalVertexToCellLocalVertexMaps );

      // add properties
      map[id]->meshBoundaryFlag_ = cell.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( cell.getGeometryMap() );
   };

   // create local and halo vertices
   PrimitiveStorage::VertexMap vtxs_ps, nbrVtxs_ps;
   for ( auto& [_, vtx] : _V )
   {
      if ( vtx.isLocal() )
      {
         add_vertex( vtxs_ps, vtx );
      }
      else if ( vtx.onHalo() )
      {
         add_vertex( nbrVtxs_ps, vtx );
      }
   }

   // create local and halo edges
   PrimitiveStorage::EdgeMap edges_ps, nbrEdges_ps;
   for ( auto& [_, edge] : edges )
   {
      WALBERLA_ASSERT( _ == edge.getPrimitiveID() );
      if ( edge.isLocal() )
      {
         add_edge( edges_ps, edge );
      }
      else if ( edge.onHalo() )
      {
         add_edge( nbrEdges_ps, edge );
      }
   }

   // create local and halo faces
   PrimitiveStorage::FaceMap faces_ps, nbrFaces_ps;
   for ( auto& [_, face] : faces )
   {
      WALBERLA_ASSERT( _ == face.getPrimitiveID() );
      if ( face.isLocal() )
      {
         add_face( faces_ps, face );
      }
      else if ( face.onHalo() )
      {
         add_face( nbrFaces_ps, face );
      }
   }

   // create local and halo cells
   PrimitiveStorage::CellMap cells_ps, nbrCells_ps;
   for ( auto& [_, cell] : cells )
   {
      WALBERLA_ASSERT( _ == cell.getPrimitiveID() );
      if ( cell.isLocal() )
      {
         add_cell( cells_ps, cell );
      }
      else if ( cell.onHalo() )
      {
         add_cell( nbrCells_ps, cell );
      }
   }

   // coordinates only required on rank 0
   if ( rank != 0 )
   {
      _coords.clear();
   }

   // ****** add neighborhood information to primitives ******

   // add neighbor edges to vertices
   for ( auto& [id, edge] : edges )
   {
      auto& v = edge.get_vertices();

      for ( auto& vtxIdx : v )
      {
         auto& vtx = _V[vtxIdx];
         WALBERLA_ASSERT( vtx.get_vertices()[0] == vtxIdx );

         if ( vtx.isLocal() )
         {
            vtxs_ps[vtx.getPrimitiveID()]->addEdge( id );
         }
         else if ( vtx.onHalo() )
         {
            nbrVtxs_ps[vtx.getPrimitiveID()]->addEdge( id );
         }
      }
   }

   // add neighbor faces to vertices and edges
   for ( auto& [id, face] : faces )
   {
      auto& v = face.get_vertices();

      for ( auto& vtxIdx : v )
      {
         auto& vtx = _V[vtxIdx];
         WALBERLA_ASSERT( vtx.get_vertices()[0] == vtxIdx );

         if ( vtx.isLocal() )
         {
            vtxs_ps[vtx.getPrimitiveID()]->addFace( id );
         }
         else if ( vtx.onHalo() )
         {
            nbrVtxs_ps[vtx.getPrimitiveID()]->addFace( id );
         }
      }

      for ( uint_t i = 0; i < 3; ++i )
      {
         Idx< 2 >  edgeVtxs{ v[i], v[( i + 1 ) % 3] };
         EdgeData* edge = vertexIDXsToEdge[edgeVtxs];

         if ( edge->isLocal() )
         {
            edges_ps[edge->getPrimitiveID()]->addFace( id );
         }
         else if ( edge->onHalo() )
         {
            nbrEdges_ps[edge->getPrimitiveID()]->addFace( id );
         }
      }
   }

   // add neighbor cells to vertices, edges and faces
   for ( auto& [id, cell] : cells )
   {
      auto& v = cell.get_vertices();

      for ( auto& vtxIdx : v )
      {
         auto& vtx = _V[vtxIdx];
         WALBERLA_ASSERT( vtx.get_vertices()[0] == vtxIdx );

         if ( vtx.isLocal() )
         {
            vtxs_ps[vtx.getPrimitiveID()]->addCell( id );
         }
         else if ( vtx.onHalo() )
         {
            nbrVtxs_ps[vtx.getPrimitiveID()]->addCell( id );
         }
      }

      for ( uint_t i = 0; i < 4; ++i )
      {
         for ( uint_t j = i + 1; j < 4; ++j )
         {
            Idx< 2 >  edgeVtxs{ v[i], v[j] };
            EdgeData* edge = vertexIDXsToEdge[edgeVtxs];

            if ( edge->isLocal() )
            {
               edges_ps[edge->getPrimitiveID()]->addCell( id );
            }
            else if ( edge->onHalo() )
            {
               nbrEdges_ps[edge->getPrimitiveID()]->addCell( id );
            }
         }
      }

      for ( uint_t j = 0; j < 4; ++j )
      {
         Idx< 3 >  faceVtxs{ v[j], v[( j + 1 ) % 4], v[( j + 2 ) % 4] };
         FaceData* face = vertexIDXsToFace[faceVtxs];

         if ( face->isLocal() )
         {
            faces_ps[face->getPrimitiveID()]->addCell( id );
         }
         else if ( face->onHalo() )
         {
            nbrFaces_ps[face->getPrimitiveID()]->addCell( id );
         }
      }
   }

   // vertex data only required on rank 0
   if ( rank != 0 )
   {
      _V.clear();
   }

   // ****** create PrimitiveStorage ******

   return std::make_shared< PrimitiveStorage >(
       vtxs_ps, edges_ps, faces_ps, cells_ps, nbrVtxs_ps, nbrEdges_ps, nbrFaces_ps, nbrCells_ps, nbrRanks, cells.size() > 0 );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::check_integrity( uint_t hanging_nodes_allowed, bool unassigned_nodes_allowed ) const
{
   WALBERLA_DEBUG_SECTION()
   {
      WALBERLA_ROOT_SECTION()
      {
         std::set< std::shared_ptr< K_Simplex > > problematic_elements;

         std::set< uint_t > v_el;
         for ( auto& el : _T )
         {
            if ( el == nullptr )
            {
               WALBERLA_LOG_WARNING( "el == nullptr" );
            }
            else
            {
               if ( hanging_nodes_allowed == 0 && el->vertices_on_edges() > 0 )
               {
                  problematic_elements.insert( el );
               }
               if ( el->has_children() )
               {
                  WALBERLA_LOG_WARNING( "el has children" );
               }
               for ( auto& idx : el->get_vertices() )
               {
                  v_el.insert( idx );
               }
               for ( auto& e : el->get_edges() )
               {
                  if ( e == nullptr )
                  {
                     WALBERLA_LOG_WARNING( "e == nullptr" );
                  }
                  else
                  {
                     if ( hanging_nodes_allowed == 0 )
                     {
                        if ( e->has_children() )
                        {
                           WALBERLA_LOG_WARNING( "e has children" );
                        }
                     }
                     else
                     {
                        if ( e->inner_vertices() > hanging_nodes_allowed )
                        {
                           WALBERLA_LOG_WARNING( "e has too many inner vertices" );
                        }
                     }
                     for ( auto& idx : e->get_vertices() )
                     {
                        if ( !el->has_vertex( idx ) )
                        {
                           WALBERLA_LOG_WARNING( "vertex mismatch element/edge" );
                        }
                     }
                  }
               }
               if constexpr ( VOL == CELL )
               {
                  for ( auto& f : el->get_faces() )
                  {
                     if ( f == nullptr )
                     {
                        WALBERLA_LOG_WARNING( "f == nullptr" );
                     }
                     else
                     {
                        if ( !hanging_nodes_allowed )
                        {
                           if ( f->has_children() )
                           {
                              WALBERLA_LOG_WARNING( "f has children" );
                           }
                        }
                        for ( auto& idx : f->get_vertices() )
                        {
                           if ( !el->has_vertex( idx ) )
                           {
                              WALBERLA_LOG_WARNING( "vertex mismatch element/face" );
                           }
                        }
                     }
                  }
               }
            }
         }

         if ( !unassigned_nodes_allowed )
         {
            std::vector< uint_t > v_msh, diff;
            for ( auto& [idx, _] : _V )
            {
               WALBERLA_UNUSED( _ );
               v_msh.push_back( idx );
            }
            std::set_difference( v_msh.begin(), v_msh.end(), v_el.begin(), v_el.end(), std::back_inserter( diff ) );
            if ( !diff.empty() )
            {
               WALBERLA_LOG_WARNING( "vertices in mesh but not in elements" );
               for ( auto& idx : diff )
               {
                  WALBERLA_LOG_WARNING( idx );
               }
            }
            diff.clear();
            std::set_difference( v_el.begin(), v_el.end(), v_msh.begin(), v_msh.end(), std::back_inserter( diff ) );
            if ( !diff.empty() )
            {
               WALBERLA_LOG_WARNING( "vertices in elements but not in mesh" );
               for ( auto& idx : diff )
               {
                  WALBERLA_LOG_WARNING( idx );
               }
            }
         }

         for ( auto& el : problematic_elements )
         {
            WALBERLA_LOG_WARNING( "hanging nodes: " << el->vertices_on_edges() );
         }
      }
   }
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::extract_data( std::map< PrimitiveID, EdgeData >& edgeData,
                                        std::map< PrimitiveID, FaceData >& faceData,
                                        std::map< PrimitiveID, CellData >& cellData ) const
{
   if ( walberla::mpi::MPIManager::instance()->rank() != 0 )
      return;

   std::set< std::shared_ptr< Simplex1 > > edges;
   std::set< std::shared_ptr< Simplex2 > > faces;
   std::set< std::shared_ptr< Simplex3 > > cells;

   // collect cells
   if constexpr ( VOL == CELL )
   {
      cells = _T;
   }
   // collect faces
   if constexpr ( VOL == FACE )
   {
      faces = _T;
   }
   for ( auto& cell : cells )
   {
      for ( auto& face : cell->get_faces() )
      {
         faces.insert( face );
      }
   }
   // collect edges
   for ( auto& face : faces )
   {
      for ( auto& edge : face->get_edges() )
      {
         edges.insert( edge );
      }
   }

   // collect edge data
   for ( auto& edge : edges )
   {
      edgeData[edge->getPrimitiveID()] = EdgeData( edge.get() );
   }
   // collect face data
   for ( auto& face : faces )
   {
      faceData[face->getPrimitiveID()] = FaceData( face.get() );
   }
   // collect cell data
   for ( auto& cell : cells )
   {
      cellData[cell->getPrimitiveID()] = CellData( cell.get() );
   }
}

template < class K_Simplex >
inline std::pair< std::set< std::shared_ptr< K_Simplex > >, std::set< std::shared_ptr< K_Simplex > > >
    K_Mesh< K_Simplex >::init_R_Pc( const std::vector< PrimitiveID >& id_r, const std::vector< PrimitiveID >& id_c ) const
{
   std::map< PrimitiveID, std::shared_ptr< K_Simplex > > idToEl;
   for ( auto& el : _T )
   {
      idToEl[el->getPrimitiveID()] = el;
   }

   std::set< std::shared_ptr< K_Simplex > > R, Pc;

   // fill R
   for ( auto& id : id_r )
   {
      auto el = idToEl[id];
      if ( el->has_green_edge() )
      {
         R.insert( el->get_parent() );
      }
      else
      {
         R.insert( el );
      }
   }

   // the children of p must be kept if any of them are marked for refinement
   auto keep_children = [&]( std::shared_ptr< K_Simplex > p ) {
      for ( auto& c : p->get_children() )
      {
         if ( R.count( c ) > 0 )
         {
            return true;
         }
      }
      return false;
   };

   // fill Pc
   for ( auto& id : id_c )
   {
      auto el = idToEl[id];
      auto p  = el->get_parent();

      // el is already a root element
      if ( p == nullptr )
      {
         continue;
      }
      // since green elements are being removed anyway, we ignore them here
      if ( el->has_green_edge() )
      {
         continue;
      }
      // el's parent can't be un-refined if any of its children are marked for refinement
      if ( keep_children( p ) )
      {
         continue;
      }

      Pc.insert( p );
   }

   return { R, Pc };
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::unrefine( const std::set< std::shared_ptr< K_Simplex > >& Pc )
{
   std::set< std::shared_ptr< Simplex1 > > edges_to_unrefine;
   std::set< std::shared_ptr< Simplex2 > > faces_to_unrefine;

   /* mark faces/edges of elements in T if they have children
      This should only be the case if an element has been green refined before
   */
   for ( auto& el : _T )
   {
      for ( auto& edge : el->get_edges() )
      {
         if ( edge->has_children() )
         {
            edges_to_unrefine.insert( edge );
         }
      }
      if constexpr ( VOL == CELL )
      {
         for ( auto& face : el->get_faces() )
         {
            if ( face->has_children() )
            {
               faces_to_unrefine.insert( face );
            }
         }
      }
   }

   // repeatedly iterate over Pc and un-refine admissable elements until no admissable elements are left
   auto repeat = true;
   auto queue  = Pc;
   while ( repeat )
   {
      repeat = false;
      for ( auto it = queue.begin(); it != queue.end(); )
      {
         auto el = *it;

         if ( el->has_grandkids() )
         {
            // if el has grandkids, don't process it until they are removed
            ++it;
         }
         else
         {
            // replace el's children in T with el
            for ( auto& child : el->get_children() )
            {
               _T.erase( child );
            }
            el->kill_children();
            _T.insert( el );

            // mark el's faces and edges for un-refinement
            for ( auto& edge : el->get_edges() )
            {
               edges_to_unrefine.insert( edge );
            }
            if constexpr ( VOL == CELL )
            {
               for ( auto& face : el->get_faces() )
               {
                  faces_to_unrefine.insert( face );
               }
            }

            // we need to check all elements again to potentially unrefine el's parent
            repeat = true;

            it = queue.erase( it );
         }
      }
   }

   // check whether marked edges/faces can be unrefined
   for ( auto& el : _T )
   {
      /*if a face/edge belongs to an ancestor of an el, some descendent of
         this face/edge belongs to el. Thus, we can't unrefine the face/edge.
      */
      for ( auto p = el->get_parent(); p != nullptr; p = p->get_parent() )
      {
         for ( auto& edge : p->get_edges() )
         {
            edges_to_unrefine.erase( edge );
         }
         if constexpr ( VOL == CELL )
         {
            for ( auto& face : p->get_faces() )
            {
               faces_to_unrefine.erase( face );
            }
         }
      }
   }

   // remove obsolete vertices, edges and faces
   // vertices must be processed recursively to make sure that all vertices are removed
   std::function< void( std::shared_ptr< Simplex1 > ) > remove_inner_vtxs = [&]( std::shared_ptr< Simplex1 > edge ) {
      for ( auto& child : edge->get_children() )
      {
         remove_inner_vtxs( child );
      }
      if ( edge->has_children() )
      {
         auto vtx = edge->get_midpoint_idx();
         _V.erase( vtx );
         _coords.erase( vtx );
      }
   };
   for ( auto& edge : edges_to_unrefine )
   {
      remove_inner_vtxs( edge );
      edge->reset_midpoint_idx();
      edge->kill_children();
   }
   for ( auto& face : faces_to_unrefine )
   {
      face->kill_children();
   }

   check_integrity( 100 );
}

template < class K_Simplex >
bool K_Mesh< K_Simplex >::refine_faces()
{
   bool faces_refined = false;
   if constexpr ( VOL == CELL )
   {
      for ( auto& el : _T )
      {
         for ( auto& face : el->get_faces() )
         {
            if ( face->vertices_on_edges() > 1 && !face->has_children() )
            {
               refine_face_red( _coords, _V, face );
               faces_refined = true;
            }
         }
      }
   }
   return faces_refined;
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::refine_red( const std::set< std::shared_ptr< K_Simplex > >& R )
{
   for ( auto& el : R )
   {
      if ( _T.erase( el ) == 0 )
      {
         // for el ∉ T: don't try to refine
         continue;
      }

      if constexpr ( VOL == FACE )
      {
         _T.merge( refine_face_red( _coords, _V, el ) );
      }
      if constexpr ( VOL == CELL )
      {
         _T.merge( refine_cell_red( _coords, _V, el ) );
      }
   }

   // In 3d, we additionally refine all faces with more than one hanging node
   bool rf = true;
   while ( rf )
   {
      rf = refine_faces();
   }
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::remove_green_edges()
{
   auto T_cpy = _T;

   for ( auto& el : T_cpy )
   {
      if ( el->has_green_edge() )
      {
         // replace el with its parent in T
         _T.erase( el );
         _T.insert( el->get_parent() );

         // remove green child-elements
         el->get_parent()->kill_children();
         if constexpr ( VOL == CELL )
         {
            // remove green edges from faces
            for ( auto& face : el->get_parent()->get_faces() )
            {
               if ( face->get_children().size() == 2 )
               {
                  face->kill_children();
               }
            }
         }
      }
   }
}

template < class K_Simplex >
std::set< std::shared_ptr< K_Simplex > > K_Mesh< K_Simplex >::find_elements_for_red_refinement() const
{
   std::set< std::shared_ptr< K_Simplex > > R;

   auto mark_for_refinement = [&]( std::shared_ptr< K_Simplex > el ) -> bool {
      if constexpr ( VOL == FACE )
      {
         // mark face for red refinement if it has more than one hanging node
         return el->vertices_on_edges() > 1;
      }
      if constexpr ( VOL == CELL )
      {
         /* mark cell for red refinement if it has more than 3 hanging nodes
            This is the case if either of the following conditions is met:
            (a) more than one of the cell's faces are refined (red)
            (b) there are hanging nodes on the interior of one of its faces
         */
         uint_t rf = 0; // number of refined faces
         for ( auto& face : el->get_faces() )
         {
            if ( face->has_children() )
            {
               ++rf;
               // check for condition (a)
               if ( rf > 1 )
               {
                  return true;
               }
               /* check for condition (b)
                  If there is a hanging node on the interior of the face, this node must be
                  on the boundary of the interior child (face->get_children()[3]).
               */
               if ( face->get_children()[3]->vertices_on_edges() > 0 )
               {
                  return true;
               }
            }
         }
      }
      return false;
   };

   for ( auto& el : _T )
   {
      if ( mark_for_refinement( el ) )
      {
         R.insert( el );
      }
   }

   return R;
}

template <>
void K_Mesh< Simplex2 >::refine_green()
{
   check_integrity( 1 );

   auto U = _T;

   for ( auto& el : U )
   {
      uint_t hanging_nodes = el->vertices_on_edges();

      if ( hanging_nodes > 0 )
      {
         WALBERLA_ASSERT( hanging_nodes == 1 );
         WALBERLA_ASSERT( !el->has_green_edge() );

         _T.erase( el );
         _T.merge( refine_face_green( el ) );
      }
   }

   check_integrity( 0 );
}

template <>
void K_Mesh< Simplex3 >::refine_green()
{
   check_integrity( 1, true );

   auto U = _T;

   for ( auto& el : U )
   {
      uint_t hanging_nodes = el->vertices_on_edges();

      if ( hanging_nodes == 0 )
      {
         continue;
      }

      std::set< std::shared_ptr< Simplex3 > > new_els;

      if ( hanging_nodes == 1 )
      {
         new_els = refine_cell_green_1( el );
      }
      else if ( hanging_nodes == 2 )
      {
         new_els = refine_cell_green_2( el );
      }
      else if ( hanging_nodes == 3 )
      {
         new_els = refine_cell_green_3( el );
      }
      else
      {
         WALBERLA_ASSERT( hanging_nodes <= 3 );
      }

      _T.erase( el );
      _T.merge( new_els );
   }
   check_integrity( 0 );
}

template < class K_Simplex >
std::pair< real_t, real_t > K_Mesh< K_Simplex >::min_max_angle() const
{
   std::pair< real_t, real_t > mm{ 10, 0 };

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         auto mm_el = el->min_max_angle( _coords );

         mm.first  = std::min( mm.first, mm_el.first );
         mm.second = std::max( mm.second, mm_el.second );
      }
   }

   walberla::mpi::broadcastObject( mm );

   return mm;
}

template < class K_Simplex >
std::pair< real_t, real_t > K_Mesh< K_Simplex >::mean_min_max_angle() const
{
   std::pair< real_t, real_t > mm{ 0, 0 };

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         auto mm_el = el->min_max_angle( _coords );

         mm.first += mm_el.first;
         mm.second += mm_el.second;
      }
   }

   mm.first /= real_t( n_elements() );
   mm.second /= real_t( n_elements() );

   walberla::mpi::broadcastObject( mm );

   return mm;
}

template < class K_Simplex >
std::pair< real_t, real_t > K_Mesh< K_Simplex >::min_max_volume() const
{
   std::pair< real_t, real_t > mm{ std::numeric_limits< real_t >::max(), 0 };

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         auto v = el->volume( _coords );

         mm.first  = std::min( mm.first, v );
         mm.second = std::max( mm.second, v );
      }
   }

   walberla::mpi::broadcastObject( mm );

   return mm;
}

template < class K_Simplex >
real_t K_Mesh< K_Simplex >::volume() const
{
   real_t v_tot = 0;

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         v_tot += el->volume( _coords );
      }
   }

   walberla::mpi::broadcastObject( v_tot );

   return v_tot;
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::exportMesh( const std::string& filename ) const
{
   WALBERLA_ROOT_SECTION()
   {
      int elType;
      if ( VOL == 2 )
         elType = 2;
      else if ( VOL == 3 )
         elType = 4;
      else
         WALBERLA_ABORT( "Only implementation of 3d and 2d meshes implemented." )

      std::ofstream file( filename );

      file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << n_vtx() << "\n";

      // due to repeated coarsening/refinement our indices may be non-contiguous
      std::map< uint_t, uint_t > vtxIdxMap;

      uint_t i = 0;
      for ( auto& [idx, coord] : _coords )
      {
         ++i;
         vtxIdxMap[idx] = i;

         file << i;
         for ( int j = 0; j < 3; ++j )
            file << " " << coord[j];
         file << "\n";
      }

      file << "$EndNodes\n$Elements\n" << n_elements() << "\n";

      uint_t k = 0;
      for ( auto& el : _T )
      {
         ++k;
         auto& v = el->get_vertices();

         file << k << " " << elType << " 2 0 0";
         for ( uint_t j = 0; j <= VOL; ++j )
            file << " " << vtxIdxMap[v[j]];
         file << "\n";
      }

      file << "$EndElements\n";
   }
}

template class K_Mesh< Simplex2 >;
template class K_Mesh< Simplex3 >;

} // namespace adaptiveRefinement
} // namespace hyteg
