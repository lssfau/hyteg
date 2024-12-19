// todo make this a test
#include <core/Environment.h>
#include <core/Format.hpp>
#include <core/config/Create.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <core/timing/Timer.h>
#include <filesystem>
#include <hyteg/polynomial/elementwise/leastSquares.hpp>

using walberla::uint_t;
using walberla::math::pi;
using walberla::math::realRandom;

template < uint8_t dim, uint8_t degree >
void lsqtst( uint_t lvl, uint_t downsampling, int print_iterator, bool use_precomputed, bool save_svd, bool delete_svd )
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "LeastSquares: dim=%d, degree=%d, lvl=%d, downsampling=%d, use_precomputed=%d",
                                                dim,
                                                degree,
                                                lvl,
                                                downsampling,
                                                use_precomputed ) );

   std::unique_ptr< hyteg::surrogate::LeastSquares< dim, degree > > lsq;
   std::string path_to_svd = walberla::format( "svd_%d_%d_%d_%d", dim, degree, lvl, downsampling );

   auto t0 = walberla::timing::getWcTime();

   if ( use_precomputed )
   {
      lsq = std::make_unique< hyteg::surrogate::LeastSquares< dim, degree > >( path_to_svd, lvl, downsampling );
   }
   else
   {
      lsq = std::make_unique< hyteg::surrogate::LeastSquares< dim, degree > >( lvl, downsampling );
   }

   auto t1 = walberla::timing::getWcTime();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "n_samples = %d; dim(P) = %d", lsq->rows, lsq->cols ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Time for computing/loading Vander+SVD: %f", t1 - t0 ) );

   if ( !use_precomputed && save_svd )
   {
      if ( !std::filesystem::exists( path_to_svd ) )
      {
         std::filesystem::create_directory( path_to_svd );
      }
      t0 = walberla::timing::getWcTime();
      lsq->write_to_file( path_to_svd );
      t1 = walberla::timing::getWcTime();
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Time for storing Vander+SVD: %f", t1 - t0 ) );
   }
   if ( use_precomputed && delete_svd )
   {
      if ( std::filesystem::exists( path_to_svd ) )
      {
         std::filesystem::remove_all( path_to_svd );
      }
   }

   // define smooth function
   auto f = []( const std::array< double, 3 >& x ) {
      return std::sin( pi / 2.0 * x[0] ) * std::cos( pi / 2.0 * x[1] ) * std::sin( pi / 2.0 * x[2] );
   };
   // coordinates on tetrahedron
   hyteg::surrogate::polynomial::Coordinates coords( lvl );
   // fill right-hand side
   auto it = lsq->samplingIterator();
   while ( it != it.end() )
   {
      auto x = coords.x( it.i(), it.j(), it.k() );
      lsq->setRHS( it(), f( x ) );
      ++it;
   }
   // solve least squares problem
   auto p = lsq->solve();
   // evaluate polynomial and compute error
   uint_t len_edge = ( uint_t( 1 ) << lvl );
   double error    = 0.0;
   uint_t k_max    = ( dim == 3 ) ? len_edge : 1;
   for ( uint_t k = 0; k < k_max; ++k )
   {
      auto z = coords.x( k );
      if ( dim == 3 )
      {
         p.fix_z( z );
      }
      for ( uint_t j = 0; j < len_edge - k; ++j )
      {
         auto y = coords.x( j );
         if ( dim >= 2 )
         {
            p.fix_y( y );
         }
         for ( uint_t i = 0; i < len_edge - k - j; ++i )
         {
            auto x = coords.x( i );
            auto e = f( { x, y, z } ) - p.eval( x );
            error += e * e;
         }
      }
   }
   auto h  = 1.0 / double( len_edge );
   auto dV = std::pow( h, dim );
   error   = std::sqrt( error * dV );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Error: ||f-p||_0 = %e", std::sqrt( error ) ) );

   if ( print_iterator == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Sample points:" );
      it = lsq->samplingIterator();

      while ( it != it.end() )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "i = %d, j = %d, k = %d", it.i(), it.j(), it.k() ) );
         ++it;
      }
   }
   if ( print_iterator == 2 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Sample points:" );
      it = lsq->samplingIterator();

      walberla::uint_t i = 0, j = 0, k = 0;
      std::cout << "\n";
      while ( it != it.end() )
      {
         if ( k < it.k() )
         {
            std::cout << "\n\n";
            for ( i = 0; i < len_edge; ++i )
            {
               std::cout << "-";
            }
            std::cout << "\n\n";
            k = it.k();
            j = 0;
            i = 0;
         }
         for ( ; j < it.j(); ++j )
         {
            for ( i = 0; i < len_edge; ++i )
            {
               std::cout << " ";
            }
            std::cout << "\n";
            i = 0;
         }
         if ( it.i() > 0 )
         {
            ++i;
         }
         for ( ; i < it.i(); ++i )
         {
            std::cout << " ";
         }
         std::cout << "X";

         ++it;
      }
      std::cout << "\n\n";
      for ( i = 0; i < len_edge; ++i )
      {
         std::cout << "-";
      }
      std::cout << "\n\n";
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

template < uint8_t D, uint8_t Q >
void polytst()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Polynomial< D=%d, Q=%d >", D, Q ) );
   // setup polynomial
   std::array< double, hyteg::surrogate::polynomial::dimP< D, Q > > c;
   for ( uint32_t i = 0; i < c.size(); ++i )
   {
      c[i] = realRandom();
   }
   hyteg::surrogate::polynomial::Polynomial< D, Q > p( c );

   // coordinates
   std::array< double, 3 > x{ realRandom(), realRandom(), realRandom() };

   // evaluate polynomial using Horner's method
   if constexpr ( D == 3 )
   {
      p.fix_z( x[2] );
   }
   if constexpr ( D >= 2 )
   {
      p.fix_y( x[1] );
   }
   auto split = p.eval( x[0] );

   // evaluate polynomial by summing up basis functions
   auto naive = p.eval_naive( x );

   // reference evaluation
   std::array< double, Q + 1 > x_pow{ 1.0 };
   std::array< double, Q + 1 > y_pow{ 1.0 };
   std::array< double, Q + 1 > z_pow{ 1.0 };
   for ( uint8_t q = 0; q < Q; ++q )
   {
      x_pow[q + 1] = x_pow[q] * x[0];
      y_pow[q + 1] = y_pow[q] * x[1];
      z_pow[q + 1] = z_pow[q] * x[2];
   }
   auto ref = 0.0;
   for ( uint32_t ijk = 0; ijk < c.size(); ++ijk )
   {
      auto [i, j, k] = p.basis( ijk );
      ref += c[ijk] * x_pow[i] * y_pow[j] * z_pow[k];
   }
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "relative error; split %e; naive %e", ( split - ref ) / ref, ( naive - ref ) / ref ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   //    walberla::shared_ptr< walberla::config::Config > cfg;

   //    if ( argc == 1 )
   //    {
   //       walberla::shared_ptr< walberla::config::Config > cfg_( new walberla::config::Config );
   //       // cfg_->readParameterFile("../../data/param/PolynomialBlending.prm");
   //       cfg_->readParameterFile( "../../hyteg/data/param/PolynomialBlending.prm" );
   //       cfg = cfg_;
   //    }
   //    else
   //    {
   //       cfg = walberla::config::create( argc, argv );
   //    }
   lsqtst< 2, 5 >( 4, 1, 0, 0, 1, 0 );
   lsqtst< 2, 5 >( 4, 1, 0, 1, 0, 1 );
   lsqtst< 3, 5 >( 4, 1, 0, 0, 1, 0 );
   lsqtst< 3, 5 >( 4, 1, 0, 1, 0, 1 );
   lsqtst< 3, 10 >( 7, 1, 0, 1, 0, 0 );

   // polytst< 3, 15 >();

   return 0;
}