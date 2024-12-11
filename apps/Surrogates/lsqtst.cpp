#include <core/Environment.h>
#include <core/Format.hpp>
#include <core/config/Create.h>
#include <core/timing/Timer.h>
#include <hyteg/polynomial/elementwise/leastSquares.hpp>

template < uint8_t dim, uint8_t degree, uint8_t lvl, bool reduced >
void lsqtst( int print_iterator )
{
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "LeastSquares< dim=%d, degree=%d, lvl=%d, reduced=%d >", dim, degree, lvl, reduced ) );

   auto t0 = walberla::timing::getWcTime();

   hyteg::surrogate::LeastSquares< dim, degree, lvl, reduced, false > lsq;

   auto t1 = walberla::timing::getWcTime();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Time for Vander+SVD: %f", t1 - t0 ) );

   // todo sample polynomial

   if ( print_iterator == 1 )
   {
      auto it = lsq.samplingIterator();

      while ( it.n < lsq.rows )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "i = %d, j = %d, k = %d", it.i, it.j, it.k ) );
         ++it;
      }
   }
   if ( print_iterator == 2 )
   {
      auto it = lsq.samplingIterator();

      walberla::uint_t i = 0, j = 0, k = 0;
      while ( it.n < lsq.rows )
      {
         if ( k < it.k )
         {
            std::cout << "\n\n";
            k = it.k;
            j = 0;
            i = 0;
         }
         for ( ; j < it.j; ++j )
         {
            for ( i = 0; i < ( 1 << lvl ); ++i )
            {
               std::cout << " ";
            }
            std::cout << "\n";
            i = 0;
         }
         if ( it.i > 0 )
         {
            ++i;
         }
         for ( ; i < it.i; ++i )
         {
            std::cout << " ";
         }
         std::cout << "X";

         ++it;
      }
   }
}

template < uint8_t D, uint8_t Q >
void polytst()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Polynomial< D=%d, Q=%d >", D, Q ) );
   // setup polynomial
   std::array< double, hyteg::surrogate::polynomial::dimP< D, Q > > c;
   for ( uint32_t i = 0; i < c.size(); ++i )
   {
      c[i] = std::sin( double( i ) );
   }
   hyteg::surrogate::polynomial::Polynomial< D, Q > p( c );

   // coordinates
   std::array< double, 3 > x{ 1e-2, 1e-1, 3.0 };

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
   auto                        ref = 0.0;
   for ( uint32_t ijk = 0; ijk < c.size(); ++ijk )
   {
      auto [i, j, k] = p.basis( ijk );
      ref += c[ijk] * x_pow[i] * y_pow[j] * z_pow[k];
   }
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "relative error; split %e; naive %e", ( split - ref ) / ref, ( naive - ref ) / ref ) );
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

   // lsqtst< 2, 2, 2, false >( 1 );
   // lsqtst< 3, 2, 2, false >( 1 );

   // lsqtst< 3, 6, 3, false >( 0 );
   // lsqtst< 3, 6, 4, false >( 0 );
   // lsqtst< 3, 6, 5, false >( 0 );

   // lsqtst< 3, 8, 3, false >( 0 );
   // lsqtst< 3, 8, 4, false >( 0 );
   // lsqtst< 3, 8, 5, false >( 0 );

   // lsqtst< 3, 10, 3, false >( 0 );
   // lsqtst< 3, 10, 4, false >( 0 );
   // lsqtst< 3, 3, 6, false >( 2 );
   // lsqtst< 3, 3, 6, true >( 2 );

   // lsqtst< 3, 10, 5, true >( 0 );
   // lsqtst< 3, 10, 8, true >( 0 );

   polytst< 1, 1 >();
   polytst< 1, 2 >();
   polytst< 1, 3 >();
   polytst< 1, 5 >();
   polytst< 1, 10 >();
   polytst< 1, 15 >();

   polytst< 2, 1 >();
   polytst< 2, 2 >();
   polytst< 2, 3 >();
   polytst< 2, 5 >();
   polytst< 2, 10 >();
   polytst< 2, 15 >();

   polytst< 3, 1 >();
   polytst< 3, 2 >();
   polytst< 3, 3 >();
   polytst< 3, 5 >();
   polytst< 3, 15 >();
   return 0;
}