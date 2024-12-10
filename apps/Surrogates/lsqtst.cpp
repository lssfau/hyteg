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

   lsqtst< 2, 2, 2, false >( 1 );
   lsqtst< 3, 2, 2, false >( 1 );

   lsqtst< 3, 6, 3, false >( 0 );
   lsqtst< 3, 6, 4, false >( 0 );
   lsqtst< 3, 6, 5, false >( 0 );

   lsqtst< 3, 8, 3, false >( 0 );
   lsqtst< 3, 8, 4, false >( 0 );
   lsqtst< 3, 8, 5, false >( 0 );

   lsqtst< 3, 10, 3, false >( 0 );
   lsqtst< 3, 10, 4, false >( 0 );
   lsqtst< 3, 3, 6, false >( 2 );
   lsqtst< 3, 3, 6, true >( 2 );

   lsqtst< 3, 10, 5, true >( 0 );
   lsqtst< 3, 10, 8, true >( 0 );

   return 0;
}