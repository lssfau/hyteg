#include <immintrin.h>

int main()
{
   // avx
   const __m256d resultd = _mm256_mul_pd( _mm256_set_pd( 4.0, 5.0, 6.0, 7.0 ), _mm256_set_pd( 5.0, 6.0, 7.0, 8.0 ) );
   // avx2
   const __m256i resulti = _mm256_and_si256( _mm256_set1_epi64x( 6 ), _mm256_set1_epi64x( 7 ) );
}
