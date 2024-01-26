#include <immintrin.h>

int main()
{
   const __m256d result = _mm256_mul_pd( _mm256_set_pd( 4.0, 5.0, 6.0, 7.0 ), _mm256_set_pd( 5.0, 6.0, 7.0, 8.0 ) );
}
