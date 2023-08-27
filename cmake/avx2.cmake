# The purpose of this if to check the availability of avx2 ob the system,
# define the macro P4EST_ENABLE_AVX2 if necesessary.
include(CheckCSourceRuns)

if(NOT CMAKE_CROSSCOMPILING)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -mavx2")
  check_c_source_runs("
    #include <immintrin.h>
    int main (void) {
      __m256i a;
      a = _mm256_set1_epi32 (1);
      _mm256_abs_epi32 (a);
      return 0;
    }"
    P4EST_ENABLE_AVX2)
    set(CMAKE_REQUIRED_FLAGS)
endif()
