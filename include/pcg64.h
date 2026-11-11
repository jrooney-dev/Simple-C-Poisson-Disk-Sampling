#ifndef PCG64_H
#define PCG64_H

#include <stdint.h>

typedef struct
{
    __uint128_t state;
    __uint128_t inc;
} pcg64_t;

/* Function prototypes */
void pcg64_seed(pcg64_t *rng, uint64_t seed, uint64_t stream);
uint64_t pcg64_random(pcg64_t *rng);
double pcg64_double(pcg64_t *rng);
uint64_t pcg64_uint64(pcg64_t *rng);
uint64_t pcg64_range(pcg64_t *rng, uint64_t min, uint64_t max);
void pcg64_generate(pcg64_t *rng, uint64_t min, uint64_t max, int count);

#endif // FUNCTION_H