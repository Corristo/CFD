#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <immintrin.h>
#include "helper.h"


void computeDensity(const double *const currentCell, double *density, const int tot_cells)
{
    *density = 0.0;
    int i;
    for (i = 0; i < PARAMQ; i++)
      *density += *(currentCell + i * tot_cells);
}

void computeDensityAVX(const double * const currentCell, double *density, const int tot_cells)
{
    const double zero = 0.0;
    __m256d vsum = _mm256_broadcast_sd(&zero);
    int i;
    for (i = 0; i < PARAMQ; i ++)
    {
        __m256d v = _mm256_loadu_pd(currentCell + i * tot_cells);
        vsum = _mm256_add_pd(vsum, v);
    }
   _mm256_storeu_pd(density, vsum);
}

void computeDensitySSE(const double * const currentCell, double *density)
{
    __m128d vsum = _mm_set1_pd(0.0);
    int i;
    for (i = 0; i < PARAMQ - 1; i += 2)
    {
      __m128d v = _mm_loadu_pd(&currentCell[i]);

        vsum = _mm_add_pd(vsum, v);
    }
    vsum = _mm_hadd_pd(vsum, vsum);
    _mm_storeh_pd(density, vsum);
    if (i < PARAMQ)
    {
      *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity, const int tot_cells)
{
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;
    for (int i = 0; i < PARAMQ; i++)
    {
        velocity[0] += *(currentCell + i * tot_cells) * LATTICEVELOCITIES[i][0];

        velocity[1] += *(currentCell + i * tot_cells) * LATTICEVELOCITIES[i][1];

        velocity[2] += *(currentCell + i * tot_cells) * LATTICEVELOCITIES[i][2];
    }

    velocity[0] = velocity[0] / (*density);
    velocity[1] = velocity[1] / (*density);
    velocity[2] = velocity[2] / (*density);

}

void computeVelocityAVX(const double * const currentCell, const double * const density, double velocity[3][4], const int tot_cells)
{
    __m256d v0, v1, v2;
    __m256d vc, vl0, vl1, vl2;
    __m128i vtemp;
    int i;
    const double zero = 0.0;
    v0 = v1 = v2 = _mm256_broadcast_sd(&zero);

    for (i = 0; i < PARAMQ; i++)
    {
        vc = _mm256_loadu_pd(currentCell + i * tot_cells);
        vtemp = _mm_set1_epi32(LATTICEVELOCITIES[i][0]);
        vl0 = _mm256_cvtepi32_pd(vtemp);
        vtemp = _mm_set1_epi32(LATTICEVELOCITIES[i][1]);
        vl1 = _mm256_cvtepi32_pd(vtemp);
        vtemp = _mm_set1_epi32(LATTICEVELOCITIES[i][2]);
        vl2 = _mm256_cvtepi32_pd(vtemp);
        v0 = _mm256_add_pd(v0, _mm256_mul_pd(vc, vl0));
        v1 = _mm256_add_pd(v1, _mm256_mul_pd(vc, vl1));
        v2 = _mm256_add_pd(v2, _mm256_mul_pd(vc, vl2));
    }
    vc = _mm256_loadu_pd(density);
    v0 = _mm256_div_pd(v0, vc);
    v1 = _mm256_div_pd(v1, vc);
    v2 = _mm256_div_pd(v2, vc);

    _mm256_storeu_pd( velocity[0], v0);
    _mm256_storeu_pd( velocity[1], v1);
    _mm256_storeu_pd( velocity[2], v2);

}
void computeVelocitySSE(const double * const currentCell, const double * const density, double *velocity)
{
    __m128d v0, v1, v2;
    int i;
    v0 = v1 = v2 = _mm_setzero_pd();
    for (i = 0; i < PARAMQ - 1; i += 2)
    {
        __m128d vc, vl0, vl1, vl2;
        __m128i vtemp;

        vc = _mm_loadu_pd(&currentCell[i]);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[0][i]);
        vl0 = _mm_cvtepi32_pd(vtemp);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[1][i]);
        vl1 = _mm_cvtepi32_pd(vtemp);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[2][i]);
        vl2 = _mm_cvtepi32_pd(vtemp);
        v0 = _mm_add_pd(v0, _mm_mul_pd(vc, vl0));
        v1 = _mm_add_pd(v1, _mm_mul_pd(vc, vl1));
        v2 = _mm_add_pd(v2, _mm_mul_pd(vc, vl2));
    }
    v0 = _mm_hadd_pd(v0, v0);
    v1 = _mm_hadd_pd(v1, v1);
    v2 = _mm_hadd_pd(v2, v2);
    _mm_store_sd (&velocity[0], v0);
    _mm_store_sd (&velocity[1], v1);
    _mm_store_sd (&velocity[2], v2);
    if (i < PARAMQ)
    {
        velocity[0] += currentCell[i] * LATTICEVELOCITIES2[0][i];
        velocity[1] += currentCell[i] * LATTICEVELOCITIES2[1][i];
        velocity[2] += currentCell[i] * LATTICEVELOCITIES2[2][i];
    }
    velocity[0] = velocity[0] / (*density);
    velocity[1] = velocity[1] / (*density);
    velocity[2] = velocity[2] / (*density);
}
void computeFeq(const double * const density, const double * const velocity, double *feq)
{
    double velocityMagnitudeSqr, dotProduct;
    double c_s_square = C_S * C_S;
    velocityMagnitudeSqr = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    for (int i = 0; i < PARAMQ; i++)
    {
        dotProduct = velocity[0] * LATTICEVELOCITIES[i][0] + velocity[1] * LATTICEVELOCITIES[i][1] + velocity[2] * LATTICEVELOCITIES[i][2];
        feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + dotProduct / c_s_square + dotProduct * dotProduct / (2  * c_s_square * c_s_square) - velocityMagnitudeSqr / (2 * c_s_square));
    }
}

void computeFeqAVXv2(double density[4], double velocity[3][4], double feq[PARAMQ][4])
{
    __m256d latticeWeightsVector, latticeVelocityX, latticeVelocityY, latticeVelocityZ, dotProductVector, result;
    __m128i tmp;
    __m256d densityVector = _mm256_loadu_pd(density);
    __m256d velocityX = _mm256_loadu_pd(velocity[0]);
    __m256d velocityY = _mm256_loadu_pd(velocity[1]);
    __m256d velocityZ = _mm256_loadu_pd(velocity[2]);
    __m256d velocityMagnitudeSqrVector = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ));

    const double c_s_square = C_S * C_S;
    __m256d c_s_squareVector = _mm256_broadcast_sd(&c_s_square);

    const double one = 1.0;
    __m256d oneVector = _mm256_broadcast_sd(&one);

    for (int i = 0; i < PARAMQ; i++)
    {
        tmp = _mm_set1_epi32(LATTICEVELOCITIES[i][0]);
        latticeVelocityX = _mm256_cvtepi32_pd(tmp);
        tmp = _mm_set1_epi32(LATTICEVELOCITIES[i][1]);
        latticeVelocityY = _mm256_cvtepi32_pd(tmp);
        tmp = _mm_set1_epi32(LATTICEVELOCITIES[i][2]);
        latticeVelocityZ = _mm256_cvtepi32_pd(tmp);
        latticeWeightsVector = _mm256_broadcast_sd(&LATTICEWEIGHTS[i]);

        dotProductVector = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, latticeVelocityX), _mm256_mul_pd(velocityY, latticeVelocityY)), _mm256_mul_pd(velocityZ, latticeVelocityZ));
//        result = _mm256_add_pd(oneVector, _mm256_div_pd(dotProductVector, c_s_squareVector));
//        result = _mm256_add_pd(result, _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(c_s_squareVector, c_s_squareVector), _mm256_mul_pd(c_s_squareVector, c_s_squareVector))));
//        result = _mm256_sub_pd(result, _mm256_div_pd(velocityMagnitudeSqrVector, _mm256_add_pd(c_s_squareVector, c_s_squareVector)));
//        result = _mm256_mul_pd(result, _mm256_mul_pd(latticeWeightsVector, densityVector));
        result = _mm256_mul_pd(latticeWeightsVector, _mm256_mul_pd(densityVector, _mm256_sub_pd(_mm256_add_pd(oneVector, _mm256_add_pd(_mm256_div_pd(dotProductVector, c_s_squareVector), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(c_s_squareVector, c_s_squareVector),_mm256_mul_pd(c_s_squareVector, c_s_squareVector))))), _mm256_div_pd(velocityMagnitudeSqrVector, _mm256_add_pd(c_s_squareVector, c_s_squareVector)))));

        _mm256_storeu_pd( feq[i], result );
    }
}
void computeFeqAVX(const double * const density, const double * const velocity, double *feq)
{
    double velocityMagnitudeSqr, dotProduct;
    double const one = 1.0;
    __m256d dotProductVector, velocityX, velocityY, velocityZ, latWeights, result;
    __m128i temp;
    __m256d vl0, vl1, vl2;
    int i;
    double c_s_square = C_S * C_S;
    velocityX = _mm256_broadcast_sd(velocity);
    velocityY = _mm256_broadcast_sd(velocity + 1);
    velocityZ = _mm256_broadcast_sd(velocity + 2);

    velocityMagnitudeSqr = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    __m256d c_s_squareVector = _mm256_broadcast_sd(&c_s_square);
    __m256d velocityMagnitudeSqrVector = _mm256_broadcast_sd(&velocityMagnitudeSqr);
    __m256d densityVector = _mm256_broadcast_sd(density);
    __m256d oneVector = _mm256_broadcast_sd(&one);
    for (i = 0; i < PARAMQ - 3; i += 4)
    {
        latWeights = _mm256_loadu_pd(&LATTICEWEIGHTS[i]);
        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[0][i]);
        vl0 = _mm256_cvtepi32_pd(temp);
        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[1][i]);
        vl1 = _mm256_cvtepi32_pd(temp);
        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[2][i]);
        vl2 = _mm256_cvtepi32_pd(temp);
        dotProductVector = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, vl0), _mm256_mul_pd(velocityY, vl1)), _mm256_mul_pd(velocityZ, vl2));
        result = _mm256_mul_pd(latWeights, _mm256_mul_pd(densityVector, _mm256_sub_pd(_mm256_add_pd(oneVector, _mm256_add_pd(_mm256_div_pd(dotProductVector, c_s_squareVector), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(c_s_squareVector, c_s_squareVector),_mm256_mul_pd(c_s_squareVector, c_s_squareVector))))), _mm256_div_pd(velocityMagnitudeSqrVector, _mm256_add_pd(c_s_squareVector, c_s_squareVector)))));
        _mm256_storeu_pd(&feq[i], result);
    }
    for (; i < PARAMQ; i++)
    {
        dotProduct = velocity[0] * LATTICEVELOCITIES[i][0] + velocity[1] * LATTICEVELOCITIES[i][1] + velocity[2] * LATTICEVELOCITIES[i][2];
        feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + dotProduct / c_s_square + dotProduct * dotProduct / (2  * c_s_square * c_s_square) - velocityMagnitudeSqr / (2 * c_s_square));
    }
}
