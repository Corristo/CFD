#include "collision.h"
#include "LBDefinitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq, int tot_cells)
{
    int i;
    for (i = 0; i < PARAMQ; i++)
    {
        *(currentCell + tot_cells * i) = *(currentCell + tot_cells * i) * (1 - 1./ (*tau )) + 1./ (*tau)  * feq[i];
    }
}
void computePostCollisionDistributionsAVX (double *currentCell, const double * const tau, double feq[PARAMQ][4], int tot_cells)
{
    int i;
    __m256d currentDistributions, feqVector, tauVector, oneVector;
    const double one = 1.0;
    oneVector = _mm256_broadcast_sd(&one);
    tauVector = _mm256_broadcast_sd(tau);
    tauVector = _mm256_div_pd(oneVector, tauVector);
    for (i = 0; i < PARAMQ; i++)
    {
        currentDistributions = _mm256_loadu_pd(currentCell + tot_cells * i);
        feqVector = _mm256_loadu_pd(feq[i]);
        currentDistributions = _mm256_add_pd(_mm256_sub_pd(currentDistributions, _mm256_mul_pd(currentDistributions, tauVector)), _mm256_mul_pd(tauVector, feqVector));
        _mm256_storeu_pd(currentCell + tot_cells * i, currentDistributions);
    }

}
void doCollision(double *collideField, int *flagField,const double * const tau,int *xlength)
{
    double velocity[3], density,  feq[PARAMQ], *currentCell;
    int x, y, z, tot_cells =  (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
    for (z = 1; z <= xlength[2]; z++)
        for (y = 1; y <= xlength[1]; y++)
            for (x = 1; x <= xlength[0]; x++)
            {

                density = 0.0;
                currentCell = collideField + (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x);
                computeDensity(currentCell, &density, tot_cells);
                if (__builtin_expect(fabs(1.0 - density) > 0.1, 0))
                    fprintf(stderr, "WARNING: Density is %.3f in cell (%d, %d, %d)\n", density, x, y, z);

                computeVelocity(currentCell, &density, velocity, tot_cells);
                computeFeqAVX(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq, tot_cells);

            }
}

void doCollisionSSE(double *collideField, int *flagField,const double  tau, int *xlength, double * densities, double ** velocities)
{
    int x, y, z, j;
    double tmpFeq0, tmpFeq1, tmpFeq2, dotProduct;
    int numberOfCells = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);

    // i = 0 and i = 18
    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[1][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[j] = (1 - 1./tau) * collideField[j] + 1./tau * tmpFeq1;
                collideField[18 * numberOfCells + j] = (1 - 1./tau) * collideField[18 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 1 and i = 17
    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[0][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[numberOfCells + j] = (1 - 1./tau) * collideField[numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[17 * numberOfCells + j] = (1 - 1./tau) * collideField[17 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 2 and i = 16
    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[2][j];
                tmpFeq1 = 1./18 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[2 * numberOfCells + j] = (1 - 1./tau) * collideField[2 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[16 * numberOfCells + j] = (1 - 1./tau) * collideField[16 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 3 and i = 15
    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = velocities[0][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[3 * numberOfCells + j] = (1 - 1./tau) * collideField[3 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[15 * numberOfCells + j] = (1 - 1./tau) * collideField[15 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 4 and i = 14
    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = velocities[1][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[4 * numberOfCells + j] = (1 - 1./tau) * collideField[4 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[14 * numberOfCells + j] = (1 - 1./tau) * collideField[14 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 5 and i = 13
    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[0][j] - velocities[1][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[5 * numberOfCells + j] = (1 - 1./tau) * collideField[5 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[13 * numberOfCells + j] = (1 - 1./tau) * collideField[13 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 6 and i = 12
    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[1][j];
                tmpFeq1 = 1./18 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[6 * numberOfCells + j] = (1 - 1./tau) * collideField[6 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[12 * numberOfCells + j] = (1 - 1./tau) * collideField[12 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 7 and i = 11
    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = velocities[0][j] - velocities[1][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[7 * numberOfCells + j] = (1 - 1./tau) * collideField[7 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[11 * numberOfCells + j] = (1 - 1./tau) * collideField[11 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 8 and i = 10
    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[0][j];
                tmpFeq1 = 1./18 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[8 * numberOfCells + j] = (1 - 1./tau) * collideField[8 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[10 * numberOfCells + j] = (1 - 1./tau) * collideField[10 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 9
    // c_9 = (0, 0, 0), w_9 = 1/3
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                collideField[9 * numberOfCells + j] = (1 - 1. / tau) * collideField[9 * numberOfCells + j] + 1. / (3 * tau) * densities[j] * tmpFeq0;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }


}

void doCollisionAVXv2(double *collideField, int *flagField,const double  tau, int *xlength, double * densities, double ** velocities)
{
    int x, y, z, j;
    double tmpFeq0, tmpFeq1, tmpFeq2, dotProduct;
    int numberOfCells = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);

    __m256d densityVector, velocityX, velocityY, velocityZ;
    __m256d tmpFeq0v, tmpFeq1v, tmpFeq2v, dotProductVector;
    const double zero = 0.0;
    const double one = 1.0;
    const double c_s_square = C_S * C_S;
    __m256d cs2v = _mm256_broadcast_sd(&c_s_square);
    __m256d oneVector = _mm256_broadcast_sd(&one);
    __m256d zeroVector = _mm256_broadcast_sd(&zero);
    __m256d tauInvVector = _mm256_div_pd(oneVector, _mm256_broadcast_sd(&tau));
    __m256d weight36 = _mm256_broadcast_sd(&LATTICEWEIGHTS[7]);
    __m256d weight18 = _mm256_broadcast_sd(&LATTICEWEIGHTS[8]);
    __m256d weight3 = _mm256_broadcast_sd(&LATTICEWEIGHTS[9]);

    // i = 0 and i = 18
    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));

                _mm256_storeu_pd(collideField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 18 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[1][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[j] = (1 - 1./tau) * collideField[j] + 1./tau * tmpFeq1;
                collideField[18 * numberOfCells + j] = (1 - 1./tau) * collideField[18 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 1 and i = 17
    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 17 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[0][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[numberOfCells + j] = (1 - 1./tau) * collideField[numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[17 * numberOfCells + j] = (1 - 1./tau) * collideField[17 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 2 and i = 16
    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 2 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 16 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[2][j];
                tmpFeq1 = 1./18 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[2 * numberOfCells + j] = (1 - 1./tau) * collideField[2 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[16 * numberOfCells + j] = (1 - 1./tau) * collideField[16 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 3 and i = 15
    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 3 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 15 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = velocities[0][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[3 * numberOfCells + j] = (1 - 1./tau) * collideField[3 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[15 * numberOfCells + j] = (1 - 1./tau) * collideField[15 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 4 and i = 14
    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 4 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 14 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = velocities[1][j] - velocities[2][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[4 * numberOfCells + j] = (1 - 1./tau) * collideField[4 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[14 * numberOfCells + j] = (1 - 1./tau) * collideField[14 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 5 and i = 13
    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 5 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 13 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[0][j] - velocities[1][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[5 * numberOfCells + j] = (1 - 1./tau) * collideField[5 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[13 * numberOfCells + j] = (1 - 1./tau) * collideField[13 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 6 and i = 12
    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 6 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 12 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[1][j];
                tmpFeq1 = 1./18 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[6 * numberOfCells + j] = (1 - 1./tau) * collideField[6 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[12 * numberOfCells + j] = (1 - 1./tau) * collideField[12 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 7 and i = 11
    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 7 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 11 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = velocities[0][j] - velocities[1][j];
                tmpFeq1 = 1./36 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[7 * numberOfCells + j] = (1 - 1./tau) * collideField[7 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[11 * numberOfCells + j] = (1 - 1./tau) * collideField[11 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 8 and i = 10
    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(collideField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 8 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(collideField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 10 * numberOfCells + j)), _mm256_mul_pd(tauInvVector, tmpFeq2v)));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                dotProduct = - velocities[0][j];
                tmpFeq1 = 1./18 * densities[j] * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * densities[j] * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                collideField[8 * numberOfCells + j] = (1 - 1./tau) * collideField[8 * numberOfCells + j] + 1./tau * tmpFeq1;
                collideField[10 * numberOfCells + j] = (1 - 1./tau) * collideField[10 * numberOfCells + j] + 1./tau * tmpFeq2;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }

    // i = 9
    // c_9 = (0, 0, 0), w_9 = 1/3
    j = (xlength[0] + 2)*(xlength[1] + 3) + 1;
    for (z = 1; z <= xlength[2]; z++)
    {
        for (y = 1; y <= xlength[1]; y++)
        {
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                velocityX = _mm256_loadu_pd(&velocities[0][j]);
                velocityY = _mm256_loadu_pd(&velocities[1][j]);
                velocityZ = _mm256_loadu_pd(&velocities[2][j]);
                densityVector = _mm256_loadu_pd(&densities[j]);
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_add_pd(_mm256_mul_pd(velocityY, velocityY), _mm256_mul_pd(velocityZ, velocityZ))), _mm256_add_pd(cs2v, cs2v)));

                _mm256_storeu_pd(collideField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), _mm256_loadu_pd(collideField + 9 * numberOfCells + j)), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            for (; x <= xlength[0]; x++)
            {
                tmpFeq0 = 1 - (velocities[0][j] * velocities[0][j] + velocities[1][j] * velocities[1][j] + velocities[2][j] * velocities[2][j]) / (2 * C_S * C_S);
                collideField[9 * numberOfCells + j] = (1 - 1. / tau) * collideField[9 * numberOfCells + j] + 1. / (3 * tau) * densities[j] * tmpFeq0;
                j ++;
            }
            j += 2;
        }
        j += 2 * xlength[0] + 4;
    }


}

void doCollisionAVX(double *collideField, int *flagField,const double * const tau,int *xlength)
{
    double velocity[3][4] __attribute__((aligned(32))), density[4] __attribute__((aligned(32))), feq[PARAMQ][4] __attribute__((aligned(32))), *currentCell;
    double velocitySequential[3], densitySequential, feqSequential[PARAMQ] __attribute__((aligned(32)));
    int i, tot_cells =  (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
    for (i = 0; i < tot_cells - 3; i += 4)
    {
        currentCell = collideField + i;
        computeDensityAVX(currentCell, density, tot_cells);
        computeVelocityAVX(currentCell, density, velocity, tot_cells);
        computeFeqAVXv2(density, velocity, feq);
        computePostCollisionDistributionsAVX(currentCell, tau, feq, tot_cells);

    }
    for (; i < tot_cells; i++)
    {
        currentCell = collideField + i;
        computeDensity(currentCell, &densitySequential, tot_cells);
//        if (__builtin_expect(fabs(1.0 - density) > 0.1, 0))
//            fprintf(stderr, "WARNING: Density is %.3f in cell (%d, %d, %d)\n", density, x, y, z);

        computeVelocity(currentCell, &densitySequential, velocitySequential, tot_cells);
        computeFeqAVX(&densitySequential, velocitySequential, feqSequential);
        computePostCollisionDistributions(currentCell, tau, feqSequential, tot_cells);
    }
}

