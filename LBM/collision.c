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
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
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

void doCollisionAVX(double *collideField, int *flagField,const double * const tau,int *xlength)
{
    double velocity[3][4], density[4], feq[PARAMQ][4], *currentCell;
    double velocitySequential[3], densitySequential, feqSequential[PARAMQ];
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

