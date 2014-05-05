#include "collision.h"
#include "LBDefinitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq)
{
    for (int i = 0; i < PARAMQ; i++)
    {
        currentCell[i] = currentCell[i] - 1. / (*tau) * (currentCell[i] - feq[i]);
     /*   if (__builtin_expect(currentCell[i] < 0.0, 0))
        {
                fprintf(stderr, "ERROR: Probability density is negative after post collision update: %.3f\n", currentCell[i]);
        } */
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength)
{
    double velocity[3], density,  feq[PARAMQ], *currentCell;
    for (int z = 1; z <= xlength; z++)
        for (int y = 1; y <= xlength ; y++)
            for (int x = 1; x <= xlength; x++)
            {
                currentCell = collideField + PARAMQ * (z * (xlength + 2) * (xlength + 2) + y * (xlength + 2) + x);
                density = 0.0;
             //   computeDensity(currentCell, &density);
                computeDensitySSE(currentCell, &density);
          /*      if (fabs(density - density2) > 0.0000001)
                {
                    fprintf(stderr, "New density fn does not work correctly\n density old: %.3f\n density new: %.3f\n", density, density);
                    abort();
                } */
           /*     if (__builtin_expect(fabs(1.0 - density) > 0.05, 0))
                    fprintf(stderr, "WARNING: Density is %.3f in cell (%d, %d, %d)\n", density, x, y, z); */
             //   computeVelocity(currentCell, &density, velocity2);
                computeVelocitySSE(currentCell, &density, velocity);
             /*  if (fabs(velocity[0] - velocity2[0]) > 0.000001)
                {
                    fprintf(stderr, "New velocity fn does not work correctly\n Velocity old: %.3f\n Velocity new: %.3f\n", velocity[0], velocity2[0]);
                    abort();
                } */
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);

            }
}

