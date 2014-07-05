#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <stdlib.h>
#include "helper.h"


void computeDensity(const double *const currentCell, double *density, const int tot_cells)
{
    *density = 0.0;
    int i;
    for (i = 0; i < PARAMQ; i++)
        *density += *(currentCell + i * tot_cells);
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

// TODO: Rework computeFeq
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
