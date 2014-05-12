#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"


/**
 *  The array bddParams has the following structure
 *  bddParams[0]   Inflow velocity in x direction
 *  bddParams[1]   Inflow velocity in y direction
 *  bddParams[2]   Inflow velocity in z direction
 *  bddParams[3]   Pressure surplus for PRESSURE_IN cells
 *  bddParams[4]   Moving wall velocity in x direction
 *  bddParams[5]   Moving wall velocity in y direction
 *  bddParams[6]   Moving wall velocity in z direction
 */
void treatBoundary(double *collideField, int* flagField, const double * const bddParams, int *xlength)
{
    int i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex;
    const double * const wallVelocity = bddParams + 4;
    const double pressureIn = bddParams[3];
    double density = 0.0;
    double velocity[3];
    double feq[PARAMQ];
    double referenceDensity = 1.0;

    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++)
            for(x = 0; x <= xlength[0] + 1; x++)
                if (__builtin_expect(!!flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z], 0))
                {
                    for (i = 0; i < PARAMQ; i++)
                    {
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (!flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ])
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
                                {
                                    case NO_SLIP: collideField[currentCellIndex] = collideField[neighbourCellIndex]; break;
                                    case MOVING_WALL:
                                    {
                                        double density = 0.0;
                                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                                        collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                          + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                                        break;
                                    }
                                    case FREE_SLIP:
                                    {

                                    }
                                    case INFLOW:
                                    {
                                        computeFeq(&referenceDensity, bddParams, feq);
                                        collideField[currentCellIndex] = feq[i];
                                        break;
                                    }
                                    case OUTFLOW:
                                    {
                                        density = 0.0;
                                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                                        computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                                        computeFeq(&referenceDensity, velocity, feq);
                                        collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                                        break;

                                    }
                                    case PRESSURE_IN:
                                    {
                                        density = 0.0;
                                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                                        computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                                        density = referenceDensity + pressureIn;
                                        computeFeq(&density, velocity, feq);
                                        collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                                        break;
                                    }

                                }
                            }
                        }
                    }




//
//
//
//
//    /* top boundary */
//    y = xlength[1] + 1;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        density = 0.0;
//                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
//                        collideField[currentCellIndex] =  collideField[neighbourCellIndex]
//                                                          + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
//                    }
//            }
//    /* back boundary */
//    z = 0;
//    for (y = 0; y <= xlength[1] + 1; y++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
//                    }
//            }
//
//
//    /* bottom boundary */
//    y = 0;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
//                    }
//            }
//
//    /* left boundary */
//    x = 0;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (y = 0; y <= xlength[1] + 1; y++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
//                    }
//            }
//
//    /* right boundary */
//    x = xlength[0] + 1;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (y = 0; y <= xlength[1] + 1; y++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
//                    }
//            }
//
//    /* front boundary, i.e. z = xlength + 1 */
//    z = xlength[2] + 1;
//    for (y = 0; y <= xlength[1] + 1; y++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
//                    }
//            }
//


}

