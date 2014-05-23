#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdlib.h>

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
    int boundaryType;
    int coordinate[3];
    int iList[5];
    int fsList[5];

    /* top boundary
     * only need to set directions 0, 5, 6, 7, 14 */
    y = xlength[1] + 1; coordinate[1] = y;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (x = 0; x <= xlength[0] + 1; x++){

            boundaryType = flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z];
            coordinate[0] = x;  coordinate[2] = z;
            *iList = 0; *(iList + 1) = 5; *(iList + 2) = 7; *(iList + 3) = 14; *(iList + 4) = 6;

            if( boundaryType != FREE_SLIP ){
                compute_boundary(collideField, bddParams, flagField, boundaryType, xlength, iList, coordinate);
            } else  {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y - 1) * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y - 1;
                        neighbourCoordZ = z;

                        fsList[0] = 4; fsList[1] = 11; fsList[2] = 13; fsList[3] = 18; fsList[4] = 12;

                        for(i = 0; i < 5; i++){
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y - 1) * (xlength[0] + 2) + x] == NO_SLIP)
                    {

                        for(i = 0; i < 4; i++){
                            neighbourCoordX = x + LATTICEVELOCITIES[iList[i]][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[iList[i]][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[iList[i]][2];
                            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                                if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                                {
                                    currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                                    neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                                    collideField[currentCellIndex] = collideField[neighbourCellIndex];
                                }
                        }
                    }
            }
        }

    /* back boundary */
    // i = 14,15,16,17,18
    z = 0; coordinate[2] = 0;
    for (y = 0; y <= xlength[1] + 1; y++)
        for (x = 0; x <= xlength[0] + 1; x++){

            boundaryType = flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z];
            coordinate[0] = x; coordinate[1] = y; 
            *iList = 14; *(iList + 1) = 15; *(iList + 2) = 17; *(iList + 3) = 18; *(iList + 4) = 16;

            if( boundaryType != FREE_SLIP ){
                compute_boundary(collideField, bddParams, flagField, boundaryType, xlength, iList, coordinate);
            } else {
                    if (!flagField[(z + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y;
                        neighbourCoordZ = z + 1;

                        fsList[0] = 0; fsList[1] = 1; fsList[2] = 3; fsList[3] = 4; fsList[4] = 2;

                        for(i = 0; i < 5; i++){
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    }
                    else if (flagField[(z + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == NO_SLIP)
                    {

                        for(i = 0; i < 4; i++){
                            neighbourCoordX = x + LATTICEVELOCITIES[iList[i]][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[iList[i]][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[iList[i]][2];
                            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                                if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                                {
                                    currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                                    neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                                    collideField[currentCellIndex] = collideField[neighbourCellIndex];
                                }
                        }
                    }
            }
        }

    /* bottom boundary */
    // i = 4, 11, 12, 13, 18
    y = 0; coordinate[1] = 0;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (x = 0; x <= xlength[0] + 1; x++){

            boundaryType = flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z];
            coordinate[0] = x; coordinate[2] = z;
            *iList = 4; *(iList + 1) = 11; *(iList + 2) = 13; *(iList + 3) = 18; *(iList + 4) = 12;

            if( boundaryType != FREE_SLIP ){
                compute_boundary(collideField, bddParams, flagField, boundaryType, xlength, iList, coordinate);

            } else {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y + 1) * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y + 1;
                        neighbourCoordZ = z;

                        fsList[0] = 0; fsList[1] = 5; fsList[2] = 7; fsList[3] = 14; fsList[4] = 6;

                        for(i = 0; i < 5; i++){
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y + 1) * (xlength[0] + 2) + x] == NO_SLIP)
                    {
                        for(i = 0; i < 4; i++){
                            neighbourCoordX = x + LATTICEVELOCITIES[iList[i]][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[iList[i]][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[iList[i]][2];
                            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                                if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                                {
                                    currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                                    neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                                    collideField[currentCellIndex] = collideField[neighbourCellIndex];
                                }
                        }
                    }
            }
        }
            
    /* left boundary */
    // i = 3, 7 ,10, 13, 17
    x = 0; coordinate[0] = 0;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++){

            boundaryType = flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z];
            coordinate[1] = y; coordinate[2] = z;
            *iList = 3; *(iList + 1) = 7; *(iList + 2) = 13; *(iList + 3) = 17; *(iList + 4) = 10;

            if( boundaryType != FREE_SLIP ){
                compute_boundary(collideField, bddParams, flagField, boundaryType, xlength, iList, coordinate);
            } else{
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x + 1])
                    {

                        neighbourCoordX = x + 1;
                        neighbourCoordY = y;
                        neighbourCoordZ = z;

                        fsList[0] = 1; fsList[1] = 5; fsList[2] = 11; fsList[3] = 15; fsList[4] = 8;

                        for(i = 0; i < 5; i++){
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x + 1] == NO_SLIP)
                    {
                        for(i = 0; i < 4; i++){
                            neighbourCoordX = x + LATTICEVELOCITIES[iList[i]][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[iList[i]][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[iList[i]][2];
                            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                                if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                                {
                                    currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                                    neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                                    collideField[currentCellIndex] = collideField[neighbourCellIndex];
                                }
                        }
                    }
            }
        }
            
    // /* right boundary */
    // // i = 1,5,8,11,15
    x = xlength[0] + 1; coordinate[0] = x;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++){

            boundaryType = flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z];
            coordinate[1] = y; coordinate[2] = z;
            *iList = 1; *(iList + 1) = 5; *(iList + 2) = 11; *(iList + 3) = 15; *(iList + 4) = 8;

            if( boundaryType != FREE_SLIP ){
                compute_boundary(collideField, bddParams, flagField, boundaryType, xlength, iList, coordinate);
            } else {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x - 1])
                    {
                        neighbourCoordX = x - 1;
                        neighbourCoordY = y;
                        neighbourCoordZ = z;

                        fsList[0] = 3; fsList[1] = 7; fsList[2] = 13; fsList[3] = 17; fsList[4] = 10;

                        for(i = 0; i < 5; i++){
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x - 1] == NO_SLIP)
                    {
                        for(i = 0; i < 4; i++){
                            neighbourCoordX = x + LATTICEVELOCITIES[iList[i]][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[iList[i]][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[iList[i]][2];
                            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                                if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                                {
                                    currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                                    neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                                    collideField[currentCellIndex] = collideField[neighbourCellIndex];
                                }
                        }
                    }
            }
        }

    // /* front boundary, i.e. z = xlength + 1 */
    // i=0,1,2,3,4
    z = xlength[2] + 1;  coordinate[2] = z;
    for (y = 0; y <= xlength[1] + 1; y++)
        for (x = 0; x <= xlength[0] + 1; x++){

            boundaryType = flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z];
            coordinate[0] = x; coordinate[1] = y;
            *iList = 0; *(iList + 1) = 1; *(iList + 2) = 3; *(iList + 3) = 4; *(iList + 4) = 2;

            if( boundaryType != FREE_SLIP ){
                compute_boundary(collideField, bddParams, flagField, boundaryType, xlength, iList, coordinate);
            } else {
                    if (!flagField[(z - 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y ;
                        neighbourCoordZ = z - 1;

                        fsList[0] = 14; fsList[1] = 15; fsList[2] = 17; fsList[3] = 18; fsList[4] = 16;

                        for(i = 0; i < 5; i++){
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    }
                    else if (flagField[(z - 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == NO_SLIP)
                    {
                        for(i = 0; i < 4; i++){
                            neighbourCoordX = x + LATTICEVELOCITIES[iList[i]][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[iList[i]][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[iList[i]][2];
                            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                                if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                                {
                                    currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + iList[i];
                                    neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                                    collideField[currentCellIndex] = collideField[neighbourCellIndex];
                                }
                        }
                    }
            }
        }

    /* inner boundary cells
     * assumes inner boundary cells can only be NO_SLIP */
    for (z = 1; z <= xlength[2]; z++)
        for (y = 1; y <= xlength[1]; y++)
            for(x = 1; x <= xlength[0]; x++)
                if (!!flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
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
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
}

/** boundary helper function */
void compute_boundary(double *collideField, const double * const bddParams, int *flagField,
                         int boundaryType, int *xlength, int *iList, int *const coordinate){

    int i;
    int neighbourCoordX, neighbourCoordY, neighbourCoordZ;
    int currentCellIndex, neighbourCellIndex;
    double density;
    const double * const wallVelocity = bddParams + 4;
    double referenceDensity = 1.0;
    double feq[PARAMQ];
    double velocity[3];
    const double pressureIn = bddParams[3];


    switch (boundaryType)
    {
        case NO_SLIP:
        {
            for(i = 0; i < 5; i++ ){
                    neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
                    neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
                    neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1 
                        && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
            }
            break;
        }

        case MOVING_WALL:
        {
            for(i = 0; i < 5; i++ ){
                    neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
                    neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
                    neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1 
                        && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[iList[i]] * (LATTICEVELOCITIES[iList[i]][0] * wallVelocity[0] 
                                                                + LATTICEVELOCITIES[iList[i]][1] * wallVelocity[1] + LATTICEVELOCITIES[iList[i]][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
            }
            break;
        }

        case INFLOW:
        {
            for(i = 0; i < 5; i++ ){
                neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
                neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
                neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1 
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                        computeFeq(&referenceDensity, bddParams, feq);
                        collideField[currentCellIndex] = feq[iList[i]];
                    }
            }

            break;
        }

        case OUTFLOW:
        {
            for(i = 0; i < 5; i++ ){                  
                neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
                neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
                neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1 
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                        density = 0.0;
                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                        computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                        computeFeq(&referenceDensity, velocity, feq);
                        collideField[currentCellIndex] = feq[PARAMQ - iList[i] - 1] + feq[iList[i]] - collideField[neighbourCellIndex];
                    }
            }
            break;
        }

        case PRESSURE_IN:
        {
            for(i = 0; i < 5; i++ ){                  
                neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
                neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
                neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1 
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - iList[i] - 1] + feq[iList[i]] - collideField[neighbourCellIndex];
                        }
            }            
            break;
        }
    }

}

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



