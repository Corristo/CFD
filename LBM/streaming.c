#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int *xlength)
{
// collision optimized version
/*
    int streamCellIndex, neighbourCellIndex, x ,y ,z, i;
    for (z = 1; z <= xlength[2] ; z++)
        for (y = 1; y <= xlength[1] ; y++)
            for (x = 1; x <= xlength[0] ; x++)
                if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    for (i = 0; i < PARAMQ; i++)
                    {
                        streamCellIndex = PARAMQ * (x + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2)) + i;
                        neighbourCellIndex = streamCellIndex - PARAMQ * (LATTICEVELOCITIES[i][0] + LATTICEVELOCITIES[i][1] * (xlength[0] + 2) + LATTICEVELOCITIES[i][2] * (xlength[0] + 2) * (xlength[1] + 2));
                        streamField[streamCellIndex] = collideField[neighbourCellIndex];
                    }
*/

// streaming/propogation optimized version
    int streamCellIndex , neighbourCellIndex, x ,y ,z, i;
    //number_cells = (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) ;
    for (i = 0; i < PARAMQ; i++)
       for (z = 1; z <= xlength[2] ; z++)
          for (y = 1; y <= xlength[1] ; y++)
             for (x = 1; x <= xlength[0] ; x++)
                if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                {
                    streamCellIndex = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2)*i + z*(xlength[0]+2)*(xlength[1]+2) + y*(xlength[0]+2) + x;
                    neighbourCellIndex = streamCellIndex - (LATTICEVELOCITIES[i][0] + LATTICEVELOCITIES[i][1] * (xlength[0] + 2) + LATTICEVELOCITIES[i][2] * (xlength[0] + 2) * (xlength[1] + 2));
//                    neighbourCellIndex = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2)*i + (z-LTTICEVELOCITIES[i][2])*(xlength[0]+2)*(xlength[1]+2) + (y-LATTICEVELOCITIES[i][1])*(xlength[0]+2) + (x-LATTICEVELOCITIES[i][0]);
                    streamField[streamCellIndex] = collideField[neighbourCellIndex];
                }
}

void doStreamingSSE(double *collideField, double *streamField,int *flagField,int *xlength)
{
// only for streaming/propogation optimized version
    int streamCellIndex , neighbourCellIndex, x ,y ,z, i;
    //number_cells = (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) ;
    for (i = 0; i < PARAMQ; i++)
       for (z = 1; z <= xlength[2] ; z++)
          for (y = 1; y <= xlength[1] ; y++)
             for (x = 1; x <= xlength[0] ; x++)
                if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                {
                    streamCellIndex = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2)*i + z*(xlength[0]+2)*(xlength[1]+2) + y*(xlength[0]+2) + x;
                    neighbourCellIndex = streamCellIndex - (LATTICEVELOCITIES[i][0] + LATTICEVELOCITIES[i][1] * (xlength[0] + 2) + LATTICEVELOCITIES[i][2] * (xlength[0] + 2) * (xlength[1] + 2));
//                    neighbourCellIndex = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2)*i + (z-LTTICEVELOCITIES[i][2])*(xlength[0]+2)*(xlength[1]+2) + (y-LATTICEVELOCITIES[i][1])*(xlength[0]+2) + (x-LATTICEVELOCITIES[i][0]);
                    streamField[streamCellIndex] = collideField[neighbourCellIndex];
                }
}

void doStreamingAVX(double *collideField, double *streamField,int *flagField,int *xlength)
{
// only for streaming/propogation optimized version
    int streamCellIndex , neighbourCellIndex, x ,y ,z, i;
    //number_cells = (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) ;
    for (i = 0; i < PARAMQ; i++)
       for (z = 1; z <= xlength[2] ; z++)
          for (y = 1; y <= xlength[1] ; y++)
             for (x = 1; x <= xlength[0] ; x++)
                if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                {
                    streamCellIndex = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2)*i + z*(xlength[0]+2)*(xlength[1]+2) + y*(xlength[0]+2) + x;
                    neighbourCellIndex = streamCellIndex - (LATTICEVELOCITIES[i][0] + LATTICEVELOCITIES[i][1] * (xlength[0] + 2) + LATTICEVELOCITIES[i][2] * (xlength[0] + 2) * (xlength[1] + 2));
//                    neighbourCellIndex = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2)*i + (z-LTTICEVELOCITIES[i][2])*(xlength[0]+2)*(xlength[1]+2) + (y-LATTICEVELOCITIES[i][1])*(xlength[0]+2) + (x-LATTICEVELOCITIES[i][0]);
                    streamField[streamCellIndex] = collideField[neighbourCellIndex];
                }
}
