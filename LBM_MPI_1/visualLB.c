#include <stdio.h>
#include "helper.h"
#include "visualLB.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"

void write_vtkHeader( FILE *fp, int *xlength)
{
    if( fp == NULL )
    {
        char szBuff[80];
        sprintf( szBuff, "Null pointer in write_vtkHeader" );
        ERROR( szBuff );
        return;
    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength[0], xlength[1], xlength[2]);
    fprintf(fp,"POINTS %i float\n", (xlength[0]) * (xlength[1]) * (xlength[2]));
    fprintf(fp,"\n");
}

void write_vtkPointCoordinates( FILE *fp, int *xlength)
{
    double originX = 0.0;
    double originY = 0.0;
    double originZ = 0.0;

    int x = 0;
    int y = 0;
    int z = 0;

    for(z = 1; z <= xlength[2]; z++)
        for(y = 1; y <= xlength[1]; y++)
            for ( x = 1; x <= xlength[0]; x++)
                fprintf(fp, "%f %f %f\n", originX+(double) x , originY + (double) y,  originZ + (double) z);


}


void writeVtkOutput(const double * const collideField, const int * const flagField, const char * filename, unsigned int t, int *xlength,int rank)
{
    int x, y, z, currentCellIndex;
    double cellVelocity[3], cellDensity;

    char szFileName[80];
    FILE *fp=NULL;
    sprintf( szFileName, "%s.%i.%i.vtk", filename, rank , t );
    fp = fopen( szFileName, "w");
    if( fp == NULL )
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to open %s", szFileName );
        ERROR( szBuff );
        return;
    }

    write_vtkHeader( fp, xlength);
    write_vtkPointCoordinates(fp, xlength);

    fprintf(fp,"POINT_DATA %i \n", (xlength[0]) * (xlength[1]) * (xlength[2]));

    fprintf(fp,"\n");
    fprintf(fp, "VECTORS velocity float\n");
    int globalcellindex;
    for(z = 1; z <= xlength[2]; z++)
        for (y = 1; y <= xlength[1]; y++)
            for (x = 1; x <= xlength[0]; x++)
            {
                currentCellIndex = PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x);
                computeDensitySSE(collideField + currentCellIndex, &cellDensity);
                computeVelocitySSE(collideField + currentCellIndex, &cellDensity, cellVelocity);
                
                fprintf(fp, "%f %f %f %d\n", cellVelocity[0], cellVelocity[1], cellVelocity[2],globalcellindex);
            }


    fprintf(fp,"\n");

    fprintf(fp, "SCALARS density float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for(z = 1; z <= xlength[2]; z++)
        for(y = 1; y <= xlength[1]; y++)
            for (x = 1; x <= xlength[0]; x++)
            {
                currentCellIndex = PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x);
                computeDensitySSE(collideField + currentCellIndex, &cellDensity);
                fprintf(fp, "%f %d\n", cellDensity,globalcellindex);
            }



    if( fclose(fp) )
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to close %s", szFileName );
        ERROR( szBuff );
    }
}

