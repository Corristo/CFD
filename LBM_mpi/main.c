#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "communication.h"
#include <time.h>
#include <unistd.h>
#include <mpi.h>
#include "parallel.h"


int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    char problem[100];
    char pgmInput[1000];
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;

    int xlength[3], local_xlength[3], timesteps, timestepsPerPlotting;
    double tau, bddParams[7];

    int rank;
    int number_of_ranks;
    int iProc, jProc, kProc;

    int neighbours[6];  // [0: left, 1: right, 2: top, 3: bottom, 4: front, 5: back]

    // send and read buffers for all possible directions :
    // [0: left, 1: right, 2: top, 3: bottom, 4: front, 5: back]
    double *sendBuffer[6];
    double *readBuffer[6];

    double * exactCollideField; // for debugging only



    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );      /* asking for the number of processes  */
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );                 /* asking for the local process id   */

    if(readParameters( xlength, &tau, bddParams, &iProc, &jProc, &kProc, &timesteps, &timestepsPerPlotting, problem, pgmInput, argc, argv ) == 0)
    {
        if (number_of_ranks != iProc * jProc * kProc)
        {
            if (rank == 0)
                printf("ERROR: number of processes started does not match the number specified in the input file!\n");
            MPI_Barrier( MPI_COMM_WORLD);
            MPI_Finalize();
            return -1;
        }
        begin = clock();
        time_t start = time(NULL);


        int iCoord, jCoord, kCoord; // position of the domain in the decomposition; needed to compute neighbouring processes and for vtk output
        computePosition( iProc, jProc, kProc, &iCoord, &jCoord, &kCoord );
        computeLocalDomainSize( xlength, local_xlength, iProc, jProc, kProc );

        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (local_xlength[0] + 2)*(local_xlength[1] + 2)*(local_xlength[2] + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (local_xlength[0] + 2)*(local_xlength[1] + 2)*(local_xlength[2] + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (local_xlength[0] + 2)*(local_xlength[1] + 2)*(local_xlength[2] + 2));

        computeNeighbours( iProc, jProc, kProc, neighbours );
        initialiseBuffers( local_xlength, sendBuffer, readBuffer, neighbours );
        initialiseFields( collideField, streamField, flagField, xlength, local_xlength, problem, pgmInput, rank, iProc, jProc, kProc );

        if (!rank)
            printf("Progress:     ");
        for(int t = 0; t < timesteps; t++)
        {
            double *swap = NULL;
            communicateBoundaryValues(local_xlength, sendBuffer, readBuffer, flagField, collideField, neighbours, 0); // communicate along x - axis (right to left - left to right)
            communicateBoundaryValues(local_xlength, sendBuffer, readBuffer, flagField, collideField, neighbours, 2); // communicate along z - axis (back to front - front to back)
            communicateBoundaryValues(local_xlength, sendBuffer, readBuffer, flagField, collideField, neighbours, 1); // communicate along y - axis (bottom to top - top to bottom)
            doStreaming(collideField, streamField, flagField, local_xlength);
            swap = collideField;
            collideField = streamField;
            streamField = swap;

            doCollision(collideField, flagField, &tau, local_xlength);
            treatBoundary(collideField, flagField, bddParams, local_xlength);

            if (t % timestepsPerPlotting == 0)
                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength, local_xlength, rank, iCoord, jCoord, kCoord, iProc, jProc, kProc);

            if (!rank)
            {
                int pct = ((float) t / timesteps) * 100;
                printf("\b\b\b%02d%%", pct);
                fflush(stdout);
            }

            /** debugging code: check collideField */
            /* check correctness of collideField with reference data */
            if (t % timestepsPerPlotting == 0)
            {

                exactCollideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));

                int x, y, z, i;
                FILE *fp = NULL;
                unsigned int line = 0;
                int error = 0;
                char szFileName[1200];
                sprintf( szFileName, "Testdata/%s/%i.dat", problem, t / timestepsPerPlotting );
                fp = fopen(szFileName,"r");
                if (fp != NULL)
                {
                    for (line = 0; line < PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2); line++)
                        fscanf(fp,"%lf",&exactCollideField[line]);

                    for (z = 1; z <= local_xlength[2]; z++)
                        for (y = 1; y <= local_xlength[1]; y++)
                            for(x = 1; x <= local_xlength[0]; x++)
                                for (i = 0; i < PARAMQ; i++)
                                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FLUID)
                                       if (fabs(collideField[PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + i] - exactCollideField[PARAMQ * ((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1]/jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc) )) + i]) > 1e-5)
                                            error = 1;
                    if (error)
                        printf("ERROR: Process %d has a different collideField in timestep %d\n",rank, t);

                    fclose(fp);

                }
                else
                    printf("ERROR: Process %d cannot read file %s\n", rank, szFileName);

                free(exactCollideField);

            }
            /** end of debugging code */


        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (!rank)
        {
            printf("\b\b\b\b100%%\n");
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

            printf("Running time (CPU time): %.2fs\n", time_spent);
            printf("Running time (Wall clock): %2.fs\n", (double)(time(NULL) - start) );
            printf("MLUPS: %.3f\n", ((double) (xlength[0]) * (xlength[1]) * (xlength[2]) * timesteps) / (1000000.0 * time_spent));
        }

        free(collideField);
        free(streamField);
        free(flagField);
        free(sendBuffer[0]);
        free(sendBuffer[1]);
        free(sendBuffer[2]);
        free(sendBuffer[3]);
        free(sendBuffer[4]);
        free(sendBuffer[5]);
        free(readBuffer[0]);
        free(readBuffer[1]);
        free(readBuffer[2]);
        free(readBuffer[3]);
        free(readBuffer[4]);
        free(readBuffer[5]);
    }

    MPI_Finalize();
    return 0;
}

#endif
