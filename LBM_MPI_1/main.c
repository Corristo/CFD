#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include <time.h>
#include <unistd.h>
#include <parallel.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    char problem[100];
    char pgmInput[1000];
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;
    int rank ;
    int number_of_ranks ;
    double *sendBuffer[6];
    double *readBuffer[6];
    int iProc,jProc,kProc ;
    int xlength_global[3],timesteps, timestepsPerPlotting;
    double tau, bddParams[7];
    MPI_Status status;



    initializeMPI(&rank,&number_of_ranks,argc,argv);
/*  
    printf("rank = %d \n" , rank);
    printf("argv[0] = %s \n" , argv[0]); 
    printf("argv[1] = %s \n" , argv[1]);
    printf("argv[2] = %s \n" , argv[2]);
*/
    int xlength[3];
    if((readParameters(xlength_global, &tau, bddParams, &timesteps, &timestepsPerPlotting, problem, pgmInput, argc, argv,&iProc,&jProc,&kProc) == 0)&&(number_of_ranks==iProc*jProc*kProc))
    {
  
        begin = clock();
        xlength[0] = (xlength_global[0]/iProc) ;xlength[1] = (xlength_global[1]/jProc) ;
        xlength[2] = (xlength_global[2]/kProc) ;

        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        
        // initialise pointers here with correct size for the domain decomposition!
        initialiseFields(collideField, streamField, flagField, xlength, problem, pgmInput,rank,number_of_ranks,iProc,jProc,kProc);

        // allocate the buffers
        initialiseBuffers(sendBuffer,readBuffer,xlength);
        
        writeVtkOutput(streamField, flagField, "./Paraview/output", 0, xlength,rank);
        //printf("Progress:     ");
        int il,ir,jb,jt,kf,kb;
        decideneighbours(&il,&ir,&jb,&jt,&kf,&kb,iProc,jProc,kProc,rank,xlength);
        
        for(int t = 0; t < 1; t++)
        {
            double *swap = NULL;
            
//          Extraction
            extract( sendBuffer , collideField , xlength );
//          Swap
                        
            swap_send_read( sendBuffer , readBuffer , xlength ,il,ir,jb,jt,kf,kb,rank,&status);            

//          Injection
            inject(readBuffer,collideField,xlength);
            doStreaming(collideField, streamField, flagField, xlength);
            swap = collideField;
            collideField = streamField;
            streamField = swap;

            doCollision(collideField, flagField, &tau, xlength,rank);

            treatBoundary(collideField, flagField, bddParams, xlength);

            if (t % timestepsPerPlotting == 0)
                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength,rank);
            int pct = ((float) t / timesteps) * 100;

            printf("\b\b\b%02d%%", pct);
            fflush(stdout);
        }
/*        
        printf("\b\b\b\b100%%\n");
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

        printf("Running time: %.2f\n", time_spent);
        printf("MLUPS: %.3f\n", ((double) (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2) * timesteps) / (1000000.0 * time_spent));
*/
        free(collideField);
        free(streamField);
        free(flagField);

    }
    finalizeMPI();
    return 0;
}

#endif

