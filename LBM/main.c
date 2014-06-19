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

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    char problem[100];
    char pgmInput[1000];
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;

    int xlength[3], timesteps, timestepsPerPlotting;
    double tau, bddParams[7];

//    double * exactCollideField;


    if(readParameters(xlength, &tau, bddParams, &timesteps, &timestepsPerPlotting, problem, pgmInput, argc, argv) == 0)
    {
        begin = clock();
        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        initialiseFields(collideField, streamField, flagField, xlength, problem, pgmInput);

        /** debugging code */
//        /* output the flagField */
//        char szFileName2[80];
//        FILE *fp2 = NULL;
//        sprintf( szFileName2, "Testdata/%s/flagField.dat", problem);
//        fp2 = fopen(szFileName2,"w");
//        for (int i = 0; i < (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2); i++)
//                    fprintf(fp2, "%d\n", flagField[i]);
        /** debugging code end */

        printf("Progress:     ");
        for(int t = 0; t < timesteps; t++)
        {
            double *swap = NULL;

            doStreaming(collideField, streamField, flagField, xlength);
            swap = collideField;
            collideField = streamField;
            streamField = swap;

            doCollision(collideField, flagField, &tau, xlength);

            treatBoundary(collideField, flagField, bddParams, xlength);

            if (t % timestepsPerPlotting == 0)
            {
                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength);
                 /** debugging code */
//                 /* create reference files */
//                FILE *fp = NULL;
//                char szFileName[80];
//                sprintf( szFileName, "Testdata/%s/%i.dat", problem, t / timestepsPerPlotting );
//                fp = fopen(szFileName,"w");
//                for (int i = 0; i < PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2); i++)
//                    fprintf(fp, "%0.7f\n", collideField[i]);


                /* check correctness */
//                exactCollideField = (double *) malloc ( ( size_t ) sizeof(double) * PARAMQ * (xlength[0] + 2) *  (xlength[1] + 2) * (xlength[2] + 2));
//                FILE *fp = NULL;
//                unsigned int line = 0;
//                int noOfReadEntries;
//                int error = 0;
//                char szFileName[80];
//                sprintf( szFileName, "Testdata/%s/%i.dat", problem, t / timestepsPerPlotting );
//                fp = fopen(szFileName,"r");
//                if (fp != NULL)
//                {
//                    for (line = 0; line < PARAMQ * (xlength[0] + 2) *  (xlength[1] + 2) * (xlength[2] + 2); line++)
//                    {
//                        noOfReadEntries = fscanf(fp,"%lf",&exactCollideField[line]);
//                        if (noOfReadEntries != 1)
//                            continue;
//                    }
//                }
//                fclose(fp);
//                for (int z = 1; z <= xlength[2]; z++)
//                    for (int y = 1; y <= xlength[1]; y++)
//                        for(int x = 1; x <= xlength[0]; x++)
//                            for (int i = 0; i < PARAMQ; i++)
//                                if (fabs(collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2 + x) + i)] - exactCollideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2 + x) + i)]) > 1e-4)
//                                    error = 1;
//                if (error)
//                    printf("ERROR: Process %d has a different collideField in timestep %d\n", 0, t);
//                free(exactCollideField);
                /** end of debugging code */
            }

            int pct = ((float) t / timesteps) * 100;

            printf("\b\b\b%02d%%", pct);
            fflush(stdout);


        }
        printf("\b\b\b\b100%%\n");
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

        printf("Running time: %.2f\n", time_spent);
        printf("MLUPS: %.3f\n", ((double) (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2) * timesteps) / (1000000.0 * time_spent));

        free(collideField);
        free(streamField);
        free(flagField);

    }
    return 0;
}

#endif

