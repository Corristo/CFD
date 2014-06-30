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
#include <papi.h>

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    double *velocities[3];
    double *densities;

    char problem[100];
    char pgmInput[1000];
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;

    int xlength[3], timesteps, timestepsPerPlotting;
    double tau, bddParams[7];

    long long counters[3];
	int PAPI_events[] = {
		PAPI_TOT_CYC,
		PAPI_L2_DCM,
		PAPI_L2_DCA };

	PAPI_library_init(PAPI_VER_CURRENT);

    double * exactCollideField;


    if(readParameters(xlength, &tau, bddParams, &timesteps, &timestepsPerPlotting, problem, pgmInput, argc, argv) == 0)
    {
        begin = clock();
        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
//        densities = (double *) malloc((size_t) sizeof(double) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
//        velocities[0] = (double *) malloc((size_t) sizeof(double) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
//        velocities[1] = (double *) malloc((size_t) sizeof(double) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
//        velocities[2] = (double *) malloc((size_t) sizeof(double) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
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

        PAPI_start_counters( PAPI_events, 3 );
        for(int t = 0; t < timesteps; t++)
        {
            double *swap = NULL;

//            doStreamingAVX(collideField, streamField, flagField, xlength);
//            doStreamingAVXv2(collideField, streamField, flagField, xlength, densities, velocities);
            doStreamingAndCollisionAVX(collideField, streamField, flagField, xlength, tau);
            swap = collideField;
            collideField = streamField;
            streamField = swap;

//            doCollisionAVX(collideField, flagField, &tau, xlength);
//            doCollisionAVXv2(collideField, flagField, tau, xlength, densities, velocities);

            treatBoundary(collideField, flagField, bddParams, xlength);

            if (t % timestepsPerPlotting == 0)
            {
//                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength);
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
//                for (int i = 0; i < PARAMQ; i++)
//                    for (int z = 1; z <= xlength[2]; z++)
//                        for (int y = 1; y <= xlength[1]; y++)
//                            for(int x = 1; x <= xlength[0]; x++)
//                                if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == FLUID)
//                                    if (fabs(collideField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2 ) + x + (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2) * i] - exactCollideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i]) > 1e-4)
//                                        error = 1;
//                if (error)
//                    printf("ERROR: Different collideField in timestep %d\n", t);
//                free(exactCollideField);
                /** end of debugging code */
            }
            PAPI_read_counters( counters, 3 );
            int pct = ((float) t / timesteps) * 100;

            printf("\b\b\b%02d%%", pct);
            fflush(stdout);


        }
        printf("\b\b\b\b100%%\n");
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

        printf("Running time: %.2f\n", time_spent);
        printf("MLUPS: %.3f\n", ((double) (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2) * timesteps) / (1000000.0 * time_spent));


        printf("%lld L2 cache misses (%.3lf%% misses) in %lld cycles\n",
		counters[1],
		(double)counters[1] / (double)counters[2] * 100,
		counters[0] );
        free(collideField);
        free(streamField);
        free(flagField);
//        free(velocities[0]);
//        free(velocities[1]);
//        free(velocities[2]);
//        free(densities);

    }
    return 0;
}

#endif

