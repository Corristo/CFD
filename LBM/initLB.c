#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *bddParams, int *timesteps, int *timestepsPerPlotting, char *problem, char *pgmInput, int argc, char *argv[])
{
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
    if (argc == 2)
    {
        read_string(argv[1], "problem", problem);
        read_int(argv[1], "xlength", xlength);
        read_int(argv[1], "ylength", xlength + 1);
        read_int(argv[1], "zlength", xlength + 2);
        READ_DOUBLE(argv[1], *tau);

        if(!strcmp(problem, "drivenCavity"))
        {
            read_double(argv[1], "velocityWallX", bddParams + 4);
            read_double(argv[1], "velocityWallY", bddParams + 5);
            read_double(argv[1], "velocityWallZ", bddParams + 6);
        }
        else
        {
            bddParams[4] = 0.0;
            bddParams[5] = 0.0;
            bddParams[6] = 0.0;
        }

        if(!strcmp(problem, "tiltedPlate"))
        {
            read_double(argv[1], "velocityInflowX", bddParams);
            read_double(argv[1], "velocityInflowY", bddParams + 1);
            read_double(argv[1], "velocityInflowZ", bddParams + 2);
            read_string(argv[1], "pgmInput", pgmInput);
        }
        else
        {
            bddParams[0] = 0.0;
            bddParams[1] = 0.0;
            bddParams[2] = 0.0;
            strcpy(pgmInput, "");
        }

        if (!strcmp(problem, "flowStep"))
        {
            read_double(argv[1], "velocityInflowX", bddParams);
            read_double(argv[1], "velocityInflowY", bddParams + 1);
            read_double(argv[1], "velocityInflowZ", bddParams + 2);
        }

        if(!strcmp(problem, "planeShearFlow"))
        {
            read_double(argv[1], "pressureIn", bddParams + 3);
        }
        else
        {
            bddParams[3] = 0.0;
        }
        READ_INT(argv[1], *timesteps);
        READ_INT(argv[1], *timestepsPerPlotting);
    }
    else

    {

        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return -1;
    }

    return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int *xlength, char *problem, char* pgmInput)
{
    int i, j, x, y, z;


    #pragma omp parallel for default(none), private(i, j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
    for (i = 0; i < PARAMQ; i++)
        for (j = 0; j < (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2); j++)
            collideField[j + i * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2)] = streamField[j+ i * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2)] = LATTICEWEIGHTS[i];
    #pragma omp parallel for default(none), private(i, j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
    for (i = 0; i < (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2); i++)
        flagField[i] = FLUID;

    if (!strcmp(problem, "drivenCavity"))
    {
        /* back boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[y * (xlength[0] + 2) + x] = NO_SLIP;

        /* bottom boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + x] = NO_SLIP;

        /* left boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)] = NO_SLIP;

        /* right boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1] = NO_SLIP;

        /* front boundary, i.e. z = xlength + 1 */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = NO_SLIP;

        /* top boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x] = MOVING_WALL;
    }
    if (!strcmp(problem, "tiltedPlate"))
    {
        int** pgmImage;
        pgmImage = read_pgm(pgmInput);
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem, pgmImage) schedule(static)
        for (z = 1; z <= xlength[2]; z++)
            for (y = 1; y <= xlength[1]; y++)
                for (x = 1; x <= xlength[0]; x++)
                    flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = !!pgmImage[x][y];
        free_imatrix(pgmImage, 0, xlength[0] + 2, 0, xlength[1] + 2);

        /* front boundary, i.e. z = xlength + 1 */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                if (flagField[xlength[2] * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) * y + x] == NO_SLIP)
                    flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = NO_SLIP;
                else
                    flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = FREE_SLIP;

        /* back boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                if (flagField[y * (xlength[0] + 2) + x] == NO_SLIP)
                    flagField[y * (xlength[0] + 2) + x] = NO_SLIP;
                else
                    flagField[y * (xlength[0] + 2) + x] = FREE_SLIP;


        /* left boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)] = INFLOW;

        /* right boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1] = OUTFLOW;


        /* top boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x] = NO_SLIP;

        /* bottom boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + x] = NO_SLIP;

    }
    if (!strcmp(problem, "flowStep"))
    {
        /* front boundary, i.e. z = xlength + 1 */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = NO_SLIP;

        /* back boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[y * (xlength[0] + 2) + x] = NO_SLIP;

        /* left boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)] = INFLOW;

        /* right boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1] = OUTFLOW;


        /* top boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x] = NO_SLIP;

        /* bottom boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + x] = NO_SLIP;


        /* step */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 1; y <= xlength[1]/2; y++) /* integer division on purpose, half of the channel is blocked by step */
                for (x = 1; x <= xlength[1]/2; x++)
                    flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = NO_SLIP;

    }
    if (!strcmp(problem, "planeShearFlow"))
    {

        /* front boundary, i.e. z = xlength + 1 */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] = FREE_SLIP;

        /* back boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[y * (xlength[0] + 2) + x] = FREE_SLIP;



        /* left boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)] = PRESSURE_IN;

        /* right boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1] = OUTFLOW;


        /* top boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x] = NO_SLIP;

        /* bottom boundary */
        #pragma omp parallel for default(none), private(i,j, x, y, z, pgmInput), shared(collideField, streamField, flagField, xlength, problem) schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
                flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + x] = NO_SLIP;

    }
    /** debugging code: checking the flagField */
//    int * exactFlagField;
//    exactFlagField = (int *) malloc( (size_t) sizeof( int ) * (xlength[0] + 2) *  (xlength[1] + 2) * (xlength[2] + 2));
//    FILE *fp2 = NULL;
//    unsigned int line2 = 0;
//    int error2 = 0;
//    char szFileName2[80];
//
//    sprintf( szFileName2, "Testdata/%s/flagField.dat", problem );
//    fp2 = fopen(szFileName2,"r");
//    if (fp2 != NULL)
//    {
//        for (line2 = 0; line2 < (xlength[0] + 2) *  (xlength[1] + 2) * (xlength[2] + 2); line2++)
//            fscanf(fp2,"%d",&exactFlagField[line2]);
//    }
//    fclose(fp2);
//    for (z = 1; z <= xlength[2]; z++)
//        for (y = 1; y <= xlength[1]; y++)
//            for(x = 1; x <= xlength[0]; x++)
//                for (i = 0; i < PARAMQ; i++)
//                    if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                        error2 = 1;
//    if (error2)
//        printf("ERROR: Different flagField at inner nodes.\n");
//
//    error2 = 0;
//    // Check global boundaries as well
//// Left boundary
//    x = 0;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (y = 0; y <= xlength[1] + 1; y++)
//            if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                error2 = 1;
//
//
//    // Right boundary
//    x = xlength[0] + 1;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (y = 0; y <= xlength[1] + 1; y++)
//            if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                error2 = 1;
//
//    // Bottom boundary
//    y = 0;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                error2 = 1;
//
//    // Top boundary
//    y = xlength[1] + 1;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                error2 = 1;
//
//    // back boundary
//    z = 0;
//    for (y = 0; y <= xlength[1] + 1; y++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                error2 = 1;
//
//    // front boundary
//    z = xlength[2] + 1;
//    for (y = 0; y <= xlength[1] + 1; y++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)])
//                error2 = 1;
//
//
//    if (error2)
//        printf("ERROR: Different flagField at global boundary.\n");
//
//    free(exactFlagField);

}

