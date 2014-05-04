#ifndef _MAIN_C_
#define _MAIN_C_

#include <stdio.h>
#include <stdlib.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){
  
    /** Initialization */
	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	int xlength;
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
	int t;
	double *swap;

    int msg = readParameters( &xlength, &tau, velocityWall, &timesteps, 
    	                  &timestepsPerPlotting, argc, argv );
    if (msg != 0){
    	fprintf(stderr, "Error in processing function readParameters.\n");
    	return -1;
    }

    /** allocate mems as D3Q19 */
    int dim = xlength + 2;
	collideField = (double*)calloc( (size_t)( 19 * pow(dim, 3) ), sizeof(double) );
	streamField = (double*)calloc( (size_t)( 19 * pow(dim, 3) ), sizeof(double) );
	flagField = (int*)calloc( (size_t)(pow(dim, 3)), sizeof(int) );
	
    initialiseFields(collideField, streamField, flagField, xlength);

    for(t = 0; t < timesteps; t++){
	    swap = NULL;
	    doStreaming(collideField, streamField, flagField, xlength);
	    swap = collideField;
	    collideField = streamField;
	    streamField = swap;

	    doCollision(collideField, flagField, &tau, xlength);
	    treatBoundary(collideField, flagField, velocityWall, xlength);

	    if(t % timestepsPerPlotting == 0){
	    	writeVtkOutput(collideField, flagField, argv[1], t, xlength);
	    }
	}

    		
    return 0;
}

#endif

