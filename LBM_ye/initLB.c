#include "initLB.h"
#include "helper.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
  
    /** Argument checking */
    if ( argc != 2 ) {
        fprintf(stderr, "Usage: ./lbsim <input_file>\n");
        return -1;
    }
  	
  	/** read file and get values */
    read_int( argv[1], "xlength", xlength);
    read_double( argv[1], "tau", tau);
    read_double( argv[1], "velocityWallx", &velocityWall[0]);
    read_double( argv[1], "velocityWally", &velocityWall[1]);
    read_double( argv[1], "velocityWallz", &velocityWall[2]);
    read_int( argv[1], "timesteps", timesteps);
    read_int( argv[1], "timestepsPerPlotting", timestepsPerPlotting);

    return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

	int Q = 19; // number of velocities in our case
	int dim = xlength + 2;
	int x, y, z, i, index;
	int dim_pow2 = dim * dim;
	int xlen_pow2 = (dim - 1) * dim;
	int xlen_pow3 = xlen_pow2 * dim;

	for( z = 0; z < dim; z++ )
		for( y = 0; y < dim; y++ )
			for( x = 0; x < dim; x++ ) 
				for( i = 0; i < Q; i++ ){
					index = Q * (z * dim_pow2 + y * dim + x) + i;
					collideField[index] = LATTICEWEIGHTS[i];
					streamField[index] = LATTICEWEIGHTS[i];
				}

	// flagField
	// bottom, z = 0
 	for( x = 0; x < dim; x++ )
 		for( y = 0; y < dim; y++ )
 			flagField[ y * dim + x ] = 1;

 	// x = 0
    for( y = 1; y < dim - 1; y++ )
 		for( z = 1; z < dim - 1; z++ )
 			flagField[ z * dim_pow2 + y * dim ] = 1;

 	// y = 0
 	for( x = 0; x < dim; x++ )
 		for( z = 1; z < dim - 1; z++ )
 			flagField[ z * dim_pow2 + x ] = 1;

 	// x = dim - 1
 	for( y = 1; y < dim - 1; y++ )
 		for( z = 1; z < dim - 1; z++ )
 			flagField[ z * dim_pow2 + y * dim + dim - 1 ] = 1;

 	// y = dim - 1
 	for( x = 0; x < dim; x++ )
 		for( z = 1; z < dim - 1; z++ )
 			flagField[ z * dim_pow2 + xlen_pow2 + x ] = 1;

 	// top z = dim - 1
 	for( x = 0; x < dim; x++ )
 		for( y = 0; y < dim; y++ )
 			flagField[ xlen_pow3 + y * dim + x ] = 2;
    
}

