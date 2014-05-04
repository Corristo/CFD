#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
  int x, y, z, i;
  int dim = xlength + 2, dim_pow2 = dim * dim;

  for( z = 1; z < xlength + 1; z++ )
  	for( y = 1; y < xlength + 1; y++ )
  		for( x = 1; x < xlength + 1; x++ )
  			for( i = 0; i < 19; i++ ){
  				streamField[ 19 * (z * dim_pow2 + y * dim + x) + i ] = 
  				 collideField[ 19 * ( (z + LATTICEVELOCITIES[i][2]) * 
  				 	                 dim_pow2 + (y + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0]) + i ];
  			}

}

