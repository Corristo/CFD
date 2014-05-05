#include "collision.h"
#include <math.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
  
  	int i;
  	double omg = *tau;

  	omg = 1.0 / omg;

  	for(i = 0; i < 19; i++)
		currentCell[i] = (1.0 - omg) * currentCell[i] + omg * feq[i];

}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
  
  int x, y, z;
  double density, velocity[3], feq[19];
  int dim = xlength + 2;
  int dim_pow2 = dim * dim;
  int index;

  for( z = 1; z < xlength + 1; z++ )
  	for( y = 1; y < xlength + 1; y++ )
  		for( x = 1; x < xlength + 1; x++ ){

  			index = 19 * (z * dim_pow2 + y * dim + x);

  			computeDensity(&collideField[index], &density);
  			computeVelocity(&collideField[index], &density, velocity);
  			computeFeq(&density, velocity, feq);
  			computePostCollisionDistributions(&collideField[index], tau, feq);
  		}
}

