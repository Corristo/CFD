#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
  int i;
  *density = 0.0;

  for(i = 0; i < 19; i++) (*density) = (*density) + currentCell[i];

}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
  int i;
  velocity[0] = velocity[1] = velocity[2] = 0.0;

  for(i = 0; i < 19; i++){
       velocity[0] = velocity[0] + LATTICEVELOCITIES[i][0] * currentCell[i];
       velocity[1] = velocity[1] + LATTICEVELOCITIES[i][1] * currentCell[i];
       velocity[2] = velocity[2] + LATTICEVELOCITIES[i][2] * currentCell[i];
  }

  velocity[0] = velocity[0] / (*density); 
  velocity[1] = velocity[1] / (*density); 
  velocity[2] = velocity[2] / (*density); 

}

void computeFeq(const double * const density, const double * const velocity, double *feq){
   int i;
   double pdt;
   double cs_2 = C_S * C_S, cs_4 = cs_2 * cs_2;

   for(i = 0; i < 19; i++){

   	    pdt = velocity[0] * LATTICEVELOCITIES[i][0] + velocity[1] * LATTICEVELOCITIES[i][1] + velocity[2] * LATTICEVELOCITIES[i][2];
   	    feq[i] = LATTICEWEIGHTS[i] * (*density) * ( 1 + pdt / (cs_2) + pdt * pdt / (double)(2 * cs_4) - 
   	    			(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]) / (double)(2 * cs_2) );

   }
}

