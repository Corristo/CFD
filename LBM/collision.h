#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"
#include "LBDefinitions.h"

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(double *currentCell, const double * const tau, const double * const feq, int tot_cells);
void computePostCollisionDistributionsAVX (double *currentCell, const double * const tau, double feq[PARAMQ][4], int tot_cells);


/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(double *collideField,int *flagField,const double * const tau,int *xlength);
void doCollisionSSE(double *collideField, int *flagField,const double tau,int *xlength, double * densities, double ** velocities);
void doCollisionAVXv2(double *collideField, int *flagField,const double  tau, int *xlength, double * densities, double ** velocities);
void doCollisionAVX(double *collideField, int *flagField, const double * const tau, int *xlength);
#endif

