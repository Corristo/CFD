#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, int* flagField, const double * const bddParams,int *xlength);

/** compute boundary helper function */
void compute_boundary(double *collideField, const double * const bddParams, int *flagField,
                         int boundaryType, int *xlength, int *iList, int *const coordinate);
#endif

