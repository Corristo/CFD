#ifndef _COMMUNICATION_H_
#define _COMMUNICATION_H_
#include <mpi.h>

void extract( int * local_xlength, double * sendBuffer, double * collideField, int axis, int step );
void inject( int * local_xlength, double * readBuffer, int * flagField, double * collideField, int axis, int step );
void communicateBoundaryValues ( int * local_xlength, MPI_Datatype * MPISendTypes, MPI_Datatype * MPIRecvTypes, double * collideField, int * neighbours, int axis);
void computeNeighbours( int iProc, int jProc, int kProc, int * neighbours );
void computePosition( int iProc, int jProc, int kProc, int * iCoord, int * jCoord, int * kCoord );
void initialiseBuffers( int * local_xlength, double ** sendBuffer, double ** readBuffer, int * neighbours );
void computeLocalDomainSize( int * xlength, int * local_xlength, int iProc, int jProc, int kProc );
void getSendRecvCount( const int * const local_xlength, const int * const neighbours, const int * const flagField, int * sendCount, int * recvCount );
void getSendRecvIndices( const int * const local_xlength, const int * const neighbours, const int * const flagField, int ** sendIndices, int ** recvIndices );
void initialiseMPITypes( const int * const local_xlength, const int * const neighbours, const int * const flagField, MPI_Datatype * MPISendTypes, MPI_Datatype * MPIRecvTypes );
#endif
