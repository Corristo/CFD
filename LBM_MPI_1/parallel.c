#include "parallel.h"
#include <mpi.h>

void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void finalizeMPI()
{
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
}

void initializeMPI(int * rank,int * number_of_ranks, int argc , char ** argv)
{
   MPI_Init( &argc, &argv );                    /* execute n processes      */
   MPI_Comm_size( MPI_COMM_WORLD, &(*number_of_ranks) );     /* asking for the number of processes  */
   MPI_Comm_rank( MPI_COMM_WORLD, &(*rank) );    /* asking for the local process id   */
}
