#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "boundary_val.h"
#include <stdio.h>
#include <time.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args)
{
    double **U, **V, **P, **R, **F, **G;
    double t = 0.0;
    double t_end, xlength, ylength, tau, omg, alpha, eps, Re, GX, GY, PI, UI, VI, dx, dy, dt, res, dt_value;
    int n, imax, jmax, itermax, it, i;

    char* paraviewOutput = "./Paraview/cavity100";
    char* inputFile = "./cavity100.dat";

    /* statistics */
    int sumIter, maxIter, minIter;
    double minTimeStep, maxTimestep;
    clock_t begin, end;
    double time_spent;
    sumIter = 0;
    maxIter = 0;
    maxTimestep = 0;


    begin = clock();

    /* Parse command line arguments */
    for (i = 1; i < argn; i++)
    {
        if (strcmp(args[i], "-i") == 0)
            inputFile = args[i + 1];

        if (strcmp(args[i], "-o") == 0)
            paraviewOutput = args[i + 1];


    }

     /* Initialization */
    read_parameters(inputFile, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value);
    U = matrix(0, imax, 0, jmax + 1);
    V = matrix(0, imax + 1, 0, jmax);
    P = matrix(0, imax + 1, 0, jmax + 1);
    R = matrix(0, imax + 1, 0, jmax + 1);
    F = matrix(0, imax, 0, jmax);
    G = matrix(0, imax, 0, jmax);
    n = 0;

    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
    minIter = itermax;
    minTimeStep = t_end;


    while (t < t_end)
    {
        calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
        if (dt < minTimeStep)
            minTimeStep = dt;
        if (dt > maxTimestep)
            maxTimestep = dt;
        boundaryvalues(imax, jmax, U, V);
        calculate_fg(Re, GX, GY, alpha, dt, dy, dy, imax, jmax, U, V, F, G);
        calculate_rs(dt, dx, dy, imax, jmax, F, G, R);
        it = 0;
        while(++it < itermax)
        {
            sor(omg, dx, dy, imax, jmax, P, R, &res);
            if (res <= eps)
                break;
        }
        if (res > eps)
            fprintf(stderr, "WARNING: SOR did not converge in timestep %d, residual: %f\n", n + 1, res);
        if (it > maxIter)
            maxIter = it;
        if (it < minIter)
            minIter = it;
        sumIter += it;
        calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
        t = t + dt;
        write_vtkFile(paraviewOutput,n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
        n++;
    }
    end = clock();
    time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

    /* print statistics */
    printf("Simulation took %f seconds\n", time_spent);
    printf("Number of timesteps: %d\n", n);
    printf("Min. timestep: %f     Max. timestep: %f    Avg. timestep: %f\n", minTimeStep, maxTimestep, t_end/n);
    printf("Min. # iterations: %d     Max. # iterations: %d     Avg: %f\n", minIter, maxIter, sumIter/(double) n);

    /* Free memory */
    free_matrix(U, 0, imax, 0, jmax + 1);
    free_matrix(V, 0, imax + 1, 0, jmax);
    free_matrix(P, 0, imax + 1, 0, jmax + 1);
    free_matrix(R, 0, imax, 0 ,jmax);
    free_matrix(F, 0, imax, 0, jmax);
    free_matrix(G, 0, imax, 0, jmax);


    return 0;
}
