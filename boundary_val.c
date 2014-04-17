/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
    int imax,
    int jmax,
    double **U,
    double **V
)
{
    int i,j;
    for (i = 0; i <= imax; i++)
    {
        for (j = 0; j <= jmax; j++)
        {
            if ((i == 0) && (j > 0))
            {
                U[0][j] = 0;
                V[0][j] = -V[1][j];
            }

            if ((i == imax) && (j > 0))
            {
                U[imax][j] = 0;
                V[imax + 1][j] = - V[imax][j];
            }

            if ((j == 0) && (i > 0))
            {
                U[i][0] = -U[i][1];
                V[i][0] = 0;
            }

            if ((j == jmax) && (i >0))
            {
                U[i][jmax + 1] = 2.0-U[i][jmax];
                V[i][jmax] = 0;
            }
        }
    }

    /* Cosmetics for a nicer structure of the matrices */
    /* The values set here are never accessed */
    U[0][0] = 0;
    U[imax][0] = 0;
    U[imax][jmax + 1] = 2.0;
    U[0][jmax + 1] = 0.0;
    V[0][0] = 0;
    V[imax + 1][0] = 0;
    V[imax + 1][jmax] = 0;
    V[0][jmax] = 0;
}
