/*
 * This is a part of Diffusion Equation MPI
 * Nathan Smith (c) 2016
 *
 * This code implements error handling. If an error is detected somewhere
 * an error can be thrown by calling error("<Error string here>").
 */

#include <mpi.h>
#include <stdio.h>

#include "pfc_header.h"

void my_error(const char *error_string)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("-------------------------------------------------------\n");
    printf("Error caught on process %d: %s\n", rank, error_string);
    printf("-------------------------------------------------------\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
}
