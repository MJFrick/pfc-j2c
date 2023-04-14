#include <mpi.h>
#include <hdf5.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <stdbool.h>

#include "pfc_header.h"

#define WISDOM_FILE "fft_wisdom"

void load_wisdom(int argc,
                 char **argv)
{
    MPI_Init(&argc, &argv);

    // Initialize fftw
    fftw_mpi_init();

    // Import wisdom from file and broadcast
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0 && access("data/plans.wisdom", F_OK) != -1)
    {
        int err = fftw_import_wisdom_from_filename("data/plans.wisdom");
        if (err == 0)
            my_error("Importing FFTW wisdom failed!");
    }
    fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    return;
}

void save_wisdom(void)
{
    // Gather wisdom from procs and save to file
    fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        int err = fftw_export_wisdom_to_filename("data/plans.wisdom");
        if (err == 0)
        {
            remove("data/plans.wisdom");
            my_error("Failed to correctly export FFTW wisdom");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Clean up fftw and finalize MPI runtime
    fftw_mpi_cleanup();
    MPI_Finalize();
    return;
}
