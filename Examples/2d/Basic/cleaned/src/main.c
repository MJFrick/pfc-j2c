/*
 *      Simulation from Input file
 *
 * Usage:
 *  >> mpiexec -np <procs> ./main <inputfile> <outputfile> <final_t> <save_every>
 */

// Standard libraries
#include <mpi.h>
#include <hdf5.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>

// Local libraries
#include "pfc_header.h"

#define PI 2.0 * acos(0.0)

int main(int argc, char **argv)
{
  // define the file id and blank state
  hid_t file_id;
  state *s;

  // Convert inputs to numbers
  int index;
  double t_stop = atof(argv[3]);
  double sv_step = atof(argv[4]);
  double next_sv = sv_step;
  // double angle;

  char string[50];
  // mpi_print("loading wisdom");

  // load any fft wisdom that may exist
  load_wisdom(argc, argv);
  mpi_print("creating state");

  // load the state from file
  s = new_state_from_file(argv[1]);
  MPI_Barrier(MPI_COMM_WORLD);
  // mpi_print("calculating k-space");

  // Calc_halfReals(s);

  if (isnan(s->dx))
  {
    mpi_print("Failed to load file...");
    return 1;
  }
  mpi_print("Loaded File");

  mpi_print("Seeded State");
  // // initialize output file
  file_id = io_init_new_file(argv[2]);

  // seed

  seed(s, 1.0, 0.3);

  // Save initial State
  save_state(s, file_id);
  io_finalize(file_id);
  file_id = io_init_from_file(argv[2]);

  // sanity check
  mpi_print("\nStarting Simulation.\n");

  // run
  while (s->tm < t_stop)
  {

    MPI_Barrier(MPI_COMM_WORLD);
    step(s);
    MPI_Barrier(MPI_COMM_WORLD);
    if (isnan(s->n[0]))
    {
      mpi_print("NaN encountered...");
      io_finalize(file_id);
      destroy_state(s);
      return 1;
    }

    if (s->tm >= next_sv)
    {
      snprintf(string, 50, "\tSaving State\tt = %g at timestep = %i", s->tm, s->step);
      mpi_print(string);

      save_state(s, file_id);
      io_finalize(file_id);
      file_id = io_init_from_file(argv[2]);
      next_sv += sv_step;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Finish up
  mpi_print("Finishing up...");
  io_finalize(file_id);
  destroy_state(s);
  save_wisdom();

  return 0;
}
