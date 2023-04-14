// Define dependencies
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "pfc_header.h"
state *create_state(Nx, Ny)
{
  ptrdiff_t local_alloc;
  state *s = malloc(sizeof(state));
  if (s == NULL)
    return NULL;
  local_alloc = fftw_mpi_local_size_2d_transposed(Nx, (Ny >> 1) + 1, MPI_COMM_WORLD,
                                                  &s->local_n0, &s->local_0_start,
                                                  &s->local_n1, &s->local_1_start);
  s->n = fftw_alloc_real(2 * local_alloc);
  s->NL = fftw_alloc_real(2 * local_alloc);

  s->nk = fftw_alloc_complex(local_alloc);
  s->NL_k = fftw_alloc_complex(local_alloc);

  s->k2 = fftw_alloc_real(local_alloc);
  s->C2 = fftw_alloc_real(local_alloc);
  s->P = fftw_alloc_real(local_alloc);
  s->Q = fftw_alloc_real(local_alloc);

  s->Nx = Nx;
  s->Ny = Ny;

  s->step = 0;

  s->dt = 0.0;
  s->dx = 0.0;
  s->dy = 0.0;
  s->Bx = 0.0;
  s->t = 0.0;
  s->v = 0.0;
  s->M = 0.0;
  s->tm = 0.0;

  if (s->n == NULL ||
      s->NL == NULL ||

      s->nk == NULL ||
      s->NL_k == NULL ||

      s->k2 == NULL ||
      s->C2 == NULL ||
      s->P == NULL ||
      s->Q == NULL ||

      false)
  {
    free(s->n);
    free(s->NL);

    free(s->nk);
    free(s->NL_k);

    free(s->k2);
    free(s->C2);
    free(s->P);
    free(s->Q);

    return NULL;
  }
  s->forward = fftw_mpi_plan_dft_r2c_2d(Nx, Ny, s->n, s->nk, MPI_COMM_WORLD,
                                        FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
  s->back = fftw_mpi_plan_dft_c2r_2d(Nx, Ny, s->nk, s->n, MPI_COMM_WORLD,
                                     FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
  // Nathan style entropy
  s->rng = gsl_rng_alloc(gsl_rng_default);
  if (s->rng == NULL)
  {
    gsl_rng_free(s->rng);
    return NULL;
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  gsl_rng_set(s->rng, (
                          (unsigned long)time(NULL) + (unsigned long)clock() + (unsigned long)getpid() + (unsigned long)getppid()) *
                          (rank + 1));

  // Barrier to avoid pipeline issues.
  MPI_Barrier(MPI_COMM_WORLD);

  return s;
}
void destroy_state(state *s)
{
  // Free memory
  if (s != NULL)
  {
    fftw_destroy_plan(s->forward);
    fftw_destroy_plan(s->back);
    free(s->n);
    free(s->NL);

    free(s->nk);
    free(s->NL_k);

    free(s->k2);
    free(s->C2);
    free(s->P);
    free(s->Q);

    free(s);
  }
}