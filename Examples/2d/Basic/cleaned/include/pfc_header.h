// dependencies
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <hdf5.h>
#include <gsl/gsl_rng.h>
#define PI 2.0*acos(0.0)
typedef struct
{
double* n;
double* NL;

fftw_complex* nk;
fftw_complex* NL_k;

double* k2;
double* C2;
double* P;
double* Q;

ptrdiff_t Nx;
ptrdiff_t Ny;

int step;

double dt;
double dx;
double dy;
double Bx;
double t;
double v;
double M;
double tm;

fftw_plan forward;
fftw_plan back;


ptrdiff_t local_n0;       // local endpoint of a processor
ptrdiff_t local_0_start;  // real start of a processor
ptrdiff_t local_n1;       // Used for transposed fftw
ptrdiff_t local_1_start;  // Same.
gsl_rng* rng;
} state;

void Calc_halfReals(state* s);
void Calc_NL(state* s);
void step(state* s);
void seed(state* s, double q, double A);

state* create_state(int Nx, int Ny);
void destroy_state(state *s);
herr_t save_state(state *s, hid_t file_id);
state* new_state_from_file(const char *filename);
state* load_state(hid_t file_id, const char *datafile);
void mpi_print(const char *str);
hid_t io_init_from_file(const char *filename);
hid_t io_init_new_file(const char *filename);
herr_t io_finalize(hid_t file_id);