#include <mpi.h>        // MPI runtime
#include <hdf5.h>       // HDF5 header
#include <stdlib.h>     // malloc etc
#include <unistd.h>     // access() to check if file exists
#include <stdbool.h>    // Booleans
#include <assert.h>     // Assertions
#include <string.h>     // String manipulation sprintf etc...
#include <math.h>       // Needed for sqrt

#include "pfc_header.h"    // local libbinary header
#include "io_convenience.c"

#define FAIL -1
#define LEN 8

const int io_verbose = true;
herr_t save_state (state *s,
   hid_t  file_id)
{
     hid_t group_id;
     herr_t status;

     /* Make Group from simulation time `t` */

     char groupname[50];
     char step_str[10];
     sprintf(step_str, "%d", s->step);
     sprintf(groupname, "%0*d%s", LEN-(int)strlen(step_str), 0, step_str);

     group_id = H5Gcreate (file_id,
                           groupname,
                           H5P_DEFAULT,
                           H5P_DEFAULT,
                           H5P_DEFAULT);status = write_array_dataset("n", group_id, s->n, s);
status = write_array_dataset("NL", group_id, s->NL, s);


status = write_Rek_array_dataset("k2", group_id, s->k2, s);
status = write_Rek_array_dataset("C2", group_id, s->C2, s);
status = write_Rek_array_dataset("P", group_id, s->P, s);
status = write_Rek_array_dataset("Q", group_id, s->Q, s);

int Nx_int = (int)s->Nx; 
 status = write_int_attribute("Nx", group_id, &s->Nx);
int Ny_int = (int)s->Ny; 
 status = write_int_attribute("Ny", group_id, &s->Ny);

status = write_int_attribute("step", group_id, &s->step);

status = write_double_attribute("dt", group_id, &s->dt);
status = write_double_attribute("dx", group_id, &s->dx);
status = write_double_attribute("dy", group_id, &s->dy);
status = write_double_attribute("Bx", group_id, &s->Bx);
status = write_double_attribute("t", group_id, &s->t);
status = write_double_attribute("v", group_id, &s->v);
status = write_double_attribute("M", group_id, &s->M);
status = write_double_attribute("tm", group_id, &s->tm);



status = H5Gclose (group_id);
return status;
}
state* load_state (hid_t       file_id,
   const char *datafile)
{
    state* s;
    hid_t group_id;
    herr_t status;

    group_id = H5Gopen2 (file_id, datafile, H5P_DEFAULT);int Nx_int;
 read_int_attribute ("Nx", group_id, &Nx_int);
 ptrdiff_t Nx = (ptrdiff_t)Nx_int;
int Ny_int;
 read_int_attribute ("Ny", group_id, &Ny_int);
 ptrdiff_t Ny = (ptrdiff_t)Ny_int;
s = create_state(Nx, Ny);
assert (s != NULL);status = read_array_dataset("n", group_id, &s->n, s);
status = read_array_dataset("NL", group_id, &s->NL, s);


status = read_array_dataset("k2", group_id, &s->k2, s);
status = read_array_dataset("C2", group_id, &s->C2, s);
status = read_array_dataset("P", group_id, &s->P, s);
status = read_array_dataset("Q", group_id, &s->Q, s);


status = read_int_attribute("step", group_id, &s->step);

status = read_double_attribute("dt", group_id, &s->dt);
status = read_double_attribute("dx", group_id, &s->dx);
status = read_double_attribute("dy", group_id, &s->dy);
status = read_double_attribute("Bx", group_id, &s->Bx);
status = read_double_attribute("t", group_id, &s->t);
status = read_double_attribute("v", group_id, &s->v);
status = read_double_attribute("M", group_id, &s->M);
status = read_double_attribute("tm", group_id, &s->tm);



status = H5Gclose (group_id);
return s;
}
