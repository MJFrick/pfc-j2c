#include <mpi.h>     // MPI runtime
#include <hdf5.h>    // HDF5 header
#include <stdlib.h>  // malloc etc
#include <unistd.h>  // access() to check if file exists
#include <stdbool.h> // Booleans
#include <assert.h>  // Assertions
#include <string.h>  // String manipulation sprintf etc...
#include <math.h>    // Needed for sqrt

#include "pfc_header.h" // local libbinary header

#define FAIL -1
#define LEN 8

const int io_verbose = true;


state *new_state_from_file(const char *filename)
{
   FILE *ifp;
   char *mode = "r";
   state *s;
   char name[20];
   double value;
   int Nx, Ny;
   double dt, dx, dy;
   double n0, c1, c2;

   ifp = fopen(filename, mode);
   assert(ifp != NULL);

   mpi_print("Getting system size");

   while (fscanf(ifp, "%s : %lf", name, &value) != EOF)
   {
      if (strcmp(name, "Nx") == 0)
         Nx = (int)value;
      else if (strcmp(name, "Ny") == 0)
         Ny = (int)value;
   }
   // // initialize output file

   s = create_state(Nx, Ny);

   mpi_print("Creating state");

   assert(s != NULL);
   rewind(ifp);

   mpi_print("Populating attributes");

   char string[50];

   while (fscanf(ifp, "%s : %lf", name, &value) != EOF)
   {
      if (strcmp(name, "dx") == 0)
         s->dx = value;
      else if (strcmp(name, "dy") == 0)
         s->dy = value;
      else if (strcmp(name, "dt") == 0)
         s->dt = value;
      else if (strcmp(name, "Bx") == 0)
         s->Bx = value;
      else if (strcmp(name, "t") == 0)
         s->t = value;
      else if (strcmp(name, "v") == 0)
         s->v = value;
      else if (strcmp(name, "M") == 0)
         s->M = value;

      else if (strcmp(name, "n0") == 0)
         n0 = value;
   }
   MPI_Barrier(MPI_COMM_WORLD);

   int index;
   for (int i = 0; i < s->local_n0; i++)
   {
      for (int j = 0; j < 2 * ((s->Ny >> 1) + 1); j++)
      {
         index = j + 2 * ((s->Ny >> 1) + 1) * i;
         s->n[index] = n0;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   s->step = 0;
   s->tm = 0.0;

   Calc_halfReals(s);
   MPI_Barrier(MPI_COMM_WORLD);
   return s;
}

herr_t save_state(state *s,
                  hid_t file_id)
{
   hid_t group_id;
   herr_t status;

   /* Make Group from simulation time `t` */

   char groupname[50];
   char step_str[10];
   sprintf(step_str, "%d", s->step);
   sprintf(groupname, "%0*d%s", LEN - (int)strlen(step_str), 0, step_str);

   group_id = H5Gcreate(file_id,
                        groupname,
                        H5P_DEFAULT,
                        H5P_DEFAULT,
                        H5P_DEFAULT);
   status = write_array_dataset("n", group_id, s->n, s);
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

   status = H5Gclose(group_id);
   return status;
}
state *load_state(hid_t file_id,
                  const char *datafile)
{
   state *s;
   hid_t group_id;
   herr_t status;

   group_id = H5Gopen2(file_id, datafile, H5P_DEFAULT);
   int Nx_int;
   read_int_attribute("Nx", group_id, &Nx_int);
   ptrdiff_t Nx = (ptrdiff_t)Nx_int;
   int Ny_int;
   read_int_attribute("Ny", group_id, &Ny_int);
   ptrdiff_t Ny = (ptrdiff_t)Ny_int;
   s = create_state(Nx, Ny);
   assert(s != NULL);
   status = read_array_dataset("n", group_id, &s->n, s);
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

   status = H5Gclose(group_id);
   return s;
}

void mpi_print (const char *str)
{
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0) printf("%s\n", str);
    MPI_Barrier (MPI_COMM_WORLD);
}

hid_t io_init_from_file (const char *filename)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t plist_id;
    hid_t file_id;
    herr_t status;

    plist_id = H5Pcreate (H5P_FILE_ACCESS);

    H5Pset_fapl_mpio (plist_id, comm, info);

    file_id = H5Fopen (filename, H5F_ACC_RDWR, plist_id);

    status = H5Pclose (plist_id);
    assert (status != FAIL);

    return file_id;
}


hid_t io_init_new_file (const char *filename)
{
    /*
    *  Create a file to save data to for this session
    */
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t plist_id;     /* Property list id */
    hid_t file_id;      /* File id */
    herr_t err;         /* Error status */

    plist_id = H5Pcreate (H5P_FILE_ACCESS);

    H5Pset_fapl_mpio (plist_id, comm, info);

    file_id = H5Fcreate (filename,
                         H5F_ACC_TRUNC,
                         H5P_DEFAULT,
                         plist_id );

    err = H5Pclose (plist_id);
    assert (err != FAIL);
    return file_id;
}

herr_t io_finalize (hid_t file_id)
{
    return H5Fclose(file_id);
}

herr_t write_array_dataset (const char *name,
                            hid_t       group_id,
                            double     *arr,
                            state      *s)
{
     hid_t dset_id, dataspace;
     hid_t memspace, plist_id;
     herr_t status;
     hsize_t size[2], count[2], offset[2];

     /* Describe the shape of the array */
     size[0] = s->Nx;
     size[1] = s->Ny;
     dataspace = H5Screate_simple(2, size, NULL);


      /* Create Dataset */
     dset_id = H5Dcreate(group_id,
                         name,
                         H5T_NATIVE_DOUBLE,
                         dataspace,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

     /* Describe memory shape, property list and data shape */
     count[0] = 1;
     count[1] = s->Ny;
     offset[1] = 0;
     memspace = H5Screate_simple (2, count, NULL);

     /* Set up some of the MPI things */
     plist_id = H5Pcreate (H5P_DATASET_XFER);
     H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_INDEPENDENT);

     /* Write data row by row in slabs */
     for (int row = 0; row < s->local_n0; row++)
     {
         offset[0] = s->local_0_start + row;

         status = H5Sselect_hyperslab (dataspace,
                                       H5S_SELECT_SET,
                                       offset,
                                       NULL,
                                       count,
                                       NULL);

         status = H5Dwrite (dset_id,
                            H5T_NATIVE_DOUBLE,
                            memspace,
                            dataspace,
                            plist_id,
                            arr + row*2*((s->Ny>>1) + 1));
     }

     /* Close everything you opened */
     status = H5Pclose (plist_id);
     status = H5Sclose (memspace);
     status = H5Dclose (dset_id);
     status = H5Sclose (dataspace);

     return status;
}

herr_t write_Rek_array_dataset (const char *name,
                            hid_t       group_id,
                            double     *arr,
                            state      *s)
{
     hid_t dset_id, dataspace;
     hid_t memspace, plist_id;
     herr_t status;
     hsize_t size[2], count[2], offset[2];

     /* Describe the shape of the array */
     size[0] = (s->Ny>>1) + 1;
     size[1] = s->Nx;

     dataspace = H5Screate_simple(2, size, NULL);

      /* Create Dataset */
     dset_id = H5Dcreate(group_id,
                         name,
                         H5T_NATIVE_DOUBLE,
                         dataspace,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

     /* Describe memory shape, property list and data shape */
     count[0] = 1;
     count[1] = s->Nx;
     offset[1] = 0;
     memspace = H5Screate_simple (2, count, NULL);

     /* Set up some of the MPI things */
     plist_id = H5Pcreate (H5P_DATASET_XFER);
     H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_INDEPENDENT);

     /* Write data row by row in slabs */
     for (int row = 0; row < s->local_n1; row++)
     {
         offset[0] = s->local_1_start + row;

         status = H5Sselect_hyperslab (dataspace,
                                       H5S_SELECT_SET,
                                       offset,
                                       NULL,
                                       count,
                                       NULL);

         status = H5Dwrite (dset_id,
                            H5T_NATIVE_DOUBLE,
                            memspace,
                            dataspace,
                            plist_id,
                            arr + row*s->Nx);
     }

     /* Close everything you opened */
     status = H5Pclose (plist_id);
     status = H5Sclose (memspace);
     status = H5Dclose (dset_id);
     status = H5Sclose (dataspace);

     return status;
}

herr_t read_array_dataset (const char *name,
                           hid_t       group_id,
                           double     *arr,
                           state      *s)
{
    hid_t dset_id, dataspace;
    hid_t memspace, plist_id;
    herr_t status;
    hsize_t count[2], offset[2];

    /* Open Dataset */
    dset_id = H5Dopen1 (group_id, name);

    /* Describe memory shape, property list and data shape */
    count[0] = 1;
    count[1] = s->Ny;
    offset[1] = 0;
    memspace = H5Screate_simple (2, count, NULL);

    /* Set up some of the MPI things */
    dataspace = H5Dget_space (dset_id);
    plist_id = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_INDEPENDENT);

    /* Write data row by row in slabs */
    for (int row = 0; row < s->local_n0; row++)
    {
        offset[0] = s->local_0_start + row;

        status = H5Sselect_hyperslab (dataspace,
                                      H5S_SELECT_SET,
                                      offset,
                                      NULL,
                                      count,
                                      NULL);

        status = H5Dread (dset_id,
                          H5T_NATIVE_DOUBLE,
                          memspace,
                          dataspace,
                          plist_id,
                          arr + row*2*(s->Ny/2 + 1));
    }

    /* Close everything you opened */
    status = H5Pclose (plist_id);
    status = H5Sclose (memspace);
    status = H5Dclose (dset_id);
    status = H5Sclose (dataspace);

    return status;
}

herr_t write_double_attribute (const char *name,
                               hid_t       group_id,
                               double     *value)
{
     hsize_t size = 1;
     herr_t status;
     hid_t attr_id, dataspace;

     dataspace = H5Screate_simple(1, &size, NULL);
     attr_id = H5Acreate2 (group_id,
                           name,
                           H5T_NATIVE_DOUBLE,
                           dataspace,
                           H5P_DEFAULT,
                           H5P_DEFAULT);

     status = H5Awrite (attr_id, H5T_NATIVE_DOUBLE, value);
     status = H5Aclose (attr_id);

     return status;
}

herr_t read_double_attribute (const char *name,
                              hid_t       group_id,
                              double     *value)
{
    /* Read integer attribute from dataset 'group_id' */
    hid_t attr_id;
    herr_t status;

    attr_id = H5Aopen (group_id, name, H5P_DEFAULT);
    status = H5Aread (attr_id, H5T_NATIVE_DOUBLE, value);

    status = H5Aclose (attr_id);
    return status;
}

herr_t write_int_attribute (const char *name,
                            hid_t       group_id,
                            int        *value)
{
     hsize_t size = 1;
     herr_t status;
     hid_t attr_id, dataspace;

     dataspace = H5Screate_simple (1, &size, NULL);
     attr_id = H5Acreate2 (group_id,
                           name,
                           H5T_NATIVE_INT,
                           dataspace,
                           H5P_DEFAULT,
                           H5P_DEFAULT);

     status = H5Awrite (attr_id, H5T_NATIVE_INT, value);
     status = H5Aclose (attr_id);

     return status;
}

herr_t read_int_attribute (const char *name,
                           hid_t       group_id,
                           int        *value)
{
    /* Read integer attribute from dataset 'group_id' */
    hid_t attr_id;
    herr_t status;

    attr_id = H5Aopen (group_id, name, H5P_DEFAULT);
    status = H5Aread (attr_id, H5T_NATIVE_INT, value);

    status = H5Aclose (attr_id);

    return status;
}