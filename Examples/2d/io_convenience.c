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