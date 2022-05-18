/************************************************************

  This example shows how to read and write data to a dataset
  using gzip compression (also called zlib or deflate).  The
  program first checks if gzip compression is available,
  then if it is it writes integers to a dataset using gzip,
  then closes the file.  Next, it reopens the file, reads
  back the data, and outputs the type of compression and the
  maximum value in the dataset to the screen.

  This file is intended for use with HDF5 Library version 1.8

 ************************************************************/
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include "hdf5.h"

#define FILE            "h5ex_d_compress.h5"
#define DATASET         "DS1"

/**** PROGRAM PARAMETERS *****/
//#define DIM0            65536 //8388608 /* data buffer size 0-dir */
//#define CHUNK0          16384 /* chunk size in 0-dir = DIM0/CHUNK0 */
#define DIM0            256 //8388608 /* data buffer size 0-dir */
#define CHUNK0           16 /* chunk size in 0-dir = DIM0/CHUNK0 */
/**** PROGRAM PARAMETERS *****/

#define RANK            1

#define DEBUG 0

/* Generate a simple, 1D sinusioidal data array with some noise */
#define TYPINT 1
#define TYPDBL 2
static int gen_data(size_t npoints, double noise, double amp, void **_buf, int typ)
{
    size_t i;
    double *pdbl = 0;
    int *pint = 0;

    /* create data buffer to write */
    if (typ == TYPINT)
        pint = (int *) malloc(npoints * sizeof(int));
    else
        pdbl = (double *) malloc(npoints * sizeof(double));
    srandom(0xDeadBeef);
    for (i = 0; i < npoints; i++)
    {
        double x = 2 * M_PI * (double) i / (double) (npoints-1);
        double n = noise * ((double) random() / ((double)(1<<31)-1) - 0.5);
        if (typ == TYPINT)
            pint[i] = (int) (amp * (1 + sin(x)) + n);
        else
            pdbl[i] = (double) (amp * (1 + sin(x)) + n);
    }
    if (typ == TYPINT)
        *_buf = pint;
    else
        *_buf = pdbl;
    return 0;
}

/* Populate the hyper-dimensional array with samples of a radially symmetric
   sinc() function but where certain sub-spaces are randomized through dimindx arrays */
static void
hyper_smooth_radial(void *b, int typ, int n, int ndims, int const *dims, int const *m,
    int const * const dimindx[10])
{
    int i;
    double hyper_radius = 0;
    const double amp = 10000;
    double val;

    for (i = ndims-1; i >= 0; i--)
    {
        int iar = n / m[i];
        iar = dimindx[i][iar]; /* allow for randomized shuffle of this axis */
        iar -= dims[i]/2;      /* ensure centering in middle of the array */
        n = n % m[i];
        hyper_radius += iar*iar;
    }
    hyper_radius = sqrt(hyper_radius);

    if (hyper_radius < 1e-15)
        val = amp;
    else
        val = amp * sin(0.4*hyper_radius) / (0.4*hyper_radius);

    if (typ == TYPINT)
    {
        int *pi = (int*) b;
        *pi = (int) val;
    }
    else
    {
        double *pd = (double*) b;
        *pd = val;
    }
}

static double func(int i, double arg)
{
    /* a random assortment of interesting, somewhat bounded, unary functions */
    double (*const funcs[])(double x) = {cos, j0, fabs, sin, cbrt, erf};
    int const nfuncs = sizeof(funcs)/sizeof(funcs[0]);
    return funcs[i%nfuncs](arg);
}

/* Populate the hyper-dimensional array with samples of set of seperable functions
   but where certain sub-spaces are randomized through dimindx arrays */
static void
hyper_smooth_separable(void *b, int typ, int n, int ndims, int const *dims, int const *m,
    int const * const dimindx[10])
{
    int i;
    double val = 1;

    for (i = ndims-1; i >= 0; i--)
    {
        int iar = n / m[i];
        iar = dimindx[i][iar]; /* allow for randomized shuffle of this axis */
        iar -= dims[i]/2;      /* ensure centering in middle of the array */
        n = n % m[i];
        val *= func(i, (double) iar);
    }

    if (typ == TYPINT)
    {
        int *pi = (int*) b;
        *pi = (int) val;
    }
    else
    {
        double *pd = (double*) b;
        *pd = val;
    }
}

/* Produce multi-dimensional array test data with the property that it is random
   in the UNcorrelated dimensions but smooth in the correlated dimensions. This
   is achieved by randomized shuffling of the array indices used in specific
   dimensional axes of the array. */
static void *
gen_random_correlated_array(int typ, int ndims, int const *dims, int nucdims, int const *ucdims)
{
    int i, n;
    int nbyt = (int) (typ == TYPINT ? sizeof(int) : sizeof(double)); 
    unsigned char *buf, *buf0;
    int m[10]; /* subspace multipliers */
    int *dimindx[10];
   
    assert(ndims <= 10);

    /* Set up total size and sub-space multipliers */
    for (i=0, n=1; i < ndims; i++)
    {
        n *= dims[i];
        m[i] = i==0?1:m[i-1]*dims[i-1];
    }

    /* allocate buffer of suitable size (doubles or ints) */
    buf0 = buf = (unsigned char*) malloc(n * nbyt);
    
    /* set up dimension identity indexing (e.g. Idx[i]==i) so that
       we can randomize those dimenions we wish to have UNcorrelated */
    for (i = 0; i < ndims; i++)
    {
        int j;
        dimindx[i] = (int*) malloc(dims[i]*sizeof(int));
        for (j = 0; j < dims[i]; j++)
            dimindx[i][j] = j;
    }

    /* Randomize selected dimension indexing */
    srandom(0xDeadBeef);
    for (i = 0; i < nucdims; i++)
    {
        int j, ucdimi = ucdims[i];
        for (j = 0; j < dims[ucdimi]-1; j++)
        {
            int tmp, k = random() % (dims[ucdimi]-j);
            if (k == j) continue;
            tmp = dimindx[ucdimi][j];
            dimindx[ucdimi][j] = k;
            dimindx[ucdimi][k] = tmp;
        }
    }

    /* populate the array data */
    for (i = 0; i < n; i++)
    {
        hyper_smooth_separable(buf, typ, i, ndims, dims, m, (int const * const *) dimindx);
        buf += nbyt;
    }

    /* free dimension indexing */
    for (i = 0; i < ndims; i++)
        free(dimindx[i]);

    return buf0;
}

static void
modulate_by_time(void *data, int typ, int ndims, int const *dims, int t)
{
    int i, n;

    for (i = 0, n = 1; i < ndims; i++)
        n *= dims[i];

    if (typ == TYPINT)
    {
        int *p = (int *) data;
        for (i = 0; i < n; i++, p++)
        {
            double val = *p;
            val *= exp(0.1*t*sin(t/9.0*2*M_PI));
            *p = val;
        }
    }
    else
    {
        double *p = (double *) data;
        for (i = 0; i < n; i++, p++)
        {
            double val = *p;
            val *= exp(0.1*t*sin(t/9.0*2*M_PI));
            *p = val;
        }
    }
}

static void
buffer_time_step(void *tbuf, void *data, int typ, int ndims, int const *dims, int t)
{
    int i, n;
    int k = t % 4;
    int nbyt = (int) (typ == TYPINT ? sizeof(int) : sizeof(double)); 

    for (i = 0, n = 1; i < ndims; i++)
        n *= dims[i];

    memcpy((char*)tbuf+k*n*nbyt, data, n*nbyt);
}

#if 0
static int read_data(char const *fname, size_t npoints, double **_buf)
{
    size_t const nbytes = npoints * sizeof(double);
    int fd;

    if (0 > (fd = open(fname, O_RDONLY))) ERROR(open);
    if (0 == (*_buf = (double *) malloc(nbytes))) ERROR(malloc);
    if (nbytes != read(fd, *_buf, nbytes)) ERROR(read);
    if (0 != close(fd)) ERROR(close);
    return 0;
}
#endif

#if 0
static hid_t setup_filter(int n, hsize_t *chunk, int zfpmode,
    double rate, double acc, uint prec,
    uint minbits, uint maxbits, uint maxprec, int minexp)
{
    hid_t cpid;
    unsigned int cd_values[10];
    int i;
    size_t cd_nelmts = 10;

    /* setup dataset creation properties */
    if (0 > (cpid = H5Pcreate(H5P_DATASET_CREATE))) ERROR(H5Pcreate);
    if (0 > H5Pset_chunk(cpid, n, chunk)) ERROR(H5Pset_chunk);

#ifdef H5Z_ZFP_USE_PLUGIN
    /* setup zfp filter via generic (cd_values) interface */
    if (zfpmode == H5Z_ZFP_MODE_RATE)
        H5Pset_zfp_rate_cdata(rate, cd_nelmts, cd_values);
    else if (zfpmode == H5Z_ZFP_MODE_PRECISION)
        H5Pset_zfp_precision_cdata(prec, cd_nelmts, cd_values);
    else if (zfpmode == H5Z_ZFP_MODE_ACCURACY)
        H5Pset_zfp_accuracy_cdata(acc, cd_nelmts, cd_values);
    else if (zfpmode == H5Z_ZFP_MODE_EXPERT)
        H5Pset_zfp_expert_cdata(minbits, maxbits, maxprec, minexp, cd_nelmts, cd_values);
    else if (zfpmode == H5Z_ZFP_MODE_REVERSIBLE)
        H5Pset_zfp_reversible_cdata(cd_nelmts, cd_values);
    else
        cd_nelmts = 0; /* causes default behavior of ZFP library */

    /* print cd-values array used for filter */
    printf("%d cd_values= ", (int) cd_nelmts);
    for (i = 0; i < (int) cd_nelmts; i++)
        printf("%u,", cd_values[i]);
    printf("\n");

    /* Add filter to the pipeline via generic interface */
    if (0 > H5Pset_filter(cpid, H5Z_FILTER_ZFP, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values)) ERROR(H5Pset_filter);

#else 

    /* When filter is used as a library, we need to init it */
    H5Z_zfp_initialize();

    /* Setup the filter using properties interface. These calls also add
       the filter to the pipeline */
    if (zfpmode == H5Z_ZFP_MODE_RATE)
        H5Pset_zfp_rate(cpid, rate);
    else if (zfpmode == H5Z_ZFP_MODE_PRECISION)
        H5Pset_zfp_precision(cpid, prec);
    else if (zfpmode == H5Z_ZFP_MODE_ACCURACY)
        H5Pset_zfp_accuracy(cpid, acc);
    else if (zfpmode == H5Z_ZFP_MODE_EXPERT)
        H5Pset_zfp_expert(cpid, minbits, maxbits, maxprec, minexp);
    else if (zfpmode == H5Z_ZFP_MODE_REVERSIBLE)
        H5Pset_zfp_reversible(cpid);

#endif

    return cpid;
}
#endif
int
main (int argc, char **argv)
{
    hid_t           file, dset, dcpl;    /* Handles */
     hid_t          filespace, memspace;      /* file and memory dataspace identifiers */
    herr_t          status;
    htri_t          avail;
    hsize_t         chunk[1];
    double          *wdata = 0; /* pointer to data buffer to write */
    double          *rdata; 
    hsize_t         dimsf[1];                 /* dataset dimensions */
    hsize_t         i;
    int             max;
    hsize_t	    count[1];	          /* hyperslab selection parameters */
    hsize_t	    offset[1];
    hid_t	    plist_id;                 /* property list identifier */
    unsigned int    filter_info;

    /*
     * MPI variables
     */
    int mpi_size, mpi_rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    /*
     * Initialize MPI
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    /* 
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    H5Pset_fclose_degree(plist_id,H5F_CLOSE_WEAK);

    if(mpi_rank == 0) {
      printf("chunk size (MB) = %d\n ", DIM0/CHUNK0*sizeof(double)/1048576);
      printf("%d %d \n",DIM0/CHUNK0*sizeof(double),1048576);
    }
    
    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
   

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = DIM0;
    filespace = H5Screate_simple(RANK, dimsf, NULL); 

    /*
     * Create the dataset creation property list, add the gzip
     * compression filter and set the chunk size.
     */
    
    chunk[0] = DIM0/CHUNK0;

    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_layout(dcpl, H5D_CHUNKED);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    //    status = H5Pset_deflate (dcpl, 9);
    status = H5Pset_chunk (dcpl, 1, chunk);

    /*
     * Create the dataset and close filespace.
     */
    dset = H5Dcreate (file, DATASET, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcpl,
                H5P_DEFAULT);
    status = H5Sclose (filespace);
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = dimsf[0]/mpi_size;
    offset[0] = mpi_rank * count[0];
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Initialize data.
     */
    wdata = gen_random_correlated_array(TYPDBL, 1, (int*)count, 0, 0);

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    /*
     * Write the data to the dataset.
     */

    //   status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, wdata);

    int nchunks = count[0]/chunk[0];
    hsize_t offset_chunk[1];
    for (i = 0; i < nchunks; i++) {
      offset_chunk[0] = count[0]*mpi_rank + chunk[0]*i;
      status = H5Dwrite_chunk(dset, plist_id,  0, offset_chunk, chunk[0]*sizeof(double), &wdata[count[0]*i]);
    }

    free(wdata);
    H5Pclose(plist_id);

    /*
     * Close and release resources.
     */
    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
    status = H5Fclose (file);
#if 0
    /*
     * Now we begin the read section of this example.
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen (FILE, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);


    dset = H5Dopen (file, DATASET, H5P_DEFAULT);

    /*
     * Retrieve dataset creation property list.
     */
    dcpl = H5Dget_create_plist (dset);

#if 0 // DEBUG
    /*
     * Retrieve and print the filter type.  Here we only retrieve the
     * first filter because we know that we only added one filter.
     */
    size_t nelmts = 0;
    unsigned int flags;
    H5Z_filter_t filter_type = H5Pget_filter (dcpl, 0, &flags, &nelmts, NULL, 0, NULL,
                &filter_info);

    if(mpi_rank == 0) {
      printf ("Filter type is: ");
      switch (filter_type) {
      case H5Z_FILTER_DEFLATE:
	printf ("H5Z_FILTER_DEFLATE\n");
	break;
      case H5Z_FILTER_SHUFFLE:
	printf ("H5Z_FILTER_SHUFFLE\n");
	break;
      case H5Z_FILTER_FLETCHER32:
	printf ("H5Z_FILTER_FLETCHER32\n");
	break;
      case H5Z_FILTER_SZIP:
	printf ("H5Z_FILTER_SZIP\n");
	break;
      case H5Z_FILTER_NBIT:
	printf ("H5Z_FILTER_NBIT\n");
	break;
      case H5Z_FILTER_SCALEOFFSET:
	printf ("H5Z_FILTER_SCALEOFFSET\n");
      }
    }
#endif
    /*
     * Initialize data.
     */
    rdata = (double *) malloc(sizeof(double)*count[0]);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    /*
     * Read the data using the default properties.
     */
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, rdata);
#if DEBUG
      hsize_t icnt = 0;
      hsize_t j;
      for (i=0; i<count[0]; i++) {
	printf("[ ");
	for (j=0; j<count[1]; j++) {
	  printf("%d ",rdata[icnt]);
	  icnt +=1;
	}
	printf("]\n");
      }
#endif
    /*
     * Find the maximum value in the dataset, to verify that it was
     * read correctly.
     */
    max = rdata[0];
    for (i=1; i < count[0]*count[1]; i++)
      if (max < rdata[i])
	max = rdata[i];
    /*
     * Print the maximum value.
     */
    if(mpi_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      printf ("Maximum value in %s is: %d\n", DATASET, max);
    } else {
      MPI_Reduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    status = H5Pclose (plist_id);
    /*
     * Close and release resources.
     */
    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
    status = H5Fclose (file);
    free(rdata); // free allocated memory
    status = H5Sclose (memspace);
    status = H5Sclose (filespace);
#endif
    MPI_Finalize();
    
    return 0;
}
