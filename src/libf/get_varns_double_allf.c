/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface -infile=../lib/pnetcdf.h -deffile=defs
 * DO NOT EDIT
 */
#include "mpinetcdf_impl.h"


#ifdef F77_NAME_UPPER
#define nfmpi_get_varns_double_all_ NFMPI_GET_VARNS_DOUBLE_ALL
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_get_varns_double_all_ nfmpi_get_varns_double_all__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_get_varns_double_all_ nfmpi_get_varns_double_all
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_get_varns_double_all_ ( int *v1, int *v2, int *v3, MPI_Offset*v4, MPI_Offset*v5, MPI_Offset*v6, double*v7 ){
    int ierr;
    int l2 = *v2 - 1;
    MPI_Offset **l4 = 0;
    MPI_Offset **l5 = 0;
    MPI_Offset **l6 = 0;

    { int ndims = ncmpixVardim(*v1,*v2-1); /* var's number dims */
    if (ndims > 0) {
        int li, lj;
        l4    = (MPI_Offset**)malloc(*v3 * sizeof(MPI_Offset*) );
        l4[0] = (MPI_Offset*) malloc(*v3 * ndims * sizeof(MPI_Offset));
        for (lj=1; lj<*v3; lj++) 
            l4[lj] = l4[lj-1] + ndims;
        for (lj=0; lj<*v3; lj++) 
            for (li=0; li<ndims; li++) 
                l4[lj][li] = v4[lj*ndims + ndims-1-li] - 1;
    }
    else if (ndims < 0) {
        /* Error return */
        ierr = ndims; 
	return ierr;
    }
    }

    { int ndims = ncmpixVardim(*v1,*v2-1); /* var's number dims */
    if (ndims > 0) {
        int li, lj;
        l5    = (MPI_Offset**)malloc(*v3 * sizeof(MPI_Offset*) );
        l5[0] = (MPI_Offset*) malloc(*v3 * ndims * sizeof(MPI_Offset));
        for (lj=1; lj<*v3; lj++) 
            l5[lj] = l5[lj-1] + ndims;
        for (lj=0; lj<*v3; lj++) 
            for (li=0; li<ndims; li++) 
                l5[lj][li] = v5[lj*ndims + ndims-1-li];
    }
    else if (ndims < 0) {
        /* Error return */
        ierr = ndims; 
	return ierr;
    }
    }

    { int ndims = ncmpixVardim(*v1,*v2-1); /* var's number dims */
    if (ndims > 0) {
        int li, lj;
        l6    = (MPI_Offset**)malloc(*v3 * sizeof(MPI_Offset*) );
        l6[0] = (MPI_Offset*) malloc(*v3 * ndims * sizeof(MPI_Offset));
        for (lj=1; lj<*v3; lj++) 
            l6[lj] = l6[lj-1] + ndims;
        for (lj=0; lj<*v3; lj++) 
            for (li=0; li<ndims; li++) 
                l6[lj][li] = v6[lj*ndims + ndims-1-li];
    }
    else if (ndims < 0) {
        /* Error return */
        ierr = ndims; 
	return ierr;
    }
    }
    ierr = ncmpi_get_varns_double_all( *v1, l2, *v3, l4, l5, l6, v7 );

    if (l4) { free(l4[0]); free(l4); }

    if (l5) { free(l5[0]); free(l5); }

    if (l6) { free(l6[0]); free(l6); }
    return ierr;
}
