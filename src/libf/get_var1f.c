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
#define nfmpi_get_var1_ NFMPI_GET_VAR1
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_get_var1_ nfmpi_get_var1__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_get_var1_ nfmpi_get_var1
/* Else leave name alone */
#endif

FORTRAN_API void FORT_CALL nfmpi_get_var1_ ( int *v1, int *v2, int v3[], void*v4, int *v5, MPI_Fint *v6, MPI_Fint *ierr ){
    *ierr = ncmpi_get_var1( *v1, *v2, v3, v4, *v5, *v6 );
}
