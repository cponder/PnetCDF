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
#define nfmpi_put_att_long_ NFMPI_PUT_ATT_LONG
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_put_att_long_ nfmpi_put_att_long__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_put_att_long_ nfmpi_put_att_long
/* Else leave name alone */
#endif

FORTRAN_API void FORT_CALL nfmpi_put_att_long_ ( int *v1, int *v2, char *v3 FORT_MIXED_LEN(d3), nc_type v4, int v5, long*v6, MPI_Fint *ierr FORT_END_LEN(d3) ){
    *ierr = ncmpi_put_att_long( *v1, *v2, v3, v4, v5, v6 );
}
