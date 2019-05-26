/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_def_var()                  : dispatcher->def_var()
 * ncmpi_inq_varid()                : dispatcher->inq_varid()
 * ncmpi_inq_var()                  : dispatcher->inq_var()
 * ncmpi_rename_var()               : dispatcher->rename_var()
 *
 * ncmpi_get_var<kind>()            : dispatcher->get_var()
 * ncmpi_put_var<kind>()            : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>()     : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>()     : dispatcher->put_var()
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 *
 * ncmpi_iget_var<kind>()           : dispatcher->iget_var()
 * ncmpi_iput_var<kind>()           : dispatcher->iput_var()
 * ncmpi_iget_var<kind>_<type>()    : dispatcher->iget_var()
 * ncmpi_iput_var<kind>_<type>()    : dispatcher->iput_var()
 *
 * ncmpi_buffer_attach()            : dispatcher->buffer_attach()
 * ncmpi_buffer_detach()            : dispatcher->buffer_detach()
 * ncmpi_bput_var<kind>_<type>()    : dispatcher->bput_var()
 *
 * ncmpi_get_varn_<type>()          : dispatcher->get_varn()
 * ncmpi_put_varn_<type>()          : dispatcher->put_varn()
 *
 * ncmpi_iget_varn_<type>()         : dispatcher->iget_varn()
 * ncmpi_iput_varn_<type>()         : dispatcher->iput_varn()
 * ncmpi_bput_varn_<type>()         : dispatcher->bput_varn()
 *
 * ncmpi_get_vard()                 : dispatcher->get_vard()
 * ncmpi_put_vard()                 : dispatcher->put_vard()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>
#include <netcdf_par.h>

#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

int
nc4io_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_def_var */
    err = nc_def_var(nc4p->ncid, name, xtype, ndims, dimids, varidp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Default mode in NetCDF is indep, set to coll if in coll mode */
    if (!(nc4p->flag & NC_MODE_INDEP)){
        err = nc_var_par_access(nc4p->ncid, *varidp, NC_COLLECTIVE);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    if (ndims > nc4p->maxndim){
        nc4p->maxndim = ndims;
    }

    return NC_NOERR;
}

int
nc4io_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err, vid;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_inq_varid */
    err = nc_inq_varid(nc4p->ncid, name, &vid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* NetCDF does not support NULL varid
     * When varid is NULL, NC_NOERR will always return even given invalid name
     */
    if (varid != NULL) *varid = vid;

    return NC_NOERR;
}

int
nc4io_inq_var(void       *ncdp,
              int         varid,
              char       *name,
              nc_type    *xtypep,
              int        *ndimsp,
              int        *dimids,
              int        *nattsp,
              MPI_Offset *offsetp,
              int        *no_fillp,
              void       *fill_valuep)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    if (offsetp != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nc_inq_var(nc4p->ncid, varid, name, xtypep, ndimsp, dimids, nattsp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    if (no_fillp != NULL || fill_valuep != NULL) {
        err = nc_inq_var_fill(nc4p->ncid, varid, no_fillp, fill_valuep);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    return NC_NOERR;
}

int
nc4io_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_rename_var */
    err = nc_rename_var(nc4p->ncid, varid, newname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

/*
nc4io_get_var is implemented in ncmpio_get_put.m4
*/

/*
nc4io_put_var is implemented in ncmpio_get_put.m4
*/

int
nc4io_iget_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    int req_id;
    NC_nc4_req req;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    err = nc4io_init_req(nc4p, &req, varid, start, count, stride, imap, buf, buftype, 0);

    // Add to req list
    nc4io_req_list_add(&(nc4p->getlist), &req_id);
    nc4p->getlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2;
    }

    return NC_NOERR;
}

int
nc4io_iput_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    int req_id;
    NC_nc4_req req;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    err = nc4io_init_req(nc4p, &req, varid, start, count, stride, imap, buf, buftype, 0);

    // Add to req list
    nc4io_req_list_add(&(nc4p->putlist), &req_id);
    nc4p->putlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2 + 1;
    }

    return NC_NOERR;
}

int
nc4io_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    return NC_NOERR;
}

int
nc4io_buffer_detach(void *ncdp)
{
    return NC_NOERR;
}

int
nc4io_bput_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    int req_id;
    NC_nc4_req req;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    err = nc4io_init_req(nc4p, &req, varid, start, count, stride, imap, buf, buftype, 1);

    // Add to req list
    nc4io_req_list_add(&(nc4p->putlist), &req_id);
    nc4p->putlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2 + 1;
    }

    return NC_NOERR;
}

int
nc4io_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int err;
    int i, j;
    int num_all;
    int ndim;
    int esize;
    MPI_Offset putsize;
    MPI_Offset *start = NULL, *count = NULL;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    // ndim
    err = nc_inq_var(nc4p->ncid, varid, NULL, NULL, &ndim, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    // element size
    MPI_Type_size(buftype, &esize);

    // sync num
    MPI_Allreduce(&num, &num_all, 1, MPI_INT, MPI_MAX, nc4p->comm);

    // We need dummy call
    if (num < num_all){
        start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim * 2);
        count = start + ndim;
        memset(start, 0, sizeof(MPI_Offset) * ndim * 2);
    }

    // convert to n vara call
    for(i = 0; i < num_all; i++){
        if (i < num){
            err |= nc4io_get_var(ncdp, varid, starts[i], counts[i], NULL, NULL, buf, bufcount, buftype, reqMode);
            // move buf ptr
            putsize = esize;
            for(j = 0; j < ndim; j++){
                putsize *= counts[i][j];
            }
            buf = (void*)(((char*)buf) + putsize);
        }
        else{
            // Dummy call
            err |= nc4io_get_var(ncdp, varid, start, count, NULL, NULL, buf, bufcount, buftype, reqMode);
        }
    }
    if (err != NC_NOERR) return err;

    // Free buffer
    if (num < num_all){
        NCI_Free(start);
    }

    return NC_NOERR;
}

int
nc4io_put_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               const void        *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int err;
    int i, j;
    int num_all;
    int ndim;
    int esize;
    MPI_Offset putsize;
    MPI_Offset *start = NULL, *count = NULL;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    // ndim
    err = nc_inq_var(nc4p->ncid, varid, NULL, NULL, &ndim, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    // element size
    MPI_Type_size(buftype, &esize);

    // sync num
    MPI_Allreduce(&num, &num_all, 1, MPI_INT, MPI_MAX, nc4p->comm);

    // We need dummy call
    if (num < num_all){
        start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim * 2);
        count = start + ndim;
        memset(start, 0, sizeof(MPI_Offset) * ndim * 2);
    }

    // convert to n vara call
    for(i = 0; i < num_all; i++){
        if (i < num){
            err |= nc4io_put_var(ncdp, varid, starts[i], counts[i], NULL, NULL, buf, bufcount, buftype, reqMode);
            // move buf ptr
            putsize = esize;
            for(j = 0; j < ndim; j++){
                putsize *= counts[i][j];
            }
            buf = (void*)(((char*)buf) + putsize);
        }
        else{
            // Dummy call
            err |= nc4io_put_var(ncdp, varid, start, count, NULL, NULL, buf, bufcount, buftype, reqMode);
        }
    }
    if (err != NC_NOERR) return err;

    // Free buffer
    if (num < num_all){
        NCI_Free(start);
    }

    return NC_NOERR;
}

int
nc4io_iget_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                void               *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    int req_id;
    NC_nc4_req req;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    if (num > 1){
        err = nc4io_init_varn_req(nc4p, &req, varid, num, starts, counts, buf, buftype, 0);
    }
    else{
        err = nc4io_init_req(nc4p, &req, varid, starts[0], counts[0], NULL, NULL, buf, buftype, 0);
    }

    // Add to req list
    nc4io_req_list_add(&(nc4p->getlist), &req_id);
    nc4p->getlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2;
    }

    return NC_NOERR;
}

int
nc4io_iput_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                const void         *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    int req_id;
    NC_nc4_req req;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    if (num > 1){
        err = nc4io_init_varn_req(nc4p, &req, varid, num, starts, counts, buf, buftype, 0);
    }
    else{
        err = nc4io_init_req(nc4p, &req, varid, starts[0], counts[0], NULL, NULL, buf, buftype, 0);
    }

    // Add to req list
    nc4io_req_list_add(&(nc4p->putlist), &req_id);
    nc4p->putlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2 + 1;
    }

    return NC_NOERR;
}

int
nc4io_bput_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                const void         *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    int req_id;
    NC_nc4_req req;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    if (num > 1){
        err = nc4io_init_varn_req(nc4p, &req, varid, num, starts, counts, buf, buftype, 1);
    }
    else{
        err = nc4io_init_req(nc4p, &req, varid, starts[0], counts[0], NULL, NULL, buf, buftype, 1);
    }

    // Add to req list
    nc4io_req_list_add(&(nc4p->putlist), &req_id);
    nc4p->putlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2 + 1;
    }

    return NC_NOERR;
}

int
nc4io_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
}

int
nc4io_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
}

