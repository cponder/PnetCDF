/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

int nc4io_init_req( NC_nc4 *nc4p,
                    NC_nc4_req *req,
                    int        varid,
                    const MPI_Offset *start,
                    const MPI_Offset *count,
                    const MPI_Offset *stride, 
                    const MPI_Offset *imap,
                    const void *buf,
                    MPI_Datatype buftype,
                    int deepcp) {
    int err;
    int i;
    int ndim;
    nc_type xtype;

    // ndim
    err = nc_inq_var(nc4p->ncid, varid, NULL, &xtype, &ndim, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    // Zero out the request
    memset(req, 0, sizeof(NC_nc4_req));

    // Record request
    req->starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*));
    req->start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim);
    req->starts[0] = req->start;
    memcpy(req->start, start, sizeof(MPI_Offset) * ndim);
    req->counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*));
    req->count = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim);
    req->counts[0] = req->count;
    memcpy(req->count, count, sizeof(MPI_Offset) * ndim);
    if (stride != NULL){
        req->stride = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim);
        memcpy(req->stride, stride, sizeof(MPI_Offset) * ndim);
    }
    if (imap != NULL){
        req->imap = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim);
        memcpy(req->imap, imap, sizeof(MPI_Offset) * ndim);
    }

    // bput
    if (deepcp){
        size_t bsize;
        
        err = nc_inq_type(nc4p->ncid, xtype, NULL, &bsize);
        if (err != NC_NOERR){
            return err;
        }
        
        for(i = 0; i < ndim; i++){
            bsize *= count[i];
        }

        req->buf = (char*)NCI_Malloc(bsize);

        memcpy(req->buf, buf, bsize);
    }
    else{
        req->buf = (void*)buf;
    }

    req->varid = varid;
    req->bufs = (char**)NCI_Malloc(sizeof(char*));
    req->bufs[0] = req->buf;
    req->nreq = 1;
    req->deepcp = deepcp;
    req->buftype = buftype;

    return NC_NOERR;
}

int nc4io_init_varn_req( NC_nc4 *nc4p,
                        NC_nc4_req *req,
                        int        varid,
                        int        nreq,
                        MPI_Offset *const*starts,
                        MPI_Offset *const*counts, 
                        const void *buf,
                        MPI_Datatype buftype,
                        int deepcp) {

    int err;
    int i, j;
    MPI_Offset rsize, boff;
    int ndim;
    size_t esize;
    nc_type xtype;

    // ndim
    err = nc_inq_var(nc4p->ncid, varid, NULL, &xtype, &ndim, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    // element size
    err = nc_inq_type(nc4p->ncid, xtype, NULL, &esize);
    if (err != NC_NOERR){
        return err;
    }

    // Zero out the request
    memset(req, 0, sizeof(NC_nc4_req));

    // Record request
    req->starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nreq);
    req->start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim * nreq);
    for(i = 0; i < nreq; i++){
        req->starts[i] = req->start + i * ndim;
        memcpy(req->starts[i], starts[i], sizeof(MPI_Offset) * ndim);
    }
    req->counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nreq);
    req->count = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndim * nreq);
    for(i = 0; i < nreq; i++){
        req->counts[i] = req->count + i * ndim;
        memcpy(req->counts[i], counts[i], sizeof(MPI_Offset) * ndim);
    }

    // bput
    if (deepcp){      
        // count buffer size
        rsize = 0;
        for(i = 0; i < nreq; i++){
            boff = esize;
            for(j = 0; j < ndim; j++){
                boff *= counts[i][j];
            }
            rsize += boff;
        }

        req->buf = (char*)NCI_Malloc(rsize);

        memcpy(req->buf, buf, rsize);
    }
    else{
        req->buf = (void*)buf;
    }

    // Calculate buffer for each individual request
    req->bufs = (char**)NCI_Malloc(sizeof(char*) * nreq);
    boff = 0;
    for(i = 0; i < nreq; i++){
        req->bufs[i] = (((char*)req->buf) + boff);

        // Advance pointer by size of the request
        rsize = esize;
        for(j = 0; j < ndim; j++){
            rsize *= counts[i][j];
        }
        boff += rsize;
    }

    req->varid = varid;
    req->buf = (void*)buf;
    req->nreq = nreq;
    req->deepcp = deepcp;
    req->buftype = buftype;

    return NC_NOERR;
}
