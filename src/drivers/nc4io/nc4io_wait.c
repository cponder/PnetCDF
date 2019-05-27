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
#include <nc4io_internal.h>

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nc4io_wait_put_reqs(NC_nc4 *nc4p, int nreq, int *reqids, int *stats){
    int err;
    int i, j;
    int nvar;
    int *nums, *nums_all;  // Number of reqs in each varn
    int *offs;  // Number of reqs in each varn
    int *idmap;
    MPI_Offset **starts, **counts, **strides;
    MPI_Offset *start_dummy, *count_dummy, *stride_dummy;
    MPI_Datatype *btypes;
    char **bufs;
    NC_nc4_req *req;

    // Create dummy stride and count
    start_dummy = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * nc4p->maxndim * 3);
    count_dummy = start_dummy + nc4p->maxndim;
    stride_dummy = count_dummy + nc4p->maxndim;
    memset(start_dummy, 0, sizeof(MPI_Offset) * nc4p->maxndim);
    memset(count_dummy, 0, sizeof(MPI_Offset) * nc4p->maxndim);
    for(i = 0; i < nc4p->maxndim; i++){
        stride_dummy[i] = 1;
    }

    /* Call nc_inq */
    err = nc_inq(nc4p->ncid, NULL, &nvar, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    // Count total number of request per variable 
    nums = (int*)NCI_Malloc(sizeof(int) * nvar * 2);
    nums_all = nums + nvar;
    offs = (int*)NCI_Malloc(sizeof(int) * (nvar + 1));
    memset(nums, 0, sizeof(int) * nvar);
    for(i = 0; i < nreq; i++){
        req = nc4p->putlist.reqs + reqids[i];
        nums[req->varid] += req->nreq;
    }

    // Start position of each var
    offs[0] = 0;
    for(i = 0; i < nvar; i++){
        offs[i + 1] = offs[i] + nums[i];
    }

    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * offs[nvar] * 3);
    counts = starts + offs[nvar];
    strides = counts + offs[nvar];
    bufs = (char**)NCI_Malloc(sizeof(char*) * offs[nvar]);
    btypes = (MPI_Datatype*)NCI_Malloc(sizeof(MPI_Datatype) * offs[nvar]);
    if (stats != NULL){
        idmap = (int*)NCI_Malloc(sizeof(int) * offs[nvar]);
        memset(stats, 0, sizeof(int) * nreq);
    }
    // Convert all call to vars call since netcdf does not support different type of API call
    memset(nums, 0, sizeof(int) * nvar);
    for(i = 0; i < nreq; i++){
        req = nc4p->putlist.reqs + reqids[i];
        memcpy(starts + offs[req->varid] + nums[req->varid], req->starts, sizeof(MPI_Offset*) * req->nreq);
        memcpy(counts + offs[req->varid] + nums[req->varid], req->counts, sizeof(MPI_Offset*) * req->nreq);
        memcpy(bufs + offs[req->varid] + nums[req->varid], req->bufs, sizeof(char*) * req->nreq);
        if (req->stride != NULL){
            strides[offs[req->varid] + nums[req->varid]] = req->stride;
        }
        else{
            strides[offs[req->varid] + nums[req->varid]] = stride_dummy;
        }
        for(j = 0; j < req->nreq; j++){
            btypes[offs[req->varid] + nums[req->varid] + j] = req->buftype;
        }
        if (stats != NULL){
            for(j = 0; j < req->nreq; j++){
                idmap[offs[req->varid] + nums[req->varid] + j] = i;
            }
        }
        nums[req->varid] += req->nreq;
    }

    // sync nums
    MPI_Allreduce(nums, nums_all, nvar, MPI_INT, MPI_MAX, nc4p->comm);

    // Call NetCDF APIs
    for(i = 0; i < nvar; i++){
        for(j = 0; j < nums_all[i]; j++){
            if (j < nums[i]){
                err = nc4io_put_var(nc4p, i, starts[offs[i] + j], counts[offs[i] + j], strides[offs[i] + j], NULL, bufs[offs[i] + j], -1, btypes[offs[i] + j], NC_REQ_COLL);
            }
            else{
                nc4io_get_var(nc4p, i, start_dummy, count_dummy, stride_dummy, NULL, &err, -1, MPI_INT, NC_REQ_COLL);
            }
            if (stats != NULL){
                stats[idmap[offs[i] + j]] |= err;
            }
        }
    }

    NCI_Free(start_dummy);
    NCI_Free(nums);
    NCI_Free(offs);
    NCI_Free(starts);
    NCI_Free(bufs);
    NCI_Free(btypes);
    if (stats != NULL){
        NCI_Free(idmap);
    }

    return NC_NOERR;
}

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nc4io_wait_get_reqs(NC_nc4 *nc4p, int nreq, int *reqids, int *stats){
    int err;
    int i, j;
    int nvar;
    int *nums, *nums_all;  // Number of reqs in each varn
    int *offs;  // Number of reqs in each varn
    int *idmap;
    MPI_Offset **starts, **counts, **strides;
    MPI_Offset *start_dummy, *count_dummy, *stride_dummy;
    MPI_Datatype *btypes;
    char **bufs;
    NC_nc4_req *req;

    // Create dummy stride and count
    start_dummy = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * nc4p->maxndim * 3);
    count_dummy = start_dummy + nc4p->maxndim;
    stride_dummy = count_dummy + nc4p->maxndim;
    memset(start_dummy, 0, sizeof(MPI_Offset) * nc4p->maxndim);
    memset(count_dummy, 0, sizeof(MPI_Offset) * nc4p->maxndim);
    for(i = 0; i < nc4p->maxndim; i++){
        stride_dummy[i] = 1;
    }

    /* Call nc_inq */
    err = nc_inq(nc4p->ncid, NULL, &nvar, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    // Count total number of request per variable 
    nums = (int*)NCI_Malloc(sizeof(int) * nvar * 2);
    nums_all = nums + nvar;
    offs = (int*)NCI_Malloc(sizeof(int) * (nvar + 1));
    memset(nums, 0, sizeof(int) * nvar);
    for(i = 0; i < nreq; i++){
        req = nc4p->putlist.reqs + reqids[i];
        nums[req->varid] += req->nreq;
    }

    // Start position of each var
    offs[0] = 0;
    for(i = 0; i < nvar; i++){
        offs[i + 1] = offs[i] + nums[i];
    }

    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * offs[nvar] * 3);
    counts = starts + offs[nvar];
    strides = counts + offs[nvar];
    bufs = (char**)NCI_Malloc(sizeof(char*) * offs[nvar]);
    btypes = (MPI_Datatype*)NCI_Malloc(sizeof(MPI_Datatype) * offs[nvar]);
    if (stats != NULL){
        idmap = (int*)NCI_Malloc(sizeof(int) * offs[nvar]);
        memset(stats, 0, sizeof(int) * nreq);
    }
    // Convert all call to vars call since netcdf does not support different type of API call
    memset(nums, 0, sizeof(int) * nvar);
    for(i = 0; i < nreq; i++){
        req = nc4p->putlist.reqs + reqids[i];
        memcpy(starts + offs[req->varid] + nums[req->varid], req->starts, sizeof(MPI_Offset*) * req->nreq);
        memcpy(counts + offs[req->varid] + nums[req->varid], req->counts, sizeof(MPI_Offset*) * req->nreq);
        memcpy(bufs + offs[req->varid] + nums[req->varid], req->bufs, sizeof(char*) * req->nreq);
        if (req->stride != NULL){
            strides[offs[req->varid] + nums[req->varid]] = req->stride;
        }
        else{
            strides[offs[req->varid] + nums[req->varid]] = stride_dummy;
        }
        for(j = 0; j < req->nreq; j++){
            btypes[offs[req->varid] + nums[req->varid] + j] = req->buftype;
        }
        if (stats != NULL){
            for(j = 0; j < req->nreq; j++){
                idmap[offs[req->varid] + nums[req->varid] + j] = i;
            }
        }
        nums[req->varid] += req->nreq;
    }

    // sync nums
    MPI_Allreduce(nums, nums_all, nvar, MPI_INT, MPI_MAX, nc4p->comm);

    // Call NetCDF APIs
    for(i = 0; i < nvar; i++){
        for(j = 0; j < nums_all[i]; j++){
            if (j < nums[i]){
                err = nc4io_put_var(nc4p, i, starts[offs[i] + j], counts[offs[i] + j], strides[offs[i] + j], NULL, bufs[offs[i] + j], -1, btypes[offs[i] + j], NC_REQ_COLL);
            }
            else{
                nc4io_get_var(nc4p, i, start_dummy, count_dummy, stride_dummy, NULL, &err, -1, MPI_INT, NC_REQ_COLL);
            }
            if (stats != NULL){
                stats[idmap[offs[i] + j]] |= err;
            }
        }
    }

    NCI_Free(start_dummy);
    NCI_Free(nums);
    NCI_Free(offs);
    NCI_Free(starts);
    NCI_Free(bufs);
    NCI_Free(btypes);
    if (stats != NULL){
        NCI_Free(idmap);
    }

    return NC_NOERR;
}