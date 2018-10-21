/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_misc()         : dispatcher->inq_misc()
 * ncmpi_wait()             : dispatcher->wait()
 * ncmpi_wait_all()         : dispatcher->wait()
 * ncmpi_cancel()           : dispatcher->cancel()
 *
 * ncmpi_set_fill()         : dispatcher->set_fill()
 * ncmpi_fill_var_rec()     : dispatcher->fill_rec()
 * ncmpi_def_var_fill()     : dispatcher->def_var_fill()
 * ncmpi_inq_var_fill()     : dispatcher->inq()
 *
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_flush()             : dispatcher->flush()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>

int ncadios_sync_header(NC_ad *ncadp) {
    int i;
    int bsize;
    int namelen;
    char *buf, *cur;

    if (ncadp->rank == 0){
        bsize = sizeof(int) * 2;   // natt and nvar
        for(i = 0; i < ncadp->dims.cnt; i++){
            bsize += strlen(ncadp->dims.data[i].name) + 1 + sizeof(int) * 2;
        }
        for(i = 0; i < ncadp->vars.cnt; i++){
            bsize += strlen(ncadp->vars.data[i].name) + 1 + sizeof(int) * 2 + ncadp->vars.data[i].ndim * sizeof(int) + sizeof(nc_type);
        }
    }
    
    MPI_Bcast(&bsize, 1, MPI_INT, 0, ncadp->comm);

    buf = NCI_Malloc(bsize);
    cur = buf;

    if (ncadp->rank == 0){
        *((int*)cur) = ncadp->dims.cnt;
        cur += sizeof(int);
        for(i = 0; i < ncadp->dims.cnt; i++){
            *((int*)cur) = ncadp->dims.data[i].len;
            cur += sizeof(int);
            namelen = strlen(ncadp->dims.data[i].name);
            *((int*)cur) = namelen;
            cur += sizeof(int);
            strcpy(cur, ncadp->dims.data[i].name);
            cur += namelen + 1;
        }
        *((int*)cur) = ncadp->vars.cnt;
        cur += sizeof(int);
        for(i = 0; i < ncadp->vars.cnt; i++){
            *((nc_type*)cur) = ncadp->vars.data[i].type;
            cur += sizeof(nc_type);
            *((int*)cur) = ncadp->vars.data[i].ndim;
            cur += sizeof(int);
            namelen = strlen(ncadp->vars.data[i].name);
            *((int*)cur) = namelen;
            cur += sizeof(int);
            memcpy(cur, ncadp->vars.data[i].dimids, ncadp->vars.data[i].ndim * sizeof(int));
            cur += ncadp->vars.data[i].ndim * 4;
            strcpy(cur, ncadp->vars.data[i].name);
            cur += namelen + 1;         
        }
    }

    MPI_Bcast(buf, bsize, MPI_BYTE, 0, ncadp->comm);


    if (ncadp->rank != 0){
        int ndim, nvar;
        int id, len;
        nc_type type;
        int *dimids;
        char *name;

        ndim = *((int*)cur);
        cur += 4;
        //printf("ndim = %d\n", ndim);
        for(i = 0; i < ndim; i++){
            len = *((int*)cur);
            cur += sizeof(int);
            namelen = *((int*)cur);
            cur += sizeof(int);
            name = cur;
            cur += namelen + 1;
            ncadiosi_def_dim(ncadp, name, len, &id);
            //printf("def_dim(%s, %d), namelen = %d\n", name, len, namelen);
        }

        
        nvar = *((int*)cur);
        cur += 4;
        for(i = 0; i < nvar; i++){
            type = *((nc_type*)cur);
            cur += sizeof(nc_type);
            ndim = *((int*)cur);
            cur += sizeof(int);
            namelen = *((int*)cur);
            cur += sizeof(int);
            dimids = (int*)cur;
            cur += ndim * sizeof(int);
            name = cur;
            cur += namelen + 1;   
            ncadiosi_def_var(ncadp, name, type, ndim, dimids, &id);
            //printf("def_var(%s, %d, %d), namelen = %d\n", name, type, ndim, namelen);
        }
        
    }

    NCI_Free(buf);
}