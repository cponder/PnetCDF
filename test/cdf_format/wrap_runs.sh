#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

outfile=`basename $1`

for j in 0 1 ; do
    export PNETCDF_SAFE_MODE=$j
    echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc
done

if [ -n "${TESTBB}" ]; then
   echo "---- testing burst buffering"
   for j in 0 1 ; do
       export PNETCDF_SAFE_MODE=$j
       echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
       export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.bb.nc
       unset PNETCDF_HINTS
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.bb.nc

       # running ncmpidiff on dim_cdf12 on one process requires more than 2 GB
       # memory. Use more processes to check. Disable the check for now.
       if test "$1" != "./dim_cdf12" ; then
          echo "--- ncmpidiff $outfile.nc $outfile.bb.nc ---"
          ${TESTSEQRUN} ${NCMPIDIFF} ${TESTOUTDIR}/$outfile.nc ${TESTOUTDIR}/$outfile.bb.nc
       fi
   done
fi
