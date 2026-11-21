#!/bin/sh
set -eu
cd "$(dirname "$0")"
export R_LIBS_USER="$PWD/environment/library"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
exec Rscript --vanilla "$@"
