#!/usr/bin/bash 

SOURCE_PATH="/u/area/ntosato/scratch/timing/time_scale_final/software/source"
INSTALL_PATH="/u/area/ntosato/scratch/timing/time_scale_final/software/programs/R"
MODULE_PATH="/u/area/ntosato/scratch/timing/time_scale_final/software/modules"
VERSION="4.3.3"
PLATFORM="h100"
mkdir -p "$INSTALL_PATH/$VERSION-$PLATFORM"

###clone

cd $SOURCE_PATH/
 rm -rf R-4.3.3*
wget https://cloud.r-project.org/src/base/R-4/R-${VERSION}.tar.gz
tar -xzf R-${VERSION}.tar.gz

### use the modules 

#source /etc/profile.d/lmod.sh

module use $MODULE_PATH
module load openBLAS/0.3.29-$PLATFORM

echo $OPENBLAS_LIB

###compile
cd $SOURCE_PATH/R-$VERSION
./configure --prefix=${INSTALL_PATH}/${VERSION}-${PLATFORM} --with-x=no \
                          	    --with-blas="-L${OPENBLAS_LIB} -lopenblas" \
                          	    --with-lapack=yes \
                          	    --with-system-valgrind-headers \
                          	    --enable-memory-profiling \
                          	    --with-valgrind-instrumentation=2 \
                          	    --with-system-valgrind-headers \
                          	    --with-jpeglib \
                          	    --with-libpng \
                          	    --with-tcltk \
                          	    --with-readline \
                          	    --with-cairo=yes \
                          	    --enable-R-profiling  \
                          	    --with-libtiff \
                          	    --enable-lto \
                                --enable-R-shlib
make -j $(nproc --all)
make install


