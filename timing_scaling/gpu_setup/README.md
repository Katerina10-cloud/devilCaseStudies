# Devil setup

Before proceeding, ensure that: 

- Your system has CUDA installed and configured properly (CUDA > 12.0 recommended).
- `Lmod` is installed for module management.

The module system is necessary in a shared enviroment such a Cluster, since you are not supposed to have root access.


We provide a sample setup with installation scripts and module file.

# OpenBLAS setup

Since by default R can use the internal BLAS as fallback if no system library are present, we provide a compile script for openBLAS and the module file to consume them, suitable to use in a cluster environment.

```bash
#!/usr/bin/bash -e

SOURCE_PATH="software/source"
INSTALL_PATH="software/programs/openBLAS"
MODULE_PATH="software/modules/openBLAS"

VERSION="0.3.29"
PLATFORM="h100"

mkdir "$INSTALL_PATH/$VERSION-$PLATFORM"

###clone

cd $SOURCE_PATH
rm -rf OpenBLAS/
git clone --depth 1 --branch v$VERSION  git@github.com:OpenMathLib/OpenBLAS.git

###compile

cd $SOURCE_PATH/OpenBLAS
mkdir -p  ${INSTALL_PATH}/$VERSION-$PLATFORM/
make clean
make distclean
make USE_OPENMP=0 -j $(nproc --all)
make PREFIX=${INSTALL_PATH}/$VERSION-$PLATFORM/ install
```

The setup that we adopted is to use pthread based libraries instead of openMP based, since our program use openMP, we wwant to avoid any kind of issue. 

After that, is possible to build a module file to load openBLAS libraries:

```lua
-- -*- lua -*-

local name      = "openBLAS"
local version   = "0.3.29-h100"

whatis("Name         : " .. name)
whatis("Version      : " .. version)

family("BLAS")

local home    = "software/programs/openBLAS/0.3.29-h100/"

prepend_path{"PATH", home .. "/bin",delim=":",priority="0"}
prepend_path{"LD_LIBRARY_PATH", home .. "/lib",delim=":",priority="0"}
prepend_path{"LIBRARY_PATH", home .. "/lib",delim=":",priority="0"}
prepend_path{"CPATH", home .. "/include",delim=":",priority="0"}
setenv("OPENBLAS_DIR", home )
setenv("OPENBLAS_ROOT", home)
setenv("OPENBLAS_LIB", home .. "/lib")
setenv("OPENBLAS_IN",home .. "/include")
```


# R setup

Now to exploit the maximum performance we compile the R software using the previous module as dependency.

```bash
#!/usr/bin/bash

SOURCE_PATH="software/source"
INSTALL_PATH="software/programs/R"
MODULE_PATH="software/modules"
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
```

This allow  us ot obtain a software optimized for the underling architecture in order to achieve the maximum performance from the hardware. 

```lua
-- -*- lua -*-

local name      = "R"
local version   = "4.3.3-h100"
whatis("Name         : " .. name)
whatis("Version      : " .. version)

family("R")
depends_on("0.3.29-h100")

--conflict()

local home    = "software/programs/R/4.3.3-h100"


prepend_path("PATH", home .. "/bin")
prepend_path("LD_LIBRARY_PATH", home .."/lib64")
prepend_path("MANPATH", home .."/share/man")
prepend_path("R_LIBS_USER","software/programs/r_packages_h100")
```

Is foundamental to specify a custom `R_LIBS_USER` in order to avoid installation clash between architecture.

# Cutensor setup

Cutensor is distributed by nvidia. Follow a script to download, extract and *install* it.

```bash
#!/usr/bin/bash

SOURCE_PATH="software/source"
INSTALL_PATH="software/programs/cutensor"
MODULE_PATH="software/modules"
VERSION="2.2.0.0"

###clone

cd $SOURCE_PATH/
rm libcutensor*

wget https://developer.download.nvidia.com/compute/cutensor/redist/libcutensor/linux-x86_64/libcutensor-linux-x86_64-$VERSION-archive.tar.xz

tar -xf libcutensor-linux-x86_64-$VERSION-archive.tar.xz
cd libcutensor-linux-x86_64-$VERSION-archive/
mkdir -p $INSTALL_PATH/cutensor/$VERSION/
cp -r ./* $INSTALL_PATH/cutensor/$VERSION/
```

This will install several version of cutensor, compatible with multiple cuda version, we will use the last one, so > 12.0

The relative module will be:

```lua
-- -*- lua -*-

local name      = "cutensor"
local version   = "2.2.0.0"
whatis("Name         : " .. name)
whatis("Version      : " .. version)

family("cutensor")
depends_on("cuda")

--conflict()

--In this section all dependencies will be specified

--binary location without architecture extension
local home    = "software/programs/cutensor/2.2.0.0"


prepend_path("CPATH", home .. "/include")
prepend_path("INCLUDE", home .. "/include")
prepend_path("LD_LIBRARY_PATH", home .. "/lib/12/")
prepend_path("LIBRARY_PATH", home .. "/lib/12/")
prepend_path("CUTENSOR_HOME", home)
```

# Note on cuda enviroment and modules

This procedure assume that you have to change the relative path on the modules to absolute path, moreover you need a working cuda installation. The environmental variable foundamental for the installation of the GPU code are `CUDA_HOME` and `CUTENSOR_HOME`.

# Devil setup

In order to install devil, you need to instruct the module system (`Lmod` in our setup) to use the new stack software.
```bash
$ module use software/modules/
$ module load R/4.3.3-h100 openBLAS/0.3.29-h100 cutensor/2.2.0.0
```
Check that the environmental variables are correctly setted (sample path):

```bash
$ echo $OPENBLAS_LIB
/u/area/ntosato/scratch/timing/time_scale_final/software/programs/openBLAS/0.3.29-h100//lib
$ echo $CUDA_HOME
/opt/programs/cuda/12.1
$ echo $CUTENSOR_HOME
/u/area/ntosato/scratch/timing/time_scale_final/software/programs/cutensor/2.2.0.0
``

And then install `devil`, using the branch `devel` that contain the accellerated code:

```
devtools::install_github("caravagnalab/devil@devel") 
```




