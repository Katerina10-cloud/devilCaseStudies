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



