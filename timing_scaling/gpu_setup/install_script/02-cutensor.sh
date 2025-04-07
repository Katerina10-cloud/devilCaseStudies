#!/usr/bin/bash 

SOURCE_PATH="/u/area/ntosato/scratch/timing/time_scale_final/software/source"
INSTALL_PATH="/u/area/ntosato/scratch/timing/time_scale_final/software/programs/R"
MODULE_PATH="/u/area/ntosato/scratch/timing/time_scale_final/software/modules"
VERSION="2.2.0.0"

###clone

cd $SOURCE_PATH/
rm libcutensor*

wget https://developer.download.nvidia.com/compute/cutensor/redist/libcutensor/linux-x86_64/libcutensor-linux-x86_64-$VERSION-archive.tar.xz

tar -xf libcutensor-linux-x86_64-$VERSION-archive.tar.xz
cd libcutensor-linux-x86_64-$VERSION-archive/
mkdir -p $INSTALL_PATH/cutensor/$VERSION/
cp -r ./* $INSTALL_PATH/cutensor/$VERSION/
