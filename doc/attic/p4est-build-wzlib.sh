#! /bin/bash

# Build and install the latest p4est develop branch including zlib.
# We first download and install a recent zlib to a local directory,
# then clone and install the current develop branch of libsc using
# that zlib, and then clone and install the current develop branch
# of p4est using that libsc and zlib.

# This results in three installation directories that any higher
# level software package may be compiled and linked against.
# The options are similar to those used in this script.

# set installation root to local subdirectory
PREFIX="$PWD/local"

# feel free to make changes to the libsc and p4est configure line.
CONFIG="--enable-mpi"

# download, build and install zlib
wget -N https://www.zlib.net/zlib-1.3.tar.gz
tar -xvzf zlib-1.3.tar.gz
cd zlib-1.3
./configure --prefix="$PREFIX/zlib"
make -j install
cd ..
rm -r zlib-1.3 zlib-1.3.tar.gz

# provide environment that links to installed zlib
export CPPFLAGS="-I$PREFIX/zlib/include"
export CFLAGS="-O2 -g -Wall -Wl,-rpath=$PREFIX/zlib/lib"
export LDFLAGS="-L$PREFIX/zlib/lib"

# clone, build and install libsc
git clone https://github.com/cburstedde/libsc.git -b develop
cd libsc
./bootstrap
mkdir build
cd build
../configure $CONFIG --prefix="$PREFIX/libsc"
make -j install V=0
cd ../../
rm -rf libsc/.git
rm -r libsc

# clone, build and install p4est
git clone https://github.com/cburstedde/p4est.git -b develop
cd p4est
./bootstrap "$PREFIX/libsc/share/aclocal"
mkdir build
cd build
../configure $CONFIG --with-sc="$PREFIX/libsc" --prefix="$PREFIX/p4est"
make -j install V=0
cd ../../
rm -rf p4est/.git
rm -r p4est
