#! /bin/bash

# This program comes with ABSOLUTELY NO WARRANTY.

# unpack under current directory
UNPACK=`pwd`

# choose names for fast and debug compilation directories
BUILD_DIR="$UNPACK/p4est-build"
BUILD_FAST="$BUILD_DIR/FAST"
BUILD_DEBUG="$BUILD_DIR/DEBUG"

# functions
busage () {
        echo "Usage: `basename $0` <p4est_tar.gz_file> [<install location>]"
}
bdie () {
        echo "Error: $@"
        exit 1
}

if test "x$CFLAGS" = x && test "x$P4EST_CFLAGS_FAST" = x ; then
        export CFLAGS_FAST="-O2"
else
        export CFLAGS_FAST="$CFLAGS $P4EST_CFLAGS_FAST"
fi
echo "CFLAGS_FAST: $CFLAGS_FAST"
if test "x$CFLAGS" = x && test "x$P4EST_CFLAGS_DEBUG" = x ; then
        export CFLAGS_DEBUG="-O0 -g"
else
        export CFLAGS_DEBUG="$CFLAGS $P4EST_CFLAGS_DEBUG"
fi
echo "CFLAGS_DEBUG: $CFLAGS_DEBUG"

TGZ="$1"
if test ! -f "$TGZ" ; then
        busage
        bdie "File not found"
fi
if ! (echo "$TGZ" | grep -q 'p4est.*.tar.gz') ; then
        busage
        bdie "File name mismatch"
fi
shift

# choose names for fast and debug installation directories
INSTALL_DIR="$1"
if test "x$INSTALL_DIR" = x ; then
        INSTALL_DIR="$UNPACK/p4est-install"
else
	shift
fi
INSTALL_FAST="$INSTALL_DIR/FAST"
INSTALL_DEBUG="$INSTALL_DIR/DEBUG"

echo
echo "This script tries to unpack, configure and build the p4est library."
echo "Build FAST: $BUILD_FAST"
echo "Build DEBUG: $BUILD_DEBUG"
echo "Install FAST: $INSTALL_FAST"
echo "Install DEBUG: $INSTALL_DEBUG"
echo "Checking environment: CFLAGS P4EST_CFLAGS_FAST P4EST_CFLAGS_DEBUG"

# remove old versions
if test -d "$BUILD_DIR" ; then
        rm -rf "$BUILD_DIR"
fi

DIR=`echo "$TGZ" | sed 's/\(p4est.*\).tar.gz/\1/'`
DIR=`basename $DIR`
echo "Unpack directory: $UNPACK/$DIR"
if test -d "$UNPACK/$DIR" ; then
        echo "Source directory found (remove it to unpack anew)"
else
        echo -n "Unpacking... "
        tar -xvz -f "$TGZ" -C "$UNPACK" >/dev/null
        echo "done"
fi
test -f "$UNPACK/$DIR/src/p4est.h" || bdie "Main header file missing"
test -f "$UNPACK/$DIR/configure" || bdie "Configure script missing"

echo "See output in files .../config.output and .../make.output"
echo
echo "Build FAST version in $BUILD_FAST"
mkdir -p "$BUILD_FAST"
cd "$BUILD_FAST"
"$UNPACK/$DIR/configure" --enable-mpi --enable-shared \
        --disable-vtk-binary --without-blas --disable-mpithread \
        --prefix="$INSTALL_FAST" CFLAGS="$CFLAGS_FAST" \
        CPPFLAGS="-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL" \
        "$@" > config.output || bdie "Error in configure"
make -C sc -j 8 > make.output || bdie "Error in make sc"
make -j 8 >> make.output || bdie "Error in make p4est"
make install >> make.output || bdie "Error in make install"
echo "FAST version installed in $INSTALL_FAST"

echo
echo "Build DEBUG version in $BUILD_DEBUG"
mkdir -p "$BUILD_DEBUG"
cd "$BUILD_DEBUG"
"$UNPACK/$DIR/configure" --enable-debug --enable-mpi --enable-shared \
        --disable-vtk-binary --without-blas --disable-mpithread \
        --prefix="$INSTALL_DEBUG" CFLAGS="$CFLAGS_DEBUG" \
        CPPFLAGS="-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL" \
        "$@" > config.output || bdie "Error in configure"
make -C sc -j 8 > make.output || bdie "Error in make sc"
make -j 8 >> make.output || bdie "Error in make p4est"
make install >> make.output || bdie "Error in make install"
echo "DEBUG version installed in $INSTALL_DEBUG"
echo
