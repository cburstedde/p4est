#! /bin/bash

# Build and install the latest p4est develop branch including zlib and
# jansson.  We first download and install a recent zlib to a local
# directory, then do the same for the jansson library.  Then we clone and
# install the current develop branch of libsc using that zlib and jansson
# installation and finally clone and install the current develop branch of
# p4est, linking against all of the above libraries.

# This results in four installation directories that any higher
# level software package may be compiled and linked against.
# The options are similar to those used in this script.
# In particular, the -rpath option may turn out useful.

# let the libsc installation script set the stage
BLSC="libsc-build-wdeps.sh"
wget -N "https://github.com/cburstedde/libsc/raw/develop/doc/$BLSC"     && \
source "$BLSC"                                          && \
rm "$BLSC"                                              || exit 1

# clone, build and install p4est
git clone --depth 1 https://github.com/cburstedde/p4est.git -b develop  && \
cd p4est                                                && \
./bootstrap "$PREFIX/libsc/share/aclocal"               && \
mkdir build                                             && \
cd build                                                && \
../configure $CONFIG --with-sc="$PREFIX/libsc" --prefix="$PREFIX/p4est" && \
make -j install V=0                                     && \
cd ../../                                               && \
rm -rf p4est/.git                                       && \
rm -r p4est                                             || bdie "p4est"
