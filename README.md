# The p4est Library

This is the README file for `p4est`.

`p4est` is a C library to manage a collection (a forest) of multiple
connected adaptive quadtrees or octrees in parallel.

    Copyright (C) 2010 The University of Texas System
    Additional copyright (C) 2011 individual authors

`p4est` is written by Carsten Burstedde, Lucas C. Wilcox, Tobin Isaac
and others and released under the GNU General Public Licence version 2.

The official web page for source code and documentation is
[p4est.org](https://www.p4est.org/).
Presently, we maintain the code on
[github](https://www.github.com/cburstedde/p4est.git).
Please feel free to participate.

## Acknowledgement and Disclaimer

The development of `p4est` was partially supported by the US National Science
Foundation (NSF Grants No. CCF-0427985, CMMI-1028889, CNS-0540372, CNS-0619838,
DMS-0724746, OCI-0749334, OPP-0941678) and the US Department of Energy (DOE
Grants No. 06ER25782, 08ER25860, SC0002710).

The authors thank the Texas Advanced
Computing Center (TACC) for providing them with access to the Ranger
supercomputer under NSF TeraGrid award MCA04N026, and the National Center for
Computational Science (NCCS) for early-user access to the Jaguar Cray XT5
supercomputer.  Any opinions, findings and conclusions or recommendations
expressed in the source code and documentation are those of the authors and do
not necessarily reflect the views of the National Science Foundation (NSF).

The authors gratefully acknowledge partial support under DARPA Cooperative
Agreement HR00112120003 via a subcontract with Embry-Riddle Aeronautical
University.  This work is approved for public release; distribution is
unlimited.  The information in this repository does not necessarily reflect
the position or the policy of the US Government.

## Purpose

Many applications in applied mathematics and numerical simulation use a mesh of
computational cells that covers the domain of interest.  The numerical solution
is then approximated by functions, each of which is associated with a small set
of cells (or even one).  Dynamic "adaptive" methods change the mesh during
the simulation by local refinement or coarsening, and "parallel" methods
distribute ("partition") the mesh between multiple processors, where each
processor ideally receives an equal share of the computational load.  `p4est`
isolates the task of parallel dynamic adaptive mesh refinement (AMR),
coarsening, and load balancing, and encapsulates algorithms that scale well to
large numbers (10^6) of MPI processes.

These algorithms are predominantly used for the numerical solution of partial
differential equations, but also support various other tasks where fast
hierarchical space subdivision is required, for example to locate particles in
space or to organize, sort, and search in large data sets.

## Geometric structure

The basic structure used in `p4est` is a "connectivity" of quadtrees (2D) or
octrees (3D) that covers the domain of interest in a conforming macro-mesh.
This includes the case of using exactly one tree for representing the
hypercube.  The trees can be arbitrarily refined and coarsened, where the
storage of quadrants/octants is distributed in parallel.  Thus, the mesh
primitives used in `p4est` are quadrilaterals in 2D and hexahedra in 3D.  The
adaptive structure allows for quadrants/octants of different sizes to neighbor
each other, which is commonly called "non-conforming".  This concept leads to
"hanging" faces and edges.

## Core `p4est` (2D) and `p8est` (3D) routines

`p?est_new`: Create an equi-partitioned, uniformly refined forest.

`p?est_refine`: Adaptively subdivide octants based on a callback function, once
for each octant or recursively.

`p?est_coarsen`: Replace families of eight child octants by their common parent
octant, once or recursively.

`p?est_partition`: Redistribute the octants in parallel, according to a given
target number of octants for each process, or weights prescribed for all
octants.

`p?est_balance`: Ensure at most 2:1 size relations between neighboring octants by
local refinement where necessary.

`p?est_checksum`: Compute a partition-independent integer "fingerprint" of a
forest.  This is useful for debugging and regression testing.
Please note that the checksum is not cryptographically secure;
we do encountere false positive matches.

## Interfacing to `p4est` from applications

`p?est_ghost`: Collect one layer of off-process octants touching the process
boundaries from the outside.  This function does *not* require a previous
call to `p?est_balance`.  This is the most generally useful function for
external applications.  By querying the ghost layer, the application can
associate degrees of freedom with the mesh which are the basis for all
numerical computation.

`p?est_lnodes`: Create a globally unique numbering of finite element degrees of
freedom for arbitrary order C0 nodal polynomials, also known as GLL basis.  The
function takes into account the classification into "independent" and
"hanging" degrees of freedom.  Numberings for standard piecewise d-linear
finite elements or lowest-order Raviart-Thomas velocity variables are created
by calling this function with the parameters degree = 1 and -1, respectively.
This function requires previous calls to `p?est_balance` and `p?est_ghost`.

`p?est_nodes`: Like `p?est_lnodes`, but only for piecewise linear finite elements.
Technically no longer required, but we keep it around for its simplicity.

## Installation from a release tarball

Our official tarballs are linked from (p4est.org)[http://www.p4est.org/] and work fine.
The unnumbered archives under [releases](https://github.com/cburstedde/p4est/releases/)
are created by github and are *not* endorsed by us.  In particular, they are
missing the subdirectory sc, the configure script, and other generated files.

`p4est` uses the GNU autoconf/automake/libtool build system, but also supports
CMake.  See the INSTALL file and the CMake section below for details.
Our official tarballs contain the configure script and all
required files.  The following configure options are commonly used.

`--enable-mpi` is necessary for parallel functionality.  Without this option,
`p4est` does not require MPI, which is replaced by a dummy implementation that
always assumes an `MPI_Comm_size` of 1 and an `MPI_Comm_rank` of 0.

`--enable-debug` enables assertions and additional code for verification.  This
is generally helpful for development, even if it is somewhat slower and
produces a lot of diagnostic log messages.

`CFLAGS="-Wall -g -O0"` is an example of setting flags for a development build
compiled with the GNU compilers.  `p4est` is written in plain C and does not
require `FFLAGS` or `CXXFLAGS`.

So, with a gcc style compiler a good debug configure line would be

    /path/to/configure CFLAGS="-O0 -g -Wall -Wextra" \
                       --enable-debug --enable-mpi

and for a production build use for example

    /path/to/configure CFLAGS="-O2 -Wall -Wno-unused-parameter" \
                       --enable-mpi --disable-shared

Both in-source and out-of-source builds are supported.  For a one-time
installation, `./configure` can be called in the root directory of the unpacked
tree.  For development or for testing multiple configure options, create an
empty directory for each configuration and invoke configure with a relative
path.
Note that the `--disable-shared` option suppresses the build of `.so` shared object lib files.
In case these are also required this option needs to be omitted or changed to `--enable-shared`.

The `configure` call is then followed by the usual `make` and `make install` commands,
the latter potentially in `sudo` mode depending on the desired installation directory.

The subpackage `sc` is contained in the tarball and used by default.  It is
possible to use a version of `libsc` that has been `make install`ed independently:

 1. Create an empty directory and call `sc/configure` with its relative path.
   Use `--prefix=path-to-install-sc` to specify the destination for the
   compiled files.  Then call `make install`.

 2. Create another empty directory and call the `p4est` `configure` with its
   relative path and the options `--prefix=path-to-install-p4est` and
   `--with-sc=path-to-install-sc`.  Make sure to use the same compiler flags
   and configure options as in (1.).  Finally, call `make install`.

This is the proposed procedure to split the `p4est` installation into two
packages, say for a linux binary distribution.  The benefit is that other
packages that might require sc do not force `p4est` to get installed.

## Installation from source

When installing from source, GNU autotools must be invoked to generate the
necessary configuration files and the configure script.  A bootstrap script is
provided for this purpose: give the shell command `./bootstrap` before
following the instructions in the `INSTALL` file.  Then proceed as above.

Note that `autoreconf` does not work in the source directory.
It is trying to grep `ACLOCAL_AMFLAGS` from the `Makefile.am`.
We use this variable for substitution, which autoreconf cannot survive.

## Installation from git repository

If you have obtained `p4est` from a git repository, such as via the shell command

    git clone https://www.github.com/cburstedde/p4est.git

then the libsc submodule, which resides in the subdirectory sc, must be
downloaded before configuring and compiling.  This can be accomplished with the
shell commands

    git submodule init; git submodule update

After the submodule has been obtained, install from source as described above.

## Building a release tarball from a git clone (one sophisticated way)

The following script may be used to build a release tarball.
It assumes an in-source build here, which makes it look relatively complex;
for a simpler procedure see below under *build tarball*.

    #!/bin/sh

    SRC_DIR=$(pwd)
    WORK_DIR=$(mktemp -d)

    cleanup() {
      rm -rf "$WORK_DIR"
      echo ""
      echo "Deleted working directory $WORK_DIR"
    }
    trap cleanup EXIT

    echo "Building p4est tarball in $WORK_DIR"
    echo ""

    # Get p4est
    cd "$WORK_DIR"
    git clone --recursive https://github.com/cburstedde/p4est
    cd p4est

    # Setup exclude files so that version number is not dirty
    ./bootstrap
    ./configure
    git ls-files --others --exclude-standard >> .git/info/exclude
    GIT_DIR=.git/modules/sc git ls-files --others --exclude-standard \
                                     >> .git/modules/sc/info/exclude
    echo ChangeLog >> .git/info/exclude
    echo ChangeLog >> .git/modules/sc/info/exclude
    git clean -xdf
    git submodule foreach git clean -xdf

    # Configure and build distribution tarball
    ./bootstrap
    ./configure
    echo ""
    echo "Building versions"
    echo ""
    git describe --abbrev=4 --match 'v*'
    git submodule foreach git describe --abbrev=4 --match 'v*'
    echo ""
    make -j3 distcheck
    mkdir -p "$SRC_DIR"
    cp p4est-*.tar.gz "$SRC_DIR"

## Building a release tarball from a git clone (the simple way)

The above simplifies a great deal when building in a dedicated empty directory;
basically it reduces to

    # [... Checkout the sc submodule as described above]
    ./bootstrap
    mkdir build && cd build && ../configure
    make -j3 distcheck

## CMake community support

We have added CMake support, which we regard as a community effort.
We are encouraging volunteers to test and possibly improve it.  CMake seems
to be working fine, but neither the pkgconfig files nor the tarball are
currently identical with the ones created by autotools.

p4est offers similar configuration options in CMake and in autoconf.
In particular, the following options are commonly used.

`--Dmpi=on` is necessary for parallel functionality, equivalent to
`--enable-mpi` in autoconf.

`-DCMAKE_BUILD_TYPE=Debug` enables assertions and additional code for
verification, equivalent to `--enable-debug` in autoconf.

A good debug configuration from the source directory would be

    cmake -B build -DCMAKE_BUILD_TYPE=Debug -Dmpi=on

and for a production build use for example

    cmake -B build -Dmpi=on

The configuration is then followed by

    cmake --build build --parallel

Optionally, `p4est` may be installed using

    cmake -B build -DCMAKE_INSTALL_PREFIX=~/local
    cmake --build build --parallel
    cmake --install build

For developers, source and binary distribution packages can be generated by:

    cmake -B build -DCMAKE_INSTALL_PREFIX=~/local
    cmake --build build --parallel
    cpack --config build/CPackSourceConfig.cmake
    cpack --config build/CPackConfig.cmake

which creates files:

 * `build/package/p4est-<version>.tar.gz` containing source code
 * `build/package/p4est-<version>-<platform>.tar.gz` containing binary libraries
    and executables suitable for computers of same operating system and compatible
    CPU arch

Note that the tarball from CPack is structured differently from the one
built using the autotools (see above).  Only the latter adheres to the GNU
file system standards and is recommended for releases.

On Windows above cpack command will yield `.zip`-Files.