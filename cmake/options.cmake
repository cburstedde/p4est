option(dev "p4est developer mode")

option(Ep4est "build 2D p4est" on)
option(Ep6est "build p6est" off)
option(Ep8est "build p8est" off)

option(mpi "use MPI library" on)

option(sc_external "force build of libsc")


if(dev)

else(dev)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_SC true)
endif(dev)
