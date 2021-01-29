option(dev "p4est developer mode")

option(Ep6est "build p6est")
option(Ep8est "build p8est")

option(vtk "VTK interface" on)

option(mpi "use MPI library" on)

option(sc_external "force build of libsc")


if(dev)

else(dev)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_SC true)
endif(dev)
