target_sources(p4est PRIVATE p4est_base.c p4est_connectivity.c p4est.c p4est_bits.c p4est_search.c p4est_build.c
p4est_algorithms.c p4est_communication.c p4est_ghost.c p4est_nodes.c p4est_points.c p4est_geometry.c p4est_iterate.c
p4est_lnodes.c p4est_mesh.c p4est_balance.c p4est_io.c p4est_connrefine.c
p4est_wrap.c p4est_plex.c p4est_empty.c p4est_vtk.c
)

if(enable_p8est)
  target_sources(p8est PRIVATE p8est_connectivity.c p8est.c p8est_bits.c p8est_search.c p8est_build.c
  p8est_algorithms.c p8est_communication.c p8est_ghost.c p8est_nodes.c p8est_vtk.c p8est_points.c p8est_geometry.c
  p8est_iterate.c p8est_lnodes.c p8est_mesh.c p8est_tets_hexes.c p8est_balance.c p8est_io.c p8est_connrefine.c
  p8est_wrap.c p8est_plex.c p8est_empty.c p8est_vtk.c
  )
endif(enable_p8est)

if(enable_p6est AND enable_p8est)
  target_sources(p6est PRIVATE p6est.c p6est_ghost.c p6est_lnodes.c p6est_profile.c p6est_vtk.c
  p6est_communication.c p6est_empty.c
  )
endif()
