
# Suggestion thanks to @pkestene and @VictorEijkhout
# gnu compatibility : https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include(GNUInstallDirs)

# --- generate pkg-config .pc
set(pc_libs_private)
set(pc_req_private "ompi ompi-c orte zlib")

set(pc_req_public "p4est libsc")

set(pc_filename p4est-${git_version}.pc)
configure_file(${CMAKE_CURRENT_LIST_DIR}/pkgconf.pc.in ${pc_filename} @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${pc_filename}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

set(pc_target ${pc_filename})
set(pc_link ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig/p4est.pc)

install(CODE "file(CREATE_LINK ${pc_target} ${pc_link} SYMBOLIC)")
