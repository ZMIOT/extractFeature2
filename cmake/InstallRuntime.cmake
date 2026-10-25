# Copy only the actual algorithm DLL closure. Failing on Qt/VTK prevents a
# future source/link change from silently bloating the cloud-only SDK.
get_filename_component(_pcl_root "${PCL_DIR}" DIRECTORY)
set(_runtime_dirs "${_pcl_root}/bin" "${_pcl_root}/dll"
  "${_pcl_root}/3rdParty/Boost/lib" "${_pcl_root}/3rdParty/FLANN/bin")
install(CODE "
  file(GET_RUNTIME_DEPENDENCIES
    LIBRARIES \"$<TARGET_FILE:seam_core>\"
    DIRECTORIES \"$<TARGET_FILE_DIR:seam_core>\" ${_runtime_dirs}
    PRE_EXCLUDE_REGEXES \"api-ms-.*\" \"ext-ms-.*\"
      \"[Vv][Cc][Rr][Uu][Nn][Tt][Ii][Mm][Ee].*\" \"[Mm][Ss][Vv][Cc][Pp].*\" \"[Cc][Oo][Nn][Cc][Rr][Tt].*\"
    POST_EXCLUDE_REGEXES \".*[Ss][Yy][Ss][Tt][Ee][Mm]32.*\" \".*[Ss][Yy][Ss][Ww][Oo][Ww]64.*\"
    RESOLVED_DEPENDENCIES_VAR resolved
    UNRESOLVED_DEPENDENCIES_VAR unresolved)
  if(unresolved)
    message(FATAL_ERROR \"Unresolved runtime DLLs: \${unresolved}\")
  endif()
  foreach(dependency IN LISTS resolved)
    get_filename_component(name \"\${dependency}\" NAME)
    if(name MATCHES \"^([Qq][Tt]|[Vv][Tt][Kk]|pcl_io|pcl_visualization)\")
      message(FATAL_ERROR \"Cloud-only SDK unexpectedly depends on \${name}\")
    endif()
  endforeach()
  file(INSTALL DESTINATION \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/bin\" TYPE SHARED_LIBRARY FILES \${resolved})
")
set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION bin)
include(InstallRequiredSystemLibraries)
