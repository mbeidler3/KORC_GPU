#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (korc_sources "")
list (APPEND korc_sources
    interp_gpu.f90
    params_gpu.f90
    main.f90
)

add_executable (xkorc ${korc_sources})
set_property(TARGET xkorc PROPERTY LINKER_LANGUAGE Fortran)

if (${USE_PSPLINE})
   target_link_libraries (xkorc ${PSPLINE_LIBRARIES})
endif()
