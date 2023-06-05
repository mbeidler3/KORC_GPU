#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (korcgpu_sources "")
list (APPEND korcgpu_sources
    interp_gpu.f90
    params_gpu.f90
    main.f90
)

add_executable (xkorcgpu ${korcgpu_sources})
set_property(TARGET xkorcgpu PROPERTY LINKER_LANGUAGE Fortran)

if (${USE_PSPLINE})
   target_link_libraries (xkorcgpu ${PSPLINE_LIBRARIES})
endif()