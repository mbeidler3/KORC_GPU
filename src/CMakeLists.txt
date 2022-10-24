#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (korc_sources "")
list (APPEND korc_sources
    main.f90
)

add_executable (xkorc ${korc_sources})
set_property(TARGET xkorc PROPERTY LINKER_LANGUAGE Fortran)
