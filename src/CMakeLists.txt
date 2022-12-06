#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (korc_sources "")
list (APPEND korc_sources
    main.f90
    io.F90
    physdat.f90
    pardata.F90
    global.F90
    input.F90
    main_nim.f90
    map_mod.f90
    exchange.f90
    lagrange_quad.F90
    get_fields.f90
    nim_stop.f90
    tblock.F90
    computation_pointers.f90
    edge.f90
    hdf5_fcapi.F90
    hdf5_psapi.f90
    local.F90
    mpi_serial.f90
    dump_ser.F90
    spline.f90
    bicube.f90
    time.f90
    memlog.F90
    lagr_disc_quad.f90
    lagrange_edge.F90
    modal_disc_quad.f90
    rblock_type_mod.F90
    tri_linear.f90
    vector_type_mod.F90
    tblock_type_mod.F90
    fields.f90
    send_rblock_mod.f90
    polynomials.f90
    math_tran.f90
    edge_type_mod.f90
    seam_storage.F90
    normalize.f90
    rt_swap.f90
    matrix_type_mod.F90
    fft2d.F90
    fft_mod.F90
    rblock.F90
    parallel.F90
)

add_executable (xkorc ${korc_sources})
set_property(TARGET xkorc PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries( xkorc ${MPI_Fortran_LIBRARIES} )
target_link_libraries( xkorc ${HDF5_Fortran_LIBRARIES} )