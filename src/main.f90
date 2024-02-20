program namac
    use utilities
    use surfacer
    implicit none
    real, allocatable       ::  coordinates(:,:), gcn(:)
    integer, allocatable    ::  coordination(:), neigh_list(:,:), is_surface(:)
      

    character(len=20)   ::  fname
    character(len=50)   ::  outname, gcname, surfname
    fname = "test.xyz"
    outname = "out_cord.xyz"
    gcname = "out_gcn.xyz"
    surfname = "out_surf.xyz"
    
    call read_xyz(fname, coordinates)
    allocate(coordination(size(coordinates,2)), neigh_list(size(coordinates,2),12))
    call coordination_calc(coordinates, 3.5, 0, coordination, neigh_list)
    call  gcn_calc(coordination, neigh_list, 12.0, gcn)

    call write_xyz_coordination(outname, coordinates, coordination)
    call write_xyz_gcn(gcname, coordinates, gcn)
    call find_surface(gcn, is_surface, 10.0)
    call write_xyz_coordination(surfname, coordinates, is_surface) 
    deallocate(coordination, coordinates, neigh_list, gcn)
end program namac
