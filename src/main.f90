program namac
    use utilities
    implicit none
    real, allocatable   ::  coordinates(:,:)
    integer, allocatable    ::  coordination(:), neigh_list(:,:)

    character(len=20)   ::  fname
    fname = "test.xyz"
    call read_xyz(fname, coordinates)
    allocate(coordination(size(coordinates,2)), neigh_list(size(coordinates,2),12))
    call coordination_calc(coordinates, 3.5, 0, coordination, neigh_list)

    deallocate(coordination, coordinates, neigh_list)
end program namac
