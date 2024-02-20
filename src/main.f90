program namac
    use utilities
    use surfacer
    use porosity
    implicit none
    real, allocatable       ::  coordinates(:,:), gcn(:)
    integer, allocatable    ::  coordination(:), neigh_list(:,:), is_surface(:)
    real                    ::  r_atom , r_probe 
    integer                 ::  n_samples, acc, pbc
    real                    ::  pore, CNmax, cutoff, gcn_cutoff
    real, allocatable       ::  insert(:,:)


    character(len=20)   ::  fname
    character(len=50)   ::  outname, gcname, surfname, porname


    fname = "large.xyz"
    namelist/parameters/fname,cutoff,CNmax,gcn_cutoff,pbc
    namelist/porosity/r_atom, r_probe, n_samples
    ! Call the subroutine to read coordinates from a file
    
    open(25, file="input.in", status="old", action="read")
    read(25, parameters)
    read(25, porosity)
    close(25)
    outname = "out_cord.xyz"
    gcname = "out_gcn.xyz"
    surfname = "out_surf.xyz"
    porname = "out_pore.xyz"
    write(*,*) "Welcome to the main NaMaC loop"
    write(*,*) "We will begin by loading the xyz file you want to analyse"
    write(*,*) ""
    call read_xyz(fname, coordinates)
    write(*,*) "Loaded ", size(coordinates, 2), " atoms"
    allocate(coordination(size(coordinates,2)), neigh_list(size(coordinates,2),12))
    
    write(*,*) "Now for the lengthiest part so far, the CN, this will take a while"
    write(*,*)
    call coordination_calc(coordinates, cutoff, pbc, coordination, neigh_list)
    write(*,*) "CN and neighbor list computed"

    write(*,*) "Now the gcn"
    call  gcn_calc(coordination, neigh_list, 12.0, gcn)
    write(*,*) "GCN computed"

    write(*,*) "Outputting all files"
    call write_xyz_coordination(outname, coordinates, coordination)
    call write_xyz_gcn(gcname, coordinates, gcn)
    call find_surface(gcn, is_surface, gcn_cutoff)
    call write_xyz_coordination(surfname, coordinates, is_surface) 
    
    write(*,*) "You thought the cn took a long time?"
    write(*,*) "Well, if you do not have a CUDA GPU then you'll be right"
    write(*,*) "Now here is the longest calculation, this bad boy can not be parallelized due to some exit clauses"
    write(*,*) "Computing the porosity"
    call find_porosity(coordinates, r_atom, r_probe, n_samples, pore, insert, acc)
    write(*,*) "Porosity = ", pore
    call write_xyz_porosity(porname, coordinates, insert, acc)
    deallocate(coordination, coordinates, neigh_list, gcn)
    write(*,*) "Seems like we are done, hopefully nothing crashed along the way"
end program namac
