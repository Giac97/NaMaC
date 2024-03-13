program namac
    use utilities
    use surfacer
    use porosity
    implicit none
    real, allocatable       ::  coordinates(:,:), gcn(:), strain(:)
    integer, allocatable    ::  coordination(:), neigh_list(:,:), is_surface(:)
    real                    ::  r_atom , r_probe 
    integer                 ::  n_samples, acc, pbc, i
    integer                 ::  N_slices, samples_per_slice
    real                    ::  pore, CNmax, cutoff, gcn_cutoff, err
    real, allocatable       ::  insert(:,:)


    character(len=20)   ::  fname, infile, arg, mode
    character(len=50)   ::  outname, gcname, surfname, porname, slice_file, net_name, strain_name
    net_name = "network.dat"
    strain_name = "out_strain.xyz"
    infile = "input.in"
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        
        if (trim(arg) .eq. "--help" .or. trim(arg) .eq. "-h") then
            print*, ""
            print*, "use -i or --infile to enter the name of the input file"
            stop 
        else if (trim(arg).eq. "-i" .or. trim(arg).eq."--infile") then
            if (i < command_argument_count()) then
                call get_command_argument(i + 1, arg)
                read(arg, *) infile
            else
                print *, "File name requires a value"
                stop
            endif
        endif

    enddo

    mode = "full"
    fname = "large.xyz"
    namelist/parameters/fname,cutoff,CNmax,gcn_cutoff,pbc, mode
    namelist/porosity/r_atom, r_probe, n_samples
    namelist/slice/N_slices, samples_per_slice
    ! Call the subroutine to read coordinates from a file
    
    open(25, file=infile, status="old", action="read")
    read(25, parameters)
    read(25, porosity)
    read(25, slice)
    close(25)
    outname = "out_cord.xyz"
    gcname = "out_gcn.xyz"
    surfname = "out_surf.xyz"
    porname = "out_pore.xyz"
    slice_file = "slice_poro.dat"
    write(*,*) "Welcome to the main NaMaC loop"
    write(*,*) "We will begin by loading the xyz file you want to analyse"
    write(*,*) ""
    call read_xyz(fname, coordinates)
    write(*,*) "Loaded ", size(coordinates, 2), " atoms"
    allocate(coordination(size(coordinates,2)), neigh_list(size(coordinates,2),12))
    if (mode.eq."full".or.mode.eq."surf") then   
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
        write(*,*) "Fraction of atoms on surface: ", real(sum(is_surface)) / real(size(coordinates, 2))
        call write_xyz_coordination(surfname, coordinates, is_surface) 
        call  output_surface_network(coordinates, is_surface, neigh_list, net_name)
    else if (mode.eq."full".or.mode.eq."porofull") then 
        write(*,*) "You thought the cn took a long time?"
        write(*,*) "Well, if you do not have a CUDA GPU then you'll be right"
        write(*,*) "Now here is the longest calculation, this bad boy can not be parallelized due to some exit clauses"
        write(*,*) "Computing the porosity"
        call find_porosity(coordinates, r_atom, r_probe, n_samples, pore, insert, acc, err)
        write(*,*) "Porosity = ", pore, " +/- ", err
        call write_xyz_porosity(porname, coordinates, insert, acc)
    else if (mode.eq."full".or.mode.eq."slice") then
        call slice_porosity(coordinates, r_atom, r_probe, samples_per_slice, N_slices, slice_file)

    else if (mode.eq."strain") then
        allocate(strain(size(coordinates, 2)))
        write(*,*) "Computation of the strain"
        write(*,*) "Fisrst step: building neighbour list"
        call coordination_calc(coordinates,  cutoff, pbc, coordination, neigh_list)
        write(*,*) "Neighbour list built"
        call strain_system(coordinates, pbc,neigh_list, coordination,r_atom, strain)
        write(*,*) "Strain computed"
        call write_xyz_strain(strain_name, coordinates, coordination, strain) 
    endif


    deallocate(coordination, coordinates, neigh_list, gcn, strain)
    write(*,*) "Seems like we are done, hopefully nothing crashed along the way"
end program namac
