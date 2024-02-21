module porosity
    use utilities
    implicit none
    private
    public  ::  find_porosity, write_xyz_porosity, slice_porosity
contains
    
    !> @brief Computes the porosity of a structure given the coordinates of the atoms
    !!
    !! Using a MOnte Carlo integration scheme, the subroutine computes the porosity of an assembly of atoms
    !! using probe atoms of a given radius
    !!
    !!@param[in]    coordinates Coordinates of the atoms forming the structure
    !!@param[in]    r_atoms     Radius of the atoms forming the structure
    !!@param[in]    r_probe     Radius of the probe atoms (can be 0.0 for ideal point like probes)
    !!@param[in]    n_samples   Total number of attempted insertions
    !!@param[out]   pore        The porosity (acc/n_samples)
    !!@param[out]   insert      Array containing the coordinates of the inserted probe atoms (the non inserted have p = (0,0,0))
    !!@param[out]   acc         Number of accepted insertion (used by write_xyz_porosity)
    subroutine find_porosity(coordinates, r_atoms, r_probe, n_samples, pore, insert, acc)
        implicit none
        real, intent(in)        ::  coordinates(:,:)
        real, intent(in)        ::  r_atoms, r_probe
        integer, intent(in)     ::  n_samples
        
        real, intent(out)       ::  pore
        real, intent(out), allocatable  ::  insert(:,:)
        integer, intent(out)    ::  acc
        real, allocatable       ::  insert_temp(:,:)
        real                    ::  test_insert(3)
        
        real                    ::  x_min, x_max 
        real                    ::  y_min, y_max
        real                    ::  z_min, z_max
        
        real                    ::  dist

        integer                 ::  i, j, N_atoms, rej, seed, idx

        N_atoms = size(coordinates, 2)
        
        x_min = minval(coordinates(1,:))
        x_max = maxval(coordinates(1,:))
        
        y_min = minval(coordinates(2,:))
        y_max = maxval(coordinates(2,:))
        
        z_min = minval(coordinates(3,:))
        z_max = maxval(coordinates(3,:))
        allocate(insert(3,n_samples))
        
        acc = 0
        seed = 176 
        idx = 0
        call random_seed()
        do i = 1, n_samples
            rej = 0
            call random_number(test_insert)
            test_insert(1) = (x_max - x_min) * test_insert(1) + x_min
            test_insert(2) = (y_max - y_min) * test_insert(2) + y_min
            test_insert(3) = (z_max - z_min) * test_insert(3) + z_min
            do j = 1, N_atoms
                call distance(test_insert, coordinates(:,j), dist)
                if (dist - (r_atoms + r_probe) .lt. 0) then
                    rej = 1
                    exit 
                endif
            enddo
            
            if (rej .eq. 0) then
                acc = acc + 1
                insert(:,i) = test_insert
                idx = idx + 1
            endif
            
            if (mod(i,1000).eq.0) then
                write(*,*) "Sample ", i, " of ", n_samples, ", porosity= ", real(acc) / real(i)
            endif
        enddo

        pore = real(acc) / real(n_samples)
    end subroutine find_porosity

    !> @brief Computes the porosity of a structure given the coordinates of the atoms in a box with fixed x, y sides
    !!
    !! Using a MOnte Carlo integration scheme, the subroutine computes the porosity of an assembly of atoms
    !! using probe atoms of a given radius in a box with fixed x, y sides. Used in computing the porosity in slices to 
    !! avoid shrinking 
    !!
    !!@param[in]    coordinates Coordinates of the atoms forming the structure
    !!@param[in]    r_atoms     Radius of the atoms forming the structure
    !!@param[in]    r_probe     Radius of the probe atoms (can be 0.0 for ideal point like probes)
    !!@param[in]    n_samples   Total number of attempted insertions
    !!@param[out]   pore        The porosity (acc/n_samples)
    !!@param[out]   insert      Array containing the coordinates of the inserted probe atoms (the non inserted have p = (0,0,0))
    !!@param[out]   acc         Number of accepted insertion (used by write_xyz_porosity)
    !!@param[in]    x_min       The lower x bound of the box
    !!@param[in]    x_max       The upper x bound of the box
    !!@param[in]    y_min       The lower y bound of the box
    !!@param[in]    y_max       The upper y bound of the box
    subroutine find_porosity_fixed_box(coordinates, r_atoms, r_probe, n_samples, pore, insert, acc, x_min, x_max, y_min, y_max)
        implicit none
        real, intent(in)        ::  coordinates(:,:)
        real, intent(in)        ::  r_atoms, r_probe
        integer, intent(in)     ::  n_samples
        
        real, intent(out)       ::  pore
        real, intent(out), allocatable  ::  insert(:,:)
        integer, intent(out)    ::  acc
        real, allocatable       ::  insert_temp(:,:)
        real                    ::  test_insert(3)
        
        real, intent(in)        ::  x_min, x_max 
        real, intent(in)        ::  y_min, y_max
        real                    ::  z_min, z_max
        
        real                    ::  dist

        integer                 ::  i, j, N_atoms, rej, seed, idx

        N_atoms = size(coordinates, 2)
        
        
        
        z_min = minval(coordinates(3,:))
        z_max = maxval(coordinates(3,:))
        allocate(insert(3,n_samples))
        
        acc = 0
        seed = 176 
        idx = 0
        call random_seed()
        do i = 1, n_samples
            rej = 0
            call random_number(test_insert)
            test_insert(1) = (x_max - x_min) * test_insert(1) + x_min
            test_insert(2) = (y_max - y_min) * test_insert(2) + y_min
            test_insert(3) = (z_max - z_min) * test_insert(3) + z_min
            do j = 1, N_atoms
                call distance(test_insert, coordinates(:,j), dist)
                if (dist - (r_atoms + r_probe) .lt. 0) then
                    rej = 1
                    exit 
                endif
            enddo
            
            if (rej .eq. 0) then
                acc = acc + 1
                insert(:,i) = test_insert
                idx = idx + 1
            endif
            
            if (mod(i,1000).eq.0) then
                write(*,*) "Sample ", i, " of ", n_samples, ", porosity= ", real(acc) / real(i)
            endif
        enddo

        pore = real(acc) / real(n_samples)
    end subroutine find_porosity_fixed_box
    subroutine write_xyz_porosity(fname, coordinates, insert, acc)
        implicit none
        character(len=50), intent(in)   ::  fname
        real, intent(in)                ::  coordinates(:,:), insert(:,:)
        integer, intent(in)             ::  acc
        integer                         ::  N_atoms,  i, N_samples

        N_atoms = size(coordinates, 2)
        N_samples = size(insert, 2)
        open(10, file=fname, status="replace", action="write")
        write(10, *) N_atoms + acc
        write(10, *) "# Writing porosity and atoms to xyz"
        do i = 1, N_atoms
            write(10, '(A2, 3F15.10)') "Au", coordinates(:,i)
        enddo

        do i = 1, N_samples
            if (abs(insert(1,i)) .ge. 0.0003 .and. abs(insert(2,i)) .ge. 0.0002 .and. abs(insert(3,i)) .ge. 0.0002) then
                write(10, '(A2, 3F15.10)') "XX", insert(:,i)
            endif
        enddo

        close(10)
            
    end subroutine
    !> @brief Uses the MC method to compute the porosity in slices of the material from the bottom to the top,
    !!
    !!@param[in]    coordinates Array of the coordinates of the atoms in the system
    !!@param[in]    r_atoms     Radius of the atoms in the system
    !!@param[in]    r_probe     Radius of the fictitious probe atoms
    !!@param[in]    n_samples   The number of samples per slice
    !!@param[in]    N_layers    The number of slices in the system
    !!@param[out]   slice_file  The output file for the porosity at each slice
    subroutine slice_porosity(coordinates, r_atoms, r_probe, n_samples, N_layers, slice_file)
        implicit none
        real, intent(in)        ::  coordinates(:,:)
        real, intent(in)        ::  r_atoms, r_probe
        integer, intent(in)     ::  n_samples
        integer, intent(in)     ::  N_layers
        character(len=50)       ::  slice_file, finalp
        
        integer     ::  i, j, N_per_slice, idx, acc, N_atoms
        real        ::  pore, Lz, dz, z_min, z_max
        real        ::  x_min, x_max, y_min, y_max
        real, allocatable   ::  insert(:,:)
        real, allocatable   ::  slice(:,:)        
        allocate(insert(3,n_samples))
        N_atoms = size(coordinates, 2)
        x_min = minval(coordinates(1,:))
        x_max = maxval(coordinates(1,:))
        
        y_min = minval(coordinates(2,:))
        y_max = maxval(coordinates(2,:))
        Lz = maxval(coordinates(3,:)) - minval(coordinates(3,:))
        dz = Lz / real(N_layers)
        open(50, file=slice_file, status="replace", action="write")
        z_min = minval(coordinates(3,:))
        do i = 1, N_layers
            write(*,*) "Evaluating porosity in slice ", i, " of ", N_layers
            pore = 0
            acc = 0
            
            z_max = z_min + dz
            N_per_slice = count(coordinates(3,:) .lt. z_max .and. coordinates(3,:) .ge. z_min)
            write(*,*) "Number of atoms in slice ", z_min ," to ", z_max, " = ", N_per_slice
            allocate(slice(3, N_per_slice))    
            idx = 1
            do j = 1, N_atoms
                if (coordinates(3,j) .lt. z_max .and. coordinates(3,j) .ge. z_min) then
                    slice(:,idx) = coordinates(:,j)
                    idx = idx + 1
                endif
            enddo
            call find_porosity_fixed_box(slice, r_atoms, r_probe, n_samples, pore, insert, acc, x_min, x_max, y_min, y_max)
            write(*,*) "finished slice"
            write(50, '(2I5, 3F10.5)')N_per_slice, i, pore, z_min, z_max
            z_min = z_min + dz
            write(*,*) z_min
            deallocate(slice)
        enddo
        close(50)
        deallocate(insert)
    end subroutine slice_porosity
        

end module porosity
