module surfacer
    use utilities
    implicit none
    private
    public  ::  find_surface,  output_surface_network
    !> @brief Determines whether an atom lies on the surface or is part of the bulk, returns as an array of 1s and 0s 
    !! depending on whether the i-th atoms is on the surface or not
    !!
    !! @param
contains
    subroutine find_surface(gcn, is_surface, gcn_cutoff)
        real, dimension(:), intent(in)          ::  gcn
        real,intent(in)                         ::  gcn_cutoff
        
        integer, allocatable, intent(out)    ::  is_surface(:)
        
        integer     ::  i, j, N_atoms
        
        N_atoms = size(gcn)
        allocate(is_surface(N_atoms))
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT
        do i = 1, N_atoms
            if (gcn(i).le.gcn_cutoff) then
                is_surface(i) = 1
            else 
                is_surface(i) = 0
            endif

        enddo
        !$ACC END KERNELS        
    end subroutine find_surface

    subroutine output_surface_network(coordinates, is_surface, neigh_list, file_name)
        implicit none
        real, intent(in)        ::  coordinates(:,:)
        integer, intent(in)     ::  is_surface(:)
        integer, intent(in)     ::  neigh_list(:,:)
        character(len=50), intent(in)    ::  file_name
        
        integer     ::  N_atoms, i, j
        integer     ::  temp_neigh(12)
        N_atoms = size(is_surface)

        open(15, file=file_name, status="replace", action="write")
        
        do i = 1, N_atoms
            temp_neigh = 0
            if (is_surface(i) .eq. 1) then
                do j = 1, 12
                    if (is_surface(neigh_list(j, i)) .eq. 1) then
                        temp_neigh(j) = neigh_list(j, i)
                    endif
                enddo
                write(15,'(13I8)') i,  temp_neigh(:)
            endif
        enddo

        close(15)
    end subroutine output_surface_network
end module surfacer


