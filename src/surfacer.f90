module surfacer
    use utilities
    implicit none
    private
    public  ::  find_surface
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
end module surfacer
