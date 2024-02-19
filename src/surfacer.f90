module surfacer
    use utilities
    implicit none
    private
    public: surfacer
    !> @brief Determines whether an atom lies on the surface or is part of the bulk, returns as an array of 1s and 0s 
    !! depending on whether the i-th atoms is on the surface or not
    !!
    !! @param
    subroutine surfacer(gcn, is_surface, gcn_cutoff)
        real, dimension(:), intent(in)    ::  gcn
        real,intent(in9                     ::  gcn_cutoff
        real, allocatable(:), intent(out)   ::  is_surface
        
        integer     ::  i, j, N_atoms
        
        N_atoms = size(gcn)
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT
        do i = 1, N_atoms
            if (gcn.le.gcn_cutoff) is_surface(i) = 1
            else
                is_surface(i) = 0
            endif

        enddo
        
    end subroutine surfacer
end module surfacer
