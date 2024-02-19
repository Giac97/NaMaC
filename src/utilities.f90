module utilities
    implicit none
    private
    public: read_xyz
    
    contains

    !> @brief Reads an xyz file, reads the coordinates and saves them to an array
    !!
    !!This code will read from a standard xyz file, doesn't support extended xyz files
    !!
    !!@param[in]    fname   the name of the input xyz file
    !!@param[out]   coordinates an array 3xN_atoms of coordintes
    subroutine read_xyz(fname, coordinates)
        implicit none
        integer :: io, N_atoms, i
        character(len=20), intent(in) :: fname
        character(len=2) :: element
        real :: x_t, y_t, z_t
        real, allocatable :: coordinates(:,:)
        open(newunit=io, file=fname, status="old", action="read")
        read(io, *) N_atoms
        read(io, *) 
        allocate(coordinates(3, N_atoms))
        
        do i = 1, N_atoms
            read(io, *) element, x_t, y_t, z_t
            coordinates(1, i) = x_t
            coordinates(2, i) = y_t
            coordinates(3, i) = z_t
        end do
        close(io)
    end subroutine

    
end module utilities
