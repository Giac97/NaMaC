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
    
    !> @brief Utility function that computes the squared distance between two points
    !!
    !! This subrouitine requires the coordinates of two points and it will compute the squared distance
    !! between them. used mostly when the actual distance is not required to be known to save time
    !! one the calculation of the squared root
    !!
    !!@param[in]    pointA  The coordinates of one of the points
    !!@param[in]    pointB  The coordinates of the other point
    !!@param[out]   distance2   The square of the distance between two points
    subroutine distance2(pointA, pointB, distance2)
        implicit none
        real, intent(in)    ::  pointA(3), pointB(3)
        real, intent(out)   ::  distance2
        real                ::  dist_vec
        
        dist_vec = pointB - pointA     
        distance2 = dist_vec(1) * dist_vec(1) + dist_vec(2) * dist_vec(2) + dist_vec(3) * dist_vec(3)  
    end subroutine

end module utilities
