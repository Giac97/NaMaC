module utilities
    implicit none
    private
    public: read_xyz, distance2, distance, coordination_calc 
    
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
    !!@param[out]   distance2   The square of the distance between the two points
    subroutine distance2(pointA, pointB, distance2)
        implicit none
        real, intent(in)    ::  pointA(3), pointB(3)
        real, intent(out)   ::  distance2
        real                ::  dist_vec
        
        dist_vec = pointB - pointA     
        distance2 = dist_vec(1) * dist_vec(1) + dist_vec(2) * dist_vec(2) + dist_vec(3) * dist_vec(3)  
    end subroutine distance
    
    !> @brief Utility function that computes the distance between two points
    !!
    !! This subroutine computes the distance between two points. If the actual distance is not essential
    !! and if it has to be computed for a large number of pairs of atoms use the distance2 subroutine
    !! as the calculation of the square root is quite expensive
    !!
    !!@param[in]    pointA  The coordinates of one of the points
    !!@param[in]    pointB  The coordinates of the other point
    !!@param[out]   distance    The distance between the two points
    subroutine distance(pointA, pointB, distance)
        implicit none
        real, intent(in)    ::  pointA(3), pointB(3)
        real, intent(out)   ::  distance
        real                ::  d2

        call distance2(pointA, pointB, d2)
        distance = sqrt(d2)
    end subroutine distance

    !> @brief Computes for each atom its coordination number
    !!
    !! For each atom in the system the coordination number is computed as
    !! the number of other atoms within a certain cutoff distance, if desired
    !! (as set by the pbc variable) it compute the distance using periodic boundary
    !! conditions. Also returns a list of neighbours for each atom as a list of 
    !! indeces pointing to the atoms neighboring the atom
    !!
    !!@param[in]    coordinates 3xN_atoms array of the cooridnates of the atoms
    !!@param[in]    cutoff  real number determining the cutoff for the neighbors 
    !!@param[out]   coordination    array containing for each atom its coordination number
    !!@param[out]   neigh_list  array containing for each atom the indeces of its neighbors
    subroutine coordination_calc(coordinates, cutoff, pbc, coordination, neigh_list)
        implicit none
        integer :: N_atoms, i, j
        real, intent(in) :: coordinates(:, :)
        real, intent(in) :: cutoff
        integer,intent(in)  ::  pbc
        integer, intent(out), allocatable :: coordination(:), neigh_list(:,:)
        real :: dist
        real :: distance_v(3)
        integer :: neighbors, co2
        real :: Lx, Ly
        if ( pbc.eq.1) then
            Lx = maxval(coordinates(1, :)) - minval(coordinates(1,:))
            Ly = maxval(coordinates(2, :)) - minval(coordinates(2,:))
        endif
        N_atoms = size(coordinates, 2)
        allocate(coordination(N_atoms), neigh_list(12,N_atoms))
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT
        do i = 1, N_atoms
            neighbors = 0
            do j = 1, N_atoms
                if (i .ne. j) then
                    distance_v(:) = coordinates(:, i) - coordinates(:, j)
                    if (pbc.eq.1.and.abs(distance_v(1)) .gt. Lx / 2) then
                        distance_v(1) = Lx - abs(distance_v(1))
                    elseif (pbc.eq.1.and.abs(distance_v(2)) .gt. Ly / 2) then
                        distance_v(2) = Ly - abs(distance_v(2))
                    endif
                    call distance(dist)
                    if (dist .lt. cutoff) then
                        neighbors = neighbors + 1
                        neigh_list(neighbors, i) = j
                    endif
                endif
            enddo
        enddo
    end subroutine coordinaion_calc            
end module utilities
