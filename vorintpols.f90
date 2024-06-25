module VorCells_Polygons

  ! Types and subroutines to deal with Voronoi cells and convex polygons.
  
  PRIVATE :: tol,SortVertCell,quicksort,EdgesNeighbors,SutherlandHodgman,edgeClipping,intersection,visible,crossProduct
  PUBLIC :: Polygon,VoronoiCell,CurvilinearPolygon,voronoi_cell,VorCellInterConvPol,VorCellInterConvPol_Collinear,interCellBall
  
  ! Tolerance used in some subroutines.
  real(kind=8) :: tol = 1000.0D+00 * epsilon ( tol )
 
  ! Type for convex polygons.
  type Polygon

    ! number of vertices: to define a polygon, the first and the last vertices have to be the same.
    integer :: n
    
    ! vertices of the polygon.
    real(kind=8), dimension(:,:), allocatable :: vertex
    
    ! A flag to indicate whether a polygon is degenerate.
    logical :: deg
    
    ! Edge identification: each edge of a polygon can be of three types: 0 => edge is contained in the boundary of
    ! A; 1 => edge is in int(A); a negative integer X in {-nBalls,...,-1} => edge is contained in some Voronoi edge with
    ! identification given by X.
    integer, dimension(:), allocatable :: edges
    
  end type Polygon
  
  ! Type for Voronoi cells. It is used when there are at least 3 noncollinear generators (ball's centers).
  type VoronoiCell
    
    ! flag to indicate whether the cell is bounded or unbounded: flag = true => bounded (cell is a convex polygon).
    logical :: flag
    
    ! number of vertices: for bounded cells, the first and the last vertices must be the same.
    integer :: n
    
    ! vertices of the cell.
    real(kind=8), dimension(:,:), allocatable :: vertex
    
    ! Rays: if flag = false, this array contains the origin and direction of two rays that compose unbounded cells.
    ! [origin, origin, direction, direction], (1--3 and 2--4). If flag = true, this array will not be used.
    real(kind=8) :: rays(2,4)
    
    ! Voronoi edge identification: each Voronoi edge L (or also a ray when the cell is unbounded) of a Voronoi cell V_i is
    ! identified with a negative integer flag {-1,...,-nBalls}, that represents the index j of the center whose 
    ! Voronoi cell V_j is neighbor of V_i and are separated by L.
    integer, dimension(:), allocatable :: edges
    
  end type VoronoiCell
  
  ! Type for curvilinear polygons, resulting from the intersection of a convex polygon with a circle.
  type CurvilinearPolygon
    
    ! number of vertices: the first and the last vertices have to be the same.
    integer :: n
    
    ! vertices
    real(kind=8), dimension(:,:), allocatable :: vertex
    
    ! vertex identification: an integer to identify the vertex type.
    ! Each vertex can be the intersection of a ball with: (i) boundary of A; (ii) an auxiliary edge from the partition of A
    ! in convex sets; (iii) another ball (a Voronoi edge). The identification will be done as:
    ! (i) 0; (ii) 1; (iii) a negative integer in {-nBalls, ..., -1}.
    ! Furthermore, the vertex can be a corner of the convex polygon. In this case, we will identify the vertex as 2.
    integer, dimension(:), allocatable :: vertex_id
    
    ! A flag to indicate whether the polygon is degenerate or empty.
    logical :: deg
    
    ! integer array to identify straight and curved edges: 0 => straight edge
    !                                                      1 => arc
    integer, dimension(:), allocatable :: edges
    
  end type CurvilinearPolygon
  
  contains
        
    ! -------------------------------------------------------- !
    subroutine SortVertCell(center, vorcell)
    ! Sort the vertices of a Voronoi cell in counter-clockwise order. For unbounded cells, it is supposed to be
    ! at least three Voronoi vertices.
        
        implicit none
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: center(2)
        
        ! TYPE ARGUMENTS
        type(VoronoiCell), intent(inout) :: vorcell
        
        ! LOCAL SCALARS
        integer :: ve,vj
        
        ! LOCAL ARRAYS
        real(kind=8) :: angles(vorcell%n - 1),angles_sorted(vorcell%n - 1),vertices_sorted(2,vorcell%n - 1),last(2)
        
        ! FUNCTIONS
        real(kind=8) angle_rad_2d
        
        ! Bounded cells
        if ( vorcell%flag ) then
            ! Sorting from the angle obtained with respect to the center.
            do ve=1,vorcell%n - 1
                angles(ve) = datan2( vorcell%vertex(2,ve) - center(2), vorcell%vertex(1,ve) - center(1) )
            end do
            angles_sorted = angles
            ! Sort the angles.
            call quicksort(angles_sorted, vorcell%n - 1, 1, vorcell%n - 1)
            do vj = 1,vorcell%n - 1
                do ve = 1,vorcell%n - 1
                    if ( angles_sorted(vj) .eq. angles(ve)  ) then
                        vertices_sorted(1:2,vj) = vorcell%vertex(1:2,ve)
                    end if
                end do
            end do
            ! Store in vorcell%vertex.
            vorcell%vertex(1:2,1:vorcell%n - 1) = vertices_sorted(1:2,1:vorcell%n - 1)
            vorcell%vertex(1:2,vorcell%n) = vorcell%vertex(1:2,1)
        else
            ! First and last vertices stored in vorcell%vertex are also origins of the rays. 
            if ( visible( vorcell%vertex(1:2,2), vorcell%vertex(1:2,1), vorcell%vertex(1:2,vorcell%n) ) ) then
                ! Origin of vector is the last sorted vertex.
                do ve=2,vorcell%n
                    ! Sort from the angle in radians swept out between the rays 
                    ! ( vorcell%vertex(:,1), vorcell%vertex(:,ve) ) and 
                    ! ( vorcell%vertex(:,1), vorcell%vertex(:,vorcell%n) ).
                    angles(ve - 1) = angle_rad_2d( vorcell%vertex(1:2,ve), vorcell%vertex(1:2,1),&
                    vorcell%vertex(1:2,vorcell%n) )
                end do
                angles_sorted = angles
                ! Sort the angles.
                call quicksort(angles_sorted, vorcell%n - 1, 1, vorcell%n - 1)
                
                do vj = 2,vorcell%n
                    do ve = 2,vorcell%n
                        if ( angles_sorted(vj - 1) .eq. angles(ve - 1)  ) then
                            vertices_sorted(1:2,vj - 1) = vorcell%vertex(1:2,ve)
                        end if
                    end do
                end do
                
                ! Store in vorcell%vertex.
                last(1:2) = vorcell%vertex(1:2,1)
                vorcell%vertex(1:2,1:vorcell%n - 1) = vertices_sorted(1:2,1:vorcell%n - 1)
                vorcell%vertex(1:2,vorcell%n) = last(1:2)
            else
                ! Last vertex already sorted.
                do ve=1,vorcell%n - 1
                    ! Sort from the angle in radians swept out between the rays 
                    ! (vorcell%vertex(:,vorcell%n), vorcell%vertex(:,ve) ) and 
                    ! (vorcell%vertex(:,vorcell%n), vorcell%vertex(:,1) ).
                    angles(ve) = angle_rad_2d( vorcell%vertex(1:2,ve), vorcell%vertex(1:2,vorcell%n),&
                    vorcell%vertex(1:2,1) )
                end do
                angles_sorted = angles
                ! Sort the angles.
                call quicksort(angles_sorted, vorcell%n - 1, 1, vorcell%n - 1)
                do vj = 1,vorcell%n - 1
                    do ve = 1,vorcell%n - 1
                        if ( angles_sorted(vj) .eq. angles(ve) ) then
                            vertices_sorted(1:2,vj) = vorcell%vertex(1:2,ve)
                        end if
                    end do
                end do
                ! Store in vorcell%vertex.
                vorcell%vertex(1:2,1:vorcell%n - 1) = vertices_sorted(1:2,1:vorcell%n - 1)
            end if
        end if
    end subroutine SortVertCell
    
    ! -------------------------------------------------------- !
    recursive subroutine quicksort(a, n, first, last)
    ! Sort elements in array a with size n from a(first) to a(last).
        
        implicit none
        
        ! SCALAR ARGUMENTS
        integer, intent(in) :: n, first, last
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(inout) :: a(n)
        
        ! LOCAL SCALARS
        integer :: i,j
        real(kind=8) :: x,t
        
        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) .lt. x)
                i=i+1
            end do
            do while (x .lt. a(j))
                j=j-1
            end do
            if (i .ge. j) exit
            t = a(i)
            a(i) = a(j)
            a(j) = t
            i = i + 1
            j= j - 1
        end do
        if (first .lt. i-1) call quicksort(a, n, first, i-1)
        if (j+1 .lt. last)  call quicksort(a, n, j+1, last)
    end subroutine quicksort
    
    ! -------------------------------------------------------- !
    subroutine voronoi_cell(ind_ball, nBalls, centers, v_num, nod_tri, tnbr, vor_xy, vor_cell, status)
        ! Obtain the Voronoi cell (vor_cell) of the generator (ball center) ind_ball.
        
        implicit none
        
        ! SCALARS ARGUMENTS
        integer, intent(in) :: ind_ball, nBalls, v_num
        integer, intent(out) :: status
        
        ! ARRAYS ARGUMENTS
        integer, intent(in) :: nod_tri(3,v_num),tnbr(3,v_num) ! Delaunay triangulation data.
        real(kind=8), intent(in) :: centers(2,nBalls), vor_xy(2,v_num) ! Ball's centers and Voronoi vertices.
        
        ! TYPE ARGUMENTS
        type(VoronoiCell), intent(out) :: vor_cell
        
        ! LOCAL SCALARS
        integer :: cur_neigh,i,iv,in,it,itu,j,k,l,piv,v,repeated_vertex
        logical :: TwoFacesHaveNbrtTriang
        real(kind=8) :: absinnprod,min_absinnprod
        
        ! LOCAL ARRAYS
        integer :: jj(2),neigh_jj(2),triangles_used(v_num),repeated_origins(2)
        real(kind=8) :: dir(2),intermediate_vertices(2,v_num)
        
        ! FUNCTIONS
        ! Euclidean distance of two points.
        real(kind=8) :: r8vec_distance
        
        status = 0
        
        ! Storing the Voronoi vertices and triangle indexes of the cell.
        it = 0
        iv = 0
        in = 0
        vor_cell%flag = .true.
        do v=1,v_num ! For each triangle v
            repeated_vertex = 0
            do j=1,3 ! For each vertex of triangle v
                ! Check if center ind_ball is a vertex of the triangle v.
                if ( ind_ball .eq. nod_tri(j,v) ) then
                    
                    ! Obtaining the indexes of the edges and neighbors with respect to the vertex ind_ball.
                    call EdgesNeighbors(j, jj, neigh_jj)
                    
                    ! Store the triangles used.
                    it = it + 1
                    triangles_used(it) = v
                    
                    ! When both edges jj(1) and jj(2) are also edges of another triangle (neighbor to v), then the Voronoi
                    ! vertex associated to the triangle v is not an origin of a ray, that is, it is an 
                    ! intermediate vertex. We will distinguish the possible cases using the variable below.
                    TwoFacesHaveNbrtTriang = .true.
                    
                    ! Prevent two intermediate vertices from being repeated (too close).
                    if ( iv .gt. 0 ) then
                        do piv = 1,iv
                            if ( r8vec_distance(2,vor_xy(:,v),intermediate_vertices(:,piv)) .lt. tol ) then
                                repeated_vertex = 1
                                exit
                            end if
                        end do
                    end if
                    
                    do l=1,2
                        ! Check whether jj(l) is also edge of another triangle.
                        k = tnbr(jj(l),v)
                        if ( k .lt. 0 ) then ! jj(l) IS NOT edge of another triangle.
                            
                            vor_cell%flag = .false.
                            TwoFacesHaveNbrtTriang = .false.
                            
                            ! Store the vertex (origin of the ray).
                            in = in + 1
                            vor_cell%rays(1:2,in) = vor_xy(1:2,v)
                            repeated_origins(in) = repeated_vertex
                            
                            ! Prevent two origins from being repeated.
                            if ( in .eq. 2 ) then
                                !if ( vor_cell%rays(1,in) .eq. vor_cell%rays(1,1) .and. &
                                !vor_cell%rays(2,in) .eq. vor_cell%rays(2,1) ) then
                                if ( r8vec_distance(2,vor_cell%rays(:,in),vor_cell%rays(:,1)) .lt. tol ) then
                                    repeated_origins(in) = 1
                                end if
                            end if
                            
                            ! Compute and store the outward normal to the triangle edge (direction of the ray).
                            if (jj(l) .lt. 3) then
                                call line_exp_normal_2d( centers(1:2, nod_tri(jj(l), v)), & 
                                centers(1:2, nod_tri(jj(l)+1, v)), vor_cell%rays(1:2,in+2) )
                                if ( jj(l) .eq. 1 .and. dot_product( vor_cell%rays(:,in+2),&
                                centers(:,nod_tri(3,v)) - centers(:,nod_tri(1,v)) ) .gt. 0.0d0 ) then
                                    vor_cell%rays(:,in+2) = - 1.0d0 * vor_cell%rays(:,in+2)
                                end if
                                if ( jj(l) .eq. 2 .and. dot_product( vor_cell%rays(:,in+2),&
                                centers(:,nod_tri(1,v)) - centers(:,nod_tri(2,v)) ) .gt. 0.0d0 ) then
                                    vor_cell%rays(:,in+2) = - 1.0d0 * vor_cell%rays(:,in+2)
                                end if
                            else
                                call line_exp_normal_2d( centers(1:2, nod_tri(3, v)), centers(1:2, nod_tri(1, v)), &
                                vor_cell%rays(1:2,in+2) )
                                if ( dot_product( vor_cell%rays(:,in+2), centers(:,nod_tri(2,v)) - centers(:,nod_tri(3,v)) )&
                                .gt. 0.0d0 ) then
                                    vor_cell%rays(:,in+2) = - 1.0d0 * vor_cell%rays(:,in+2)
                                end if
                            end if
                        end if
                    end do

                    ! Prevent an intermediate vertex and an origin from being repeated.
                    if ( v .gt. 1 .and. in .gt. 0 ) then
                        if ( in .eq. 1 ) then
                            !if ( vor_xy(1,v) .eq. vor_cell%rays(1,1) .and. vor_xy(2,v) .eq. vor_cell%rays(2,1) ) then
                            if ( r8vec_distance(2,vor_xy(:,v),vor_cell%rays(:,1)) .lt. tol ) then
                                repeated_vertex = 1
                            end if
                        end if
                        if ( in .eq. 2 ) then
                            !if ( vor_xy(1,v) .eq. vor_cell%rays(1,1) .and. vor_xy(2,v) .eq. vor_cell%rays(2,1) ) then
                            if ( r8vec_distance(2,vor_xy(:,v),vor_cell%rays(:,1)) .lt. tol ) then
                                repeated_vertex = 1
                            end if
                            !if ( vor_xy(1,v) .eq. vor_cell%rays(1,2) .and. vor_xy(2,v) .eq. vor_cell%rays(2,2) ) then
                            if ( r8vec_distance(2,vor_xy(:,v),vor_cell%rays(:,2)) .lt. tol ) then
                                repeated_vertex = 1
                            end if
                        end if
                    end if
                    
                    ! Adding the intermediate vertices.
                    if ( TwoFacesHaveNbrtTriang .and. repeated_vertex .eq. 0 ) then
                        iv = iv + 1
                        intermediate_vertices(1:2,iv) = vor_xy(1:2,v)
                    end if
                end if
            end do
        end do
        
        ! Put the vor_cell%n vertices in vor_cell%vertex.
        if ( vor_cell%flag ) then
            vor_cell%n = iv + 1
            allocate(vor_cell%vertex(2,vor_cell%n))
            vor_cell%vertex(1:2,1:vor_cell%n - 1) = intermediate_vertices(1:2,1:iv)
            ! Sorting the vertices.
            call SortVertCell(centers(1:2,ind_ball), vor_cell)
        else
            ! Adding the endpoints (origin of the rays). Sort the vertices.
            if ( repeated_origins(1) .eq. 0 .and. repeated_origins(2) .eq. 0 ) then
                vor_cell%n = iv + 2
                allocate(vor_cell%vertex(2,vor_cell%n))
                vor_cell%vertex(1:2,1) = vor_cell%rays(1:2,1)
                vor_cell%vertex(1:2,vor_cell%n) = vor_cell%rays(1:2,2)
                if ( iv .gt. 0 ) then
                    vor_cell%vertex(1:2,2:vor_cell%n-1) = intermediate_vertices(1:2,1:iv)
                end if
            else if ( repeated_origins(1) .eq. 0 .and. repeated_origins(2) .eq. 1 ) then
                vor_cell%n = iv + 1
                allocate(vor_cell%vertex(2,vor_cell%n))
                vor_cell%vertex(1:2,1) = vor_cell%rays(1:2,1)
                if ( iv .gt. 0 ) then
                    vor_cell%vertex(1:2,2:vor_cell%n) = intermediate_vertices(1:2,1:iv)
                end if
            else if ( repeated_origins(1) .eq. 1 .and. repeated_origins(2) .eq. 0 ) then
                vor_cell%n = iv + 1
                allocate(vor_cell%vertex(2,vor_cell%n))
                if ( iv .gt. 0 ) then
                    vor_cell%vertex(1:2,1:vor_cell%n-1) = intermediate_vertices(1:2,1:iv)
                end if
                vor_cell%vertex(1:2,vor_cell%n) = vor_cell%rays(1:2,2)
            else
                vor_cell%n = iv
                allocate(vor_cell%vertex(2,vor_cell%n))
                vor_cell%vertex(1:2,1:vor_cell%n) = intermediate_vertices(1:2,1:vor_cell%n)
            end if
            if ( vor_cell%n .gt. 2 ) then
                call SortVertCell(centers(1:2,ind_ball), vor_cell)
            end if
        end if

        
        ! Identifying the edges of the Voronoi cell.
        if ( vor_cell%flag ) then
            
            ! Bounded cells
            allocate(vor_cell%edges(vor_cell%n - 1))
            
            do i = 1,vor_cell%n - 1
                min_absinnprod = 1.0d+20
                do itu = 1,it
                    if ( r8vec_distance(2,vor_xy(:,triangles_used(itu)),vor_cell%vertex(:,i)) .lt. tol ) then
                        do j=1,3
                            if ( ind_ball .eq. nod_tri(j,triangles_used(itu)) ) then
                                call EdgesNeighbors(j, jj, neigh_jj)
                                exit
                            end if
                        end do
                        do l=1,2
                            absinnprod = dabs(dot_product(centers(:,nod_tri(neigh_jj(l),triangles_used(itu)))-centers(:,ind_ball),&
                            vor_cell%vertex(:,i+1)-vor_cell%vertex(:,i)))
                            if ( absinnprod .lt. min_absinnprod ) then
                                min_absinnprod = absinnprod
                                cur_neigh = nod_tri(neigh_jj(l),triangles_used(itu))
                            end if
                        end do
                    end if
                end do
                
                vor_cell%edges(i) = - cur_neigh 
                
                ! Check whether the index is correct.
                if ( .not. (1 .le. iabs(vor_cell%edges(i)) .and. iabs(vor_cell%edges(i)) .le. nBalls) ) then
                    write(*,*) 'The edge ID is not between 1 and nBalls.'
                    status = -21
                    return
                end if
                
            end do
        else
            ! Unbounded cells
            allocate(vor_cell%edges(vor_cell%n + 1))
            
            ! First ray
            min_absinnprod = 1.0d+20
            do itu = 1,it
                if ( r8vec_distance(2,vor_xy(:,triangles_used(itu)),vor_cell%vertex(:,1)) .lt. tol ) then
                    do j=1,3
                        if ( ind_ball .eq. nod_tri(j,triangles_used(itu)) ) then
                            call EdgesNeighbors(j, jj, neigh_jj)
                            exit
                        end if
                    end do
                    if ( vor_cell%vertex(1,1) .eq. vor_cell%rays(1,1) .and. &
                    vor_cell%vertex(2,1) .eq. vor_cell%rays(2,1) ) then
                        dir(:) = vor_cell%rays(:,3)
                    else
                        dir(:) = vor_cell%rays(:,4)
                    end if
                    do l=1,2
                        absinnprod = dabs(dot_product(centers(:,nod_tri(neigh_jj(l),triangles_used(itu)))-centers(:,ind_ball),&
                        dir))
                        if ( absinnprod .lt. min_absinnprod ) then
                            min_absinnprod = absinnprod
                            cur_neigh = nod_tri(neigh_jj(l),triangles_used(itu))
                        end if
                    end do
                end if
            end do
            
            vor_cell%edges(1) = - cur_neigh
            
            ! Check whether the index is correct.
            if ( .not. (1 .le. iabs(vor_cell%edges(1)) .and. iabs(vor_cell%edges(1)) .le. nBalls) ) then
                write(*,*) 'The edge ID is not between 1 and nBalls.'
                status = -22
                return
            end if
            
            ! Second ray
            min_absinnprod = 1.0d+20
            do itu = 1,it
                if ( r8vec_distance(2,vor_xy(:,triangles_used(itu)),vor_cell%vertex(:,vor_cell%n)) .lt. tol ) then
                    do j=1,3
                        if ( ind_ball .eq. nod_tri(j,triangles_used(itu)) ) then
                            call EdgesNeighbors(j, jj, neigh_jj)
                            exit
                        end if
                    end do
                    if ( vor_cell%vertex(1,vor_cell%n) .eq. vor_cell%rays(1,2) .and. &
                    vor_cell%vertex(2,vor_cell%n) .eq. vor_cell%rays(2,2) ) then
                        dir(:) = vor_cell%rays(:,4)
                    else
                        dir(:) = vor_cell%rays(:,3)
                    end if
                    do l=1,2
                        absinnprod = dabs(dot_product(centers(:,nod_tri(neigh_jj(l),triangles_used(itu)))-centers(:,ind_ball),&
                        dir))
                        if ( absinnprod .lt. min_absinnprod ) then
                            min_absinnprod = absinnprod
                            cur_neigh = nod_tri(neigh_jj(l),triangles_used(itu))
                        end if
                    end do
                end if
            end do
            
            vor_cell%edges(vor_cell%n + 1) = - cur_neigh
            
            ! Check whether the index is correct.
            if (.not. (1 .le. iabs(vor_cell%edges(vor_cell%n + 1)) .and. iabs(vor_cell%edges(vor_cell%n + 1)) .le. nBalls)) then
                write(*,*) 'The edge ID is not between 1 and nBalls.'
                status = -23
                return
            end if
            
            ! Intermediate edges
            if ( vor_cell%n .ge. 2 ) then
                do i = 1,vor_cell%n - 1
                    min_absinnprod = 1.0d+20
                    do itu = 1,it
                        if ( r8vec_distance(2,vor_xy(:,triangles_used(itu)),vor_cell%vertex(:,i)) .lt. tol ) then
                            do j=1,3
                                if ( ind_ball .eq. nod_tri(j,triangles_used(itu)) ) then
                                    call EdgesNeighbors(j, jj, neigh_jj)
                                    exit
                                end if
                            end do
                            do l=1,2
                                absinnprod = dabs(dot_product(centers(:,nod_tri(neigh_jj(l),triangles_used(itu)))&
                                -centers(:,ind_ball),vor_cell%vertex(:,i+1)-vor_cell%vertex(:,i)))
                                if ( absinnprod .lt. min_absinnprod ) then
                                    min_absinnprod = absinnprod
                                    cur_neigh = nod_tri(neigh_jj(l),triangles_used(itu))
                                end if
                            end do
                        end if
                    end do
                    
                    vor_cell%edges(i+1) = - cur_neigh
                    
                    ! Check whether the index is correct.
                    if ( .not. (1 .le. iabs(vor_cell%edges(i+1)) .and. iabs(vor_cell%edges(i+1)) .le. nBalls) ) then
                        write(*,*) 'The edge ID is not between 1 and nBalls.'
                        status = -24
                        return
                    end if
                    
                end do
            end if
        end if
        
    end subroutine voronoi_cell
    
    ! -------------------------------------------------------- !
    subroutine EdgesNeighbors(j, jj, neigh_jj)
        ! Given any triangle with vertices nod_tri(j,:), j=1,2, or 3, this routine computes the indices jj of the edges
        ! that contain the vertex with index j. Also, computes the indexes of each neighbor vertex of nod_tri(j,:) in jj.
        implicit none
        
        ! SCALAR ARGUMENTS
        integer, intent(in) :: j
        
        ! ARRAY ARGUMENTS
        integer, intent(out) :: jj(2),neigh_jj(2)
        
        if ( j .eq. 1 ) then
            jj(1) = 1
            neigh_jj(1) = 2
            jj(2) = 3
            neigh_jj(2) = 3
        end if
        if ( j .eq. 2 ) then
            jj(1) = 2
            neigh_jj(1) = 3
            jj(2) = 1
            neigh_jj(2) = 1
        end if
        if ( j .eq. 3 ) then
            jj(1) = 3
            neigh_jj(1) = 1
            jj(2) = 2
            neigh_jj(2) = 2
        end if
    end subroutine EdgesNeighbors
    
    ! -------------------------------------------------------- !
    subroutine VorCellInterConvPol(center, convPol, vor_cell, cell)
        ! This subroutine computes the intersection of a convex polygon (convPol) with a Voronoi cell (vor_cell), 
        ! returning a convex polygon (cell).
        implicit none
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: center(2) ! generator of the Voronoi cell.
        
        ! TYPE ARGUMENTS
        type(VoronoiCell), intent(in) :: vor_cell
        type(Polygon), intent(in) :: convPol
        type(Polygon), intent(out) :: cell

        ! LOCAL SCALARS
        integer :: aux,aux_edge
        
        ! LOCAL ARRAYS
        real(kind=8) :: aux_vec(2)
        
        ! LOCAL TYPES
        type(VoronoiCell) :: clipping_window_cell
        
        ! FUNCTIONS
        real(kind=8) :: r8vec_distance
        
        ! Allocate the maximum number of vertices and edges in the output cell.
        allocate(cell%vertex(2, convPol%n + vor_cell%n + 2))
        allocate(cell%edges(convPol%n + vor_cell%n + 1))
        
        ! Bounded cells
        if ( vor_cell%flag ) then
            ! Compute the intersection.
            call SutherlandHodgman( vor_cell, convPol, cell )
        else ! Unbounded cells
            aux = vor_cell%n + 2
            ! To pass the unbounded cell as an argument of the SH algorithm, we need to add an auxiliary vertex for each half-line.
            clipping_window_cell%flag = .false.
            allocate( clipping_window_cell%vertex(2,aux), clipping_window_cell%edges(vor_cell%n + 1) )
            clipping_window_cell%vertex(:,2:aux - 1) = vor_cell%vertex(:,1:vor_cell%n)
            clipping_window_cell%edges(1:vor_cell%n + 1) = vor_cell%edges(1:vor_cell%n + 1)
            clipping_window_cell%rays(:,1:4) = vor_cell%rays(:,1:4)
            clipping_window_cell%n = aux
            if ( vor_cell%n .eq. 1 ) then ! One vertex => two half-lines.
                ! Adding the vertices in the correct order: center must be visible with respect to the line passing through
                ! the points, moving from the first to second point added, in this order.
                if ( visible( center, clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,3),&
                clipping_window_cell%vertex(:,2) ) ) then
                    clipping_window_cell%vertex(:,1) = clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,3)
                    clipping_window_cell%vertex(:,aux) = clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,4)
                else
                    clipping_window_cell%vertex(:,1) = clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,4)
                    clipping_window_cell%vertex(:,aux) = clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,3)
                    ! Correct the order of edge identification.
                    aux_edge = clipping_window_cell%edges(1)
                    clipping_window_cell%edges(1) = clipping_window_cell%edges(2)
                    clipping_window_cell%edges(2) = aux_edge
                end if
                ! Compute the intersection.
                call SutherlandHodgman( clipping_window_cell, convPol, cell )
            else if ( vor_cell%n .eq. 2 ) then ! Two vertices => one edge and two half-lines.
                ! Adding the vertices in the correct order: center must be visible with respect to the line passing through
                ! the points, moving from the first to second point added, in this order.
                if ( visible( center, clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,3),&
                clipping_window_cell%vertex(:,2) ) ) then
                    clipping_window_cell%vertex(:,1) = clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,3)
                    clipping_window_cell%vertex(:,aux) = clipping_window_cell%vertex(:,aux - 1) + clipping_window_cell%rays(:,4)
                else
                    aux_vec(1:2) = clipping_window_cell%vertex(:,2)
                    clipping_window_cell%vertex(:,1) = clipping_window_cell%vertex(:,aux - 1) + clipping_window_cell%rays(:,4)
                    clipping_window_cell%vertex(:,aux) = clipping_window_cell%vertex(:,2) + clipping_window_cell%rays(:,3)
                    clipping_window_cell%vertex(:,2) = clipping_window_cell%vertex(:,aux - 1)
                    clipping_window_cell%vertex(:,aux - 1) = aux_vec(1:2)
                    ! Correct the order of edge identification.
                    aux_edge = clipping_window_cell%edges(1)
                    clipping_window_cell%edges(1) = clipping_window_cell%edges(3)
                    clipping_window_cell%edges(3) = aux_edge
                end if
                call SutherlandHodgman( clipping_window_cell, convPol, cell )
            else ! Two half-lines and clipping_window_cell%n - 1 edges.
                ! Choose the origin of the rays correctly, since, in this case, all vertices were previously sorted.
                if ( r8vec_distance(2,clipping_window_cell%vertex(:,2),clipping_window_cell%rays(:,1)) .le. tol ) then
                    clipping_window_cell%vertex(:,1) = clipping_window_cell%rays(:,1) + clipping_window_cell%rays(:,3)
                    clipping_window_cell%vertex(:,aux) = clipping_window_cell%rays(:,2) + clipping_window_cell%rays(:,4)
                else
                    clipping_window_cell%vertex(:,1) = clipping_window_cell%rays(:,2) + clipping_window_cell%rays(:,4)
                    clipping_window_cell%vertex(:,aux) = clipping_window_cell%rays(:,1) + clipping_window_cell%rays(:,3)
                end if
                ! Compute the intersection.
                call SutherlandHodgman( clipping_window_cell, convPol, cell )
            end if
            deallocate( clipping_window_cell%vertex, clipping_window_cell%edges )
        end if
    end subroutine VorCellInterConvPol
    
    ! -------------------------------------------------------- !
    subroutine VorCellInterConvPol_Collinear(ind_ball, nBalls, centers, convPol, cell)
        ! This subroutine computes the intersection of a convex polygon (convPol) with the Voronoi cell ind_ball
        ! when the generators (ball's centers) are collinear. Return: a convex polygon (cell).
        implicit none
        
        ! SCALAR ARGUMENTS
        integer, intent(in) :: ind_ball,nBalls
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: centers(2,nBalls)
        
        ! TYPE ARGUMENTS
        type(Polygon), intent(in) :: convPol
        type(Polygon), intent(out) :: cell

        ! LOCAL SCALARS
        integer :: i,neighbor,neighbor2
        real(kind=8) :: d1,d2,mind1,mind2
        
        ! LOCAL ARRAYS
        real(kind=8) :: midpoint(2),normal(2)
        
        ! LOCAL TYPES
        type(Polygon) :: aux_pol
        type(VoronoiCell) :: clipping_window_cell
        
        ! FUNCTIONS
        real(kind=8) :: r8vec_distance
        
        neighbor = 0
        neighbor2 = 0
        mind1 = 1.0d+15
        mind2 = 1.0d+15
        
        ! Check that the centers do not form a vertical line.
        if ( dabs( dot_product(centers(:,2) - centers(:,1), (/ 1.0d0, 0.0d0 /)) ) .gt. tol ) then
            ! Compute the neighbors of center ind_ball.
            do i = 1,nBalls
                if ( i .ne. ind_ball ) then
                    ! Look at the points on the left side of ind_ball.
                    if ( centers(1,i) - centers(1,ind_ball) .lt. 0.0d0 ) then
                        ! Find the point closest to ind_ball.
                        d1 = r8vec_distance(2,centers(:,i),centers(:,ind_ball))
                        if ( d1 .lt. mind1 ) then
                            mind1 = d1
                            ! Set the neighbor.
                            neighbor = i
                        end if
                    end if
                    ! Look at the points on the right side of ind_ball.
                    if ( centers(1,i) - centers(1,ind_ball) .gt. 0.0d0 ) then
                        ! Find the point closest to ind_ball.
                        d2 = r8vec_distance(2,centers(:,i),centers(:,ind_ball))
                        if ( d2 .lt. mind2 ) then
                            mind2 = d2
                            ! Set the neighbor.
                            neighbor2 = i
                        end if
                    end if
                    ! Check whether there are overlapped balls.
                    if ( dabs( centers(1,i) - centers(1,ind_ball) ) .lt. tol ) then
                        write(*,*) 'There are centers too close!'
                        stop
                    end if
                end if
            end do
        else ! The centers form a vertical line.
            do i = 1,nBalls
                if ( i .ne. ind_ball ) then
                    ! Look at the points below ind_ball.
                    if ( centers(2,i) - centers(2,ind_ball) .lt. 0.0d0 ) then
                        ! Find the point closest to ind_ball.
                        d1 = r8vec_distance(2,centers(:,i),centers(:,ind_ball))
                        if ( d1 .lt. mind1 ) then
                            mind1 = d1
                            ! Set the neighbor.
                            neighbor = i
                        end if
                    end if
                    ! Look at the points above ind_ball.
                    if ( centers(2,i) - centers(2,ind_ball) .gt. 0.0d0 ) then
                        ! Find the point closest to ind_ball.
                        d2 = r8vec_distance(2,centers(:,i),centers(:,ind_ball))
                        if ( d2 .lt. mind2 ) then
                            mind2 = d2
                            ! Set the neighbor.
                            neighbor2 = i
                        end if
                    end if
                    ! Check whether there are overlapped balls.
                    if ( dabs( centers(2,i) - centers(2,ind_ball) ) .lt. tol ) then
                        write(*,*) 'There are centers too close!'
                        stop
                    end if
                end if
            end do
        end if
        
        ! Defining a half-plane as an unbounded cell.
        clipping_window_cell%n = 2 ! Two points on the boundary (line).
        clipping_window_cell%flag = .false.
        allocate(clipping_window_cell%vertex(2,clipping_window_cell%n))
        
        ! Checking if center ind_ball is at one end. In this case, there is only one neighbor.
        if ( neighbor .eq. 0 .or. neighbor2 .eq. 0 ) then
            if ( neighbor2 .ne. 0 ) then
                neighbor = neighbor2
            end if
            ! Calculate the points on the line that delimit the half-plane.
            midpoint(1) = ( centers(1,ind_ball) + centers(1,neighbor) ) / 2.0d0
            midpoint(2) = ( centers(2,ind_ball) + centers(2,neighbor) ) / 2.0d0
            call line_exp_normal_2d ( centers(:,ind_ball), centers(:,neighbor), normal )
            ! Adding the two points to compose the boundary of the half-plane. They must be visible to the center ind_ball.
            if ( visible( centers(:,ind_ball), midpoint, midpoint(1:2) + normal(1:2) ) ) then
                clipping_window_cell%vertex(:,1) = midpoint(1:2)
                clipping_window_cell%vertex(:,2) = midpoint(1:2) + normal(1:2)
            else
                clipping_window_cell%vertex(:,1) = midpoint(1:2) + normal(1:2)
                clipping_window_cell%vertex(:,2) = midpoint(1:2)
            end if
            ! Adding the edge (line) identification.
            allocate(clipping_window_cell%edges(1))
            clipping_window_cell%edges(1) = - neighbor
            ! Allocate the maximum number of vertices and edges in cell.
            allocate(cell%vertex(2, convPol%n + clipping_window_cell%n))
            allocate(cell%edges(convPol%n + clipping_window_cell%n - 1))
            ! Intersection of convPol with the half-plane.
            call SutherlandHodgman(clipping_window_cell, convPol, cell)
            deallocate(clipping_window_cell%vertex, clipping_window_cell%edges)
        else
            ! In this case, there are two neighbors.
            ! Calculate the points on the line that delimit the first half-plane.
            midpoint(1) = ( centers(1,ind_ball) + centers(1,neighbor) ) / 2.0d0
            midpoint(2) = ( centers(2,ind_ball) + centers(2,neighbor) ) / 2.0d0
            call line_exp_normal_2d ( centers(:,ind_ball), centers(:,neighbor), normal )
            ! Adding the two points to compose the boundary of the first half-plane. They must be visible to the center ind_ball.
            if ( visible( centers(:,ind_ball), midpoint, midpoint(1:2) + normal(1:2) ) ) then
                clipping_window_cell%vertex(:,1) = midpoint(1:2)
                clipping_window_cell%vertex(:,2) = midpoint(1:2) + normal(1:2)
            else
                clipping_window_cell%vertex(:,1) = midpoint(1:2) + normal(1:2)
                clipping_window_cell%vertex(:,2) = midpoint(1:2)
            end if
            ! Adding the edge (line) identification.
            allocate(clipping_window_cell%edges(2))
            clipping_window_cell%edges(1) = - neighbor
            ! Allocate the maximum number of vertices and edges in aux_pol.
            allocate(aux_pol%vertex(2, convPol%n + clipping_window_cell%n))
            allocate(aux_pol%edges(convPol%n + clipping_window_cell%n - 1))
            ! Intersection of convPol with the first half-plane.
            call SutherlandHodgman(clipping_window_cell, convPol, aux_pol)
            ! Calculate the points on the line that delimit the second half-plane.
            midpoint(1) = ( centers(1,ind_ball) + centers(1,neighbor2) ) / 2.0d0
            midpoint(2) = ( centers(2,ind_ball) + centers(2,neighbor2) ) / 2.0d0
            call line_exp_normal_2d ( centers(:,ind_ball), centers(:,neighbor2), normal )
            ! Adding the two points to compose the boundary of the second half-plane. They must be visible to the center ind_ball.
            if ( visible( centers(:,ind_ball), midpoint, midpoint(1:2) + normal(1:2) ) ) then
                clipping_window_cell%vertex(:,1) = midpoint(1:2)
                clipping_window_cell%vertex(:,2) = midpoint(1:2) + normal(1:2)
            else
                clipping_window_cell%vertex(:,1) = midpoint(1:2) + normal(1:2)
                clipping_window_cell%vertex(:,2) = midpoint(1:2)
            end if
            ! Adding the edge (line) identification.
            clipping_window_cell%edges(2) = - neighbor2
            ! Allocate the maximum number of vertices and edges in cell.
            allocate(cell%vertex(2, aux_pol%n + clipping_window_cell%n))
            allocate(cell%edges(convPol%n + clipping_window_cell%n - 1))
            ! Intersection of aux_pol with the second half-plane.
            call SutherlandHodgman(clipping_window_cell, aux_pol, cell)
            deallocate(clipping_window_cell%vertex, clipping_window_cell%edges, aux_pol%vertex)
        end if
    end subroutine VorCellInterConvPol_Collinear
    
    ! -------------------------------------------------------- !
    subroutine interCellBall(center, r, cell, curv_pol)
        ! This subroutine computes the intersection of a convex polygon (cell) with the ball B(center, r).
        ! Return: a curvilinear polygon (curv_pol).
        implicit none
        
        ! SCALAR ARGUMENTS
        real(kind=8), intent(in) :: r
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: center(2)
        
        ! TYPE ARGUMENTS
        type(Polygon), intent(in) :: cell
        type(CurvilinearPolygon), intent(out) :: curv_pol

        ! LOCAL SCALARS
        integer :: i,iv,ie,int_num
        integer :: last_added ! is equal to 1 if the last vertex added is an intersection and 0 otherwise.
        real(kind=8) :: d1,d2,d
        
        ! LOCAL ARRAYS
        real(kind=8) :: p(2,2)
        
        ! FUNCTIONS
        real(kind=8) :: r8vec_distance
        
        ! Allocating the maximum number of vertices for the intersection.
        allocate( curv_pol%vertex( 2, 2 * cell%n ) )
        allocate( curv_pol%vertex_id( 2 * cell%n ) )
        allocate( curv_pol%edges( 2 * (cell%n - 1) ) )
        
        iv = 0
        ie = 0
        last_added = 0
        curv_pol%n = 0
        curv_pol%deg = .false.
        
        do i = 1,cell%n - 1
            ! Check if two consecutive points are too close.
            d = r8vec_distance(2, cell%vertex(1:2,i), cell%vertex(1:2,i+1))
            if ( d .lt. tol) then
                go to 10
            end if
            
            ! Intersection between circle and line segment with endpoints cell%vertex(1:2,i) and cell%vertex(1:2,i+1).
            call circle_imp_segment_intersect ( r, center, cell%vertex(1:2,i), cell%vertex(1:2,i+1), int_num, p )
            
            ! Two intersections.
            if ( int_num .eq. 2 ) then
                if ( iv .gt. 0 ) then
                    ie = ie + 1
                    curv_pol%edges(ie) = 1
                end if
                d1 = r8vec_distance ( 2, cell%vertex(:,i), p(:,1) )
                d2 = r8vec_distance ( 2, cell%vertex(:,i), p(:,2) )
                if ( d1 .lt. d2 ) then
                    iv = iv + 1
                    curv_pol%vertex(1:2,iv) = p(1:2,1)
                    curv_pol%vertex_id(iv) = cell%edges(i)
                    iv = iv + 1
                    curv_pol%vertex(1:2,iv) = p(1:2,2)
                    curv_pol%vertex_id(iv) = cell%edges(i)
                else
                    iv = iv + 1
                    curv_pol%vertex(1:2,iv) = p(1:2,2)
                    curv_pol%vertex_id(iv) = cell%edges(i)
                    iv = iv + 1
                    curv_pol%vertex(1:2,iv) = p(1:2,1)
                    curv_pol%vertex_id(iv) = cell%edges(i)
                end if
                last_added = 1
                ie = ie + 1
                curv_pol%edges(ie) = 0
            ! One intersection.
            else if ( int_num .eq. 1 ) then
                if ( iv .gt. 0 ) then
                    if ( last_added .eq. 1 ) then
                        ie = ie + 1
                        curv_pol%edges(ie) = 1
                    else
                        ie = ie + 1
                        curv_pol%edges(ie) = 0
                    end if
                end if
                iv = iv + 1
                curv_pol%vertex(1:2,iv) = p(1:2,1)
                curv_pol%vertex_id(iv) = cell%edges(i)
                last_added = 1
                if ( (cell%vertex(1,i+1) - center(1))**2 + (cell%vertex(2,i+1) - center(2))**2 .lt. r**2 ) then
                    iv = iv + 1
                    curv_pol%vertex(1:2,iv) = cell%vertex(1:2,i+1)
                    curv_pol%vertex_id(iv) = 2
                    last_added = 0
                    ie = ie + 1
                    curv_pol%edges(ie) = 0
                end if
            ! No intersection.
            else
                if ( (cell%vertex(1,i+1) - center(1))**2 + (cell%vertex(2,i+1) - center(2))**2 .lt. r**2 ) then
                    if ( iv .gt. 0 ) then
                        ie = ie + 1
                        curv_pol%edges(ie) = 0
                    end if
                    iv = iv + 1
                    curv_pol%vertex(1:2,iv) = cell%vertex(1:2,i+1)
                    curv_pol%vertex_id(iv) = 2
                    last_added = 0
                end if
            end if
10          continue  
        end do
        
        ! Number of vertices of curv_pol.
        curv_pol%n = iv
        
        ! Check whether is degenerate.
        if (curv_pol%n .eq. 0) then
            curv_pol%deg = .true.
        else
            ! adding the last vertex (equal to the first) and the corresponding edge type.
            iv = iv + 1
            curv_pol%n = iv
            curv_pol%vertex(1:2,iv) = curv_pol%vertex(1:2,1)
            curv_pol%vertex_id(iv) = curv_pol%vertex_id(1)
            ie = ie + 1
            if ( last_added .eq. 1 ) then
                curv_pol%edges(ie) = 1
            else
                curv_pol%edges(ie) = 0
            end if
        end if            
    end subroutine interCellBall
    
    ! -------------------------------------------------------- !
    subroutine SutherlandHodgman( W, P, Q )
        ! Sutherland-Hodgman algorithm for 2d clipping.
        
        implicit none
        
        ! W is the clipping window (a Voronoi cell).
        ! P is the input polygon.
        ! Q is the output polygon, i.e., intersection of W and P.
        
        ! TYPE ARGUMENTS
        type(VoronoiCell), intent(in) :: W
        type(Polygon), intent(in) :: P
        type(Polygon), intent(inout) :: Q
    
        ! LOCAL SCALARS
        integer :: i
        real(kind=8) :: res
        
        ! LOCAL ARRAYS
        real(kind=8) :: y1(2), y2(2)    ! vertices of the edge to clip workPolygon.
        
        ! LOCAL TYPE
        type(Polygon) :: workPolygon    ! polygon clipped step by step.
        
        ! allocate workPolygon with the maximal possible size: the sum of the number of vertices in W and P.
        allocate( workPolygon%vertex( 2, W%n + P%n ), workPolygon%edges(W%n + P%n - 1) )
        
        ! initialize the work polygon with P.
        workPolygon%n = P%n
        workPolygon%vertex(:,1:workPolygon%n) = P%vertex(:,1:workPolygon%n)
        workPolygon%edges(1:W%n + P%n - 1) = 111
        workPolygon%edges(1:workPolygon%n - 1) = P%edges(1:workPolygon%n - 1)
        workPolygon%deg = .false.

        do i=1,W%n - 1
            y1(:) = W%vertex(:,i)
            y2(:) = W%vertex(:,i+1)
            ! Compute the intersection of workPolygon with the half-plane delimited by the line passing through (y1,y2).
            call edgeClipping( workPolygon, y1, y2, W%edges(i), Q )
            workPolygon%n = Q%n
            workPolygon%vertex(:,1:workPolygon%n) = Q%vertex(:,1:workPolygon%n)
            if ( workPolygon%n .gt. 0 ) then
                workPolygon%edges(1:workPolygon%n - 1) = Q%edges(1:workPolygon%n - 1)
            end if
        end do
        
        ! Check whether Q is degenerate.
        Q%deg = .false.
        if ( Q%n .le. 3 ) then ! Remember that first and last vertices are equal.
            Q%deg = .true.
        end if
        if ( Q%n .gt. 3 ) then
            ! Compute the area of the polygon.
            call polygon_1_2d ( Q%n - 1, Q%vertex(1:2,1:Q%n - 1), res )
            if ( abs(res) .le. tol ) then
                Q%deg = .true.
            end if
        end if
        
        deallocate(workPolygon%vertex)
    end subroutine SutherlandHodgman
    
    ! -------------------------------------------------------- !
    subroutine edgeClipping( poly, y1, y2, id_vor_edge, outputPoly )
        ! make the clipping of the polygon by the line (y1--y2)
        
        implicit none
        
        ! TYPE ARGUMENTS
        type(Polygon), intent(in) :: poly
        type(Polygon), intent(inout) :: outputPoly
        
        ! SCALAR ARGUMENTS
        integer, intent(in) :: id_vor_edge
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: y1(2), y2(2)
        
        ! LOCAL ARRAYS
        real(kind=8) :: x1(2), x2(2), intersecPoint(2)
        
        ! LOCAL SCALARS
        integer :: i,c,last_edge
        real(kind=8) :: d
        logical :: intsct
        
        ! FUNCTIONS
        real(kind=8) :: r8vec_distance
        
        ! counter for the output polygon
        c = 0
        
        ! Check if line (y1--y2) is degenerate
        d = r8vec_distance(2,y1,y2)
        if ( d .le. tol ) then
            outputPoly%n = poly%n
            outputPoly%vertex(:,1:outputPoly%n) = poly%vertex(:,1:outputPoly%n)
            outputPoly%edges(1:outputPoly%n - 1) = poly%edges(1:outputPoly%n - 1)
            return
        end if
    
        do i = 1,poly%n-1 ! for each edge i of poly
            x1(:) = poly%vertex(:,i)   ! vertex 1 of edge i
            x2(:) = poly%vertex(:,i+1) ! vertex 2 of edge i
            if ( visible(x1, y1, y2) ) then ! if x1 is visible
                if ( visible(x2, y1, y2) ) then ! if x2 is visible
                    ! add x2 to the output polygon
                    if ( c .eq. 0 .or. (c .gt. 0 .and. r8vec_distance(2,outputPoly%vertex(:,c),x2) .gt. tol) ) then
                        c = c+1
                        outputPoly%vertex(:,c) = x2(:)
                        outputPoly%edges(c) = poly%edges(i)
                    end if
                else ! x2 is not visible
                    call intersection( x1, x2, y1, y2, intersecPoint, intsct )
                    if ( .not. intsct ) then
                        ! add x2 to the output polygon
                        if ( c .eq. 0 .or. (c .gt. 0 .and. r8vec_distance(2,outputPoly%vertex(:,c),x2) .gt. tol) ) then
                            c = c+1
                            outputPoly%vertex(:,c) = x2(:)
                            outputPoly%edges(c) = poly%edges(i)
                        end if
                    end if
                    if ( c .eq. 0 .or. (c .gt. 0 .and. r8vec_distance(2,outputPoly%vertex(:,c),intersecPoint) .gt. tol) ) then
                        c = c+1
                        outputPoly%vertex(:,c) = intersecPoint(:)
                        outputPoly%edges(c) = poly%edges(i)
                    end if
                end if
            else ! x1 is not visible
                if ( visible(x2, y1, y2) ) then ! x2 is visible
                    call intersection( x1, x2, y1, y2, intersecPoint, intsct )
                    if ( .not. intsct ) then
                        if ( c .eq. 0 .or. (c .gt. 0 .and. r8vec_distance(2,outputPoly%vertex(:,c),x2) .gt. tol) ) then
                            c = c+1
                            outputPoly%vertex(:,c) = x2(:)
                            outputPoly%edges(c) = poly%edges(i)
                        end if
                    else
                        if ( c .eq. 0 .or. (c .gt. 0 .and. r8vec_distance(2,outputPoly%vertex(:,c),intersecPoint) .gt. tol) ) then
                            c = c+1
                            outputPoly%vertex(:,c) = intersecPoint(:)
                            outputPoly%edges(c) = id_vor_edge
                        end if
                        if ( r8vec_distance(2,outputPoly%vertex(:,c),x2) .gt. tol ) then
                            c = c+1
                            outputPoly%vertex(:,c) = x2(:)
                            outputPoly%edges(c) = poly%edges(i)
                        end if
                    end if
                end if
            end if
        end do
        
        ! Correct the order of the edges.
        if (c .gt. 0) then
            last_edge = outputPoly%edges(1)
            outputPoly%edges(1:c-1) = outputPoly%edges(2:c)
            outputPoly%edges(c) = last_edge
        end if
        
        ! Last and first vertices must be the same.
        if (c .gt. 0) then
            d = r8vec_distance(2,outputPoly%vertex(:,1),outputPoly%vertex(:,c))
            if ( d .gt. tol ) then
                c = c+1
                outputPoly%vertex(:,c) = outputPoly%vertex(:,1)
            end if
        end if
        
        ! set the size of the outputPolygon
        outputPoly%n = c

    end subroutine edgeClipping
    
    ! -------------------------------------------------------- !
    subroutine intersection( x1, x2, y1, y2, inter, intsct )
        ! It computes the intersection (if it exists) between a segment [x1,x2] and the line passing through y1 and y2.
        implicit none
        
        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x1(2), x2(2), y1(2), y2(2)
        real(kind=8), intent(out) :: inter(2)
        
        ! SCALAR ARGUMENTS
        logical, intent(out) :: intsct
        
        ! LOCAL ARRAYS
        real(kind=8) :: p(2)
        
        call xedge(1, y1(1), y1(2), y2(1), y2(2), x1(1), x1(2), x2(1), x2(2), p(1), p(2), intsct )
        if ( intsct ) then
            inter(1:2) = p(1:2)
        else
            call xedge(1, y2(1), y2(2), y1(1), y1(2), x1(1), x1(2), x2(1), x2(2), p(1), p(2), intsct )
            if ( intsct ) then
                inter(1:2) = p(1:2)
            end if
        end if
    
    end subroutine intersection
    
    ! -------------------------------------------------------- !
    function visible( p, y1, y2 )
        ! This function checks whether p belongs to the half-plane to the left of line (y1--y2).
        ! (i.e., p is "visible"?).
        
        implicit none
        
        ! ARRAY ARGUMENTS
        real(kind=8) :: p(2), y1(2), y2(2)
        
        ! LOCAL SCALAR
        logical :: visible
        
        ! LOCAL ARRAYS
        real(kind=8) :: v1(2), v2(2)
        
        v1(:) = y2(:) -  y1(:)
        v2(:) = p(:)  -  y1(:)
        
        if ( crossProduct(v1,v2) .ge. 0.0d0) then
            visible = .true.
        else 
            visible = .false.
        end if
        
    end function visible
    
    ! -------------------------------------------------------- !
    function crossProduct( v1, v2 )
        ! It compute the crossproduct of vectors v1 and v2.
        
        implicit none
        
        ! ARRAY ARGUMENTS
        real(kind=8) :: v1(2), v2(2)
        
        ! LOCAL SCALAR
        real(kind=8) :: crossProduct
        
        crossProduct = v1(1)*v2(2) - v1(2)*v2(1)
        
    end function crossProduct
 
end module VorCells_Polygons
