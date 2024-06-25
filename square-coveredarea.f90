program squarecoveredarea

  use iso_c_binding, only: c_ptr, c_loc,c_f_pointer
  use VorCells_Polygons

  implicit none
  
  ! PARAMETERS
  integer, parameter :: n = 5
  
  ! PROBLEM DATA TYPE
  type :: pdata_type
     integer :: npols,maxnvpols
     integer, dimension(:), allocatable :: nvpols
     integer, dimension(:,:), allocatable :: poledges
     real(kind=8), dimension(:,:), allocatable :: xpol,ypol
  end type pdata_type

  ! LOCAL SCALARS
  integer :: i,j,allocerr
  type(pdata_type), target :: pdata
  real(kind=8) :: finish,covarea,h,start
  
  ! LOCAL ARRAYS
  real(kind=8) :: x(n),hvec(5)
  
  ! COMMON SCALARS
  integer :: nballs

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdatafirst/ dbl,dtr,nballs
  
  ! Problem data:
  
  nballs = 2
  
  hvec(1:5) = (/ 1.0d-01, 1.0d-02, 1.0d-03, 1.0d-04, 1.0d-05 /)
  
  dbl(1) = 0.0d0
  dbl(2) = 0.0d0
  dtr(1) = 3.0d0
  dtr(2) = 3.0d0
  
  ! Convex polygons
  pdata%npols = 1
  pdata%maxnvpols = 5
  
  allocate( pdata%nvpols(pdata%npols), pdata%poledges(pdata%npols, pdata%maxnvpols), &
            pdata%xpol(pdata%npols, pdata%maxnvpols), pdata%ypol(pdata%npols, pdata%maxnvpols), stat=allocerr )
  
  ! Number of vertices for each polygon
  pdata%nvpols(1:pdata%npols) = (/ 5 /)
  
  ! Vertices
  pdata%xpol(1, 1:pdata%nvpols(1)) = (/ 0.0d0, 3.0d0, 3.0d0, 0.0d0, 0.0d0 /)
  pdata%ypol(1, 1:pdata%nvpols(1)) = (/ 0.0d0, 0.0d0, 3.0d0, 3.0d0, 0.0d0 /)

  ! Edges ID
  pdata%poledges(1, 1:pdata%nvpols(1) - 1) = (/ 0, 0, 0, 0 /)
  
  ! Config
  
  x(1:2) = (/ 0.0d0, 3.0d0 /)
  x(3:4) = (/ 1.2d0, 1.7d0 /)
  x(5) = 1.0d0
  
  call cpu_time(start)
  do j = 1,100000
    call evalccovsec(x,n,nballs,covarea,c_loc(pdata))
  end do
  call cpu_time(finish)
  write(*,*) 'covered area = ',covarea
  write(*,*) 'CPU time  = ',(finish-start)/100000
  
  do i = 1,5
    h = hvec(i)
    write(*,*) 'h = ',h
    if ( i .eq. 1 .or. i .eq. 2 ) then
        call cpu_time(start)
        do j = 1,10000
            call evalccovfirst(n,x,covarea,h)
        end do
        call cpu_time(finish)
        write(*,*) 'covered area = ',covarea
        write(*,*) 'CPU time = ',(finish-start)/10000
    else
        call cpu_time(start)
        call evalccovfirst(n,x,covarea,h)
        call cpu_time(finish)
        write(*,*) 'covered area = ',covarea
        write(*,*) 'CPU time = ',finish-start
    end if
  end do
  
  deallocate(pdata%nvpols, pdata%poledges, pdata%xpol, pdata%ypol, stat=allocerr)
  
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if
  
  stop
  
contains

  subroutine evalccovsec(x,n,nballs,covarea,pdataptr)
  
    use VorCells_Polygons

    implicit none
    
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nballs
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: covarea

    ! LOCAL SCALARS
    logical :: inside
    integer :: allocerr,i,indpol,ierror,k,status,v,v_num
    real(kind=8) :: area,convPolArea,r,theta1,theta2,totalCoveredArea
    type(pdata_type), pointer :: pdata
    
    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol
    
    ! LOCAL ARRAYS
    integer :: indexes(nballs),nod_tri(3,2*nballs),tnbr(3,2*nballs)
    real(kind=8) :: centers(2,nballs)
  
    ! ALLOCATABLE ARRAYS
    real(kind=8), dimension (:,:), allocatable :: vor_xy
  
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d
    
    call c_f_pointer(pdataptr,pdata)
    
    do i = 1,nballs
       indexes(i) = i
       centers(1,i) = x(2*i-1)
       centers(2,i) = x(2*i)
    end do
    
    ! Radius
    
    r = x(2*nballs+1)
    
    ! Delaunay triangulation
    
    if ( nballs .ge. 3 ) then
      
      call dtris2(nballs, centers, indexes, v_num, nod_tri, tnbr, ierror)
      
      if ( ierror .eq. 0 ) then
        ! Obtaining the Voronoi vertices.
        allocate ( vor_xy(2,v_num), stat=allocerr )
        
        if ( allocerr .ne. 0 ) then
            write(*,*) 'Allocation error.'
            stop
        end if
        
        do k=1,v_num
            call triangle_circumcenter_2d ( centers(1:2, nod_tri(1:3, k)), vor_xy(1:2, k) )
        end do
      
      end if
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Error in dtris2'
        write ( *, '(a,i6)' ) '  IERROR = ', ierror
        return
      end if
      
    end if
    
    totalCoveredArea = 0.0d0
    
    do i = 1,nballs
    
        if ( nballs .ge. 3 .and. ierror .eq. 0 ) then
        
            call voronoi_cell(i, nballs, centers, v_num, nod_tri, tnbr, vor_xy, vorCell, status)
            
            if ( status .ne. 0 ) then
                write(*,*) 'Some problem with subroutine voronoi_cell.'
                return
            end if
            
        end if
        
        do indpol=1,pdata%npols
            
            area = 0.0d0
            convPol%n = pdata%nvpols(indpol)
            convPol%deg = .false.
            
            allocate ( convPol%vertex(2,convPol%n), convPol%edges(convPol%n - 1), stat=allocerr )
            
            if ( allocerr .ne. 0 ) then
               write(*,*) 'Allocation error.'
               stop
            end if
            
            convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
            convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
            convPol%edges(1:convPol%n - 1) = pdata%poledges(indpol, 1:convPol%n - 1)
            
            if ( nballs .ge. 3 .and. ierror .eq. 0 ) then
                call VorCellInterConvPol(centers(:,i), convPol, vorCell, curCell)
            end if
            
            if ( nballs .eq. 2 .or. (nballs .ge. 3 .and. ierror .eq. 225) ) then
                call VorCellInterConvPol_Collinear(i, nballs, centers, convPol, curCell)
            end if
            
            if ( (.not. curCell%deg) .or. nballs .eq. 1 ) then
                
                if ( nballs .eq. 1 ) then
                    call interCellBall(centers(1:2,i), r, convPol, curvPol)
                else
                    call interCellBall(centers(1:2,i), r, curCell, curvPol)
                end if
                
                if ( curvPol%deg ) then
                    
                    if ( nballs .eq. 1 ) then
                        call polygon_contains_point_2d_convex(convPol%n - 1,convPol%vertex(1:2,1:convPol%n - 1),centers(1:2,i),&
                        inside)
                    else
                        call polygon_contains_point_2d_convex(curCell%n - 1,curCell%vertex(1:2,1:curCell%n - 1),centers(1:2,i),&
                        inside)
                    end if
                    
                    if ( inside ) then
                        area = area + (PI * r**2)
                    end if
                else
                    do v = 1,curvPol%n - 1
                        
                        if ( curvPol%edges(v) .eq. 0 ) then
                            area = area + 0.5d0 * ( curvPol%vertex(2,v+1) - curvPol%vertex(2,v) )&
                            * ( curvPol%vertex(1,v) + curvPol%vertex(1,v+1) )
                        else
                            theta1 = angle_rad_2d(curvPol%vertex(1:2,v), centers(1:2,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /))
                            
                            theta2 = angle_rad_2d(curvPol%vertex(1:2,v+1), centers(1:2,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /))
                            
                            if ( theta2 .lt. theta1) then
                                theta2 = theta2 + 2.0d0 * PI
                            end if
                            
                            area = area + 0.5d0 * r**2 * (  (theta2 - theta1) + &
                            (dsin(theta2) * dcos(theta2)) - (dsin(theta1) * dcos(theta1)) ) + centers(1,i) * r * ( &
                            dsin(theta2) - dsin(theta1) )
                        end if
                        
                    end do
                end if
            end if
            
            totalCoveredArea = totalCoveredArea + area
            
            deallocate(convPol%vertex,convPol%edges, stat=allocerr)
            
            if ( allocerr .ne. 0 ) then
               write(*,*) 'Deallocation error.'
               stop
            end if
            
            if ( allocated(curCell%vertex) ) deallocate(curCell%vertex)
            if ( allocated(curCell%edges) ) deallocate(curCell%edges)
            if ( allocated(curvPol%vertex) ) deallocate(curvPol%vertex)
            if ( allocated(curvPol%vertex_id) ) deallocate(curvPol%vertex_id)
            if ( allocated(curvPol%edges) ) deallocate(curvPol%edges)
            
        end do
        
        if ( allocated(vorCell%vertex) ) deallocate(vorCell%vertex)
        if ( allocated(vorCell%edges) ) deallocate(vorCell%edges)
        
    end do
    
    if ( allocated(vor_xy) ) deallocate(vor_xy)

    covarea = totalCoveredArea
    
  end subroutine evalccovsec
  
  !*********************************
  !*********************************
  subroutine evalccovfirst(n,x,covarea,h)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: n
    real(kind=8), intent(in) :: h
    real(kind=8), intent(out) :: covarea

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)

    ! COMMON SCALARS
    integer :: nballs

    ! COMMON ARRAYS
    real(kind=8) :: dbl(2),dtr(2)

    ! COMMON BLOCKS
    common /pdatafirst/ dbl,dtr,nballs

    ! LOCAL SCALARS
    logical :: inA,insomeball
    integer :: i,j,k,nx,ny
    real(kind=8) :: hx,hy,count

    ! LOCAL ARRAYS
    real(kind=8) :: p(2)
    
    nx = ceiling( ( dtr(1) - dbl(1) ) / h )
    ny = ceiling( ( dtr(2) - dbl(2) ) / h )
        
    hx = ( dtr(1) - dbl(1) ) / nx
    hy = ( dtr(2) - dbl(2) ) / ny
        
    count = 0.0d0
    do i = 1,nx
        p(1) = dbl(1) + ( i - 0.5d0 ) * hx
        do j = 1,ny
            p(2) = dbl(2) + ( j - 0.5d0 ) * hy

            inA = .true.

            if ( inA ) then
                k = 0
                insomeball = .false.
                do while ( .not. insomeball .and. k .lt. nballs )
                    k = k + 1
                    if ( ( x(2*k-1) - p(1) ) ** 2 + ( x(2*k) - p(2) ) ** 2 .le. x(2*nballs+1) ** 2 ) then
                        insomeball = .true.
                    end if
                end do
                if ( insomeball ) then
                    count = count + 1.0d0
                end if
            end if
        end do
    end do
        
    covarea = hx * hy * count
    
  end subroutine evalccovfirst

end program squarecoveredarea
