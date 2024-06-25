! *****************************************************************
! *****************************************************************

program algencama

  use bmalgencan, only: algencan
  use iso_c_binding, only: c_ptr, c_loc,c_f_pointer
  use VorCells_Polygons

  implicit none
  
  ! PROBLEM DATA TYPE
  type :: pdata_type
     real(kind=8) :: dbl(2),dtr(2)
     integer :: counters(5)
     integer :: npols,maxnvpols
     integer, dimension(:), allocatable :: nvpols
     integer, dimension(:,:), allocatable :: poledges
     real(kind=8), dimension(:,:), allocatable :: xpol,ypol
  end type pdata_type

  ! LOCAL SCALARS
  logical :: corrin,extallowed,inside,rhoauto,scale
  integer :: allocerr,hlnnzmax,i,ierr,indpol,inform,istop,itrial,jnnzmax,k,m,maxoutit,n,nwcalls,&
             nwtotit,outiter,outside,p,totiter
  real(kind=8) :: bestrad,bdsvio,csupn,epsfeas,epscompl,epsopt,f,finish,nlpsupn,rhoini,seed,ssupn,start
  type(pdata_type), target :: pdata
  
  ! LOCAL ARRAYS
  logical, allocatable :: lind(:),uind(:)
  real(kind=8) :: rp(2)
  real(kind=8), allocatable :: c(:),lbnd(:),ubnd(:),lambda(:),x(:)
  real(kind=8), dimension(:,:), allocatable :: polpart
  
  ! FUNCTIONS
  real(kind=8) :: drand
  
  ! COMMON SCALARS
  integer :: nballs
  
  ! COMMON BLOCKS
  common /data/ nballs

  ! Problem data:
  ! Convex polygons
  pdata%npols = 14
  pdata%maxnvpols = 7
  
  allocate( pdata%nvpols(pdata%npols), pdata%poledges(pdata%npols, pdata%maxnvpols), &
            pdata%xpol(pdata%npols, pdata%maxnvpols), pdata%ypol(pdata%npols, pdata%maxnvpols), stat=allocerr )
  
  ! Number of vertices for each polygon
  pdata%nvpols(1:pdata%npols) = (/ 5, 6, 5, 5, 5, 5, 6, 5, 6, 5, 6, 7, 5, 4 /)
  
  ! Vertices
  pdata%xpol(1, 1:pdata%nvpols(1)) = (/ 0.0d0, 0.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.0d0 /)
  pdata%ypol(1, 1:pdata%nvpols(1)) = (/ 2.0d0/3.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0 /)
  
  pdata%xpol(2, 1:pdata%nvpols(2)) = (/ 35.0d0/150.0d0, 2.0d0/15.0d0, 0.0d0, 0.0d0, 2.0d0/15.0d0, 35.0d0/150.0d0 /)
  pdata%ypol(2, 1:pdata%nvpols(2)) = (/ 1.0d0/3.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 0.2d0, 0.2d0, 1.0d0/3.0d0 /)
  
  pdata%xpol(3, 1:pdata%nvpols(3)) = (/ 0.0d0, 0.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.0d0 /)
  pdata%ypol(3, 1:pdata%nvpols(3)) = (/ 0.2d0, -4.0d0/15.0d0, -4.0d0/15.0d0, 0.2d0, 0.2d0 /)
  
  pdata%xpol(4, 1:pdata%nvpols(4)) = (/ 0.0d0, 0.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.0d0 /)
  pdata%ypol(4, 1:pdata%nvpols(4)) = (/ -4.0d0/15.0d0, -1.0d0/3.0d0, -1.0d0/3.0d0, -4.0d0/15.0d0, -4.0d0/15.0d0 /)
  
  pdata%xpol(5, 1:pdata%nvpols(5)) = (/ 2.0d0/15.0d0, 2.0d0/15.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/15.0d0 /)
  pdata%ypol(5, 1:pdata%nvpols(5)) = (/ 2.0d0/3.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0 /)
  
  pdata%xpol(6, 1:pdata%nvpols(6)) = (/ 7.0d0/15.0d0, 8.0d0/15.0d0, 0.6d0, 7.0d0/15.0d0, 7.0d0/15.0d0 /)
  pdata%ypol(6, 1:pdata%nvpols(6)) = (/ 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0, 7.0d0/15.0d0 /)
  
  pdata%xpol(7, 1:pdata%nvpols(7)) = (/ 8.0d0/15.0d0, 7.0d0/15.0d0, 0.3d0, 7.0d0/15.0d0, 0.6d0, 8.0d0/15.0d0 /)
  pdata%ypol(7, 1:pdata%nvpols(7)) = (/ 7.0d0/15.0d0, 7.0d0/15.0d0, 1.0d0/3.0d0, 0.2d0, 1.0d0/3.0d0, 7.0d0/15.0d0 /)
  
  pdata%xpol(8, 1:pdata%nvpols(8)) = (/ 7.0d0/15.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.4d0, 7.0d0/15.0d0 /)
  pdata%ypol(8, 1:pdata%nvpols(8)) = (/ 0.2d0, 0.2d0, -4.0d0/15.0d0, 0.0d0, 0.2d0 /)
  
  pdata%xpol(9, 1:pdata%nvpols(9)) = (/ 11.0d0/15.0d0, 0.6d0, 7.0d0/15.0d0, 0.4d0, 0.6d0, 11.0d0/15.0d0 /)
  pdata%ypol(9, 1:pdata%nvpols(9)) = (/ 2.0d0/15.0d0, 1.0d0/3.0d0, 0.2d0, 0.0d0, 0.0d0, 2.0d0/15.0d0 /)
  
  pdata%xpol(10, 1:pdata%nvpols(10)) = (/ 13.0d0/15.0d0, 13.0d0/15.0d0, 1.0d0, 1.0d0, 13.0d0/15.0d0 /)
  pdata%ypol(10, 1:pdata%nvpols(10)) = (/ -4.0d0/15.0d0, -1.0d0/3.0d0, -1.0d0/3.0d0, -4.0d0/15.0d0, -4.0d0/15.0d0 /)
  
  pdata%xpol(11, 1:pdata%nvpols(11)) = (/ 13.0d0/15.0d0, 11.0d0/15.0d0, 13.0d0/15.0d0, 1.0d0, 1.0d0, 13.0d0/15.0d0 /)
  pdata%ypol(11, 1:pdata%nvpols(11)) = (/ 1.0d0/3.0d0, 2.0d0/15.0d0, -4.0d0/15.0d0, -4.0d0/15.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0 /)
  
  pdata%xpol(12, 1:pdata%nvpols(12)) = (/ 13.0d0/15.0d0, 0.8d0, 11.0d0/15.0d0, 13.0d0/15.0d0, 1.0d0, 1.0d0, 13.0d0/15.0d0 /)
  pdata%ypol(12, 1:pdata%nvpols(12)) = (/ 2.0d0/3.0d0, 2.0d0/3.0d0, 8.0d0/15.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 8.0d0/15.0d0, &
  2.0d0/3.0d0 /)
  
  pdata%xpol(13, 1:pdata%nvpols(13)) = (/ 11.0d0/15.0d0, 0.8d0, 8.0d0/15.0d0, 0.6d0, 11.0d0/15.0d0 /)
  pdata%ypol(13, 1:pdata%nvpols(13)) = (/ 8.0d0/15.0d0, 2.0d0/3.0d0, 7.0d0/15.0d0, 1.0d0/3.0d0, 8.0d0/15.0d0 /)
  
  pdata%xpol(14, 1:pdata%nvpols(14)) = (/ 11.0d0/15.0d0, 0.6d0, 13.0d0/15.0d0, 11.0d0/15.0d0 /)
  pdata%ypol(14, 1:pdata%nvpols(14)) = (/ 2.0d0/15.0d0, 0.0d0, -4.0d0/15.0d0, 2.0d0/15.0d0 /)
  
  ! Edges ID
  pdata%poledges(1, 1:pdata%nvpols(1) - 1) = (/ 0, 1, 1, 0 /)
  pdata%poledges(2, 1:pdata%nvpols(2) - 1) = (/ 0, 1, 0, 1, 0 /)
  pdata%poledges(3, 1:pdata%nvpols(3) - 1) = (/ 0, 1, 1, 1 /)
  pdata%poledges(4, 1:pdata%nvpols(4) - 1) = (/ 0, 0, 0, 1 /)
  pdata%poledges(5, 1:pdata%nvpols(5) - 1) = (/ 1, 0, 1, 0 /)
  pdata%poledges(6, 1:pdata%nvpols(6) - 1) = (/ 1, 0, 0, 1 /)
  pdata%poledges(7, 1:pdata%nvpols(7) - 1) = (/ 1, 0, 0, 1, 1 /)
  pdata%poledges(8, 1:pdata%nvpols(8) - 1) = (/ 0, 1, 0, 1 /)
  pdata%poledges(9, 1:pdata%nvpols(9) - 1) = (/ 0, 1, 1, 0, 1 /)
  pdata%poledges(10, 1:pdata%nvpols(10) - 1) = (/ 0, 0, 0, 1 /)
  pdata%poledges(11, 1:pdata%nvpols(11) - 1) = (/ 0, 1, 1, 0, 1 /)
  pdata%poledges(12, 1:pdata%nvpols(12) - 1) = (/ 0, 1, 0, 1, 0, 0 /)
  pdata%poledges(13, 1:pdata%nvpols(13) - 1) = (/ 1, 0, 1, 0 /)
  pdata%poledges(14, 1:pdata%nvpols(14) - 1) = (/ 1, 0, 1 /)
  
  ! Rectangle containing A
  pdata%dbl(1:2) = (/ -1.0d0/6.0d0, -0.5d0 /)
  pdata%dtr(1:2) = (/ 7.0d0/6.0d0, 5.0d0/6.0d0 /)
  
  ! Constraints
  
  m = 1
  p = 0

  allocate(lambda(m+p),c(m+p),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if
  
  ! Parameters setting
  
  epsfeas  = 1.0d-08
  epscompl = 1.0d-08
  epsopt   = 1.0d-08

  maxoutit = 50

  rhoauto = .true.

  if ( .not. rhoauto ) then
      rhoini = 1.0d-08
  end if

  scale = .false.

  extallowed = .true.

  corrin = .true.
  
  ! Reading parameters
  write(*,*) 'Enter nballs, itrial and bestrad: '
  read(*,*) nballs,itrial,bestrad
    
  ! Number of variables (x = (c_1,...,c_nballs,r))

  n = 2 * nballs + 1
    
  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),lind(n),lbnd(n),uind(n),ubnd(n),stat=allocerr)

  if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
  end if
    
  lind(1:2*nballs) = .false.
  uind(1:2*nballs) = .false.
  lind(n) = .true.
  lbnd(n) = 0.0d0
  uind(n) = .true.
  ubnd(n) = sqrt(0.89)

  ! Upper bounds on the number of sparse-matrices non-null elements
  
  jnnzmax = n
  hlnnzmax = huge( 1 )
    
  ! Initialize counters
        
  pdata%counters(1:5) = 0
        
  ! Initial guess and Lagrange multipliers
        
  seed = 654321.0d0 * itrial
        
  do i = 1,nballs
10  continue

    rp(1) = pdata%dbl(1) + ( pdata%dtr(1) - pdata%dbl(1) ) * drand(seed)
    rp(2) = pdata%dbl(2) + ( pdata%dtr(2) - pdata%dbl(2) ) * drand(seed)
        
    outside = 0
    do indpol=1,pdata%npols
            
        allocate( polpart(2, pdata%nvpols(indpol)-1), stat=allocerr )
        
        if ( allocerr .ne. 0 ) then
            write(*,*) 'Allocation error.'
            stop
        end if
            
        polpart(1,1:pdata%nvpols(indpol)-1) = pdata%xpol(indpol,1:pdata%nvpols(indpol)-1)
        polpart(2,1:pdata%nvpols(indpol)-1) = pdata%ypol(indpol,1:pdata%nvpols(indpol)-1)
        call polygon_contains_point_2d_convex(pdata%nvpols(indpol)-1,polpart,rp,inside)
            
        deallocate( polpart, stat=allocerr )
            
        if ( allocerr .ne. 0 ) then
            write(*,*) 'Deallocation error.'
            stop
        end if
            
        if ( inside ) then
            exit
        else
            outside = outside + 1
        end if
            
    end do
    if ( outside .eq. pdata%npols ) then
        go to 10
    end if
    x(2*i-1) = rp(1)
    x(2*i)   = rp(2)
  end do
    
  x(2*nballs+1) = ( 0.5d0 + 1.0d0 * drand(seed) ) / nballs

  lambda(1:m+p) = 0.0d0
    
  ! Optimize
    
  call cpu_time(start)
  call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
      n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
      scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
      outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))
  call cpu_time(finish)
  
  if ( ierr .eq. 0  ) then
    if ( csupn .le. epsfeas .and. f .lt. bestrad ) then
      open(10,file='bestrad.txt')
      write(10,*) f
      close(10)
      call drawsol(n,x,c_loc(pdata))
      open(20,file='tabline.txt')
      write(20,9000) n,m+p,f,csupn,ssupn,nlpsupn,bdsvio,outiter,totiter,&
      nwcalls,nwtotit,pdata%counters(1),pdata%counters(2),pdata%counters(3),pdata%counters(4),&
      pdata%counters(5),finish - start,itrial,istop,ierr
      close(20)
    end if
  end if
  
  deallocate(lind,lbnd,uind,ubnd,x,stat=allocerr)

  if ( allocerr .ne. 0 ) then
      write(*,*) 'Deallocation error.'
      stop
  end if
  
  deallocate(lambda,c,stat=allocerr )
  
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if
  
  deallocate(pdata%nvpols, pdata%poledges, pdata%xpol, pdata%ypol, stat=allocerr)
  
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if
  
  stop

9000 format(2(1X,I6),1X,1P,D24.16,4(1X,1P,D7.1),4(1X,I8),5(1X,I7),1X,D8.2,3(1X,I4))

contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf(n,x,f,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    
    ! COMMON SCALARS
    integer :: nballs
  
    ! COMMON BLOCKS
    common /data/ nballs
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(1) = pdata%counters(1) + 1
    
    f = x(n)
    
  end subroutine evalf

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg(n,x,g,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr
  
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    
    ! COMMON SCALARS
    integer :: nballs
  
    ! COMMON BLOCKS
    common /data/ nballs
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1
    
    g(1:n-1) = 0.0d0
    g(n)     = 1.0d0

  end subroutine evalg

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalc(n,x,m,p,c,inform,pdataptr)
  
    use VorCells_Polygons

    implicit none
    
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m+p)
    
    ! COMMON SCALARS
    integer :: nballs
  
    ! COMMON BLOCKS
    common /data/ nballs

    ! LOCAL SCALARS
    logical :: inside
    integer :: allocerr,i,indpol,ierror,k,status,v,v_num
    real(kind=8) :: area,convPolArea,dist,eps,r,seed,theta1,theta2,totalCoveredArea,totalPolsArea
    type(pdata_type), pointer :: pdata
    
    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol
    
    ! LOCAL ARRAYS
    integer :: indexes(nballs),nod_tri(3,2*nballs),tnbr(3,2*nballs)
    real(kind=8) :: centers(2,nballs),pD(2,4)
  
    ! ALLOCATABLE ARRAYS
    real(kind=8), dimension (:,:), allocatable :: vor_xy
  
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(3) = pdata%counters(3) + 1
    
    eps = 1.0d-10
    seed = 654321.0
    
    ! Indexes of the balls and centers
    
    do i = 1,nballs
       indexes(i) = i
       centers(1,i) = x(2*i-1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
       centers(2,i) = x(2*i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
    end do
  
    ! Radius
    
    r = x(2*nballs+1)
    
    ! Points of D
    
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)
    
    ! Check whether r is positive
    if ( r .le. 0.0d0 ) then
        c(1) = 1.0d+09
        write(*,*) '*** IN EVALC SUBROUTINE: RADIUS LESS OR EQUAL THAN ZERO. ***'
        return
    else
        ! Check whether centers are close enough to A
        do i = 1,nballs
            call quad_point_dist_signed_2d( pD, centers(:,i), dist )
            if ( dist .gt. r ) then
                c(1) = 1.0d+09
                write(*,*) '*** IN EVALC SUBROUTINE: THERE IS A CENTER TOO FAR FROM D. ***'
                return
            end if
        end do
    end if
    
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
                inform = -51
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

    totalPolsArea = 0.0d0
    
    do indpol=1,pdata%npols
        
        convPol%n = pdata%nvpols(indpol)
        convPol%deg = .false.
        
        allocate ( convPol%vertex(2,convPol%n), stat=allocerr )
        
        if ( allocerr .ne. 0 ) then
            write(*,*) 'Allocation error.'
            stop
        end if
        
        convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
        convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
        
        call polygon_area_2d ( convPol%n - 1, convPol%vertex(1:2,1:convPol%n - 1), convPolArea )
        
        totalPolsArea = totalPolsArea + convPolArea
        
        deallocate(convPol%vertex, stat=allocerr)
        
        if ( allocerr .ne. 0 ) then
            write(*,*) 'Deallocation error.'
            stop
        end if
        
    end do
    
    c(1) = totalPolsArea - totalCoveredArea
    inform = 0
    
  end subroutine evalc

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)
  
    use VorCells_Polygons

    implicit none
    
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m+p)
    real(kind=8), intent(in) :: x(n)
    logical, intent(out) :: sorted(m+p)
    integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
    real(kind=8), intent(out) :: jval(lim)
    
    ! COMMON SCALARS
    integer :: nballs
  
    ! COMMON BLOCKS
    common /data/ nballs
    
    ! LOCAL SCALARS
    integer :: i,ierror,indpol,k,status,v,v_num
    real(kind=8) :: eps,r,seed,theta1,theta2,dist
    logical :: inside
    type(pdata_type), pointer :: pdata

    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol

    ! LOCAL ARRAYS
    integer :: indexes(nballs),nod_tri(3,2*nballs),tnbr(3, 2*nballs)
    real(kind=8) :: centers(2,nballs),pD(2,4)
    
    ! ALLOCATABLE ARRAYS
    real ( kind = 8 ), dimension (:,:), allocatable :: vor_xy
    
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(4) = pdata%counters(4) + 1
    
    ! Points of D
    
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)

    ! Only gradients of constraints j such that ind(j) = .true. need
    ! to be computed.
    
    if ( ind(1) ) then
       if ( lim .lt. n ) then
          inform = -94
          return
       end if
       
       jsta(1) = 1
       jlen(1) = n
       
       jvar(1:n) = (/ (i,i=1,n) /)

       jval(1:n) = 0.0d0
       
       eps = 1.0d-10
       seed = 654321.0
    
       ! Indexes of the balls and centers
    
       do i = 1,nballs
          indexes(i) = i
          centers(1,i) = x(2*i-1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
          centers(2,i) = x(2*i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
       end do
    
       ! Radius
       r = x(2*nballs+1)
       
       ! Check whether r is positive
       if ( r .le. 0.0d0 ) then
           write(*,*) '*** IN EVALJ SUBROUTINE: RADIUS LESS OR EQUAL THAN ZERO. ***'
           return
       else
           ! Check whether centers are close enough to A
           do i = 1,nballs
               call quad_point_dist_signed_2d( pD, centers(:,i), dist )
               if ( dist .gt. r ) then
                   write(*,*) '*** IN EVALJ SUBROUTINE: THERE IS A CENTER TOO FAR FROM D. ***'
                   return
               end if
           end do
       end if
       
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
       
       do i = 1,nballs
          
          if ( nballs .ge. 3 .and. ierror .eq. 0 ) then
          
             call voronoi_cell(i, nballs, centers, v_num, nod_tri, tnbr, vor_xy, vorCell, status)
             
             if ( status .ne. 0 ) then
                inform = -52
                return
             end if
             
          end if
          
          do indpol=1,pdata%npols
             
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
             
             if ( .not. curCell%deg .or. nballs .eq. 1 ) then
                
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
                      jval(2*nballs+1) = jval(2*nballs+1) + 2.0d0 * PI * r
                   end if
                   
                else
                   
                   do v = 1,curvPol%n - 1
                      if ( curvPol%edges(v) .eq. 1 ) then
                         
                         theta1 = angle_rad_2d(curvPol%vertex(1:2,v), centers(1:2,i), &
                         (/ centers(1,i) + 1.0d0, centers(2,i) /))
                         theta2 = angle_rad_2d(curvPol%vertex(1:2,v+1), centers(1:2,i), &
                         (/ centers(1,i) + 1.0d0, centers(2,i) /))
                         
                         if ( theta2 .lt. theta1) then
                            theta2 = theta2 + 2.0d0 * PI
                         end if
                         
                         jval(2*nballs+1) = jval(2*nballs+1) + r * (theta2 - theta1)
                         jval(2*i-1)      = jval(2*i-1) + r * (dsin(theta1) - dsin(theta2))
                         jval(2*i)        = jval(2*i)   + r * (dcos(theta2) - dcos(theta1))
                         
                      end if
                   end do
                   
                end if
                
             end if
             
             deallocate( convPol%vertex,convPol%edges, stat=allocerr )
             
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
       jval(2*nballs+1) = -1.0d0 * jval(2*nballs+1)
       if ( allocated(vor_xy) ) deallocate(vor_xy)

       ! Says whether the variables' indices in jvar (related to this
       ! constraint) are in increasing order. In case they are,
       ! Algencan takes advantage of this. Implement sorted gradients
       ! of constraints if you can do this in a natural (cheap)
       ! way. Under no circumnstance use a sorting algorithm. (It is
       ! better to set sorted(1) = .false. in this case.)
       
       sorted(1) = .true.
       inform = 0
    end if
    
  end subroutine evalj

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)
    
    use VorCells_Polygons

    implicit none
    
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: inclf
    integer, intent(in) :: m,n,lim,p
    integer, intent(out) :: hlnnz
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: lambda(m+p),x(n)
    integer, intent(out) :: hlrow(lim),hlcol(lim)
    real(kind=8), intent(out) :: hlval(lim)
    
    ! COMMON SCALARS
    integer :: nballs
  
    ! COMMON BLOCKS
    common /data/ nballs
    
    ! LOCAL SCALARS
    integer :: i,ierror,indpol,j,k,status,v,v_num
    real(kind=8) :: dist,eps,inner_prod_div,r,seed,theta1,theta2,theta_set_inter
    logical :: inside
    type(pdata_type), pointer :: pdata

    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol

    ! LOCAL ARRAYS
    integer :: block_diag_ind(3),block_row_vector_ind(2),indexes(nballs),nod_tri(3,2*nballs),tnbr(3, 2*nballs)
    real(kind=8) :: block_diag_val(3),block_row_vector_val(2),centers(2,nballs),normal(2),pbA(2),pD(2,4)
    
    ! ALLOCATABLE ARRAYS
    integer, dimension (:,:), allocatable :: blocks_below_diag_ind
    integer, dimension (:), allocatable :: intersections
    real(kind=8), dimension (:,:), allocatable :: blocks_below_diag_vals,vor_xy
    
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(5) = pdata%counters(5) + 1
    
    ! Points of D
    
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)

    ! If .not. inclf then the Hessian of the objective function must not be included
    ! The Hessian of the objective function is zero.
    
    eps = 1.0d-10
    seed = 654321.0
    
    ! Indexes of the balls and centers
    
    do i = 1,nballs
       indexes(i) = i
       centers(1,i) = x(2*i-1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
       centers(2,i) = x(2*i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
    end do
  
    ! Radius
    r = x(2*nballs+1)
    
    hlnnz = 0
    
    ! Check whether r is positive
    if ( r .le. 0.0d0 ) then
        write(*,*) '*** IN EVALHL SUBROUTINE: RADIUS LESS OR EQUAL THAN ZERO. ***'
        write(*,*) '*** RETURNING THE HESSIAN AS ZERO. ***'
        return
    else
        ! Check whether centers are close enough to A
        do i = 1,nballs
            call quad_point_dist_signed_2d( pD, centers(:,i), dist )
            if ( dist .gt. r ) then
                write(*,*) '*** IN EVALHL SUBROUTINE: THERE IS A CENTER TOO FAR FROM D. ***'
                write(*,*) '*** RETURNING THE HESSIAN AS ZERO. ***'
                return
            end if
        end do
    end if
    
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
    
    ! Note that entries of the Hessian of the Lagrangian can be
    ! repeated. If this is case, them sum of repeated entrances is
    ! considered. This feature simplifies the construction of the
    ! Hessian of the Lagrangian.
    
    if ( hlnnz + 1 .gt. lim ) then
       inform = -95
       return
    end if
    
    ! Second order derivative with respect to the radius.
    hlnnz = hlnnz + 1 
    hlrow(hlnnz) = n
    hlcol(hlnnz) = n
    hlval(hlnnz) = 0.0d0
    
    do i = 1,nballs
        
        ! These arrays are useful for computing blocks below the diagonal.
        if ( i .ne. nballs ) then
            allocate( intersections(nballs - i),blocks_below_diag_ind(4,nballs - i),blocks_below_diag_vals(4,nballs - i), &
            stat=allocerr )
            
            if ( allocerr .ne. 0 ) then
               write(*,*) 'Allocation error.'
               stop
            end if
            
            intersections(1:nballs - i) = 0
            blocks_below_diag_ind(1:4,1:nballs - i) = 0
            blocks_below_diag_vals(1:4,1:nballs - i) = 0.0d0
        end if
        
        if ( hlnnz + 5 .gt. lim ) then
           inform = -95
           return
        end if
        
        ! Diagonal and row-vector blocks.
        
        hlnnz = hlnnz + 1
        block_diag_ind(1) = hlnnz
        hlrow(hlnnz) = 2*i - 1
        hlcol(hlnnz) = 2*i - 1
        
        hlnnz = hlnnz + 1
        block_diag_ind(2) = hlnnz
        hlrow(hlnnz) = 2*i
        hlcol(hlnnz) = 2*i - 1
        
        hlnnz = hlnnz + 1
        block_diag_ind(3) = hlnnz
        hlrow(hlnnz) = 2*i
        hlcol(hlnnz) = 2*i
        
        block_diag_val(1:3) = 0.0d0
        
        hlnnz = hlnnz + 1
        block_row_vector_ind(1) = hlnnz
        hlrow(hlnnz) = n
        hlcol(hlnnz) = 2*i - 1
        
        hlnnz = hlnnz + 1
        block_row_vector_ind(2) = hlnnz
        hlrow(hlnnz) = n
        hlcol(hlnnz) = 2*i
        
        block_row_vector_val(1:2) = 0.0d0
        
        if ( nballs .ge. 3 .and. ierror .eq. 0 ) then
        
            call voronoi_cell(i, nballs, centers, v_num, nod_tri, tnbr, vor_xy, vorCell, status)
            
            if ( status .ne. 0 ) then
                inform = -53
                return
            end if
            
        end if
        
        do indpol=1,pdata%npols
            
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
            
            if ( .not. curCell%deg .or. nballs .eq. 1 ) then
                
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
                        hlval(1) = hlval(1) - 2.0d0 * PI
                    end if
                    
                else
                
                    do v = 1,curvPol%n - 1
                        
                        if ( curvPol%edges(v) .eq. 1 ) then ! Connection is an arc
                            
                            theta1 = angle_rad_2d( curvPol%vertex(:,v), centers(:,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /) )
                            
                            theta2 = angle_rad_2d( curvPol%vertex(:,v+1), centers(:,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /) )
                            
                            if ( theta2 .lt. theta1) then
                                theta2 = theta2 + 2.0d0 * PI
                            end if
                            
                            ! Sec. order derivative with respect to the radius (perimeter).
                            hlval(1) = hlval(1) + theta1 - theta2
                            
                            ! Integrals.
                            
                            ! Diagonal blocks:
                            block_diag_val(1) = block_diag_val(1) + dsin(theta1 - theta2)*dcos(theta1 + theta2)
                            block_diag_val(2) = block_diag_val(2) + dcos(theta2)**2 - dcos(theta1)**2
                            block_diag_val(3) = block_diag_val(3) + dsin(theta2 - theta1)*dcos(theta1 + theta2)

                            ! Two-dimensional row vector:
                            block_row_vector_val(1) = block_row_vector_val(1) + dsin(theta1) - dsin(theta2)
                            block_row_vector_val(2) = block_row_vector_val(2) + dcos(theta2) - dcos(theta1)
                            
                            ! Starting points
                            if ( curvPol%vertex_id(v) .lt. 0 ) then ! intersection with another ball
                                
                                theta_set_inter = angle_rad_2d( curvPol%vertex(:,v), centers(:,iabs(curvPol%vertex_id(v))), &
                                (/ centers(1,iabs(curvPol%vertex_id(v))) + 1.0d0, centers(2,iabs(curvPol%vertex_id(v))) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) - ((dcos(theta_set_inter - theta1) - 1.0d0) / dsin(theta_set_inter - theta1))
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) - (1.0d0 / dtan(theta_set_inter - theta1)) &
                                                    * (dcos(theta1)**2)
                                block_diag_val(2) = block_diag_val(2) - (1.0d0 / dtan(theta_set_inter - theta1)) &
                                                    * dsin(theta1) * dcos(theta1)
                                block_diag_val(3) = block_diag_val(3) - (1.0d0 / dtan(theta_set_inter - theta1)) &
                                                    * (dsin(theta1)**2)
                                
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) &
                                                            - (1.0d0 / dtan(theta_set_inter - theta1)) * dcos(theta1) &
                                                            + (dcos(theta1) / dsin(theta_set_inter - theta1))
                                block_row_vector_val(2) = block_row_vector_val(2) &
                                                            - (1.0d0 / dtan(theta_set_inter - theta1)) * dsin(theta1) &
                                                            + (dsin(theta1) / dsin(theta_set_inter - theta1))
                                
                                ! Derivative with respect to the intersecting balls centers.
                                if ( iabs(curvPol%vertex_id(v)) .gt. i ) then
                                    if ( intersections(iabs(curvPol%vertex_id(v)) - i) .eq. 0 ) then
                                    
                                        intersections(iabs(curvPol%vertex_id(v)) - i) = 1
                                        
                                        if ( hlnnz + 4 .gt. lim ) then
                                            inform = -95
                                            return
                                        end if
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(1,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v)) - 1
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(2,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v))
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(3,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v)) - 1
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(4,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v))
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                    else
                                        
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                    end if
                                end if
                            end if
                            
                            if ( curvPol%vertex_id(v) .eq. 0 ) then ! intersection with the boundary of A
                                
                                ! Computing the normal vector to the boundary of A.
                                if ( v .eq. 1 ) then
                                    pbA(1:2) = curvPol%vertex(1:2,curvPol%n - 1)
                                else
                                    pbA(1:2) = curvPol%vertex(1:2,v - 1)
                                end if
                                call line_exp_normal_2d ( curvPol%vertex(:,v), pbA, normal )
                                
                                inner_prod_div = dot_product( normal, (/ dcos(theta1), dsin(theta1) /) ) / &
                                                dot_product( normal, (/ - dsin(theta1), dcos(theta1) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) - inner_prod_div
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) - inner_prod_div * (dcos(theta1)**2)
                                block_diag_val(2) = block_diag_val(2) - inner_prod_div * dsin(theta1) * dcos(theta1)
                                block_diag_val(3) = block_diag_val(3) - inner_prod_div * (dsin(theta1)**2)
                            
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) - inner_prod_div * dcos(theta1)
                                block_row_vector_val(2) = block_row_vector_val(2) - inner_prod_div * dsin(theta1)
                            end if
                            
                            ! Ending points
                            if ( curvPol%vertex_id(v+1) .lt. 0 ) then ! intersection with another ball
                                
                                theta_set_inter = angle_rad_2d( curvPol%vertex(:,v+1), centers(:,iabs(curvPol%vertex_id(v+1))), &
                                (/ centers(1,iabs(curvPol%vertex_id(v+1))) + 1.0d0, centers(2,iabs(curvPol%vertex_id(v+1))) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) - ((1.0d0 - dcos(theta_set_inter - theta2)) / dsin(theta_set_inter - theta2))
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) + (1.0d0 / dtan(theta_set_inter - theta2)) &
                                                    * (dcos(theta2)**2)
                                block_diag_val(2) = block_diag_val(2) + (1.0d0 / dtan(theta_set_inter - theta2)) &
                                                    * dsin(theta2) * dcos(theta2)
                                block_diag_val(3) = block_diag_val(3) + (1.0d0 / dtan(theta_set_inter - theta2)) &
                                                    * (dsin(theta2)**2)
                                
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) &
                                                            + (1.0d0 / dtan(theta_set_inter - theta2)) * dcos(theta2) &
                                                            - (dcos(theta2)/dsin(theta_set_inter - theta2))
                                block_row_vector_val(2) = block_row_vector_val(2) &
                                                            + (1.0d0 / dtan(theta_set_inter - theta2)) * dsin(theta2) &
                                                            - (dsin(theta2)/dsin(theta_set_inter - theta2))
                                
                                ! Derivative with respect to the intersecting balls centers.
                                if ( iabs(curvPol%vertex_id(v+1)) .gt. i ) then
                                    if ( intersections(iabs(curvPol%vertex_id(v+1)) - i) .eq. 0 ) then
                                    
                                        intersections(iabs(curvPol%vertex_id(v+1)) - i) = 1
                                        
                                        if ( hlnnz + 4 .gt. lim ) then
                                            inform = -95
                                            return
                                        end if
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(1,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1)) - 1
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(2,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1))
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )

                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(3,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1)) - 1
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(4,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1))
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                    else
                                        
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                    end if
                                end if
                            end if
                            
                            if ( curvPol%vertex_id(v+1) .eq. 0 ) then ! intersection with the boundary of A
                                
                                ! Computing the normal vector to the boundary of A.
                                if ( v+1 .eq. curvPol%n ) then
                                    pbA(1:2) = curvPol%vertex(1:2,2)
                                else
                                    pbA(1:2) = curvPol%vertex(1:2,v + 2)
                                end if
                                call line_exp_normal_2d ( pbA, curvPol%vertex(:,v+1), normal )
                                
                                inner_prod_div = dot_product( normal, (/ dcos(theta2), dsin(theta2) /) ) / &
                                                dot_product( normal, (/ - dsin(theta2), dcos(theta2) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) + inner_prod_div
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) + inner_prod_div * (dcos(theta2)**2)
                                block_diag_val(2) = block_diag_val(2) + inner_prod_div * dsin(theta2) * dcos(theta2)
                                block_diag_val(3) = block_diag_val(3) + inner_prod_div * (dsin(theta2)**2)
                                
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) + inner_prod_div * dcos(theta2)
                                block_row_vector_val(2) = block_row_vector_val(2) + inner_prod_div * dsin(theta2)
                            end if
                        end if
                    end do
                end if
            end if
            
            deallocate( convPol%vertex,convPol%edges,stat=allocerr )
            
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
        
        hlval(block_diag_ind(1)) = block_diag_val(1)
        hlval(block_diag_ind(2)) = block_diag_val(2)
        hlval(block_diag_ind(3)) = block_diag_val(3)
        
        do j = 1,nballs - i
            if ( intersections(j) .eq. 1 ) then
                hlval(blocks_below_diag_ind(1,j)) = blocks_below_diag_vals(1,j)
                hlval(blocks_below_diag_ind(2,j)) = blocks_below_diag_vals(2,j)
                hlval(blocks_below_diag_ind(3,j)) = blocks_below_diag_vals(3,j)
                hlval(blocks_below_diag_ind(4,j)) = blocks_below_diag_vals(4,j)
            end if
        end do
        
        hlval(block_row_vector_ind(1)) = block_row_vector_val(1)
        hlval(block_row_vector_ind(2)) = block_row_vector_val(2)
        
        if ( allocated(intersections) ) deallocate(intersections)
        if ( allocated(blocks_below_diag_ind) ) deallocate(blocks_below_diag_ind)
        if ( allocated(blocks_below_diag_vals) ) deallocate(blocks_below_diag_vals)
    
    end do
    if ( allocated(vor_xy) ) deallocate(vor_xy)
    
    hlval(1:hlnnz) = hlval(1:hlnnz) * lambda(1)
    inform = 0
    
  end subroutine evalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine drawsol(n,x,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! COMMON SCALARS
    integer :: nballs
  
    ! COMMON BLOCKS
    common /data/ nballs

    ! LOCAL SCALARS
    integer :: i,j,k
    type(pdata_type), pointer :: pdata

    ! LOCAL ARRAYS
    character(len=80) :: filename
    
    call c_f_pointer(pdataptr,pdata)

    write(filename,"(A15,I0,A3)") 'covering-stoyan',nballs,'.mp'
    
    open(unit=10,file=trim(filename))

    ! BEGINING
    write(10,10) nballs

    ! BALLS FILL
    do i = 1,nballs
        write(10,29) 2.0d0*x(2*nballs+1),x(2*i-1),x(2*i)
    end do

    ! REGION A. 
    ! Define convex polygons as paths
    do j = 1,pdata%npols
        write(10,50) j
        do k = 1,pdata%nvpols(j)
            write(10,100) pdata%xpol(j,k),pdata%ypol(j,k)
        end do
        write(10,80)
    end do
    
    ! Fill the region
    do j = 1,pdata%npols
        write(10,150) j
    end do
    
    ! Draw boundary using the partitions
    do j = 1,pdata%npols
        do k = 1,pdata%nvpols(j)
            if ( k .lt. pdata%nvpols(j) ) then
                if ( pdata%poledges(j,k) .eq. 1 ) then
                    write(10,125) pdata%xpol(j,k),pdata%ypol(j,k),pdata%xpol(j,k+1),pdata%ypol(j,k+1)
                end if
            end if
        end do
    end do
    
    do j = 1,pdata%npols
        do k = 1,pdata%nvpols(j)
            if ( k .lt. pdata%nvpols(j) ) then
                if ( pdata%poledges(j,k) .eq. 0 ) then
                    write(10,130) pdata%xpol(j,k),pdata%ypol(j,k),pdata%xpol(j,k+1),pdata%ypol(j,k+1)
                end if
            end if
        end do
    end do
    
    ! BALLS DRAW
    do i = 1,nballs
        write(10,30) 2.0d0*x(2*nballs+1),x(2*i-1),x(2*i)
    end do

    write(10,40) x(2*nballs+1)

    close(10)

    ! NON-EXECUTABLE STATEMENTS

    10 format('prologues := 3;',/, &
            'outputtemplate := "%j-%c.mps";',/, &
            'input mpcolornames;',/, &
            'beginfig(',i3,');',/, &
            'u = 5cm;',/, &
            'path pols[];')
    50 format('pols[',I0,'] = ')
    80 format('cycle;')
    100 format('(',f20.10,'u,',f20.10,'u)--')
    125 format('draw (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.7 withcolor 0.9Dandelion;')
    130 format('draw (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.7;')
    150 format('fill pols[',I0,'] withcolor 0.9Dandelion;')
    29 format('fill fullcircle scaled ',f20.10,'u shifted (',f20.10,'u,',f20.10,'u) withcolor 0.9white;')
    30 format('draw fullcircle scaled ',f20.10,'u shifted (',f20.10,'u,',f20.10,'u) withcolor RoyalBlue withpen',/, &
    'pencircle scaled 0.7;')
    40 format('% r = ',f20.10,/,'endfig;',/,'end;') 

  end subroutine drawsol
  
end program algencama
