module modcesaro

  PUBLIC :: includes_triang_segs
  
  contains
  
    ! ==============================================
    ! ==============================================
    recursive subroutine includes_triang_segs(vxb,vyb,nsegs,iter)

        implicit none
        
        ! PARAMETERS
        real(kind=8), parameter :: PI = dacos( - 1.0d0 )
        
        ! SCALAR ARGUMENTS
        integer, intent(in) :: iter,nsegs

        ! ARRAY ARGUMENTS
        real(kind=8), intent(inout) :: vxb(64),vyb(64)

        ! LOCAL SCALARS
        integer :: i,ind
        real(kind=8) :: rotv(2)

        ! LOCAL ARRAYS
        real(kind=8) :: mthirdx(2),mthirdy(2),auxvxb(64),auxvyb(64)
        
        auxvxb(1:64) = vxb(1:64)
        auxvyb(1:64) = vyb(1:64)
        
        ind = 1
        do i = 1,nsegs
            if ( i .lt. nsegs ) then
                mthirdx(1) = auxvxb(i) + (1.0d0/3.0d0)*(auxvxb(i+1) - auxvxb(i))
                mthirdx(2) = auxvxb(i) + (2.0d0/3.0d0)*(auxvxb(i+1) - auxvxb(i))
                mthirdy(1) = auxvyb(i) + (1.0d0/3.0d0)*(auxvyb(i+1) - auxvyb(i))
                mthirdy(2) = auxvyb(i) + (2.0d0/3.0d0)*(auxvyb(i+1) - auxvyb(i))
                
                call vector_rotate_base_2d( (/ mthirdx(2), mthirdy(2) /), (/ mthirdx(1), mthirdy(1) /), PI / 3.0d0, rotv )
                
            else
                mthirdx(1) = auxvxb(i) + (1.0d0/3.0d0)*(auxvxb(1) - auxvxb(i))
                mthirdx(2) = auxvxb(i) + (2.0d0/3.0d0)*(auxvxb(1) - auxvxb(i))
                mthirdy(1) = auxvyb(i) + (1.0d0/3.0d0)*(auxvyb(1) - auxvyb(i))
                mthirdy(2) = auxvyb(i) + (2.0d0/3.0d0)*(auxvyb(1) - auxvyb(i))
                
                call vector_rotate_base_2d( (/ mthirdx(2), mthirdy(2) /), (/ mthirdx(1), mthirdy(1) /), PI / 3.0d0, rotv )
            end if
            
            vxb(ind) = auxvxb(i)
            vyb(ind) = auxvyb(i)
            ind = ind + 1
            vxb(ind) = mthirdx(1)
            vyb(ind) = mthirdy(1)
            ind = ind + 1
            vxb(ind) = rotv(1)
            vyb(ind) = rotv(2)
            ind = ind + 1
            vxb(ind) = mthirdx(2)
            vyb(ind) = mthirdy(2)
            
            ind = ind + 1
        end do
        
        if ( iter .lt. 2 ) then
        call includes_triang_segs(vxb, vyb, nsegs * 4, iter + 1)
        end if
        
    end subroutine includes_triang_segs

end module modcesaro
