module model_spin_boson
    implicit none

    ! constants
    real*8, parameter :: pi     = dacos(-1.d0)
    complex*16,parameter:: iota = (0.d0, 1.d0)
    real*8, parameter :: hbar   = 1.d0
    real*8, parameter :: kb     = 1.d0
    real*8, parameter :: clight = 2.99792458d10
    real*8, parameter :: mp     = 1836.d0

    ! unit conversions
    real*8, parameter :: sec2au = 4.13414d16
    real*8, parameter :: omg2au = clight*2*pi/sec2au
    real*8, parameter :: eng2au = omg2au    !4.55634d-6 
    real*8, parameter :: tk2au  = 3.16681d-6
    
    ! model parameters
    real*8, parameter :: en     = -120.d0*eng2au
    real*8, parameter :: en1    = en/2  !0.d0*eng2au      !410.d0*eng2au
    real*8, parameter :: en2    = -en/2 !120.d0*eng2au    !530.d0*eng2au   
    real*8, parameter :: vcoup  = 87.7d0*eng2au
    real*8, parameter :: omega  = 1060.d0*omg2au
    real*8, parameter :: eta    = 10.d0*omega
    real*8, parameter :: gamma_B= eta 
    real*8, parameter :: g      = 2.61d-3       ! a.u.
    real*8, parameter :: tempK  = 77.d0*tk2au
    real*8, parameter :: beta   = 1.d0/(kb*tempK)
    

    ! model dimensions
    integer, parameter :: nquant = 2
    integer, parameter :: nclass = 1
    integer, parameter :: nbasis = 2

    real*8, parameter, dimension(nbasis)  :: mass = mp
    real*8, dimension(nbasis, nbasis) :: H, H_diab, H_new   !, H_adiab
    real*8, dimension(nbasis, nbasis, nclass) :: dV_dq      !, dij

    contains

    function ham_diab(q)    result(H_diab)
        implicit none
        real*8, intent(in) :: q(:)
        real*8  H_diab(nbasis,nbasis)

        H_diab  = 0.d0
        H_diab(1,1) = en1 + (0.5* mass(1) * omega**2 * q(1)**2) + g*q(1)  ! sum{mass,q : nclass}
        H_diab(2,2) = en2 + (0.5* mass(1) * omega**2 * q(1)**2) - g*q(1)  ! sum{mass,q : nclass}
        H_diab(1,2) = vcoup
        H_diab(2,1) = H_diab(1,2)
    end function ham_diab

    function delV_delq(q)   result(dV_dq)
        implicit none
        real*8, intent(in) :: q(:)

        real*8  dV_dq(nbasis,nbasis,nclass)
        
        dV_dq   = 0.d0
        dV_dq (1,1,1) = (mass(1) * omega**2 * q(1)) + g    ! sum{mass,q : nclass}
        dV_dq (2,2,1) = (mass(1) * omega**2 * q(1)) - g    ! sum{mass,q : nclass}
    end function delV_delq

    !function d_ij(ci, dV_dq, E) result(dij)
    !     implicit none
    !     real*8, intent(in) :: ci, dV_dq, E, dij(nbasis,nbasis,nclass)
        
    !     d_ij = dij 
    ! end function d_ij

    ! vdotd
    ! ham_adiab

    !-------------------------------------------------------------------------!
    subroutine pes(ipes)
        implicit none
        integer, intent(in) :: ipes
        integer i 
        real*8, dimension(nclass) :: iarr, x
        real*8, dimension(nbasis,nbasis) :: y1
        real*8, dimension(nbasis,nbasis, nclass) :: y2

        iarr = 1.d0
        ! 1-D PES : nclass = 1
        do i = 1, 1000
            x= (-0.25 + 0.5*i/1000)*iarr
            if (ipes==1) then 
                y1 = ham_diab(x)
                write(200,*) x, y1(1,1), y1(1,2), y1(2,1), y1(2,2)
            endif
            if (ipes==2) then 
                y2 = delV_delq(x)
                write(200,*) x, y2(1,1,1), y2(1,2,1), y2(2,1,1), y2(2,2,1)
            endif
        enddo
    end subroutine pes
    !-------------------------------------------------------------------------!
end module model_spin_boson
