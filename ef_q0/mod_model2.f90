module model_spin_boson_q0
    implicit none

    !! constants
    real*8, parameter :: pi     = dacos(-1.d0)
    complex*16,parameter:: iota = (0.d0, 1.d0)
    real*8, parameter :: hbar   = 1.d0
    real*8, parameter :: kb     = 1.d0
    real*8, parameter :: clight = 2.99792458d10
    real*8, parameter :: mp     = 1836.d0

    !! unit conversions
    real*8, parameter :: sec2au = 4.13414d16
    real*8, parameter :: omg2au = clight*2*pi/sec2au
    real*8, parameter :: eng2au = omg2au    !4.55634d-6 
    real*8, parameter :: tk2au  = 3.16681d-6
    real*8, parameter :: ang2au = 1.8897259886
    
    !! model parameters
    real*8 :: en     = -900.d0*eng2au
    real*8 :: en1    = 0.d0*eng2au     
    real*8 :: en2    = -900.d0*eng2au   
    real*8 :: vcoup  = 150.d0*eng2au
    real*8, parameter :: omg1   = 1200.d0*omg2au
    real*8, parameter :: omg2   = 400.d0*omg2au
    real*8, parameter :: eta    = 400.d0*omg2au
    real*8, parameter :: gamma_B= eta
    real*8 :: lmda1  = 12000.d0*eng2au    
    real*8 :: tempK  = 400.d0*tk2au
    real*8 :: beta   != 1.d0/(kb*tempK)
    

    !! model dimensions
    integer, parameter :: nquant = 2
    integer, parameter :: nclass = 2    ! q1, q2
    integer, parameter :: nbasis_e = 2
    !integer, parameter :: nbasis = nbasis_e 

    real*8, parameter, dimension(nclass)  :: mass = mp
    ! real*8, dimension(nbasis_e, nbasis_e) :: H_el, H_el_diab, H, H_new
    ! real*8, dimension(nbasis_ev, nbasis_ev) :: H_vc, H_vc_diab, H_vc_new   !, H_adiab
    ! real*8, dimension(nbasis_e, nbasis_e, nclass) :: dV_dq      !, dij

    real*8, dimension(nclass), parameter :: omega = (/omg1, omg2/) 
    real*8 :: g1   != dsqrt(0.5*lmda1*mass(1)*omg1**2)      
    real*8 :: g2   != 6.d-3       ! a.u.
    real*8, dimension(nclass) :: g     != (/g1, g2/) 

    contains

    subroutine setup_param()
        implicit none
        logical :: tk_file_exists, vc_file_exists, l1_file_exists, et_file_exists

        ! read simulation parameters
        INQUIRE(FILE="inparam_tk.txt", EXIST=tk_file_exists)    ! Temperature
        INQUIRE(FILE="inparam_vc.txt", EXIST=vc_file_exists)    ! Vcoupling
        INQUIRE(FILE="inparam_l1.txt", EXIST=l1_file_exists)    ! Lambda_1
        INQUIRE(FILE="inparam_et.txt", EXIST=et_file_exists)    ! Exothermicity

        if (tk_file_exists .eqv. .True.) then
             open(1, file="inparam_tk.txt")
             read(1,*) tempK
             tempK = tempK*tk2au
        endif

        if (vc_file_exists .eqv. .True.) then
            open(2, file="inparam_vc.txt")
            read(2,*) vcoup
            vcoup = vcoup*eng2au
        endif

        if (l1_file_exists .eqv. .True.) then
            open(3, file="inparam_l1.txt")
            read(3,*) lmda1
            lmda1 = lmda1*eng2au
        endif

        if (et_file_exists .eqv. .True.) then
            open(4, file="inparam_et.txt")
            read(4,*) en2 !en, en1, en2
            !en = en*eng2au
            !en1= en1*eng2au
            en2= en2*eng2au
        endif

        print *, "tempK: ", tempK
        print *, "vcoup: ", vcoup
        print *, "lmda1: ", lmda1
        print *, "en, en1, en2: ", en, en1, en2

        beta = 1.d0/(kb*tempK)
        g1   = dsqrt(0.5*lmda1*mass(1)*omg1**2)
        g2   = 6.d-3
        g    = (/g1, g2/) 

    end subroutine setup_param


    function ham_diab(q)    result(H_el_diab)
        implicit none
        real*8, intent(in) :: q(nclass)
        real*8  H_el_diab(nbasis_e,nbasis_e)
        real*8 ho_1, ho_2     

        H_el_diab  = 0.d0
        ho_1 = (0.5* mass(1) * omg1**2 * q(1)**2)
        ho_2 = (0.5* mass(2) * omg2**2) * (q(2) - g2*q(1)/(mass(2)*omg2**2))**2
        
        H_el_diab(1,1) = en1 + ho_1 + g1*q(1) + ho_2  
        H_el_diab(2,2) = en2 + ho_1 - g1*q(1) + ho_2
        H_el_diab(1,2) = vcoup
        H_el_diab(2,1) = H_el_diab(1,2)
    end function ham_diab

    function delV_delq(q)   result(dV_dq)
        implicit none
        real*8, intent(in) :: q(nclass)
        real*8  dV_dq(nbasis_e,nbasis_e,nclass)       

        dV_dq   = 0.d0
        dV_dq (1,1,1) = mass(1) * omg1**2 * q(1) + g1 - g2*q(2) + g2**2*q(1)/(mass(2)*omg2**2)
        dV_dq (1,1,2) = mass(2) * omg2**2 *q(2) - g2*q(1)       
        dV_dq (2,2,1) = mass(1) * omg1**2 * q(1) - g1 - g2*q(2) + g2**2*q(1)/(mass(2)*omg2**2)
        dV_dq (2,2,2) = dV_dq (1,1,2)  
    end function delV_delq

    subroutine  ham_adiab(q, nbasis, eval)
        implicit none
        real*8, intent(in)  :: q(:)
        integer, intent(in) :: nbasis
        real*8, intent(out) :: eval(nbasis)
        real*8 evec(nbasis,nbasis)

        evec = ham_diab(q)
        call diasym(evec, eval, nbasis, 'N') ! only eigenvalues calculated     
    
    end subroutine ham_adiab

    !-------------------------------------------------------------------------!

    !-----------------------------------------------------------!
    !Calculates EXPECTATION VALUE of the operator               !
    !input:  ci(n)  = complex coefficients of the wavefunction  !
    !        M(n,n) = Operator matrix                           !
    !output: expval = expectation value of the operator         !
    !NOTE : If the vectors are COMPLEX(*), dot_product(X,Y)     !
    !       returns SUM(CONJG(X)*Y).                            !
    !      -Since ci is complex, it's complex conjugate is taken!
    !      -For consistency, take real part of calculated expval!         
    !-----------------------------------------------------------!
    function expec_val(ci_b, M, ci_k) result(expval)
        implicit none
        complex*16, intent(in)  :: ci_b(:), ci_k(:)
        real*8,     intent(in)  :: M(:,:)
        real*8  expval
    
        expval = real( dot_product(ci_b, matmul(M, ci_k)) ) 
      end function expec_val

    !-------------------!

    !---------------------------------------------------------!
    !Calls the LAPACK diagonalization subroutine DSYEV        !
    !input:  a(n,n) = real symmetric matrix to be diagonalized!
    !            n  = size of a                               !
    !output: a(n,n) = orthonormal eigenvectors of a           !
    !        eig(n) = eigenvalues of a in ascending order     !
    !*add at complile time : -L/usr/local/lib -llapack -lblas !
    !---------------------------------------------------------!
    subroutine diasym(a,eig,n, iev)
    implicit none
    
    character iev 
    integer n,l,inf
    real*8  a(n,n),eig(n),work(n*(3+n/2))
    
    l=n*(3+n/2)
    call dsyev(iev,'U',n,a,n,eig,work,l,inf)
    end subroutine diasym
    !---------------------!

    !-------------------------------------------------------------------------!

    subroutine pes(ipes, wrt_au)
        implicit none
        integer, intent(in) :: ipes, wrt_au
        integer i
        real*8, dimension(nclass) :: iarr, x
        real*8 xlim
        real*8, dimension(nbasis_e,nbasis_e) :: y1
        real*8, dimension(nbasis_e) :: y2
        real*8, dimension(nbasis_e,nbasis_e, nclass) :: y3
        !! without vibrational quantization : electronic hamiltonian
        
        !  1-D PES : along x(1) keeping x(2) at equilibrium
        iarr = (/1.d0,0.d0/)!1.d0
        xlim = 2*ang2au     !! au
        do i = 1, 1000
            x= (-xlim + 2*xlim*i/1000)*iarr
            x(2) = g2*x(1)/(mass(2)*omg1**2) 

            !! diabatic pes
            if (ipes==1) then 
                y1 = ham_diab(x)
                if (wrt_au==0)  then 
                    x   = x/ang2au 
                    y1  = y1/eng2au
                endif
                write(200,*) x, y1(1,1), y1(1,2), y1(2,1), y1(2,2)
            endif
            
            !! adiabatic pes
            if (ipes==2) then
                call ham_adiab(x, nbasis_e, y2)
                if (wrt_au==0)  then 
                    x   = x/ang2au 
                    y2  = y2/eng2au
                endif
                write(201,*) x, y2(1), y2(2)
            endif

            !! derivative of pes
            if (ipes==3) then 
                y3 = delV_delq(x)
                if (wrt_au==0)  then 
                    x   = x/ang2au 
                    y3  = y3/eng2au
                endif
                write(203,*) x, y3(1,1,1), y3(1,2,1), y3(2,1,1), y3(2,2,1)
                write(204,*) x, y3(1,1,2), y3(1,2,2), y3(2,1,2), y3(2,2,2)
            endif
        enddo
    end subroutine pes
    !-------------------------------------------------------------------------!
end module model_spin_boson_q0
