module model_tully_1
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
	real*8 :: A = 0.01
	real*8 :: B = 1.6
	real*8 :: C = 0.005   
	real*8 :: D = 1 
  real*8 :: ngamma = 0.366
  real*8 :: Rmax = 10.0  ! a.u.
  real*8 :: Rmin = -10.0 ! a.u.

	!! model dimensions
	! integer, parameter :: nquant = 2
	integer, parameter :: nclass = 1    ! q1, q2
	integer, parameter :: nbasis_e = 2
	!integer, parameter :: nbasis = nbasis_e

  integer, parameter, dimension(nbasis_e) :: N_s1 = [1,0]
  integer, parameter, dimension(nbasis_e) :: N_s2 = [0,1]
	real*8, parameter  :: mass = 2000.0 !mp !dimension(nclass)  :: mass = mp
	! real*8, dimension(nbasis_e, nbasis_e) :: H_el, H_el_diab, H, H_new
	! real*8, dimension(nbasis_ev, nbasis_ev) :: H_vc, H_vc_diab, H_vc_new   !, H_adiab
	! real*8, dimension(nbasis_e, nbasis_e, nclass) :: dV_dq      !, dij

  ! real*8, dimension(nclass) :: Rn, Pn 
	contains
  
  function ham_eDiab(Rn) result(H_e)
    implicit none

    real*8, intent(in) :: Rn(nclass)
    real*8  H_e(nbasis_e, nbasis_e)
    real*8  R
    ! real*8 A,B,C,D 

    R = Rn(1)
    ! Electronic Hamiltonian
    H_e   = 0.d0
    if (R>=0.0) then
      H_e(1,1) = A*(1.0-exp(-B*R))
    else
      H_e(1,1) = -A*(1.0-exp(B*R))
    endif
    H_e(2,2) = -H_e(1,1)
    H_e(1,2) = C*exp(-D*R**2)
    H_e(2,1) = H_e(1,2) 
  end function ham_eDiab

  function ham_ne(Pn, Rn, pe, xe) result(H_ne)
    implicit none

    real*8, intent(in), dimension(nclass)   :: Pn, Rn
    real*8, intent(in), dimension(nbasis_e) :: pe, xe
    real*8  H_e(nbasis_e, nbasis_e)
    real*8  H_el, H_ne  ! Hamiltonian function
    real*8  P

    P = Pn(1)
    H_e = ham_eDiab(Rn)
    ! Electronic Hamiltonian
    H_el  = 0.5*(H_e(1,1)+H_e(2,2))
    H_el  = H_el+ 0.25*(pe(1)**2+xe(1)**2 - pe(2)**2-xe(2)**2)*(H_e(1,1)-H_e(2,2))
    H_el  = H_el+ (pe(1)*pe(2) + xe(1)*xe(2))*H_e(1,2) 

    ! Nuclear-Electronic Hamiltonian
    H_ne  = P**2/(2*mass) + H_el

    ! xe_dot(1) = 0.5*pe(1)*(H_e(1,1)-H_e(2,2)) + pe(2)*H_e(1,2)
    ! xe_dot(2) = -0.5*pe(2)*(H_e(1,1)-H_e(2,2)) + pe(1)*H_e(1,2)
    !pe_dot(1) = -(0.5*xe(1)*(H_e(1,1)-H_e(2,2)) + xe(2)*H_e(1,2))
    ! pe_dot(2) = -(-0.5*xe(2)*(H_e(1,1)-H_e(2,2)) + xe(1)*H_e(1,2))

  end function ham_ne

  function delHeDiab_delR(Rn)   result(dV_dR)
    implicit none
    
    real*8, intent(in) :: Rn(nclass)
    real*8  dV_dR(nbasis_e,nbasis_e,nclass)
    real*8  R       

    R = Rn(1)
    dV_dR   = 0.d0
    if (R>=0.0) then
      dV_dR (1,1,1) = A*B*exp(-B*R)
    else
      dV_dR (1,1,1) = A*B*exp(B*R)
    endif
    dV_dR (2,2,1) = -dV_dR (1,1,1)
    dV_dR (1,2,1) = -2.0*R*C*D*exp(-D*R**2)
    dV_dR (2,1,1) = dV_dR (1,2,1)
  end function delHeDiab_delR

	subroutine  ham_eAdiab(Rn, nbasis, eval)
		implicit none

		real*8, intent(in)  :: Rn(nclass)
		integer, intent(in) :: nbasis
		real*8, intent(out) :: eval(nbasis)
		real*8 evec(nbasis,nbasis)

		evec = ham_eDiab(Rn)
		call diasym(evec, eval, nbasis, 'N') ! only eigenvalues calculated     

	end subroutine ham_eAdiab

    !-------------------!
  subroutine nq_xp(nk,qk, xe, pe)
    implicit none

    real*8, intent(in), dimension(nbasis_e) :: nk, qk
    real*8, intent(out), dimension(nbasis_e):: xe, pe
    xe = sqrt(2*(nk+ngamma))*cos(qk)
    pe = -sqrt(2*(nk+ngamma))*sin(qk)
  end subroutine nq_xp 

  function ck_nq(nk,qk) result(ck)
    implicit none
    
    real*8, intent(in), dimension(nbasis_e) :: nk, qk
    complex*16  ck(nbasis_e)
    ck = sqrt(nk+ngamma)*exp(-iota*qk)
  end function ck_nq

  function ck_xp(xe,pe) result(ck)
    implicit none

    real*8, intent(in), dimension(nbasis_e) :: xe, pe 
    complex*16  ck(nbasis_e)
    ck = (xe + iota*pe)/sqrt(2.0)
  end function ck_xp

  function state_vec(k) result(vec)
    implicit none
    integer, intent(in) :: k 
    integer, dimension(nbasis_e):: vec
    vec = 0
    vec(k) = 1
  end function state_vec

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
    !*add at compile time : -L/usr/local/lib -llapack -lblas  !
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
      integer i, ngrids
      real*8, dimension(nclass) :: iarr, x
      real*8 xlim
      real*8, dimension(nbasis_e,nbasis_e) :: y1
      real*8, dimension(nbasis_e) :: y2
      real*8, dimension(nbasis_e,nbasis_e, nclass) :: y3
      !! without vibrational quantization : electronic hamiltonian
      
      !  1-D PES : along x(1) keeping x(2) at equilibrium
      iarr = 1.d0 !(/1.d0,0.d0/)!1.d0
      xlim = Rmax !5 !! au
      ngrids = 1000
      do i = 1, ngrids
        x= (-xlim + 2*xlim*i/ngrids)*iarr
        ! x(2) = g2*x(1)/(mass(2)*omg1**2) 

        !! diabatic pes
        if (ipes==1) then 
          y1 = ham_eDiab(x)
          if (wrt_au==0)  then 
              x   = x/ang2au 
              y1  = y1/eng2au
          endif
          write(200,*) x, y1(1,1), y1(1,2), y1(2,1), y1(2,2)
        endif
        
        !! adiabatic pes
        if (ipes==2) then
          call ham_eAdiab(x, nbasis_e, y2)
          if (wrt_au==0)  then 
              x   = x/ang2au 
              y2  = y2/eng2au
          endif
          write(201,*) x, y2(1), y2(2)
        endif

        !! derivative of pes
        if (ipes==3) then 
          y3 = delHeDiab_delR(x)
          if (wrt_au==0)  then 
              x   = x/ang2au 
              y3  = y3/eng2au
          endif
          write(203,*) x, y3(1,1,1), y3(1,2,1), y3(2,1,1), y3(2,2,1)
          !write(204,*) x, y3(1,1,2), y3(1,2,2), y3(2,1,2), y3(2,2,2)
        endif
      enddo
    end subroutine pes
    !-------------------------------------------------------------------------!
end module model_tully_1
