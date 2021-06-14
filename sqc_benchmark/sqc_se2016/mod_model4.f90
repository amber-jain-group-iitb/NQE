module model_spinboson_b
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
	
	!! model parameters: spin-boson
  ! specify spectral density of continuous ohmic distribution  
  ! using Kondo param (alpha) & characteristic frequency (omgc)
  real*8, parameter :: vcoup  = 10.d0*eng2au
  real*8, parameter :: omgc   = 2.5d0*vcoup
  real*8, parameter :: alpha  = 0.09d0 
  real*8, parameter :: eps    = 0.d0

  ! real*8, parameter :: tempK  = !77.d0*tk2au
  real*8, parameter :: beta   = 0.1d0/vcoup !1.d0/(kb*tempK)
  real*8, parameter :: ngamma = 0.366
  real*8, parameter :: tempK  = 1/(kb*beta)/tk2au

	!! model dimensions
	integer, parameter :: nbath = 100
	integer, parameter :: nclass = nbath !1    ! q1, q2
	integer, parameter :: nbasis_e = 2
	!integer, parameter :: nbasis = nbasis_e 

	real*8, dimension(nclass)  :: mass
  real*8, dimension(nclass)  :: omega, Johm, cpk
  real*8 omg_max

  integer, parameter, dimension(nbasis_e) :: N_s1 = [1,0]
  integer, parameter, dimension(nbasis_e) :: N_s2 = [0,1]

	! real*8, dimension(nbasis_e, nbasis_e) :: H_el, H_el_diab, H, H_new
	! real*8, dimension(nbasis_ev, nbasis_ev) :: H_vc, H_vc_diab, H_vc_new   !, H_adiab
	! real*8, dimension(nbasis_e, nbasis_e, nclass) :: dV_dq      !, dij

  ! real*8, dimension(nclass) :: Rn, Pn 
	contains

  subroutine setup_model()
    implicit none

    real*8 domg
    integer i

    mass  = mp
    omg_max = 3*omgc
    domg  = omg_max/nbath

    omega = (/(i*domg, i=1,nbath)/)
    Johm  = (pi/2.0*alpha)*omega*exp(-omega/omgc) ! Ohmic spectral density
    cpk   = sqrt(2.0/pi*domg*mass*omega*Johm)     ! coupling parameter
  end subroutine setup_model
  
  function ham_eDiab(Rn) result(H_e)
    implicit none

    real*8, intent(in) :: Rn(nclass)
    real*8  V(nbasis_e)
    real*8, dimension(nbasis_e, nbasis_e) ::  H_e
    
    V(1) = sum(0.5*mass*omega**2*Rn**2)
    V(2) = sum(cpk*Rn)

    ! Electronic Hamiltonian
    H_e   = 0.d0
    H_e(1,1) = V(1)+V(2)
    H_e(2,2) = V(1)-V(2)
    H_e(1,2) = vcoup
    H_e(2,1) = H_e(1,2)
  end function ham_eDiab

  function ham_ne(Pn, Rn, pe, xe) result(H_ne)
    implicit none

    real*8, intent(in), dimension(nclass)   :: Pn, Rn
    real*8, intent(in), dimension(nbasis_e) :: pe, xe
    real*8  H_e(nbasis_e, nbasis_e)
    real*8  H_el, H_ne  ! Hamiltonian function

    H_e = ham_eDiab(Rn)
    ! Electronic Hamiltonian
    H_el  = 0.5*(H_e(1,1)+H_e(2,2))
    H_el  = H_el+ 0.25*(pe(1)**2+xe(1)**2 - pe(2)**2-xe(2)**2)*(H_e(1,1)-H_e(2,2))
    H_el  = H_el+ (pe(1)*pe(2) + xe(1)*xe(2))*H_e(1,2) 

    ! Nuclear-Electronic Hamiltonian
    H_ne  = sum(Pn**2/(2*mass)) + H_el

  end function ham_ne

  function delHeDiab_delR(Rn)   result(dV_dR)
    implicit none
    
    real*8, intent(in) :: Rn(nclass)
    real*8  dV_dR(nbasis_e,nbasis_e,nclass)

    dV_dR   = 0.d0
    dV_dR (1,1,:) = mass*omega**2*Rn + cpk
    dV_dR (2,2,:) = mass*omega**2*Rn - cpk
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
    !call nq_xp
    ! qk = arctan(-pe/xe)
    ! nk = sum(0.5*(xe)**2 + 0.5*(pe)**2) + ngamma
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
      xlim = 0.25 !Rmax !5 !! au
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
end module model_spinboson_b
