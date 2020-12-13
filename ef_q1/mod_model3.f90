module model_spin_boson_q1
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
  real*8, parameter :: en     = -900.d0*eng2au
  real*8, parameter :: en1    = 0.d0*eng2au     
  real*8, parameter :: en2    = -900.d0*eng2au   
  real*8, parameter :: vcoup  = 150.d0*eng2au
  real*8, parameter :: omg1   = 1200.d0*omg2au
  real*8, parameter :: omg2   = 400.d0*omg2au
  real*8, parameter :: eta    = 400.d0*omg2au
  real*8, parameter :: gamma_B= eta
  real*8, parameter :: lmda1  = 12000.d0*eng2au    
  real*8, parameter :: tempK  = 400.d0*tk2au
  real*8, parameter :: beta   = 1.d0/(kb*tempK)
  

  !! model dimensions
  integer, parameter :: nquant  = 2
  integer, parameter :: nclass  = 2     ! q1, q2
  integer, parameter :: nclass_q= 1     ! q1 mode quantised
  integer, parameter :: nclass_c= 1     ! q2 mode classical
  integer, parameter :: nbasis_e= 2
  integer, parameter :: nbasis_v= 25
  integer, parameter :: nbasis_ev = nbasis_e*nbasis_v
  integer, parameter :: ngrid   = 501   ! size of DVR grid
   

  ! real*8, dimension(nbasis_e, nbasis_e) :: H_el, H_el_diab
  !! real*8, dimension(nbasis_ev, nbasis_ev) :: H_vc_diab   !, H_adiab
  !! real*8, dimension(nbasis_ev, nbasis_ev, nclass) :: dV_dq      !, dij
  real*8, dimension(ngrid)  :: x_grid 
  real*8, dimension(ngrid, ngrid)  :: H_dvr_R, H_dvr_P, evec_R, evec_P
  real*8, dimension(ngrid)  :: eval_R, eval_P
  real*8, dimension(nbasis_ev, nbasis_ev) :: fc_mat, q_op, qsq_op

  real*8, parameter, dimension(nclass) :: mass = mp
  real*8, parameter, dimension(2) ::  xgrid_lim = [-5.0, 5.0]*ang2au     !! position limit of DVR basis

  real*8, dimension(nclass), parameter :: omega = (/omg1, omg2/) 
  real*8, parameter :: g1   = dsqrt(0.5*lmda1*mass(1)*omg1**2)      
  real*8, parameter :: g2   = 6.d-3       ! a.u.
  real*8, dimension(nclass), parameter :: g     = (/g1, g2/) 
  real*8, parameter :: delxg= (xgrid_lim(2) - xgrid_lim(1))/(ngrid - 1)     !! step size in position of DVR basis

  contains

!-------------------------------------------------------------------------!

  subroutine setup_ham_dvr()
    implicit none
    real*8  KE_dvr(ngrid,ngrid)
    real*8  H_dvr_R(ngrid, ngrid), H_dvr_P(ngrid, ngrid) 
    real*8 x  !(nclass_c)
    integer i,j
    print *, "nbasis_ev: ", nbasis_ev
    KE_dvr = compute_KE_matrix_dvr(mass(1))   !! subroutine call
    x_grid = (/(xgrid_lim(1) +delxg*i, i=1,ngrid)/)

    do i=1,ngrid
      x = x_grid(i)

      do j=1,ngrid
        if (i==j) then
          H_dvr_R(i,i) = KE_dvr(i,j) + 0.5*mass(1)*omg1**2*x**2 + g1*x
          H_dvr_P(i,i) = KE_dvr(i,j) + 0.5*mass(1)*omg1**2*x**2 - g1*x
        else 
          H_dvr_R(i,j) =  KE_dvr(i,j)
          H_dvr_P(i,j) =  KE_dvr(i,j)
        endif
      enddo
    enddo

    call diag_dvr(H_dvr_R, nbasis_v, eval_R, evec_R)
    call diag_dvr(H_dvr_P, nbasis_v, eval_P, evec_P)

    do i=1,ngrid   ! nbasis_v
      write(300, *) eval_R(i), eval_P(i)
    enddo

    call compute_overlap_mat(evec_R, evec_P, fc_mat, q_op, qsq_op)

    call write_matrix(fc_mat, nbasis_ev, 301)
    call write_matrix(q_op, nbasis_ev, 302) 
    call write_matrix(qsq_op, nbasis_ev, 303)

  end subroutine setup_ham_dvr

  !--------------------------!

  function compute_KE_matrix_dvr(m) result(KE_matrix)     !! m
    !! computes KE matrix in DVR basis
    !! Appendix A of JCP 96, 1982 (1991)
    implicit none
    real*8 KE_matrix(ngrid,ngrid)
 
    real*8,intent(in) :: m         !! mass
    integer i,j 
   
    KE_matrix=hbar**2/(2*m*delxg**2)
    do i=1,ngrid
      do j=1,ngrid
        KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
        if(i==j) then
          KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
          else
          KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
        endif
      enddo
    enddo
  end function compute_KE_matrix_dvr
  
  subroutine  diag_dvr(H_dvr, nbasis_v, evl, evc) 
    implicit none
    integer,  intent(in)  :: nbasis_v
    real*8,   intent(in)  :: H_dvr(ngrid,ngrid)
    real*8,   intent(out) :: evl(ngrid)
    real*8,   intent(out) :: evc(ngrid, ngrid)

    evc = H_dvr
    call diasym(evc, evl, ngrid, 'V') 
  end subroutine diag_dvr

  subroutine compute_overlap_mat(evec_R, evec_P, fc_mat, q_op, qsq_op)
    implicit none
    real*8, dimension(ngrid, ngrid),  intent(in)  :: evec_R(:,:), evec_P(:,:)
    real*8, dimension(nbasis_ev, nbasis_ev),intent(out) :: fc_mat(:,:)
    real*8, dimension(nbasis_ev, nbasis_ev),intent(out) :: q_op(:,:), qsq_op(:,:)
    integer i,j,k,l

    fc_mat  = 0.d0
    q_op    = 0.d0
    qsq_op  = 0.d0

    do i= 1, nbasis_v
      do j= 1, nbasis_v
        k = i + nbasis_v
        l = j + nbasis_v
        ! print *, shape(evec_R(:,i)*x_grid*evec_R(:,j))

        fc_mat(i,j) = sum(evec_R(:,i)*evec_R(:,j))    ! 1
        fc_mat(k,l) = sum(evec_P(:,i)*evec_P(:,j))    ! 1
        fc_mat(i,l) = sum(evec_P(:,i)*evec_R(:,j))    
        fc_mat(k,j) = sum(evec_R(:,i)*evec_P(:,j)) 

        q_op(i,j)=sum(evec_R(:,i)*x_grid*evec_R(:,j))
        q_op(k,l)=sum(evec_P(:,i)*x_grid*evec_P(:,j))

        qsq_op(i,j)=sum(evec_R(:,i)*x_grid**2*evec_R(:,j))
        qsq_op(k,l)=sum(evec_P(:,i)*x_grid**2*evec_P(:,j))
      enddo
    enddo

    ! call write_matrix(fc_mat, nbasis_ev, 301)
    ! call write_matrix(q_op, nbasis_ev, 302) 
    ! call write_matrix(qsq_op, nbasis_ev, 303)
  end subroutine compute_overlap_mat
!-------------------------------------------------------------------------!

  function ham_vibronic_diab(R) result(H_vc_diab)
      implicit none
      real*8, intent(in)  :: R(nclass_c)
      real*8  H_vc_diab(nbasis_ev,nbasis_ev)
      integer i,j,k,l
      real*8 rc

      !!!call setup_ham_dvr()    !!!
      ! do i=1,ngrid   ! nbasis_v
      !   write(333, *) eval_R(i), eval_P(i)
      ! enddo

      rc = R(1)      ! reaction coordinate
      !print *, rc

      do i= 1, nbasis_v
        do j= 1, nbasis_v
          k = i + nbasis_v
          l = j + nbasis_v
          !! (i,j), (k,l)...(m,n) more than 2 level electronic systems

          !! PP, RR
          H_vc_diab(i,j) = -g2*rc*q_op(i,j) + g2**2/(2*mass(2)*(omg2**2))*qsq_op(i,j) 
          H_vc_diab(k,l) = -g2*rc*q_op(k,l) + g2**2/(2*mass(2)*(omg2**2))*qsq_op(k,l)
          if (i==j) then
            H_vc_diab(i,i) = en1 + H_vc_diab(i,i) + eval_R(i) + 0.5*mass(2)*(omg2*rc)**2 
            H_vc_diab(k,k) = en2 + H_vc_diab(k,k) + eval_P(i) + 0.5*mass(2)*(omg2*rc)**2
          endif

          !! PR, RP 
          H_vc_diab(i,l) = vcoup*fc_mat(i,l)
          H_vc_diab(k,j) = vcoup*fc_mat(k,j)
        enddo
      enddo
  end function ham_vibronic_diab

  function delV_delq(R)   result(dVdr)
    implicit none
    real*8, intent(in) :: R(nclass_c)
    real*8  dVdr(nbasis_ev,nbasis_ev,nclass_c)
    integer i,j,k,l
    real*8 rc

    rc = R(1)      ! reaction coordinate
    !print *, rc       

    dVdr   = 0.d0
    do i= 1, nbasis_v
      do j= 1, nbasis_v
        k = i + nbasis_v
        l = j + nbasis_v

        !! PP, RR
        dVdr(i,j, 1) = -g2*q_op(i,j) 
        dVdr(k,l, 1) = -g2*q_op(k,l) 
        if (i==j) then
          dVdr(i,i, 1) = dVdr(i,i, 1) + mass(2)*(omg2**2)*rc 
          dVdr(k,k, 1) = dVdr(k,k, 1) + mass(2)*(omg2**2)*rc
        endif
      enddo
    enddo 
  end function delV_delq

  subroutine  ham_vibronic_adiab(R, nbasis, eval)
    implicit none
    real*8, intent(in)  :: R(nclass_c)
    integer, intent(in) :: nbasis
    real*8, intent(out) :: eval(nbasis)
    real*8 evec(nbasis,nbasis)

    evec = ham_vibronic_diab(R)
    !call write_matrix(evec, nbasis, 205) ! vibrionic hamiltonina
    call diasym(evec, eval, nbasis, 'N') ! only eigenvalues calculated     
  end subroutine ham_vibronic_adiab

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
  !           iev = 'N' for just eigenvalue                 !
  !                 'V' for both eigenvalue and eigenvectors!
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
    integer i, npts
    real*8, dimension(nclass_c) :: x  ! iarr
    real*8 xlim
    real*8, dimension(nbasis_ev,nbasis_ev) :: H
    real*8, dimension(nbasis_ev) :: E
    !! with vibrational quantization : vibronic hamiltonian
    
    ! 1-D PES : nclass_c = 1
    ! iarr = (/1.d0,0.d0/)!1.d0

    xlim = 4.d0*ang2au    !! au
    
    npts =1001
    do i = 1, npts
      x= (-xlim + 2*xlim*i/npts)!*iarr

      ! diabatic pes
      if (ipes==1) then 
        H = ham_vibronic_diab(x)
        write(200,*) x, H(1,1), H(1,2), H(2,1), H(2,2), H(3,1)
      endif

      ! adiabatic pes
      if (ipes==2) then 
        call  ham_vibronic_adiab(x, nbasis_ev, E)
        
        if (wrt_au==0)  then 
          x   = x/ang2au 
          E   = E/eng2au
        endif

        write(202,*) x, E(1), E(2), E(3), E(4), E(5)    ! eigenvalues : adiabatic PES
      endif
    enddo
  end subroutine pes

  subroutine write_matrix(M, n, fname)
    implicit none
    real*8, dimension(n,n), intent(in) :: M(:,:)
    integer, intent(in) :: n , fname
    integer i,j 

    do i = 1,n
      write(fname,*) (M(i,j), j = 1,n)
    enddo
  end subroutine write_matrix
  !-------------------------------------------------------------------------!
end module model_spin_boson_q1
