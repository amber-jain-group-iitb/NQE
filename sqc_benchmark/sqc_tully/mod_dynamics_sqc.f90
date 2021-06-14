module dynamics
  use model_tully_1
  implicit none
  
  integer     :: iflow, ncore, ntraj, seed
  real*8      :: dtc, ttime
  integer     :: ifriction, istate
  integer     :: iprint, istr_sci, istr_snk, istr_eng, wrt_qs, nwrt_traj
  integer     :: wrt_tstate
  integer     :: nsteps, nreadouts
  integer     :: itraj_exit

  real*8,     allocatable :: H(:,:), H_new(:,:)
  real*8,     allocatable :: dV_dq(:,:,:)

  complex*16, allocatable :: ci(:) 
  real*8,     allocatable :: xe(:), pe(:)
  integer,    allocatable :: Nnk(:)
  real*8,     allocatable :: q(:), qdot(:), acc(:)
  real*8,     allocatable :: force(:)
  real*8 Pn !nclass

  ! real*8,     allocatable :: tsteps(:)
  real*8,     allocatable :: sum_ci_sq(:,:), tot_energy(:,:), sum_nk(:,:)
  real*8,     allocatable :: pop(:), pop_sum_xp(:) ! pop_sum(:), pop_avg(:)
  integer,    allocatable :: initial_states(:), final_states(:)
  integer,    allocatable :: count_sts(:)
      
  contains
  subroutine setup_dynamics()
    implicit none

    ! read runtime parameters
    open(10, file="input0.txt")
    read(10,*) iflow
    read(10,*) ncore
    read(10,*) ntraj
    
    read(10,*) dtc 
    read(10,*) ttime
    
    read(10,*) ifriction ! 0
    read(10,*) istate

    read(10,*) iprint
    read(10,*) nwrt_traj
    read(10,*) istr_sci
    read(10,*) istr_snk
    read(10,*) istr_eng
    read(10,*) wrt_qs
    read(10,*) wrt_tstate
    if (iflow == 1) read(10,*) seed
    if (iflow == 5) read(10,*) Pn
    close(10)

    nsteps      = ceiling(ttime/dtc)
    nreadouts   = ceiling(real(nsteps)/nwrt_traj)

    ! print *, iflow, ncore, ntraj 
    ! print *, dtc, ttime, nsteps, ifriction, istate
    ! print *, iprint, istr_sci, istr_eng, wrt_qs
    ! print *, nsteps, nreadouts, seed

    ! allocate arrays
    allocate( H(nbasis_e,nbasis_e), H_new(nbasis_e,nbasis_e) )
    allocate( dV_dq(nbasis_e,nbasis_e,nclass) )
    allocate( ci(nbasis_e) )
    allocate( xe(nbasis_e), pe(nbasis_e) )
    allocate( Nnk(nbasis_e) )
    allocate( q(nclass), qdot(nclass), acc(nclass) )
    allocate( force(nclass))

    ! allocate( tsteps(nsteps) )
    allocate( sum_ci_sq(1,nreadouts), sum_nk(1,nreadouts) ) !(ntraj,nsteps)
    allocate( tot_energy(1,nreadouts) )
    allocate( pop(nreadouts), pop_sum_xp(nbasis_e) ) !pop_sum(nreadouts)
    allocate( initial_states(ntraj), final_states(ntraj) ) 
    allocate( count_sts(-nbasis_e:nbasis_e) ) !count_sts(nbasis_e) 
  end subroutine setup_dynamics

  !----------------------------------------------------!

  subroutine init_condition(P)
    implicit none
    ! integer :: i
    ! real*8, dimension(nclass) :: sig_q, sig_qdot
    real*8 :: rnd(2*nbasis_e) 
    real*8, intent(in) :: P !nclass
    real*8 nk(nbasis_e), qk(nbasis_e)    
    
    ! do i = 1, nclass
    !   sig_q(i)    = 1.d0/dsqrt(beta*mass(i)*omega(i)**2) 
    !   sig_qdot(i) = 1.d0/dsqrt(beta*mass(i))

    !   call gauss_rand_n(rnd1)
    !   call gauss_rand_n(rnd2)
    !   ! print *, rnd1, rnd2
      
    !   q(i)       = -g(i)/(mass(i)*omega(i)**2) + rnd1*sig_q(i)
    !   qdot(i)    = rnd2*sig_qdot(i)
    ! enddo
    itraj_exit = 0
    q  = -9.99
    qdot = P/mass
    ! write(*,*) P, qdot

    H       = ham_eDiab(q)
    dV_dq   = delHeDiab_delR(q)

    Nnk = N_s1  !state_vec(istate)   ! state 1
    call random_number(rnd)
    nk = Nnk + ngamma*[2*rnd(1)-1, 2*rnd(2)-1]
    qk = [2*pi*rnd(3), 2*pi*rnd(4)]
    call nq_xp(nk,qk, xe, pe)

    call compute_force(xe, pe, dV_dq, force)
    write(*,*) "rnd",   rnd
    write(*,*) "init",  q, qdot, force
    write(*,*) "eng",   calc_tot_eng_Hn(qdot, q, pe, xe)
    ! write(*,*) nk, qk
    ! write(*,*) "ck^2: ", calc_ci_ss(ck_nq(nk,qk))
    ! write(*,*) xe, pe
    ! write(*,*) "ck^2: ", calc_ci_ss(ck_xp(xe,pe))
  end subroutine init_condition

  !----------------------------------------------------------!

  subroutine compute_force(xe, pe, dV_dq, force)
    implicit none
    integer i
    real*8,intent(in), dimension(nbasis_e)  :: xe, pe
    real*8,intent(in) :: dV_dq(nbasis_e,nbasis_e,nclass)
    real*8,intent(out):: force(nclass)
    ! complex*16 ci((nbasis_e))
    
    ! ci = ck_xp(xe,pe)
    do i = 1, nclass
      !force(i) =- real(sum(conjg(ci)*matmul(delV_diab(:,:,i), ci)))
      ! force(i) = - expec_val(ci, dV_dq(:,:,i), ci)
      force(i) = 0.25*(pe(1)**2+xe(1)**2 - pe(2)**2-xe(2)**2)*(dV_dq(1,1,i)-dV_dq(2,2,i))
      force(i) = force(i)+ (pe(1)*pe(2) + xe(1)*xe(2))*dV_dq(1,2,i)
      force(i) = force(i)+ 0.5*(dV_dq (1,1,1) + dV_dq (2,2,1))
      force(i) = -force(i)
    enddo
  end subroutine compute_force

  subroutine evolve_classical(q, qdot, force)
    ! Velocity Verlet
    implicit none
    real*8, intent(inout), dimension(nclass) :: q, qdot, force
    ! real*8, dimension(nclass) :: delr,delv
    ! real*8 c0, c1, c2, gdt

    if (ifriction == 0) then 
        acc = force/mass
        q   = q + qdot*dtc + 0.5*acc*(dtc**2)
        dV_dq   = delHeDiab_delR(q)
        call compute_force(xe, pe, dV_dq, force) !compute_force(ci,dV_dq,force)

        acc     = 0.5*(acc + force/mass)
        qdot    = qdot + acc*dtc   
    endif
    
  end subroutine evolve_classical

  ! function xp_eqnM(xe, pe, H_e) result(m1m2)
  !   implicit none
  !   real*8, intent(in), dimension(nbasis_e)  :: xe, pe
  !   real*8, intent(in)  :: H_e(nbasis_e, nbasis_e)
  !   real*8 m1(2*nbasis_e,2), m2(2)
  !   real*8  m1m2(2*nbasis_e)

  !   ! do i=1,nbasis_e
  !   !   do j=i,nbasis_e
  !   !   k = j+nbasis_e
  !   !     m1(j,1) = (1.0/nbasis_e)*pe(i) 
  !   !     m1(j,2) = pe(j)
  !   !     m1(k,1) = (1.0/nbasis_e)*xe(i)
  !   !     m1(k,2) = xe(j)

  !   !     m2(1) = H_e(i,i)-H_e(j,j)
  !   !     m2(2) = H_e(i,j)
  !   !   enddo
  !   ! enddo   

  !   m1(1,:) = [0.5*pe(1), pe(2)]
  !   m1(2,:) = [-0.5*pe(2), pe(1)]
  !   m1(3,:) = -[0.5*xe(1), xe(2)]
  !   m1(4,:) = -[-0.5*xe(2), xe(1)]

  !   m2(1) = H_e(1,1)-H_e(2,2)
  !   m2(2) = H_e(1,2)

  !   m1m2 = matmul(m1, m2) !
  ! end function

  ! subroutine evolve_quantumM(xe, pe, Hel_k1, Hel_k4)
  !   ! Runge-Kutta 4th order integrator
  !   ! J. Chem. Phys. 139, 234112 (2013)
  !   ! Evolving H_el(R,xe,pe) part of eq.10
  !   implicit none

  !   real*8, intent(inout), dimension(nbasis_e)  :: xe, pe
  !   real*8, intent(in), dimension(nbasis_e, nbasis_e) :: Hel_k1, Hel_k4  
  !   real*8, dimension(1:2*nbasis_e) :: xepe, k1, k2, k3, k4
  !   real*8  Hel_k23(nbasis_e, nbasis_e)
  !   integer, dimension(nbasis_e)   :: xe_idx, pe_idx
  !   integer i

  !   xe_idx = (/(i, i=1,nbasis_e)/)            ! 1, 2
  !   pe_idx = (/(i, i=nbasis_e+1,2*nbasis_e)/) ! 3, 4
  !   xepe = (/xe,pe/)

  !   Hel_k23 = (Hel_k1 + Hel_k4)/2.d0    ! H_el{R_t}, H_el{R_t+1} 

  !   k1 = xp_eqn(xe, pe, Hel_k1)
  !   k2 = xp_eqn(xe+dtc*k1(xe_idx)/2.d0, pe+dtc*k1(pe_idx)/2.d0, Hel_k23)
  !   k3 = xp_eqn(xe+dtc*k2(xe_idx)/2.d0, pe+dtc*k2(pe_idx)/2.d0, Hel_k23)
  !   k4 = xp_eqn(xe+dtc*k3(xe_idx), pe+dtc*k3(pe_idx), Hel_k4)

  !   ! k1p = xp_eqn(xe, pe, Hel_k1, 2)
  !   ! k2p = xp_eqn(xe, pe+dtc*k1p/2.d0, Hel_k23, 2)
  !   ! k3p = xp_eqn(xe, pe+dtc*k2p/2.d0, Hel_k23, 2)
  !   ! k4p = xp_eqn(xe, pe+dtc*k3p, Hel_k4, 2)

  !   ! xe = xe + dtc*(k1+2*k2+2*k3+k4)/6.d0
  !   ! pe = pe + dtc*(k1p+2*k2p+2*k3p+k4p)/6.d0

  !   xepe = xepe + dtc*(k1+2*k2+2*k3+k4)/6.d0
  !   xe = xepe(xe_idx)
  !   pe = xepe(pe_idx) 
  ! end subroutine evolve_quantum

  function xedot(H_e, pe) result(xe_dot)
    implicit none
    real*8, intent(in)  :: H_e(nbasis_e, nbasis_e)
    real*8, intent(in), dimension(nbasis_e)  :: pe
    real*8  xe_dot(nbasis_e)

    xe_dot(1) = 0.5*pe(1)*(H_e(1,1)-H_e(2,2)) + pe(2)*H_e(1,2)
    xe_dot(2) = -0.5*pe(2)*(H_e(1,1)-H_e(2,2)) + pe(1)*H_e(1,2)
  end function xedot

  function pedot(H_e, xe) result(pe_dot)
    implicit none
    real*8, intent(in)  :: H_e(nbasis_e, nbasis_e)
    real*8, intent(in), dimension(nbasis_e)  :: xe
    real*8  pe_dot(nbasis_e)

    pe_dot(1) = -(0.5*xe(1)*(H_e(1,1)-H_e(2,2)) + xe(2)*H_e(1,2))
    pe_dot(2) = -(-0.5*xe(2)*(H_e(1,1)-H_e(2,2)) + xe(1)*H_e(1,2))
  end function pedot

  subroutine evolve_quantum(xe, pe, Hel_k1, Hel_k4)
    ! Runge-Kutta 4th order integrator
    ! J. Chem. Phys. 139, 234112 (2013)
    ! Evolving H_el(R,xe,pe) part of eq.10
    implicit none

    real*8, intent(inout), dimension(nbasis_e)  :: xe, pe
    real*8, intent(in), dimension(nbasis_e, nbasis_e) :: Hel_k1, Hel_k4  
    real*8, dimension(nbasis_e) :: k1x, k2x, k3x, k4x
    real*8, dimension(nbasis_e) :: k1p, k2p, k3p, k4p
    real*8  Hel_k23(nbasis_e, nbasis_e)


    Hel_k23 = (Hel_k1 + Hel_k4)/2.d0    ! H_el{R_t}, H_el{R_t+1} 

    k1x = xedot(Hel_k1, pe)
    k1p = pedot(Hel_k1, xe)

    k2x = xedot(Hel_k23, pe+dtc*k1p/2)
    k2p = pedot(Hel_k23, xe+dtc*k1x/2)

    k3x = xedot(Hel_k23, pe+dtc*k2p/2)
    k3p = pedot(Hel_k23, xe+dtc*k2x/2)

    k4x = xedot(Hel_k4, pe+dtc*k3p)
    k4p = pedot(Hel_k4, xe+dtc*k3x)

    xe = xe + dtc*(k1x+2*k2x+2*k3x+k4x)/6.d0
    pe = pe + dtc*(k1p+2*k2p+2*k3p+k4p)/6.d0

  end subroutine evolve_quantum

  subroutine run_dynamics()
    implicit none

    integer ns, nt, i, j
    integer fs, ipos
    ! integer :: s(35)    ! seeder

    ! !call setup()   !! setup dynamics in main program
    ! s = 10
    ! if(iflow==1) s = s*seed
    ! print *, s
    ! call random_seed(put=s)

    ! initialize trajectory observables
    pop_sum_xp  = 0.d0
    sum_ci_sq   = 0.d0
    sum_nk      = 0.d0
    tot_energy  = 0.d0
    count_sts   = 0

    all_traj: do nt = 1, ntraj

      ! Pn = 10.0
      call init_condition(Pn)
      initial_states(nt) = calc_state(xe, pe) !initial state of each trajectory
      j = 1   ! variable to keep track of readout

      traj_evolve: do ns = 1, nsteps
        ! write(*,*) "obs1", ns, nwrt_traj  ,ntraj 
        !----------write observables----------!
        if (iprint==1)  print *, nt, q
        if ( (mod(ns, nwrt_traj)==1) )  then !(mod(ns, nwrt_traj)==1 )
          if (ntraj==1) then 
            if (istr_sci==1) sum_ci_sq(nt,j)  = calc_ci_ss(ck_xp(xe,pe))  ! ci
            if (istr_snk==1) sum_nk(nt,j)     = calc_sum_nk(xe,pe)
            if (istr_eng==1) tot_energy(nt,j) = calc_tot_eng_Hn(qdot,q,pe,xe) 
            if (wrt_qs==1)   write(105,*)  nt, ns, q
          endif 
            ! pop(j) = calc_state(xe, pe) !calc_pop(ci, istate)
            j = j+1
        endif
        !-------------------------------------!

        !----------evolve trajectory----------!
        call evolve_classical(q, qdot, force)
        H_new   = ham_eDiab(q)
        call evolve_quantum(xe, pe, H, H_new)
        H       = H_new
        !-------------------------------------!
        
        call traj_exit_check(q, ns)
        if (itraj_exit == 1) then
          ipos = int(sign(1.d0,q(1)))
          final_states(nt) = calc_state(xe, pe)*ipos
          fs = final_states(nt)
          count_sts(fs) = count_sts(fs) + 1
          exit traj_evolve
        endif
      enddo traj_evolve

      if (iprint==2) print *, nt, ns
      ! pop_sum = pop_sum + pop
      ! call calc_state_pop(xe, pe, pop_sum_xp)
    enddo all_traj
    call writefiles(ntraj)
  end subroutine run_dynamics

  !----------------------------------------------------------!

  function calc_pop_ci(ci, istate) result(popn)
    implicit none
    complex*16, intent(in) :: ci(nbasis_e)
    integer,    intent(in) :: istate
    real*8 popn(nbasis_e)
    popn = abs(ci(istate))**2
  end function calc_pop_ci

  function calc_state(xe, pe) result(k_state)
    implicit none
    real*8, intent(in), dimension(nbasis_e) :: xe, pe
    integer k_state
    logical isState1, isState2
    real*8 var(nbasis_e)
    ! write(*,*) "state (xe,pe): ", xe, pe
    ! ci_sq = (xe)**2 + (pe)**2
    var = 0.5*((xe)**2 + (pe)**2) - ngamma

    !isStatek(k) = all(abs(var-state_vec(k)) <= ngamma, 1)
    isState1 = all(abs(var-N_s1) <= ngamma, 1)
    isState2 = all(abs(var-N_s2) <= ngamma, 1)
    
    ! write(*,*) "iss", isState1, isState2 
    if (isState1)  then 
      k_state = 1
    elseif (isState2) then
      k_state = 2
    else 
      k_state = 0
    endif
  end function calc_state

  subroutine calc_state_pop(xe, pe, pop_sum_xp) 
    implicit none
    real*8, intent(in), dimension(nbasis_e)  :: xe, pe
    real*8,intent(out) :: pop_sum_xp(nbasis_e)
    integer state 

    state = calc_state(xe,pe)
    pop_sum_xp(state) = pop_sum_xp(state) + 1.0
    ! if (calc_state(xe,pe) == 1) then
    !   pop_sum_xp(1) = pop_sum_xp(1) + 1.0
    ! else
    !   pop_sum_xp(2) = pop_sum_xp(2) + 1.0
    ! endif
  end subroutine calc_state_pop

  function calc_ci_ss(ci) result(ci_ss)
    implicit none
    complex*16, intent(in)  :: ci(nbasis_e)
    real*8 ci_ss
    ci_ss   = sum(abs(ci)**2)
  end function calc_ci_ss

  function calc_sum_nk(xe,pe) result(sum_nk)
    implicit none
    real*8, intent(in), dimension(nbasis_e) :: xe, pe
    real*8 sum_nk
    sum_nk   = sum(0.5*(xe)**2 + 0.5*(pe)**2) - nbasis_e*ngamma
  end function calc_sum_nk

  function calc_tot_eng_ci(ci, H) result(tot_eng)
    implicit none
    complex*16, intent(in)  :: ci(nbasis_e)
    real*8, intent(in)  :: H(nbasis_e,nbasis_e) ! H_e = ham_eDiab(Rn)
    real*8 tot_eng      !! electronic energy(PE) + KE
    tot_eng = expec_val(ci, H, ci) + sum(0.5*mass*qdot**2)
  end function calc_tot_eng_ci

  function calc_tot_eng_Hn(qdot, q, pe, xe) result(tot_eng)
    implicit none
    real*8, intent(in), dimension(nclass)   :: qdot, q 
    real*8, intent(in), dimension(nbasis_e) :: pe, xe 
    real*8 tot_eng      !! electronic energy(PE) + KE

    tot_energy = ham_ne(mass*qdot, q, pe, xe)
  end function calc_tot_eng_Hn
  
  !-------------------------------------!  
  subroutine gauss_rand_n(rnd)    !! function
    implicit none
    real*8,intent(out)::rnd ! gaussian around 0??
    real*8 rnd1,rnd2
    
    call random_number(rnd1)
    call random_number(rnd2)
    rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)
  end subroutine gauss_rand_n

  subroutine traj_exit_check(q, ns)
    implicit none

    real*8, intent(in) :: q(nclass)
    integer, intent(in):: ns
    if (q(1)>Rmax .or. q(1)<Rmin) then
      itraj_exit = 1
      write(*,*) 'traj exit R: ', q(1)
    endif

    if (ns > nsteps) then
      itraj_exit = 1
      write(*,*) 'traj exit ns: ', ns
    endif
  end subroutine traj_exit_check    
  !----------------------------------------------------------!

  subroutine writefiles(flag)
    implicit none
    !-----------------------------------------------!
    ! flag  = ntraj = 1 write a traj + pop data     !
    !               > 1 write only population data  !
    !-----------------------------------------------!
    integer, intent(in) :: flag
    real*8 tsteps     
    integer i, j, fs
    j = 1
    ! tsteps = (/(i*dtc, i = 1, nsteps)/)
    do i = 0, nsteps, nwrt_traj
      tsteps = i*dtc
      
      if (flag==1) then
        if (istr_sci==1) write(102,*) tsteps, sum_ci_sq(flag,j)
        if (istr_snk==1) write(103,*) tsteps, sum_nk(flag,j)
        if (istr_eng==1) write(104,*) tsteps, tot_energy(flag,j)
      endif

      ! if (iflow==0)   write(100,*) tsteps, pop_sum_xp(j)/ntraj   ! average population
      ! if (iflow==1)   write(100,*) tsteps, pop_sum_xp(j)         ! population_sum
      j = j+1
    enddo

    do i = 1, ntraj
      if (wrt_tstate==1) write(111,*) i, initial_states(i), final_states(i)
    enddo

    open(112,file="trajStates.dat", status="new")
    do i =-nbasis_e, nbasis_e
      write(112,*) Pn, i, count_sts(i)
    enddo
    close(112)
  end subroutine writefiles
  !----------------------------------------------------------!
end module dynamics
