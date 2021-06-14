module dynamics
  use model_spin_boson_q0
  implicit none
  
  integer     :: iflow, ncore, ntraj, seed
  real*8      :: dtc, ttime
  integer     :: ifriction, istate
  integer     :: iprint, istr_sci, istr_snk, istr_eng, wrt_qs, nwrt_traj
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
  ! real*8,     allocatable :: pop_sum_xp(:) ! pop_sum(:), pop_avg(:)
  integer,    allocatable :: initial_states(:)!, final_states(:) 
  integer,    allocatable :: traj_states(:,:)
  real*8,     allocatable :: pop(:,:)
  integer,    allocatable :: summ_all(:), summ_nz(:)
  ! integer,    allocatable :: count_sts(:)
      
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
    if (iflow == 1) read(10,*) seed
    close(10)

    nsteps      = nint(ttime/dtc)
    nreadouts   = nint(real(nsteps)/nwrt_traj)

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
    ! allocate( pop_sum_xp(nbasis_e) ) !pop_sum(nreadouts)
    !! [nbasis_e+1 x nreadouts+1]
    allocate( traj_states(0:nbasis_e,0:nreadouts), pop(0:nbasis_e,0:nreadouts) )
    allocate( initial_states(ntraj) ) 
    allocate( summ_all(nreadouts+1), summ_nz(nreadouts+1)  )
    ! allocate( count_sts(0:nbasis_e) ) !count_sts(nbasis_e) 
  end subroutine setup_dynamics

  !----------------------------------------------------!

  subroutine init_condition()
    implicit none
    ! integer :: i
    real*8, dimension(nclass) :: sig_q, sig_qdot
    real*8 :: rnd_xp(2*nbasis_e), rnd1, rnd2 
    real*8 nk(nbasis_e), qk(nbasis_e)
    integer i    
    
    do i = 1, nclass
      sig_q(i)    = 1.d0/dsqrt(beta*mass(i)*omega(i)**2) !(i)
      sig_qdot(i) = 1.d0/dsqrt(beta*mass(i))

      call gauss_rand_n(rnd1)
      call gauss_rand_n(rnd2)
      ! print *, rnd1, rnd2
      
      q(i)       = rnd1*sig_q(i) -g(i)/(mass(i)*omega(i)**2)
      qdot(i)    = rnd2*sig_qdot(i)
    enddo
    itraj_exit = 0

    H       = ham_eDiab(q)
    dV_dq   = delHeDiab_delR(q)

    Nnk = N_s1  !state_vec(istate)   ! state 1
    call random_number(rnd_xp)
    nk = Nnk + ngamma*[2*rnd_xp(1)-1, 2*rnd_xp(2)-1]
    qk = [2*pi*rnd_xp(3), 2*pi*rnd_xp(4)]
    call nq_xp(nk,qk, xe, pe)

    call compute_force(xe, pe, dV_dq, force)
    ! write(*,*) "rnd",   rnd_xp
    write(*,*) "init",  q(1), qdot(1), force(1)
    write(*,*) "eng",   calc_tot_eng_Hn(qdot, q, pe, xe)
    ! write(*,*) nk, qk
    ! write(*,*) "ck^2: ", calc_ci_ss(ck_nq(nk,qk))
    ! write(*,*) xe, pe
    ! write(*,*) "ck^2: ", calc_ci_ss(ck_xp(xe,pe))

    if (ntraj==1) then
      write(102,*) 0.d0, calc_ci_ss(ck_xp(xe,pe))
      write(103,*) 0.d0, calc_sum_nk(xe,pe)
      write(104,*) 0.d0, calc_tot_eng_Hn(qdot, q, pe, xe)
    endif
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
      force(i) = force(i)+ 0.5*(dV_dq (1,1,i) + dV_dq (2,2,i))
      force(i) = -force(i)
    enddo
  end subroutine compute_force

  subroutine evolve_classical(q, qdot, force)
    ! Velocity Verlet
    implicit none
    real*8, intent(inout), dimension(nclass) :: q, qdot, force
    real*8, dimension(nclass) :: delr,delv
    real*8 c0, c1, c2, gdt

    if (ifriction == 0) then 
      acc = force/mass
      q   = q + qdot*dtc + 0.5*acc*(dtc**2)

      dV_dq   = delHeDiab_delR(q)
      call compute_force(xe, pe, dV_dq, force)

      acc     = 0.5*(acc + force/mass)
      qdot    = qdot + acc*dtc   
    endif
    
    if (ifriction == 1) then
      gdt = gamma_B*dtc

      c0 = exp(-gdt)
      c1 = 1.d0/gdt*(1-c0)
      c2 = 1.d0/gdt*(1-c1)

      acc = force/mass
      call stochastic_force(delr,delv)

      q(1)    = q(1) + qdot(1)*dtc + 0.5*acc(1)*(dtc**2)
      q(2)    = q(2) + c1*dtc*qdot(2) + c2*(dtc**2)*acc(2) + delr(2)

      dV_dq   = delHeDiab_delR(q)
      call compute_force(xe, pe, dV_dq, force)

      acc(1)  = 0.5*(acc(1) + force(1)/mass(1))
      acc(2)  = (c1-c2)*acc(2) + c2*(force(2)/mass(2))
      qdot(1) = qdot(1) + acc(1)*dtc
      qdot(2) = c0*qdot(2) + acc(2)*dtc + delv(2)    
    endif
  end subroutine evolve_classical

  subroutine stochastic_force(delr,delv)
    ! stochastic forces for langevin equation
    implicit none
    integer i
    real*8 :: rnd1,rnd2,sig_r,sig_v,sig_rv, gdt
    real*8, dimension(nclass), intent(out) :: delr(:),delv(:)
  
    gdt=gamma_B*dtc
  
    do i=1,nclass
      sig_r=dtc*dsqrt(kb*tempK/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
      sig_v=dsqrt(kb*tempK/mass(i)*(1-dexp(-2*gdt)))
      sig_rv=(dtc*kb*tempK/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient
  
      call gauss_rand_n(rnd1)
      call gauss_rand_n(rnd2)
      delr(i)=sig_r*rnd1
      delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
    enddo      
  end subroutine stochastic_force

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
    integer state_val
    integer :: s(35)    ! seeder

    !call setup()   !! setup dynamics in main program
    s = 10
    if(iflow==1) s = s*seed
    print *, s
    call random_seed(put=s)

    ! initialize trajectory observables
    sum_ci_sq   = 0.d0
    sum_nk      = 0.d0
    tot_energy  = 0.d0
    traj_states = 0
    pop         = 0.d0

    ! pop_sum_xp  = 0.d0
    ! count_sts   = 0
    
    all_traj: do nt = 1, ntraj

      call init_condition()
      !initial state of each trajectory
      state_val               = calc_state(xe, pe)
      initial_states(nt)      = state_val 
      traj_states(state_val,0)= traj_states(state_val,0) + 1
      j = 1   ! variable to keep track of readout

      traj_evolve: do ns = 1, nsteps
        ! write(*,*) "obs1", ns, nwrt_traj  ,ntraj 
        !----------write observables----------!
        if (iprint==1)  print *, nt, ns, q(100)
        if ( (mod(ns, nwrt_traj)==1) )  then !(mod(ns, nwrt_traj)==1 )
          if (ntraj==1) then 
            if (istr_sci==1) sum_ci_sq(nt,j)  = calc_ci_ss(ck_xp(xe,pe))  ! ci
            if (istr_snk==1) sum_nk(nt,j)     = calc_sum_nk(xe,pe)
            if (istr_eng==1) tot_energy(nt,j) = calc_tot_eng_Hn(qdot,q,pe,xe) 
            if (wrt_qs==1)   write(105,*)  nt, ns, q(100)
          endif
          state_val = calc_state(xe, pe)
          traj_states(state_val,j) = traj_states(state_val,j) + 1
          j = j+1
        endif
        !-------------------------------------!

        !----------evolve trajectory----------!
        call evolve_classical(q, qdot, force)
        H_new   = ham_eDiab(q)
        call evolve_quantum(xe, pe, H, H_new)
        H       = H_new
        !-------------------------------------!
        
        ! ! call traj_exit_check(q, ns)
        ! if (ns == nsteps) then !(itraj_exit == 1)
        !   final_states(nt) = calc_state(xe, pe)
        !   fs = final_states(nt)
        !   count_sts(fs) = count_sts(fs) + 1
        !   exit traj_evolve
        ! endif
      enddo traj_evolve

      if (iprint==2) print *, nt, ns
    enddo all_traj
    write(*,*) "-----All trajectories evolved-----"
    call calc_state_pop(traj_states)
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
    ! write(*,*) k_state
  end function calc_state

  ! subroutine store_states(state_val, ns) 
  !   implicit none
  !   integer, intent(in) :: state_val, ns

  !   traj_states(state_val,ns) = traj_states(state_val,ns) + 1
  ! end subroutine store_states

  subroutine calc_state_pop(traj_states)
    implicit none

    integer, intent(in) :: traj_states(0:nbasis_e,0:nreadouts)
    integer i ,j

    summ_all = sum(traj_states,DIM=1)
    summ_nz  = summ_all - traj_states(0,:)
    do i = 1, nreadouts
      ! write(121,*) i, summ_all(i) 
      do j = 1,nbasis_e
        if (j==0) then
          pop(j,i) =  traj_states(0,i)/real(summ_all(i))
        else
          pop(j,i) =  traj_states(j,i)/real(summ_nz(i))
        endif
      enddo
    enddo
    
    ! write(*,*) pop(0,:)
             
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

    tot_eng = ham_ne(mass*qdot, q, pe, xe)
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

  ! subroutine traj_exit_check(q, ns)
  !   implicit none

  !   real*8, intent(in) :: q(nclass)
  !   integer, intent(in):: ns
  !   if (q(1)>Rmax .or. q(1)<Rmin) then
  !     itraj_exit = 1
  !     write(*,*) 'traj exit R: ', q(1)
  !   endif

  !   if (ns > nsteps) then
  !     itraj_exit = 1
  !     write(*,*) 'traj exit ns: ', ns
  !   endif
  ! end subroutine traj_exit_check    
  !----------------------------------------------------------!

  subroutine writefiles(flag)
    implicit none
    !-----------------------------------------------!
    ! flag  = ntraj = 1 write a traj + pop data     !
    !               > 1 write only population data  !
    !-----------------------------------------------!
    integer, intent(in) :: flag
    real*8 tsteps     
    integer i, j, st
    j = 1
    ! tsteps = (/(i*dtc, i = 1, nsteps)/)
    
    open(112,file="trajStates.dat", status="replace")
    do i = 0, nsteps, nwrt_traj
      tsteps = i*dtc
      
      if ((flag==1) .and. (i .ne. 0)) then
        if (istr_sci==1) write(102,*) tsteps, sum_ci_sq(flag,j)
        if (istr_snk==1) write(103,*) tsteps, sum_nk(flag,j)
        if (istr_eng==1) write(104,*) tsteps, tot_energy(flag,j)
      endif

      ! write(112,*) tsteps, (/(pop(st,j-1), st=0,nbasis_e)/)
      if (iflow==0) then ! average population
        write(112,*)tsteps, summ_all(j-1),summ_nz(j-1),(/(pop(st,j-1), st=0,nbasis_e)/)   
      elseif (iflow==1) then ! population_sum
        write(112,*)tsteps, summ_all(j-1),summ_nz(j-1),(/(traj_states(st,j-1), st=0,nbasis_e)/)
      endif

      j = j+1
    enddo
    close(112)
  end subroutine writefiles
  !----------------------------------------------------------!
end module dynamics