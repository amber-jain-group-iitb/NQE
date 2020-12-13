module dynamics
    use model_spin_boson_q1
    implicit none

    integer     :: iflow, ncore, ntraj, seed
    real*8      :: dtc, ttime
    integer     :: ifriction, istate
    integer     :: iprint, istr_sci, istr_eng, wrt_qs, nwrt_traj
    integer     :: nsteps, nreadouts

    real*8,     allocatable :: H(:,:), H_new(:,:)
    real*8,     allocatable :: dV_dq(:,:,:)

    complex*16, allocatable :: ci(:)
    real*8,     allocatable :: q(:), qdot(:)
    real*8,     allocatable :: R(:), Rdot(:)
    real*8,     allocatable :: foRce(:)

    ! real*8,     allocatable :: tsteps(:)
    real*8,     allocatable :: sum_ci_sq(:,:), tot_energy(:,:)
    real*8,     allocatable :: pop(:), pop_sum(:), pop_avg(:)
     
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
        
        read(10,*) ifriction
        read(10,*) istate

        read(10,*) iprint
        read(10,*) nwrt_traj
        read(10,*) istr_sci
        read(10,*) istr_eng
        read(10,*) wrt_qs
        if (iflow == 1) read(10,*) seed
        close(10)

        nsteps      = ceiling(ttime/dtc)
        nreadouts   = ceiling(real(nsteps)/nwrt_traj)

        ! print *, iflow, ncore, ntraj
        ! print *, dtc, ttime, nsteps, ifriction, istate
        ! print *, iprint, istr_sci, istr_eng
        print *, nsteps, nreadouts

        ! allocate arrays
        allocate( H(nbasis_ev, nbasis_ev), H_new(nbasis_ev, nbasis_ev) )
        allocate( dV_dq(nbasis_ev, nbasis_ev, nclass_c) )

        allocate( ci(nbasis_ev) )
        allocate( q(nclass), qdot(nclass) )
        allocate( R(nclass_c), Rdot(nclass_c) )
        allocate( foRce(nclass_c) )

        ! allocate( tsteps(nsteps) )
        allocate( sum_ci_sq(1,nreadouts), tot_energy(1,nreadouts) )    ! (ntraj,nsteps)
        allocate( pop(nreadouts), pop_sum(nreadouts) )
    end subroutine setup_dynamics

!----------------------------------------------------!

    subroutine init_condition()
        implicit none
        integer :: i
        real*8, dimension(nclass) :: sig_q, sig_qdot
        real*8 :: rnd1, rnd2   
        
        do i = 1, nclass
            sig_q(i)    = 1.d0/sqrt(beta*mass(i)*omega(i)**2) 
            sig_qdot(i) = 1.d0/sqrt(beta*mass(i))

            call gauss_rand_n(rnd1)
            call gauss_rand_n(rnd2)
            ! print *, rnd1, rnd2
            
            q(i)       = -g(i)/(mass(i)*omega(i)**2) + rnd1*sig_q(i)
            qdot(i)    = rnd2*sig_qdot(i)
        enddo

        R(1)    = q(2)
        Rdot(1) = qdot(2)

        H       = ham_vibronic_diab(R)
        dV_dq   = delV_delq(R)

        ci          = (0.d0,0.d0)
        ci(istate)  = (1.d0,0.d0)
        call compute_force(ci, dV_dq, foRce)
    end subroutine init_condition

    !----------------------------------------------------------!

    subroutine evolve_classical(q, qdot, force) !! dummy variables for (R, Rdot, foRce; accR)
        ! Velocity Verlet
        implicit none
        real*8, intent(inout) :: q(:), qdot(:), force(:)
        real*8, dimension(nclass_c) :: acc 
        real*8, dimension(nclass) :: delr,delv
        real*8 c0, c1, c2, gdt

        if (ifriction == 0) then 
            acc = force/mass
            q   = q + qdot*dtc + 0.5*acc*(dtc**2)

            dV_dq   = delV_delq(q)
            call compute_force(ci, dV_dq, force)
    
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

            q   = q + c1*dtc*qdot + c2*(dtc**2)*acc + delr(2)

            dV_dq   = delV_delq(q)
            call compute_force(ci, dV_dq, force)

            acc     = (c1-c2)*acc + c2*(force/mass)
            qdot    = c0*qdot + acc*dtc + delv(2)
        endif

    end subroutine evolve_classical

    subroutine evolve_quantum(ci, H_k1, H_k4)
        ! Runge-Kutta 4th order integrator
        implicit none
        complex*16, intent(inout)   :: ci(:)
        real*8,     intent(in)      :: H_k1(:,:), H_k4(:,:) 
        complex*16, dimension(nbasis_ev)         :: k1, k2, k3, k4
        real*8,     dimension(nbasis_ev, nbasis_ev)   :: H_k23

        k1 = -iota/hbar*matmul(H_k1, ci)
        
        H_k23 = (H_k1 + H_k4)/2.d0
        k2 = -iota/hbar*matmul(H_k23, ci + dtc*k1/2.d0)
        k3 = -iota/hbar*matmul(H_k23, ci + dtc*k2/2.d0)
        
        k4 = -iota/hbar*matmul(H_k4, ci + dtc*k3)

        ci = ci + dtc*(k1+2*k2+2*k3+k4)/6.d0
    end subroutine evolve_quantum

    subroutine run_dynamics()
        implicit none
        integer ns, nt, i, j
        integer :: s(35)    ! seeder

        !call setup()   !! setup dynamics in main program
        s = 10
        if(iflow==1) s = s*seed
        print *, s
        call random_seed(put=s)
        ! initialize trajectory observables
        pop_sum     = 0.d0
        sum_ci_sq   = 0.d0
        tot_energy  = 0.d0
        do nt = 1, ntraj
            call init_condition()
            j = 1
            do ns = 1, nsteps
                if (iprint==1)  print *, nt, ns
                if (mod(ns, nwrt_traj)==1) then    
                    if (ntraj==1) then 
                        if (istr_sci==1) sum_ci_sq(nt,j)   = calc_ci_ss(ci)
                        if (istr_eng==1) tot_energy(nt,j)  = calc_tot_eng(ci, H)
                        if (wrt_qs==1)   write(101,*)   nt, ns, q(1), q(2), R
                    endif 
                    pop(j) = calc_pop(ci, istate)
                    j = j+1
                endif

                !----------evolve trajectory----------!
                call evolve_classical(R, Rdot, foRce)
                H_new   = ham_vibronic_diab(R)
                call evolve_quantum(ci, H, H_new)
                H       = H_new
                !-------------------------------------!
            enddo
            if (iprint==2) print *, nt, ns
            pop_sum = pop_sum + pop
        enddo
        call writefiles(ntraj)
    end subroutine run_dynamics

    !----------------------------------------------------------!
    subroutine compute_force(ci, dV_dq, force)      !! dummy variables for ( , , foRce)
        implicit none
        integer i
        complex*16, intent(in)      :: ci(:)
        real*8,     intent(in)      :: dV_dq(:,:,:)
        real*8,dimension(nclass_c), intent(out) ::  force(:)
        
        do i = 1, nclass_c
            !force(i) =- real(sum(conjg(ci)*matmul(delV_diab(:,:,i), ci)))
            force(i) = - expec_val(ci, dV_dq(:,:,i), ci)
        enddo
    end subroutine compute_force
    
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

    !----------------------------------------------------------!

    function calc_pop(ci, istate) result(popn)
        implicit none
        complex*16, intent(in) :: ci(:)
        integer,    intent(in) :: istate
        real*8 popn
        if (istate==1) popn = sum(abs(ci(1 : nbasis_v))**2)
        if (istate==2) popn = sum(abs(ci(nbasis_v+1 : 2*nbasis_v))**2)
    end function calc_pop

    ! subroutine avg_pop(pop)
    !     implicit none
    !     real*8,     intent(in):: pop(:,:)
    !     !real*8 pop_sum(nsteps), pop_avg(nsteps)
    !     pop_avg = pop_sum/ntraj
    ! end subroutine avg_pop

    function calc_ci_ss(ci) result(ci_ss)
        implicit none
        complex*16, intent(in)  :: ci(:)
        real*8 ci_ss
        ci_ss   = sum(abs(ci)**2)
    end function calc_ci_ss

    function calc_tot_eng(ci, H) result(tot_eng) 
        implicit none
        complex*16, intent(in)  :: ci(:)
        real*8,     intent(in)  :: H(:,:)
        real*8 tot_eng      !! electronic energy(PE) + KE
        tot_eng = expec_val(ci, H, ci) + (0.5*mass(2)*Rdot(1)**2)
    end function calc_tot_eng
    
    subroutine gauss_rand_n(rnd)    !! function
        implicit none
        real*8,intent(out)::rnd
        real*8 rnd1,rnd2
        
        call random_number(rnd1)
        call random_number(rnd2)
        rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)
    end subroutine gauss_rand_n

    !----------------------------------------------------------!

    subroutine writefiles(flag)
        implicit none
        !-----------------------------------------------!
        ! flag  = ntraj = 1 write a traj + pop data     !
        !               > 1 write only population data  !
        !-----------------------------------------------!
        integer, intent(in) :: flag
        real*8 tsteps     
        integer i, j
        j = 1
        ! tsteps = (/(i*dtc, i = 1, nsteps)/)
        do i = 0, nsteps, nwrt_traj
            tsteps = i*dtc
            
            if (flag==1) then
                if (istr_sci==1) write(102,*) tsteps, sum_ci_sq(flag,j)
                if (istr_eng==1) write(103,*) tsteps, tot_energy(flag,j)
            endif

            if (iflow==0)   write(100,*) tsteps, pop_sum(j)/ntraj   ! average population
            if (iflow==1)   write(100,*) tsteps, pop_sum(j)         ! population_sum
            j = j+1
        enddo
    end subroutine writefiles
    !----------------------------------------------------------!
!--------------------------------------------------------------------!
end module dynamics





