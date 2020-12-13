module dynamics
    use model_spin_boson_q0
    implicit none

    integer     :: iflow, ncore, ntraj, seed
    real*8      :: dtc, ttime
    integer     :: ifriction, istate
    integer     :: iprint, istr_sci, istr_eng
    integer     :: nsteps

    complex*16, allocatable :: ci(:)
    real*8,     allocatable :: q(:), qdot(:), acc(:)
    real*8,     allocatable :: force(:)

    real*8,     allocatable :: tsteps(:)
    real*8,     allocatable :: sum_ci_sq(:,:), tot_energy(:,:), pop(:,:)
    real*8,     allocatable :: pop_sum(:), pop_avg(:)
     
    contains
    subroutine setup()
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
        read(10,*) istr_sci
        read(10,*) istr_eng
        if (iflow == 1) read(10,*) seed
        close(10)

        nsteps = ceiling(ttime/dtc)

        ! print *, iflow, ncore, ntraj
        ! print *, dtc, ttime, nsteps, ifriction, istate
        ! print *, iprint, istr_sci, istr_eng

        ! allocate arrays
        allocate( ci(nquant) )
        allocate( q(nclass), qdot(nclass), acc(nclass) )
        allocate( force(nclass))

        allocate( tsteps(nsteps) )
        allocate( sum_ci_sq(ntraj,nsteps), tot_energy(ntraj,nsteps), pop(ntraj,nsteps) )
        allocate( pop_sum(nsteps), pop_avg(nsteps))
    end subroutine setup

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

        H       = ham_diab(q)
        dV_dq   = delV_delq(q)

        ci          = complex(0d0,0d0)
        ci(istate)  = complex(1d0,0d0)
        call compute_force(ci, dV_dq, force)
    end subroutine init_condition

    !----------------------------------------------------------!

    subroutine evolve_classical(q, qdot, force)
        ! Velocity Verlet
        implicit none
        real*8, intent(inout) :: q(:), qdot(:), force(:)
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

            q   = q + c1*dtc*qdot + c2*(dtc**2)*acc + delr

            dV_dq   = delV_delq(q)
            call compute_force(ci, dV_dq, force)

            acc     = (c1-c2)*acc + c2*(force/mass)
            qdot    = c0*qdot + acc*dtc + delv    
        endif

    end subroutine evolve_classical

    subroutine evolve_quantum(ci, H_k1, H_k4)
        ! Runge-Kutta 4th order integrator
        implicit none
        complex*16, intent(inout)   :: ci(:)
        real*8,     intent(in)      :: H_k1(:,:), H_k4(:,:) 
        complex*16, dimension(1:nquant)         :: k1, k2, k3, k4
        real*8,     dimension(nbasis, nbasis)   :: H_k23

        k1 = -iota/hbar*matmul(H_k1, ci)
        
        H_k23 = (H_k1 + H_k4)/2.d0
        k2 = -iota/hbar*matmul(H_k23, ci + dtc*k1/2.d0)
        k3 = -iota/hbar*matmul(H_k23, ci + dtc*k2/2.d0)
        
        k4 = -iota/hbar*matmul(H_k4, ci + dtc*k3)

        ci = ci + dtc*(k1+2*k2+2*k3+k4)/6.d0
    end subroutine evolve_quantum

    subroutine run_dynamics()
        implicit none
        integer ns, nt, i

        !call setup()
        ! initialize trajectory observables
        pop         = 0.d0
        sum_ci_sq   = 0.d0
        tot_energy  = 0.d0
        ! print *, (pop(1,i), i=1,10)

        do nt = 1, 1
            call init_condition()

            do ns = 1, nsteps
                if (iprint==1)      print *, nt, ns        ! print(q,q_dot)
                if (istr_sci==1)    call check_ci(ci, nt, ns)
                if (istr_eng==1)    call check_energy(ci, H, nt, ns)
                call calculate_population(ci, istate, nt, ns)

                ! evolve trajectory
                call evolve_classical(q, qdot, force)
                H_new   = ham_diab(q)
                call evolve_quantum(ci, H, H_new)
                H       = H_new 
                write(10,*) ns,  q(1), q(2)
            enddo
            if (iprint==2) print *, nt, ns
            call avg_pop(pop)
        enddo
        call writefiles()
    end subroutine run_dynamics

    !----------------------------------------------------------!
    subroutine compute_force(ci, dV_dq, force)
        implicit none
        integer i
        complex*16, intent(in)      :: ci(:)
        real*8,     intent(in)      :: dV_dq(:,:,:)
        real*8,dimension(nclass), intent(out) ::  force(:)
        
        do i = 1, nclass
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

    subroutine calculate_population(ci,  istate, nt, ns)
        implicit none
        complex*16, intent(in) :: ci(:)
        integer,    intent(in) :: istate, nt, ns
        !real*8 pop(ntraj,nsteps)

        pop(nt,ns) = abs(ci(istate))**2
    end subroutine calculate_population

    subroutine avg_pop(pop)
        implicit none
        real*8,     intent(in):: pop(:,:)
        !real*8 pop_sum(nsteps), pop_avg(nsteps)

        pop_sum = sum(pop, dim=1)
        pop_avg = pop_sum/ntraj
    end subroutine avg_pop

    subroutine check_ci(ci, nt, ns)
        implicit none
        complex*16, intent(in)  :: ci(:)
        integer,    intent(in)  :: nt, ns
        
        sum_ci_sq(nt,ns)  = sum(abs(ci)**2)
    end subroutine check_ci

    subroutine check_energy(ci, H, nt, ns) !! electronic energy(PE) + KE
        implicit none
        complex*16, intent(in)  :: ci(:)
        real*8,     intent(in)  :: H(:,:)
        integer,    intent(in)  :: nt, ns
        
        tot_energy(nt,ns) = expec_val(ci, H, ci) + sum(0.5*mass*qdot**2)
    end subroutine check_energy
    
    subroutine gauss_rand_n(rnd)    !! function
        implicit none
        real*8,intent(out)::rnd
        real*8 rnd1,rnd2
        
        call random_number(rnd1)
        call random_number(rnd2)
        rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)
    end subroutine gauss_rand_n

    !----------------------------------------------------------!

    subroutine writefiles()
        implicit none
        integer i, itraj
        itraj = 1   ! Analyse a trajectory
        
        ! tsteps = (/(i*dtc, i = 1, nsteps)/)
        do i = 1, nsteps
            tsteps(i) = i*dtc

            if (iflow==0)       write(100,*) tsteps(i), pop_avg(i)      ! average population
            if (iflow==1)       write(100,*) tsteps(i), pop_sum(i)      ! population_sum
            
            if (istr_sci==1)    write(102,*) tsteps(i), sum_ci_sq(itraj,i)
            if (istr_eng==1)    write(103,*) tsteps(i), tot_energy(itraj,i)
        enddo

    end subroutine writefiles
    !----------------------------------------------------------!
!--------------------------------------------------------------------!
end module dynamics





