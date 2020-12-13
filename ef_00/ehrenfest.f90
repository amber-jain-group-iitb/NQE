program ehrenfest_dynamics

    !! generate pes
    !-------------------------------!
    ! SUBROUTINE : pes(ipes, wrt_au)!
    ! write position vs energy value!
    !       in fort.200             !
    ! ipes  : 1 - diabatic pes      !
    !        11 - adiabatic pes     !
    !        12 - adiabatic pes (vc)!
    !         2 - pes gradient plot !
    ! wrt_au: 0 - angstrom vs cm_1  !
    !         1 - atomic units      !
    !-------------------------------!
    use model_spin_boson_q1 
    call setup_ham_dvr()  
    call pes(1, 0)

    ! !! run dynamics
    ! use dynamics    ! specify model to use in mod_dynamics_ef
    ! ! call setup_ham_dvr() 
    ! call setup()
    ! call run_dynamics()
    
end program ehrenfest_dynamics
