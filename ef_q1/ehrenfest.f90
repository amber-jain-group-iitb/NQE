program ehrenfest_dynamics

    !! generate pes
    !-------------------------------!
    ! SUBROUTINE : pes(ipes, wrt_au)!
    ! write position vs energy value!
    !       in fort.200             !
    ! ipes  : 1 - diabatic pes      !
    !         2 - adiabatic pes     !
    !         3 - pes gradient plot !
    ! wrt_au: 0 - angstrom vs cm_1  !
    !         1 - atomic units      !
    !-------------------------------!
    ! use model_spin_boson_q0 
    ! call setup_param() 
    ! call pes(2, 0)

    !! run dynamics : specify model to use in mod_dynamics_ef
    use dynamics
    call setup_param()     
    call setup_dynamics()
    call run_dynamics()
    
end program ehrenfest_dynamics
