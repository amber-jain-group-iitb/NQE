program sqc_dynamics
    use model_spin_boson_q1
    use dynamics
    implicit none
    
    call setup_model
    call setup_ham_dvr() 

    !---------generate pes----------!
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
    ! call pes(1, 1)

    !---------run dynamics----------!  
    ! specify model to use in mod_dynamics_xxx
    call setup_dynamics()
    call run_dynamics()
    
end program sqc_dynamics
