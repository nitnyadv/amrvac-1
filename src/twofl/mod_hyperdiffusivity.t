!> Subroutines for Roe-type Riemann solver for HD
module mod_hyperdiffusivity


  implicit none
  private

  public :: hyperdiffusivity_init

contains

  subroutine hyperdiffusivity_init()
    use mod_global_parameters
    use mod_physics, only: phys_req_diagonal

    nghostcells = 3
    phys_req_diagonal = .true.

  end subroutine hyperdiffusivity_init


end module mod_hyperdiffusivity
