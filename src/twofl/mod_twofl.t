module mod_twofl
  use mod_twofl_phys

  use mod_amrvac

  implicit none
  public

contains

  subroutine twofl_activate()
    call twofl_phys_init()
  end subroutine twofl_activate

end module mod_twofl
