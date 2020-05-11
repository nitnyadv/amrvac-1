!> Module for including collisions
module mod_collisions
  implicit none



  double precision :: sigma_en

  contains


  !> Initialize the module
  subroutine collisions_init(SI_unit)
  !> Collision cross sections in both cgs and si
  double precision, parameter :: sigma_en_cgs = 1d-17 ! cm
  double precision, parameter :: sigma_en_si = 1d-19 ! m
    
    logical, intent(in) :: SI_unit

    write(*,*) "collisions init"
    if(SI_unit) then
      sigma_en = sigma_en_si

    else
      sigma_en = sigma_en_cgs

    end if
  end subroutine collisions_init


end module mod_collisions
