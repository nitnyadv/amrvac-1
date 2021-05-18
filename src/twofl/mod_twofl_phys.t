!> Magneto-hydrodynamics module
module mod_twofl_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len
  use mod_physics
  implicit none
  private





  !! E_c = E_kin + E_mag + E_int
  !! E_n = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_TOT=2
  !! E_c = E_int
  !! E_n = E_int
  integer, public, parameter              :: EQ_ENERGY_INT=1
  !! E_n, E_c are calculated from density as c_adiab rho^gamma
  !! No energy equation => no variable assigned for it
  integer, public, parameter              :: EQ_ENERGY_NONE=0
  !! E_c = E_kin + E_int
  !! E_n = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_KI=3
  !! additional variable for the charges energy at index eaux_
  !! E_c (index e_) = E_kin + E_mag + E_int, E_c (index eaux_) = E_int
  !! E_n (index e_) = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_TOT2=4

  integer, public, protected              :: twofl_eq_energy = EQ_ENERGY_KI


  !> Whether thermal conduction is used
  logical, public, protected              :: twofl_thermal_conduction_c = .false.


  !> Whether radiative cooling is added
  logical, public, protected              :: twofl_radiative_cooling = .false.

  !> Whether viscosity is added
  logical, public, protected              :: twofl_viscosity = .false.

  !> Whether gravity is added: common flag for charges and neutrals
  logical, public, protected              :: twofl_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: twofl_Hall = .false.
#if defined(ONE_FLUID) && ONE_FLUID==1
  !> Whether Ambipolar term is used
  logical, public, protected              :: twofl_ambipolar = .false.

  !> Whether Ambipolar term is implemented using supertimestepping
  logical, public, protected              :: twofl_ambipolar_sts = .false.

  !> Whether Ambipolar term is implemented explicitly
  logical, public, protected              :: twofl_ambipolar_exp = .false.

  !> The MHD ambipolar coefficient
  double precision, public                :: twofl_eta_ambi = 0.0d0
#endif

  !> Whether TRAC method is used
  logical, public, protected              :: twofl_trac = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: twofl_glm = .false.

  !> Which TRAC method is used            
  integer, public, protected              :: twofl_trac_type=1

  !> Height of the mask used in the TRAC method
  double precision, public, protected              :: twofl_trac_mask = 0.d0


  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: twofl_glm_alpha = 0.5d0

  !> Use Boris approximation
  character(len=20) :: twofl_boris_method = "none"

  integer, parameter :: boris_none           = 0
  integer, parameter :: boris_reduced_force  = 1
  integer, parameter :: boris_simplification = 2
  integer            :: twofl_boris_type       = boris_none

  !> Speed of light for Boris' approximation. If negative, test changes to the
  !> momentum equation with gamma_A = 1
  double precision                        :: twofl_boris_c = 0.0d0

  !> MHD fourth order
  logical, public, protected              :: twofl_4th_order = .false.

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_c_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom_c(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_c_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_c_
  integer, public, protected              :: Tweight_c_

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of auxiliary internal energy
  integer, public, protected :: eaux_c_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> equi vars flags
  logical, public, protected :: has_equi_rho_c0 = .false.  
  logical, public, protected :: has_equi_pe_c0 = .false.  

  !> equi vars indices in the state%equi_vars array
  integer, public, protected :: equi_rho_c0_ = -1
  integer, public, protected :: equi_pe_c0_ = -1

#if !defined(ONE_FLUID) || ONE_FLUID==0
  !neutrals:
  logical, public, protected              :: twofl_thermal_conduction_n = .false.

  integer, public, protected              :: rho_n_
  integer, allocatable, public, protected :: mom_n(:)
  integer, public, protected              :: e_n_
  integer, public, protected              :: Tcoff_n_
  integer, public, protected              :: Tweight_n_
  logical, public, protected :: has_equi_rho_n0 = .false. 
  logical, public, protected :: has_equi_pe_n0 = .false.  
  integer, public, protected :: equi_rho_n0_ = -1
  integer, public, protected :: equi_pe_n0_ = -1

  ! related to collisions:
  !> collisional alpha
  double precision, public                :: twofl_alpha_coll = 0d0
  !> whether include thermal exchange collisional terms
  logical, public                         :: twofl_coll_inc_te = .true.
  logical, public                         :: twofl_implicit_coll_terms = .true.
  double precision, public                :: dtcollpar = -1d0 !negative value does not impose restriction on the timestep

#endif

  ! Eq of state:
  !> Helium abundance over Hydrogen
  !> He_abundance = (He2+ + He+ + He)/(H+ + H)
  double precision, public, protected  :: He_abundance=0.1d0
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  !> H_ion_fr is not allowed to be 0
  double precision, public, protected  :: H_ion_fr=1d0
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected  :: He_ion_fr=1d0
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected  :: He_ion_fr2=1d0
  double precision, public, protected  :: Rc  ! defined for compat with the new eq of state, it is set to 1 for ONE_FLUID
#if !defined(ONE_FLUID) || ONE_FLUID==0
  double precision, public, protected  :: Rn,rho_nc_fr
#endif
  ! Eq of state end

  !> The adiabatic index
  double precision, public                :: twofl_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: twofl_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: twofl_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: twofl_eta_hyper = 0.0d0

  !> The MHD Hall coefficient
  double precision, public                :: twofl_etah = 0.0d0


  !> The small_est allowed energy
  double precision, protected             :: small_e


  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: twofl_divb_4thorder = .false.

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'


  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.



  !> added from modules: gravity
  !> source split or not
  logical :: grav_split= .false.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_multigrid     = -1
  integer, parameter :: divb_glm           = 1
  integer, parameter :: divb_powel         = 2
  integer, parameter :: divb_janhunen      = 3
  integer, parameter :: divb_linde         = 4
  integer, parameter :: divb_lindejanhunen = 5
  integer, parameter :: divb_lindepowel    = 6
  integer, parameter :: divb_lindeglm      = 7
  integer, parameter :: divb_ct            = 8


  ! Public methods
  public :: twofl_phys_init
  public :: twofl_to_conserved
  public :: twofl_to_primitive
  public :: get_divb
  public :: get_rhoc_tot
#if !defined(ONE_FLUID) || ONE_FLUID==0
  public :: get_rhon_tot
#endif
  public :: get_current
  public :: twofl_get_pthermal_c
  public :: twofl_get_csound2
  public :: get_normalized_divb
  public :: b_from_vector_potential
  {^NOONED
  public :: twofl_clean_divb_multigrid
  }

#if defined(ONE_FLUID) && ONE_FLUID==1
  !define the subroutine interface for the ambipolar mask
  abstract interface

    subroutine mask_subroutine(ixI^L,ixO^L,w,x,res)
       use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(inout) :: res(ixI^S)
    end subroutine mask_subroutine

  end interface

   procedure (mask_subroutine), pointer :: usr_mask_ambipolar => null()
   public :: usr_mask_ambipolar 

#endif
contains

  !> Read this module"s parameters from a file
  subroutine twofl_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /twofl_list/ twofl_eq_energy, twofl_gamma, twofl_adiab,&
      twofl_eta, twofl_eta_hyper, twofl_etah, twofl_glm_alpha,& 
      twofl_thermal_conduction_c, twofl_radiative_cooling, twofl_Hall, twofl_gravity,&
      twofl_viscosity, twofl_4th_order, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, type_ct, divbwave,He_abundance, SI_unit, B0field,&
      B0field_forcefree, Bdip, Bquad, Boct, Busr,&
      !added:
      has_equi_rho_c0, has_equi_pe_c0,&
      H_ion_fr, He_ion_fr, He_ion_fr2,&
#if !defined(ONE_FLUID) || ONE_FLUID==0
      has_equi_pe_n0, has_equi_rho_n0, twofl_thermal_conduction_n,  &
      twofl_alpha_coll, twofl_implicit_coll_terms, twofl_coll_inc_te, dtcollpar,&
#else
      twofl_ambipolar, twofl_ambipolar_sts, twofl_eta_ambi,&
#endif
      !added end
      has_equi_rho_c0, has_equi_pe_c0,&
      boundary_divbfix, boundary_divbfix_skip, twofl_divb_4thorder, &
      twofl_boris_method, twofl_boris_c, clean_initial_divb,  &
      twofl_trac, twofl_trac_type, twofl_trac_mask

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, twofl_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine twofl_read_params

  !> Write this module's parameters to a snapsoht
  subroutine twofl_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = twofl_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine twofl_write_info

  subroutine twofl_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim),  wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
    integer                            :: hxO^L, kxC^L, iw
    double precision                   :: inv_volume(ixI^S)

    call mpistop("to do")

  end subroutine twofl_angmomfix

  subroutine twofl_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    !use mod_gravity, only: gravity_init
    use mod_supertimestepping, only: sts_init, add_sts_method

    {^NOONED
    use mod_multigrid_coupling
    }

    integer :: itr, idir

    call twofl_read_params(par_files)
    if(mype .eq. 0) then
#if defined(ONE_FLUID) && ONE_FLUID==1
      print*, "ONE FLUID SET IN amrvac.h"
#else
      print*, "TWO FLUID VERSION"
#endif
    endif
    physics_type = "twofl"
    phys_energy=.true.
  !> Solve total energy equation or not
  ! for the two fluid the true value means 
  ! E_charges = E_mag + E_kin_charges + E_int_charges
  ! E_neutrals =  E_kin_neutrals + E_int_neutrals
    phys_total_energy=.false.

  !> Solve internal enery instead of total energy
  ! for the two fluid the true vale means 
  ! E_charges = E_int_charges
  ! E_neutrals = E_int_neutrals
    phys_internal_e=.false.

  ! For the two fluid phys_energy=.true. and phys_internal_e=.false. means
  ! E_charges = E_kin_charges + E_int_charges
  ! E_neutrals =  E_kin_neutrals + E_int_neutrals




  !> Solve internal energy and total energy equations
  ! this implies two equations of energy solved
    phys_solve_eaux=.false.


    if(twofl_eq_energy == EQ_ENERGY_INT) then
      phys_internal_e = .true.
    elseif(twofl_eq_energy == EQ_ENERGY_TOT .or. twofl_eq_energy == EQ_ENERGY_TOT2) then
      phys_total_energy = .true.
      if(twofl_eq_energy == EQ_ENERGY_TOT2) then
        phys_solve_eaux = .true.
      endif

    elseif(twofl_eq_energy == EQ_ENERGY_NONE) then
      phys_energy = .false.
    endif

    phys_trac=twofl_trac
    phys_trac_type=twofl_trac_type



    if(.not. phys_energy) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
      if(twofl_thermal_conduction_n) then
        twofl_thermal_conduction_n=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_thermal_conduction_n=F when twofl_energy=F'
      end if
#endif
      if(twofl_thermal_conduction_c) then
        twofl_thermal_conduction_c=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_thermal_conduction_c=F when twofl_energy=F'
      end if
      if(twofl_radiative_cooling) then
        twofl_radiative_cooling=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_radiative_cooling=F when twofl_energy=F'
      end if
      if(twofl_trac) then
        twofl_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_trac=F when twofl_energy=F'
      end if
    end if
    {^IFONED
      if(twofl_trac .and. twofl_trac_type .gt. 1) then
        twofl_trac_type=1
        if(mype==0) write(*,*) 'WARNING: set twofl_trac_type=1 for 1D simulation'
      end if
    }
    if(twofl_trac .and. twofl_trac_type .le. 3) then
      twofl_trac_mask=bigdouble
      if(mype==0) write(*,*) 'WARNING: set twofl_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=twofl_trac_mask

    if(phys_solve_eaux) prolongprimitive=.true.

    ! set default gamma for polytropic/isothermal process
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    {^NOONED
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => twofl_clean_divb_multigrid
    }
    case ('glm')
      twofl_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm
    case ('powel', 'powell')
      type_divb = divb_powel
    case ('janhunen')
      type_divb = divb_janhunen
    case ('linde')
      type_divb = divb_linde
    case ('lindejanhunen')
      type_divb = divb_lindejanhunen
    case ('lindepowel')
      type_divb = divb_lindepowel
    case ('lindeglm')
      twofl_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    select case (twofl_boris_method)
    case ("none")
      twofl_boris_type = boris_none
    case ("reduced_force")
      twofl_boris_type = boris_reduced_force
    case ("simplification")
      twofl_boris_type = boris_simplification
    case default
      call mpistop("Unknown twofl_boris_method (none, reduced_force, simplification)")
    end select

    !allocate charges first and the same order as in mhd module
    
    rho_c_ = var_set_fluxvar("rho_c", "rho_c")
    !set variables from mod_variables to point to charges vars
    iw_rho = rho_c_

    allocate(mom_c(ndir))
    do idir=1,ndir
      mom_c(idir) = var_set_fluxvar("m_c","v_c",idir)
    enddo

    allocate(iw_mom(ndir))
    iw_mom(1:ndir) = mom_c(1:ndir)

    ! Set index of energy variable
    if (phys_energy) then
      e_c_ = var_set_fluxvar("e_c", "p_c")
      iw_e = e_c_
    else
      e_c_     = -1
    end if


    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

  !now allocate neutrals
  ! ambipolar sts assumes mag and energy charges are continuous

#if !defined(ONE_FLUID) || ONE_FLUID==0

    ! Determine flux variables
    rho_n_ = var_set_fluxvar("rho_n", "rho_n")
    allocate(mom_n(ndir))
    do idir=1,ndir
      mom_n(idir) = var_set_fluxvar("m_n","v_n",idir)
    enddo
    if (phys_energy) then
      e_n_ = var_set_fluxvar("e_n", "p_n")
    else
      e_n_     = -1
    end if

    ! check dtcoll par
    if(.not. twofl_implicit_coll_terms .and. (dtcollpar .le. 0d0 .or. dtcollpar .ge. 1d0)) then
      if (mype .eq. 0) print*, "Explicit update of coll terms requires 0<dtcollpar<1, dtcollpar set to 0.8."
      dtcollpar = 0.8
    endif 

#endif



    if (twofl_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    !  set auxiliary internal energy variable
    if(phys_energy .and. phys_solve_eaux) then
      eaux_c_ = var_set_fluxvar("eaux_c", "paux_c",need_bc=.false.)
      iw_eaux = eaux_c_
    else
      eaux_c_ = -1
    end if


    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_c_ = -1
#if !defined(ONE_FLUID) || ONE_FLUID==0
    Tweight_n_ = -1
#endif
    if(twofl_trac) then
      Tcoff_c_ = var_set_fluxvar('Tcoff_c', 'TCoff_c', need_bc=.false.)
      iw_tcoff = Tcoff_c_
#if !defined(ONE_FLUID) || ONE_FLUID==0
      Tcoff_n_ = var_set_fluxvar('Tcoff_n', 'TCoff_n', need_bc=.false.)
#endif
      if(twofl_trac_type .ge. 2) then
        Tweight_c_ = var_set_extravar('Tweight_c', 'Tweight_c')
#if !defined(ONE_FLUID) || ONE_FLUID==0
        Tweight_n_ = var_set_extravar('Tweight_n', 'Tweight_n')
#endif
      endif
    else
#if !defined(ONE_FLUID) || ONE_FLUID==0
      Tcoff_n_ = -1
#endif
      Tcoff_c_ = -1
    end if

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
#if !defined(ONE_FLUID) || ONE_FLUID==0
    if(has_equi_rho_n0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho_n0_ = number_equi_vars
    endif  
    if(has_equi_pe_n0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe_n0_ = number_equi_vars
    endif  
#endif
    if(has_equi_rho_c0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho_c0_ = number_equi_vars
    endif  
    if(has_equi_pe_c0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe_c0_ = number_equi_vars
    endif  
      


    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! determine number of stagger variables
    if(stagger_grid) nws=ndim


    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(ndim>1) then
      if(twofl_glm) then
        flux_type(:,psi_)=flux_special
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_special
        end do
      else
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_tvdlf
        end do
      end if
    end if

    select case (twofl_boris_method)
    case ("none")
      twofl_boris_type = boris_none
    case ("reduced_force")
      twofl_boris_type = boris_reduced_force
    case ("simplification")
      twofl_boris_type = boris_simplification
      do idir = 1, ndir
        phys_iw_methods(mom_c(idir))%inv_capacity => twofl_gamma2_alfven
      end do
    case default
      call mpistop("Unknown twofl_boris_method (none, reduced_force, simplification)")
    end select

    phys_get_dt              => twofl_get_dt
    phys_get_cmax            => twofl_get_cmax
    phys_get_a2max           => twofl_get_a2max
    !phys_get_tcutoff         => twofl_get_tcutoff
    phys_get_cbounds         => twofl_get_cbounds
    phys_get_flux            => twofl_get_flux
    phys_add_source_geom     => twofl_add_source_geom
    phys_add_source          => twofl_add_source
    phys_to_conserved        => twofl_to_conserved
    phys_to_primitive        => twofl_to_primitive
    !phys_ei_to_e             => twofl_ei_to_e
    !phys_e_to_ei             => twofl_e_to_ei
    phys_check_params        => twofl_check_params
    phys_check_w             => twofl_check_w
    !phys_get_pthermal        => twofl_get_pthermal
    phys_write_info          => twofl_write_info
    phys_angmomfix           => twofl_angmomfix
    phys_handle_small_values => twofl_handle_small_values
    phys_energy_synchro      => twofl_energy_synchro
    ! implicit collisional terms update
#if !defined(ONE_FLUID) || ONE_FLUID==0
    if(twofl_implicit_coll_terms .and. twofl_alpha_coll>0d0) then
      phys_implicit_update => twofl_implicit_coll_terms_update
    endif
#endif
    !set equilibrium variables for the new grid
    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif

    if(type_divb==divb_glm) then
      phys_modify_wLR => twofl_modify_wLR
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_update_faces => twofl_update_faces
      phys_face_to_center => twofl_face_to_center
      phys_modify_wLR => twofl_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => twofl_boundary_adjust
    end if

    {^NOONED
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => twofl_clean_divb_multigrid
    }

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call twofl_physical_units()

    if(.not. phys_energy .and. (twofl_thermal_conduction_c& 
#if !defined(ONE_FLUID) || ONE_FLUID==0
        .or. twofl_thermal_conduction_n&
#endif
    )) then
      call mpistop("thermal conduction needs twofl_energy=T")
    end if
    if(.not. phys_energy .and. twofl_radiative_cooling) then
      call mpistop("radiative cooling needs twofl_energy=T")
    end if

    ! initialize thermal conduction module
!!TODO
!    if (twofl_thermal_conduction) then
!      phys_req_diagonal = .true.
!      if(use_twofl_tc .eq. MHD_TC) then
!
!        if(twofl_internal_e) then
!          call tc_init_twofl_for_internal_energy(twofl_gamma,[rho_,e_,mag(1)],twofl_get_temperature_from_eint)
!        else
!          if(twofl_solve_eaux) then
!            call tc_init_twofl_for_total_energy(twofl_gamma,[rho_,e_,mag(1),eaux_],twofl_get_temperature_from_etot, twofl_get_temperature_from_eint)
!          else
!            call tc_init_twofl_for_total_energy(twofl_gamma,[rho_,e_,mag(1)],twofl_get_temperature_from_etot, twofl_get_temperature_from_eint)
!          endif
!        endif
!
!      else if(use_twofl_tc .eq. HD_TC) then
!        if(twofl_internal_e) then
!          call tc_init_hd_for_internal_energy(twofl_gamma,[rho_,e_],twofl_get_temperature_from_eint)
!        else
!          if(twofl_solve_eaux) then
!            call tc_init_hd_for_total_energy(twofl_gamma,[rho_,e_,eaux_],twofl_get_temperature_from_etot, twofl_get_temperature_from_eint)
!          else
!            call tc_init_hd_for_total_energy(twofl_gamma,[rho_,e_],twofl_get_temperature_from_etot, twofl_get_temperature_from_eint)
!          endif
!        endif
!      endif
!    end if

    ! Initialize radiative cooling module
    !TODO
    !if (twofl_radiative_cooling) then
    !  call radiative_cooling_init(twofl_gamma,He_abundance)
    !end if

    ! Initialize viscosity module
    !!TODO
    !if (twofl_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(twofl_gravity) then
    !  call gravity_init()
       call grav_params_read(par_files)
    end if

    ! Initialize particles module


    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (twofl_hall) then
       phys_req_diagonal = .true.
       if (twofl_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

#if defined(ONE_FLUID) && ONE_FLUID==1
    if(twofl_ambipolar) then
      phys_req_diagonal = .true.
      if(twofl_ambipolar_sts) then
        call sts_init()
        !!ADDED  
!        if(twofl_4th_order) then
!          phys_wider_stencil = 2
!        else
!          phys_wider_stencil = 1
!        end if
        !!ADDED end 
        if(phys_internal_e) then
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mag(1),&
               ndir,mag(1),ndir,.true.)
        else
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mom_c(ndir)+1,&
               mag(ndir)-mom_c(ndir),mag(1),ndir,.true.)
        end if
      else
        twofl_ambipolar_exp=.true.
        ! For flux ambipolar term, we need one more reconstructed layer since currents are computed
        ! in mhd_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
        ! added in nghostcells.
        if(twofl_4th_order) then
          phys_wider_stencil = 2
        else
          phys_wider_stencil = 1
        end if
      end if
    end if
#endif
  end subroutine twofl_phys_init

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid_faces(igrid,x,ixI^L,ixO^L)
    use mod_global_parameters
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in) :: igrid, ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)

    double precision :: delx(ixI^S,1:ndim)
    double precision :: xC(ixI^S,1:ndim),xshift^D
    integer :: idims, ixC^L, hxO^L, ix, idims2

    if(slab_uniform)then
      ^D&delx(ixI^S,^D)=rnode(rpdx^D_,igrid)\
    else
      ! for all non-cartesian and stretched cartesian coordinates
      delx(ixI^S,1:ndim)=ps(igrid)%dx(ixI^S,1:ndim)
    endif
  
  
    do idims=1,ndim
      hxO^L=ixO^L-kr(idims,^D);
      if(stagger_grid) then
        ! ct needs all transverse cells
        ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
      else
        ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
        ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
      end if
      ! always xshift=0 or 1/2
      xshift^D=half*(one-kr(^D,idims));
      do idims2=1,ndim
        select case(idims2)
        {case(^D)
          do ix = ixC^LIM^D
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ix^D%ixC^S,^D)=x(ix^D%ixC^S,^D)+(half-xshift^D)*delx(ix^D%ixC^S,^D)
          end do\}
        end select
      end do
      call usr_set_equi_vars(ixI^L,ixC^L,xC,ps(igrid)%equi_vars(ixI^S,1:number_equi_vars,idims))
    end do

  end subroutine set_equi_vars_grid_faces


  !> sets the equilibrium variables
  subroutine set_equi_vars_grid(igrid)
    use mod_global_parameters
    use mod_usr_methods
  
    integer, intent(in) :: igrid

    !values at the center
    call usr_set_equi_vars(ixG^LL,ixG^LL,ps(igrid)%x,ps(igrid)%equi_vars(ixG^T,1:number_equi_vars,0))

    !values at the interfaces
    call set_equi_vars_grid_faces(igrid,ps(igrid)%x,ixG^LL,ixM^LL)
  
 
 end subroutine set_equi_vars_grid


  !> copied from mod_gravity
  subroutine grav_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grav_list/ grav_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grav_list, end=111)
111    close(unitpar)
    end do

  end subroutine grav_params_read

  subroutine twofl_check_params
    use mod_global_parameters
    use mod_usr_methods

    ! after user parameter setting
    gamma_1=twofl_gamma-1.d0
    if (.not. phys_energy) then
       if (twofl_gamma <= 0.0d0) call mpistop ("Error: twofl_gamma <= 0")
       if (twofl_adiab < 0.0d0) call mpistop ("Error: twofl_adiab < 0")
       small_pressure = twofl_adiab*small_density**twofl_gamma
    else
       if (twofl_gamma <= 0.0d0 .or. twofl_gamma == 1.0d0) &
            call mpistop ("Error: twofl_gamma <= 0 or twofl_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    if (twofl_boris_type > 0 .and. abs(twofl_boris_c) <= 0.0d0) then
      call mpistop("You have not specified twofl_boris_c")
    end if

!    if(H_ion_fr == 0d0 .and. He_ion_fr == 0d0) then
!      call mpistop("H_ion_fr or He_ion_fr must be > 0 or use hd module")
!    endif
!    if(H_ion_fr == 1d0 .and. He_ion_fr == 1d0) then
!      call mpistop("H_ion_fr or He_ion_fr must be < 1 or use mhd module")
!    endif
    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
  end subroutine twofl_check_params

  subroutine twofl_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0
    double precision :: a,b,c,d
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
    end if


#if !defined(ONE_FLUID) || ONE_FLUID==0
    a =  4d0  * He_abundance * He_ion_fr/H_ion_fr + 1d0 !rho_c 
    b = (2d0 + He_ion_fr2) * He_abundance * He_ion_fr/H_ion_fr + 2d0 !pe_c
    c = (1d0 - H_ion_fr + 4*He_abundance*(1d0 - He_ion_fr))/H_ion_fr !rho_n
    d = (1d0 - H_ion_fr + He_abundance*(1d0 - He_ion_fr))/H_ion_fr !pe_n
    Rc=b/a
    Rn=d/c
    rho_nc_fr = c/a
    if(mype .eq.0) then
      print*, "eq state Rn=", Rn, " Rc=",Rc
    endif
#else
    !a = a+c from above
    !b = b+d
    a = 1d0 + He_abundance*(1d0 + 3d0 * He_ion_fr)
    b = 1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0)
    Rc = 1d0 !for compatibility between single fluid eq state defined by units only and  
    ! twofl eq of state p = rho R T
#endif

    if(unit_velocity==0) then
#if defined(ONE_FLUID) && ONE_FLUID==1
      unit_density=a*mp*unit_numberdensity
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
#else
      unit_density=mp*unit_numberdensity
      unit_pressure=unit_numberdensity*kB*unit_temperature
#endif
      unit_velocity=sqrt(unit_pressure/unit_density)
    else

#if defined(ONE_FLUID) && ONE_FLUID==1
      unit_density=a*mp*unit_numberdensity
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
#else
      unit_density=mp*unit_numberdensity
      unit_temperature=unit_pressure/(unit_numberdensity*kB)
#endif
      unit_pressure=unit_density*unit_velocity**2
    end if
    unit_magneticfield=sqrt(miu0*unit_pressure)
    unit_time=unit_length/unit_velocity

  end subroutine twofl_physical_units

  subroutine twofl_check_w(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision :: tmp(ixI^S)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    flag=.false.
    !TODO 
    return
#if !defined(ONE_FLUID) || ONE_FLUID==0
        
    call get_rhon_tot(w,ixI^L,ixO^L,tmp)
    where(tmp(ixO^S) < small_density) flag(ixO^S,rho_n_) = .true.
#endif
    call get_rhoc_tot(w,ixI^L,ixO^L,tmp)
    where(tmp(ixO^S) < small_density) flag(ixO^S,rho_c_) = .true.
    if(phys_energy) then
      if(primitive) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
        tmp(ixO^S) = w(ixO^S,e_n_)
        if(has_equi_pe_n0) then
          tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_n0_,0)*gamma_1
        endif
        where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_n_) = .true.
#endif
        tmp(ixO^S) = w(ixO^S,e_c_)
        if(has_equi_pe_c0) then
          tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_c0_,0)*gamma_1
        endif
        where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_c_) = .true.
        if(twofl_eq_energy == EQ_ENERGY_TOT2) then 
          where(w(ixO^S,eaux_c_) < small_pressure) flag(ixO^S,e_c_) = .true.
        endif
      else
        if(phys_internal_e) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
          where(w(ixO^S,e_n_) < small_e) flag(ixO^S,e_n_) = .true.
#endif
          where(w(ixO^S,e_c_) < small_e) flag(ixO^S,e_c_) = .true.
        else
#if !defined(ONE_FLUID) || ONE_FLUID==0
          !neutrals
          tmp(ixO^S)=w(ixO^S,e_n_)-&
                twofl_kin_en_n(w,ixI^L,ixO^L)
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_n_) = .true.
#endif
          if(phys_total_energy) then
            tmp(ixO^S)=w(ixO^S,e_c_)-&
                twofl_kin_en_c(w,ixI^L,ixO^L)-twofl_mag_en(w,ixI^L,ixO^L)
          else
            tmp(ixO^S)=w(ixO^S,e_c_)-&
                twofl_kin_en_c(w,ixI^L,ixO^L)
          end if
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_c_) = .true.
          if(twofl_eq_energy == EQ_ENERGY_TOT2) then 
            where(w(ixO^S,eaux_c_) < small_e) flag(ixO^S,e_c_) = .true.
          endif
        end if
      endif
    end if

  end subroutine twofl_check_w

  !> Transform primitive variables into conservative ones
  subroutine twofl_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir
    double precision                :: rhoc(ixI^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision                :: rhon(ixI^S)
#endif

    if (fix_small_values) then
      call twofl_handle_small_values(.true., w, x, ixI^L, ixO^L, 'twofl_to_conserved')
    end if

#if !defined(ONE_FLUID) || ONE_FLUID==0
    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
#endif
    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(phys_energy) then
      if(phys_internal_e) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*inv_gamma_1
#endif
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*inv_gamma_1
      else
#if !defined(ONE_FLUID) || ONE_FLUID==0
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom_n(:))**2,dim=ndim+1)*rhon(ixO^S)
#endif
        if(phys_total_energy) then
          w(ixO^S,e_c_)=w(ixO^S,e_c_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom_c(:))**2,dim=ndim+1)*rhoc(ixO^S)&
                   +twofl_mag_en(w, ixI^L, ixO^L)
          if(twofl_eq_energy == EQ_ENERGY_TOT2) then
            w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)*inv_gamma_1
          endif
        else
          ! kinetic energy + internal energy is evolved   
          w(ixO^S,e_c_)=w(ixO^S,e_c_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom_c(:))**2,dim=ndim+1)*rhoc(ixO^S)
        endif
      end if
      !print*, "TOCONS ec ", w(1:10,e_c_)
      !print*, "TOCONS en ", w(1:10,e_n_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
#if !defined(ONE_FLUID) || ONE_FLUID==0
       w(ixO^S, mom_n(idir)) = rhon(ixO^S) * w(ixO^S, mom_n(idir))
#endif
       w(ixO^S, mom_c(idir)) = rhoc(ixO^S) * w(ixO^S, mom_c(idir))
    end do
  end subroutine twofl_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine twofl_to_primitive(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir
    double precision                :: rhoc(ixI^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision                :: rhon(ixI^S)
#endif


    if (fix_small_values) then
      call twofl_handle_small_values(.false., w, x, ixI^L, ixO^L, 'twofl_to_primitive')
    end if

#if !defined(ONE_FLUID) || ONE_FLUID==0
    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
#endif
    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)

    if(phys_energy) then
      if(phys_internal_e) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
#endif
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
      else
        ! neutrals evolved energy = ke + e_int 
#if !defined(ONE_FLUID) || ONE_FLUID==0
        w(ixO^S,e_n_)=gamma_1*(w(ixO^S,e_n_)&
                    -twofl_kin_en_n(w,ixI^L,ixO^L))
#endif
        ! charges
        if(phys_total_energy) then
         ! evolved energy = ke + e_int + e_mag 
          w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                    -twofl_kin_en_c(w,ixI^L,ixO^L)&
                    -twofl_mag_en(w,ixI^L,ixO^L))
          if(twofl_eq_energy == EQ_ENERGY_TOT2) then
            w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)*gamma_1
          endif
        else
         ! evolved energy = ke + e_int 
          w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                    -twofl_kin_en_c(w,ixI^L,ixO^L))
        end if
      end if
    end if
    
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/rhoc(ixO^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
       w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/rhon(ixO^S)
#endif
    end do

  end subroutine twofl_to_primitive



!!USED IN TC
!  !> Transform internal energy to total energy
!  subroutine twofl_ei_to_e(ixI^L,ixO^L,w,x)
!    use mod_global_parameters
!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision, intent(inout) :: w(ixI^S, nw)
!    double precision, intent(in)    :: x(ixI^S, 1:ndim)
! 
!    ! Calculate total energy from internal, kinetic and magnetic energy
!    if(twofl_solve_eaux) w(ixI^S,eaux_)=w(ixI^S,e_)
!    w(ixO^S,e_)=w(ixO^S,e_)&
!               +twofl_kin_en(w,ixI^L,ixO^L)&
!               +twofl_mag_en(w,ixI^L,ixO^L)
!
!  end subroutine twofl_ei_to_e
!
!  !> Transform total energy to internal energy
!  subroutine twofl_e_to_ei(ixI^L,ixO^L,w,x)
!    use mod_global_parameters
!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision, intent(inout) :: w(ixI^S, nw)
!    double precision, intent(in)    :: x(ixI^S, 1:ndim)
!
!    ! Calculate ei = e - ek - eb
!    w(ixO^S,e_)=w(ixO^S,e_)&
!                -twofl_kin_en(w,ixI^L,ixO^L)&
!                -twofl_mag_en(w,ixI^L,ixO^L)
!
!  end subroutine twofl_e_to_ei
!!USED IN TC END

  subroutine twofl_energy_synchro(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth1(ixI^S),pth2(ixI^S),alfa(ixI^S),beta(ixI^S)
    double precision, parameter :: beta_low=0.005d0,beta_high=0.05d0

!    double precision :: vtot(ixI^S),cs2(ixI^S),mach(ixI^S)
!    double precision, parameter :: mach_low=20.d0,mach_high=200.d0

    ! get magnetic energy
    alfa(ixO^S)=twofl_mag_en(w,ixI^L,ixO^L)
    pth1(ixO^S)=gamma_1*(w(ixO^S,e_c_)-twofl_kin_en_c(w,ixI^L,ixO^L)-alfa(ixO^S))
    pth2(ixO^S)=w(ixO^S,eaux_c_)*gamma_1
    ! get plasma beta
    beta(ixO^S)=min(pth1(ixO^S),pth2(ixO^S))/alfa(ixO^S)

    ! whether Mach number should be another criterion ?
!    vtot(ixO^S)=sum(w(ixO^S,mom(:))**2,dim=ndim+1)
!    call twofl_get_csound2(w,x,ixI^L,ixO^L,cs2)
!    mach(ixO^S)=sqrt(vtot(ixO^S)/cs2(ixO^S))/w(ixO^S,rho_)
    where(beta(ixO^S) .ge. beta_high)
!    where(beta(ixO^S) .ge. beta_high .and. mach(ixO^S) .le. mach_low)
      w(ixO^S,eaux_c_)=pth1(ixO^S)*inv_gamma_1
    else where(beta(ixO^S) .le. beta_low)
!    else where(beta(ixO^S) .le. beta_low .or. mach(ixO^S) .ge. mach_high)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)-pth1(ixO^S)*inv_gamma_1+w(ixO^S,eaux_c_)
    else where
      alfa(ixO^S)=dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low)
!      alfa(ixO^S)=min(dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low),
!                      dlog(mach_high(ixO^S)/mach(ixO^S))/dlog(mach_high/mach_low))
      w(ixO^S,eaux_c_)=(pth2(ixO^S)*(one-alfa(ixO^S))&
                     +pth1(ixO^S)*alfa(ixO^S))*inv_gamma_1
      w(ixO^S,e_c_)=w(ixO^S,e_c_)-pth1(ixO^S)*inv_gamma_1+w(ixO^S,eaux_c_)
    end where
  end subroutine twofl_energy_synchro

  subroutine twofl_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision :: rhoc(ixI^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision :: rhon(ixI^S)
#endif

    !TODO 
    return

    if(small_values_method == "ignore") return

    call twofl_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,rho_c_)) w(ixO^S,rho_c_) = small_density
#if !defined(ONE_FLUID) || ONE_FLUID==0
        where(flag(ixO^S,rho_n_)) w(ixO^S,rho_n_) = small_density
#endif
        do idir = 1, ndir
#if !defined(ONE_FLUID) || ONE_FLUID==0
          if(small_values_fix_iw(mom_n(idir))) then
            where(flag(ixO^S,rho_n_)) w(ixO^S, mom_n(idir)) = 0.0d0
          end if
#endif
          if(small_values_fix_iw(mom_c(idir))) then
            where(flag(ixO^S,rho_c_)) w(ixO^S, mom_c(idir)) = 0.0d0
          end if
        end do

        if(phys_energy) then
          if(primitive) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
            where(flag(ixO^S,e_n_)) w(ixO^S,e_n_) = small_pressure
#endif
            where(flag(ixO^S,e_c_)) w(ixO^S,e_c_) = small_pressure
          else if(phys_internal_e) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
            where(flag(ixO^S,e_n_))
              w(ixO^S,e_n_)=small_e
            end where
#endif
            where(flag(ixO^S,e_c_))
              w(ixO^S,e_c_)=small_e
            end where
          else
#if !defined(ONE_FLUID) || ONE_FLUID==0
            where(flag(ixO^S,e_n_))
              w(ixO^S,e_n_) = small_e+&
                 twofl_kin_en_n(w,ixI^L,ixO^L)
            end where
#endif
            if(phys_total_energy) then
              where(flag(ixO^S,e_c_))
                w(ixO^S,e_c_) = small_e+&
                   twofl_kin_en_c(w,ixI^L,ixO^L)+&
                   twofl_mag_en(w,ixI^L,ixO^L)
              end where
            else
              where(flag(ixO^S,e_c_))
                w(ixO^S,e_c_) = small_e+&
                   twofl_kin_en_c(w,ixI^L,ixO^L)
              end where
            endif
            if(phys_solve_eaux) then
              where(flag(ixO^S,e_c_))
                w(ixO^S,eaux_c_)=small_e
              end where
            end if
          end if
        end if
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag)
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(phys_energy) then
            if(phys_internal_e) then
              w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
#if !defined(ONE_FLUID) || ONE_FLUID==0
              w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
#endif
            else
#if !defined(ONE_FLUID) || ONE_FLUID==0
             w(ixO^S,e_n_)=gamma_1*(w(ixO^S,e_n_)&
                         -twofl_kin_en_n(w,ixI^L,ixO^L))
#endif
              if(phys_total_energy) then
                w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                          -twofl_kin_en_c(w,ixI^L,ixO^L)&
                          -twofl_mag_en(w,ixI^L,ixO^L))
               else
                 w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                          -twofl_kin_en_c(w,ixI^L,ixO^L))

               endif  
              if(phys_solve_eaux) w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)*gamma_1
            end if
          end if
          ! Convert momentum to velocity
#if !defined(ONE_FLUID) || ONE_FLUID==0
          if(has_equi_rho_n0) then
            rhon(ixG^T) = w(ixG^T,rho_n_) + block%equi_vars(ixG^T,equi_rho_n0_,0)
          else  
            rhon(ixG^T) = w(ixG^T,rho_n_) 
          endif
#endif
    
          if(has_equi_rho_c0) then
            rhoc(ixG^T) = w(ixG^T,rho_c_) + block%equi_vars(ixG^T,equi_rho_c0_,0)
          else  
            rhoc(ixG^T) = w(ixG^T,rho_c_) 
          endif
          do idir = 1, ndir
#if !defined(ONE_FLUID) || ONE_FLUID==0
             w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/rhon(ixO^S)
#endif
             w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/rhoc(ixO^S)
          end do
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine twofl_handle_small_values


  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine twofl_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision                :: vn(ixI^S),vc(ixI^S)

    call twofl_get_csound(w,x,ixI^L,ixO^L,idim,cmax)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    call twofl_get_v_n_idim(w,x,ixI^L,ixO^L,idim,vn)
#endif
    call twofl_get_v_c_idim(w,x,ixI^L,ixO^L,idim,vc)
    cmax(ixO^S)=&
#if !defined(ONE_FLUID) || ONE_FLUID==0
        max(abs(vn(ixO^S)),& 
#endif
            abs(vc(ixO^S))&
#if !defined(ONE_FLUID) || ONE_FLUID==0
            )&
#endif
        +cmax(ixO^S)

  end subroutine twofl_get_cmax

  subroutine twofl_get_a2max(w,x,ixI^L,ixO^L,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixI^S,ndim,nw)
    integer :: gxO^L,hxO^L,jxO^L,kxO^L,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxO^L=ixO^L-kr(i,^D);
      gxO^L=hxO^L-kr(i,^D);
      jxO^L=ixO^L+kr(i,^D);
      kxO^L=jxO^L+kr(i,^D);
      a2(ixO^S,i,1:nw)=abs(-w(kxO^S,1:nw)+16.d0*w(jxO^S,1:nw)&
         -30.d0*w(ixO^S,1:nw)+16.d0*w(hxO^S,1:nw)-w(gxO^S,1:nw))
      a2max(i)=maxval(a2(ixO^S,i,1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine twofl_get_a2max

#if !defined(ONE_FLUID) || ONE_FLUID==0
  ! COPIED from hd/moh_hd_phys
  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine twofl_get_tcutoff_n(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),w(ixI^S,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S)
    integer :: jxO^L,hxO^L
    logical :: lrlt(ixI^S)

    {^IFONED
    ! reuse lts as rhon
    call get_rhon_tot(w,ixI^L,ixO^L,lts)
    tmp1(ixI^S)=w(ixI^S,e_n_)-0.5d0*sum(w(ixI^S,mom_n(:))**2,dim=ndim+1)/lts(ixI^S)
    Te(ixI^S)=tmp1(ixI^S)/lts(ixI^S)*(twofl_gamma-1.d0)

    Tmax_local=maxval(Te(ixO^S))

    hxO^L=ixO^L-1;
    jxO^L=ixO^L+1;
    lts(ixO^S)=0.5d0*abs(Te(jxO^S)-Te(hxO^S))/Te(ixO^S)
    lrlt=.false.
    where(lts(ixO^S) > delta)
      lrlt(ixO^S)=.true.
    end where
    tco_local=zero
    if(any(lrlt(ixO^S))) then
      tco_local=maxval(Te(ixO^S), mask=lrlt(ixO^S))
    end if
    }
  end subroutine twofl_get_tcutoff_n
#endif

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine twofl_get_tcutoff_c(ixI^L,ixO^L,w,x,Tco_local,Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),w(ixI^S,1:nw)
    double precision, intent(out) :: Tco_local,Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: bunitvec
    double precision, dimension(ixI^S,1:ndim) :: gradT
    double precision :: Bdir(ndim)
    integer :: idims,jxO^L,hxO^L,ixA^D,ixB^D
    logical :: lrlt(ixI^S)

    ! reuse lts as rhoc
    call get_rhoc_tot(w,ixI^L,ixO^L,lts)
    if(phys_internal_e) then
      tmp1(ixI^S)=w(ixI^S,e_c_)
    else
      tmp1(ixI^S)=w(ixI^S,e_c_)-0.5d0*(sum(w(ixI^S,mom_c(:))**2,dim=ndim+1)/&
                       lts(ixI^S)+sum(w(ixI^S,mag(:))**2,dim=ndim+1))
    end if
    Te(ixI^S)=tmp1(ixI^S)/lts(ixI^S)*(twofl_gamma-1.d0)
    Tmax_local=maxval(Te(ixO^S))

    {^IFONED
    hxO^L=ixO^L-1;
    jxO^L=ixO^L+1;
    lts(ixO^S)=0.5d0*abs(Te(jxO^S)-Te(hxO^S))/Te(ixO^S)
    lrlt=.false.
    where(lts(ixO^S) > trac_delta)
      lrlt(ixO^S)=.true.
    end where
    Tco_local=zero
    block%special_values(1)=zero
    if(any(lrlt(ixO^S))) then
      Tco_local=zero
      block%special_values(1)=maxval(Te(ixO^S), mask=lrlt(ixO^S))
    end if
    }
    {^NOONED
    if(mod(twofl_trac_type,2) .eq. 1) then
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixI^L,ixO^L,idims,tmp1)
        gradT(ixO^S,idims)=tmp1(ixO^S)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixO^S,:)=w(ixO^S,iw_mag(:))+block%B0(ixO^S,:,0)
      else
        bunitvec(ixO^S,:)=w(ixO^S,iw_mag(:))
      end if
      if(twofl_trac_type>1) then
        ! B direction at cell center
        Bdir=zero
        {do ixA^D=0,1\}
          ixB^D=(ixOmin^D+ixOmax^D-1)/2+ixA^D;
          Bdir(1:ndim)=Bdir(1:ndim)+bunitvec(ixB^D,1:ndim)
        {end do\}
        if(sum(Bdir(:)**2) .gt. zero) then
          Bdir(1:ndim)=Bdir(1:ndim)/dsqrt(sum(Bdir(:)**2))
        end if
        block%special_values(3:ndim+2)=Bdir(1:ndim)
      end if
      tmp1(ixO^S)=dsqrt(sum(bunitvec(ixO^S,:)**2,dim=ndim+1))
      where(tmp1(ixO^S)/=0.d0)
        tmp1(ixO^S)=1.d0/tmp1(ixO^S)
      elsewhere
        tmp1(ixO^S)=bigdouble
      end where
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixO^S,idims)=bunitvec(ixO^S,idims)*tmp1(ixO^S)
      end do
      ! temperature length scale inversed
      lts(ixO^S)=abs(sum(gradT(ixO^S,1:ndim)*bunitvec(ixO^S,1:ndim),dim=ndim+1))/Te(ixO^S)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixO^S)=minval(dxlevel)*lts(ixO^S)
      else
        lts(ixO^S)=minval(block%ds(ixO^S,:),dim=ndim+1)*lts(ixO^S)
      end if
      lrlt=.false.
      where(lts(ixO^S) > trac_delta)
        lrlt(ixO^S)=.true.
      end where
      if(any(lrlt(ixO^S))) then
        block%special_values(1)=maxval(Te(ixO^S), mask=lrlt(ixO^S))
      else
        block%special_values(1)=zero
      end if
      block%special_values(2)=Tmax_local
    end if
    }
  end subroutine twofl_get_tcutoff_c

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine twofl_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision :: wmean(ixI^S,nw)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision :: rhon(ixI^S)
#endif
    double precision :: rhoc(ixI^S)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer                            :: idimE,idimN

    if (typeboundspeed=='cmaxmean') then
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
#if !defined(ONE_FLUID) || ONE_FLUID==0
      call get_rhon_tot(wmean,ixI^L,ixO^L,rhon)
      tmp2(ixO^S)=wmean(ixO^S,mom_n(idim))/rhon(ixO^S)
#endif
      call get_rhoc_tot(wmean,ixI^L,ixO^L,rhoc)
      tmp1(ixO^S)=wmean(ixO^S,mom_c(idim))/rhoc(ixO^S)
      call twofl_get_csound(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
        cmax(ixO^S)=max(max(abs(tmp2(ixO^S)), abs(tmp1(ixO^S)) ) +csoundR(ixO^S),zero)
        cmin(ixO^S)=min(min(abs(tmp2(ixO^S)),  abs(tmp1(ixO^S)) ) -csoundR(ixO^S),zero)
#else
        cmax(ixO^S)=max(abs(tmp1(ixO^S))+csoundR(ixO^S),zero)
        cmin(ixO^S)=min(abs(tmp1(ixO^S))-csoundR(ixO^S),zero)
#endif
      else
#if !defined(ONE_FLUID) || ONE_FLUID==0
        cmax(ixO^S)= max(abs(tmp2(ixO^S)),abs(tmp1(ixO^S)))+csoundR(ixO^S)
#else
        cmax(ixO^S)= abs(tmp1(ixO^S))+csoundR(ixO^S)

#endif

      end if
    else
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      call get_rhoc_tot(wLP,ixI^L,ixO^L,rhoc)
#if !defined(ONE_FLUID) || ONE_FLUID==0
      call get_rhon_tot(wLP,ixI^L,ixO^L,rhon)
      tmp1(ixO^S)=sqrt(abs(rhoc(ixO^S)  +rhon(ixO^S)))
#else
      tmp1(ixO^S)=sqrt(abs(rhoc(ixO^S)))
#endif

      call get_rhoc_tot(wRP,ixI^L,ixO^L,rhoc)
#if !defined(ONE_FLUID) || ONE_FLUID==0
      call get_rhon_tot(wRP,ixI^L,ixO^L,rhon)
      tmp2(ixO^S)=sqrt(abs(rhoc(ixO^S) +rhon(ixO^S)))
#else
      tmp2(ixO^S)=sqrt(abs(rhoc(ixO^S)))
#endif

      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
#if !defined(ONE_FLUID) || ONE_FLUID==0
      umean(ixO^S)=(0.5*(wLp(ixO^S,mom_n(idim))+wLp(ixO^S,mom_c(idim)))*tmp1(ixO^S) + &
                    0.5*(wRp(ixO^S,mom_n(idim))+wRp(ixO^S,mom_c(idim)))*tmp2(ixO^S))*tmp3(ixO^S)
#else
      umean(ixO^S)=(wLp(ixO^S,mom_c(idim))*tmp1(ixO^S)+wRp(ixO^S,mom_c(idim))*tmp2(ixO^S))*tmp3(ixO^S)
#endif
      call twofl_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call twofl_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)


#if !defined(ONE_FLUID) || ONE_FLUID==0
      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*(&
       0.5*(wRp(ixO^S,mom_n(idim))+wRp(ixO^S,mom_c(idim)))- & 
       0.5*(wLp(ixO^S,mom_n(idim))+wLp(ixO^S,mom_c(idim))))**2
#else
      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
       (wRp(ixO^S,mom_c(idim)) - wLp(ixO^S,mom_c(idim)))**2
#endif
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    end if

    !!TODO
!    if(stagger_grid) then
!      ! calculate velocities related to different UCT schemes
!      select case(type_ct)
!      case('average')
!      case('uct_contact')
!        if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixI^S,1:ndim))
!        ! get average normal velocity at cell faces
!        vcts%vnorm(ixO^S,idim)=0.5d0*(wLp(ixO^S,mom(idim))+wRp(ixO^S,mom(idim)))
!      case('uct_hll')
!        if(.not.allocated(vcts%vbarC)) then
!          allocate(vcts%vbarC(ixI^S,1:ndir,2),vcts%vbarLC(ixI^S,1:ndir,2),vcts%vbarRC(ixI^S,1:ndir,2))
!          allocate(vcts%cbarmin(ixI^S,1:ndim),vcts%cbarmax(ixI^S,1:ndim)) 
!        end if
!        ! Store magnitude of characteristics
!        if(present(cmin)) then
!          vcts%cbarmin(ixO^S,idim)=max(-cmin(ixO^S),zero)
!          vcts%cbarmax(ixO^S,idim)=max( cmax(ixO^S),zero)
!        else
!          vcts%cbarmax(ixO^S,idim)=max( cmax(ixO^S),zero)
!          vcts%cbarmin(ixO^S,idim)=vcts%cbarmax(ixO^S,idim)
!        end if
!
!        idimN=mod(idim,ndir)+1 ! 'Next' direction
!        idimE=mod(idim+1,ndir)+1 ! Electric field direction
!        ! Store velocities
!        vcts%vbarLC(ixO^S,idim,1)=wLp(ixO^S,mom(idimN))
!        vcts%vbarRC(ixO^S,idim,1)=wRp(ixO^S,mom(idimN))
!        vcts%vbarC(ixO^S,idim,1)=(vcts%cbarmax(ixO^S,idim)*vcts%vbarLC(ixO^S,idim,1) &
!             +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
!            /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))
!
!        vcts%vbarLC(ixO^S,idim,2)=wLp(ixO^S,mom(idimE))
!        vcts%vbarRC(ixO^S,idim,2)=wRp(ixO^S,mom(idimE))
!        vcts%vbarC(ixO^S,idim,2)=(vcts%cbarmax(ixO^S,idim)*vcts%vbarLC(ixO^S,idim,2) &
!             +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
!            /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))
!      case default
!        call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
!      end select
!    end if

  end subroutine twofl_get_cbounds

  !> Calculate fast magnetosonic wave speed
  subroutine twofl_get_csound(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S), gamma2(ixO^S)
    double precision :: rhoc(ixI^S)
!    double precision :: rhon(ixI^S)
    !TODO csound
    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
!#if !defined(ONE_FLUID) || ONE_FLUID==0
!    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
!    inv_rho(ixO^S) = 1d0/(rhon(ixO^S)+rhoc(ixO^S)) 
!#else
    inv_rho=1.d0/rhoc(ixO^S)
!#endif
    if (twofl_boris_type == boris_reduced_force) then
      call twofl_gamma2_alfven(ixI^L, ixO^L, w, gamma2)
    else
      gamma2 = 1.0d0
    end if

    call twofl_get_csound2(w,x,ixI^L,ixO^L,csound)

    ! store |B|^2 in v
    b2(ixO^S) = twofl_mag_en_all(w,ixI^L,ixO^L) * gamma2

    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * twofl_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho * gamma2

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. twofl_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
       if (twofl_boris_type == boris_simplification) then
          csound(ixO^S) = twofl_gamma_alfven(w, ixI^L,ixO^L) * csound(ixO^S)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            twofl_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
    end if

  end subroutine twofl_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine twofl_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S), gamma_A2(ixO^S)

    if(has_equi_rho_c0) then
      inv_rho=1.d0/(w(ixO^S,rho_c_)+block%equi_vars(ixO^S,equi_rho_c0_,0))
    else  
      inv_rho=1.d0/w(ixO^S,rho_c_)
    end if

    if (twofl_boris_type == boris_reduced_force) then
      call twofl_gamma2_alfven(ixI^L, ixO^L, w, gamma_A2)
    else
      gamma_A2 = 1.0d0
    end if

    if(phys_energy) then
#if !defined(ONE_FLUID) || ONE_FLUID==0
       call twofl_get_csound2_from_pe(w,x,ixI^L,ixO^L,w(ixI^S,e_c_),w(ixI^S,e_n_),csound)
#else
       call twofl_get_csound2_from_pe(w,x,ixI^L,ixO^L,w(ixI^S,e_c_),csound)
#endif
            
    else
       call twofl_get_csound2_adiab(w,x,ixI^L,ixO^L,csound)
    end if
    ! store |B|^2 in v
    b2(ixO^S)        = twofl_mag_en_all(w,ixI^L,ixO^L) * gamma_A2
    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * twofl_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho * gamma_A2

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. twofl_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
       if (twofl_boris_type == boris_simplification) then
          csound(ixO^S) = twofl_gamma_alfven(w, ixI^L,ixO^L) * csound(ixO^S)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            twofl_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
    end if

  end subroutine twofl_get_csound_prim

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine twofl_get_pthermal_c(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    if(phys_energy) then
      if(phys_internal_e) then
        pth(ixO^S)=gamma_1*w(ixO^S,e_c_)
      elseif(phys_total_energy) then
        pth(ixO^S)=gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)&
           - twofl_mag_en(w,ixI^L,ixO^L))
      else
        pth(ixO^S)=gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L))
      end if
    else
      pth(ixO^S)=twofl_adiab*w(ixO^S,rho_c_)**twofl_gamma
    end if

    ! no check when equi
    if(.not. has_equi_pe_c0) then
      if (check_small_values) then
        {do ix^DB= ixO^LIM^DB\}
           if(pth(ix^D)<small_pressure) then
             write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                  " encountered when call twofl_get_pthermal_c"
             write(*,*) "Iteration: ", it, " Time: ", global_time
             write(*,*) "Location: ", x(ix^D,:)
             write(*,*) "Cell number: ", ix^D
             do iw=1,nw
               write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
             end do
             ! use erroneous arithmetic operation to crash the run
             if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
             write(*,*) "Saving status at the previous time step"
             crash=.true.
           end if
        {enddo^D&\}
      end if
  
      if (fix_small_values) then
        {do ix^DB= ixO^LIM^DB\}
           if(pth(ix^D)<small_pressure) then
              pth(ix^D)=small_pressure
           end if
        {enddo^D&\}
      end if
    end if

  end subroutine twofl_get_pthermal_c

#if !defined(ONE_FLUID) || ONE_FLUID==0
  !> Calculate thermal pressure of neutrals =(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine twofl_get_pthermal_n(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    if(phys_energy) then
      if(phys_internal_e) then
        pth(ixO^S)=gamma_1*w(ixO^S,e_n_)
      else
        pth(ixO^S)=gamma_1*(w(ixO^S,e_n_)&
           - twofl_kin_en_n(w,ixI^L,ixO^L))
      end if
    else
      pth(ixO^S)=twofl_adiab*w(ixO^S,rho_n_)**twofl_gamma
    end if

    if(.not. has_equi_pe_n0) then
      if (check_small_values) then
        {do ix^DB= ixO^LIM^DB\}
           if(pth(ix^D)<small_pressure) then
             write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                  " encountered when call twofl_get_pthermal_n"
             write(*,*) "Iteration: ", it, " Time: ", global_time
             write(*,*) "Location: ", x(ix^D,:)
             write(*,*) "Cell number: ", ix^D
             do iw=1,nw
               write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
             end do
             ! use erroneous arithmetic operation to crash the run
             if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
             write(*,*) "Saving status at the previous time step"
             crash=.true.
           end if
        {enddo^D&\}
      end if
  
      if (fix_small_values) then
        {do ix^DB= ixO^LIM^DB\}
           if(pth(ix^D)<small_pressure) then
              pth(ix^D)=small_pressure
           end if
        {enddo^D&\}
      end if
    end if
  end subroutine twofl_get_pthermal_n
#endif

!  !> Calculate temperature=p/rho when in e_ the internal energy is stored
!  subroutine twofl_get_temperature_from_eint_n(w, x, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: w(ixI^S, 1:nw)
!    double precision, intent(in) :: x(ixI^S, 1:ndim)
!    double precision, intent(out):: res(ixI^S)
!    res(ixO^S) = gamma_1 * w(ixO^S, e_n_) /w(ixO^S,rho_n_)
!  end subroutine twofl_get_temperature_from_eint_n
!
!  !> Calculate temperature=p/rho when in e_ the internal energy is stored
!  subroutine twofl_get_temperature_from_eint_c(w, x, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: w(ixI^S, 1:nw)
!    double precision, intent(in) :: x(ixI^S, 1:ndim)
!    double precision, intent(out):: res(ixI^S)
!    res(ixO^S) = gamma_1 * w(ixO^S, e_c_) /w(ixO^S,rho_c_)
!  end subroutine twofl_get_temperature_from_eint_c
!
!  !> Calculate temperature=p/rho when in e_ the total energy is stored
!  !> this does not check the values of twofl_energy and twofl_internal_e, 
!  !>  twofl_energy = .true. and twofl_internal_e = .false.
!  !> also check small_values is avoided
!  subroutine twofl_get_temperature_from_etot_n(w, x, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: w(ixI^S, 1:nw)
!    double precision, intent(in) :: x(ixI^S, 1:ndim)
!    double precision, intent(out):: res(ixI^S)
!    res(ixO^S)=(gamma_1*(w(ixO^S,e_)&
!           - twofl_kin_en(w,ixI^L,ixO^L)&
!           - twofl_mag_en(w,ixI^L,ixO^L)))/w(ixO^S,rho_)
!  end subroutine twofl_get_temperature_from_etot_n
!  !> Calculate temperature=p/rho when in e_ the total energy is stored
!  !> this does not check the values of twofl_energy and twofl_internal_e, 
!  !>  twofl_energy = .true. and twofl_internal_e = .false.
!  !> also check small_values is avoided
!  subroutine twofl_get_temperature_from_etot_n(w, x, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: w(ixI^S, 1:nw)
!    double precision, intent(in) :: x(ixI^S, 1:ndim)
!    double precision, intent(out):: res(ixI^S)
!    res(ixO^S)=(gamma_1*(w(ixO^S,e_)&
!           - twofl_kin_en(w,ixI^L,ixO^L)&
!  end subroutine twofl_get_temperature_from_etot_n
!  !> Calculate temperature=p/rho when in e_ the total energy is stored
!  !> this does not check the values of twofl_energy and twofl_internal_e, 
!  !>  twofl_energy = .true. and twofl_internal_e = .false.
!  !> also check small_values is avoided
!  subroutine twofl_get_temperature_from_etot_c(w, x, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: w(ixI^S, 1:nw)
!    double precision, intent(in) :: x(ixI^S, 1:ndim)
!    double precision, intent(out):: res(ixI^S)
!    res(ixO^S)=(gamma_1*(w(ixO^S,e_)&
!           - twofl_kin_en(w,ixI^L,ixO^L)&
!           - twofl_mag_en(w,ixI^L,ixO^L)))/w(ixO^S,rho_)
!  end subroutine twofl_get_temperature_from_etot_c

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p_tot/rho_tot
  subroutine twofl_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: pe_c1(ixI^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision  :: pe_n1(ixI^S)
#endif

    if(phys_energy) then
      call twofl_get_pthermal_c(w,x,ixI^L,ixO^L,pe_c1)
#if !defined(ONE_FLUID) || ONE_FLUID==0
      call twofl_get_pthermal_n(w,x,ixI^L,ixO^L,pe_n1)
      call twofl_get_csound2_from_pe(w,x,ixI^L,ixO^L,pe_c1,pe_n1,csound2)
#else
      call twofl_get_csound2_from_pe(w,x,ixI^L,ixO^L,pe_c1,csound2)
#endif
    else
      call twofl_get_csound2_adiab(w,x,ixI^L,ixO^L,csound2)
    endif
  end subroutine twofl_get_csound2


  subroutine twofl_get_csound2_adiab(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: rhoc(ixI^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision  :: rhon(ixI^S)
#endif

    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
    csound2(ixO^S)=twofl_gamma*twofl_adiab*&
                  max((rhoc(ixO^S)**twofl_gamma + rhon(ixO^S)**twofl_gamma)/(rhoc(ixO^S)+ rhon(ixO^S)),&
                  rhon(ixO^S)**gamma_1,rhoc(ixO^S)**gamma_1) 
#else
    csound2(ixO^S)=twofl_gamma*twofl_adiab* rhoc(ixO^S)**gamma_1
#endif
  end subroutine twofl_get_csound2_adiab



#if !defined(ONE_FLUID) || ONE_FLUID==0
  subroutine twofl_get_csound2_from_pe(w,x,ixI^L,ixO^L,pe_c1,pe_n1,csound2)
                                      &
                                      
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: pe_c1(ixI^S)
    double precision, intent(in)    :: pe_n1(ixI^S)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: csound1(ixI^S),rhon(ixI^S),rhoc(ixI^S)

    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
    if(has_equi_pe_n0) then
      csound1(ixO^S) = pe_n1(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,block%iw0)
    else
      csound1(ixO^S) = pe_n1(ixO^S) 
    endif

    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
    if(has_equi_pe_c0) then
      csound2(ixO^S) = pe_c1(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,block%iw0)
    else
      csound2(ixO^S) = pe_c1(ixO^S) 
    endif
     !TODO csound 
!    csound2(ixO^S)=twofl_gamma*(csound2(ixO^S) + csound1(ixO^S))/(rhoc(ixO^S) + rhon(ixO^S))
    csound2(ixO^S)=twofl_gamma*max((csound2(ixO^S) + csound1(ixO^S))/(rhoc(ixO^S) + rhon(ixO^S)),&
                      csound1(ixO^S)/rhon(ixO^S), csound2(ixO^S)/rhoc(ixO^S))
  end subroutine twofl_get_csound2_from_pe

#else
  subroutine twofl_get_csound2_from_pe(w,x,ixI^L,ixO^L,pe_c1,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: pe_c1(ixI^S)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision                :: rhoc(ixI^S)

    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
    if(has_equi_pe_c0) then
      csound2(ixO^S) = pe_c1(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,block%iw0)
    else
      csound2(ixO^S) = pe_c1(ixO^S) 
    endif
    csound2(ixO^S)=twofl_gamma* &
                      csound2(ixO^S)/rhoc(ixO^S)
  end subroutine twofl_get_csound2_from_pe
#endif


  !> Calculate total pressure within ixO^L including magnetic pressure
!  subroutine twofl_get_p_total(w,x,ixI^L,ixO^L,p)
!    use mod_global_parameters
!
!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision, intent(in)    :: w(ixI^S,nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(out)   :: p(ixI^S)
!
!    call twofl_get_pthermal(w,x,ixI^L,ixO^L,p)
!
!    p(ixO^S) = p(ixO^S) + 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
!
!  end subroutine twofl_get_p_total

  !> Calculate fluxes within ixO^L.
  subroutine twofl_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: pgas(ixO^S), ptotal(ixO^S),tmp(ixI^S)
    double precision, allocatable:: vHall(:^D&,:)
    integer                      :: idirmin, iw, idir, jdir, kdir
!ambipolar
#if defined(ONE_FLUID) && ONE_FLUID==1
    double precision, allocatable, dimension(:^D&,:) :: Jambi, btot
    double precision, allocatable, dimension(:^D&) :: tmp2, tmp3, tmp4
#endif

    ! value at the interfaces, idim =  block%iw0 
    ! reuse tmp, used afterwards
    ! value at the interface so we can't put momentum
    call get_rhoc_tot(w,ixI^L,ixO^L,tmp)
    ! Get flux of density
    f(ixO^S,rho_c_)=w(ixO^S,mom_c(idim))*tmp(ixO^S)

    if (twofl_Hall) then
      allocate(vHall(ixI^S,1:ndir))
      call twofl_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    end if

    if(B0field) tmp(ixO^S)=sum(block%B0(ixO^S,:,idim)*w(ixO^S,mag(:)),dim=ndim+1)

    if(phys_energy) then
      pgas(ixO^S)=w(ixO^S,e_c_)
    else
      pgas(ixO^S)=twofl_adiab*w(ixO^S,rho_c_)**twofl_gamma
    end if

    ptotal(ixO^S) = pgas(ixO^S) + 0.5d0*sum(w(ixO^S, mag(:))**2, dim=ndim+1)



    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    if (twofl_boris_type == boris_reduced_force) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixO^S,mom_c(idir)) = pgas(ixO^S)
        else
          f(ixO^S,mom_c(idir)) = 0.0d0
        end if
        f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+w(ixO^S,mom_c(idim))*wC(ixO^S,mom_c(idir))
      end do
    else
      ! Normal case (no Boris approximation)
      do idir=1,ndir
        if(idim==idir) then
          f(ixO^S,mom_c(idir))=ptotal(ixO^S)-w(ixO^S,mag(idim))*w(ixO^S,mag(idir))
          if(B0field) f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+tmp(ixO^S)
        else
          f(ixO^S,mom_c(idir))= -w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
        end if
        if (B0field) then
          f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))&
               -w(ixO^S,mag(idir))*block%B0(ixO^S,idim,idim)&
               -w(ixO^S,mag(idim))*block%B0(ixO^S,idir,idim)
        end if
        f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+w(ixO^S,mom_c(idim))*wC(ixO^S,mom_c(idir))
      end do
    end if



    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(phys_energy) then
      if (phys_internal_e) then
         f(ixO^S,e_c_)=w(ixO^S,mom_c(idim))*w(ixO^S,e_c_)
         if (twofl_Hall) then
            call mpistop("solve internal energy not implemented for Hall MHD")
         endif
      else if(twofl_eq_energy == EQ_ENERGY_KI) then

        f(ixO^S,e_c_)=w(ixO^S,mom_c(idim))*(wC(ixO^S,e_c_)+pgas(ixO^S))
      else
        f(ixO^S,e_c_)=w(ixO^S,mom_c(idim))*(wC(ixO^S,e_c_)+ptotal(ixO^S))&
           -w(ixO^S,mag(idim))*sum(w(ixO^S,mag(:))*w(ixO^S,mom_c(:)),dim=ndim+1)
        !if(phys_solve_eaux) f(ixO^S,eaux_)=w(ixO^S,mom(idim))*wC(ixO^S,eaux_)

        if (B0field) then
           f(ixO^S,e_c_) = f(ixO^S,e_c_) &
              + w(ixO^S,mom_c(idim)) * tmp(ixO^S) &
              - sum(w(ixO^S,mom_c(:))*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
        end if

        if (twofl_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (twofl_etah>zero) then
              f(ixO^S,e_c_) = f(ixO^S,e_c_) + vHall(ixO^S,idim) * &
                 sum(w(ixO^S, mag(:))**2,dim=ndim+1) &
                 - w(ixO^S,mag(idim)) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
              if (B0field) then
                 f(ixO^S,e_c_) = f(ixO^S,e_c_) &
                    + vHall(ixO^S,idim) * tmp(ixO^S) &
                    - sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
              end if
           end if
        end if
      end if !total_energy
      ! add flux of equilibrium internal energy corresponding to pe_c0
      if(has_equi_pe_c0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
        f(ixO^S,e_c_)=  f(ixO^S,e_c_) &
          + w(ixO^S,mom_c(idim)) * block%equi_vars(ixO^S,equi_pe_c0_,idim) * inv_gamma_1
#else
        f(ixO^S,e_c_)=  f(ixO^S,e_c_) &
          + w(ixO^S,mom_c(idim)) * block%equi_vars(ixO^S,equi_pe_c0_,idim) * twofl_gamma * inv_gamma_1
#endif
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (twofl_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=w(ixO^S,mom_c(idim))*w(ixO^S,mag(idir))-w(ixO^S,mag(idim))*w(ixO^S,mom_c(idir))

        if (B0field) then
          f(ixO^S,mag(idir))=f(ixO^S,mag(idir))&
                +w(ixO^S,mom_c(idim))*block%B0(ixO^S,idir,idim)&
                -w(ixO^S,mom_c(idir))*block%B0(ixO^S,idim,idim)
        end if

        if (twofl_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (twofl_etah>zero) then
            if (B0field) then
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*(w(ixO^S,mag(idim))+block%B0(ixO^S,idim,idim)) &
                   + vHall(ixO^S,idim)*(w(ixO^S,mag(idir))+block%B0(ixO^S,idir,idim))
            else
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*w(ixO^S,mag(idim)) &
                   + vHall(ixO^S,idim)*w(ixO^S,mag(idir))
            end if
          end if
        end if

      end if
    end do

    if (twofl_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

    if (twofl_Hall) then
      deallocate(vHall)
    end if

#if !defined(ONE_FLUID) || ONE_FLUID==0
    !!neutrals
    call get_rhon_tot(w,ixI^L,ixO^L,tmp)
    f(ixO^S,rho_n_)=w(ixO^S,mom_n(idim))*tmp(ixO^S)
    if(phys_energy) then
      pgas(ixO^S) = w(ixO^S, e_n_)
    else
      pgas(ixO^S)=twofl_adiab*w(ixO^S,rho_n_)**twofl_gamma
    endif
    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
       f(ixO^S, mom_n(idir)) = w(ixO^S,mom_n(idim)) * wC(ixO^S, mom_n(idir))
    end do

    f(ixO^S, mom_n(idim)) = f(ixO^S, mom_n(idim)) + pgas(ixO^S)

    if(phys_energy) then
      !reuse pgas for storing a in the term: div (u_n * a) and make multiplication at the end
      pgas(ixO^S) = wC(ixO^S,e_n_)
      if(.not. phys_internal_e) then
        ! add pressure
        pgas(ixO^S) = pgas(ixO^S) + w(ixO^S,e_n_)
      endif
      ! add flux of equilibrium internal energy corresponding to pe_n0
      if(has_equi_pe_n0) then
        pgas(ixO^S) = pgas(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,idim) * inv_gamma_1
      endif
      ! add u_n * a in the flux
      f(ixO^S, e_n_) = w(ixO^S,mom_n(idim)) * pgas(ixO^S)

      ! Viscosity fluxes - viscInDiv
      !if (hd_viscosity) then
      !  call visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, phys_energy)
      !endif
    end if
#else
    ! Contributions of ambipolar term in explicit scheme
    if(twofl_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixI^S,1:3))
      call twofl_get_Jambi(w,x,ixI^L,ixO^L,Jambi)
      allocate(btot(ixO^S,1:3))
      if(B0field) then
        do idir=1,3
          btot(ixO^S, idir) = w(ixO^S,mag(idir)) + block%B0(ixO^S,idir,idim)
        enddo
      else
        btot(ixO^S,1:3) = w(ixO^S,mag(1:3))
      endif
      allocate(tmp2(ixO^S),tmp3(ixO^S))
      !tmp2 = Btot^2
      tmp2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixO^S) = sum(Jambi(ixO^S,:)*btot(ixO^S,:),dim=ndim+1)




      !v->Jambi
      !Exb
      !{
      !bx*(by1*bz - by*bz1)*vx - (by*by1 + bz*bz1)*(-bz*vy + by*vz) + bx^2*(bz1*vy - by1*vz), 
      !-bz^2*bz1*vx - bx1*by*bz*vy + bx^2*bx1*vz + by^2*(-bz1*vx + bx1*vz) + bx*(-bx1*bz*vx + by*bz1*vy + bz*bz1*vz), 
      !bx*bx1*by*vx + by^2*by1*vx - bx^2*bx1*vy + bz^2*(by1*vx - bx1*vy) + bx1*by*bz*vz - bx*by1*(by*vy + bz*vz)
      !}

      !b2 = bx^2 + by^2 + bz^2
      !t1 = bz1 vy - by1 vz
      !t2 = bx1 vz - bz1 vx
      !t3 = by1 vx - bx1 vy
      !
      !Print["Energy x remainder"]
      !Print[Simplify[exb[[1]] - b2 * t1  ]]
      !Print["Energy y remainder"]
      !Print[Simplify[exb[[2]] - b2 * t2  ]]
      !Print["Energy z remainder"]
      !Print[Simplify[exb[[3]] - b2 * t3  ]]
      !Energy x remainder
      !(by1*bz - by*bz1)*(bx*vx + by*vy + bz*vz)
      !Energy y remainder
      !(-(bx1*bz) + bx*bz1)*(bx*vx + by*vy + bz*vz)
      !Energy z remainder
      !(bx1*by - bx*by1)*(bx*vx + by*vy + bz*vz)


      !mag1
      !{0, -(bz*(bx*vx + by*vy)) + (bx^2 + by^2)*vz, bx*by*vx - bx^2*vy + bz*(-(bz*vy) + by*vz)}
      !mag2
      !{bz*(bx*vx + by*vy) - (bx^2 + by^2)*vz, 0, by^2*vx - bx*by*vy - bz*(-(bz*vx) + bx*vz)}
      !mag3
      !{-(bx*by*vx) + bx^2*vy - bz*(-(bz*vy) + by*vz), -(by^2*vx) + bx*by*vy + bz*(-(bz*vx) + bx*vz), 0}


      if (B0field) allocate(tmp4(ixO^S))

      select case(idim)
        case(1)
          tmp(ixO^S)=w(ixO^S,mag(3)) *Jambi(ixO^S,2) - w(ixO^S,mag(2)) * Jambi(ixO^S,3)
          if(B0field) tmp4(ixO^S) = w(ixO^S,mag(2)) * btot(ixO^S,3) - w(ixO^S,mag(3)) * btot(ixO^S,2)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) - tmp2(ixO^S) * Jambi(ixO^S,3) + tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) + tmp2(ixO^S) * Jambi(ixO^S,2) - tmp3(ixO^S) * btot(ixO^S,2)
        case(2)
          tmp(ixO^S)=w(ixO^S,mag(1)) *Jambi(ixO^S,3) - w(ixO^S,mag(3)) * Jambi(ixO^S,1)
          if(B0field) tmp4(ixO^S) = w(ixO^S,mag(3)) * btot(ixO^S,1) - w(ixO^S,mag(1)) * btot(ixO^S,3)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) + tmp2(ixO^S) * Jambi(ixO^S,3) - tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) - tmp2(ixO^S) * Jambi(ixO^S,1) + tmp3(ixO^S) * btot(ixO^S,1)
        case(3)
          tmp(ixO^S)=w(ixO^S,mag(2)) *Jambi(ixO^S,1) - w(ixO^S,mag(1)) * Jambi(ixO^S,2)
          if(B0field) tmp4(ixO^S) = w(ixO^S,mag(1)) * btot(ixO^S,2) - w(ixO^S,mag(2)) * btot(ixO^S,1)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) - tmp2(ixO^S) * Jambi(ixO^S,2) + tmp3(ixO^S) * btot(ixO^S,2)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) + tmp2(ixO^S) * Jambi(ixO^S,1) - tmp3(ixO^S) * btot(ixO^S,1)
      endselect

      if(phys_total_energy) then
        f(ixO^S,e_c_) = f(ixO^S,e_c_) + tmp2(ixO^S) *  tmp(ixO^S)
        if(B0field) f(ixO^S,e_c_) = f(ixO^S,e_c_) +  tmp3(ixO^S) *  tmp4(ixO^S)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
      if (B0field) deallocate(tmp4)
    endif

#endif
!    print*, "2GETFLUX mom_c", f(ixOmin1:ixOmin1+5,mom_c(1))
!    print*, "2GETFLUX mom_c2", f(ixOmin1:ixOmin1+5,mom_c(2))
!    print*, "GETFLUX rho_c", f(ixOmin1:ixOmin1+5,rho_c_)
!    print*, "GETFLUX e_c", f(ixOmin1:ixOmin1+5,e_c_)
!    print*, "GETFLUX b", f(ixOmin1:ixOmin1+5,mag(1:3))
!    print*, "2GETFLUX mom_n", f(ixOmin1:ixOmin1+5,mom_n(1))
!    print*, "2GETFLUX mom_n2", f(ixOmin1:ixOmin1+5,mom_n(2))
!    print*, "GETFLUX rho_n", f(ixOmin1:ixOmin1+5,rho_n_)

  end subroutine twofl_get_flux


#if defined(ONE_FLUID) && ONE_FLUID==1
  !> Source terms J.E in internal energy. 
  !> For the ambipolar term E = ambiCoef * JxBxB=ambiCoef * B^2(-J_perpB) 
  !=> the source term J.E = ambiCoef * B^2 * J_perpB^2 = ambiCoef * JxBxB^2/B^2
  !> ambiCoef is calculated as mhd_ambi_coef/rho^2,  see also the subroutine mhd_get_Jambi
  subroutine add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x,ie)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: tmp(ixI^S)
    double precision :: jxbxb(ixI^S,1:3)

    call twofl_get_jxbxb(wCT,x,ixI^L,ixO^L,jxbxb)
    tmp(ixO^S) = sum(jxbxb(ixO^S,1:3)**2,dim=ndim+1) / twofl_mag_en_all(wCT, ixI^L, ixO^L)
    call multiplyAmbiCoef(ixI^L,ixO^L,tmp,wCT,x)   
    w(ixO^S,ie)=w(ixO^S,ie)+qdt * tmp

  end subroutine add_source_ambipolar_internal_energy

  subroutine twofl_get_jxbxb(w,x,ixI^L,ixO^L,res)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: res(:^D&,:)

    double precision  :: btot(ixI^S,1:3)
    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: tmp(ixI^S),b2(ixI^S)

    res=0.d0
    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    !!!here we know that current has nonzero values only for components in the range idirmin, 3
 
    if(B0field) then
      do idir=1,3
        btot(ixO^S, idir) = w(ixO^S,mag(idir)) + block%B0(ixO^S,idir,block%iw0)
      enddo
    else
      btot(ixO^S,1:3) = w(ixO^S,mag(1:3))
    endif
    tmp(ixO^S) = sum(current(ixO^S,idirmin:3)*btot(ixO^S,idirmin:3),dim=ndim+1) !J.B
    b2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1) !B^2
    do idir=1,idirmin-1
      res(ixO^S,idir) = btot(ixO^S,idir) * tmp(ixO^S)
    enddo
    do idir=idirmin,3
      res(ixO^S,idir) = btot(ixO^S,idir) * tmp(ixO^S) - current(ixO^S,idir) * b2(ixO^S)
    enddo
  end subroutine twofl_get_jxbxb

  !> Sets the sources for the ambipolar
  !> this is used for the STS method
  ! The sources are added directly (instead of fluxes as in the explicit) 
  !> at the corresponding indices
  !>  store_flux_var is explicitly called for each of the fluxes one by one
  subroutine sts_set_source_ambipolar(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixI^L, ixO^L,igrid,nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step

    double precision, dimension(ixI^S,1:3) :: tmp,ff
    double precision, allocatable, dimension(:^D&,:,:) :: fluxall
    double precision, allocatable :: fE(:^D&,:)
    double precision  :: btot(ixI^S,1:3),tmp2(ixI^S)
    integer :: i, ixA^L, ie_

    ixA^L=ixO^L^LADD1;

    call twofl_get_jxbxb(w,x,ixI^L,ixA^L,tmp)

    ! set electric field in tmp: E=nuA * jxbxb, where nuA=-etaA/rho^2
    do i=1,3
      call multiplyAmbiCoef(ixI^L,ixA^L,tmp(ixI^S,i),w,x)   
    enddo

    if(fix_conserve_at_step) then
      allocate(fluxall(ixI^S,1:nflux,1:ndim))
      fluxall=0.d0
    end if

    if(phys_total_energy ) then
      ! tmp is E
      ! btot is B1
      ! energy term is  -div(ExB1) 
      btot(ixA^S,1:3)=0.d0
      !TODO this has to be mag pert only! CHECK!!!!
      !if(B0field) then
      !  do i=1,ndir
      !    btot(ixA^S, i) = w(ixA^S,mag(i)) + block%B0(ixA^S,i,0)
      !  enddo
      !else
        btot(ixA^S,1:ndir) = w(ixA^S,mag(1:ndir))
      !endif
      call cross_product(ixI^L,ixA^L,tmp,btot,ff)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,1,1:ndim)=ff(ixI^S,1:ndim)
      !- sign comes from the fact that the flux divergence is a source now
      wres(ixO^S,e_c_)=-tmp2(ixO^S)
    endif

    if(stagger_grid) then
      if(ndir>ndim) then
        !!!Bz
        ff(ixA^S,1) = tmp(ixA^S,2)
        ff(ixA^S,2) = -tmp(ixA^S,1)
        ff(ixA^S,3) = 0.d0
        call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixI^S,1+ndir,1:ndim)=ff(ixI^S,1:ndim)
        wres(ixO^S,mag(ndir))=-tmp2(ixO^S)
      end if
      allocate(fE(ixI^S,7-2*ndim:3))
      call update_faces_ambipolar(ixI^L,ixO^L,w,x,tmp,fE,btot)
      ixAmax^D=ixOmax^D;
      ixAmin^D=ixOmin^D-1;
      wres(ixA^S,mag(1:ndim))=-btot(ixA^S,1:ndim)
    else
      !write curl(ele) as the divergence
      !m1={0,ele[[3]],-ele[[2]]}
      !m2={-ele[[3]],0,ele[[1]]}
      !m3={ele[[2]],-ele[[1]],0}

      !!!Bx
      ff(ixA^S,1) = 0.d0
      ff(ixA^S,2) = tmp(ixA^S,3)
      ff(ixA^S,3) = -tmp(ixA^S,2)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,2,1:ndim)=ff(ixI^S,1:ndim)
      !flux divergence is a source now
      wres(ixO^S,mag(1))=-tmp2(ixO^S)
      !!!By
      ff(ixA^S,1) = -tmp(ixA^S,3)
      ff(ixA^S,2) = 0.d0
      ff(ixA^S,3) = tmp(ixA^S,1)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,3,1:ndim)=ff(ixI^S,1:ndim)
      wres(ixO^S,mag(2))=-tmp2(ixO^S)

      if(ndir==3) then
        !!!Bz
        ff(ixA^S,1) = tmp(ixA^S,2)
        ff(ixA^S,2) = -tmp(ixA^S,1)
        ff(ixA^S,3) = 0.d0
        call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixI^S,1+ndir,1:ndim)=ff(ixI^S,1:ndim)
        wres(ixO^S,mag(ndir))=-tmp2(ixO^S)
      end if

    end if

    if(fix_conserve_at_step) then
      fluxall=my_dt*fluxall
      call store_flux(igrid,fluxall,1,ndim,nflux)
      if(stagger_grid) then
        call store_edge(igrid,ixI^L,my_dt*fE,1,ndim)
        deallocate(fE)
      end if
      deallocate(fluxall)
    end if

  end subroutine sts_set_source_ambipolar

  !> get ambipolar electric field and the integrals around cell faces
  subroutine update_faces_ambipolar(ixI^L,ixO^L,w,x,ECC,fE,circ)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    ! amibipolar electric field at cell centers
    double precision, intent(in)       :: ECC(ixI^S,1:3)
    double precision, intent(out)      :: fE(ixI^S,7-2*ndim:3)
    double precision, intent(out)      :: circ(ixI^S,1:ndim)

    integer                            :: hxC^L,ixC^L,ixA^L
    integer                            :: idim1,idim2,idir,ix^D

    fE=zero
    ! calcuate ambipolar electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D+kr(idir,^D)-1;
     {do ix^DB=0,1\}
        if({ ix^D==1 .and. ^D==idir | .or.}) cycle
        ixAmin^D=ixCmin^D+ix^D;
        ixAmax^D=ixCmax^D+ix^D;
        fE(ixC^S,idir)=fE(ixC^S,idir)+ECC(ixA^S,idir)
     {end do\}
      fE(ixC^S,idir)=fE(ixC^S,idir)*0.25d0*block%dsC(ixC^S,idir)
    end do

    ! Calculate circulation on each face to get value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    circ=zero

    do idim1=1,ndim ! Coordinate perpendicular to face 
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
      circ(ixC^S,idim1)=circ(ixC^S,idim1)/block%surfaceC(ixC^S,idim1)
    end do

  end subroutine update_faces_ambipolar

  !> use cell-center flux to get cell-face flux
  !> and get the source term as the divergence of the flux
  subroutine get_flux_on_cell_face(ixI^L,ixO^L,ff,src)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, dimension(:^D&,:), intent(inout) :: ff
    double precision, intent(out) :: src(ixI^S)

    double precision :: ffc(ixI^S,1:ndim)
    double precision :: dxinv(ndim)
    integer :: idims, ix^D, ixA^L, ixB^L, ixC^L

    ixA^L=ixO^L^LADD1;
    dxinv=1.d0/dxlevel
    ! cell corner flux in ffc
    ffc=0.d0
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
    {do ix^DB=0,1\}
      ixBmin^D=ixCmin^D+ix^D;
      ixBmax^D=ixCmax^D+ix^D;
      ffc(ixC^S,1:ndim)=ffc(ixC^S,1:ndim)+ff(ixB^S,1:ndim)
    {end do\}
    ffc(ixC^S,1:ndim)=0.5d0**ndim*ffc(ixC^S,1:ndim)
    ! flux at cell face
    ff(ixI^S,1:ndim)=0.d0
    do idims=1,ndim
      ixB^L=ixO^L-kr(idims,^D);
      ixCmax^D=ixOmax^D; ixCmin^D=ixBmin^D;
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixCmin^D-ix^D;
           ixBmax^D=ixCmax^D-ix^D;
           ff(ixC^S,idims)=ff(ixC^S,idims)+ffc(ixB^S,idims)
         end if
      {end do\}
      ff(ixC^S,idims)=ff(ixC^S,idims)*0.5d0**(ndim-1)
    end do
    src=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        ff(ixA^S,idims)=dxinv(idims)*ff(ixA^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        src(ixO^S)=src(ixO^S)+ff(ixO^S,idims)-ff(ixB^S,idims)
      end do
    else
      do idims=1,ndim
        ff(ixA^S,idims)=ff(ixA^S,idims)*block%surfaceC(ixA^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        src(ixO^S)=src(ixO^S)+ff(ixO^S,idims)-ff(ixB^S,idims)
      end do
      src(ixO^S)=src(ixO^S)/block%dvolume(ixO^S)
    end if
  end subroutine get_flux_on_cell_face
  !> Calculates the explicit dt for the ambipokar term
  !> This function is used by both explicit scheme and STS method
  function get_ambipolar_dt(w,ixI^L,ixO^L,dx^D,x)  result(dtnew)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    double precision              :: coef
    double precision              :: dxarr(ndim)
    double precision              :: tmp(ixI^S)
    ^D&dxarr(^D)=dx^D;

    tmp(ixO^S) = twofl_mag_en_all(w, ixI^L, ixO^L)
    call multiplyAmbiCoef(ixI^L,ixO^L,tmp,w,x) 
    coef = maxval(abs(tmp(ixO^S)))
    if(slab_uniform) then
      dtnew=minval(dxarr(1:ndim))**2.0d0/coef
    else
      dtnew=minval(block%ds(ixO^S,1:ndim))**2.0d0/coef
    end if

  end function get_ambipolar_dt

  !> multiply res by the ambipolar coefficient
  !> The ambipolar coefficient is calculated as -mhd_eta_ambi/rho^2
  !> The user may mask its value in the user file
  !> by implemneting usr_mask_ambipolar subroutine
  subroutine multiplyAmbiCoef(ixI^L,ixO^L,res,w,x)   
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: res(ixI^S)
    double precision :: tmp(ixI^S)
    double precision :: rhoc(ixI^S)

    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
    ! pu density (split or not) in tmp
    !print* , "MULTAMB TOTRHO ", tmp(ixOmin1:ixOmin1+10)
    tmp(ixO^S) = -(twofl_eta_ambi/rhoc(ixO^S)**2) 
    !print* , "MULTAMB NUA ", tmp(ixOmin1:ixOmin1+10)
    if (associated(usr_mask_ambipolar)) then
      call usr_mask_ambipolar(ixI^L,ixO^L,w,x,tmp)
    endif

    res(ixO^S) = tmp(ixO^S) * res(ixO^S)
  end subroutine multiplyAmbiCoef
#endif

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine twofl_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    !use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if(phys_internal_e) then
        active = .true.
#if !defined(ONE_FLUID) || ONE_FLUID==0
        call internal_energy_add_source_n(qdt,ixI^L,ixO^L,wCT,w,x)
#endif
        call internal_energy_add_source_c(qdt,ixI^L,ixO^L,wCT,w,x,e_c_)
        if(phys_solve_eaux) then
          call internal_energy_add_source_c(qdt,ixI^L,ixO^L,wCT,w,x,eaux_c_)
        endif
      else 
#if !defined(E_RM_W0) || E_RM_W0==1
        ! add -p0 div v source terms when equi are present
#if !defined(ONE_FLUID) || ONE_FLUID==0
        if(has_equi_pe_n0) then
          call add_pe_n0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
        endif
#endif
        if(has_equi_pe_c0) then
          call add_pe_c0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
        endif
#endif
        if(twofl_eq_energy == EQ_ENERGY_KI) then
          active = .true.
          call add_source_lorentz_work(qdt,ixI^L,ixO^L,w,wCT,x)
        endif
      endif
#if defined(ONE_FLUID) && ONE_FLUID==1
      if((twofl_eq_energy == EQ_ENERGY_KI .or. twofl_eq_energy == EQ_ENERGY_INT)&
           .and. twofl_ambipolar) then
        call add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x,e_c_)
      endif
#endif
      

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(twofl_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if (twofl_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
      end if
#if !defined(ONE_FLUID) || ONE_FLUID==0
      !it is not added in a split manner
      if(.not. twofl_implicit_coll_terms.and. twofl_alpha_coll>0d0) then
        active = .true.
        call  twofl_explicit_coll_terms_update(qdt,ixI^L,ixO^L,w,wCT,x)
      endif
#endif
    end if

      {^NOONED
    if(.not.source_split_divb .and. .not.qsourcesplit .and. istep==nstep) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
        call add_source_janhunen(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
        call add_source_powel(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
        call add_source_glm(dt,ixI^L,ixO^L,pso(saveigrid)%w,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    else if(source_split_divb .and. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
    }

!    if(twofl_radiative_cooling) then
!      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,qsourcesplit,active)
!    end if
!
!    if(twofl_viscosity) then
!      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,phys_energy,qsourcesplit,active)
!    end if
!
    if(twofl_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,twofl_eq_energy .eq. EQ_ENERGY_KI .or. phys_total_energy,qsourcesplit,active)
    end if

    if (twofl_boris_type == boris_reduced_force) then
      call boris_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,qsourcesplit,active)
    end if


  end subroutine twofl_add_source

#if !defined(ONE_FLUID) || ONE_FLUID==0
  subroutine add_pe_n0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v(ixI^S,1:ndir)

    call twofl_get_v_n(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-block%equi_vars(ixI^S,equi_pe_n0_,0),w,x,e_n_)

  end subroutine add_pe_n0_divv
#endif

  subroutine add_pe_c0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v(ixI^S,1:ndir)

    call twofl_get_v_c(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-block%equi_vars(ixI^S,equi_pe_c0_,0),w,x,e_c_)

  end subroutine add_pe_c0_divv


  subroutine add_geom_PdivV(qdt,ixI^L,ixO^L,v,p,w,x,ind)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L,ind
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: p(ixI^S), v(ixI^S,1:ndir), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divv(ixI^S)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixI^L,ixO^L,divv,sixthorder=.true.)
      else
        call divvector(v,ixI^L,ixO^L,divv,fourthorder=.true.)
      end if
    else
     call divvector(v,ixI^L,ixO^L,divv)
    end if
    w(ixO^S,ind)=w(ixO^S,ind)+qdt*p(ixO^S)*divv(ixO^S)
  end subroutine add_geom_PdivV



  subroutine boris_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
       qsourcesplit,active)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active

    double precision :: JxB(ixI^S,3)
    double precision :: gamma_A2(ixO^S)
    integer          :: idir

    ! Boris source term is always unsplit
    if (qsourcesplit) return

    call get_lorentz(ixI^L,ixO^L,w,JxB)
    call twofl_gamma2_alfven(ixI^L, ixO^L, wCT, gamma_A2)

    do idir = 1, ndir
      w(ixO^S,mom_c(idir)) = w(ixO^S,mom_c(idir)) + &
           qdt * gamma_A2 * JxB(ixO^S, idir)
    end do

  end subroutine boris_add_source

  !> Compute the Lorentz force (JxB)
  subroutine get_lorentz(ixI^L,ixO^L,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: JxB(ixI^S,3)
    double precision                :: a(ixI^S,3), b(ixI^S,3), tmp(ixI^S,3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixO^S, idir) = twofl_mag_i_all(w, ixI^L, ixO^L,idir)
    end do

    ! store J current in a
    call get_current(w,ixI^L,ixO^L,idirmin,current)

    a=0.0d0
    do idir=7-2*ndir,3
      a(ixO^S,idir)=current(ixO^S,idir)
    end do

    call cross_product(ixI^L,ixO^L,a,b,JxB)

    ! for B0  splitting add contribution of J0

    if(B0field) then
      a=0.0d0
      do idir=7-2*ndir,3
        a(ixO^S,idir)=block%J0(ixO^S,idir)
      end do
      if(B0field_forcefree) then
        ! store in b only b1
        b=0.0d0
        do idir = 1, ndir
          b(ixO^S, idir) = w(ixO^S,mag(idir))
        end do
      endif
      call cross_product(ixI^L,ixO^L,a,b,tmp)
      do idir = 1, ndir
        JxB(ixO^S, idir) = JxB(ixO^S,idir) + tmp(ixO^S,idir)
      end do
      
    endif


  end subroutine get_lorentz

  subroutine  add_source_lorentz_work(qdt,ixI^L,ixO^L,w,wCT,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: a(ixI^S,3), b(ixI^S,3)
    
    call get_lorentz(ixI^L, ixO^L,wCT,a)
    call twofl_get_v_c(wCT,x,ixI^L,ixO^L,b)
    w(ixO^S,e_c_)=w(ixO^S,e_c_)+qdt*sum(a(ixO^S,1:3)*b(ixO^S,1:3),dim=ndim+1)


  end subroutine  add_source_lorentz_work


  !> Compute 1/(1+v_A^2/c^2) for Boris' approximation, where v_A is the Alfven
  !> velocity
  subroutine twofl_gamma2_alfven(ixI^L, ixO^L, w, gamma_A2)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision, intent(out) :: gamma_A2(ixO^S)
    double precision              :: rhoc(ixI^S) 
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision              :: rhon(ixI^S)
#endif
    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)

#if !defined(ONE_FLUID) || ONE_FLUID==0
    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
    rhoc(ixO^S) = rhon(ixO^S) + rhoc(ixO^S)
#endif
    if (twofl_boris_c < 0.0d0) then
      ! Good for testing the non-conservative momentum treatment
      gamma_A2 = 1.0d0
    else
      ! Compute the inverse of 1 + B^2/(rho * c^2)
      gamma_A2 = 1.0d0 / (1.0d0 + twofl_mag_en_all(w, ixI^L, ixO^L) / (rhoc(ixO^S) * twofl_boris_c**2))
    end if
  end subroutine twofl_gamma2_alfven

  !> Compute 1/sqrt(1+v_A^2/c^2) for Boris simplification, where v_A is the
  !> Alfven velocity
  function twofl_gamma_alfven(w, ixI^L, ixO^L) result(gamma_A)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: gamma_A(ixO^S)

    call twofl_gamma2_alfven(ixI^L, ixO^L, w, gamma_A)
    gamma_A = sqrt(gamma_A)
  end function twofl_gamma_alfven

#if !defined(ONE_FLUID) || ONE_FLUID==0
  !> Calculate v_n vector
  subroutine twofl_get_v_n(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)
    double precision              :: rhon(ixI^S)
    integer :: idir

    call get_rhon_tot(w,ixI^L,ixO^L,rhon)

    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom_n(idir)) / rhon(ixO^S)
    end do

  end subroutine twofl_get_v_n

  subroutine get_rhon_tot(w,ixI^L,ixO^L,rhon)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(out) :: rhon(ixI^S)
    if(has_equi_rho_n0) then
      rhon(ixO^S) = w(ixO^S,rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,block%iw0)
    else  
      rhon(ixO^S) = w(ixO^S,rho_n_) 
    endif

  end subroutine get_rhon_tot



  !> Calculate v component
  subroutine twofl_get_v_n_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)
    double precision              :: rhon(ixI^S)

    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
    v(ixO^S) = w(ixO^S, mom_n(idim)) / rhon(ixO^S)

  end subroutine twofl_get_v_n_idim

  subroutine internal_energy_add_source_n(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call twofl_get_pthermal_n(wCT,x,ixI^L,ixO^L,pth)
    if(has_equi_pe_n0) then
      ! usually block%iw0 should be used
      ! here for sure it is the value at the center
      pth(ixI^S) = pth(ixI^S) + block%equi_vars(ixI^S,equi_pe_n0_,0)
    endif
    call twofl_get_v_n(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-pth,w,x,e_n_)

    if(fix_small_values .and. .not. has_equi_pe_n0) then
      call twofl_handle_small_ei(w,x,ixI^L,ixO^L,e_n_,'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source_n
#endif

  !> Calculate v_c vector
  subroutine twofl_get_v_c(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)
    double precision              :: rhoc(ixI^S)
    integer :: idir

    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom_c(idir)) / rhoc(ixO^S)
    end do

  end subroutine twofl_get_v_c


  subroutine get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(out) :: rhoc(ixI^S)
    if(has_equi_rho_c0) then
      rhoc(ixO^S) = w(ixO^S,rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,block%iw0)
    else  
      rhoc(ixO^S) = w(ixO^S,rho_c_) 
    endif

  end subroutine get_rhoc_tot


  !> Calculate v_c component
  subroutine twofl_get_v_c_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)
    double precision              :: rhoc(ixI^S)

    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
    v(ixO^S) = w(ixO^S, mom_c(idim)) / rhoc(ixO^S)

  end subroutine twofl_get_v_c_idim



  subroutine internal_energy_add_source_c(qdt,ixI^L,ixO^L,wCT,w,x,ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call twofl_get_pthermal_c(wCT,x,ixI^L,ixO^L,pth)
    if(has_equi_pe_c0) then
      pth(ixI^S) = pth(ixI^S) + block%equi_vars(ixI^S,equi_pe_c0_,0)
    endif
    call twofl_get_v_c(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-pth,w,x,ie)
    if(fix_small_values .and. .not. has_equi_pe_c0) then
      call twofl_handle_small_ei(w,x,ixI^L,ixO^L,ie,'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source_c


  !> handle small or negative internal energy
  subroutine twofl_handle_small_ei(w, x, ixI^L, ixO^L, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L, ie
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision              :: rhoc(ixI^S)
#if !defined(ONE_FLUID) || ONE_FLUID==0
    double precision              :: rhon(ixI^S)
#endif

    !TODO
    return

    flag=.false.
    where(w(ixO^S,ie)<small_e) flag(ixO^S,ie)=.true.
    if(any(flag(ixO^S,ie))) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, ie)
      case default
        ! small values error shows primitive variables
        ! TODO is this what was meant to be?
#if !defined(ONE_FLUID) || ONE_FLUID==0
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
        call get_rhon_tot(w,ixI^L,ixO^L,rhon)
#endif
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
        call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)
        do idir = 1, ndir
#if !defined(ONE_FLUID) || ONE_FLUID==0
           w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/rhon(ixO^S)
#endif
           w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/rhoc(ixO^S)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine twofl_handle_small_ei

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a(ixI^S,3), b(ixI^S,3), axb(ixI^S,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixO^S,1:ndir)=block%B0(ixO^S,1:ndir,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixO^S,idir)=block%J0(ixO^S,idir)
      end do
      call cross_product(ixI^L,ixO^L,a,b,axb)
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
      ! add J0xB0 source term in momentum equations
      w(ixO^S,mom_c(1:ndir))=w(ixO^S,mom_c(1:ndir))+axb(ixO^S,1:ndir)
    end if

    if(phys_total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixO^S,:)=wCT(ixO^S,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixO^S,:)=b(ixO^S,:)+block%B0(ixO^S,:,0)
      ! store velocity in a
      do idir=1,ndir
        call twofl_get_v_c_idim(wCT,x,ixI^L,ixO^L,idir,a(ixI^S,idir))
      end do
      call cross_product(ixI^L,ixO^L,a,b,axb)
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixO^S,e_c_)=w(ixO^S,e_c_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
      end do
#if defined(ONE_FLUID) && ONE_FLUID == 1
      if(twofl_ambipolar) then
        !reuse axb
        call twofl_get_jxbxb(wCT,x,ixI^L,ixO^L,axb)
        ! calcuate electric field on cell edges from cell centers
        do idir=7-2*ndim,3
          !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
          !jxbxb(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
          call multiplyAmbiCoef(ixI^L,ixO^L,axb(ixI^S,idir),wCT,x)   
          w(ixO^S,e_c_)=w(ixO^S,e_c_)+axb(ixO^S,idir)*block%J0(ixO^S,idir)
        enddo
      endif
#endif
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_B0')

  end subroutine add_source_B0split

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ixA^L,idir,jdir,kdir,idirmin,idim,jxO^L,hxO^L,ix
    integer :: lxO^L, kxO^L

    double precision :: tmp(ixI^S),tmp2(ixI^S)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    double precision :: gradeta(ixI^S,1:ndim), Bf(ixI^S,1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (twofl_4th_order) then
      ixA^L=ixO^L^LADD2;
    else
      ixA^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    if (twofl_eta>zero)then
       eta(ixA^S)=twofl_eta
       gradeta(ixO^S,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixI^L,ixO^L,idim,tmp)
          gradeta(ixO^S,idim)=tmp(ixO^S)
       end do
    end if

    if(B0field) then
      Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (twofl_4th_order) then
         tmp(ixO^S)=zero
         tmp2(ixI^S)=Bf(ixI^S,idir)
         do idim=1,ndim
            lxO^L=ixO^L+2*kr(idim,^D);
            jxO^L=ixO^L+kr(idim,^D);
            hxO^L=ixO^L-kr(idim,^D);
            kxO^L=ixO^L-2*kr(idim,^D);
            tmp(ixO^S)=tmp(ixO^S)+&
                 (-tmp2(lxO^S)+16.0d0*tmp2(jxO^S)-30.0d0*tmp2(ixO^S)+16.0d0*tmp2(hxO^S)-tmp2(kxO^S)) &
                 /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixO^S)=zero
         tmp2(ixI^S)=Bf(ixI^S,idir)
         do idim=1,ndim
            jxO^L=ixO^L+kr(idim,^D);
            hxO^L=ixO^L-kr(idim,^D);
            tmp(ixO^S)=tmp(ixO^S)+&
                 (tmp2(jxO^S)-2.0d0*tmp2(ixO^S)+tmp2(hxO^S))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixO^S)=tmp(ixO^S)*eta(ixO^S)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (twofl_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixO^S)=tmp(ixO^S)-gradeta(ixO^S,jdir)*current(ixO^S,kdir)
                else
                   tmp(ixO^S)=tmp(ixO^S)+gradeta(ixO^S,jdir)*current(ixO^S,kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixO^S,mag(idir))=w(ixO^S,mag(idir))+qdt*tmp(ixO^S)
       if (phys_energy) then
          w(ixO^S,e_c_)=w(ixO^S,e_c_)+qdt*tmp(ixO^S)*Bf(ixO^S,idir)
          if(phys_solve_eaux) then
            w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)+qdt*tmp(ixO^S)*Bf(ixO^S,idir)
          end if
       end if
    end do ! idir

    if (phys_energy) then
       ! de/dt+=eta*J**2
      tmp(ixO^S)=qdt*eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)+tmp(ixO^S)
      if(phys_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')

  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S),curlj(ixI^S,1:3)
    double precision :: tmpvec(ixI^S,1:3),tmp(ixO^S)
    integer :: ixA^L,idir,idirmin,idirmin1

    ixA^L=ixO^L^LADD2;

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res2: Non-conforming input limits")

    ixA^L=ixO^L^LADD1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixI^L,ixA^L,idirmin,current)

    if (twofl_eta>zero)then
       eta(ixA^S)=twofl_eta
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixA^S,1:ndir)=zero
    do idir=idirmin,3
       tmpvec(ixA^S,idir)=current(ixA^S,idir)*eta(ixA^S)
    end do
    curlj=0.d0
    call curlvector(tmpvec,ixI^L,ixO^L,curlj,idirmin1,1,3)
    if(stagger_grid.and.ndim==2.and.ndir==3) then
      ! if 2.5D
      w(ixO^S,mag(ndir)) = w(ixO^S,mag(ndir))-qdt*curlj(ixO^S,ndir)
    else
      w(ixO^S,mag(1:ndir)) = w(ixO^S,mag(1:ndir))-qdt*curlj(ixO^S,1:ndir)
    end if

    if(phys_energy) then
      ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
      ! de1/dt= eta J^2 - B1 dot curl(eta J)
      tmp(ixO^S)=eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)+qdt*(tmp(ixO^S)-&
        sum(wCT(ixO^S,mag(1:ndir))*curlj(ixO^S,1:ndir),dim=ndim+1))
      if(phys_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res2')
  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    double precision                :: current(ixI^S,7-2*ndir:3)
    double precision                :: tmpvec(ixI^S,1:3),tmpvec2(ixI^S,1:3),tmp(ixI^S),ehyper(ixI^S,1:3)
    integer                         :: ixA^L,idir,jdir,kdir,idirmin,idirmin1

    ixA^L=ixO^L^LADD3;
    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixI^L,ixA^L,idirmin,current)
    tmpvec(ixA^S,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixA^S,jdir)=current(ixA^S,jdir)
    end do

    ixA^L=ixO^L^LADD2;
    call curlvector(tmpvec,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    ixA^L=ixO^L^LADD1;
    tmpvec(ixA^S,1:ndir)=zero
    call curlvector(tmpvec2,ixI^L,ixA^L,tmpvec,idirmin1,1,3)
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*twofl_eta_hyper

    ixA^L=ixO^L;
    tmpvec2(ixA^S,1:ndir)=zero
    call curlvector(ehyper,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-tmpvec2(ixO^S,idir)*qdt
    end do

    if (phys_energy) then
      ! de/dt= +div(B x Ehyper)
      ixA^L=ixO^L^LADD1;
      tmpvec2(ixA^S,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixA^S,idir) = tmpvec(ixA^S,idir)&
        + lvc(idir,jdir,kdir)*wCT(ixA^S,mag(jdir))*ehyper(ixA^S,kdir)
      end do; end do; end do
      tmp(ixO^S)=zero
      call divvector(tmpvec2,ixI^L,ixO^L,tmp)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)+tmp(ixO^S)*qdt
    end if

    if (fix_small_values)  call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S)

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, twofl_divb_4thorder)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (twofl_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(twofl_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*twofl_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*twofl_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
      end if
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       case("limited")
          call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       end select
       if (phys_total_energy) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixO^S,e_c_) = w(ixO^S,e_c_)-qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom_c(idir))=w(ixO^S,mom_c(idir))-qdt*twofl_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm')

  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S),v(ixI^S,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, twofl_divb_4thorder)

    ! calculate velocity
    call twofl_get_v_c(wCT,x,ixI^L,ixO^L,v)

    if (phys_total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixO^S,e_c_)=w(ixO^S,e_c_)-&
           qdt*sum(v(ixO^S,:)*wCT(ixO^S,mag(:)),dim=ndim+1)*divb(ixO^S)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom_c(idir))=w(ixO^S,mom_c(idir))-qdt*twofl_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, twofl_divb_4thorder)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*wCT(ixO^S,mom_c(idir))/wCT(ixO^S,rho_c_)*divb(ixO^S)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: idim, idir, ixp^L, i^D, iside
    double precision :: divb(ixI^S),graddivb(ixI^S)
    logical, dimension(-1:1^D&) :: leveljump

    ! Calculate div B
    ixp^L=ixO^L^LADD1;
    call get_divb(wCT,ixI^L,ixp^L,divb, twofl_divb_4thorder)

    ! for AMR stability, retreat one cell layer from the boarders of level jump
    {do i^DB=-1,1\}
      if(i^D==0|.and.) cycle
      if(neighbor_type(i^D,saveigrid)==2 .or. neighbor_type(i^D,saveigrid)==4) then
        leveljump(i^D)=.true.
      else
        leveljump(i^D)=.false.
      end if
    {end do\}

    ixp^L=ixO^L;
    do idim=1,ndim
      select case(idim)
       {case(^D)
          do iside=1,2
            i^DD=kr(^DD,^D)*(2*iside-3);
            if (leveljump(i^DD)) then
              if (iside==1) then
                ixpmin^D=ixOmin^D-i^D
              else
                ixpmax^D=ixOmax^D-i^D
              end if
            end if
          end do
       \}
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixI^L,ixp^L,idim,graddivb)
       case("limited")
         call gradientS(divb,ixI^L,ixp^L,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
       else
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff &
                          /(^D&1.0d0/block%ds(ixp^S,^D)**2+)
       end if

       w(ixp^S,mag(idim))=w(ixp^S,mag(idim))+graddivb(ixp^S)

       if (typedivbdiff=='all' .and. phys_total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_c_)=w(ixp^S,e_c_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixI^L,ixO^L,divb, fourthorder)

    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: divb(ixI^S)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: bvec(ixI^S,1:ndir)
    double precision                   :: divb_corner(ixI^S), sign
    double precision                   :: aux_vol(ixI^S)
    integer                            :: ixC^L, idir, ic^D, ix^L

    if(stagger_grid) then
      divb=0.d0
      do idir=1,ndim
        ixC^L=ixO^L-kr(idir,^D);
        divb(ixO^S)=divb(ixO^S)+block%ws(ixO^S,idir)*block%surfaceC(ixO^S,idir)-&
                                block%ws(ixC^S,idir)*block%surfaceC(ixC^S,idir)
      end do
      divb(ixO^S)=divb(ixO^S)/block%dvolume(ixO^S)
    else
      bvec(ixI^S,:)=w(ixI^S,mag(:))
      select case(typediv)
      case("central")
        call divvector(bvec,ixI^L,ixO^L,divb,fourthorder)
      case("limited")
        call divvectorS(bvec,ixI^L,ixO^L,divb)
      end select
    end if

  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixI^L,ixO^L,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision                   :: divb(ixI^S), dsurface(ixI^S)

    integer :: ixA^L,idims

    call get_divb(w,ixI^L,ixO^L,divb)
    if(slab_uniform) then
      divb(ixO^S)=0.5d0*abs(divb(ixO^S))/sqrt(twofl_mag_en_all(w,ixI^L,ixO^L))/sum(1.d0/dxlevel(:))
    else
      ixAmin^D=ixOmin^D-1;
      ixAmax^D=ixOmax^D-1;
      dsurface(ixO^S)= sum(block%surfaceC(ixO^S,:),dim=ndim+1)
      do idims=1,ndim
        ixA^L=ixO^L-kr(idims,^D);
        dsurface(ixO^S)=dsurface(ixO^S)+block%surfaceC(ixA^S,idims)
      end do
      divb(ixO^S)=abs(divb(ixO^S))/sqrt(twofl_mag_en_all(w,ixI^L,ixO^L))*&
      block%dvolume(ixO^S)/dsurface(ixO^S)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixO^L, ixI^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    integer, intent(out) :: idirmin
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),bvec(ixI^S,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))

    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixO^S,idirmin0:3)=current(ixO^S,idirmin0:3)+&
        block%J0(ixO^S,idirmin0:3)

  end subroutine get_current

  ! copied from gravity
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    double precision       :: vel(ixI^S)
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    if(qsourcesplit .eqv. grav_split) then
      active = .true.

      if (.not. associated(usr_gravity)) then
        write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
        write(*,*) "like the phys_gravity in mod_usr_methods.t"
        call mpistop("gravity_add_source: usr_gravity not defined")
      else
        call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
      end if
  
      do idim = 1, ndim
#if !defined(ONE_FLUID) || ONE_FLUID==0
        w(ixO^S,mom_n(idim)) = w(ixO^S,mom_n(idim)) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,rho_n_)
#endif
        w(ixO^S,mom_c(idim)) = w(ixO^S,mom_c(idim)) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,rho_c_)
        if(energy) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
#if !defined(ONE_FLUID) || ONE_FLUID==0
          call twofl_get_v_n_idim(wCT,x,ixI^L,ixO^L,idim,vel)
          w(ixO^S,e_n_)=w(ixO^S,e_n_) &
              + qdt * gravity_field(ixO^S,idim) * vel(ixO^S) * wCT(ixO^S,rho_n_)
#endif
          call twofl_get_v_c_idim(wCT,x,ixI^L,ixO^L,idim,vel)
          w(ixO^S,e_c_)=w(ixO^S,e_c_) &
              + qdt * gravity_field(ixO^S,idim) * vel(ixO^S) * wCT(ixO^S,rho_c_)
#else
#if !defined(ONE_FLUID) || ONE_FLUID==0
          call twofl_get_v_n_idim(wCT,x,ixI^L,ixO^L,idim,vel)
          w(ixO^S,e_n_)=w(ixO^S,e_n_) &
              + qdt * gravity_field(ixO^S,idim) *  wCT(ixO^S,mom_n(idim))
#endif
          call twofl_get_v_c_idim(wCT,x,ixI^L,ixO^L,idim,vel)
          w(ixO^S,e_c_)=w(ixO^S,e_c_) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,mom_c(idim))
#endif


        end if
      end do
    end if

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: dxinv(1:ndim), max_grav
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    ^D&dxinv(^D)=one/dx^D;

    if(.not. associated(usr_gravity)) then
      write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
      write(*,*) "like the phys_gravity in mod_usr_methods.t"
      call mpistop("gravity_get_dt: usr_gravity not defined")
    else
      call usr_gravity(ixI^L,ixO^L,w,x,gravity_field)
    end if

    do idim = 1, ndim
      max_grav = maxval(abs(gravity_field(ixO^S,idim)))
      max_grav = max(max_grav, epsilon(1.0d0))
      dtnew = min(dtnew, 1.0d0 / sqrt(max_grav * dxinv(idim)))
    end do

  end subroutine gravity_get_dt


  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine twofl_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    !use mod_radiative_cooling, only: cooling_get_dt
    !use mod_viscosity, only: viscosity_get_dt
    !use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)

    dtnew = bigdouble

    ^D&dxarr(^D)=dx^D;
    if (twofl_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/twofl_eta
    else if (twofl_eta<zero)then
       call get_current(w,ixI^L,ixO^L,idirmin,current)
       call usr_special_resistivity(w,ixI^L,ixO^L,idirmin,x,current,eta)
       dtnew=bigdouble
       do idim=1,ndim
         if(slab_uniform) then
           dtnew=min(dtnew,&
                dtdiffpar/(smalldouble+maxval(eta(ixO^S)/dxarr(idim)**2)))
         else
           dtnew=min(dtnew,&
                dtdiffpar/(smalldouble+maxval(eta(ixO^S)/block%ds(ixO^S,idim)**2)))
         end if
       end do
    end if

    if(twofl_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/twofl_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/twofl_eta_hyper,dtnew)
      end if
    end if


#if !defined(ONE_FLUID) || ONE_FLUID==0
    ! the timestep related to coll terms: 1/(rho_n rho_c alpha)
    if(dtcollpar>0d0 .and. twofl_alpha_coll >0d0) then
        call coll_get_dt(w,x,ixI^L,ixO^L,dtnew)
    endif
#else
    if(twofl_ambipolar_exp) then
      dtnew=min(dtdiffpar*get_ambipolar_dt(w,ixI^L,ixO^L,dx^D,x),dtnew)
    endif
#endif

!    if(twofl_radiative_cooling) then
!      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
!    end if
!
!    if(twofl_viscosity) then
!      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
!    end if
!
    if(twofl_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

  end subroutine twofl_get_dt


#if !defined(ONE_FLUID) || ONE_FLUID==0
  subroutine coll_get_dt(w,x,ixI^L,ixO^L,dtnew)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: dtnew

    double precision :: rhon(ixI^S), rhoc(ixI^S)

    call get_rhon_tot(w,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,ixI^L,ixO^L,rhoc)

    dtnew = min(dtcollpar/(twofl_alpha_coll * maxval(max(rhon(ixO^S), rhoc(ixO^S)))), dtnew)


  end subroutine coll_get_dt
#endif
  

  ! Add geometrical source terms to w
  subroutine twofl_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

!    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
!    br_=mag(1); bphi_=mag(1)-1+phi_
!
!    select case (coordinate)
!    case (cylindrical)
!      if (angmomfix) then
!        call mpistop("angmomfix not implemented yet in MHD")
!      endif
!      call twofl_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
!      if(phi_>0) then
!        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(tmp(ixO^S)-&
!                  wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_))
!        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
!                 -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_) &
!                 +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
!        if(.not.stagger_grid) then
!          w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt/x(ixO^S,1)*&
!                   (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
!                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
!                   /wCT(ixO^S,rho_)
!        end if
!      else
!        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*tmp(ixO^S)
!      end if
!      if(twofl_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
!    case (spherical)
!       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
!       call twofl_get_p_total(wCT,x,ixI^L,ixO^L,tmp1)
!       tmp(ixO^S)=tmp1(ixO^S)
!       if(B0field) then
!         tmp2(ixO^S)=sum(block%B0(ixO^S,:,0)*wCT(ixO^S,mag(:)),dim=ndim+1)
!         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
!       end if
!       ! m1
!       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
!                  *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
!       if(ndir>1) then
!         do idir=2,ndir
!           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2/wCT(ixO^S,rho_)-wCT(ixO^S,mag(idir))**2
!           if(B0field) tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,idir,0)*wCT(ixO^S,mag(idir))
!         end do
!       end if
!       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)/x(ixO^S,1)
!       ! b1
!       if(twofl_glm) then
!         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt/x(ixO^S,1)*2.0d0*wCT(ixO^S,psi_)
!       end if
!
!       {^NOONED
!       ! m2
!       tmp(ixO^S)=tmp1(ixO^S)
!       if(B0field) then
!         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
!       end if
!       ! This will make hydrostatic p=const an exact solution
!       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S) &
!            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
!            /block%dvolume(ixO^S)
!       tmp(ixO^S)=-(wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))/wCT(ixO^S,rho_) &
!            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
!       if (B0field) then
!          tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(2)) &
!               +wCT(ixO^S,mag(1))*block%B0(ixO^S,2,0)
!       end if
!       if(ndir==3) then
!         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2/wCT(ixO^S,rho_) &
!              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
!         if (B0field) then
!            tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,3,0)*wCT(ixO^S,mag(3))&
!                 *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
!         end if
!       end if
!       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
!       ! b2
!       if(.not.stagger_grid) then
!         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
!              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))/wCT(ixO^S,rho_)
!         if(B0field) then
!           tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,2,0) &
!                -wCT(ixO^S,mom(2))*block%B0(ixO^S,1,0))/wCT(ixO^S,rho_)
!         end if
!         if(twofl_glm) then
!           tmp(ixO^S)=tmp(ixO^S) &
!                + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
!         end if
!         w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
!       end if
!       }
!
!       if(ndir==3) then
!         ! m3
!         if(.not.angmomfix) then
!           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))/wCT(ixO^S,rho_) &
!                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
!                -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))/wCT(ixO^S,rho_) &
!                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
!                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
!           if (B0field) then
!              tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(3)) &
!                   +wCT(ixO^S,mag(1))*block%B0(ixO^S,3,0) {^NOONED &
!                   +(block%B0(ixO^S,2,0)*wCT(ixO^S,mag(3)) &
!                   +wCT(ixO^S,mag(2))*block%B0(ixO^S,3,0)) &
!                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
!           end if
!           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
!         else
!           call mpistop("angmomfix not implemented yet in MHD")
!         end if
!         ! b3
!         if(.not.stagger_grid) then
!           tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
!                -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))/wCT(ixO^S,rho_) {^NOONED &
!                -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
!                -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
!                /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
!           if (B0field) then
!              tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,3,0) &
!                   -wCT(ixO^S,mom(3))*block%B0(ixO^S,1,0))/wCT(ixO^S,rho_){^NOONED &
!                   -(wCT(ixO^S,mom(3))*block%B0(ixO^S,2,0) &
!                   -wCT(ixO^S,mom(2))*block%B0(ixO^S,3,0))*dcos(x(ixO^S,2)) &
!                   /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
!           end if
!           w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
!         end if
!       end if
!    end select
  end subroutine twofl_add_source_geom

  !> Compute 2 times total magnetic energy
  function twofl_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,block%iw0))**2, dim=ndim+1)
    else
      mge = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if
  end function twofl_mag_en_all

  !> Compute full magnetic field by direction
  function twofl_mag_i_all(w, ixI^L, ixO^L,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idir
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mgf(ixO^S)

    if (B0field) then
      mgf = w(ixO^S, mag(idir))+block%B0(ixO^S,idir,block%iw0)
    else
      mgf = w(ixO^S, mag(idir))
    end if
  end function twofl_mag_i_all

  !> Compute evolving magnetic energy
  function twofl_mag_en(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    mge = 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
  end function twofl_mag_en

#if !defined(ONE_FLUID) || ONE_FLUID==0
  !> compute kinetic energy of neutrals
  function twofl_kin_en_n(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    if(has_equi_rho_n0) then
      ke = 0.5d0 * sum(w(ixO^S, mom_n(:))**2, dim=ndim+1) / (w(ixO^S, rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,0))
    else
      ke = 0.5d0 * sum(w(ixO^S, mom_n(:))**2, dim=ndim+1) / w(ixO^S, rho_n_)
    endif

  end function twofl_kin_en_n
#endif

  !> compute kinetic energy of charges
  !> w are conserved variables
  function twofl_kin_en_c(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    if(has_equi_rho_c0) then
      ke = 0.5d0 * sum(w(ixO^S, mom_c(:))**2, dim=ndim+1) / (w(ixO^S, rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,0))
    else
      ke = 0.5d0 * sum(w(ixO^S, mom_c(:))**2, dim=ndim+1) / w(ixO^S, rho_c_)
    endif
  end function twofl_kin_en_c

  subroutine twofl_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: vHall(ixI^S,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: rho(ixI^S)

    if(has_equi_rho_c0) then
      rho(ixO^S) = w(ixO^S,rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,0)
    else
      rho(ixO^S) = w(ixO^S,rho_c_) 
    endif


    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    vHall(ixO^S,1:3) = zero
    vHall(ixO^S,idirmin:3) = - twofl_etah*current(ixO^S,idirmin:3)
    do idir = idirmin, 3
       vHall(ixO^S,idir) = vHall(ixO^S,idir)/w(ixO^S,rho_c_)
    end do

  end subroutine twofl_getv_Hall

#if defined(ONE_FLUID) && ONE_FLUID == 1
  subroutine twofl_get_Jambi(w,x,ixI^L,ixO^L,res)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, allocatable, intent(inout) :: res(:^D&,:)


    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)

    res = 0d0

    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
 
    res(ixO^S,idirmin:3)=-current(ixO^S,idirmin:3)
    do idir = idirmin, 3
      call multiplyAmbiCoef(ixI^L,ixO^L,res(ixI^S,idir),w,x)   
    enddo

  end subroutine twofl_get_Jambi
#endif

! the following not used
!  subroutine twofl_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
!    use mod_global_parameters
!
!    integer, intent(in) :: ixI^L, ixO^L
!    double precision, intent(in)    :: dx^D
!    double precision, intent(in)    :: w(ixI^S,1:nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(out)   :: dthall
!    !.. local ..
!    double precision :: dxarr(ndim)
!    double precision :: bmag(ixI^S)
!
!    dthall=bigdouble
!
!    ! because we have that in cmax now:
!    return
!
!    ^D&dxarr(^D)=dx^D;
!
!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,block%iw0))**2))
!    end if
!
!    if(slab_uniform) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(twofl_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_c_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(twofl_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_c_)))
!    end if
!
!  end subroutine twofl_getdt_Hall

  subroutine twofl_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixI^S), dPsi(ixI^S)

    if(stagger_grid) then
      wLC(ixO^S,mag(idir))=s%ws(ixO^S,idir)
      wRC(ixO^S,mag(idir))=s%ws(ixO^S,idir)
      wLp(ixO^S,mag(idir))=s%ws(ixO^S,idir)
      wRp(ixO^S,mag(idir))=s%ws(ixO^S,idir)
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-MHD system.
      ! 23/04/2013 Oliver Porth
      dB(ixO^S)   = wRp(ixO^S,mag(idir)) - wLp(ixO^S,mag(idir))
      dPsi(ixO^S) = wRp(ixO^S,psi_) - wLp(ixO^S,psi_)

      wLp(ixO^S,mag(idir))   = 0.5d0 * (wRp(ixO^S,mag(idir)) + wLp(ixO^S,mag(idir))) &
           - 0.5d0/cmax_global * dPsi(ixO^S)
      wLp(ixO^S,psi_)       = 0.5d0 * (wRp(ixO^S,psi_) + wLp(ixO^S,psi_)) &
           - 0.5d0*cmax_global * dB(ixO^S)

      wRp(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wRp(ixO^S,psi_) = wLp(ixO^S,psi_)
      wRC(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wRC(ixO^S,psi_) = wLp(ixO^S,psi_)
      wLC(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wLC(ixO^S,psi_) = wLp(ixO^S,psi_)
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine twofl_modify_wLR

  subroutine twofl_boundary_adjust
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid
    integer :: ixG^L, ixO^L, i^D

    ixG^L=ixG^LL;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>ps(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i^D=kr(^D,idim)*(2*iside-3);
              if (neighbor_type(i^D,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                ! MF nonlinear force-free B field extrapolation and data driven
                ! require normal B of the first ghost cell layer to be untouched by
                ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
                select case (idim)
                {case (^D)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin^DD=ixGmax^D+1-nghostcells+boundary_divbfix_skip(2*^D)^D%ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmax^DD;
                   else
                      ! minimal boundary
                      ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmin^D-1+nghostcells-boundary_divbfix_skip(2*^D-1)^D%ixOmax^DD=ixGmax^DD;
                   end if \}
                end select
                call fixdivB_boundary(ixG^L,ixO^L,ps(igrid)%w,ps(igrid)%x,iB)
              end if
           end do
        end do
     end do

  end subroutine twofl_boundary_adjust

  subroutine fixdivB_boundary(ixG^L,ixO^L,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixG^L,ixO^L,iB
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision, intent(in) :: x(ixG^S,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix^D,ixF^L

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=w(ix1+1,ixFmin2:ixFmax2,mag(1)) &
            +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-&
                    w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=( (w(ix1+1,ixFmin2:ixFmax2,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1,ixFmin2:ixFmax2,1)&
           +(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,ixFmin2:ixFmax2,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,2)&
           -(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,2) )&
            /block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,mag(1))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
                     w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) &
             +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-&
                     w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) &
             +dx1x3*(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
          ( (w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)&
           +(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,2)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,3)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)-&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
       }
     case(2)
       {^IFTWOD
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=w(ix1-1,ixFmin2:ixFmax2,mag(1)) &
            -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-&
                    w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=( (w(ix1-1,ixFmin2:ixFmax2,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)&
           -(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,ixFmin2:ixFmax2,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,2)&
           +(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,2) )&
            /block%surfaceC(ix1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,mag(1))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
                     w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) &
             -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-&
                     w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) &
             -dx1x3*(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
          ( (w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)&
           -(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,2)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,3)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)-&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
       }
     case(3)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=w(ixFmin1:ixFmax1,ix2+1,mag(2)) &
            +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))-&
                    w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         enddo
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,2)&
           +(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,1)&
           -(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,1) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)-w(ixFmin1:ixFmax1,ix2,mag(2))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
             ix2+1,ixFmin3:ixFmax3,mag(2)) &
             +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1))) &
             +dx2x3*(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=&
          ( (w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,2)&
           +(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,1)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,1)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,3)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,2)-&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
       }
     case(4)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=w(ixFmin1:ixFmax1,ix2-1,mag(2)) &
            -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))-&
                    w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)&
           -(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,1)&
           +(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,1) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2,2)-w(ixFmin1:ixFmax1,ix2,mag(2))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
             ix2-1,ixFmin3:ixFmax3,mag(2)) &
             -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1))) &
             -dx2x3*(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=&
          ( (w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,2)&
           -(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,1)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,1)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,3)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,2)-&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
       }
     {^IFTHREED
     case(5)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3+1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=w(ixFmin1:ixFmax1,&
             ixFmin2:ixFmax2,ix3+1,mag(3)) &
             +dx3x1*(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1))) &
             +dx3x2*(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))-&
                     w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=&
          ( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,3)&
           +(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,1)&
           -(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,1)&
           +(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,2)&
           -(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,2) )&
            /block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,3)-&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
     case(6)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3-1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=w(ixFmin1:ixFmax1,&
             ixFmin2:ixFmax2,ix3-1,mag(3)) &
             -dx3x1*(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1))) &
             -dx3x2*(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))-&
                     w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=&
          ( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,3)&
           -(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,1)&
           +(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,1)&
           -(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,2)&
           +(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,2) )&
            /block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,3)-&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  {^NOONED
  subroutine twofl_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active
    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ix^L, ixC^L, idim
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixG^T), grad(ixG^T, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb

    mg%operator_type = mg_laplacian

    ! Set boundary conditions
    do n = 1, 2*ndim
       idim = (n+1)/2
       select case (typeboundary(mag(idim), n))
       case ('symm')
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('special')
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case default
          print *, "divb_multigrid warning: unknown b.c.: ", &
               trim(typeboundary(mag(idim), n))
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    ix^L=ixM^LL^LADD1;
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       call get_divb(ps(igrid)%w(ixG^T, 1:nw), ixG^LL, ixM^LL, tmp, &
            twofl_divb_4thorder)
       mg%boxes(id)%cc({1:nc}, mg_irhs) = tmp(ixM^T)
       max_divb = max(max_divb, maxval(abs(tmp(ixM^T))))
    end do

    ! Solve laplacian(phi) = divB
    if(stagger_grid) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION, &
           MPI_MAX, icomm, ierrmpi)

      if (mype == 0) print *, "Performing multigrid divB cleaning"
      if (mype == 0) print *, "iteration vs residual"
      ! Solve laplacian(phi) = divB
      do n = 1, max_its
         call mg_fas_fmg(mg, n>1, max_res=residual_it(n))
         if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
         if (residual_it(n) < residual_reduction * max_divb) exit
      end do
      if (mype == 0 .and. n > max_its) then
         print *, "divb_multigrid warning: not fully converged"
         print *, "current amplitude of divb: ", residual_it(max_its)
         print *, "multigrid smallest grid: ", &
              mg%domain_size_lvl(:, mg%lowest_lvl)
         print *, "note: smallest grid ideally has <= 8 cells"
         print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
         print *, "note: dx/dy/dz should be similar"
      end if
    else
      do n = 1, max_its
         call mg_fas_vcycle(mg, max_res=res)
         if (res < max_residual) exit
      end do
      if (res > max_residual) call mpistop("divb_multigrid: no convergence")
    end if


    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! Compute the gradient of phi
       tmp(ix^S) = mg%boxes(id)%cc({:,}, mg_iphi)

       if(stagger_grid) then
         do idim =1, ndim
           ixCmin^D=ixMlo^D-kr(idim,^D);
           ixCmax^D=ixMhi^D;
           call gradientx(tmp,ps(igrid)%x,ixG^LL,ixC^L,idim,grad(ixG^T,idim),.false.)
           ! Apply the correction B* = B - gradient(phi)
           ps(igrid)%ws(ixC^S,idim)=ps(igrid)%ws(ixC^S,idim)-grad(ixC^S,idim)
         end do
         ! store cell-center magnetic energy
         tmp(ixM^T) = sum(ps(igrid)%w(ixM^T, mag(1:ndim))**2, dim=ndim+1)
         ! change cell-center magnetic field
         call twofl_face_to_center(ixM^LL,ps(igrid))
       else
         do idim = 1, ndim
            call gradient(tmp,ixG^LL,ixM^LL,idim,grad(ixG^T, idim))
         end do
         ! store cell-center magnetic energy
         tmp(ixM^T) = sum(ps(igrid)%w(ixM^T, mag(1:ndim))**2, dim=ndim+1)
         ! Apply the correction B* = B - gradient(phi)
         ps(igrid)%w(ixM^T, mag(1:ndim)) = &
              ps(igrid)%w(ixM^T, mag(1:ndim)) - grad(ixM^T, :)
       end if

       if(phys_total_energy) then
         ! Determine magnetic energy difference
         tmp(ixM^T) = 0.5_dp * (sum(ps(igrid)%w(ixM^T, &
              mag(1:ndim))**2, dim=ndim+1) - tmp(ixM^T))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixM^T, e_c_) = ps(igrid)%w(ixM^T, e_c_) + tmp(ixM^T)
       end if
    end do

    active = .true.

  end subroutine twofl_clean_divb_multigrid
  }

  subroutine twofl_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixI^S,1:nw)
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    case('uct_contact')
      call update_faces_contact(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s)
    case('uct_hll')
      call update_faces_hll(ixI^L,ixO^L,qt,qdt,fE,sCT,s)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine twofl_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer                            :: hxC^L,ixC^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi
#if defined(ONE_FLUID) && ONE_FLUID==1
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_ambi
#endif

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

#if defined(ONE_FLUID) && ONE_FLUID==1
    ! if there is ambipolar diffusion, get E_ambi
    if(twofl_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)
#endif
    fE=zero

    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            hxC^L=ixC^L+kr(idim2,^D);
            ! Interpolate to edges
            fE(ixC^S,idir)=quarter*(fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)&
                                   -fC(ixC^S,iwdim2,idim1)-fC(hxC^S,iwdim2,idim1))

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(twofl_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
#if defined(ONE_FLUID) && ONE_FLUID==1
            ! add ambipolar electric field
            if(twofl_ambipolar_exp) fE(ixC^S,idir)=fE(ixC^S,idir)+E_ambi(ixC^S,idir)
#endif
            fE(ixC^S,idir)=qdt*s%dsC(ixC^S,idir)*fE(ixC^S,idir)

            if (.not.slab) then
              where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixC^S,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) &
      call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate

  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    double precision                   :: circ(ixI^S,1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixI^S,7-2*ndim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixI^S),ER(ixI^S)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixI^S),ERC(ixI^S)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi, E_ambi
    ! total magnetic field at cell centers
    double precision                   :: Btot(ixI^S,1:ndim)
    integer                            :: hxC^L,ixC^L,jxC^L,ixA^L,ixB^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    if(B0field) then
      Btot(ixI^S,1:ndim)=wp(ixI^S,mag(1:ndim))+block%B0(ixI^S,1:ndim,0)
    else
      Btot(ixI^S,1:ndim)=wp(ixI^S,mag(1:ndim))
    end if
    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=7-2*ndim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)+Btot(ixI^S,idim1)*wp(ixI^S,mom_c(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)-Btot(ixI^S,idim1)*wp(ixI^S,mom_c(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)
#if defined(ONE_FLUID) && ONE_FLUID==1
    if(twofl_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)
#endif
    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    fE=zero
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax^D=ixOmax^D;
            ixCmin^D=ixOmin^D+kr(idir,^D)-1;
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            hxC^L=ixC^L+kr(idim2,^D);
            ! average cell-face electric field to cell edges
            fE(ixC^S,idir)=quarter*&
            (fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)&
            -fC(ixC^S,iwdim2,idim1)-fC(hxC^S,iwdim2,idim1))

            ! add slope in idim2 direction from equation (50)
            ixAmin^D=ixCmin^D;
            ixAmax^D=ixCmax^D+kr(idim1,^D);
            EL(ixA^S)=fC(ixA^S,iwdim1,idim2)-ECC(ixA^S,idir)
            hxC^L=ixA^L+kr(idim2,^D);
            ER(ixA^S)=fC(ixA^S,iwdim1,idim2)-ECC(hxC^S,idir)
            where(vnorm(ixC^S,idim1)>0.d0)
              ELC(ixC^S)=EL(ixC^S)
            else where(vnorm(ixC^S,idim1)<0.d0)
              ELC(ixC^S)=EL(jxC^S)
            else where
              ELC(ixC^S)=0.5d0*(EL(ixC^S)+EL(jxC^S))
            end where
            hxC^L=ixC^L+kr(idim2,^D);
            where(vnorm(hxC^S,idim1)>0.d0)
              ERC(ixC^S)=ER(ixC^S)
            else where(vnorm(hxC^S,idim1)<0.d0)
              ERC(ixC^S)=ER(jxC^S)
            else where
              ERC(ixC^S)=0.5d0*(ER(ixC^S)+ER(jxC^S))
            end where
            fE(ixC^S,idir)=fE(ixC^S,idir)+0.25d0*(ELC(ixC^S)+ERC(ixC^S))

            ! add slope in idim1 direction from equation (50)
            jxC^L=ixC^L+kr(idim2,^D);
            ixAmin^D=ixCmin^D;
            ixAmax^D=ixCmax^D+kr(idim2,^D);
            EL(ixA^S)=-fC(ixA^S,iwdim2,idim1)-ECC(ixA^S,idir)
            hxC^L=ixA^L+kr(idim1,^D);
            ER(ixA^S)=-fC(ixA^S,iwdim2,idim1)-ECC(hxC^S,idir)
            where(vnorm(ixC^S,idim2)>0.d0)
              ELC(ixC^S)=EL(ixC^S)
            else where(vnorm(ixC^S,idim2)<0.d0)
              ELC(ixC^S)=EL(jxC^S)
            else where
              ELC(ixC^S)=0.5d0*(EL(ixC^S)+EL(jxC^S))
            end where
            hxC^L=ixC^L+kr(idim1,^D);
            where(vnorm(hxC^S,idim2)>0.d0)
              ERC(ixC^S)=ER(ixC^S)
            else where(vnorm(hxC^S,idim2)<0.d0)
              ERC(ixC^S)=ER(jxC^S)
            else where
              ERC(ixC^S)=0.5d0*(ER(ixC^S)+ER(jxC^S))
            end where
            fE(ixC^S,idir)=fE(ixC^S,idir)+0.25d0*(ELC(ixC^S)+ERC(ixC^S))

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(twofl_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
#if defined(ONE_FLUID) && ONE_FLUID==1
            ! add ambipolar electric field
            if(twofl_ambipolar_exp) fE(ixC^S,idir)=fE(ixC^S,idir)+E_ambi(ixC^S,idir)
#endif
            ! times time step and edge length 
            fE(ixC^S,idir)=fE(ixC^S,idir)*qdt*s%dsC(ixC^S,idir)
            if (.not.slab) then
              where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixC^S,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) &
      call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate

  end subroutine update_faces_contact

  !> update faces
  subroutine update_faces_hll(ixI^L,ixO^L,qt,qdt,fE,sCT,s)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
    type(state)                        :: sCT, s

    double precision                   :: vtilL(ixI^S,2)
    double precision                   :: vtilR(ixI^S,2)
    double precision                   :: bfacetot(ixI^S,ndim)
    double precision                   :: btilL(s%ixGs^S,ndim)
    double precision                   :: btilR(s%ixGs^S,ndim)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi, E_ambi
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir

    associate(bfaces=>s%ws,bfacesCT=>sCT%ws,x=>s%x,vbarC=>vcts%vbarC,cbarmin=>vcts%cbarmin,&
      cbarmax=>vcts%cbarmax)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)
#if defined(ONE_FLUID) && ONE_FLUID==1
    ! if there is ambipolar diffusion, get E_ambi
    if(twofl_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)
#endif
    fE=zero

    do idir=7-2*ndim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);

      ! Set indices and directions
      idim1=mod(idir,3)+1
      idim2=mod(idir+1,3)+1

      jxC^L=ixC^L+kr(idim1,^D);
      ixCp^L=ixC^L+kr(idim2,^D);

      ! Reconstruct transverse transport velocities
      call reconstruct(ixI^L,ixC^L,idim2,vbarC(ixI^S,idim1,1),&
               vtilL(ixI^S,2),vtilR(ixI^S,2))

      call reconstruct(ixI^L,ixC^L,idim1,vbarC(ixI^S,idim2,2),&
               vtilL(ixI^S,1),vtilR(ixI^S,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      if(B0field) then
        bfacetot(ixI^S,idim1)=bfacesCT(ixI^S,idim1)+block%B0(ixI^S,idim1,idim1)
        bfacetot(ixI^S,idim2)=bfacesCT(ixI^S,idim2)+block%B0(ixI^S,idim2,idim2)
      else
        bfacetot(ixI^S,idim1)=bfacesCT(ixI^S,idim1)
        bfacetot(ixI^S,idim2)=bfacesCT(ixI^S,idim2)
      end if
      call reconstruct(ixI^L,ixC^L,idim2,bfacetot(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))

      call reconstruct(ixI^L,ixC^L,idim1,bfacetot(ixI^S,idim2),&
               btilL(ixI^S,idim2),btilR(ixI^S,idim2))

      ! Take the maximum characteristic

      cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
      cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

      cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
      cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
     

      ! Calculate eletric field
      fE(ixC^S,idir)=-(cp(ixC^S,1)*vtilL(ixC^S,1)*btilL(ixC^S,idim2) &
                     + cm(ixC^S,1)*vtilR(ixC^S,1)*btilR(ixC^S,idim2) &
                     - cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))&
                     /(cp(ixC^S,1)+cm(ixC^S,1)) &
                     +(cp(ixC^S,2)*vtilL(ixC^S,2)*btilL(ixC^S,idim1) &
                     + cm(ixC^S,2)*vtilR(ixC^S,2)*btilR(ixC^S,idim1) &
                     - cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))&
                     /(cp(ixC^S,2)+cm(ixC^S,2))

      ! add current component of electric field at cell edges E=-vxB+eta J
      if(twofl_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
#if defined(ONE_FLUID) && ONE_FLUID==1
      ! add ambipolar electric field
      if(twofl_ambipolar_exp) fE(ixC^S,idir)=fE(ixC^S,idir)+E_ambi(ixC^S,idir)
#endif
      fE(ixC^S,idir)=qdt*s%dsC(ixC^S,idir)*fE(ixC^S,idir)

      if (.not.slab) then
        where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixC^S,idir)=zero
        end where
      end if

    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) &
      call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face: interal(fE dot dl)

    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate
  end subroutine update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixI^L,ixO^L,sCT,s,jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixI^L, ixO^L
    type(state), intent(in)            :: sCT, s
    ! current on cell edges
    double precision :: jce(ixI^S,7-2*ndim:3)

    ! current on cell centers
    double precision :: jcc(ixI^S,7-2*ndim:3)
    ! location at cell faces
    double precision :: xs(ixGs^T,1:ndim)
    ! resistivity
    double precision :: eta(ixI^S)
    double precision :: gradi(ixGs^T)
    integer :: ix^D,ixC^L,ixA^L,ixB^L,idir,idirmin,idim1,idim2

    associate(x=>s%x,dx=>s%dx,w=>s%w,wCT=>sCT%w,wCTs=>sCT%ws)
    ! calculate current density at cell edges
    jce=0.d0
    do idim1=1,ndim 
      do idim2=1,ndim
        do idir=7-2*ndim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax^D=ixOmax^D;
          ixCmin^D=ixOmin^D+kr(idir,^D)-1;
          ixBmax^D=ixCmax^D-kr(idir,^D)+1;
          ixBmin^D=ixCmin^D;
          ! current at transverse faces
          xs(ixB^S,:)=x(ixB^S,:)
          xs(ixB^S,idim2)=x(ixB^S,idim2)+half*dx(ixB^S,idim2)
          call gradientx(wCTs(ixGs^T,idim2),xs,ixGs^LL,ixC^L,idim1,gradi,.true.)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixC^S,idir)=jce(ixC^S,idir)+gradi(ixC^S)
          else
            jce(ixC^S,idir)=jce(ixC^S,idir)-gradi(ixC^S)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(twofl_eta>zero)then
      jce(ixI^S,:)=jce(ixI^S,:)*twofl_eta
    else
      ixA^L=ixO^L^LADD1;
      call get_current(wCT,ixI^L,ixO^L,idirmin,jcc)
      call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,jcc,eta)
      ! calcuate eta on cell edges
      do idir=7-2*ndim,3
        ixCmax^D=ixOmax^D;
        ixCmin^D=ixOmin^D+kr(idir,^D)-1;
        jcc(ixC^S,idir)=0.d0
       {do ix^DB=0,1\}
          if({ ix^D==1 .and. ^D==idir | .or.}) cycle
          ixAmin^D=ixCmin^D+ix^D;
          ixAmax^D=ixCmax^D+ix^D;
          jcc(ixC^S,idir)=jcc(ixC^S,idir)+eta(ixA^S)
       {end do\}
        jcc(ixC^S,idir)=jcc(ixC^S,idir)*0.25d0
        jce(ixC^S,idir)=jce(ixC^S,idir)*jcc(ixC^S,idir)
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

#if defined(ONE_FLUID) && ONE_FLUID==1
  !> get ambipolar electric field on cell edges
  subroutine get_ambipolar_electric_field(ixI^L,ixO^L,w,x,fE)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: fE(ixI^S,7-2*ndim:3)

    double precision :: jxbxb(ixI^S,1:3)
    integer :: idir,ixA^L,ixC^L,ix^D

    ixA^L=ixO^L^LADD1;
    call twofl_get_jxbxb(w,x,ixI^L,ixA^L,jxbxb)
    ! calcuate electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
      !jxbxb(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
      call multiplyAmbiCoef(ixI^L,ixA^L,jxbxb(ixI^S,idir),w,x)   
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D+kr(idir,^D)-1;
      fE(ixC^S,idir)=0.d0
     {do ix^DB=0,1\}
        if({ ix^D==1 .and. ^D==idir | .or.}) cycle
        ixAmin^D=ixCmin^D+ix^D;
        ixAmax^D=ixCmax^D+ix^D;
        fE(ixC^S,idir)=fE(ixC^S,idir)+jxbxb(ixA^S,idir)
     {end do\}
      fE(ixC^S,idir)=fE(ixC^S,idir)*0.25d0
    end do

  end subroutine get_ambipolar_electric_field
#endif

  !> calculate cell-center values from face-center values
  subroutine twofl_face_to_center(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim

    associate(w=>s%w, ws=>s%ws)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxO^L=ixO^L-kr(idim,^D);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixO^S,mag(idim))=half/s%surface(ixO^S,idim)*&
        (ws(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
        +ws(hxO^S,idim)*s%surfaceC(hxO^S,idim))
    end do

    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do

    end associate

  end subroutine twofl_face_to_center

  !> calculate magnetic field from vector potential
  subroutine b_from_vector_potential(ixIs^L, ixI^L, ixO^L, ws, x)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws)
    double precision, intent(in)       :: x(ixI^S,1:ndim)

    double precision                   :: Adummy(ixIs^S,1:3)

    call b_from_vector_potentialA(ixIs^L, ixI^L, ixO^L, ws, x, Adummy)

  end subroutine b_from_vector_potential


#if !defined(ONE_FLUID) || ONE_FLUID==0

    !> Implicit solve of psb=psa+dtfactor*dt*F_im(psb)
  subroutine twofl_implicit_coll_terms_update(dtfactor,qdt,qtC,psb,psa)
    use mod_global_parameters
    use mod_ghostcells_update

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    double precision :: tmp(ixG^T),tmp1(ixG^T),tmp2(ixG^T),tmp3(ixG^T),tmp4(ixG^T),tmp5(ixG^T)
    double precision :: v_c(ixG^T,ndir), v_n(ixG^T,ndir)
    integer :: idir
    integer :: iigrid, igrid
    double precision :: rhon(ixG^T), rhoc(ixG^T)
    !print*, "IMPL call ", it

    !$OMP PARALLEL DO PRIVATE(igrid)
    call getbc(global_time,0.d0,psa,1,nw)
    do iigrid=1,igridstail; igrid=igrids(iigrid);

      ! first copy psa into psb
      !psb(igrid)%w(ixG^T,:) = psa(igrid)%w(ixG^T,:)
      ! only the indices which are not updated here: densities and mag. field
      psb(igrid)%w(ixG^T,mag(:)) = psa(igrid)%w(ixG^T,mag(:))
      psb(igrid)%w(ixG^T,rho_n_) = psa(igrid)%w(ixG^T,rho_n_)
      psb(igrid)%w(ixG^T,rho_c_) = psa(igrid)%w(ixG^T,rho_c_)

      !print*, "PSA momn", psa(igrid)%w(ixGlo1:ixGlo1+5,mom_n(1))
      !print*, "PSA momc", psa(igrid)%w(ixGlo1:ixGlo1+5,mom_c(1))

      call get_rhon_tot(psa(igrid)%w,ixG^LL,ixG^LL,rhon)
      call get_rhoc_tot(psa(igrid)%w,ixG^LL,ixG^LL,rhoc)

      tmp3(ixG^T) =  1d0 + dtfactor * dt * twofl_alpha_coll * (rhon(ixG^T) +  rhoc(ixG^T))     
      ! momentum update
      do idir=1,ndir
        tmp(ixG^T) = twofl_alpha_coll * dtfactor * dt *(-rhoc(ixG^T) * psa(igrid)%w(ixG^T,mom_n(idir)) + &
                                          rhon(ixG^T) * psa(igrid)%w(ixG^T,mom_c(idir)))/tmp3(ixG^T)

        psb(igrid)%w(ixG^T,mom_n(idir)) = psa(igrid)%w(ixG^T,mom_n(idir)) + tmp(ixG^T)
        psb(igrid)%w(ixG^T,mom_c(idir)) = psa(igrid)%w(ixG^T,mom_c(idir)) - tmp(ixG^T)
      enddo

!      !print*, "ec ", psa(igrid)%w(1:10,e_c_) 
!      !print*, "en ", psa(igrid)%w(1:10,e_n_) 
!      !energy update
!      ! tmp4 = e int n
!      ! tmp5 = e int c
!      ! contribution from velocity
      if(.not. phys_internal_e) then
        ! kinetic energy update
        tmp1(ixG^T) =  twofl_kin_en_n(psa(igrid)%w,ixG^LL,ixG^LL) 
        tmp2(ixG^T) =  twofl_kin_en_c(psa(igrid)%w,ixG^LL,ixG^LL) 
        tmp4(ixG^T) = psa(igrid)%w(ixG^T,e_n_) - tmp1(ixG^T)
        tmp5(ixG^T) = psa(igrid)%w(ixG^T,e_c_) - tmp2(ixG^T)
      !print*, "ec int ", tmp4(1:10) 
      !print*, "en int ", tmp5(1:10) 
        if(phys_total_energy) then
          tmp5(ixG^T) = tmp5(ixG^T) - twofl_mag_en(psa(igrid)%w,ixG^LL,ixG^LL)
        endif

        !!implicit update
        tmp(ixG^T) = twofl_alpha_coll * dtfactor * dt * &
                           (-rhoc(ixG^T) * tmp1(ixG^T) + rhon(ixG^T) * tmp2(ixG^T))/tmp3(ixG^T)
        psb(igrid)%w(ixG^T,e_n_) = psa(igrid)%w(ixG^T,e_n_) + tmp(ixG^T)
        psb(igrid)%w(ixG^T,e_c_) = psa(igrid)%w(ixG^T,e_c_) - tmp(ixG^T)

       !explicit update 
!        ! calculate velocities, using the already updated variables
!        call twofl_get_v_n(psb(igrid)%w,psb(igrid)%x,ixG^LL,ixG^LL,v_n)
!        call twofl_get_v_c(psb(igrid)%w,psb(igrid)%x,ixG^LL,ixG^LL,v_c)
!        tmp(ixG^T) = 0.5d0 * (sum(v_c(ixG^T,1:ndir), dim=ndim+1) - sum(v_n(ixG^T,1:ndir)**2, dim=ndim+1)) &
!                     * dtfactor * dt * twofl_alpha_coll * rhoc(ixG^T) * rhon(ixG^T)
!        
!        psb(igrid)%w(ixG^T,e_n_) = psa(igrid)%w(ixG^T,e_n_) + tmp(ixG^T)
!        psb(igrid)%w(ixG^T,e_c_) = psa(igrid)%w(ixG^T,e_c_) - tmp(ixG^T)

       else 
        tmp4(ixG^T) = psa(igrid)%w(ixG^T,e_n_) 
        tmp5(ixG^T) = psa(igrid)%w(ixG^T,e_c_) 
        ! calculate velocities, using the already updated variables
        call twofl_get_v_n(psb(igrid)%w,psb(igrid)%x,ixG^LL,ixG^LL,v_n)
        call twofl_get_v_c(psb(igrid)%w,psb(igrid)%x,ixG^LL,ixG^LL,v_c)
        tmp(ixG^T) = 0.5d0 * sum((v_c(ixG^T,1:ndir) - v_n(ixG^T,1:ndir))**2, dim=ndim+1) &
                     * dtfactor * dt * twofl_alpha_coll * rhoc(ixG^T) * rhon(ixG^T)
        psb(igrid)%w(ixG^T,e_n_) = psa(igrid)%w(ixG^T,e_n_) + tmp(ixG^T)
        psb(igrid)%w(ixG^T,e_c_) = psa(igrid)%w(ixG^T,e_c_) + tmp(ixG^T)
       endif

      !update internal energy
      if(twofl_coll_inc_te) then
!        if(has_equi_pe_n0) then
!          tmp4(ixG^T) = tmp4(ixG^T) + psa(igrid)%equi_vars(ixG^T,equi_pe_n0_)*inv_gamma_1  
!        endif
!        if(has_equi_pe_c0) then
!          tmp5(ixG^T) = tmp5(ixG^T) + psa(igrid)%equi_vars(ixG^T,equi_pe_c0_)*inv_gamma_1 
!        endif
        tmp3(ixG^T) =  1d0 + dtfactor * dt * twofl_alpha_coll * (rhon(ixG^T)/Rc +  rhoc(ixG^T)/Rn) !2 from braginskii    
        tmp(ixG^T) = dtfactor * dt * twofl_alpha_coll *(-rhoc(ixG^T)/Rn * tmp4(ixG^T) + rhon(ixG^T)/Rc * tmp5(ixG^T))/tmp3(ixG^T)
        psb(igrid)%w(ixG^T,e_n_) = psb(igrid)%w(ixG^T,e_n_)+tmp(ixG^T)
        psb(igrid)%w(ixG^T,e_c_) = psb(igrid)%w(ixG^T,e_c_)-tmp(ixG^T)
      endif
    end do
    !$OMP END PARALLEL DO

   end subroutine twofl_implicit_coll_terms_update 






  subroutine twofl_explicit_coll_terms_update(qdt,ixI^L,ixO^L,w,wCT,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: wCT(ixI^S,1:nw)

    double precision :: tmp(ixI^S),tmp2(ixI^S) 
    double precision, dimension(:^D&,:), allocatable :: vc, vn
    double precision :: rhon(ixI^S), rhoc(ixI^S)
    integer :: idir
!    print*, "expl call ", it

!    print*, " RHOC ", wCT(ixI^S,rho_c_) 
!    print*, " RHON ", wCT(ixI^S,rho_n_) 
!    print*, " MOM N ", wCT(ixI^S,mom_n(idir))
!    print*, " MOM C ", wCT(ixO^S,mom_c(idir))
    if(has_equi_rho_n0) then
      rhon(ixO^S) = wCT(ixO^S,rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,0)
    else  
      rhon(ixO^S) = wCT(ixO^S,rho_n_) 
    endif

    if(has_equi_rho_c0) then
      rhoc(ixO^S) = wCT(ixO^S,rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_)
    else  
      rhoc(ixO^S) = wCT(ixO^S,rho_c_) 
    endif

    ! momentum update
    do idir=1,ndir
      ! the coll. term in the neutrals momentum eq. multiplied by qdt
      tmp(ixO^S) = qdt * twofl_alpha_coll * (-rhoc(ixO^S) * wCT(ixO^S,mom_n(idir)) + rhon(ixO^S) * wCT(ixO^S,mom_c(idir)))
      w(ixO^S,mom_n(idir)) = w(ixO^S,mom_n(idir)) + tmp(ixO^S) 
      w(ixO^S,mom_c(idir)) = w(ixO^S,mom_c(idir)) - tmp(ixO^S) 

    enddo
    ! energy update
    ! calculate velocities
    allocate(vc(ixO^S,1:ndir),vn(ixO^S,1:ndir))

    !use the updated vars
    call twofl_get_v_c(w,x,ixO^L,ixO^L,vc)
    call twofl_get_v_n(w,x,ixO^L,ixO^L,vn)
    if(twofl_eq_energy == EQ_ENERGY_INT) then
      ! tmp = qdt * (v_c-v_n)^2 * alpha * rho_n * rho_c = FRICTIONAL HEATING multiplied by qdt
      tmp(ixO^S) = 0.5d0 * sum((vc(ixO^S, 1:ndir) - vn(ixO^S, 1:ndir))**2, dim=ndim+1)&
         * qdt * twofl_alpha_coll * rhoc(ixO^S) * rhon(ixO^S)
    endif
    if(.not. twofl_coll_inc_te) then
      if(phys_internal_e) then
        w(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp(ixO^S)
        w(ixO^S,e_c_) = w(ixO^S,e_c_) + tmp(ixO^S)
      else  
        ! tmp = qdt * 1/2 * (v_c^2-v_n^2) * alpha * rho_n * rho_c = FRICTIONAL HEATING + WORK done by coll terms in NEUTRALS mom
        ! eq. MULTIPLIED by qdt
        tmp(ixO^S) = 0.5d0*(sum(vc(ixO^S, 1:ndir)**2, dim=ndim+1) - sum(vn(ixO^S, 1:ndir)**2,dim=ndim+1))&
                     * qdt * twofl_alpha_coll * rhoc(ixO^S) * rhon(ixO^S)
        w(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp(ixO^S)
        w(ixO^S,e_c_) = w(ixO^S,e_c_) - tmp(ixO^S)
      endif  
    else
      ! coll term THERMAL EXCHANGE in neutrals energy eq. multiplied by qdti
      if(phys_total_energy) then
        tmp2(ixO^S) = qdt * twofl_alpha_coll * (-rhoc(ixO^S) * wCT(ixO^S,e_n_)*2d0/Rn +& 
                    rhon(ixO^S) * (wCT(ixO^S,e_c_) - twofl_mag_en(wCT,ixO^L,ixO^L)) *2d0/Rc )
      else
        tmp2(ixO^S) = qdt * twofl_alpha_coll * (-rhoc(ixO^S) * wCT(ixO^S,e_n_)*2d0/Rn +& 
                    rhon(ixO^S) * wCT(ixO^S,e_c_)*2d0/Rc )
      endif
      if(.not. phys_internal_e) then
          tmp2(ixO^S) = tmp2(ixO^S) + qdt * twofl_alpha_coll * (-rhoc(ixO^S) * twofl_kin_en_n(wCT,ixO^L,ixO^L)*(1d0 -2d0/Rn) +&
                                                            rhon(ixO^S) * twofl_kin_en_c(wCT,ixO^L,ixO^L)*(1d0 - 2d0/Rc))           
      endif
      if(phys_internal_e) then
        w(ixO^S,e_n_) = w(ixO^S,e_n_) + (tmp(ixO^S) + tmp2(ixO^S))
        w(ixO^S,e_c_) = w(ixO^S,e_c_) + (tmp(ixO^S) - tmp2(ixO^S))

      else
        w(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp2(ixO^S)
        w(ixO^S,e_c_) = w(ixO^S,e_c_) - tmp2(ixO^S)
      endif
    endif
    deallocate(vc,vn)


  end subroutine twofl_explicit_coll_terms_update
#endif
end module mod_twofl_phys
