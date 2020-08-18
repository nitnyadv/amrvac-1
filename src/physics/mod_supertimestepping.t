!> Generic supertimestepping method
!> adapted from mod_thermal_conduction
!>
!> 
!> PURPOSE: 
!> IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(KAPPA_i,j . GRAD_j T)
!> where KAPPA_i,j = tc_k_para b_i b_j + tc_k_perp (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(tc_k_para . GRAD T)
!> USAGE:
!> 1. in mod_usr.t -> subroutine usr_init(), add 
!>        unit_length=<your length unit>
!>        unit_numberdensity=<your number density unit>
!>        unit_velocity=<your velocity unit>
!>        unit_temperature=<your temperature unit>
!>    before call (m)hd_activate()
!> 2. to switch on thermal conduction in the (m)hd_list of amrvac.par add:
!>    (m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_perpendicular=.true.  ! (default .false.) turn on thermal conduction perpendicular to magnetic field 
!>    tc_saturate=.false.  ! (default .true. ) turn off thermal conduction saturate effect
!>    tc_dtpar=0.9/0.45/0.3 ! stable time step coefficient for 1D/2D/3D, decrease it for more stable run
!>    tc_slope_limiter='MC' ! choose limiter for slope-limited anisotropic thermal conduction in MHD

module mod_supertimestepping

  use mod_geometry
  implicit none

  private

  public :: is_sts_initialized
  public :: sts_init
  public :: add_sts_method,add_sts_method2
  public :: sts_add_source
  public :: set_dt_sts_ncycles
  public :: sourcetype_sts, sourcetype_sts_prior, sourcetype_sts_after, sourcetype_sts_split


  !> input parameters
  double precision :: sts_dtpar=0.9d0 !the coefficient that multiplies the sts dt
  integer :: sts_ncycles=1000 !the maximum number of subcycles

  !!method 2 only, not input
  double precision,parameter :: nu_sts = 0.5
  
  integer :: method_sts = 1
  integer, parameter :: sourcetype_sts_prior =0
  integer, parameter :: sourcetype_sts_after =1
  integer, parameter :: sourcetype_sts_split =2
  integer :: sourcetype_sts = sourcetype_sts_after

  !> Whether to conserve fluxes at the current partial step
  logical :: fix_conserve_at_step = .true.

  logical :: sts_initialized = .false.
  logical :: first = .true.

  abstract interface

  !this is used for setting sources
    subroutine subr1(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,indexChangeStart, indexChangeN, indexChangeFixC)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L,igrid
      double precision, intent(in) ::  x(ixI^S,1:ndim)
      double precision, intent(inout) :: wres(ixI^S,1:nw), w(ixI^S,1:nw)
      double precision, intent(in) :: my_dt
      logical, intent(in) :: fix_conserve_at_step
      integer, intent(in), dimension(:) :: indexChangeStart, indexChangeN
      logical, intent(in), dimension(:) :: indexChangeFixC

    end subroutine subr1

 
  !this is used for getting the timestep in dtnew
   function subr2(w,ixG^L,ix^L,dx^D,x) result(dtnew)
      use mod_global_parameters
      
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      double precision, intent(in) :: w(ixG^S,1:nw) 
      double precision :: dtnew
    end function subr2

  !!for the subroutines in this module:
  !sts_add_source
   subroutine subr3(dt)
    double precision,intent(in) :: dt
   end subroutine subr3
 

 !!sts_get_ncycles 
 function subr4(dt,dtnew,dt_modified) result(s)
    double precision,intent(in) :: dtnew
    double precision,intent(inout) :: dt
    logical,intent(inout) :: dt_modified
    integer :: s
 end function subr4

 !!for before_first_cycle and after_last_cycle
 subroutine subr5(ixI^L, ixO^L, w, x)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw) 
 end subroutine subr5

  end interface

  type sts_term

    procedure (subr2), pointer, nopass :: sts_getdt
    procedure (subr1), pointer, nopass :: sts_set_sources
    procedure (subr5), pointer, nopass :: sts_before_first_cycle, sts_after_last_cycle
    type(sts_term), pointer :: next
    double precision :: dt_expl
    integer, public :: s

    !!types used for send/recv ghosts, see mod_ghostcells_update
    integer, dimension(-1:2^D&,-1:1^D&) :: type_send_srl_sts, type_recv_srl_sts
    integer, dimension(-1:1^D&,-1:1^D&) :: type_send_r_sts
    integer, dimension(-1:1^D&, 0:3^D&) :: type_recv_r_sts
    integer, dimension(-1:1^D&, 0:3^D&) :: type_recv_p_sts, type_send_p_sts

    integer :: startVar
    integer :: endVar

    integer, dimension(:), allocatable ::ixChangeStart
    integer, dimension(:), allocatable ::ixChangeN
    logical, dimension(:), allocatable ::ixChangeFixC

  end type sts_term

  type(sts_term), pointer :: head_sts_terms
  !!The following two subroutine/function  pointers 
  !make the difference between the two STS methods implemented
  procedure (subr3), pointer :: sts_add_source
  procedure (subr4), pointer :: sts_get_ncycles

contains
  !> Read this module"s parameters from a file
  subroutine sts_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /sts_list/ sts_dtpar,sts_ncycles,method_sts,sourcetype_sts

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, sts_list, end=111)
111    close(unitpar)
    end do


  end subroutine sts_params_read

  !!TODO combine the following two subroutines
  subroutine add_sts_method2(sts_getdt, sts_set_sources, &
                             startVar, endVar, ixChangeStart, ixChangeN, ixChangeFixC)
    use mod_global_parameters
    use mod_ghostcells_update


    integer, intent(in) :: startVar, endVar 
    integer, intent(in) :: ixChangeStart(:) 
    integer, intent(in) :: ixChangeN(:) 
    logical, intent(in) :: ixChangeFixC(:) 

    interface
      subroutine sts_set_sources(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step, my_dt, igrid,indexChangeStart, indexChangeN, indexChangeFixC)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L,igrid
      double precision, intent(in) ::  x(ixI^S,1:ndim)
      double precision, intent(inout) :: wres(ixI^S,1:nw), w(ixI^S,1:nw)
      double precision, intent(in) :: my_dt
      logical, intent(in) :: fix_conserve_at_step
      integer, intent(in), dimension(:) :: indexChangeStart, indexChangeN
      logical, intent(in), dimension(:) :: indexChangeFixC

      end subroutine sts_set_sources
  
      function sts_getdt(w,ixG^L,ix^L,dx^D,x) result(dtnew)
        use mod_global_parameters
        
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      double precision, intent(in) :: w(ixG^S,1:nw) 
      double precision :: dtnew
      end function sts_getdt

    end interface

    type(sts_term), pointer  :: temp
    allocate(temp)
    !!types cannot be initialized here TODO see when init in mhd_phys is called
!    type_send_srl=>temp%type_send_srl_sts
!    type_recv_srl=>temp%type_recv_srl_sts
!    type_send_r=>temp%type_send_r_sts
!    type_recv_r=>temp%type_recv_r_sts
!    type_send_p=>temp%type_send_p_sts
!    type_recv_p=>temp%type_recv_p_sts
!    call create_bc_mpi_datatype(startVar,endVar-startVar+1)
!    ! point bc mpi data type back to full type for (M)HD
!    type_send_srl=>type_send_srl_f
!    type_recv_srl=>type_recv_srl_f
!    type_send_r=>type_send_r_f
!    type_recv_r=>type_recv_r_f
!    type_send_p=>type_send_p_f
!    type_recv_p=>type_recv_p_f

    temp%sts_getdt => sts_getdt
    temp%sts_set_sources => sts_set_sources
    temp%sts_before_first_cycle => null()
    temp%sts_after_last_cycle => null()

    temp%startVar = startVar
    temp%endVar = endVar
    if(size(ixChangeStart) .ne. size(ixChangeN) .or. size(ixChangeStart) .ne. size(ixChangeFixC) ) then
      if(mype .eq. 0) print*, "sizes are not equal ", size(ixChangeStart), size(ixChangeN), size(ixChangeFixC)
      return
    endif  
    allocate(temp%ixChangeStart(size(ixChangeStart)))
    allocate(temp%ixChangeN(size(ixChangeStart)))
    allocate(temp%ixChangeFixC(size(ixChangeStart)))
    temp%ixChangeStart = ixChangeStart
    temp%ixChangeN = ixChangeN
    temp%ixChangeFixC = ixChangeFixC

    temp%next => head_sts_terms
    head_sts_terms => temp

  end subroutine add_sts_method2

  !!the same as above without the sts_before_first_cycle and sts_after_last_cycle
  subroutine add_sts_method(sts_getdt, sts_set_sources, sts_before_first_cycle, sts_after_last_cycle,&
                             startVar, endVar, ixChangeStart, ixChangeN, ixChangeFixC)
    use mod_global_parameters
    use mod_ghostcells_update


    integer, intent(in) :: startVar, endVar 
    integer, intent(in) :: ixChangeStart(:) 
    integer, intent(in) :: ixChangeN(:) 
    logical, intent(in) :: ixChangeFixC(:) 

    interface
      subroutine sts_set_sources(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step, my_dt, igrid,indexChangeStart, indexChangeN, indexChangeFixC)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L,igrid
      double precision, intent(in) ::  x(ixI^S,1:ndim)
      double precision, intent(inout) :: wres(ixI^S,1:nw), w(ixI^S,1:nw)
      double precision, intent(in) :: my_dt
      logical, intent(in) :: fix_conserve_at_step
      integer, intent(in), dimension(:) :: indexChangeStart, indexChangeN
      logical, intent(in), dimension(:) :: indexChangeFixC

      end subroutine sts_set_sources
  
      function sts_getdt(w,ixG^L,ix^L,dx^D,x) result(dtnew)
        use mod_global_parameters
        
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      double precision, intent(in) :: w(ixG^S,1:nw) 
      double precision :: dtnew
      end function sts_getdt

    subroutine sts_before_first_cycle(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw) 
    end subroutine sts_before_first_cycle

    subroutine sts_after_last_cycle(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw) 
    end subroutine sts_after_last_cycle

    end interface

    type(sts_term), pointer  :: temp
    allocate(temp)
    !!types cannot be initialized here TODO see when init in mhd_phys is called
!    type_send_srl=>temp%type_send_srl_sts
!    type_recv_srl=>temp%type_recv_srl_sts
!    type_send_r=>temp%type_send_r_sts
!    type_recv_r=>temp%type_recv_r_sts
!    type_send_p=>temp%type_send_p_sts
!    type_recv_p=>temp%type_recv_p_sts
!    call create_bc_mpi_datatype(startVar,endVar-startVar+1)
!    ! point bc mpi data type back to full type for (M)HD
!    type_send_srl=>type_send_srl_f
!    type_recv_srl=>type_recv_srl_f
!    type_send_r=>type_send_r_f
!    type_recv_r=>type_recv_r_f
!    type_send_p=>type_send_p_f
!    type_recv_p=>type_recv_p_f

    temp%sts_getdt => sts_getdt
    temp%sts_set_sources => sts_set_sources
    temp%sts_before_first_cycle => sts_before_first_cycle
    temp%sts_after_last_cycle => sts_after_last_cycle

    temp%startVar = startVar
    temp%endVar = endVar
    if(size(ixChangeStart) .ne. size(ixChangeN) .or. size(ixChangeStart) .ne. size(ixChangeFixC) ) then
      if(mype .eq. 0) print*, "sizes are not equal ", size(ixChangeStart), size(ixChangeN), size(ixChangeFixC)
      return
    endif  
    allocate(temp%ixChangeStart(size(ixChangeStart)))
    allocate(temp%ixChangeN(size(ixChangeStart)))
    allocate(temp%ixChangeFixC(size(ixChangeStart)))
    temp%ixChangeStart = ixChangeStart
    temp%ixChangeN = ixChangeN
    temp%ixChangeFixC = ixChangeFixC

    temp%next => head_sts_terms
    head_sts_terms => temp

  end subroutine add_sts_method




  !> Initialize the module
  subroutine sts_init()
    use mod_global_parameters
    use mod_physics
    if(.not. sts_initialized) then  
      nullify(head_sts_terms)
      call sts_params_read(par_files)
      sts_dtpar=sts_dtpar/dble(ndim)
      sts_initialized = .true.
      if(method_sts .eq. 1) then
        sts_add_source => sts_add_source1
        sts_get_ncycles => sts_get_ncycles1
        if(mype .eq. 0) print*, "Method 1 STS"
      else if(method_sts .eq. 2) then
        sts_add_source => sts_add_source2
        sts_get_ncycles => sts_get_ncycles2
        if(mype .eq. 0) print*, "Method 2 STS"
      else
        call mpistop("Unknown sts method")
      endif
    endif
    
  end subroutine sts_init


  pure function is_sts_initialized() result(res)
  logical :: res
   if (sts_initialized) then
      res = associated(head_sts_terms)
    else
      res = .false.
    endif
  end function is_sts_initialized


  function sts_get_ncycles1(dt,dtnew,dt_modified) result (is)
    double precision,intent(in) :: dtnew
    double precision,intent(inout) :: dt
    logical,intent(inout) :: dt_modified
    integer :: is

    double precision    :: ss
    integer:: ncycles

    !!ss is now limit of dt because of sts_ncycles
    ss = dtnew * ((2d0 * sts_ncycles + 1)**2 -9.0)/16d0
    if(dt>ss) then
       dt_modified = .true.
       dt = ss
       is = sts_ncycles 
     else  
       ss = dt/dtnew
       ! get number of sub-steps of supertime stepping (Meyer 2012 MNRAS 422,2102)
       if(ss .le. 1d0) then
         is=1
       else
         is=ceiling((dsqrt(9.d0+16.d0*ss)-1.d0)/2.d0)
         is=is/2*2+1
      endif
    endif
    !print*, dt, " --DTEXPL-- ", dtnew, ", ncycle1 ",is
  end  function sts_get_ncycles1


  function sts_get_ncycles2(dt,dtnew,dt_modified) result(is)
    double precision,intent(in) :: dtnew
    double precision,intent(inout) :: dt
    logical,intent(inout) :: dt_modified
    integer :: is

    double precision    :: ss,rr
    integer:: ncycles
      rr = dt/dtnew
      !print*, dt, " --DTEXPL-- ", dtnew, ", rr ",rr
      ncycles = sts_ncycles 
      !print*, "NCYCLES BEFORE ",ncycles
      ss=sum_chev(nu_sts,ncycles,rr)
      !print*, "NCYCLES AFTER ",ncycles
      is = ncycles
      !print*, "SUMCHEV ", ss, " NCYCLES ", is
      if(ss < rr) then
        dt_modified = .true.
        dt = ss *  dtnew
      endif

  end  function sts_get_ncycles2


  function set_dt_sts_ncycles(my_dt) result(dt_modified)
    use mod_global_parameters

    double precision,intent(inout) :: my_dt
    logical :: dt_modified

    integer:: iigrid, igrid, ncycles
    double precision :: dtnew,dtmin_mype
    double precision    :: dx^D, ss
    type(sts_term), pointer  :: temp,oldTemp
    nullify(oldTemp)
    temp => head_sts_terms
    dt_modified = .false.
    do while(associated(temp))
      dt_modified = .false.
      dtmin_mype=bigdouble
      !$OMP PARALLEL DO PRIVATE(igrid,dtnew,&
      !$OMP& dx^D) REDUCTION(min:dtmin_mype)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          !!maybe the following global variables are neeeded in get_dt!!!!!!!
          ! next few lines ensure correct usage of routines like divvector etc
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            typelimiter=type_limiter(node(plevel_,igrid))
            typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
            !!end maybe the following global variables are neeeded in get_dt!!!!!!!
            dx^D=rnode(rpdx^D_,igrid);
            dtmin_mype=min(dtmin_mype, sts_dtpar * temp%sts_getdt(ps(igrid)%w,ixG^LL,ixM^LL,dx^D,ps(igrid)%x))
      end do
      !$OMP END PARALLEL DO
      call MPI_ALLREDUCE(dtmin_mype,dtnew,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                               icomm,ierrmpi)
      temp%s = sts_get_ncycles(my_dt,dtnew,dt_modified)  
      temp%dt_expl = dtnew
 
      !if(mype==0) write(*,*) 'supertime steps:',temp%s, ", dt is ", dt,&
      !   " MODI ",dt_modified, " dt expl ", dtnew, ", my_dt=",my_dt
       if(dt_modified) then
         temp => head_sts_terms
         oldTemp=>temp
       else
         temp=>temp%next
       endif
       if(associated(oldTemp,temp)) then
        temp=>temp%next
       endif

    enddo

    !if(mype==0) write(*,*) 'EXIT'
  end function set_dt_sts_ncycles


!!!> IMPLEMENTATION2
   pure FUNCTION chev(j,nu,N)
    use mod_constants

    double precision, INTENT(IN) :: nu
    INTEGER, INTENT(IN)       :: j, N
    double precision             :: chev

    chev = 1d0 / ((-1d0 + nu)*cos(((2d0*j - 1d0) / N)* (dpi/2d0)) + 1d0 + nu)

  END FUNCTION chev

  FUNCTION sum_chev(nu,N,limMax)
    double precision, intent(in) :: nu,limmax
    integer, intent(inout)       ::  N
    double precision             :: sum_chev, tmp

    integer :: j

    j=1
    sum_chev = 0d0
    !print*, N," LIMMAX ", limMax
    do while (j < N .and. sum_chev < limMax)
      sum_chev = sum_chev + chev(j,nu,N)
      j=j+1
      !print*, " SS j ", limMax,sum_chev,j
    enddo
   N=j-1 
   !print*, " RES ", N,sum_chev
  END FUNCTION sum_chev



  PURE FUNCTION total_chev(nu,N)
    double precision, INTENT(IN) :: nu
    INTEGER, INTENT(IN)       :: N
    double precision             :: total_chev

    total_chev = N/(2d0*dsqrt(nu)) * ( (1d0 + dsqrt(nu))**(2d0*N) - (1d0 - dsqrt(nu))**(2d0*N) ) / &
         ( (1d0 + dsqrt(nu))**(2d0*N) + (1d0 - dsqrt(nu))**(2d0*N) )

  END FUNCTION total_chev
!!!> IMPLEMENTATION2 end




  subroutine sts_add_source2(my_dt)
  ! Turlough Downes 2006,2007
    use mod_ghostcells_update
    use mod_global_parameters
    use mod_fix_conserve
    
    double precision, intent(in) :: my_dt
    double precision, allocatable :: bj(:)
    double precision :: sumbj,dtj  

    integer:: iigrid, igrid,j,i,ii,ii2
    logical :: stagger_flag, prolong_flag, coarsen_flag
    type(sts_term), pointer  :: temp

    ! not do fix conserve and getbc for staggered values if stagger is used
    stagger_flag=stagger_grid
    prolong_flag = prolongprimitive
    coarsen_flag = coarsenprimitive
    stagger_grid=.false.
    prolongprimitive = .false.
    coarsenprimitive = .false.
    bcphys=.false.


    call init_comm_fix_conserve(1,ndim,1)

    fix_conserve_at_step = time_advance .and. levmax>levmin

    temp => head_sts_terms
    do while(associated(temp))

      if(associated(temp%sts_before_first_cycle)) then
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call temp%sts_before_first_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)  
        enddo 
      endif
      


      allocate(bj(1:temp%s))
      do j=1,temp%s
        bj(j) = chev(j,nu_sts,sts_ncycles)
      enddo

      type_send_srl=>temp%type_send_srl_sts
      type_recv_srl=>temp%type_recv_srl_sts
      type_send_r=>temp%type_send_r_sts
      type_recv_r=>temp%type_recv_r_sts
      type_send_p=>temp%type_send_p_sts
      type_recv_p=>temp%type_recv_p_sts
        
       if(first) then 
        call create_bc_mpi_datatype(temp%startVar,temp%endVar-temp%startVar+1)
        first = .false.
       endif 
      !!first step
      !call getbc(global_time,0.d0,ps,temp%startVar,temp%endVar-temp%startVar+1)

      sumbj=0d0
      do j=1,temp%s
        if(j .eq. temp%s .and. (sumbj + bj(j)) * temp%dt_expl > my_dt) then
          dtj = my_dt - sumbj * temp%dt_expl
        else
          dtj = bj(j)* temp%dt_expl
        endif  
        sumbj = sumbj + bj(j)
        !print*, "Substep ",j, ", sumbj=", sumbj, ", dtj=", dtj
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);

          !!maybe the following global variables are neeeded in set_sources!!!!!!!
          ! next few lines ensure correct usage of routines like divvector etc
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            typelimiter=type_limiter(node(plevel_,igrid))
            typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
            !!end maybe the following global variables are neeeded in set_sources!!!!!!!

            !print*, "ID_sts ",igrid
            call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w, ps(igrid)%x,ps1(igrid)%w, fix_conserve_at_step, &
                                    dtj,igrid,temp%ixChangeStart, temp%ixChangeN, temp%ixChangeFixC)
            do i = 1,size(temp%ixChangeStart)
              ii=temp%ixChangeStart(i)
              ii2 = ii + temp%ixChangeN(i) -1
              ps(igrid)%w(ixM^T,ii:ii2)=ps(igrid)%w(ixM^T,ii:ii2)+dtj*ps1(igrid)%w(ixM^T,ii:ii2)
            end do
            !if( igrid .eq. 1) print*, " ps1Bx " , ps1(igrid)%w(1:10,8)
            !if( igrid .eq. 1) print*, " psBx " , ps(igrid)%w(1:10,8)
            ! if( igrid .eq. 16) print*," psBxEND " ,  ps(igrid)%w(size(ps(igrid)%w,1)-9:size(ps(igrid)%w,1),8)

        end do
      !$OMP END PARALLEL DO
        if (fix_conserve_at_step) call fix_conserve_vars(ps, temp%ixChangeStart, temp%ixChangeN, temp%ixChangeFixC)
        call getbc(global_time,0.d0,ps,temp%startVar,temp%endVar-temp%startVar+1)
      end do

      if(associated(temp%sts_after_last_cycle)) then 
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call temp%sts_after_last_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do 
      endif
      deallocate(bj)

      temp=>temp%next
    enddo
    if(associated(head_sts_terms)) then
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
      bcphys=.true.
  
      ! restore stagger_grid value
      stagger_grid=stagger_flag
      prolongprimitive  = prolong_flag
      coarsenprimitive = coarsen_flag


    endif

  end subroutine sts_add_source2

  subroutine sts_add_source1(my_dt)
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    
    double precision, intent(in) :: my_dt
    double precision :: dtj
    double precision :: omega1,cmu,cmut,cnu,cnut
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid,j,i,ii,ii2
    logical :: evenstep, stagger_flag, prolong_flag, coarsen_flag
    type(sts_term), pointer  :: temp
    type(state), dimension(:), pointer :: tmpPs1, tmpPs2

    ! not do fix conserve and getbc for staggered values if stagger is used
    stagger_flag=stagger_grid
    prolong_flag = prolongprimitive
    coarsen_flag = coarsenprimitive
    stagger_grid=.false.
    prolongprimitive = .false.
    coarsenprimitive = .false.
    bcphys=.false.


    call init_comm_fix_conserve(1,ndim,1)

    fix_conserve_at_step = time_advance .and. levmax>levmin

    temp => head_sts_terms
    do while(associated(temp))

      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(associated(temp%sts_before_first_cycle)) then
          call temp%sts_before_first_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)  
        endif 
        if(.not. allocated(ps2(igrid)%w)) then
              allocate(ps2(igrid)%w(ixG^T,1:nw))
        endif  
        if(.not. allocated(ps3(igrid)%w)) then
              allocate(ps3(igrid)%w(ixG^T,1:nw))
        endif  
        if(.not. allocated(ps4(igrid)%w)) then
              allocate(ps4(igrid)%w(ixG^T,1:nw))
        endif 
          !we need the full set of variables
          !for the calculation of the sources and not
          !only those which change
        ps1(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
        ps2(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
        ps3(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
        ps4(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
      enddo
      !$OMP END PARALLEL DO
  

      allocate(bj(0:temp%s))
      bj(0)=1.d0/3.d0
      bj(1)=bj(0)
      if(temp%s>1) then
        omega1=4.d0/dble(temp%s**2+temp%s-2)
        cmut=omega1/3.d0
      else
        omega1=0.d0
        cmut=1.d0
      endif

      type_send_srl=>temp%type_send_srl_sts
      type_recv_srl=>temp%type_recv_srl_sts
      type_send_r=>temp%type_send_r_sts
      type_recv_r=>temp%type_recv_r_sts
      type_send_p=>temp%type_send_p_sts
      type_recv_p=>temp%type_recv_p_sts
        
       if(first) then 
        call create_bc_mpi_datatype(temp%startVar,temp%endVar-temp%startVar+1)
        first = .false.
       endif 
      !!first step
      !call getbc(global_time,0.d0,ps,temp%startVar,temp%endVar-temp%startVar+1)
      dtj = cmut*my_dt
      !print*, "STS1 1 dtj = ", dtj
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          !!maybe the following global variables are neeeded in set_sources!!!!!!!
          ! next few lines ensure correct usage of routines like divvector etc
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            typelimiter=type_limiter(node(plevel_,igrid))
            typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
            !!end maybe the following global variables are neeeded in set_sources!!!!!!!
!        print*, "1sID_sts ",igrid
        call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,ps4(igrid)%w, fix_conserve_at_step, dtj, &
                                    igrid, temp%ixChangeStart, temp%ixChangeN, temp%ixChangeFixC)

!        print*, "FIRST STep Bx ", igrid, minval(ps4(igrid)%w(ixM^T,8)), maxval(ps4(igrid)%w(ixM^T,8))
!        print*, "FIRST STep e ", igrid, minval(ps4(igrid)%w(ixM^T,5)), maxval(ps4(igrid)%w(ixM^T,5))
!        print*, "BEFORE ADD FIRST STep minmax wBx ", igrid, &
!                minval(ps1(igrid)%w(ixM^T,8)), maxval(ps1(igrid)%w(ixM^T,8))
!        print*, "BEFORE ADD FIRST STep minmax wE ", igrid, &
!                minval(ps1(igrid)%w(ixM^T,5)), maxval(ps1(igrid)%w(ixM^T,5))


        !!!eq solved: dU/dt = S
        !!!In ps3 is stored S^n
        do i = 1,size(temp%ixChangeStart)
          ii=temp%ixChangeStart(i)
          ii2 = ii + temp%ixChangeN(i) -1
          ps3(igrid)%w(ixM^T,ii:ii2) = my_dt * ps4(igrid)%w(ixM^T,ii:ii2)
          ps1(igrid)%w(ixM^T,ii:ii2) = ps1(igrid)%w(ixM^T,ii:ii2) + cmut * ps3(igrid)%w(ixM^T,ii:ii2)
        enddo
!        if( igrid .eq. 1) print*, " 1sperpps1 bx " , ps1(igrid)%w(1:10,8)
!        if( igrid .eq. 1)   print*," 1stepps1 e " , ps1(igrid)%w(1:10,5)
!        if( igrid .eq. 1) print*, " 1stepps4Bx " , ps4(igrid)%w(1:10,8)

!        print*, "AFTER ADD FIRST STep minmax wBx ", igrid, &
!                minval(ps1(igrid)%w(ixM^T,8)), maxval(ps1(igrid)%w(ixM^T,8))
!        print*, "AFTER ADD FIRST STep minmax wBx  first 10 ", igrid, &
!                ps1(igrid)%w(1:10,8)
!        print*, "AFTER ADD FIRST STep minmax wE ", igrid, &
!                minval(ps1(igrid)%w(ixM^T,5)), maxval(ps1(igrid)%w(ixM^T,5))
!        print*, "AFTER ADD FIRST STep minmax wE  first 10 ", igrid, &
!                ps1(igrid)%w(1:10,5)
        !ps3(igrid)%w(ixM^T,1:nw) = dt * ps1(igrid)%w(ixM^T,1:nw)
        !ps1(igrid)%w(ixM^T,1:nw) = ps1(igrid)%w(ixM^T,1:nw) + cmut * ps3(igrid)%w(ixM^T,1:nw)

      end do
      !$OMP END PARALLEL DO
      ! fix conservation of AMR grid by replacing flux from finer neighbors
      !if (fix_conserve_at_step) then
      !  call recvflux(1,ndim)
      !  call sendflux(1,ndim)
      !  call fix_conserve(ps1,1,ndim,temp%startVar,temp%endVar-temp%startVar+1)
      !end if
      if (fix_conserve_at_step) call fix_conserve_vars(ps1, temp%ixChangeStart, temp%ixChangeN, temp%ixChangeFixC)
      call getbc(global_time,0.d0,ps1,temp%startVar,temp%endVar-temp%startVar+1)
      !!first step end
      
      evenstep=.true.

!      print*, "ps1 ",loc(ps1)
!      print*, "ps2 ",loc(ps2)
      tmpPs2=>ps1

      do j=2,temp%s
        bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
        cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
        cmut=omega1*cmu
        cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
        cnut=(bj(j-1)-1.d0)*cmut
        if(evenstep) then
          tmpPs1=>ps1
          tmpPs2=>ps2
        else
          tmpPs1=>ps2
          tmpPs2=>ps1
        endif

      !print*, "-------------------------------------------------------------------"
      !print*,"IGRIDSTAIL ACTIVE ", igridstail_active
      dtj = cmut*my_dt
      !print*, "STS1 ", j, " dtj = ", dtj
      !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            !print*, "ID_sts ",igrid
          !!maybe the following global variables are neeeded in set_sources!!!!!!!
          ! next few lines ensure correct usage of routines like divvector etc
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            block=>ps(igrid)
            typelimiter=type_limiter(node(plevel_,igrid))
            typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
            !!end maybe the following global variables are neeeded in set_sources!!!!!!!
            call temp%sts_set_sources(ixG^LL,ixM^LL,tmpPs1(igrid)%w, ps(igrid)%x,ps4(igrid)%w, fix_conserve_at_step,dtj, &
                                    igrid, temp%ixChangeStart, temp%ixChangeN, temp%ixChangeFixC)
            !print*, "** tmpPs2 ",igrid, loc(tmpPs2(igrid)), " tmpPs1 ", loc(tmpPs1(igrid))
            !print*, " STep ",j , " e ",evenstep, " tmpPs1 ", minval(tmpPs1(igrid)%w),maxval(tmpPs1(igrid)%w)
!            print*, " BEFORE ADD STep ",j ,  " tmpPs2 ", minval(tmpPs1(igrid)%w),maxval(tmpPs1(igrid)%w)
!            print*, "LOCTmpps1,2 ", loc(tmpPs1), loc(tmpPs2)
!            print*, " tmpPs1bx " , tmpPs1(igrid)%w(1:10,8)
!            print*, " tmpPs1e " , tmpPs1(igrid)%w(1:10,5)
            !print*, " STep ",j ,  " ps3 ", minval(ps3(igrid)%w),maxval(ps3(igrid)%w)
            !print*, " STep ",j ,  " ps4 ", minval(ps4(igrid)%w),maxval(ps4(igrid)%w)
            do i = 1,size(temp%ixChangeStart)
              ii=temp%ixChangeStart(i)
              ii2 = ii + temp%ixChangeN(i) - 1
              tmpPs2(igrid)%w(ixM^T,ii:ii2)=cmu*tmpPs1(igrid)%w(ixM^T,ii:ii2)+cnu*tmpPs2(igrid)%w(ixM^T,ii:ii2)+(1.d0-cmu-cnu)*ps(igrid)%w(ixM^T,ii:ii2)&
                        +dtj*ps4(igrid)%w(ixM^T,ii:ii2)+cnut*ps3(igrid)%w(ixM^T,ii:ii2)
            end do
            !tmpPs2(igrid)%w(ixM^T,1:nw)=cmu*tmpPs1(igrid)%w(ixM^T,1:nw)+cnu*tmpPs2(igrid)%w(ixM^T,1:nw)+(1.d0-cmu-cnu)*ps(igrid)%w(ixM^T,1:nw)&
            !            +cmut*dt*ps4(igrid)%w(ixM^T,1:nw)+cnut*ps3(igrid)%w(ixM^T,1:nw)
            !print*, " STep ",j , " after_set ", minval(tmpPs2(igrid)%w),maxval(tmpPs2(igrid)%w)
!            if( igrid .eq. 1) print*, " ps4dtBx " , dt*ps4(igrid)%w(1:10,8)
!            if( igrid .eq. 1) print*, " ps3Bx " , dt*ps3(igrid)%w(1:10,8)
!            if( igrid .eq. 1) print*, " cmn " , cmu,cnu,cmut,cnut,cmut*dt,cnut*dt
!            if( igrid .eq. 1) print*, " psBx " , ps(igrid)%w(1:10,8)
!
!            if( igrid .eq. 1) print*, " tmpPs1Bx " , tmpPs1(igrid)%w(1:10,8)
!            if( igrid .eq. 1) print*, " tmpPs2Bx " , tmpPs2(igrid)%w(1:10,8)
!            if( igrid .eq. 1)   print*," tmpPs2e " , tmpPs2(igrid)%w(1:10,5)
        end do
      !$OMP END PARALLEL DO
        if (fix_conserve_at_step) call fix_conserve_vars(tmpPs2, temp%ixChangeStart, temp%ixChangeN, temp%ixChangeFixC)
        call getbc(global_time,0.d0,tmpPs2,temp%startVar,temp%endVar-temp%startVar+1)
        evenstep=.not. evenstep
!        print*, "tmpPs2INLOOP ",loc(tmpPs2)
      end do

!      print*, "IXCHANGE ", temp%ixChange
!      print*, "SIZE IXCHANGE ", size(temp%ixChange)
!      print*, "tmpPs2 ",loc(tmpPs2)
      !print*, "IGRIDSTAIL  ", igridstail
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        do i = 1,size(temp%ixChangeStart)
          ii=temp%ixChangeStart(i)
          ii2 = ii + temp%ixChangeN(i) - 1
          ps(igrid)%w(ixG^T,ii:ii2)=tmpPs2(igrid)%w(ixG^T,ii:ii2)
        end do 
        !ps(igrid)%w(ixG^T,1:nw)=tmpPs2(igrid)%w(ixG^T,1:nw)
!        if(igrid .eq. 1) then
!          print*, igrid, " MINMAX_EnsSTS Bx ", minval(ps(igrid)%w(ixG^T,8)), maxval(ps(igrid)%w(ixG^T,8))
!          print*, igrid,  " MINMAX_EnsSTS Bx first 10 ", ps(igrid)%w(1:10,8)
!          print*, igrid,  " MINMAX_EnsSTS E ", minval(ps(igrid)%w(ixG^T,5)), maxval(ps(igrid)%w(ixG^T,5))
!          print*, igrid,  " MINMAX_EnsSTS E firat 10 ", ps(igrid)%w(1:10,5)
!        endif
        !print*, " MINMAXloc  ", minloc(ps(igrid)%w(ixG^T,1:nw)), maxloc(ps(igrid)%w(ixG^T,1:nw))
      end do 
      !$OMP END PARALLEL DO

      deallocate(bj)
      if(associated(temp%sts_after_last_cycle)) then 
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call temp%sts_after_last_cycle(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do 
      endif

      temp=>temp%next
    enddo
    if(associated(head_sts_terms)) then
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
      bcphys=.true.
  
      ! restore stagger_grid value
      stagger_grid=stagger_flag
      prolongprimitive  = prolong_flag
      coarsenprimitive = coarsen_flag


    endif

  end subroutine sts_add_source1


end module mod_supertimestepping
