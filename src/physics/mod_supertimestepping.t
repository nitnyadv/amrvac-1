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
  public :: add_sts_method
  public :: sts_add_source


  !> input parameter
  double precision :: sts_dtpar=0.9d0

  !> Whether to conserve fluxes at the current partial step
  logical :: fix_conserve_at_step = .true.

  logical :: sts_initialized = .false.

  abstract interface

  !this is used for setting sources
    subroutine subr1(ixI^L,ixO^L,w,x,wres)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
      double precision, intent(inout) :: wres(ixI^S,1:nw)

    end subroutine subr1

 
  !this is used for getting the timestep in dtnew
   function subr2(w,ixG^L,ix^L,dx^D,x) result(dtnew)
      use mod_global_parameters
      
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      double precision, intent(in) :: w(ixG^S,1:nw) 
      double precision :: dtnew
    end function subr2
  end interface

  type sts_term

    procedure (subr2), pointer, nopass :: sts_getdt
    procedure (subr1), pointer, nopass :: sts_set_sources
    type(sts_term), pointer :: next
    double precision :: dt_sts
    integer, public :: s

    !!types used for send/recv ghosts, see mod_ghostcells_update
    integer, dimension(-1:2^D&,-1:1^D&) :: type_send_srl_sts, type_recv_srl_sts
    integer, dimension(-1:1^D&,-1:1^D&) :: type_send_r_sts
    integer, dimension(-1:1^D&, 0:3^D&) :: type_recv_r_sts
    integer, dimension(-1:1^D&, 0:3^D&) :: type_recv_p_sts, type_send_p_sts

    integer :: startVar
    integer :: endVar

    integer, dimension(:), allocatable ::ixChange

  end type sts_term

  type(sts_term), pointer :: head_sts_terms

contains
  !> Read this module"s parameters from a file
  subroutine sts_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /sts_list/ sts_dtpar

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, sts_list, end=111)
111    close(unitpar)
    end do

  end subroutine sts_params_read



  subroutine add_sts_method(sts_getdt, sts_set_sources, startVar, endVar, ixChange)
    use mod_global_parameters
    use mod_ghostcells_update


    integer, intent(in) :: startVar, endVar 
    integer, intent(in) :: ixChange(:) 

    interface
      subroutine sts_set_sources(ixI^L,ixO^L,w,x,wres)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
      double precision, intent(inout) :: wres(ixI^S,1:nw)
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


    type_send_srl=>temp%type_send_srl_sts
    type_recv_srl=>temp%type_recv_srl_sts
    type_send_r=>temp%type_send_r_sts
    type_recv_r=>temp%type_recv_r_sts
    type_send_p=>temp%type_send_p_sts
    type_recv_p=>temp%type_recv_p_sts
    call create_bc_mpi_datatype(startVar,endVar-startVar+1)
    ! point bc mpi data type back to full type for (M)HD
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f

    temp%sts_getdt => sts_getdt
    temp%sts_set_sources => sts_set_sources
    

    temp%startVar = startVar
    temp%endVar = endVar
    allocate(temp%ixChange(size(ixChange)))
    temp%ixChange = ixChange

    temp%next => head_sts_terms
    head_sts_terms => temp

  end subroutine add_sts_method




  !> Initialize the module
  subroutine sts_init()
    use mod_global_parameters
    use mod_physics
    if(.not. sts_initialized) then  
      nullify(head_sts_terms)
      sts_dtpar=sts_dtpar/dble(ndim)
      call sts_params_read(par_files)
      sts_initialized = .true.
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


  function get_sts_ncycles(temp) result(s)
    use mod_global_parameters

    type(sts_term), pointer, intent(in)  :: temp
    integer :: s

    integer:: iigrid, igrid
    double precision :: dtnew,dtmin_mype
    double precision    :: dx^D


      dtmin_mype=bigdouble
    !$OMP PARALLEL DO PRIVATE(igrid,qdtnew,&
    !$OMP& dx^D) REDUCTION(min:dtmin_mype)
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          dx^D=rnode(rpdx^D_,igrid);
          dtmin_mype=min(dtmin_mype, sts_dtpar * temp%sts_getdt(ps(igrid)%w,ixG^LL,ixM^LL,dx^D,ps(igrid)%x))
       end do
    !$OMP END PARALLEL DO
       call MPI_ALLREDUCE(dtmin_mype,dtnew,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                             icomm,ierrmpi)
      ! get number of sub-steps of supertime stepping (Meyer 2012 MNRAS 422,2102)
       if(dt/dtnew< 0.5d0) then
         s=1
       else if(dt/dtnew< 2.d0) then
         s=2
       else
         s=ceiling((dsqrt(9.d0+8.d0*dt/dtnew)-1.d0)/2.d0)
         ! only use odd s number
         s=s/2*2+1
       endif
       if(mype==0) write(*,*) 'supertime steps:',s,' normal subcycles:',s


   end  function get_sts_ncycles


  subroutine sts_add_source()
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    
    double precision :: omega1,cmu,cmut,cnu,cnut
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid,j,i,ii,s
    logical :: evenstep, stagger_flag
    type(sts_term), pointer  :: temp
    type(state), dimension(:), pointer :: tmpPs1, tmpPs2

    ! not do fix conserve and getbc for staggered values if stagger is used
    stagger_flag=stagger_grid
    stagger_grid=.false.
    bcphys=.false.


    call init_comm_fix_conserve(1,ndim,1)

    !!tod pass this to sts_set_sources
    fix_conserve_at_step = time_advance .and. levmax>levmin

    temp => head_sts_terms
    do while(associated(temp))

  
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not. allocated(ps2(igrid)%w)) then
              allocate(ps2(igrid)%w(ixG^T,1:nw))
        endif  
        if(.not. allocated(ps3(igrid)%w)) then
              allocate(ps3(igrid)%w(ixG^T,1:nw))
        endif  
        if(.not. allocated(ps4(igrid)%w)) then
              allocate(ps4(igrid)%w(ixG^T,1:nw))
        endif  
        do i = 1,size(temp%ixChange)
          j=temp%ixChange(i)
          ps1(igrid)%w(ixG^T,j)=ps(igrid)%w(ixG^T,j)
          ps2(igrid)%w(ixG^T,j)=ps(igrid)%w(ixG^T,j)
          ps3(igrid)%w(ixG^T,j)=ps(igrid)%w(ixG^T,j)
        enddo
      enddo
      !$OMP END PARALLEL DO
  
      s=get_sts_ncycles(temp)

      allocate(bj(0:s))
      bj(0)=1.d0/3.d0
      bj(1)=bj(0)
      if(temp%s>1) then
        omega1=4.d0/dble(s**2+s-2)
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

      !!first step

      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        call temp%sts_set_sources(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,ps1(igrid)%w)
        do i = 1,size(temp%ixChange)
          ii=temp%ixChange(i)
          ps3(igrid)%w(ixM^T,ii) = dt * ps1(igrid)%w(ixM^T,ii)
          ps1(igrid)%w(ixM^T,ii) = ps1(igrid)%w(ixM^T,ii) + cmut * ps3(igrid)%w(ixM^T,ii)
        enddo

      end do
      !$OMP END PARALLEL DO
      ! fix conservation of AMR grid by replacing flux from finer neighbors
      !if (fix_conserve_at_step) then
      !  call recvflux(1,ndim)
      !  call sendflux(1,ndim)
      !  call fix_conserve(ps1,1,ndim,temp%startVar,temp%endVar-temp%startVar+1)
      !end if
      call getbc(global_time,0.d0,ps1,temp%startVar,temp%endVar-temp%startVar+1)
      !!first step end
      
      evenstep=.true.
      do j=2,s
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

      !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            call temp%sts_set_sources(ixG^LL,ixM^LL,ps1(igrid)%w, ps(igrid)%x,ps4(igrid)%w)
            do i = 1,size(temp%ixChange)
              ii=temp%ixChange(i)
              tmpPs2(igrid)%w(ixM^T,ii)=cmu*tmpPs1(igrid)%w(ixM^T,ii)+cnu*tmpPs2(igrid)%w(ixM^T,ii)+(1.d0-cmu-cnu)*ps(igrid)%w(ixM^T,ii)&
                        +cmut*dt*ps4(igrid)%w(ixM^T,ii)+cnut*ps3(igrid)%w(ixM^T,ii)
            end do
        end do
      !$OMP END PARALLEL DO
        call getbc(global_time,0.d0,tmpPs2,temp%startVar,temp%endVar-temp%startVar+1)
        evenstep=.not. evenstep
      end do
  
      if(evenstep) then
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          do i = 1,size(temp%ixChange)
            ii=temp%ixChange(i)
            ps(igrid)%w(ixG^T,ii)=ps1(igrid)%w(ixG^T,ii)
          end do 
        end do 
      else
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          do i = 1,size(temp%ixChange)
            ii=temp%ixChange(i)
            ps(igrid)%w(ixG^T,ii)=ps2(igrid)%w(ixG^T,ii)
          end do 
        end do 
      end if
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


    endif

  end subroutine sts_add_source


end module mod_supertimestepping
