#include "global.h"

module scdm_m
! module list is copied from hamiltonian
  use batch_m
  use batch_ops_m
  use blas_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use cmplxscl_m
  use cube_m
  use cube_function_m
  use datasets_m
  use derivatives_m
  use energy_m
!  use hamiltonian_base_m
!  use epot_m
!  use gauge_field_m
  use fft_m
  use nfft_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use index_m
  use io_m
  use io_function_m
  use kpoints_m
  use lalg_basic_m
!  use lasers_m
  use math_m
  use mesh_m
  use mesh_cube_map_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use opencl_m
  use ob_interface_m
!  use ob_lead_m
  use opencl_m
  use par_vec_m
  use parser_m
  use poisson_m
  use poisson_fft_m
  use profiling_m
!  use projector_m
  use simul_box_m
  use smear_m
  use states_m
  use states_calc_m
  use states_dim_m
  use types_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)

  implicit none

  private
  public ::               &
       scdm_t,            &
       scdm_init,         &
       dscdm_localize,    &
       zscdm_localize

  type scdm_t
     type(states_t)   :: st          ! localized orthogonal states
     type(poisson_t)  :: poisson1    ! solver (not used, only for testing)
     type(cube_t)     :: cube        ! mesh cube for fft
     FLOAT, pointer   :: center(:,:) ! coordinates of centers of states (in same units as mesh%x)
     FLOAT            :: rcut        ! orbital cutoff radius (box size) NOTE: this could be dynamic and state dependent
     integer          :: box_size    ! number of mesh points in the dimension of local box around scdm states NOTE: this could be dynamic and state dependent
     integer          :: full_box    ! = (2*(2*box_size+1))**3, i.e. number of points in box
     type(mesh_t)     :: boxmesh     ! mesh describing the small box
     type(cube_t)     :: boxcube     ! cube of the small box (used for fft in poisson solver, has doubled size for truncation)
     integer, pointer :: box(:,:,:,:)  ! indices of global points that are contained in the local box for each state
     FLOAT, pointer   :: dpsi(:,:)   ! scdm states in their local box
     CMPLX, pointer   :: zpsi(:,:)   ! ^
     type(poisson_t)  :: poisson     ! solver used to compute exchange with localized scdm states
     type(poisson_fft_t) :: poisson_fft ! used for above poisson solver
     type(cmplxscl_t)    :: cmplxscl
     logical, pointer :: periodic(:) ! tracks whcih scdm states are split by the periodic boundary conditions
     !
     logical          :: re_ortho_normalize=.false. ! orthonormalize the scdm states
     logical          :: verbose     ! write infor about SCDM procedure
     !
     ! parallelization of scdm states
     logical          :: root        ! this is a redundat flag equal to mesh%vp%rank==0
     integer          :: nst
     integer          :: st_start    ! the distributed index
     integer          :: st_end      ! .
     integer          :: lnst        ! .
     !
  end type scdm_t
  !
  logical,public    :: scdm_is_init=.false.  ! is initialized
  logical,public    :: scdm_is_local=.false.  ! is localized
  !
  type(scdm_t), public       :: scdm
! debug stuff
  type(geometry_t), public   :: scdm_geo
  !
  contains

  subroutine scdm_init(st,der,scdm)
    !
    type(states_t), intent(in)  :: st ! this contains the KS set (for now from hm%hf_st which is confusing)
    type(derivatives_t) :: der
    type(scdm_t) :: scdm
    type(cmplxscl_t) :: cmplxscl
    !
    integer :: i,j,k, ip, rank
    !debug
    integer :: temp(3)
    !
    integer,  allocatable:: istart(:)
    integer,  allocatable:: iend(:)
    integer,  allocatable:: ilsize(:)
    integer :: box(3)
    FLOAT :: dummy, enlarge
    !
    ! check if already initialized
    if(scdm_is_init) return
    !
    if(st%lnst.ne.st%nst) call messages_not_implemented("SCDM with state parallelization")
    if(st%d%nik.gt.1) call messages_not_implemented("SCDM with k-point sampling")
    if(der%mesh%sb%periodic_dim.gt.0.and.der%mesh%sb%periodic_dim.ne.3) &
                   call messages_not_implemented("SCDM with mixed-periodicity")  
    !
!    scdm%root = (der%mesh%vp%rank == 0)
    call MPI_Comm_Rank( der%mesh%mpi_grp%comm, rank, mpi_err)
    scdm%root = (rank ==0)
    !
    ! inherit some indices from st
    !scdm%st%d%dim = st%d%dim
    scdm%st%nst   = st%nst
    scdm%nst   = st%nst
    scdm%st%d%nik = st%d%nik
    scdm%st%d     = st%d
    scdm%cmplxscl = st%cmplxscl
    call parse_logical(datasets_check('SCDM_reorthonormalize'), .false., scdm%re_ortho_normalize)
    if(scdm%re_ortho_normalize) scdm%st%d%orth_method = ORTH_CHOLESKY_SERIAL
    !
    call parse_logical(datasets_check('SCDM_verbose'), .false., scdm%verbose)
    !
    ! allocate centers
    SAFE_ALLOCATE(scdm%center(3,scdm%st%nst))
    !
    ! make a cube around the center points
    ! with side length NOTE: this should be dynamic
    call parse_float(datasets_check('SCDMCutoffRadius'), 3._8, scdm%rcut, units_inp%length)
    if(scdm%root.and.scdm%verbose) call messages_print_var_value(stdout, 'SCDM cutoff', scdm%rcut)
    ! box_size is half the size of the  box
    scdm%box_size = 0
    do i=1,3
       scdm%box_size = max(scdm%box_size,ceiling(scdm%rcut/der%mesh%spacing(i)))
    enddo
    !
    if(scdm%root.and.scdm%verbose) then
       call messages_print_var_value(stdout,'SCDM box_size', scdm%box_size)
       call messages_print_var_value(stdout,'SCDM box_size[Ang]', scdm%box_size*der%mesh%spacing(1)*0.529177249)
    endif
!    scdm%full_box = (2*(2*scdm%box_size+1))**3
scdm%full_box = (2*scdm%box_size+1)**3
    !check if scdm is not bigger than fft-grid of full simualtion cell  
    if(scdm%full_box.gt.der%mesh%np_global) then
       message(1) = 'SCDM box larger than mesh, no point in using it'
       call messages_fatal(1,only_root_writes = .true.)
    endif
    dummy = 2*(2*scdm%box_size+1)*der%mesh%spacing(1)*0.529177249
    if(scdm%root.and.scdm%verbose) call messages_print_var_value(stdout, 'SCDM fullbox[Ang]', dummy)
    SAFE_ALLOCATE(scdm%box(scdm%box_size*2+1,scdm%box_size*2+1,scdm%box_size*2+1,scdm%st%nst))
    !
    ! the localzied states defined in the box are distributed over state index
    SAFE_ALLOCATE(istart(der%mesh%mpi_grp%size))
    SAFE_ALLOCATE(iend(der%mesh%mpi_grp%size))
    SAFE_ALLOCATE(ilsize(der%mesh%mpi_grp%size))
    !
    call multicomm_divide_range(st%nst,der%mesh%mpi_grp%size, istart, iend, lsize=ilsize)
    scdm%st_start = istart(der%mesh%vp%rank+1)
    scdm%st_end = iend(der%mesh%vp%rank+1)
    scdm%lnst = ilsize(der%mesh%vp%rank+1)
    !
    ! allocate local chunk of states
    ! localized SCDM states in full box (not really needed, but convenient for distribution)
    ! root process holds all states NOTE: this is not great... clearly ... but will go away when SCDM procedure is parallel
    if(.not.states_are_real(st)) then
       if(scdm%root) then
          SAFE_ALLOCATE(scdm%st%zpsi(der%mesh%np_global, scdm%st%d%dim, scdm%st%nst, scdm%st%d%nik))
       else
          SAFE_ALLOCATE(scdm%st%zpsi(der%mesh%np_global, scdm%st%d%dim, scdm%lnst, scdm%st%d%nik))
       endif
       ! localized SCDM states defined on box twice their size for coulomb truncation
       SAFE_ALLOCATE(scdm%zpsi(1:scdm%full_box,scdm%lnst))
    else ! real
       if(scdm%root) then
          SAFE_ALLOCATE(scdm%st%dpsi(der%mesh%np_global, scdm%st%d%dim, scdm%st%nst, scdm%st%d%nik))
       else
          SAFE_ALLOCATE(scdm%st%dpsi(der%mesh%np_global, scdm%st%d%dim, scdm%lnst, scdm%st%d%nik))
       endif
       ! localized SCDM states defined on box twice their size for coulomb truncation
       SAFE_ALLOCATE(scdm%dpsi(1:scdm%full_box,scdm%lnst))
    endif
    !
    SAFE_ALLOCATE(scdm%periodic(scdm%lnst))
    !
    ! create a mesh object for the small box (for now each scdm state is in the same box, should be dynamic)
    ! only initialize values needed in the following (e.g. by poisson_fft_init)
    scdm%boxmesh%spacing(:) = minval(der%mesh%spacing(:))
    SAFE_ALLOCATE(scdm%boxmesh%sb)
    scdm%boxmesh%sb%periodic_dim = 0
    scdm%boxmesh%sb%dim = 3
    scdm%boxmesh%sb%klattice_primitive(:,:) = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
    scdm%boxmesh%sb%rlattice_primitive(:,:) = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
    !
    !set mesh points with double size for coulomb truncation
    scdm%boxmesh%np = scdm%full_box
    scdm%boxmesh%np_global = scdm%boxmesh%np
    scdm%boxmesh%np_part = scdm%boxmesh%np_global
    scdm%boxmesh%np_part_global = scdm%boxmesh%np_global
    ! set index type of mesh
    scdm%boxmesh%idx%is_hypercube = .false.
    ! mesh has to be centered around zero with left overhang otherwise mesh_cub_map doesn't seem to work
!    scdm%boxmesh%idx%nr(1,:) = -(scdm%box_size*2+1)
!    scdm%boxmesh%idx%nr(2,:) =  (scdm%box_size*2+1)-1
scdm%boxmesh%idx%nr(1,:) = -(scdm%box_size)
scdm%boxmesh%idx%nr(2,:) =  (scdm%box_size) 
    !
    scdm%boxmesh%idx%dim = 3
    scdm%boxmesh%idx%ll(:) = scdm%boxmesh%idx%nr(2,:) - scdm%boxmesh%idx%nr(1,:) + 1
!???
    scdm%boxmesh%idx%enlarge(:) = 0
    SAFE_ALLOCATE(scdm%boxmesh%idx%lxyz(scdm%boxmesh%np,scdm%boxmesh%idx%dim))
    ! need to copy indices because otherwise line gets too long (precompiler???)
!    i=-(scdm%box_size*2+1)
!    j=(scdm%box_size*2+1)-1
i=-(scdm%box_size)
j=(scdm%box_size)
    SAFE_ALLOCATE(scdm%boxmesh%idx%lxyz_inv(i:j,i:j,i:j))
    !
    ip = 0
    do i=scdm%boxmesh%idx%nr(1,1),scdm%boxmesh%idx%nr(2,1)
       do j=scdm%boxmesh%idx%nr(1,2),scdm%boxmesh%idx%nr(2,2)
          do k=scdm%boxmesh%idx%nr(1,3),scdm%boxmesh%idx%nr(2,3)
             ip = ip +1
             scdm%boxmesh%idx%lxyz(ip,1) = i
             scdm%boxmesh%idx%lxyz(ip,2) = j
             scdm%boxmesh%idx%lxyz(ip,3) = k
             scdm%boxmesh%idx%lxyz_inv(i,j,k) = ip
          enddo
       enddo
    enddo
    call mesh_cube_map_init(scdm%boxmesh%cube_map, scdm%boxmesh%idx, scdm%boxmesh%np_global)
    !
    ! create a cube object for the small box, with double size for coulomb truncation
    ! instead,this should be used to enlarge the box, but then need to keep track of dimension:
    !call mesh_double_box(scdm%boxmesh%sb, scdm%boxmesh, 2._8, temp)
    !
!    call cube_init(scdm%boxcube, scdm%boxmesh%idx%ll, scdm%boxmesh%sb, &
!               fft_type=FFT_REAL, fft_library=FFTLIB_FFTW, dont_optimize = .true.)
 
    ! without nfft we have to double the box 
#ifndef HAVE_NFFT
   if(der%mesh%sb%periodic_dim.gt.0) call messages_not_implemented("periodic SSCDM  without NFFT library")  
   box(:) = scdm%boxmesh%idx%ll(:)*2
   call cube_init(scdm%boxcube, box, scdm%boxmesh%sb,fft_type=FFT_REAL, fft_library=FFTLIB_FFTW)
# else ! nfft case
    box(:) = scdm%boxmesh%idx%ll(:) +2
    if(der%mesh%sb%periodic_dim.eq.3) then
       !enlargement factor to fit he simulationbox boundary
! ??? not sure
       enlarge = der%mesh%sb%lsize(1)/(2*scdm%box_size+1)
       !
    else ! non-periodic case
       enlarge = M_TWO
    endif
    call cube_init(scdm%boxcube, box, scdm%boxmesh%sb, &
               fft_type=FFT_COMPLEX, fft_library=FFTLIB_NFFT, &
               tp_enlarge=enlarge,spacing=der%mesh%spacing)
!call cube_init(scdm%boxcube, box*2, scdm%boxmesh%sb,fft_type=FFT_REAL, fft_library=FFTLIB_FFTW)
#endif 
    !
    ! Joseba recommends including this
    !if (der%mesh%parallel_in_domains .and. this%cube%parallel_in_domains) then
    !    call mesh_cube_parallel_map_init(this%mesh_cube_map, der%mesh, this%cube)
    !end if
    !
    ! set up poisson solver used for the exchange operator with scdm states
    ! this replictaes poisson_kernel_init()
    scdm%poisson%poisson_soft_coulomb_param = M_ZERO
    call poisson_fft_init(scdm%poisson_fft, scdm%boxmesh, scdm%boxcube, kernel=POISSON_FFT_KERNEL_SPH)
!call poisson_fft_init(scdm%poisson_fft, scdm%boxmesh, scdm%boxcube, kernel=POISSON_FFT_KERNEL_NOCUT)
    !
    ! create poisson object
    SAFE_ALLOCATE(scdm%poisson%der)
    SAFE_ALLOCATE(scdm%poisson%der%mesh)
    scdm%poisson%der%mesh = scdm%boxmesh
    scdm%poisson%method = POISSON_FFT
    scdm%poisson%kernel = POISSON_FFT_KERNEL_SPH
    scdm%poisson%cube = scdm%boxcube
    scdm%poisson%fft_solver = scdm%poisson_fft
    !
    ! set flag to do this only once
    scdm_is_init = .true.
    !
    call messages_write('done SCDM init')
    !
  end subroutine scdm_init

  subroutine dRRQR(n,np,KSt,JPVT)
    ! wrapper routine for real rank-revealing QR decompisition
    ! of the n*np matrix KSt, returning the pivot vector JPVT
    !
    integer, intent(in)  :: np,n
    FLOAT, intent(inout) :: KSt(:,:)
    integer, intent(out) :: JPVT(np)
    !
    integer            :: lwork,INFO
    FLOAT              :: TAU(n)
    FLOAT, allocatable ::  work(:)
    !
    ! dummy call to obtain dimension of work
    allocate(work(1))
    call DGEQP3( n,np, KSt, n, JPVT, TAU, WORK, -1, INFO )
    if(INFO.ne.0) then
       print *, 'Illegal argument in DGEQP3: ', INFO
       stop
    endif
    ! Note: scalapack routine is called P?GEQPF()
    !
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    !
    JPVT(:) = 0
    TAU(:) = 0.
    ! actual call
    call DGEQP3( n,np, KSt, n, JPVT, TAU, WORK, LWORK, INFO )
    if(INFO.ne.0)then
       print *, 'Illegal argument in DGEQP3: ', INFO
       stop
    endif
    !
  end subroutine dRRQR

  subroutine zRRQR(n,np,KSt,JPVT)
    ! wrapper routine for complex rank-revealing QR decompisition
    ! of the n*np matrix KSt, returning the pivot vector JPVT
    !
    integer, intent(in)  :: np,n
    CMPLX, intent(inout) :: KSt(:,:)
    integer, intent(out) :: JPVT(np)
    !
    integer            :: lwork,INFO
    CMPLX              :: TAU(n)
    CMPLX, allocatable :: work(:)
    FLOAT              :: rwork(2*np)
    !
    ! dummy call to obtain dimension of work
    allocate(work(1))
    call ZGEQP3( n,np, KSt, n, JPVT, TAU, WORK, -1, RWORK, INFO )
    if(INFO.ne.0) then
       print *, 'Illegal argument in ZGEQP3: ', INFO
       stop
    endif
    ! Note: scalapack routine is called P?GEQPF()
    !
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    !
    JPVT(:) = 0
    TAU(:) = 0.
    ! actual call
    call ZGEQP3( n,np, KSt, n, JPVT, TAU, WORK, LWORK, RWORK, INFO )
    if(INFO.ne.0)then
       print *, 'Illegal argument in ZGEQP3: ', INFO
       stop
    endif
    !
  end subroutine zRRQR

  subroutine check_periodic_box(idx,center,size,periodic)
    ! check if there are points outside simulation cell (index range of idx) by
    ! checking the corners only. This is inteded for rectangular cells
    ! should be generalized to arbitrary shapes (but then need to check the faces)
    !
    type(index_t),    intent(in)  :: idx
    integer, intent(in)  :: center(:)
    integer, intent(in)  :: size
    logical, intent(out) :: periodic
    !
    ! internal
    integer :: ix(3), corner(3,8), i1, idim
    !
    periodic = .false.
    !
    ! make the sign pattern for corners
    corner(:,1) = (/1,1,1/)
    corner(:,2) = (/1,1,-1/)
    corner(:,3) = (/1,-1,1/)
    corner(:,4) = (/-1,1,1/)
    corner(:,5) = (/1,-1,-1/)
    corner(:,6) = (/-1,1,-1/)
    corner(:,7) = (/-1,-1,1/)
    corner(:,8) = (/-1,-1,-1/)
    !
    do i1=1,8
       !
       ix(:)=center(:) + size*corner(:,i1)
       !
       do idim=1,3
          if(ix(idim).lt.idx%nr(1,idim).or.ix(idim).gt.idx%nr(2,idim)) then
             periodic = .true. 
             return
          endif
       enddo
       !
    enddo
    !
  end subroutine check_periodic_box



#include "undef.F90"
#include "real.F90"
#include "scdm_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "scdm_inc.F90"

end module scdm_m
