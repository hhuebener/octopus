!! Copyright (C) 2002-2010 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: density.F90 12493 2014-09-11 15:02:47Z dstrubbe $

#include "global.h"

module density_m
  use blas_m
  use batch_m
  use c_pointer_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use comm_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use io_m
  use kpoints_m
  use loct_m
  use math_m
  use mesh_m
  use messages_m
  use multigrid_m
  use multicomm_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use mpi_lib_m
  use opencl_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_m
  use states_dim_m
  use symmetrizer_m
  use types_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  private

  public ::                           &
    density_calc_t,                   &
    density_calc_init,                &
    density_calc_accumulate,          &
    density_calc_end,                 &
    density_calc,                     &
    states_freeze_orbitals,           &
    states_total_density

  type density_calc_t
    FLOAT,          pointer :: density(:, :)
    FLOAT,          pointer :: Imdensity(:, :)
    type(states_t), pointer :: st
    type(grid_t),   pointer :: gr
    type(opencl_mem_t)      :: buff_density
    integer                 :: pnp
    logical                 :: packed
  end type density_calc_t

contains
  
  subroutine density_calc_init(this, st, gr, density, Imdensity)
    type(density_calc_t),           intent(out)   :: this
    type(states_t),       target,   intent(in)    :: st
    type(grid_t),         target,   intent(in)    :: gr
    FLOAT,                target,   intent(out)   :: density(:, :)
    FLOAT, optional,      target,   intent(out)   :: Imdensity(:, :)

    PUSH_SUB(density_calc_init)

    this%density => density
    this%st => st
    this%gr => gr

    this%density = M_ZERO

    if(present(Imdensity)) then
      this%Imdensity => Imdensity
      this%Imdensity = M_ZERO
    else 
      nullify(this%Imdensity)
    end if      

    this%packed = .false.

    POP_SUB(density_calc_init)
  end subroutine density_calc_init

  ! ---------------------------------------------------

  subroutine density_calc_pack(this)
    type(density_calc_t),           intent(out)   :: this

    PUSH_SUB(density_calc_pack)

    this%packed = .true.
#ifdef HAVE_OPENCL    
    this%pnp = opencl_padded_size(this%gr%mesh%np)
    call opencl_create_buffer(this%buff_density, CL_MEM_READ_WRITE, TYPE_FLOAT, this%pnp*this%st%d%nspin)
    
    ! set to zero
    call opencl_set_buffer_to_zero(this%buff_density, TYPE_FLOAT, this%pnp*this%st%d%nspin)
#endif
    POP_SUB(density_calc_pack)
  end subroutine density_calc_pack

  ! ---------------------------------------------------

  subroutine density_calc_accumulate(this, ik, psib, psibL)
    type(density_calc_t), target, intent(inout) :: this
    integer,                      intent(in)    :: ik
    type(batch_t),                intent(inout) :: psib
    type(batch_t), optional,      intent(inout) :: psibL !< Left states

    integer :: ist, ip, ispin
    CMPLX   :: term, psi1, psi2
    FLOAT, pointer :: crho(:)
    FLOAT, pointer :: Imcrho(:)  
    FLOAT, allocatable :: frho(:), weight(:)
    type(profile_t), save :: prof
    logical :: correct_size, cmplxscl
#ifdef HAVE_OPENCL
    integer            :: wgsize
    type(opencl_mem_t) :: buff_weight
    type(cl_kernel)    :: kernel
#endif

    PUSH_SUB(density_calc_accumulate)
    call profiling_in(prof, "CALC_DENSITY")

    correct_size = ubound(this%density, dim = 1) == this%gr%fine%mesh%np .or. &
         ubound(this%density, dim = 1) == this%gr%fine%mesh%np_part
    ASSERT(correct_size)

    cmplxscl = associated(this%Imdensity)
    
    ispin = states_dim_get_spin_index(this%st%d, ik)

    SAFE_ALLOCATE(weight(1:psib%nst))
    forall(ist = 1:psib%nst) weight(ist) = this%st%d%kweights(ik)*this%st%occ(psib%states(ist)%ist, ik)

    if(this%st%d%ispin /= SPINORS) then 

      if(this%gr%have_fine_mesh) then
        SAFE_ALLOCATE(crho(1:this%gr%mesh%np_part))
        crho = M_ZERO
        if(cmplxscl) then
          SAFE_ALLOCATE(Imcrho(1:this%gr%mesh%np_part))
          Imcrho = M_ZERO
        end if
      else
        crho => this%density(:, ispin)
        if (cmplxscl) Imcrho => this%Imdensity(:, ispin)
      end if

      select case(batch_status(psib))
      case(BATCH_NOT_PACKED)
        if(states_are_real(this%st)) then
          do ist = 1, psib%nst
            forall(ip = 1:this%gr%mesh%np)
              crho(ip) = crho(ip) + weight(ist)*psib%states(ist)%dpsi(ip, 1)**2
            end forall
          end do
        else
          if(cmplxscl) then
            do ist = 1, psib%nst
              forall(ip = 1:this%gr%mesh%np)
                crho(ip)   = crho(ip)   + weight(ist) * real( psibL%states(ist)%zpsi(ip, 1) * psib%states(ist)%zpsi(ip, 1))
                Imcrho(ip) = Imcrho(ip) + weight(ist) * aimag(psibL%states(ist)%zpsi(ip, 1) * psib%states(ist)%zpsi(ip, 1))
              end forall
            end do
          else
            do ist = 1, psib%nst
              forall(ip = 1:this%gr%mesh%np)
                crho(ip) = crho(ip) + weight(ist)* &
                  (real(psib%states(ist)%zpsi(ip, 1), REAL_PRECISION)**2 + aimag(psib%states(ist)%zpsi(ip, 1))**2)
              end forall
            end do
          end if
        end if
      case(BATCH_PACKED)
        if(states_are_real(this%st)) then
          do ip = 1, this%gr%mesh%np
            do ist = 1, psib%nst
              crho(ip) = crho(ip) + weight(ist)*psib%pack%dpsi(ist, ip)**2
            end do
          end do
        else
          if(cmplxscl) then
            do ip = 1, this%gr%mesh%np
              do ist = 1, psib%nst
                crho(ip)   = crho(ip)   + weight(ist) * real( psibL%pack%zpsi(ist, ip) * psib%pack%zpsi(ist, ip))
                Imcrho(ip) = Imcrho(ip) + weight(ist) * aimag(psibL%pack%zpsi(ist, ip) * psib%pack%zpsi(ist, ip))
              end do
            end do
          else  
            do ip = 1, this%gr%mesh%np
              do ist = 1, psib%nst
                crho(ip) = crho(ip) + weight(ist)* &
                  (real(psib%pack%zpsi(ist, ip), REAL_PRECISION)**2 + aimag(psib%pack%zpsi(ist, ip))**2)
              end do
            end do
          end if
        end if
      case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
        if(.not. this%packed) call density_calc_pack(this)

        if(states_are_real(this%st)) then
          kernel = kernel_density_real
        else
          kernel = kernel_density_complex
        end if

        call opencl_create_buffer(buff_weight, CL_MEM_READ_ONLY, TYPE_FLOAT, psib%nst)
        call opencl_write_buffer(buff_weight, psib%nst, weight)

        call opencl_set_kernel_arg(kernel, 0, psib%nst)
        call opencl_set_kernel_arg(kernel, 1, this%gr%mesh%np)
        call opencl_set_kernel_arg(kernel, 2, this%pnp*(ispin - 1))
        call opencl_set_kernel_arg(kernel, 3, buff_weight)
        call opencl_set_kernel_arg(kernel, 4, psib%pack%buffer)
        call opencl_set_kernel_arg(kernel, 5, log2(psib%pack%size(1)))
        call opencl_set_kernel_arg(kernel, 6, this%buff_density)

        wgsize = opencl_kernel_workgroup_size(kernel)
        
        call opencl_kernel_run(kernel, (/pad(this%gr%mesh%np, wgsize)/), (/wgsize/))

        call opencl_finish()
        
        call opencl_release_buffer(buff_weight)
#endif
      end select

      if(this%gr%have_fine_mesh) then
        if(cmplxscl) then
          SAFE_ALLOCATE(frho(1:this%gr%fine%mesh%np))
          call dmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, crho, frho, order = 2)
          forall(ip = 1:this%gr%fine%mesh%np) this%density(ip, ispin) = this%density(ip, ispin) + frho(ip)    
          call dmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, Imcrho, frho, order = 2)
          forall(ip = 1:this%gr%fine%mesh%np) this%density(ip, ispin) = this%Imdensity(ip, ispin) + frho(ip)
          SAFE_DEALLOCATE_P(crho)
          SAFE_DEALLOCATE_P(Imcrho)
          SAFE_DEALLOCATE_A(frho)
        else  
          SAFE_ALLOCATE(frho(1:this%gr%fine%mesh%np))
          call dmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, crho, frho, order = 2)
          ! some debugging output that I will keep here for the moment, XA
          !      call dio_function_output(1, "./", "n_fine", this%gr%fine%mesh, frho, unit_one, ierr)
          !      call dio_function_output(1, "./", "n_coarse", this%gr%mesh, crho, unit_one, ierr)
          forall(ip = 1:this%gr%fine%mesh%np) this%density(ip, ispin) = this%density(ip, ispin) + frho(ip)
          SAFE_DEALLOCATE_P(crho)
          SAFE_DEALLOCATE_A(frho)
        end if
      end if

    else !SPINORS

      ! in this case wavefunctions are always complex
      ASSERT(.not. this%gr%have_fine_mesh)
      call batch_sync(psib)

      do ist = 1, psib%nst
        do ip = 1, this%gr%fine%mesh%np

          psi1 = psib%states(ist)%zpsi(ip, 1)
          psi2 = psib%states(ist)%zpsi(ip, 2)

          this%density(ip, 1) = this%density(ip, 1) + weight(ist)*(real(psi1, REAL_PRECISION)**2 + aimag(psi1)**2)
          this%density(ip, 2) = this%density(ip, 2) + weight(ist)*(real(psi2, REAL_PRECISION)**2 + aimag(psi2)**2)

          term = weight(ist)*psi1*conjg(psi2)
          this%density(ip, 3) = this%density(ip, 3) + real(term, REAL_PRECISION)
          this%density(ip, 4) = this%density(ip, 4) + aimag(term)

        end do
      end do
      
    end if

    SAFE_DEALLOCATE_A(weight)

    call profiling_out(prof)

    POP_SUB(density_calc_accumulate)
  end subroutine density_calc_accumulate

  ! ---------------------------------------------------

  subroutine density_calc_end(this)
    type(density_calc_t), intent(inout) :: this

    type(symmetrizer_t) :: symmetrizer
    FLOAT, allocatable :: tmpdensity(:)
    integer :: ispin
    type(profile_t), save :: reduce_prof
#ifdef HAVE_OPENCL
    integer :: ip
    FLOAT, allocatable :: fdensity(:)
#endif

    PUSH_SUB(density_calc_end)

    if(this%packed) then
#ifdef HAVE_OPENCL
      SAFE_ALLOCATE(tmpdensity(1:this%gr%mesh%np_part))

      ! the density is in device memory
      do ispin = 1, this%st%d%nspin
        call opencl_read_buffer(this%buff_density, this%gr%mesh%np, tmpdensity, offset = (ispin - 1)*this%pnp)

        if(this%gr%have_fine_mesh) then
           SAFE_ALLOCATE(fdensity(1:this%gr%fine%mesh%np))
           call dmultigrid_coarse2fine(this%gr%fine%tt, this%gr%der, this%gr%fine%mesh, tmpdensity, fdensity, order = 2)
           forall(ip = 1:this%gr%fine%mesh%np) this%density(ip, ispin) = this%density(ip, ispin) + fdensity(ip)
           SAFE_DEALLOCATE_A(fdensity)
        else
           forall(ip = 1:this%gr%mesh%np) this%density(ip, ispin) = this%density(ip, ispin) + tmpdensity(ip)
        end if

      end do

      this%packed = .false.
      call opencl_release_buffer(this%buff_density)
      SAFE_DEALLOCATE_A(tmpdensity)
#endif
    end if

    ! reduce over states and k-points
    if(this%st%parallel_in_states .or. this%st%d%kpt%parallel) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      call comm_allreduce(this%st%st_kpt_mpi_grp%comm, this%density, dim = (/this%gr%fine%mesh%np, this%st%d%nspin/))
      call profiling_out(reduce_prof)
    end if

    if(this%st%symmetrize_density) then
      SAFE_ALLOCATE(tmpdensity(1:this%gr%fine%mesh%np))
      call symmetrizer_init(symmetrizer, this%gr%fine%mesh)

      do ispin = 1, this%st%d%nspin
        call dsymmetrizer_apply(symmetrizer, field = this%density(:, ispin), symmfield = tmpdensity)
        this%density(1:this%gr%fine%mesh%np, ispin) = tmpdensity(1:this%gr%fine%mesh%np)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(tmpdensity)
    end if

    POP_SUB(density_calc_end)
  end subroutine density_calc_end


  ! ---------------------------------------------------------
  !> Computes the density from the orbitals in st. 
  subroutine density_calc(st, gr, density, Imdensity)
    type(states_t),          intent(inout)  :: st
    type(grid_t),            intent(in)     :: gr
    FLOAT,                   intent(out)    :: density(:, :)
    FLOAT, optional,         intent(out)    :: Imdensity(:, :)

    integer :: ik, ib
    type(density_calc_t) :: dens_calc
    logical :: cmplxscl

    PUSH_SUB(density_calc)

    ASSERT(ubound(density, dim = 1) == gr%fine%mesh%np .or. ubound(density, dim = 1) == gr%fine%mesh%np_part)

    cmplxscl = present(Imdensity)
    
    if (cmplxscl) then
      call density_calc_init(dens_calc, st, gr, density, Imdensity) 
    else
      call density_calc_init(dens_calc, st, gr, density)
    end if
    
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        if(cmplxscl) then
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik), st%psibL(ib, ik))
        else
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
        end if
      end do
    end do

    call density_calc_end(dens_calc)

    POP_SUB(density_calc)
  end subroutine density_calc

  ! ---------------------------------------------------------

  subroutine states_freeze_orbitals(st, gr, mc, n)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(multicomm_t), intent(in)    :: mc
    integer,           intent(in)    :: n

    integer :: ist, ik, ib, nblock
    type(states_t) :: staux
    CMPLX, allocatable :: psi(:, :, :)
    type(batch_t)  :: psib
    type(density_calc_t) :: dens_calc
#ifdef HAVE_MPI
    integer :: status(MPI_STATUS_SIZE)
#endif

    PUSH_SUB(states_freeze_orbitals)

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call messages_fatal(2)
    end if

    ASSERT(.not. st%parallel_in_states)

    if(.not.associated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np, 1:st%d%dim))
    end if

    call density_calc_init(dens_calc, st, gr, st%frozen_rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      if(n < st%st_start .or. n > st%st_end) cycle

      do ib =  st%group%block_start, st%group%block_end
        if(states_block_max(st, ib) <= n) then

          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
          if(states_block_max(st, ib) == n) exit

        else 

          nblock = n - states_block_min(st, ib) + 1

          SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim, 1:nblock))

          do ist = 1, nblock
            call states_get_state(st, gr%mesh, states_block_min(st, ib) + ist - 1, ik, psi(:, :, ist)) 
          end do

          call batch_init(psib, st%d%dim, states_block_min(st, ib), n, psi)
          call density_calc_accumulate(dens_calc, ik, psib)
          call batch_end(psib)
          exit

        end if

      end do
    end do

    call density_calc_end(dens_calc)

    call states_copy(staux, st)

    st%nst = st%nst - n

    call states_deallocate_wfns(st)
    call states_distribute_nodes(st, mc)
    call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX)

#if defined(HAVE_MPI) 

    if(staux%parallel_in_states) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = staux%st_start, staux%st_end
          if(ist <= n) cycle
          if(.not.state_is_local(st, ist-n)) then
            call mpi_send(staux%zpsi(1, 1, ist, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, staux%node(ist), &
              ist, st%mpi_grp%comm, mpi_err)

            call mpi_recv(st%zpsi(1, 1, ist-n, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, st%node(ist-n), &
              ist, st%mpi_grp%comm, status, mpi_err)
          else
            st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
          end if
   
        end do
      end do
   else
     do ik = st%d%kpt%start, st%d%kpt%end
       do ist = staux%st_start, staux%st_end
         if(ist <= n) cycle
         st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
       end do
     end do
   end if

#else

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
      end do
    end do

#endif

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)
    st%occ      = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%occ(ist, ik) = staux%occ(n+ist, ik)
        st%eigenval(ist, ik) = staux%eigenval(n+ist, ik)
      end do
    end do

    call states_end(staux)
    POP_SUB(states_freeze_orbitals)
  end subroutine states_freeze_orbitals


  ! ---------------------------------------------------------
  !> this routine calculates the total electronic density,
  !! which is the sum of the part coming from the orbitals, the
  !! non-linear core corrections and the frozen orbitals
  subroutine states_total_density(st, mesh, rho, Imrho)
    type(states_t),  intent(in)  :: st
    type(mesh_t),    intent(in)  :: mesh
    FLOAT,           intent(out) :: rho(:,:)
    FLOAT, optional, pointer, intent(out) :: Imrho(:,:)

    integer :: is, ip
    logical :: cmplxscl
    
    PUSH_SUB(states_total_density)

    cmplxscl = .false.
    if(present(Imrho)) then
      ASSERT(associated(Imrho))
      cmplxscl = .true.
    end if

    if(.not. cmplxscl) then
      forall(ip = 1:mesh%np, is = 1:st%d%nspin)
        rho(ip, is) = st%rho(ip, is)
      end forall

      if(associated(st%rho_core)) then
        forall(ip = 1:mesh%np, is = 1:st%d%spin_channels)
          rho(ip, is) = rho(ip, is) + st%rho_core(ip)/st%d%nspin
        end forall
      end if

      ! Add, if it exists, the frozen density from the inner orbitals.
      if(associated(st%frozen_rho)) then
        forall(ip = 1:mesh%np, is = 1:st%d%spin_channels)
          rho(ip, is) = rho(ip, is) + st%frozen_rho(ip, is)
        end forall
      end if
  
    else

      forall(ip = 1:mesh%np, is = 1:st%d%nspin)
        rho(ip, is)   = st%zrho%Re(ip, is)
        Imrho(ip, is) = st%zrho%Im(ip, is)
      end forall

      if(associated(st%rho_core)) then
        forall(ip = 1:mesh%np, is = 1:st%d%spin_channels)
          rho(ip, is)   = rho(ip, is)   + st%rho_core(ip)/st%d%nspin
          Imrho(ip, is) = Imrho(ip, is) + st%Imrho_core(ip)/st%d%nspin          
        end forall
      end if

      ! Add, if it exists, the frozen density from the inner orbitals.
      if(associated(st%frozen_rho)) then
        forall(ip = 1:mesh%np, is = 1:st%d%spin_channels)
          rho(ip, is) = rho(ip, is) + st%frozen_rho(ip, is)
          Imrho(ip, is) = Imrho(ip, is) + st%Imfrozen_rho(ip, is)
        end forall
      end if
      
    end if

    POP_SUB(states_total_density)
  end subroutine states_total_density

end module density_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
