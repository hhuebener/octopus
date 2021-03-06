!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: v_ks.F90 12633 2014-11-23 21:51:11Z adelgado $

#include "global.h"
 
module v_ks_m
  use berry_m
  use current_m
  use datasets_m
  use density_m
  use derivatives_m
  use energy_m
  use energy_calc_m
  use epot_m 
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_function_m
  use lalg_basic_m
  use magnetic_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use multigrid_m
  use parser_m
  use poisson_m
  use profiling_m
  use pcm_m 
  use simul_box_m
  use states_m
  use states_dim_m
  use unit_system_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)
  use xc_functl_m
  use xc_ks_inversion_m
  use xc_OEP_m

  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    v_ks_calc,          &
    v_ks_calc_t,        &
    v_ks_calc_start,    &
    v_ks_calc_finish,   &
    v_ks_freeze_hxc,    &
    v_ks_calculate_current

  integer, parameter, public :: &
    SIC_NONE   = 1,     &  !< no self-interaction correction
    SIC_PZ     = 2,     &  !< Perdew-Zunger SIC (OEP way)
    SIC_AMALDI = 3,     &  !< Amaldi correction term
    SIC_ADSIC  = 4         !< Averaged density SIC

  type v_ks_calc_t
    private
    logical                       :: calculating
    logical                       :: time_present
    FLOAT                         :: time
    FLOAT,                pointer :: density(:, :)
    logical                       :: total_density_alloc
    FLOAT,                pointer :: total_density(:)
    FLOAT                         :: amaldi_factor
    type(energy_t),       pointer :: energy
    type(states_t),       pointer :: hf_st
    FLOAT,                pointer :: vxc(:, :)
    FLOAT,                pointer :: vtau(:, :)
    FLOAT,                pointer :: axc(:, :, :)
    FLOAT,                pointer :: vberry(:, :)
    FLOAT,                pointer :: a_ind(:, :)
    FLOAT,                pointer :: b_ind(:, :)
    logical                       :: calc_energy
    !cmplxscl
    FLOAT,                pointer :: Imdensity(:, :)
    FLOAT,                pointer :: Imtotal_density(:)
    FLOAT,                pointer :: Imvxc(:, :)
    FLOAT,                pointer :: Imvtau(:, :)
    FLOAT,                pointer :: Imaxc(:, :, :)
    FLOAT,                pointer :: Imvberry(:, :)
    
  end type v_ks_calc_t

  type v_ks_t
    integer :: theory_level

    logical :: frozen_hxc !< For RPA and SAE calculations.

    integer                  :: xc_family  !< the XC stuff
    integer                  :: sic_type   !< what kind of self-interaction correction to apply
    type(xc_t)               :: xc
    type(xc_OEP_t)           :: oep
    type(xc_ks_inversion_t)  :: ks_inversion
    type(poisson_t), pointer :: hartree_solver
    logical                  :: new_hartree
    type(grid_t), pointer    :: gr
    type(v_ks_calc_t)        :: calc
    logical                  :: calculate_current
    type(current_t)          :: current_calculator
  end type v_ks_t

contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, gr, st, geo, mc)
    type(v_ks_t),         intent(out)   :: ks
    type(grid_t), target, intent(inout) :: gr
    type(states_t),       intent(in)    :: st
    type(geometry_t),     intent(inout) :: geo
    type(multicomm_t),    intent(in)    :: mc  

    PUSH_SUB(v_ks_init)

    !%Variable TheoryLevel
    !%Type integer
    !%Default dft
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with different "theory levels":
    !%Option independent_particles 2
    !% Particles will be considered as independent, <i>i.e.</i> as non-interacting.
    !% This mode is mainly used for testing purposes, as the code is usually 
    !% much faster with <tt>independent_particles</tt>.
    !%Option hartree 1
    !% Calculation within the Hartree method (experimental). Note that, contrary to popular
    !% belief, the Hartree potential is self-interaction-free. Therefore, this run 
    !% mode will not yield the same result as <tt>dft</tt> without exchange-correlation.
    !%Option hartree_fock 3
    !% This is the traditional Hartree-Fock scheme. Like the Hartree scheme, it is fully
    !% self-interaction-free. This mode is extremely slow. It is often more convenient
    !% to use <tt>dft</tt> within the OEP scheme to get similar (but not the same) results.
    !% Note that within this scheme you can use a correlation functional, or a hybrid
    !% functional (see <tt>XCFunctional</tt>). In the latter case, you will be following the
    !% quantum-chemistry recipe to use hybrids.
    !%Option dft 4
    !% This is the default density-functional theory scheme. Note that you can also use 
    !% hybrids in this scheme, but they will be handled the "DFT" way, <i>i.e.</i>, solving the
    !% OEP equation.
    !%Option classical 5
    !% (Experimental) Only the classical interaction between ions is
    !% considered. This is mainly for testing.
    !%Option rdmft 7 
    !% (Not fully implemented) Reduced Density Matrix functional theory
    !%End
    call parse_integer(datasets_check('TheoryLevel'), KOHN_SHAM_DFT, ks%theory_level)
    if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call input_error('TheoryLevel')

    call messages_obsolete_variable('NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable('HartreeFock', 'TheoryLevel')
    if(ks%theory_level == CLASSICAL) call messages_experimental('Classical theory level')
    if(ks%theory_level == RDMFT ) call messages_experimental('RDMFT theory level')
    
    ks%xc_family = XC_FAMILY_NONE
    ks%sic_type  = SIC_NONE
    
    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      ks%sic_type = SIC_NONE
    case(HARTREE)
      call messages_experimental("Hartree theory level")
      if(gr%mesh%sb%periodic_dim == gr%mesh%sb%dim) &
        call messages_experimental("Hartree in fully periodic system")
      if(gr%mesh%sb%kpoints%full%npoints > 1) &
        call messages_not_implemented("Hartree with k-points")

    case(HARTREE_FOCK)
      if(gr%mesh%sb%kpoints%full%npoints > 1) &
        call messages_not_implemented("Hartree-Fock with k-points")

      ! initialize XC modules
      call xc_init(ks%xc, gr%mesh%sb%dim, gr%mesh%sb%periodic_dim, st%qtot, hartree_fock=.true.)
      ks%xc_family = ks%xc%family
      ks%sic_type = SIC_NONE

    case(KOHN_SHAM_DFT)
      ! initialize XC modules
      call xc_init(ks%xc, gr%mesh%sb%dim, gr%mesh%sb%periodic_dim, st%qtot, hartree_fock=.false.)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA) /= 0) then

        !%Variable SICCorrection
        !%Type integer
        !%Default sic_none
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls which form of self-interaction correction to use. Note that
        !% this correction will be applied to the functional chosen by <tt>XCFunctional</tt>.
        !%Option sic_none 1
        !% No self-interaction correction.
        !%Option sic_pz 2
        !% Perdew-Zunger SIC, handled by the OEP technique.
        !%Option sic_amaldi 3
        !% Amaldi correction term.
        !%Option sic_adsic 4
        !% Average-density SIC.
        !% C. Legrand et al. J. Phys. B 35, 1115 (2002). 
        !%End
        call parse_integer(datasets_check('SICCorrection'), sic_none, ks%sic_type)
        if(.not. varinfo_valid_option('SICCorrection', ks%sic_type)) call input_error('SICCorrection')

        ! Perdew-Zunger corrections
        if(ks%sic_type == SIC_PZ) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)

      else
        ks%sic_type = SIC_NONE
      end if

      if(iand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        call xc_oep_init(ks%oep, ks%xc_family, gr, st)
      endif
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_init(ks%ks_inversion, ks%xc_family, gr, geo, mc)
      endif
    end select

    ks%frozen_hxc = .false.

    call v_ks_write_info(ks, stdout)

    ks%new_hartree = .false.
    nullify(ks%hartree_solver)
    if(ks%theory_level /= INDEPENDENT_PARTICLES) then
      if(gr%have_fine_mesh) then
        ks%new_hartree = .true.
        SAFE_ALLOCATE(ks%hartree_solver)
        call poisson_init(ks%hartree_solver, gr%fine%der, mc, &
          label = " (fine mesh)", theta = st%cmplxscl%theta)
      else
        ks%hartree_solver => psolver
      end if
    end if

    ks%gr => gr
    ks%calc%calculating = .false.

    ks%calculate_current = .false.
    call current_init(ks%current_calculator)

    POP_SUB(v_ks_init)
  end subroutine v_ks_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks, gr)
    type(v_ks_t),     intent(inout) :: ks
    type(grid_t),     intent(inout) :: gr

    PUSH_SUB(v_ks_end)

    call current_end(ks%current_calculator)

    select case(ks%theory_level)
    case(KOHN_SHAM_DFT)
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_end(ks%ks_inversion, gr)
      endif
      if(iand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        call xc_oep_end(ks%oep)
      endif
      call xc_end(ks%xc)
    end select

    if(ks%new_hartree) then
      call poisson_end(ks%hartree_solver)
      SAFE_DEALLOCATE_P(ks%hartree_solver)
    end if

    POP_SUB(v_ks_end)
  end subroutine v_ks_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_t), intent(in) :: ks
    integer,      intent(in) :: iunit

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(v_ks_write_info)

    call messages_print_stress(iunit, "Theory Level")
    call messages_print_var_option(iunit, "TheoryLevel", ks%theory_level)

    select case(ks%theory_level)
    case(HARTREE_FOCK)
      write(iunit, '(1x)')
      call xc_write_info(ks%xc, iunit)

    case(KOHN_SHAM_DFT)
      write(iunit, '(1x)')
      call xc_write_info(ks%xc, iunit)

      write(iunit, '(1x)')
      call messages_print_var_option(iunit, 'SICCorrection', ks%sic_type)

      if(iand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        call xc_oep_write_info(ks%oep, iunit)
      end if
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_write_info(ks%ks_inversion, iunit)
      end if
        

    end select

    call messages_print_stress(iunit)

    POP_SUB(v_ks_write_info)
  end subroutine v_ks_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_calc(ks, hm, st, geo, calc_eigenval, time, calc_berry, calc_energy)
    type(v_ks_t),               intent(inout) :: ks
    type(hamiltonian_t),        intent(inout) :: hm
    type(states_t),             intent(inout) :: st
    type(geometry_t),           intent(in)    :: geo
    logical,          optional, intent(in)    :: calc_eigenval
    FLOAT,            optional, intent(in)    :: time
    logical,          optional, intent(in)    :: calc_berry !< use this before wfns initialized
    logical,          optional, intent(in)    :: calc_energy

    call v_ks_calc_start(ks, hm, st, geo, time, calc_berry, calc_energy)
    call v_ks_calc_finish(ks, hm)

    if(optional_default(calc_eigenval, .false.)) then

      if(ks%gr%ob_grid%open_boundaries .and. .not. present(time)) then
        ! We know the eigenvalues.
        st%eigenval(1:st%nst, 1:st%d%nik) = st%ob_eigenval(1:st%nst, 1:st%d%nik)
      else
        call energy_calc_eigenvalues(hm, ks%gr%der, st)
      end if
      
    end if
  end subroutine v_ks_calc

  ! --------------------------------------------------------- 

  !> This routine starts the calculation of the Kohn-Sham
  !! potential. The routine v_ks_calc_finish must be called to finish
  !! the calculation. The argument hm is not modified. The argument st
  !! can be modified after the function have been used.
  subroutine v_ks_calc_start(ks, hm, st, geo, time, calc_berry, calc_energy) 
    type(v_ks_t),            target,   intent(inout) :: ks 
    type(hamiltonian_t),     target,   intent(in)    :: hm !< This MUST be intent(in), changes to hm are done in v_ks_calc_finish.
    type(states_t),                    intent(inout) :: st
    type(geometry_t) ,                 intent(in)    :: geo
    FLOAT,                   optional, intent(in)    :: time 
    logical,                 optional, intent(in)    :: calc_berry !< Use this before wfns initialized.
    logical,                 optional, intent(in)    :: calc_energy

    type(profile_t), save :: prof
    type(energy_t), pointer :: energy
    logical  :: cmplxscl

    PUSH_SUB(v_ks_calc_start)
    call profiling_in(prof, "KOHN_SHAM_CALC")
    cmplxscl = hm%cmplxscl%space

    ASSERT(.not. ks%calc%calculating)
    ks%calc%calculating = .true.

    if(in_debug_mode) then
      write(message(1), '(a)') 'Debug: Calculating Kohn-Sham potential.'
      call messages_info(1)
    end if

    if(present(time)) then
      ks%calc%time = time
    else
      ks%calc%time = M_ZERO
    endif
    ks%calc%calc_energy = optional_default(calc_energy, .true.)

    nullify(ks%calc%vberry)
    nullify(ks%calc%Imvberry) !cmplxscl
    if(associated(hm%vberry)) then
      SAFE_ALLOCATE(ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin))
      SAFE_ALLOCATE(ks%calc%Imvberry(1:ks%gr%mesh%np, 1:hm%d%nspin)) !cmplxscl
      if(optional_default(calc_berry, .true.)) then
        call berry_potential(st, ks%gr%mesh, hm%ep%E_field, ks%calc%vberry)
      else
        ! before wfns are initialized, cannot calculate this term
        ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = M_ZERO
        ks%calc%Imvberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = M_ZERO !cmplxscl
      endif
    endif

    ! If the Hxc term is frozen, there is nothing more to do (WARNING: MISSING ks%calc%energy%intnvxc)
    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_start)
      return
    end if

    SAFE_ALLOCATE(ks%calc%energy)
    energy => ks%calc%energy

    call energy_copy(hm%energy, ks%calc%energy)

    energy%intnvxc = M_ZERO
    energy%Imintnvxc = M_ZERO !cmplxscl

    ! check whether we should introduce the Amaldi SIC correction
    ks%calc%amaldi_factor = M_ONE
    if(ks%sic_type == SIC_AMALDI) ks%calc%amaldi_factor = (st%qtot - M_ONE)/st%qtot

    nullify(ks%calc%density, ks%calc%total_density)
    nullify(ks%calc%vxc, ks%calc%vtau, ks%calc%axc)
    !cmplxscl
    nullify(ks%calc%Imdensity, ks%calc%Imtotal_density)
    nullify(ks%calc%Imvxc, ks%calc%Imvtau, ks%calc%Imaxc)


    if(ks%theory_level /= INDEPENDENT_PARTICLES .and. ks%calc%amaldi_factor /= M_ZERO) then

      call calculate_density()

      if(poisson_is_async(ks%hartree_solver)) then
        if(.not. cmplxscl) then
          call dpoisson_solve_start(ks%hartree_solver, ks%calc%total_density)
        else
          call zpoisson_solve_start(ks%hartree_solver, ks%calc%total_density + M_zI * ks%calc%Imtotal_density)
        end if
      end if

      if(ks%theory_level /= HARTREE .and. ks%theory_level /= RDMFT) call v_a_xc(geo, hm)
    else
      ks%calc%total_density_alloc = .false.
    end if

    if(ks%calculate_current) then
      call states_allocate_current(st, ks%gr)
      call current_calculate(ks%current_calculator, ks%gr, hm, geo, st, st%current)
    end if

    nullify(ks%calc%hf_st) 
    if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK .or. ks%theory_level == RDMFT) then
      SAFE_ALLOCATE(ks%calc%hf_st)
      call states_copy(ks%calc%hf_st, st)
    end if

    ! Calculate the vector potential induced by the electronic current.
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    nullify(ks%calc%a_ind, ks%calc%b_ind)
    if(hm%self_induced_magnetic) then
      SAFE_ALLOCATE(ks%calc%a_ind(1:ks%gr%mesh%np_part, 1:ks%gr%sb%dim))
      SAFE_ALLOCATE(ks%calc%b_ind(1:ks%gr%mesh%np_part, 1:ks%gr%sb%dim))
      call magnetic_induced(ks%gr%der, st, ks%calc%a_ind, ks%calc%b_ind)
    end if

    call profiling_out(prof)
    POP_SUB(v_ks_calc_start)

  contains

    subroutine calculate_density()
      integer :: ip

      PUSH_SUB(v_ks_calc_start.calculate_density)

      ! get density taking into account non-linear core corrections
      SAFE_ALLOCATE(ks%calc%density(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      if (.not. cmplxscl) then
        call states_total_density(st, ks%gr%fine%mesh, ks%calc%density)
      else 
        SAFE_ALLOCATE(ks%calc%Imdensity(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        call states_total_density(st, ks%gr%fine%mesh, ks%calc%density, ks%calc%Imdensity)
      end if

      ! Amaldi correction
      if(ks%sic_type == SIC_AMALDI) &
        ks%calc%density = ks%calc%amaldi_factor*ks%calc%density

      nullify(ks%calc%total_density)
      if(cmplxscl) nullify(ks%calc%Imtotal_density)
      if(associated(st%rho_core) .or. hm%d%spin_channels > 1) then
        ks%calc%total_density_alloc = .true.

        SAFE_ALLOCATE(ks%calc%total_density(1:ks%gr%fine%mesh%np))
        if(cmplxscl) then
          SAFE_ALLOCATE(ks%calc%Imtotal_density(1:ks%gr%fine%mesh%np))
        endif

        forall(ip = 1:ks%gr%fine%mesh%np)
          ks%calc%total_density(ip) = sum(ks%calc%density(ip, 1:hm%d%spin_channels))
        end forall
        if(cmplxscl) then
          forall(ip = 1:ks%gr%fine%mesh%np)
            ks%calc%Imtotal_density(ip) = sum(ks%calc%Imdensity(ip, 1:hm%d%spin_channels))
          end forall
        end if

        ! remove non-local core corrections
        if(associated(st%rho_core)) then
          forall(ip = 1:ks%gr%fine%mesh%np)
            ks%calc%total_density(ip) = ks%calc%total_density(ip) - st%rho_core(ip)*ks%calc%amaldi_factor
          end forall
        end if
      else
        ks%calc%total_density_alloc = .false.
        ks%calc%total_density => ks%calc%density(:, 1)
        if (cmplxscl) ks%calc%Imtotal_density => ks%calc%Imdensity(:, 1)
      end if

      POP_SUB(v_ks_calc_start.calculate_density)
    end subroutine calculate_density

    !ADSIC potential is:
    !V_ADSIC[n] = V_ks[n] - (V_h[n/N] - V_xc[n/N])
    subroutine add_adsic(cmplxscl, hm)
      logical, intent(in)                :: cmplxscl
      type(hamiltonian_t), intent(in)    :: hm

      integer        :: ip, ispin, ist
      FLOAT, pointer :: vxc_sic(:,:),  Imvxc_sic(:, :), vh_sic(:), rho(:, :), Imrho(:, :), qsp(:)
      CMPLX, pointer :: zrho_total(:), zvh_sic(:)
      
      PUSH_SUB(add_adsic)
      
      if(iand(hm%xc_family, XC_FAMILY_MGGA) /= 0) then
        call messages_not_implemented('ADSIC with MGGAs')
      end if
      
      SAFE_ALLOCATE(vxc_sic(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(vh_sic(1:ks%gr%mesh%np))
      SAFE_ALLOCATE(rho(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(qsp(1:st%d%nspin))
      
      vxc_sic = M_ZERO
      vh_sic = M_ZERO
      qsp = M_ZERO
      do ist = 1, st%nst
        qsp(:) = qsp(:)+ st%occ(ist, :) * st%d%kweights(:)
      end do

      do ispin = 1, st%d%nspin
        rho(:, ispin) = ks%calc%density(:, ispin) / qsp(ispin)
      end do

      if(cmplxscl) then
        SAFE_ALLOCATE(Imrho(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(Imvxc_sic(1:ks%gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(zvh_sic(1:ks%gr%mesh%np))
        
        do ispin = 1, st%d%nspin
          Imrho(:, ispin) = ks%calc%Imdensity(:, ispin) / qsp(ispin)
        end do

        Imvxc_sic = M_ZERO
        call xc_get_vxc_cmplx(ks%gr%fine%der, ks%xc, st%d%nspin, rho, Imrho, &
          vxc_sic, Imvxc_sic, st%cmplxscl%theta)

        SAFE_ALLOCATE(zrho_total(1:ks%gr%mesh%np))

        zrho_total(:) = sum(rho, 2) + M_zI * sum(Imrho, 2)
        zrho_total(:) = zrho_total(:) / st%qtot

        zvh_sic = M_ZERO
        call zpoisson_solve(ks%hartree_solver, zvh_sic, zrho_total)

        ks%calc%vxc = ks%calc%vxc - vxc_sic
        ks%calc%Imvxc = ks%calc%Imvxc - Imvxc_sic

        do ip=1, ks%gr%mesh%np
          ks%calc%vxc(ip, :) = ks%calc%vxc(ip, :) - real(zvh_sic(ip), REAL_PRECISION)
          ks%calc%Imvxc(ip, :) = ks%calc%Imvxc(ip, :) - aimag(zvh_sic(ip))
        end do

        SAFE_DEALLOCATE_P(Imrho)
        SAFE_DEALLOCATE_P(Imvxc_sic)
        SAFE_DEALLOCATE_P(zvh_sic)
        SAFE_DEALLOCATE_P(zrho_total)
      else
        call xc_get_vxc(ks%gr%fine%der, ks%xc, &
          st, rho, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, &
          vxc_sic)
        rho(:, 1) = ks%calc%total_density / st%qtot
        call dpoisson_solve(ks%hartree_solver, vh_sic, rho(:,1))
        ks%calc%vxc = ks%calc%vxc - vxc_sic
        forall(ip = 1:ks%gr%mesh%np) ks%calc%vxc(ip,:) = ks%calc%vxc(ip,:) - vh_sic(ip)
      end if

      SAFE_DEALLOCATE_P(vxc_sic)
      SAFE_DEALLOCATE_P(vh_sic)                                
      SAFE_DEALLOCATE_P(rho)
      SAFE_DEALLOCATE_P(qsp)

      POP_SUB(add_adsic)
    end subroutine add_adsic


    ! ---------------------------------------------------------
    subroutine v_a_xc(geo, hm)
      type(geometry_t),     intent(in) :: geo
      type(hamiltonian_t),  intent(in) :: hm 

      type(profile_t), save :: prof
      logical :: cmplxscl
      FLOAT :: factor
      CMPLX :: ctmp
      integer :: ispin

      PUSH_SUB(v_ks_calc_start.v_a_xc)
      call profiling_in(prof, "XC")

      cmplxscl = hm%cmplxscl%space

      energy%exchange = M_ZERO
      energy%correlation = M_ZERO
      energy%xc_j = M_ZERO
      !cmplxscl
      energy%Imexchange = M_ZERO
      energy%Imcorrelation = M_ZERO
      energy%Imxc_j = M_ZERO


      SAFE_ALLOCATE(ks%calc%vxc(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      ks%calc%vxc = M_ZERO
      if(cmplxscl) then
        SAFE_ALLOCATE(ks%calc%Imvxc(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        ks%calc%Imvxc = M_ZERO
      end if

      nullify(ks%calc%vtau)
      if(iand(hm%xc_family, XC_FAMILY_MGGA) /= 0) then
        SAFE_ALLOCATE(ks%calc%vtau(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        ks%calc%vtau = M_ZERO
        if(cmplxscl) then
          SAFE_ALLOCATE(ks%calc%Imvtau(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
          ks%calc%Imvtau = M_ZERO
        end if
      end if

      ! Get the *local* XC term
      if(hm%d%cdft) then
        call messages_not_implemented('Current-DFT')
      else if(ks%calc%calc_energy) then
        if(iand(hm%xc_family, XC_FAMILY_MGGA) /= 0) then
          if (cmplxscl) call messages_not_implemented('Complex Scaling with XC_FAMILY_MGGA')
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, &
            ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, ks%calc%vxc, &
            ex = energy%exchange, ec = energy%correlation, deltaxc = energy%delta_xc, vtau = ks%calc%vtau)
        else
          if(.not. cmplxscl) then
            call xc_get_vxc(ks%gr%fine%der, ks%xc, &
              st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, ks%calc%vxc, &
              ex = energy%exchange, ec = energy%correlation, deltaxc = energy%delta_xc)
          else
            call xc_get_vxc_cmplx(ks%gr%fine%der, ks%xc, st%d%ispin, ks%calc%density, ks%calc%Imdensity, &
              ks%calc%vxc, ks%calc%Imvxc, hm%cmplxscl%theta, ex = energy%exchange, ec = energy%correlation, &
              Imex = energy%Imexchange, Imec = energy%Imcorrelation)
          end if
        end if
      else
        if(iand(hm%xc_family, XC_FAMILY_MGGA) /= 0) then
          if (cmplxscl) call messages_not_implemented('Complex Scaling with XC_FAMILY_MGGA')
          call xc_get_vxc(ks%gr%fine%der, ks%xc, &
            st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, &
            ks%calc%vxc, vtau = ks%calc%vtau)
        else
          if(.not. cmplxscl) then
            call xc_get_vxc(ks%gr%fine%der, ks%xc, &
              st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, &
              ks%calc%vxc)
          else
            call xc_get_vxc_cmplx(ks%gr%fine%der, ks%xc, st%d%ispin, ks%calc%density, ks%calc%Imdensity, &
              ks%calc%vxc, ks%calc%Imvxc, hm%cmplxscl%theta)
          end if
        end if
      end if

      if (ks%sic_type == SIC_ADSIC) then
        call add_adsic(cmplxscl, hm)
      end if

      if(ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to be handled specially
        if(iand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
          if (cmplxscl) call messages_not_implemented('Complex Scaling with XC_FAMILY_OEP')
          if (states_are_real(st)) then
!HH                                                                                                     
if(.not.hm%EXX)   call dxc_oep_calc(ks%oep, ks%xc, (ks%sic_type == SIC_PZ),  &
                   ks%gr, hm, st, energy%exchange, energy%correlation, vxc = ks%calc%vxc)
            call dxc_oep_calc(ks%oep, ks%xc, (ks%sic_type == SIC_PZ),  &
              ks%gr, hm, st, energy%exchange, energy%correlation, vxc = ks%calc%vxc)
          else
!HH                                                                                                     
if(.not.hm%EXX)  call zxc_oep_calc(ks%oep, ks%xc, (ks%sic_type == SIC_PZ),  &
              ks%gr, hm, st, energy%exchange, energy%correlation, vxc = ks%calc%vxc)
          end if
        endif

        if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
          if (cmplxscl) call messages_not_implemented('Complex Scaling with XC_FAMILY_KS_INVERSION')
          ! Also treat KS inversion separately (not part of libxc)
          call xc_ks_inversion_calc(ks%ks_inversion, ks%gr, hm, st, vxc = ks%calc%vxc, time = ks%calc%time)
        endif
      end if

      if(ks%calc%calc_energy) then
        ! Now we calculate Int[n vxc] = energy%intnvxc
        energy%intnvxc = M_ZERO
        energy%Imintnvxc = M_ZERO !cmplxscl

        if(hm%d%ispin == SPINORS .and. cmplxscl) &
          call messages_not_implemented('Complex Scaling with SPINORS')
        do ispin = 1, hm%d%nspin
          if(ispin <= 2) then
            factor = M_ONE
          else
            factor = M_TWO
          endif
          if (.not. cmplxscl) then
            energy%intnvxc = energy%intnvxc + factor * dmf_dotp(ks%gr%fine%mesh, st%rho(:, ispin), ks%calc%vxc(:, ispin))
          else
            ctmp = factor * zmf_dotp(ks%gr%fine%mesh, st%zrho%Re(:, ispin) + M_zI * st%zrho%Im(:, ispin), &
              ks%calc%vxc(:, ispin) + M_zI * ks%calc%Imvxc(:, ispin), dotu = .true.)
            energy%intnvxc = energy%intnvxc + real(ctmp)
            energy%Imintnvxc = energy%Imintnvxc + aimag(ctmp)          
          end if
        enddo
      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start.v_a_xc)
    end subroutine v_a_xc

  end subroutine v_ks_calc_start
  ! ---------------------------------------------------------

  subroutine v_ks_calc_finish(ks, hm)
    type(v_ks_t), target, intent(inout) :: ks
    type(hamiltonian_t),  intent(inout) :: hm

    integer :: ip, ispin

    PUSH_SUB(v_ks_calc_finish)

    ASSERT(ks%calc%calculating)
    ks%calc%calculating = .false.

    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_finish)
      return
    end if

    !change the pointer to the energy object
    SAFE_DEALLOCATE_P(hm%energy)
    hm%energy => ks%calc%energy

    if(associated(hm%vberry)) then
      hm%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin)
      SAFE_DEALLOCATE_P(ks%calc%vberry)
    endif

    if(hm%self_induced_magnetic) then
      hm%a_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim) = ks%calc%a_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim)
      hm%b_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim) = ks%calc%b_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim)

      SAFE_DEALLOCATE_P(ks%calc%a_ind)
      SAFE_DEALLOCATE_P(ks%calc%b_ind)
    end if

    if(ks%theory_level == INDEPENDENT_PARTICLES .or. ks%calc%amaldi_factor == M_ZERO) then

      hm%vhxc = M_ZERO
      hm%energy%intnvxc     = M_ZERO
      hm%energy%hartree     = M_ZERO
      hm%energy%exchange    = M_ZERO
      hm%energy%correlation = M_ZERO
      !cmplxscl
      if(hm%cmplxscl%space) then
        hm%Imvhxc = M_ZERO
        hm%energy%Imintnvxc     = M_ZERO
        hm%energy%Imhartree     = M_ZERO
        hm%energy%Imexchange    = M_ZERO
        hm%energy%Imcorrelation = M_ZERO
      end if
    else

      if(ks%theory_level /= HARTREE .and. ks%theory_level /= RDMFT ) then 
        if(ks%gr%have_fine_mesh) then
          do ispin = 1, hm%d%nspin
            call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, &
              ks%calc%vxc(:, ispin), hm%vxc(:, ispin), INJECTION)
            ! some debugging output that I will keep here for the moment, XA
            !          call dio_function_output(1, "./", "vxc_fine", ks%gr%fine%mesh, vxc(:, ispin), unit_one, ierr)
            !          call dio_function_output(1, "./", "vxc_coarse", ks%gr%mesh, hm%vxc(:, ispin), unit_one, ierr)
          end do
          if(hm%cmplxscl%space) then
            do ispin = 1, hm%d%nspin
              call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, &
                ks%calc%Imvxc(:, ispin), hm%Imvxc(:, ispin), INJECTION)
            end do            
          end if
          SAFE_DEALLOCATE_P(ks%calc%vxc)
          SAFE_DEALLOCATE_P(ks%calc%Imvxc) !cmplxscl
        else
          ! just change the pointer to avoid the copy
          SAFE_DEALLOCATE_P(hm%vxc)
          hm%vxc => ks%calc%vxc
          if(hm%cmplxscl%space) then
            SAFE_DEALLOCATE_P(hm%Imvxc)
            hm%Imvxc => ks%calc%Imvxc
          end if
        end if

        if(iand(hm%xc_family, XC_FAMILY_MGGA) /= 0) then
          do ispin = 1, hm%d%nspin
            call lalg_copy(ks%gr%fine%mesh%np, ks%calc%vtau(:, ispin), hm%vtau(:, ispin))
            if(hm%cmplxscl%space) call lalg_copy(ks%gr%fine%mesh%np, ks%calc%Imvtau(:, ispin), hm%Imvtau(:, ispin))
          end do
          SAFE_DEALLOCATE_P(ks%calc%vtau)
          SAFE_DEALLOCATE_P(ks%calc%Imvtau)          
        end if

      else
        hm%vxc = M_ZERO
        if(hm%cmplxscl%space) hm%Imvxc = M_ZERO
      end if

      hm%energy%hartree = M_ZERO
      hm%energy%Imhartree = M_ZERO
      call v_ks_hartree(ks, hm)


      ! Build Hartree + XC potential
     
      forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 1) = hm%vxc(ip, 1) + hm%vhartree(ip)
      if (hm%cmplxscl%space) forall(ip = 1:ks%gr%mesh%np) hm%Imvhxc(ip, 1) = hm%Imvxc(ip, 1) + hm%Imvhartree(ip)
      if(associated(hm%vberry)) then
        forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 1) = hm%vhxc(ip, 1) + hm%vberry(ip, 1)
      endif
      
      if(hm%d%ispin > UNPOLARIZED) then
        forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 2) = hm%vxc(ip, 2) + hm%vhartree(ip)
        if (hm%cmplxscl%space) forall(ip = 1:ks%gr%mesh%np) hm%Imvhxc(ip, 2) = hm%Imvxc(ip, 2) + hm%Imvhartree(ip)
        if(associated(hm%vberry)) then
          forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 2) = hm%vhxc(ip, 2) + hm%vberry(ip, 2)
        endif
      end if
      
      if(hm%d%ispin == SPINORS) then
        forall(ispin = 3:4, ip = 1:ks%gr%mesh%np) hm%vhxc(ip, ispin) = hm%vxc(ip, ispin)
        if (hm%cmplxscl%space) forall(ispin = 3:4, ip = 1:ks%gr%mesh%np) hm%Imvhxc(ip, ispin) = hm%Imvxc(ip, ispin)
      end if

      if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK .or. ks%theory_level == RDMFT) then

        ! swap the states object
        call states_end(hm%hf_st)
        SAFE_DEALLOCATE_P(hm%hf_st)
        hm%hf_st => ks%calc%hf_st

        select case(ks%theory_level)
        case(HARTREE_FOCK)
          hm%exx_coef = ks%xc%exx_coef
        case(HARTREE)
          hm%exx_coef = M_ONE
        case(RDMFT) 
          hm%exx_coef = M_ONE
        end select
      end if
      
    end if

    call hamiltonian_update(hm, ks%gr%mesh, time = ks%calc%time)

    SAFE_DEALLOCATE_P(ks%calc%density)
    SAFE_DEALLOCATE_P(ks%calc%Imdensity)
    if(ks%calc%total_density_alloc) then
      SAFE_DEALLOCATE_P(ks%calc%total_density)
      SAFE_DEALLOCATE_P(ks%calc%Imtotal_density)
    end if
    nullify(ks%calc%total_density)
    nullify(ks%calc%Imtotal_density)

    POP_SUB(v_ks_calc_finish)
  end subroutine v_ks_calc_finish

  ! --------------------------------------------------------- 
  !
  !> Hartree contribution to the KS potential. This function is
  !! designed to be used by v_ks_calc_finish and it cannot be called
  !! directly.
  !
  subroutine v_ks_hartree(ks, hm)
    type(v_ks_t),                intent(inout) :: ks
    type(hamiltonian_t), target, intent(inout) :: hm

    FLOAT, pointer :: pot(:), Impot(:), aux(:)
    CMPLX, pointer :: zpot(:)
    CMPLX :: ztmp

    PUSH_SUB(v_ks_hartree)

    ASSERT(associated(ks%hartree_solver))

    if(.not. ks%gr%have_fine_mesh) then
      pot => hm%vhartree
      if (hm%cmplxscl%space) then 
        Impot => hm%Imvhartree
        SAFE_ALLOCATE(zpot(1:size(Impot,1)))
      end if
    else
      if(.not. hm%cmplxscl%space) then
        SAFE_ALLOCATE(pot(1:ks%gr%fine%mesh%np_part))
        pot = M_ZERO
      else
        SAFE_ALLOCATE(aux(1:ks%gr%fine%mesh%np_part))
        SAFE_ALLOCATE(zpot(1:ks%gr%fine%mesh%np_part))
        zpot = M_z0
      end if
    end if

    if(.not. poisson_is_async(ks%hartree_solver)) then
      if (.not. hm%cmplxscl%space) then 
        ! solve the Poisson equation
        call dpoisson_solve(ks%hartree_solver, pot, ks%calc%total_density)
      else
        ! Solve the Poisson equation for the scaled density and coulomb potential
        call zpoisson_solve(ks%hartree_solver, zpot,&
          ks%calc%total_density + M_zI * ks%calc%Imtotal_density)
        pot   =   real(zpot)
        Impot =  aimag(zpot)
      end if
    else
      ! The calculation was started by v_ks_calc_start.
      if(.not. hm%cmplxscl%space) then
        call dpoisson_solve_finish(ks%hartree_solver, pot)
      else
        call zpoisson_solve_finish(ks%hartree_solver, zpot)
        pot   =   real(zpot)
        Impot =  aimag(zpot)
      end if
    end if

    !> PCM reaction field due to the electronic density
    if (hm%pcm%run_pcm) then
    !> Generates the real-space PCM potential due to electrons during the SCF calculation.
        call v_electrons_cav_li(hm%pcm%v_e, pot, hm%pcm)
        call pcm_charges(hm%pcm%q_e, hm%pcm%qtot_e, hm%pcm%v_e, hm%pcm%matrix, hm%pcm%n_tesserae) 
        call pcm_pot_rs( hm%pcm%v_e_rs, hm%pcm%q_e, hm%pcm%tess, hm%pcm%n_tesserae, ks%gr%mesh, hm%pcm%gaussian_width )

        ! Calculating the PCM term renormalizing the sum of the single-particle energies
        hm%energy%pcm_corr = dmf_dotp( ks%gr%fine%mesh, ks%calc%total_density, hm%pcm%v_e_rs + hm%pcm%v_n_rs )
    endif

    if(ks%calc%calc_energy) then
      ! Get the Hartree energy
      if(.not. hm%cmplxscl%space) then
        hm%energy%hartree = M_HALF*dmf_dotp(ks%gr%fine%mesh, ks%calc%total_density, pot)
      else
        ztmp = M_HALF*zmf_dotp(ks%gr%fine%mesh,&
           ks%calc%total_density + M_zI * ks%calc%Imtotal_density, zpot, dotu = .true.)
        hm%energy%hartree   = real(ztmp)
        hm%energy%Imhartree = aimag(ztmp)
      end if
    end if

    if(ks%gr%have_fine_mesh) then
      ! we use injection to transfer to the fine grid, we cannot use
      ! restriction since the boundary conditions are not zero for the
      ! Hartree potential (and for some XC functionals).
      if(.not. hm%cmplxscl%space) then
        call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, pot, hm%vhartree, INJECTION)
      else
        aux = real(zpot)
        call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, aux, hm%vhartree, INJECTION)
        aux = aimag(zpot)
        call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, aux, hm%Imvhartree, INJECTION)
        SAFE_DEALLOCATE_P(aux)
      end if
      ! some debugging output that I will keep here for the moment, XA
      !      call dio_function_output(1, "./", "vh_fine", ks%gr%fine%mesh, pot, unit_one, is)
      !      call dio_function_output(1, "./", "vh_coarse", ks%gr%mesh, hm%vhartree, unit_one, is)
      SAFE_DEALLOCATE_P(pot)
    end if

    if (hm%cmplxscl%space) then
      SAFE_DEALLOCATE_P(zpot)
    end if
    
    POP_SUB(v_ks_hartree)
  end subroutine v_ks_hartree
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_freeze_hxc(ks)
    type(v_ks_t), intent(inout) :: ks

    PUSH_SUB(v_ks_freeze_hxc)

    ks%frozen_hxc = .true.
    
    POP_SUB(v_ks_freeze_hxc)
  end subroutine v_ks_freeze_hxc
  ! ---------------------------------------------------------

  subroutine v_ks_calculate_current(this, calc_cur)
    type(v_ks_t), intent(inout) :: this
    logical,      intent(in)    :: calc_cur

    PUSH_SUB(v_ks_calculate_current)
    
    this%calculate_current = calc_cur

    POP_SUB(v_ks_calculate_current)
  end subroutine v_ks_calculate_current
  
end module v_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
