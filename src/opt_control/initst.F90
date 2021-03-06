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
!! $Id: initst.F90 12387 2014-08-11 18:58:28Z dstrubbe $

#include "global.h"

module initst_m
  use datasets_m
  use density_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use messages_m
  use opt_control_state_m
  use parser_m
  use profiling_m
  use restart_m
  use states_m
  use states_calc_m
  use states_restart_m
  use string_m
  use system_m
  use v_ks_m
  use varinfo_m
  use types_m

  implicit none

  private
  public :: initial_state_init

  integer, parameter ::  &
    oct_is_groundstate      = 1,      &
    oct_is_excited          = 2,      &
    oct_is_gstransformation = 3,      &
    oct_is_userdefined      = 4         


contains


  ! ---------------------------------------------------------
  subroutine initial_state_init(sys, hm, qcstate)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    type(opt_control_state_t), target, intent(inout) :: qcstate

    integer           :: ist, jst, ik, ib, idim, inst, inik, id, is, ip, ierr, &
                         no_states, istype, freeze_orbitals, ncols
    type(block_t)     :: blk
    type(states_t)    :: tmp_st 
    FLOAT             :: xx(MAX_DIM), rr, psi_re, psi_im
    CMPLX, allocatable :: rotation_matrix(:, :)
    type(restart_t) :: restart

    type(states_t), pointer :: psi

    PUSH_SUB(initial_state_init)

    call opt_control_state_init(qcstate, sys%st, sys%geo)
    psi => opt_control_point_qs(qcstate)
    call states_deallocate_wfns(psi)
    call states_allocate_wfns(psi, sys%gr%mesh, TYPE_CMPLX)

    !%Variable OCTInitialState
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 1
    !%Description
    !% Describes the initial state of the quantum system.
    !% Possible arguments are:
    !%Option oct_is_groundstate 1
    !% Start in the ground state.
    !%Option oct_is_excited 2
    !% Currently not in use.
    !%Option oct_is_gstransformation 3
    !% Start in a transformation of the ground-state orbitals, as defined in the
    !% block <tt>OCTInitialTransformStates</tt>.
    !%Option oct_is_userdefined 4
    !% Start in a userdefined state.
    !%End
    call parse_integer(datasets_check('OCTInitialState'), oct_is_groundstate, istype)
    if(.not.varinfo_valid_option('OCTInitialState', istype)) call input_error('OCTInitialState')    

    select case(istype)
    case(oct_is_groundstate) 
      message(1) =  'Info: Using ground state for initial state.'
      call messages_info(1)
      call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, psi%dom_st_kpt_mpi_grp, &
                         ierr, mesh=sys%gr%mesh, exact=.true.)
      if(ierr == 0) call states_load(restart, psi, sys%gr, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to read wavefunctions."
        call messages_fatal(1)
      end if
      call restart_end(restart)

    case(oct_is_excited)  
      message(1) = 'Using an excited state as the starting state for an '
      message(2) = 'optimal-control run is not possible yet.'
      message(3) = 'Try using "OCTInitialState = oct_is_transformation" instead.'
      call messages_fatal(3)

    case(oct_is_gstransformation)   
      message(1) =  'Info: Using superposition of states for initial state.'
      call messages_info(1)


      !%Variable OCTInitialTransformStates
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTInitialState = oct_is_gstransformation</tt>, you must specify an
      !% <tt>OCTInitialTransformStates</tt> block, in order to specify which linear
      !% combination of the states present in <tt>restart/gs</tt> is used to
      !% create the initial state.
      !% 
      !% The syntax is the same as the <tt>TransformStates</tt> block.
      !%End
      if(parse_isdef(datasets_check('OCTInitialTransformStates')) /= 0) then
        if(parse_block(datasets_check('OCTInitialTransformStates'), blk) == 0) then
          call states_copy(tmp_st, psi)
          call states_deallocate_wfns(tmp_st)
          call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, tmp_st%dom_st_kpt_mpi_grp, &
                            ierr, mesh=sys%gr%mesh, exact=.true.)
          if(ierr == 0) then
            call states_look_and_load(restart, tmp_st, sys%gr)
          else
            message(1) = "Could not read states for OCTInitialTransformStates."
            call messages_fatal(1)
          endif          
          call restart_end(restart)

          SAFE_ALLOCATE(rotation_matrix(1:psi%nst, 1:tmp_st%nst))
          rotation_matrix = M_z0
          do ist = 1, psi%nst
            ncols = parse_block_cols(blk, ist - 1)
            if(ncols /= tmp_st%nst) then
              write(message(1),'(a,i6)') "Wrong number of columns in OCTInitialTransformStates block, row ", ist
              call messages_fatal(1)
            endif
            do jst = 1, ncols
              call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(ist, jst))
            end do
          end do
          ! FIXME: rotation matrix should be R_TYPE
          if(states_are_real(psi)) then
            call dstates_rotate(sys%gr%mesh, psi, tmp_st, real(rotation_matrix, REAL_PRECISION))
          else
            call zstates_rotate(sys%gr%mesh, psi, tmp_st, rotation_matrix)
          endif
          SAFE_DEALLOCATE_A(rotation_matrix)
          call states_end(tmp_st)
        else
          message(1) = '"OCTInitialTransformStates" has to be specified as block.'
          call messages_info(1)
          call input_error('OCTInitialTransformStates')
        end if
      else
        message(1) = 'Error: if "OCTInitialState = oct_is_gstransformation", then you must'
        message(2) = 'supply an "OCTInitialTransformStates" block to define the transformation.'
        call messages_info(2)
        call input_error('OCTInitialTransformStates')
      end if


    case(oct_is_userdefined) 
      message(1) =  'Info: Building user-defined initial state.'
      call messages_info(1)
      
      !%Variable OCTInitialUserdefined
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% Define an initial state. Syntax follows the one of the <tt>UserDefinedStates</tt> block.
      !% Example:
      !%
      !% <tt>%OCTInitialUserdefined
      !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%  
      !%End
      if(parse_block(datasets_check('OCTInitialUserdefined'), blk) == 0) then
        
        no_states = parse_block_n(blk)
        do ib = 1, no_states
          call parse_block_integer(blk, ib - 1, 0, idim)
          call parse_block_integer(blk, ib - 1, 1, inst)
          call parse_block_integer(blk, ib - 1, 2, inik)

          ! read formula strings and convert to C strings
          do id = 1, psi%d%dim
            do is = 1, psi%nst
              do ik = 1, psi%d%nik   
                
                ! does the block entry match and is this node responsible?
                if(.not. (id  ==  idim .and. is  ==  inst .and. ik  ==  inik    &
                  .and. psi%st_start  <=  is .and. psi%st_end >= is) ) cycle
                
                ! parse formula string
                call parse_block_string(                            &
                  blk, ib - 1, 3, psi%user_def_states(id, is, ik))
                ! convert to C string
                call conv_to_C_string(psi%user_def_states(id, is, ik))
                
                do ip = 1, sys%gr%mesh%np
                  xx = sys%gr%mesh%x(ip, :)
                  rr = sqrt(sum(xx(:)**2))
                  
                  ! parse user-defined expressions
                  call parse_expression(psi_re, psi_im, &
                    sys%gr%sb%dim, xx, rr, M_ZERO, psi%user_def_states(id, is, ik))
                  ! fill state
                  psi%zpsi(ip, id, is, ik) = psi_re + M_zI * psi_im
                end do
                ! normalize orbital
                call zstates_normalize_orbital(sys%gr%mesh, psi%d%dim, &
                  psi%zpsi(:,:, is, ik))
              end do
            end do
          enddo
        end do
        call parse_block_end(blk)
      else
        message(1) = '"OCTInitialUserdefined" has to be specified as block.'
        call messages_fatal(1)
      end if
      
    case default
      write(message(1),'(a)') "No valid initial state defined."
      write(message(2),'(a)') "Choosing the ground state."
      call messages_info(2)
    end select

    ! Check whether we want to freeze some of the deeper orbitals.
    call parse_integer(datasets_check('TDFreezeOrbitals'), 0, freeze_orbitals)
    if(freeze_orbitals > 0) then
      ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
      call states_freeze_orbitals(psi, sys%gr, sys%mc, freeze_orbitals)
      write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
        ' orbitals have been frozen.', psi%nst, ' will be propagated.'
      call messages_info(1)
      call density_calc(psi, sys%gr, psi%rho)
      call v_ks_calc(sys%ks, hm, psi, sys%geo, calc_eigenval = .true.)
    elseif(freeze_orbitals < 0) then
      ! This means SAE approximation. We calculate the Hxc first, then freeze all
      ! orbitals minus one.
      write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
      call messages_info(1)
      call density_calc(psi, sys%gr, psi%rho)
      call v_ks_calc(sys%ks, hm, psi, sys%geo, calc_eigenval = .true.)
      call states_freeze_orbitals(psi, sys%gr, sys%mc, n = psi%nst - 1)
      call v_ks_freeze_hxc(sys%ks)
      call density_calc(psi, sys%gr, psi%rho)
    else
      ! Normal run.
      call density_calc(psi, sys%gr, psi%rho)
      call v_ks_calc(sys%ks, hm, psi, sys%geo, calc_eigenval = .true.)
    end if
    
    POP_SUB(initial_state_init)
  end subroutine initial_state_init

end module initst_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
