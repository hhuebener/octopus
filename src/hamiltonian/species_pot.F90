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
!! $Id: species_pot.F90 12558 2014-10-15 16:07:52Z dstrubbe $

#include "global.h"

module species_pot_m
  use atom_m
  use curvilinear_m
  use datasets_m
  use double_grid_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use io_function_m
  use loct_math_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use periodic_copy_m
  use poisson_m
  use profiling_m
  use ps_m
  use root_solver_m
  use simul_box_m
  use species_m
  use splines_m
  use submesh_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                       &
    species_get_density,          &
    species_get_nlcc,             &
    dspecies_get_orbital,         &
    zspecies_get_orbital,         &
    dspecies_get_orbital_submesh, &
    zspecies_get_orbital_submesh, &
    species_get_local,            &
    species_atom_density

  type(mesh_t), pointer :: mesh_p
  FLOAT, allocatable :: rho_p(:)
  FLOAT, allocatable :: grho_p(:, :)
  FLOAT :: alpha_p
  FLOAT :: pos_p(MAX_DIM)

contains


  ! ---------------------------------------------------------
  subroutine species_atom_density(mesh, sb, atom, spin_channels, rho)
    type(mesh_t),         intent(in)    :: mesh
    type(simul_box_t),    intent(in)    :: sb
    type(atom_t), target, intent(in)    :: atom
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: rho(:, :) !< (mesh%np, spin_channels)

    integer :: isp, ip, in_points, nn, icell
    FLOAT :: rr, x, pos(1:MAX_DIM)
    FLOAT :: xx(MAX_DIM), yy(MAX_DIM), rerho, imrho
    type(species_t), pointer :: spec
    type(ps_t), pointer :: ps

#if defined(HAVE_MPI)
    integer :: in_points_red
#endif
    type(periodic_copy_t) :: pp

    PUSH_SUB(species_atom_density)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    spec => atom%spec
    rho = M_ZERO

    ! build density ...
    select case (species_type(spec))
    case (SPEC_FROM_FILE, SPEC_USDEF, SPEC_SOFT_COULOMB, &
          SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN, SPEC_PS_CPI, SPEC_PS_FHI) ! ... from userdef
      do isp = 1, spin_channels
        rho(1:mesh%np, isp) = M_ONE
        x = (species_zval(spec)/real(spin_channels, REAL_PRECISION)) / dmf_integrate(mesh, rho(:, isp))
        rho(1:mesh%np, isp) = x * rho(1:mesh%np, isp)
      end do

    case (SPEC_CHARGE_DENSITY)
      ! We put, for the electron density, the same as the positive density that 
      ! creates the external potential.

      call periodic_copy_init(pp, sb, spread(M_ZERO, dim=1, ncopies = sb%dim), &
        range = M_TWO * maxval(sb%lsize(1:sb%dim)))

      rho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy(1:sb%dim) = periodic_copy_position(pp, sb, icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = atom%x, coords = xx)
          xx(1:sb%dim) = xx(1:sb%dim) + yy(1:sb%dim)
          rr = sqrt(dot_product(xx(1:sb%dim), xx(1:sb%dim)))
          call parse_expression(rerho, imrho, sb%dim, xx, rr, M_ZERO, trim(species_rho_string(spec)))
          rho(ip, 1) = rho(ip, 1) + rerho
        end do
      end do
      call periodic_copy_end(pp)
      if(spin_channels > 1) then
        rho(:, 1) = M_HALF*rho(:, 1)
        rho(:, 2) = rho(:, 1)
      end if
      ! rescale to match the valence charge
      do isp = 1, spin_channels
        x = species_zval(spec) / dmf_integrate(mesh, rho(:, isp))
        rho(1:mesh%np, isp) = x * rho(1:mesh%np, isp)
      end do

    case (SPEC_JELLI) ! ... from jellium
      in_points = 0
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = atom%x)
        if(rr <= species_jradius(spec)) then
          in_points = in_points + 1
        end if
      end do

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(in_points, in_points_red, 1, MPI_INTEGER, MPI_SUM, mesh%vp%comm, mpi_err)
        in_points = in_points_red
      end if
#endif

      if(in_points > 0) then
        ! This probably should be done inside the mesh_function_m module.
 
        if (mesh%use_curvilinear) then
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = atom%x)
            if(rr <= species_jradius(spec)) then
              rho(ip, 1:spin_channels) = species_zval(spec) /   &
                (mesh%vol_pp(ip)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        else
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = atom%x)
            if(rr <= species_jradius(spec)) then
              rho(ip, 1:spin_channels) = species_zval(spec) /   &
                (mesh%vol_pp(1)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        end if
      end if

    case (SPEC_JELLI_SLAB) ! ... from jellium slab
      in_points = 0
      do ip = 1, mesh%np
        rr = abs( mesh%x( ip, 3 ) )
        if( rr <= species_jthick(spec)/M_TWO ) then
          in_points = in_points + 1
        end if
      end do

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(in_points, in_points_red, 1, MPI_INTEGER, MPI_SUM, mesh%vp%comm, mpi_err)
        in_points = in_points_red
      end if
#endif

      if(in_points > 0) then
        ! This probably should be done inside the mesh_function_m module.

        if (mesh%use_curvilinear) then
          do ip = 1, mesh%np
            rr = abs( mesh%x( ip, 3 ) )
            if( rr <= species_jthick(spec)/M_TWO ) then
              rho(ip, 1:spin_channels) = species_zval(spec) /   &
                (mesh%vol_pp(ip)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        else
          do ip = 1, mesh%np
            rr = abs( mesh%x( ip, 3 ) )
            if( rr <= species_jthick(spec)/M_TWO ) then
              rho(ip, 1:spin_channels) = species_zval(spec) /   &
                (mesh%vol_pp(1)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        end if
      end if

    case (SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_UPF, SPEC_PSPIO) ! ...from pseudopotential

      ! the outer loop sums densities over atoms in neighbour cells
      pos(1:MAX_DIM) = M_ZERO
      ps => species_ps(spec)

      call periodic_copy_init(pp, sb, atom%x, &
        range = spline_cutoff_radius(ps%Ur(1, 1), ps%projectors_sphere_threshold))

      do icell = 1, periodic_copy_num(pp)
        pos(1:sb%dim) = periodic_copy_position(pp, sb, icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos)
          rr = max(rr, r_small)
          nn_loop: do nn = 1, ps%conf%p
            do isp = 1, spin_channels
              if(rr >= spline_range_max(ps%Ur(nn, isp))) cycle nn_loop
            enddo
            do isp = 1, spin_channels
              rho(ip, isp) = rho(ip, isp) + ps%conf%occ(nn, isp) * spline_eval(ps%Ur(nn, isp), rr)**2 /(M_FOUR*M_PI)
            enddo
          end do nn_loop
        end do
      end do
  
      call periodic_copy_end(pp)
      nullify(ps)

    end select

    POP_SUB(species_atom_density)
  end subroutine species_atom_density
  ! ---------------------------------------------------------

  subroutine species_get_density(spec, pos, mesh, rho, Imrho)
    type(species_t),    target, intent(in)  :: spec
    FLOAT,                      intent(in)  :: pos(:)
    type(mesh_t),       target, intent(in)  :: mesh
    FLOAT,                      intent(out) :: rho(:)
    FLOAT, optional,            intent(out) :: Imrho(:)

    type(root_solver_t) :: rs
    logical :: conv
    integer :: dim
    FLOAT   :: x(1:MAX_DIM+1), chi0(MAX_DIM), startval(MAX_DIM + 1)
    FLOAT   :: delta, alpha, beta, xx(MAX_DIM), yy(MAX_DIM), rr, imrho1, rerho
    FLOAT   :: dist2, dist2_min
    integer :: icell, ipos, ip
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps
    logical :: have_point
#ifdef HAVE_MPI
    real(8) :: local_min(2), global_min(2)
#endif
    type(submesh_t)       :: sphere
    type(profile_t), save :: prof
    FLOAT,    allocatable :: rho_sphere(:)
    FLOAT, parameter      :: threshold = CNST(1e-6)
    FLOAT                 :: norm_factor
    logical               :: cmplxscl
    
    PUSH_SUB(species_get_density)

    call profiling_in(prof, "SPECIES_DENSITY")

    cmplxscl = present(Imrho)
    
    select case(species_type(spec))

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF, SPEC_PSPIO)
      ps => species_ps(spec)

      call submesh_init(sphere, mesh%sb, mesh, pos, spline_cutoff_radius(ps%nlr, threshold))
      SAFE_ALLOCATE(rho_sphere(1:sphere%np))
      
      forall(ip = 1:sphere%np) rho_sphere(ip) = sphere%x(ip, 0)
      if(sphere%np > 0) call spline_eval_vec(ps%nlr, sphere%np, rho_sphere)

      rho(1:mesh%np) = M_ZERO
      ASSERT(.not.cmplxscl) ! not implemented

      ! A small amount of charge is missing with the cutoff, we
      ! renormalize so that the long range potential is exact
      norm_factor = abs(species_zval(spec)/dsm_integrate(mesh, sphere, rho_sphere))
      
      do ip = 1, sphere%np
        rho(sphere%map(ip)) = rho(sphere%map(ip)) + norm_factor*rho_sphere(ip)
      end do

      SAFE_DEALLOCATE_A(rho_sphere)
      call submesh_end(sphere)
      nullify(ps)

    case(SPEC_FULL_DELTA)

      dist2_min = huge(delta)
      ipos = 0

      do ip = 1, mesh%np

        rho(ip) = M_ZERO

        dist2 = sum((mesh%x(ip, 1:mesh%sb%dim) - pos(1:mesh%sb%dim))**2)
        if (dist2 < dist2_min) then
          ipos = ip
          dist2_min = dist2
        end if

      end do

      write(message(1), '(3a,f5.2,3a)') &
        "Info: spec_full_delta species ", trim(species_label(spec)), &
        " atom displaced ", units_from_atomic(units_out%length, sqrt(dist2_min)), &
        " [ ", trim(units_abbrev(units_out%length)), " ]"
      call messages_info(1)

      have_point = .true.
#ifdef HAVE_MPI
      ! in parallel we have to find the minimum of the whole grid
      if(mesh%parallel_in_domains) then

        local_min = (/ dist2_min, dble(mesh%mpi_grp%rank)/)
        call MPI_Allreduce(local_min, global_min, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mesh%mpi_grp%comm, mpi_err)

        if(mesh%mpi_grp%rank /= nint(global_min(2))) have_point = .false.

      end if
#endif
      if(have_point) then
        if(mesh%use_curvilinear) then
          rho(ipos) = -species_z(spec)/mesh%vol_pp(ipos)
        else
          rho(ipos) = -species_z(spec)/mesh%vol_pp(1)
        end if
      end if

    case(SPEC_FULL_GAUSSIAN)

      ! periodic copies are not considered in this routine
      if(simul_box_is_periodic(mesh%sb)) then
        call messages_experimental("spec_full_gaussian for periodic systems")
      end if

      ! --------------------------------------------------------------
      ! Constructs density for an all-electron atom with the procedure
      ! sketched in Modine et al. [Phys. Rev. B 55, 10289 (1997)],
      ! section II.B
      ! --------------------------------------------------------------
      dim = mesh%sb%dim

      SAFE_ALLOCATE(rho_p(1:mesh%np))
      SAFE_ALLOCATE(grho_p(1:mesh%np, 1:dim+1))

      mesh_p => mesh
      pos_p = pos

      ! Initial guess.
      call curvilinear_x2chi(mesh%sb, mesh%cv, pos, chi0)
      delta   = mesh%spacing(1)
      alpha   = sqrt(M_TWO)*species_sigma(spec)*delta
      alpha_p = alpha  ! global copy of alpha
      beta    = M_ONE

      ! the first dim variables are the position of the delta function
      startval(1:dim) = CNST(1.0)

      ! the dim+1 variable is the normalization of the delta function
      startval(dim+1) = beta

      ! get a better estimate for beta
      call getrho(startval)
      beta = M_ONE / dmf_integrate(mesh, rho_p)
      startval(dim+1) = beta

      ! solve equation
      call root_solver_init(rs, dim+1, &
        solver_type=ROOT_NEWTON, maxiter=500, abs_tolerance=CNST(1.0e-10))
      call droot_solver_run(rs, func, x, conv, startval=startval)

      if(.not.conv) then
        write(message(1),'(a)') 'Internal error in species_get_density.'
        call messages_fatal(1)
      end if

      ! we want a charge of -Z
      rho = -species_z(spec)*rho_p

      nullify(mesh_p)
      SAFE_DEALLOCATE_A(grho_p)
      SAFE_DEALLOCATE_A(rho_p)


    case(SPEC_CHARGE_DENSITY)

      call periodic_copy_init(pp, mesh%sb, spread(M_ZERO, dim=1, ncopies = mesh%sb%dim), &
        range = M_TWO * maxval(mesh%sb%lsize(1:mesh%sb%dim)))

      rho = M_ZERO
      if (cmplxscl) Imrho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy(1:mesh%sb%dim) = periodic_copy_position(pp, mesh%sb, icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos, coords = xx)
          xx(1:mesh%sb%dim) = xx(1:mesh%sb%dim) + yy(1:mesh%sb%dim)
          rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
          call parse_expression(rerho, imrho1, mesh%sb%dim, xx, rr, M_ZERO, trim(species_rho_string(spec)))
          rho(ip) = rho(ip) - rerho
          if (cmplxscl) Imrho(ip) = Imrho(ip) - imrho1
        end do
      end do
      if (cmplxscl) then 
        rr = M_ONE ! XXXXX normalization
        rho(1:mesh%np) = rr * rho(1:mesh%np)
        Imrho(1:mesh%np) = rr * Imrho(1:mesh%np)
      else
        rr = species_zval(spec) / abs(dmf_integrate(mesh, rho(:)))
        rho(1:mesh%np) = rr * rho(1:mesh%np)
      end if

      call periodic_copy_end(pp)

    end select

    call profiling_out(prof)
    POP_SUB(species_get_density)
  end subroutine species_get_density


  ! ---------------------------------------------------------
  subroutine func(xin, ff, jacobian)
    FLOAT, intent(in)  :: xin(:)
    FLOAT, intent(out) :: ff(:), jacobian(:,:)

    FLOAT, allocatable :: xrho(:)
    integer :: idir, jdir, dim

    PUSH_SUB(func)

    dim = mesh_p%sb%dim

    call getrho(xin)
    SAFE_ALLOCATE(xrho(1:mesh_p%np))

    ! First, we calculate the function ff.
    do idir = 1, dim
      xrho(1:mesh_p%np) = rho_p(1:mesh_p%np) * mesh_p%x(1:mesh_p%np, idir)
      ff(idir) = dmf_integrate(mesh_p, xrho) - pos_p(idir)
    end do
    ff(dim+1) = dmf_integrate(mesh_p, rho_p) - M_ONE

    ! Now the jacobian.
    do idir = 1, dim
      do jdir = 1, dim+1
        xrho(1:mesh_p%np) = grho_p(1:mesh_p%np, jdir) * mesh_p%x(1:mesh_p%np, idir)
        jacobian(idir, jdir) = dmf_integrate(mesh_p, xrho)
      end do
    end do
    do jdir = 1, dim+1
      xrho(1:mesh_p%np) = grho_p(1:mesh_p%np, jdir)
      jacobian(dim+1, jdir) = dmf_integrate(mesh_p, xrho)
    end do

    SAFE_DEALLOCATE_A(xrho)
    POP_SUB(func)
  end subroutine func

  ! ---------------------------------------------------------
  subroutine species_get_nlcc(spec, pos, mesh, rho_core)
    type(species_t), target, intent(in)  :: spec
    FLOAT,                   intent(in)  :: pos(MAX_DIM)
    type(mesh_t),            intent(in)  :: mesh
    FLOAT,                   intent(out) :: rho_core(:)

    FLOAT :: center(MAX_DIM), rr
    integer :: icell, ip
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps

    PUSH_SUB(species_get_nlcc)

    ! only for 3D pseudopotentials, please
    if(species_is_ps(spec)) then
      ps => species_ps(spec)
      rho_core = M_ZERO
      call periodic_copy_init(pp, mesh%sb, pos, range = spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold))
      do icell = 1, periodic_copy_num(pp)
        center(1:mesh%sb%dim) = periodic_copy_position(pp, mesh%sb, icell)
        do ip = 1, mesh%np
          rr = sqrt(sum((mesh%x(ip, 1:mesh%sb%dim) - center(1:mesh%sb%dim))**2))
          if(rr < spline_range_max(ps%core)) then
            rho_core(ip) = rho_core(ip) + spline_eval(ps%core, rr)
          end if
        end do
      end do
      call periodic_copy_end(pp)
    else
      rho_core = M_ZERO
    end if

    POP_SUB(species_get_nlcc)
  end subroutine species_get_nlcc

  ! ---------------------------------------------------------
  subroutine getrho(xin)
    FLOAT, intent(in) :: xin(:)

    integer :: ip, jp, idir, dim
    FLOAT   :: r, chi(MAX_DIM)

    PUSH_SUB(getrho)

    dim = mesh_p%sb%dim
    rho_p = M_ZERO
    do ip = 1, mesh_p%np

      jp = ip
      if(mesh_p%parallel_in_domains) &
        jp = mesh_p%vp%local(mesh_p%vp%xlocal+ip-1)

      chi(1:dim) = mesh_p%idx%lxyz(jp, 1:dim) * mesh_p%spacing(1:dim) + mesh_p%sb%box_offset(1:dim) 

      r = sqrt( sum( (chi(1:dim) - xin(1:dim))**2 ) )

      if( (r/alpha_p)**2 < CNST(10.0)) then
        grho_p(ip, dim+1) = exp(-(r/alpha_p)**2)
        rho_p(ip)         = xin(dim+1) * grho_p(ip, dim+1)
      else
        grho_p(ip, dim+1) = M_ZERO
        rho_p(ip)         = M_ZERO
      end if

      do idir = 1, dim
        grho_p(ip, idir) = (M_TWO/alpha_p**2) * (chi(idir)-xin(idir)) * rho_p(ip)
      end do
    end do

    POP_SUB(getrho)
  end subroutine getrho 


  ! ---------------------------------------------------------
  subroutine species_get_local(spec, mesh, x_atom, vl, Imvl)
    type(species_t), target, intent(in)  :: spec
    type(mesh_t),            intent(in)  :: mesh
    FLOAT,                   intent(in)  :: x_atom(:)
    FLOAT,                   intent(out) :: vl(:)
    FLOAT,         optional, intent(out) :: Imvl(:) !< cmplxscl: imaginary part of the potential

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(MAX_DIM), r, r2
    integer :: ip, err, idim
    type(ps_t), pointer :: ps
    CMPLX :: zpot

    type(profile_t), save :: prof

    PUSH_SUB(species_get_local)
    call profiling_in(prof, "SPECIES_GET_LOCAL")

      select case(species_type(spec))

      case(SPEC_SOFT_COULOMB)

        do ip = 1, mesh%np
          xx(1:mesh%sb%dim) = mesh%x(ip,1:mesh%sb%dim) - x_atom(1:mesh%sb%dim)
          r2 = sum(xx(1:mesh%sb%dim)**2)
          vl(ip) = -species_zval(spec)/sqrt(r2+species_sc_alpha(spec))
        end do

      case(SPEC_USDEF)

        do ip = 1, mesh%np
          
          xx = M_ZERO
          xx(1:mesh%sb%dim) = mesh%x(ip,1:mesh%sb%dim) - x_atom(1:mesh%sb%dim)
          r = sqrt(sum(xx(1:mesh%sb%dim)**2))
          
          ! Note that as the spec%user_def is in input units, we have to convert
          ! the units back and forth
          forall(idim = 1:mesh%sb%dim) xx(idim) = units_from_atomic(units_inp%length, xx(idim))
          r = units_from_atomic(units_inp%length, r)
          zpot = species_userdef_pot(spec, mesh%sb%dim, xx, r)
          vl(ip)   = units_to_atomic(units_inp%energy, real(zpot))
          if(present(Imvl)) then!cmplxscl
            Imvl(ip) = units_to_atomic(units_inp%energy, aimag(zpot))            
          end if

        end do


      case(SPEC_FROM_FILE)

        call dio_function_input(trim(species_filename(spec)), mesh, vl, err)
        if(err /= 0) then
          write(message(1), '(a)')    'Error loading file '//trim(species_filename(spec))//'.'
          write(message(2), '(a,i4)') 'Error code returned = ', err
          call messages_fatal(2)
        end if

      case(SPEC_JELLI)
        a1 = species_z(spec)/(M_TWO*species_jradius(spec)**3)
        a2 = species_z(spec)/species_jradius(spec)
        Rb2= species_jradius(spec)**2
        
        do ip = 1, mesh%np
          
          xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim) - x_atom(1:mesh%sb%dim)
          r = sqrt(sum(xx(1:mesh%sb%dim)**2))
          
          if(r <= species_jradius(spec)) then
            vl(ip) = (a1*(r*r - Rb2) - a2)
          else
            vl(ip) = -species_z(spec)/r
          end if
          
        end do
      
      case(SPEC_JELLI_SLAB)
        a1 = M_TWO *M_PI * species_z(spec)/ (M_FOUR *mesh%sb%lsize(1) *mesh%sb%lsize(2) )

        do ip = 1, mesh%np

          r = abs( mesh%x(ip, 3 ) )

          if(r <= species_jthick(spec)/M_TWO ) then
            vl(ip) = a1 *( r*r/species_jthick(spec) + species_jthick(spec)/M_FOUR )
          else
            vl(ip) = a1 *r
          end if

        end do

      case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF, SPEC_PSPIO)
       
        ps => species_ps(spec)

        do ip = 1, mesh%np
          r2 = sum((mesh%x(ip, 1:mesh%sb%dim) - x_atom(1:mesh%sb%dim))**2)
          if(r2 < spline_range_max(ps%vlr_sq)) then
            vl(ip) = spline_eval(ps%vlr_sq, r2)
          else
            vl(ip) = P_PROTON_CHARGE*species_zval(spec)/sqrt(r2)
          end if
        end do

        nullify(ps)
        
      case(SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN, SPEC_CHARGE_DENSITY)
        vl(1:mesh%np) = M_ZERO
        
      end select

      call profiling_out(prof)
    POP_SUB(species_get_local)
  end subroutine species_get_local

#include "complex.F90"
#include "species_pot_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "species_pot_inc.F90"

end module species_pot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
