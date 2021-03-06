!! Copyright (C) 2011 U. De Giovannini
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
!! $Id: photoelectron_spectrum.F90 12443 2014-08-28 10:36:33Z umberto $

#include "global.h"

program photoelectron_spectrum
  use command_line_m
  use datasets_m
  use geometry_m
  use global_m
  use io_binary_m
  use io_function_m
  use io_m
  use messages_m
  use parser_m
  use pes_m  
  use pes_mask_m  
  use profiling_m
  use simul_box_m
  use space_m
  use string_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m
  
  implicit none

  integer              :: ierr, mode, interp, integrate

  integer              :: dim, ll(MAX_DIM), ii, dir, how
  FLOAT                :: Emax, Emin,Estep, uEstep,uEspan(2), pol(3)
  FLOAT                :: uThstep,uThspan(2),uPhstep,uPhspan(2), pvec(3)
  FLOAT                :: center(3)
  FLOAT, pointer       :: lk(:),RR(:)
  FLOAT, allocatable   :: PESK(:,:,:)
  logical              :: interpol
  
  type(space_t)     :: space
  type(geometry_t)  :: geo
  type(simul_box_t) :: sb
  
  character(len=512) :: filename

  !Initial values
  ll = 1 
  mode = 1
  interpol = .true. 

  call global_init(is_serial = .true.)
  
  call datasets_init(1)

  call messages_init()
  
  call io_init()


  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if

  !set default values
  mode = 1
  interp = 1
  integrate = -1
  uEstep = -1
  uEspan = (/-1,-1/)
  uThstep = -1
  uThspan = (/-1,-1/)
  uPhstep = -1
  uPhspan = (/-1,-1/)
  center = (/0,0,0/)
  pvec = (/1,0,0/)
  Emin = M_ZERO
  Emax = M_ZERO
  
  
  call get_laser_polarization(pol)
  
  
  call getopt_photoelectron_spectrum(mode,interp,uEstep, uEspan,&
                                     uThstep, uThspan, uPhstep, &
                                     uPhspan, pol, center, pvec, integrate)
  if(interp  ==  0) interpol = .false.

  call pes_mask_read_info("td.general/", dim, Emax, Estep, ll(1), Lk,RR)

  write(message(1), '(a)') 'Read PES info file.'
  call messages_info(1)

  do ii=2, dim
    ll(ii) = ll(1)
  end do    
  
  SAFE_ALLOCATE(pesk(1:ll(1),1:ll(2),1:ll(3)))

  filename=io_workpath('td.general/PESM_map.obf')
  call io_binary_read(trim(filename),ll(1)**dim,pesk, ierr) 
  if(ierr > 0) then
    message(1) = "Failed to read file "//trim(filename)
    call messages_fatal(1)
  end if


  write(message(1), '(a)') 'Read PES restart file.'
  call messages_info(1)

  !! set user values
  if(uEstep >  0 .and. uEstep > Estep)    Estep = uEstep
  if(uEspan(1) > 0 ) Emin = uEspan(1)
  if(uEspan(2) > 0 ) Emax = uEspan(2)


  call unit_system_init()
 
  write(message(1),'(a,f10.2,a2,f10.2,a2,f10.2,a1)') &
                   "Zenith axis: (",pol(1),", ",pol(2),", ",pol(3),")"
  call messages_info(1)


  ! choose what to calculate
  ! these functions are defined in pes_mask_out_inc.F90
  select case(mode)
  case(1) ! Energy-resolved
    write(message(1), '(a)') 'Compute energy-resolved PES'
    call messages_info(1)
    call pes_mask_output_power_totalM(pesk,'./PES_power.sum', Lk, dim, Emax, Estep, interpol)
 
 
  case(2) ! Angle and energy resolved
    write(message(1), '(a)') 'Compute angle- and energy-resolved PES'
    call messages_info(1)
    call pes_mask_output_ar_polar_M(pesk,'./PES_angle_energy.map', Lk, dim, pol, Emax, Estep)


  case(3) ! On a plane
    
    dir = -1
    if(sum((pvec-(/1 ,0 ,0/))**2)  <= 1E-14  )  dir = 1
    if(sum((pvec-(/0 ,1 ,0/))**2)  <= 1E-14  )  dir = 2
    if(sum((pvec-(/0 ,0 ,1/))**2)  <= 1E-14  )  dir = 3

    filename = "PES_velocity.map."//index2axis(dir)//"=0"


    if (dir == -1) then
        write(message(1), '(a)') 'Unrecognized plane. Use -u to change.'
        call messages_fatal(1)
      else
        write(message(1), '(a)') 'Compute velocity map on plane: '//index2axis(dir)//" = 0"
        call messages_info(1)
    end if 
    
    if(integrate /= INTEGRATE_NONE) then
      write(message(1), '(a)') 'Integrate on: '//index2var(integrate)
      call messages_info(1)      
      filename = "PES_velocity.map.i_"//trim(index2var(integrate))//"."//index2axis(dir)//"=0"
    end if
    
    call pes_mask_output_full_mapM_cut(pesk, filename, Lk, dim, pol, dir, integrate)    

  case(4) ! Angle energy resolved on plane 
    write(message(1), '(a)') 'Compute angle and energy-resolved PES'
    call messages_info(1)
    if(uEstep >  0 .and. uEstep > Estep) then
      Estep = uEstep
    else
      Estep = Emax/size(Lk,1)
    end if

    call pes_mask_output_ar_plane_M(pesk,'./PES_energy.map', Lk, dim, pol, Emax, Estep)

  case(5) ! Angular-resolved  


    write(message(1), '(a,es19.12,a2,es19.12,2x,a19)') &
          'Compute PES on a spherical cut at E= ',Emin,", ",Emax, & 
           str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19) 
    call messages_info(1)

    if(uEstep >  0 .and. uEstep > Estep) then
     Estep = uEstep
    else
     Estep = Emax/size(Lk,1)
    end if
 
    call pes_mask_output_ar_spherical_cut_M(pesk,'./PES_sphere.map', Lk, dim, pol, Emin, Emax, Estep)       

  case(6) ! Full momentum resolved matrix 
 
    call space_init(space)
    call geometry_init(geo, space)
    call simul_box_init(sb, geo, space)
 
    call io_function_read_how(sb, how, ignore_error = .true.)
 
    write(message(1), '(a)') 'Compute full momentum-resolved PES'
    call messages_info(1)

    call pes_mask_output_full_mapM(pesk, './PES_fullmap', Lk, how, sb)        

    call simul_box_end(sb)
    call geometry_end(geo)
    call space_end(space)

  end select


  write(message(1), '(a)') 'Done'
  call messages_info(1)

  call io_end()
  call datasets_end()
  call messages_end()
  call global_end()
  
  SAFE_DEALLOCATE_A(pesk)    

  contains

    subroutine get_laser_polarization(lPol)
       FLOAT,   intent(out) :: lPol(:) 
       
        type(block_t)       :: blk
        integer             :: no_l
        
        PUSH_SUB(get_laser_polarization)
        
        no_l = 0
        if(parse_block('TDExternalFields', blk) == 0) then
          no_l = parse_block_n(blk)

          call parse_block_float(blk, 0, 1, lPol(1))
          call parse_block_float(blk, 0, 2, lPol(2))
          call parse_block_float(blk, 0, 3, lPol(3))


          call parse_block_end(blk)
        end if
        
        if(no_l > 1) then
          message(1)="There is more than one external field. Polarization will be selected"
          message(2)="from the first field. Use -V to change axis."
          call messages_info(2)
        end if

        POP_SUB(get_laser_polarization)
    end subroutine get_laser_polarization

    character(5) pure function index2var(ivar) result(ch)
      integer, intent(in) :: ivar

      select case(ivar)
        case(INTEGRATE_PHI)
          ch = 'phi'
        case(INTEGRATE_THETA)
          ch = 'theta'
        case(INTEGRATE_R)
          ch = 'r'
        case(INTEGRATE_KX)
          ch = 'kx'
        case(INTEGRATE_KY)
          ch = 'ky'
        case(INTEGRATE_KZ)
          ch = 'kz'
        case default
          write(ch,'(i1)') ivar
      end select
    end function index2var

end program photoelectron_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
