!! Copyright (C) 2005-2009 Heiko Appel, Florian Lorenzen, Xavier Andrade
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
!! $Id: octcl_kernel.F90 11591 2013-12-17 09:10:32Z joseba $

#include "global.h"

module octcl_kernel_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use c_pointer_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use opencl_m
  use varinfo_m

  implicit none
  private

  public ::                             &
    octcl_kernel_t,                     &
    octcl_kernel_global_init,           &
    octcl_kernel_global_end,            &
    octcl_kernel_start_call,            &
    octcl_kernel_build

#ifdef HAVE_OPENCL
  public ::                             &
    octcl_kernel_get_ref
#endif

  type octcl_kernel_t
    private
#ifdef HAVE_OPENCL
    type(cl_kernel)               :: kernel
#endif
    logical                       :: initialized = .false.
    type(octcl_kernel_t), pointer :: next
    integer                       :: arg_count
  end type octcl_kernel_t

  type(octcl_kernel_t), pointer :: head

contains

  !------------------------------------------------------------

  subroutine octcl_kernel_global_init()
    
    PUSH_SUB(octcl_kernel_global_init)

    nullify(head)

    POP_SUB(octcl_kernel_global_init)
  end subroutine octcl_kernel_global_init

  !------------------------------------------------------------
  
  subroutine octcl_kernel_global_end()
    type(octcl_kernel_t), pointer :: next_head

    PUSH_SUB(octcl_kernel_global_end)

    do
      if(.not. associated(head)) exit
      next_head => head%next
      call octcl_kernel_end(head)
      head => next_head
    end do

    POP_SUB(octcl_kernel_global_end)
  end subroutine octcl_kernel_global_end

  !------------------------------------------------------------

  subroutine octcl_kernel_build(this, file_name, kernel_name, flags)
    type(octcl_kernel_t),        intent(inout) :: this
    character(len=*),            intent(in)    :: file_name
    character(len=*),            intent(in)    :: kernel_name
    character(len=*), optional,  intent(in)    :: flags

#ifdef HAVE_OPENCL
    type(cl_program) :: prog

    PUSH_SUB(octcl_kernel_build)

    call opencl_build_program(prog, trim(conf%share)//'/opencl/'//trim(file_name), flags = flags)
    call opencl_create_kernel(this%kernel, prog, trim(kernel_name))
    call opencl_release_program(prog)
    this%initialized = .true.

    POP_SUB(octcl_kernel_build)
#endif
  end subroutine octcl_kernel_build

  !------------------------------------------------------------

  subroutine octcl_kernel_end(this)
    type(octcl_kernel_t), intent(inout) :: this
#ifdef HAVE_OPENCL
    integer :: ierr
#endif

      PUSH_SUB(octcl_kernel_end)

#ifdef HAVE_OPENCL
      call clReleaseKernel(this%kernel, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "release_kernel")
#endif
      this%initialized = .false.

      POP_SUB(octcl_kernel_end)
  end subroutine octcl_kernel_end

  !------------------------------------------------------------

  subroutine octcl_kernel_start_call(this, file_name, kernel_name, flags)
    type(octcl_kernel_t), target, intent(inout) :: this
    character(len=*),             intent(in)    :: file_name
    character(len=*),             intent(in)    :: kernel_name
    character(len=*), optional,   intent(in)    :: flags

    PUSH_SUB(octcl_kernel_start_call)

    if(.not. this%initialized) then
      call octcl_kernel_build(this, file_name, kernel_name, flags)
      this%next => head
      head => this
    end if

    POP_SUB(octcl_kernel_start_call)
  end subroutine octcl_kernel_start_call

  !--------------------------------------------------------------
#ifdef HAVE_OPENCL
  type(cl_kernel) function octcl_kernel_get_ref(this) result(ref)
    type(octcl_kernel_t), intent(in) :: this
    
    ref = this%kernel
  end function octcl_kernel_get_ref
#endif
  !--------------------------------------------------------------

end module octcl_kernel_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
