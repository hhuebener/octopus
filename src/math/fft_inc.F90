!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011 J. Alberdi-Rodriguez, P. Garcia Risueño, M. Oliveira
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
!! $Id: fft_inc.F90 12164 2014-06-02 14:50:50Z joseba $


! ---------------------------------------------------------
subroutine X(fft_forward)(fft, in, out, norm)
    type(fft_t),     intent(in)  :: fft
    R_TYPE,          intent(in)  :: in(:,:,:)
    CMPLX,           intent(out) :: out(:,:,:)
    FLOAT, optional, intent(out) :: norm 

    integer :: ii, jj, kk, slot, n1, n2, n3
    type(profile_t), save :: prof_fw
#ifdef HAVE_CLAMDFFT
    CMPLX, allocatable :: cin(:, :, :)
    type(opencl_mem_t) :: rsbuffer, fsbuffer
#endif

    PUSH_SUB(X(fft_forward))

    call profiling_in(prof_fw, "FFT_FORWARD")

    
    slot = fft%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
      ii = min(1, fft_array(slot)%rs_n(1))
      jj = min(1, fft_array(slot)%rs_n(2))
      kk = min(1, fft_array(slot)%rs_n(3))
      call fftw_execute_dft(fft_array(slot)%planf, in(ii,jj,kk), out(ii,jj,kk))
    case (FFTLIB_NFFT)
#ifdef HAVE_NFFT
      call X(nfft_forward)(fft_array(slot)%nfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%nfft%norm
#endif
    case (FFTLIB_PNFFT)
#ifdef HAVE_PNFFT
      call X(pnfft_forward)(fft_array(slot)%pnfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%pnfft%norm
#endif
    case (FFTLIB_PFFT)
      if (all(fft_array(slot)%rs_n /= 0)) then
        ASSERT(fft_array(slot)%X(rs_data)(1,1,1) == in(1,1,1))
      end if
      if (all(fft_array(slot)%fs_n /= 0)) then
        ASSERT(fft_array(slot)%fs_data(1,1,1) == out(1,1,1))
      end if
#ifdef R_TREAL
      if (fft_array(slot)%rs_n(1) > fft_array(slot)%rs_n_global(1)) then
        do kk = 1, n3
          do jj = 1, n2
            fft_array(slot)%drs_data(n1, jj, kk) = M_ZERO
          end do
        end do
      end if
#endif 
#ifdef HAVE_PFFT
      call pfft_execute(fft_array(slot)%pfft_planf)
#endif
    case(FFTLIB_CLAMD)
#ifdef HAVE_CLAMDFFT

      SAFE_ALLOCATE(cin(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)))

      cin(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)) = &
        in(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3))

      call opencl_create_buffer(rsbuffer, CL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%rs_n(1:3)))
      call opencl_create_buffer(fsbuffer, CL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%fs_n(1:3)))

      call opencl_write_buffer(rsbuffer, product(fft_array(slot)%rs_n(1:3)), cin)

      call opencl_finish()

      call clAmdFftEnqueueTransform(fft_array(slot)%cl_plan_fw, CLFFT_FORWARD, opencl%command_queue, &
        rsbuffer%mem, fsbuffer%mem, cl_status)
      if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftEnqueueTransform')

      call opencl_finish()

      call opencl_read_buffer(fsbuffer, product(fft_array(slot)%fs_n(1:3)), out)

      call opencl_finish()

      call opencl_release_buffer(rsbuffer)
      call opencl_release_buffer(fsbuffer)
#endif
    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select

    call fft_operation_count(fft)

    call profiling_out(prof_fw)

    POP_SUB(X(fft_forward))
  end subroutine X(fft_forward)

! ---------------------------------------------------------
  subroutine X(fft_forward_cl)(fft, in, out)
    type(fft_t),        intent(in)    :: fft
    type(opencl_mem_t), intent(in)    :: in
    type(opencl_mem_t), intent(inout) :: out

    integer :: slot
    type(profile_t), save :: prof_fw
#ifdef HAVE_CLAMDFFT
    type(opencl_mem_t)         :: tmp_buf
    integer                    :: bsize
    integer(8)                 :: tmp_buf_size
#endif

    PUSH_SUB(X(fft_forward_cl))

    call profiling_in(prof_fw, "FFT_FORWARD_CL")

    slot = fft%slot
    ASSERT(fft_array(slot)%library == FFTLIB_CLAMD)

#ifdef HAVE_CLAMDFFT

    call clAmdFftGetTmpBufSize(fft_array(slot)%cl_plan_bw, tmp_buf_size, cl_status)
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftGetTmpBufSize')

    if(tmp_buf_size > 0) then
      call opencl_create_buffer(tmp_buf, CL_MEM_READ_WRITE, TYPE_BYTE, int(tmp_buf_size, 4))
    end if

    if(tmp_buf_size > 0) then
      call clAmdFftEnqueueTransform(fft_array(slot)%cl_plan_fw, CLFFT_FORWARD, opencl%command_queue, &
        in%mem, out%mem, tmp_buf%mem, cl_status)
    else
      call clAmdFftEnqueueTransform(fft_array(slot)%cl_plan_fw, CLFFT_FORWARD, opencl%command_queue, &
        in%mem, out%mem, cl_status)
    end if
    
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftEnqueueTransform')
    
    call fft_operation_count(fft)

    call opencl_finish()

    if(tmp_buf_size > 0) call opencl_release_buffer(tmp_buf)

#endif

    call profiling_out(prof_fw)

    POP_SUB(X(fft_forward_cl))
  end subroutine X(fft_forward_cl)

  ! ---------------------------------------------------------

  subroutine X(fft_forward1)(fft, in, out)
    type(fft_t), intent(in)  :: fft
    R_TYPE,      intent(in)  :: in(:)
    CMPLX,       intent(out) :: out(:)

    PUSH_SUB(X(fft_forward1))

    call fftw_execute_dft(fft_array(fft%slot)%planf, in(1), out(1))
    call fft_operation_count(fft)

    POP_SUB(X(fft_forward1))
  end subroutine X(fft_forward1)

  ! ---------------------------------------------------------
  subroutine X(fft_backward)(fft, in, out, norm)
    type(fft_t), intent(in)  :: fft
    CMPLX,       intent(in)  :: in(:,:,:)
    R_TYPE,      intent(out) :: out(:,:,:)
    FLOAT, optional, intent(out) :: norm 
    
    integer :: ii, jj, kk, slot
    FLOAT :: scaling_factor
    type(profile_t), save :: prof_bw
    logical :: scale
#ifdef HAVE_CLAMDFFT
    CMPLX, allocatable :: cout(:, :, :)
    type(opencl_mem_t) :: rsbuffer, fsbuffer
#endif

    PUSH_SUB(X(fft_backward))
    
    call profiling_in(prof_bw,"FFT_BACKWARD")

    scale = .true.

    slot = fft%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
      call fftw_execute_dft(fft_array(slot)%planb, in(1,1,1), out(1,1,1))
    case (FFTLIB_NFFT)
      scale = .false. ! the result is already scaled
#ifdef HAVE_NFFT    
      call X(nfft_backward)(fft_array(slot)%nfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%nfft%norm
#endif
    case (FFTLIB_PNFFT)
      scale = .false. ! the result is already scaled
#ifdef HAVE_PNFFT    
      call X(pnfft_backward)(fft_array(slot)%pnfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%pnfft%norm
#endif
    case (FFTLIB_PFFT)
      if (all(fft_array(slot)%fs_n /= 0)) then
        ASSERT(fft_array(slot)%fs_data(1,1,1) == in(1,1,1))
      end if  
      if (all(fft_array(slot)%rs_n /= 0)) then
        ASSERT(fft_array(slot)%X(rs_data)(1,1,1) == out(1,1,1))
      end if
#ifdef HAVE_PFFT
      call pfft_execute(fft_array(slot)%pfft_planb)
#endif
    case(FFTLIB_CLAMD)
#ifdef HAVE_CLAMDFFT

      call opencl_create_buffer(rsbuffer, CL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%rs_n(1:3)))
      call opencl_create_buffer(fsbuffer, CL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%fs_n(1:3)))

      call opencl_write_buffer(fsbuffer, product(fft_array(slot)%fs_n(1:3)), in)

      call opencl_finish()

      call clAmdFftEnqueueTransform(fft_array(slot)%cl_plan_bw, CLFFT_FORWARD, opencl%command_queue, &
        fsbuffer%mem, rsbuffer%mem, cl_status)
      if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftEnqueueTransform')

      call opencl_finish()

      SAFE_ALLOCATE(cout(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)))

      call opencl_read_buffer(rsbuffer, product(fft_array(slot)%rs_n(1:3)), cout)

      call opencl_finish()

      out(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)) = &
        cout(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3))

      call opencl_release_buffer(rsbuffer)
      call opencl_release_buffer(fsbuffer)
      
      scale = .false. ! scaling is done by the library
#endif
    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select

    if(scale) then
      ! multiply by 1/(N1*N2*N2)
      scaling_factor = M_ONE/(fft_array(slot)%rs_n_global(1)*fft_array(slot)%rs_n_global(2)*fft_array(slot)%rs_n_global(3))
      !$omp parallel do
      do kk = 1, fft_array(slot)%rs_n(3)
        do jj = 1, fft_array(slot)%rs_n(2)
          do ii = 1, fft_array(slot)%rs_n(1)
            out(ii, jj, kk) = out(ii, jj, kk)*scaling_factor
          end do
        end do
      end do
      !$omp end parallel do
    end if

    call fft_operation_count(fft)

    call profiling_out(prof_bw)

    POP_SUB(X(fft_backward))
  end subroutine X(fft_backward)

  ! ---------------------------------------------------------

  subroutine X(fft_backward_cl)(fft, in, out)
    type(fft_t),        intent(in)    :: fft
    type(opencl_mem_t), intent(in)    :: in
    type(opencl_mem_t), intent(inout) :: out

    integer :: slot
    type(profile_t), save :: prof_bw
#ifdef HAVE_CLAMDFFT
    integer                    :: bsize
    integer(8)                 :: tmp_buf_size
    type(opencl_mem_t)         :: tmp_buf
#endif

    PUSH_SUB(X(fft_backward_cl))
    
    call profiling_in(prof_bw,"FFT_BACKWARD_CL")

    slot = fft%slot
    ASSERT(fft_array(slot)%library == FFTLIB_CLAMD)

#ifdef HAVE_CLAMDFFT
    
    call clAmdFftGetTmpBufSize(fft_array(slot)%cl_plan_bw, tmp_buf_size, cl_status)
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftGetTmpBufSize')

    if(tmp_buf_size > 0) then
      call opencl_create_buffer(tmp_buf, CL_MEM_READ_WRITE, TYPE_BYTE, int(tmp_buf_size, 4))
    end if

    if(tmp_buf_size > 0) then
      call clAmdFftEnqueueTransform(fft_array(slot)%cl_plan_bw, CLFFT_FORWARD, opencl%command_queue, &
        in%mem, out%mem, tmp_buf%mem, cl_status)
    else
      call clAmdFftEnqueueTransform(fft_array(slot)%cl_plan_bw, CLFFT_FORWARD, opencl%command_queue, &
        in%mem, out%mem, cl_status)
    end if
    
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftEnqueueTransform')
    
    call fft_operation_count(fft)

    call opencl_finish()

    if(tmp_buf_size > 0) call opencl_release_buffer(tmp_buf)

#endif

    call profiling_out(prof_bw)

    POP_SUB(X(fft_backward_cl))
  end subroutine X(fft_backward_cl)

  ! ---------------------------------------------------------
  subroutine X(fft_backward1)(fft, in, out)
    type(fft_t), intent(in)  :: fft
    CMPLX,       intent(in)  :: in(:)
    R_TYPE,      intent(out) :: out(:)
 
    PUSH_SUB(X(fft_backward1))
    
    call fftw_execute_dft(fft_array(fft%slot)%planb, in(1), out(1))

    call fft_operation_count(fft)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft_array(fft%slot)%rs_n_global(1), R_TOTYPE(M_ONE) / fft_array(fft%slot)%rs_n_global(1), out)

    POP_SUB(X(fft_backward1))
  end subroutine X(fft_backward1)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
