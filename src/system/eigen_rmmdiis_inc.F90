!! Copyright (C) 2009 X. Andrade
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
!! $Id: eigen_rmmdiis_inc.F90 12308 2014-07-23 20:32:09Z dstrubbe $

! ---------------------------------------------------------
!> See http://prola.aps.org/abstract/PRB/v54/i16/p11169_1
subroutine X(eigensolver_rmmdiis) (gr, st, hm, pre, tol, niter, converged, ik, diff, save_mem)
  type(grid_t),           intent(in)    :: gr
  type(states_t), target, intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  FLOAT,                  intent(out)   :: diff(:) !< (1:st%nst)
  logical,                intent(in)    :: save_mem

  R_TYPE, allocatable :: mm(:, :, :, :), evec(:, :, :), finalpsi(:)
  R_TYPE, allocatable :: eigen(:), nrmsq(:)
  FLOAT,  allocatable :: eval(:, :)
  FLOAT,  allocatable :: lambda(:)
  integer :: ist, minst, idim, ii, iter, nops, maxst, jj, bsize, ib, jter, kter, prog
  R_TYPE :: ca, cb, cc
  R_TYPE, allocatable :: fr(:, :)
  type(batch_pointer_t), allocatable :: psib(:), resb(:)
  integer, allocatable :: done(:), last(:)
  logical, allocatable :: failed(:)
  logical :: pack, save_pack_mem

  PUSH_SUB(X(eigensolver_rmmdiis))

  pack = hamiltonian_apply_packed(hm, gr%mesh)
  save_pack_mem = pack .and. save_mem

  SAFE_ALLOCATE(lambda(1:st%nst))
  SAFE_ALLOCATE(psib(1:niter))
  SAFE_ALLOCATE(resb(1:niter))
  SAFE_ALLOCATE(done(1:st%d%block_size))
  SAFE_ALLOCATE(last(1:st%d%block_size))
  SAFE_ALLOCATE(failed(1:st%d%block_size))
  SAFE_ALLOCATE(fr(1:4, 1:st%d%block_size))
  SAFE_ALLOCATE(nrmsq(1:st%d%block_size))
  SAFE_ALLOCATE(eigen(1:st%d%block_size))

  do iter = 1, niter
    if(iter /= 1) then
      SAFE_ALLOCATE(psib(iter)%batch)
    end if
    SAFE_ALLOCATE(resb(iter)%batch)
  end do

  nops = 0

  call profiling_in(prof, "RMMDIIS")

  failed = .false.
  prog = 0

  do ib = st%group%block_start, st%group%block_end
    minst = states_block_min(st, ib)
    maxst = states_block_max(st, ib)
    bsize = maxst - minst + 1

    psib(1)%batch => st%group%psib(ib, ik)

    if(pack) call batch_pack(psib(1)%batch)

    call batch_copy(psib(1)%batch, resb(1)%batch, reference = .false.)

    call X(hamiltonian_apply_batch)(hm, gr%der, psib(1)%batch, resb(1)%batch, ik)
    nops = nops + bsize

    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1)%batch, resb(1)%batch, eigen)

    st%eigenval(minst:maxst, ik) = R_REAL(eigen(1:bsize))

    call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(1)%batch, resb(1)%batch)

    done = 0

    call X(mesh_batch_dotp_vector)(gr%mesh, resb(1)%batch, resb(1)%batch, nrmsq)

    do ii = 1, bsize
      if(sqrt(abs(nrmsq(ii))) < tol) done(ii) = 1
    end do

    if(all(done(1:bsize) /= 0)) then
      if(pack) then
        call batch_unpack(st%group%psib(ib, ik))
      end if
      call batch_end(resb(1)%batch, copy = .false.)
      cycle
    end if

    call batch_copy(psib(1)%batch, psib(2)%batch, reference = .false.)

    ! get lambda 
    call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(1)%batch, psib(2)%batch)

    if(save_pack_mem) call batch_unpack(resb(1)%batch)

    call batch_copy(psib(1)%batch, resb(2)%batch, reference = .false.)

    call X(hamiltonian_apply_batch)(hm, gr%der, psib(2)%batch, resb(2)%batch, ik)
    nops = nops + bsize

    call X(mesh_batch_dotp_vector)(gr%mesh, psib(2)%batch, psib(2)%batch, fr(1, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1)%batch, psib(2)%batch, fr(2, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(gr%mesh, psib(2)%batch, resb(2)%batch, fr(3, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1)%batch, resb(2)%batch, fr(4, :), reduce = .false.)

    if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, fr)

    do ist = minst, maxst
      ii = ist - minst + 1

      ca = fr(1, ii)*fr(4, ii) - fr(3, ii)*fr(2, ii)
      cb = fr(3, ii) - st%eigenval(ist, ik)*fr(1, ii)
      cc = st%eigenval(ist, ik)*fr(2, ii) - fr(4, ii)

      lambda(ist) = R_REAL(CNST(2.0)*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc)))

      ! restrict the value of lambda to be between 0.1 and 1.0
      if(abs(lambda(ist)) > CNST(1.0)) lambda(ist) = lambda(ist)/abs(lambda(ist))
      if(abs(lambda(ist)) < CNST(0.1)) lambda(ist) = CNST(0.1)*lambda(ist)/abs(lambda(ist))
    end do

    SAFE_ALLOCATE(mm(1:niter, 1:niter, 1:2, 1:bsize))

    do iter = 2, niter

      ! for iter == 2 the preconditioning was done already
      if(iter > 2) then
        if(save_pack_mem) call batch_unpack(psib(iter - 2)%batch)
        call batch_copy(psib(iter - 1)%batch, psib(iter)%batch, reference = .false.)
        call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(iter - 1)%batch, psib(iter)%batch)
      end if

      ! predict by jacobi
      call batch_xpay(gr%mesh%np, psib(iter - 1)%batch, lambda, psib(iter)%batch)

      if(iter > 2) then
        if(save_pack_mem) call batch_unpack(resb(iter - 1)%batch)

        call batch_copy(psib(iter)%batch, resb(iter)%batch, reference = .false.)
      end if

      ! calculate the residual
      call X(hamiltonian_apply_batch)(hm, gr%der, psib(iter)%batch, resb(iter)%batch, ik)
      nops = nops + bsize

      call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(iter)%batch, resb(iter)%batch)

      call profiling_in(prof_sync, "RMMDIIS_SYNC")
      if(save_pack_mem) then
        call batch_sync(psib(iter)%batch)
        call batch_sync(resb(iter)%batch)
      end if
      call profiling_out(prof_sync)

      call profiling_in(prof_iter, 'RMMDIIS_MATRIX')
      ! calculate the matrix elements between iterations
      do jter = 1, iter
        do kter = 1, jter
          
          if(jter < iter - 1 .and. kter < iter - 1) then
            ! it was calculated on the previous iteration
            ! in parallel this was already reduced, so we set it to zero in non-root ranks
            if(gr%mesh%parallel_in_domains .and. gr%mesh%mpi_grp%rank /= 0) mm(jter, kter, 1:2, 1:bsize) = CNST(0.0)
            cycle
          endif

          call X(mesh_batch_dotp_vector)(gr%mesh, resb(jter)%batch, resb(kter)%batch, mm(jter, kter, 1, :), reduce = .false.)
          call X(mesh_batch_dotp_vector)(gr%mesh, psib(jter)%batch, psib(kter)%batch, mm(jter, kter, 2, :), reduce = .false.)

        end do
      end do
      call profiling_out(prof_iter)
            
      ! symmetrize
      do jter = 1, iter
        do kter = jter + 1, iter
          mm(jter, kter, 1:2, 1:bsize) = mm(kter, jter, 1:2, 1:bsize)
        end do
      end do

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, mm)

      SAFE_ALLOCATE(evec(1:iter, 1:1, 1:bsize))
      SAFE_ALLOCATE(eval(1:iter, 1:bsize))

      do ist = minst, maxst
        ii = ist - minst + 1

        failed(ii) = .false.
        call lalg_lowest_geneigensolve(1, iter, mm(:, :, 1, ii), mm(:, :, 2, ii), eval(:, ii), evec(:, :, ii), bof = failed(ii))
        if(failed(ii)) then
          last(ii) = iter - 1
          ! we shouldn`t do anything
          evec(1:iter - 1, 1, ii) = CNST(0.0)
          evec(iter, 1, ii) = CNST(1.0)
          cycle
        end if
      end do

      call batch_end(resb(iter)%batch, copy = .false.)

      call profiling_in(prof_lc, "RMMDIIS_LC")

      call batch_scal(gr%mesh%np, evec(iter, 1, :), psib(iter)%batch, a_start = minst)

      do jj = 1, iter - 1
        if(pack) call batch_pack(psib(jj)%batch)
        call batch_axpy(gr%mesh%np, evec(jj, 1, :), psib(jj)%batch, psib(iter)%batch, a_start = minst)
        if(pack) call batch_unpack(psib(jj)%batch, copy = .false.)
      end do

      call profiling_out(prof_lc)

      call batch_copy(psib(iter)%batch, resb(iter)%batch, reference = .false.)

      ! re-calculate the residual
      call X(hamiltonian_apply_batch)(hm, gr%der, psib(iter)%batch, resb(iter)%batch, ik)
      nops = nops + bsize
      call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(iter)%batch, resb(iter)%batch)

      ! why not allocate these outside the loop?
      SAFE_DEALLOCATE_A(eval)      
      SAFE_DEALLOCATE_A(evec)

      if(in_debug_mode) then
        call X(mesh_batch_dotp_vector)(gr%der%mesh, resb(iter)%batch, resb(iter)%batch, eigen)

        do ist = minst, maxst
          write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') &
            'Debug: RMMDIIS Eigensolver - ik', ik, ' ist ', ist, ' iter ', iter, ' res ', sqrt(abs(eigen(ist - minst + 1)))
          call messages_info(1)
        enddo
      end if

    end do ! iter

    if(save_pack_mem) call batch_unpack(psib(niter - 1)%batch, copy = any(failed(1:bsize)))

    SAFE_DEALLOCATE_A(mm)

    if(save_pack_mem) then
      call batch_pack(resb(niter - 1)%batch, copy = .false.)
    end if

    ! end with a trial move
    call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(niter)%batch, resb(niter - 1)%batch)

    call batch_xpay(gr%mesh%np, psib(niter)%batch, lambda, resb(niter - 1)%batch)

    if(any(failed(1:bsize))) then 
      SAFE_ALLOCATE(finalpsi(1:gr%mesh%np))

      do ist = minst, maxst
        ii = ist - minst + 1
        
        if(failed(ii)) then
          do idim = 1, st%d%dim
            call batch_get_state(psib(last(ii))%batch, (/ist, idim/), gr%mesh%np, finalpsi)
            call batch_set_state(resb(niter - 1)%batch, (/ist, idim/), gr%mesh%np, finalpsi)
          end do
        end if
      end do

      SAFE_DEALLOCATE_A(finalpsi)
    end if

    ! we can remove most of the batches
    do iter = 1, niter
      if(iter /= 1) call batch_end(psib(iter)%batch, copy = .false.)
      if(iter /= niter - 1) call batch_end(resb(iter)%batch, copy = .false.)
    end do

    if(save_pack_mem) call batch_pack(psib(1)%batch, copy = .false.)

    call batch_copy_data(gr%mesh%np, resb(niter - 1)%batch, st%group%psib(ib, ik))

    call batch_end(resb(niter - 1)%batch, copy = .false.)

    if(pack) call batch_unpack(st%group%psib(ib, ik))

    prog = prog + bsize
    if(mpi_grp_is_root(mpi_world) .and. .not. in_debug_mode) then
      call loct_progress_bar(st%nst*(ik - 1) + prog, st%nst*st%d%nik)
    end if
    
  end do ! ib

  call profiling_out(prof)

  call X(states_orthogonalization_full)(st, gr%mesh, ik)

  ! recalculate the eigenvalues and residuals
  do ib = st%group%block_start, st%group%block_end
    minst = states_block_min(st, ib)
    maxst = states_block_max(st, ib)
    
    call batch_copy(st%group%psib(ib, ik), resb(1)%batch, reference = .false.)
    
    call X(hamiltonian_apply_batch)(hm, gr%der, st%group%psib(ib, ik), resb(1)%batch, ik)
    call X(mesh_batch_dotp_vector)(gr%der%mesh, st%group%psib(ib, ik), resb(1)%batch, eigen)
    
    st%eigenval(minst:maxst, ik) = R_REAL(eigen(1:maxst - minst + 1))
    
    call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), st%group%psib(ib, ik), resb(1)%batch)
    
    call X(mesh_batch_dotp_vector)(gr%der%mesh, resb(1)%batch, resb(1)%batch, eigen)
    
    diff(minst:maxst) = sqrt(abs(eigen(1:maxst - minst + 1)))
    
    call batch_end(resb(1)%batch, copy = .false.)
    
    nops = nops + maxst - minst + 1
  end do
  
  converged = converged + count(diff(:) <= tol)

  do iter = 1, niter
    if(iter /= 1) then
      SAFE_DEALLOCATE_P(psib(iter)%batch)
    end if
    SAFE_DEALLOCATE_P(resb(iter)%batch)
  end do

  niter = nops

  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(psib)
  SAFE_DEALLOCATE_A(resb)
  SAFE_DEALLOCATE_A(done)
  SAFE_DEALLOCATE_A(last)
  SAFE_DEALLOCATE_A(failed)
  SAFE_DEALLOCATE_A(fr)
  SAFE_DEALLOCATE_A(nrmsq)

  POP_SUB(X(eigensolver_rmmdiis))

end subroutine X(eigensolver_rmmdiis)

! ---------------------------------------------------------

subroutine X(eigensolver_rmmdiis_min) (gr, st, hm, pre, niter, converged, ik)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik

  integer, parameter :: sweeps = 5
  integer, parameter :: sd_steps = 2

  integer :: isd, ist, minst, maxst, ib, ii
  R_TYPE  :: ca, cb, cc
  R_TYPE, allocatable :: lambda(:), diff(:)
  R_TYPE, allocatable :: me1(:, :), me2(:, :)
  logical :: pack
  type(batch_t) :: resb, kresb

  PUSH_SUB(X(eigensolver_rmmdiis_min))

  pack = hamiltonian_apply_packed(hm, gr%mesh)

  SAFE_ALLOCATE(me1(1:2, 1:st%d%block_size))
  SAFE_ALLOCATE(me2(1:4, 1:st%d%block_size))
  SAFE_ALLOCATE(lambda(1:st%nst))

  niter = 0

  if(in_debug_mode) then
    SAFE_ALLOCATE(diff(1:st%d%block_size))
  endif

  do ib = st%group%block_start, st%group%block_end
    minst = states_block_min(st, ib)
    maxst = states_block_max(st, ib)

    if(pack) call batch_pack(st%group%psib(ib, ik))

    call batch_copy(st%group%psib(ib, ik), resb, reference = .false.)
    call batch_copy(st%group%psib(ib, ik), kresb, reference = .false.)

    do isd = 1, sd_steps

      call X(hamiltonian_apply_batch)(hm, gr%der, st%group%psib(ib, ik), resb, ik)

      call X(mesh_batch_dotp_vector)(gr%mesh, st%group%psib(ib, ik), resb, me1(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%group%psib(ib, ik), st%group%psib(ib, ik), me1(2, :), reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, me1)

      forall(ist = minst:maxst) st%eigenval(ist, ik) = R_REAL(me1(1, ist - minst + 1)/me1(2, ist - minst + 1))
 
      call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), st%group%psib(ib, ik), resb)

      if(in_debug_mode) then
        call X(mesh_batch_dotp_vector)(gr%der%mesh, resb, resb, diff)

        do ist = minst, maxst
          write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') &
            'Debug: RMMDIIS MIN Eigensolver - ik', ik, ' ist ', ist, ' iter ', isd, ' res ', sqrt(abs(diff(ist - minst + 1)))
          call messages_info(1)
        enddo
      end if

      call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb, kresb)

      call X(hamiltonian_apply_batch)(hm, gr%der, kresb, resb, ik)

      niter = niter + 2*(maxst - minst + 1)

      call X(mesh_batch_dotp_vector)(gr%mesh, kresb, kresb, me2(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%group%psib(ib, ik),  kresb, me2(2, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, kresb, resb,  me2(3, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%group%psib(ib, ik),  resb,  me2(4, :), reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, me2)

      do ist = minst, maxst
        ii = ist - minst + 1

        ca = me2(1, ii)*me2(4, ii) - me2(3, ii)*me2(2, ii)
        cb = me1(2, ii)*me2(3, ii) - me1(1, ii)*me2(1, ii)
        cc = me1(1, ii)*me2(2, ii) - me1(2, ii)*me2(4, ii)

        lambda(ist) = CNST(2.0)*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

      end do

      call batch_axpy(gr%mesh%np, lambda, kresb, st%group%psib(ib, ik))
      
    end do

    if(pack) call batch_unpack(st%group%psib(ib, ik))

    call batch_end(resb, copy = .false.)
    call batch_end(kresb, copy = .false.)

    if(mpi_grp_is_root(mpi_world) .and. .not. in_debug_mode) then
      call loct_progress_bar(st%nst*(ik - 1) +  maxst, st%nst*st%d%nik)
    end if

  end do

  if(in_debug_mode) then
    SAFE_DEALLOCATE_A(diff)
  endif

  call X(states_orthogonalization_full)(st, gr%mesh, ik)

  converged = 0

  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(me1)
  SAFE_DEALLOCATE_A(me2)

  POP_SUB(X(eigensolver_rmmdiis_min))

end subroutine X(eigensolver_rmmdiis_min)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
