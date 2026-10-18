! This file is part of crest.
! SPDX-Identifier: LGPL-3.0-or-later

module mlip_api
  use iso_fortran_env,only:wp => real64,stderr => error_unit
  use ieee_arithmetic,only:ieee_is_finite
  use strucrd,only:coord
  use crest_parameters,only:autoaa
  use calc_type,only:calculation_settings
  use api_helpers,only:api_handle_output,api_print_e_grd
  use fmlip_relay_client,only:mlip_init,mlip_compute,mlip_finalize_all, &
  &                           MLIP_MAX_INSTANCES,MLIP_OK,GET_PID
!$ use omp_lib,only:omp_get_thread_num
  implicit none
  private

  real(wp),parameter :: hartree_to_ev = 27.211386245988_wp
  integer,parameter :: server_timeout_seconds = 120

  logical,save :: server_active(MLIP_MAX_INSTANCES) = .false.
  character(len=1024),save :: server_model(MLIP_MAX_INSTANCES) = ''

  public :: mlip_engrad
  public :: finalize_mlip

contains

  subroutine mlip_engrad(mol,calc,energy,gradient,iostatus)
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus

    integer :: thread_id,instance_id,port,spin,ierr
    integer :: pbc(3)
    real(wp) :: energy_ev
    real(wp) :: forces(3,mol%nat),cell(3,3),stress(3,3)
    logical :: exists,pr
    character(len=32) :: port_text
    character(len=:),allocatable :: server_command

    iostatus = 1
    energy = 0.0_wp
    gradient = 0.0_wp

    if (.not.allocated(calc%mlip_modelpath)) then
      write (stderr,'(a)') 'MLIP calculation has no model path'
      return
    end if
    if (.not.allocated(calc%mlip_backend)) then
      write (stderr,'(a)') 'MLIP calculation has no backend'
      return
    end if
    if (trim(calc%mlip_backend) /= 'mace') then
      write (stderr,'(a,a)') 'Unsupported MLIP backend: ',trim(calc%mlip_backend)
      return
    end if
    if (calc%chrg /= 0 .or. calc%uhf < 0 .or. calc%uhf > 1) then
      write (stderr,'(a)') 'The native MACE backend supports only neutral singlet systems'
      return
    end if
    inquire (file=trim(calc%mlip_modelpath),exist=exists)
    if (.not.exists) then
      write (stderr,'(a,a)') 'MACE model file not found: ',trim(calc%mlip_modelpath)
      return
    end if

    thread_id = 0
!$  thread_id = omp_get_thread_num()
    if (thread_id >= MLIP_MAX_INSTANCES) then
      write (stderr,'(a,i0,a,i0)') 'MLIP thread ',thread_id, &
      & ' exceeds the relay limit of ',MLIP_MAX_INSTANCES
      return
    end if
    instance_id = thread_id+1

    if (server_active(instance_id)) then
      if (trim(server_model(instance_id)) /= trim(calc%mlip_modelpath)) then
        write (stderr,'(a)') 'A CREST thread cannot switch MLIP model after server startup'
        return
      end if
    else
      port = 20000+modulo(GET_PID(),40000)+instance_id-1
      write (port_text,'(i0)') port
      server_command = 'fmlip-relay-server --port '//trim(port_text)// &
      & ' --backend mace --model '//shell_quote(trim(calc%mlip_modelpath))// &
      & ' --device cpu --dtype float64'
      if (len(server_command) > 1020) then
        write (stderr,'(a)') 'MLIP server command exceeds the relay client limit'
        return
      end if
      call mlip_init(instance_id,port,server_command,server_timeout_seconds,ierr)
      if (ierr /= MLIP_OK) then
        write (stderr,'(a,i0)') 'Failed to initialize native MACE relay, status ',ierr
        return
      end if
      server_active(instance_id) = .true.
      server_model(instance_id) = trim(calc%mlip_modelpath)
    end if

    pbc = 0
    cell = 0.0_wp
    cell(1,1) = 1.0_wp
    cell(2,2) = 1.0_wp
    cell(3,3) = 1.0_wp
    spin = max(1,calc%uhf)
    call mlip_compute(instance_id,mol%nat,mol%at,mol%xyz*autoaa,cell,pbc,0, &
    &                 calc%chrg,spin,energy_ev,forces,stress,ierr)
    if (ierr /= MLIP_OK) then
      write (stderr,'(a,i0)') 'Native MACE energy/force evaluation failed, status ',ierr
      return
    end if

    energy = energy_ev/hartree_to_ev
    gradient = -forces*(autoaa/hartree_to_ev)
    if (.not.ieee_is_finite(energy) .or. .not.all(ieee_is_finite(gradient))) then
      energy = 0.0_wp
      gradient = 0.0_wp
      write (stderr,'(a)') 'Native MACE returned a non-finite energy or gradient'
      return
    end if

    pr = .false.
    call api_handle_output(calc,'mlip.out',mol,pr)
    if (pr) call api_print_e_grd(pr,calc%prch,mol,energy,gradient)
    iostatus = 0
  end subroutine mlip_engrad

  subroutine finalize_mlip()
    integer :: ierr

    call mlip_finalize_all(ierr)
    server_active = .false.
    server_model = ''
    if (ierr /= MLIP_OK) then
      write (stderr,'(a,i0)') 'Failed to finalize native MLIP relay instances, status ',ierr
    end if
  end subroutine finalize_mlip

  function shell_quote(argument) result(quoted)
    character(len=*),intent(in) :: argument
    character(len=:),allocatable :: quoted
    integer :: i

    quoted = "'"
    do i = 1,len_trim(argument)
      if (argument(i:i) == "'") then
        quoted = quoted//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
      else
        quoted = quoted//argument(i:i)
      end if
    end do
    quoted = quoted//"'"
  end function shell_quote

end module mlip_api
