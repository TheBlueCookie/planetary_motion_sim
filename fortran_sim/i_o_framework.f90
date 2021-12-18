module io_manager
    use vectors
    implicit none

    ! unix path seperator
    character(len = 1) :: path_sep = '/'
    ! temporary string for messages
    character(len = 500) :: temp_dump

contains
    ! prints n x n matrix to terminal
    subroutine prettyprint_real(mat, n)
        real(rk), dimension(n, n), intent(in) :: mat
        integer(ik), intent(in) :: n
        integer(ik) :: i

        do i = 1, n
            write(6, *) mat(i, :)
        end do
        write(6, *) ""
    end subroutine prettyprint_real

    ! prints an n x n matrix of vectors
    subroutine prettyprint_vec_mat(vecs, n)
        type(vec), dimension(n, n), intent(in) :: vecs
        integer(ik), intent(in) :: n
        integer(ik) :: i, j

        do i = 1, n
            do j = 1, n
                write(6, *) i, j, ': (', vecs(i, j)%x, vecs(i, j)%y, vecs(i, j)%z, ')'
            end do
        end do
        write(6, *) ""
    end subroutine prettyprint_vec_mat

    ! prints a vector array of length n
    subroutine prettyprint_vec_arr(vecs, n)
        type(vec), dimension(n), intent(in) :: vecs
        integer(ik), intent(in) :: n
        integer(ik) :: i

        do i = 1, n
            write(6, *) i, ': (', vecs(i)%x, vecs(i)%y, vecs(i)%z, ')'
        end do
        write(6, *) ""
    end subroutine prettyprint_vec_arr

    ! prints a single vector
    subroutine prettyprint_vec_single(a)
        type(vec), intent(in) :: a

        write(6, *) '(', a%x, a%y, a%z, ')'
        write(6, *) ""
    end subroutine prettyprint_vec_single

    ! creates directory and exports ALL simulation results. only to be called after running the simulation
    subroutine export_all_sim_results(lpos, lvel, lacc, lforces, nbody, nstep, savestep)
        type(vec), dimension(nbody, nstep), intent(in) :: lpos, lvel, lacc, lforces
        integer(ik), intent(in) :: nbody, nstep, savestep
        integer(ik) :: i, j, t, ios
        character(len = 100) :: dirname
        character(len = 120) :: cmd, saveloc

        dirname = prepare_subdir()

        do i = 1, nbody
            write(saveloc, '(a, a, a, i0, a)') trim(dirname), path_sep, 'body_', i, '.txt'
            write(6, *) saveloc
            open(12, file = saveloc, iostat = ios, action = 'write', status = 'new')
            if (ios /= 0) then
                write(temp_dump, '(a, a, a, i0)') 'Error while handling ', saveloc, ' Errorcode: ', ios
                call error(temp_dump)
            end if
            do j = 1, nstep, savestep
                write(12, *, iostat = ios) lpos(i, j), lvel(i, j), lacc(i, j), lforces(i, j)
                if (ios /= 0) then
                    write(temp_dump, '(a, a, a, i0)') 'Error while handling ', saveloc, ' Errorcode: ', ios
                    call error(temp_dump)
                end if
            end do
            close(12)
        end do
    end subroutine export_all_sim_results

    ! creates a subdirectory of data using current epoch time as unique name
    function prepare_subdir() result(dirpath)
        character(len = 200) :: dirpath, cmd
        integer(ik) :: t

        t = time()
        write(dirpath, '(a, a, i0)') 'data', path_sep, t
        write(cmd, '(a, a)') 'mkdir ', trim(dirpath)
        call system(cmd)
        write(temp_dump, *) 'Created save directory: ', trim(dirpath)
        call info(temp_dump)
    end function prepare_subdir

    ! saves a single simulation step to file by creating a new file or appending to the existing one
    subroutine save_sim_step(last_pos, last_vel, last_acc, last_forces, epot, ekin, time, nbody, subdir_path)
        character(len = 100), intent(in) :: subdir_path
        character(len = 200) :: saveloc
        type(vec), dimension(nbody), intent(in) :: last_pos, last_vel, last_acc, last_forces
        real(rk), intent(in) :: epot, ekin, time
        integer(ik), intent(in) :: nbody
        integer(ik) :: i, j, t, ios
        logical :: fstat

        do i = 1, nbody
            write(saveloc, '(a, a, a, i0, a)') trim(subdir_path), path_sep, 'body_', i, '.txt'
            inquire(file = saveloc, exist = fstat)
            if (fstat) then
                open(1, file = saveloc, iostat = ios, action = 'write', status = 'old', position = 'append')
            else
                open(1, file = saveloc, iostat = ios, action = 'write', status = 'new')
            end if
            if (ios /= 0) then
                write(6, '(a, a, a, i0)') 'Error while handling ', saveloc, ' Errorcode: ', ios
                stop
            end if
            write(1, *) last_pos(i), last_vel(i), last_acc(i), last_forces(i), time
            close(1)
        end do

        write(saveloc, '(a, a, a)') trim(subdir_path), path_sep, 'energies.txt'
        inquire(file = saveloc, exist = fstat)
        if (fstat) then
            open(1, file = saveloc, iostat = ios, action = 'write', status = 'old', position = 'append')
        else
            open(1, file = saveloc, iostat = ios, action = 'write', status = 'new')
        end if
        if (ios /= 0) then
            write(6, '(a, a, a, i0)') 'Error while handling ', saveloc, ' Errorcode: ', ios
            stop
        end if
        write(1, *) epot, ekin, time
        close(1)

    end subroutine save_sim_step

    ! prints current simulation step info to terminal
    subroutine print_sim_info(last_pos, nbody, nstep, cstep, dt, savestep)
        type(vec), dimension(nbody), intent(in) :: last_pos
        integer(ik), intent(in) :: nbody, nstep, cstep, savestep
        real(rk), intent(in) :: dt

        write(6, '(a, i0, a, i0, a, i0)') 'SIMULATION STATUS AT STEP ', cstep, ' / ', nstep, '. Simulation time: Day ', &
                nint(cstep * dt / 86400)
        write(6, '(a, i0, a)') 'Simulating ', nbody, ' objects.'
        write(6, '(a, i0, a)') 'Saved ', cstep / savestep, ' steps to file.'
        write(6, *) 'Current object positions: '
        call prettyprint_vec_arr(last_pos, nbody)
    end subroutine print_sim_info

    ! reads parameter file and prints confirmation, number of bodies and steps have to be known before
    subroutine read_sim_params(fname, n, mass, ipos, ivel, ddur, fstep, sstep, nstep)
        integer(ik), intent(in) :: n, nstep
        integer(ik), intent(inout) :: fstep, sstep
        type(vec), dimension(n), intent(inout) :: ipos, ivel
        real(rk), dimension(n), intent(inout) :: mass
        real(rk), intent(inout) :: ddur
        real(rk), dimension(6) :: temp_arr
        character(len = 120), intent(in) :: fname
        character(len = 1000) :: line
        integer(ik) :: ios, dump1, dump2, i

        open(1, file = fname, iostat = ios, action = 'read')
        if (ios /= 0) then
            write(6, '(a, a, a, i0)') 'Error while handling ', fname, ' Errorcode: ', ios
            stop
        end if
        read(1, *, iostat = ios) dump1, dump2, ddur
        read(1, *, iostat = ios) mass
        read(1, *, iostat = ios) fstep, sstep

        if (fstep >= nstep) then
            write(temp_dump, *) 'Feedback printing step larger than total steps. No feedback will be printed!'
            call warning(temp_dump)
        end if

        if (sstep >= nstep) then
            write(temp_dump, *) 'Saving step larger than total steps. No data will be saved.'
            call warning(temp_dump)
        end if

        do i = 1, n
            read(1, *, iostat = ios) temp_arr
            if (ios < 0) then
                write(temp_dump, '(a, i0, a, i0)') 'Error while reading initial values. Expected ', n, ', found ', i - 1
                call error(temp_dump)
            end if
            ipos(i) = temp_arr(1:3)
            ivel(i) = temp_arr(4:6)
        end do

        close(1)
        write(6, '(a, a)') 'Initial parameters read from ', fname
        write(6, *) '   Number of objects: ', n
        write(6, *) '   Simulation steps: ', nstep
        write(6, *) '   Simulation duration (days): ', ddur
        write(6, *) '   Masses of objects (kg): ', mass
        write(6, *) '   Feedback and saving done every ... step: ', fstep, sstep
        write(6, *) '   Initial positions (cartesian, m):'
        call prettyprint_vec_arr(ipos, n)
        write(6, *) '   Initial velocities (cartesian, m/s): '
        call prettyprint_vec_arr(ivel, n)

    end subroutine read_sim_params

    ! reads number of bodies and steps from param file
    subroutine get_arr_dims(fname, n, m)
        character(len = 120), intent(in) :: fname
        integer(ik), intent(inout) :: n, m
        integer(ik) :: ios

        open(1, file = fname, iostat = ios, action = 'read')
        if (ios /= 0) then
            write(temp_dump, '(a, a, a, i0)') 'Error while handling ', trim(fname), ', Errorcode: ', ios
            call error(temp_dump)
        end if
        read(1, *, iostat = ios) n, m

        close(1)
    end subroutine get_arr_dims

    ! generic warning
    subroutine warning(msg)
        character(len = 500) :: msg

        write(6, *) ''
        write(6, '(a, a)') 'WARNING: ', trim(msg)
        write(6, *)
    end subroutine warning

    ! generic error
    subroutine error(msg)
        character(len = 500) :: msg

        write(6, *) ''
        write(6, '(a, a)') 'ERROR: ', trim(msg)
        write(6, *)
        stop
    end subroutine error

    ! generic info
    subroutine info(msg)
        character(len = 500) :: msg

        write(6, *) ''
        write(6, '(a, a)') 'INFO: ', trim(msg)
        write(6, *)
    end subroutine info
end module io_manager