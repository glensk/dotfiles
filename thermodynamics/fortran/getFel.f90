program getFel
    use qsort_c_module

    implicit none
    character(len=40) :: specFile, str, weightsFile
    real(8) :: dummy, nel, kB, mu, neps, nel1, step, d, eps, U, S, U0, T, T1, T2, Td, damping, f
    real(8), allocatable :: spectrum(:), w(:), Tmesh(:)
    integer i, i1, i2, n, cn, Maxi, spin, nelec, nT, j
    type(list) l(2)
    logical found

    ! parameters for mu search
    eps = 1e-8    ! eps in eV for determining the search direction
    Maxi = 50000  ! max number of steps for mu search

    kB = 1/11604.50613781094 ! Boltzmann in eV/K

    if ( iargc().NE.7 ) then
      write (*,*) "provide input:"
      write (*,*) "getFel.x T1 T2 Td nelec spin damping neps"
      call exit
    end if

    call getarg(1,str); read (str,*) T1       ! start T in K
    call getarg(2,str); read (str,*) T2       ! end T in K
    call getarg(3,str); read (str,*) Td       ! step T in K
    call getarg(4,str); read (str,*) nelec    ! number of electrons
    call getarg(5,str); read (str,*) spin     ! 1=unpolarized, 2=polarized
    call getarg(6,str); read (str,*) damping  ! damping factor for the linear minimization
    call getarg(7,str); read (str,*) neps     ! eps in nElec for stopping the mu search

    specFile = "_tmp_EIGENVALUES"
    weightsFile = "_tmp_weights"
    open (unit=10,file=specFile,status='old',err=99)
    open (unit=11,file=weightsFile,status='old',err=99)

    ! read in specFile once to determine n=number_of_eigenvalues
    n = 0
 10 read (10,*,end=20) dummy
    n = n+1
    go to 10
 20 continue
    if (n.eq.0) go to 99
    rewind (10)

    ! read in eigenvalues and weights
    ! eigenvalues in spectrum file exptected to be in eV
    allocate (spectrum(n),w(n))
    do i=1,n
      read(10,*) spectrum(i)
      read(11,*) w(i)
    end do
    close(10)
    close(11)

    ! QsortC routine produces SegFault for very large lists (somewhere about 2.000.000 to 3.000.000 entries)
    if (n>2000000) then
      write (*,*) ""
      write (*,*) "eigenvalues file contains more than 2.000.000 eigenvalues"
      write (*,*) "this is to large to be handled by quicksort routine"
      write (*,*) "EXITING"
      stop
    endif

    ! sort eigenvalues
    allocate(l(1)%A(n),l(2)%A(n))
    l(1)%A=spectrum
    l(2)%A=w
    call QsortC(l)
    spectrum=l(1)%A
    w=l(2)%A
    deallocate(l(1)%A)
    deallocate(l(2)%A)

    ! we write the files only for more than one temperature, otherwise on screen
    if (T1.NE.T2) then
      open(unit=10,file="Fermi_level_shift_meV")
      open(unit=11,file="Electronic_internal_energy_meV")
      open(unit=12,file="Electronic_entropy_kB")
      open(unit=13,file="Electronic_free_energy_meV")
    end if

    ! get U0 reference at T=0K
    i = 1
    U0 = 0.
    nel1 = 0.
    found = .false.
    do while (.NOT.found)     ! we fill up with electons until we match the required nr of electrons (nelec)
      nel = nel1 + w(i) * (3-spin)
      if ((nel-nelec)**2>(nel1-nelec)**2) then
        found = .true.
      else
        U0 = U0 + spectrum(i) * w(i)
        nel1 = nel
        i = i + 1
      end if
    end do

    ! get T mesh, add T=1K for U0 reference (0K difficult to converge)
    nT = int((T2-T1)/Td)
    allocate(Tmesh(nt+1))
    do i=1,nT+1
      Tmesh(i) = T1+(i-1)*Td
    end do

    do i=1,size(Tmesh)
      T = Tmesh(i)

      ! ========= START: search for the Fermi level: mu
      cn=0  ! counter for mu search
      mu=0.
      found=.false.
      do while (.NOT.found) ! loop for mu search

        ! calculate nel=number_of_electrons for current mu
        call FermiDiracSum(spectrum,w,n,mu,kB*T,nel)
        nel = nel * (3-spin)

        ! if neps is satisfied we found a good mu, otherwise estimate new mu search direction
        if (abs(nel-nelec)<neps) then
          found=.true.
        else
          ! calculate nel1=number_of_electrons for mu+eps and get new mu from linear extrapolation
          call FermiDiracSum(spectrum,w,n,mu+eps,kB*T,nel1)
          nel1 = nel1 * (3-spin)
          d = (nel1-nel)/eps
          step = (nelec-nel1)/d
          mu = mu + damping*step
        endif
        
        cn=cn+1
        if (cn>Maxi) then; 
          write (*,*)
          write (*,'(A,I0,A,F7.1,A)') "  cannot converge Fermi level search within ",Maxi," steps at T=",T," K";
          write(*,*);write (*,*) " goal neps: ",neps,"   current eps: ",nel-nelec
          write(*,*); write (*,*) "you can try a smaller damping with -d or a larger neps with -n"; write (*,*)
          stop
        endif
      enddo
      ! ========= END: search for the Fermi level: mu

      call InternalEnergy(spectrum,w,n, mu, kB*T, U)
      call        Entropy(spectrum,w,n, mu, kB*T, S)

      ! we write the files only for more than one temperature, otherwise on screen
      if (T1.NE.T2) then
        write (10,*) T,1000*mu                       ! Fermi level shift in meV
        write (11,*) T,(3-spin)*1000*(U-U0)          ! internal energy in meV
        write (12,*) T,(3-spin)*S                    ! entropy in kB
        write (13,*) T,(3-spin)*1000*((U-U0)-kB*T*S) ! free energy in meV
      else
        write (*,*)
        write (*,*) "Temperature in K:                 ", T
        write (*,*) "Number of electrons (given):      ", nelec
        write (*,*) "Number of electrons (integration):", nel
        write (*,*) "Converged in steps:               ", cn
        write (*,*) "Fermi level shift in meV:         ", 1000*mu
        write (*,*) "Internal energy in meV:           ", (3-spin)*1000*(U-U0)
        write (*,*) "Entropy in kB:                    ", (3-spin)*S
        write (*,*) "Free energy in meV:               ", (3-spin)*1000*((U-U0)-kB*T*S)
      end if

    end do ! T loop
    
    if (T1.NE.T2) then
      close(10)
      close(11)
      close(12)
      close(13)
      write (*,*) "files written:"
      write (*,*) "Fermi_level_shift_meV       Electronic_entropy_kB"
      write (*,*) "Electronic_free_energy_meV  Electronic_internal_energy_meV"
    end if

    f = 1/(exp((spectrum(n)-mu)/(kB*T))+1)
    write (*,*); write (*,'(A,F10.5)') " Occupancy at highest level:       ",f
    if (f>0.0001) then
      write (*,*); write (*,*) "WARNING: Occupancy is high, maybe too little NBANDS for the temperature??"
    end if

    stop

 99 continue
    write (*,*) "error: cannot open or empty file ",trim(specFile)," or ",trim(weightsFile)
    stop

end program


subroutine FermiDiracSum(spectrum,w,n,mu,kT,nel)
    implicit none
    integer i, n
    real(8) e, f, mu, kT, nel
    real(8) spectrum(n) ,w(n)

    nel = 0.
    do i=1,n
      e = spectrum(i)-mu
      f = 1/(exp(e/kT)+1)
      nel = nel + w(i) * f
    end do
end subroutine


subroutine InternalEnergy(spectrum,w,n,mu,kT,U)
    implicit none
    integer i, n
    real(8) e, f, mu, kT, U
    real(8) spectrum(n) ,w(n)

    U = 0.
    do i=1,n
      e = spectrum(i)-mu
      f = 1/(exp(e/kT)+1)
      U = U + w(i) * f * spectrum(i)   ! WARNING: in this line the UNSHIFTED eigenvalues are needed
    end do
end subroutine


subroutine Entropy(spectrum,w,n,mu,kT,S)
    implicit none
    integer i, n
    real(8) e, f, mu, kT, S
    real(8) spectrum(n) ,w(n)

    S = 0.
    do i=1,n
      e = spectrum(i)-mu
      f = 1/(exp(e/kT)+1)
      if (f.NE.1.AND.f.NE.0) S = S - w(i) * (f*log(f) + (1-f)*log(1-f))
    end do
end subroutine


