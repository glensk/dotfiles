module free_energy_surface
  type names
    character(20) name
  end type

  type fitV
    integer nV
    real*8,pointer :: fit(:)
  end type

  type fitTV
    integer nV
    real*8,pointer :: fit(:,:)
  end type

  type Fsurf
    integer nT,nF, nM
    real*8 Veq,E0
    character(20) Etype
    real*8,pointer :: T(:)
    type(fitV) :: E, cBya
    type(fitTV),pointer :: M(:)
    type(fitTV),pointer :: F(:)
  end type

  type GDsurf
    integer nT,nG
    real*8,pointer :: T(:)
    type(names),pointer :: Gnames(:)
    type(fitTV),pointer :: G(:)
  end type

  contains

  subroutine import(dat,ni,nj,fileName)
    use utils
    implicit none
    real*8,pointer :: dat(:,:)
    integer i,ni,nj
    character(*) fileName

    open(unit=11,file=fileName,status="old",err=99)
    do i=1,ni; read(11,*) dat(i,:); enddo
    close(11)
    return

 99 continue
    write (*,*) "error opening "//fileName//" for reading"; stop
  end subroutine


  subroutine expt(dat1,dat2,fileName)
    implicit none
    real*8 dat1(:),dat2(:)
    character(*) fileName
    integer i,n

    n = size(dat1)
    open(unit=11,file=fileName,status="new",err=99)
    do i=1,n; write (11,*) dat1(i),dat2(i); enddo
    close(11)
    return

 99 continue
    write(*,*) "error opening "//fileName//" for writing"; stop
  end subroutine


  subroutine export(dat1,dat2,fileName)
    real*8 dat1(:),dat2(:)
    character(*) fileName
    call expt(dat1,dat2,"_output/"//fileName)
  end subroutine


  subroutine readInp(base,strType,P,Fs,Gs)
    use utils
    implicit none
    type(Fsurf) Fs
    type(GDsurf) Gs
    character(2) str
    character(20) strType
    character(40) fileName
    character(*) base
    real*8 P
    integer i,j

    write (fileName,*) base//"input"
    open(unit=10,file=fileName,status="old",err=99)
    read(10,*) strType,P,Fs%nT,Fs%Etype,Fs%E%nV,Fs%nF,Gs%nG,Fs%nM

    allocate(Fs%E%fit(Fs%E%nV))
    read(10,*) Fs%E%fit(:)
    Fs%E0  = Fs%E%fit(1)
    Fs%Veq = Fs%E%fit(2)
    
    if (strType.EQ."hcp") then
      read(10,*) Fs%cBya%nV
      allocate(Fs%cBya%fit(Fs%cBya%nV))
      read(10,*) Fs%cBya%fit(:)
    endif

    if (Fs%nF>0) then
      allocate(Fs%T(Fs%nT))
      open(unit=12,file=base//"T",status="old",err=99)
      do i=1,Fs%nT; read(12,*) Fs%T(i); enddo
      close(12)

      allocate(Fs%F(Fs%nF))
      do j=1,Fs%nF
        read(10,*) Fs%F(j)%nV
        allocate(Fs%F(j)%fit(Fs%nT,Fs%F(j)%nV))
        write (str,'(i0)') j
        call import(Fs%F(j)%fit,Fs%nT,Fs%F(j)%nV,base//"F__"//str)
      enddo
    endif

    if (Gs%nG>0) then
      Gs%nT=Fs%nT
      allocate(Gs%T(Gs%nT))
      Gs%T=Fs%T
      allocate(Gs%G(Gs%nG),Gs%Gnames(Gs%nG))
      do j=1,Gs%nG
        read(10,*) Gs%G(j)%nV,Gs%Gnames(j)%name
        allocate(Gs%G(j)%fit(Gs%nT,Gs%G(j)%nV))
        write (str,'(I0)') j
        call import(Gs%G(j)%fit,Gs%nT,Gs%G(j)%nV,base//"G__"//str)
      enddo
    endif

    if (Fs%nM>0) then
      allocate(Fs%M(Fs%nM))
      do j=1,Fs%nM
        read(10,*) Fs%M(j)%nV
        allocate(Fs%M(j)%fit(Fs%nT,Fs%M(j)%nV))
        write (str,'(i0)') j
        call import(Fs%M(j)%fit,Fs%nT,Fs%M(j)%nV,base//"M__"//str)
      enddo
    endif

    close(10)
    return

 99 continue
    write (*,*) "error when opening "//base//"input"
    stop

  end subroutine


  subroutine cByaRatio (Fs,V,cBya)
    implicit none
    type(Fsurf) Fs
    real*8 V,cBya
    integer k

    cBya = 0.
    do k=1,Fs%cBya%nV
      cBya = cBya + Fs%cBya%fit(k) * V**(k-1.)
    enddo

  end subroutine


  subroutine EOS(Etype,V,V0,BM,BMder,E)
    implicit none
    character(20) Etype
    real*8 V,V0,BM,BMder,E

    select case (Etype)
     case ("EVinet")
       E = (4*BM*V0)/(BMder-1)**2. - 2*V0*BM*(BMder-1)**(-2.) &
           * (5+3*BMder*((V/V0)**(1./3.)-1) - 3*(V/V0)**(1./3.))*exp(-(3./2.)*(BMder-1)*((V/V0)**(1./3.)-1))
     case ("EMurn")
       E = (BM*V)/(BMder*(BMder-1))*(BMder*(1-V0/V)+(V0/V)**BMder-1)
     case ("EBirch")
       E = (9*V0*BM)/16*(((V0/V)**(2./3.)-1)**3.*BMder + ((V0/V)**(2./3.) - 1)**2.*(6 - 4*(V0/V)**(2./3.)))
     case ("ECubic")
       write (*,*) "ECubic not yet implemented"; stop
    end select
  end subroutine


  subroutine freeEnergyPerfect(Fs,V,Tind,F)
    use constants
    implicit none
    type(Fsurf) Fs
    integer Tind,k,j
    real*8 V,V0,BM,BMder,F,E,T,M,Fmag
    
    if (Tind>0) then
      ! T=0K contribution
      V0 = Fs%E%fit(2)
      BM = Fs%E%fit(3)*GPaTomeVAng3
      BMder = Fs%E%fit(4)
      call EOS(Fs%Etype,V,V0,BM,BMder,E)
      F = E
    else ! if Tind negative then do not include T=0K contribution
      F = 0.
      Tind = -Tind
    endif
    
    ! finite T contributions
    do j=1,Fs%nF
        do k=1,Fs%F(j)%nV
          F = F + Fs%F(j)%fit(Tind,k) * V**(k-1.)
        enddo
    enddo

    if (Fs%nM>0) then
      ! magnetic free energy
      T=Fs%T(Tind)
      Fmag=0
      do j=1,Fs%nM
          M = 0.
          do k=1,Fs%M(j)%nV; M = M + Fs%M(j)%fit(Tind,k) * V**(k-1.); enddo
          Fmag = Fmag - kB*T * log(M + 1)
      enddo
      F = F + Fmag
    endif

  end subroutine

  subroutine defectGibbsEnergy(Gs,P,Tind,Gdef)
    use constants
    implicit none
    type(GDsurf) Gs
    integer Tind,j,k
    real*8 Gf,Gdef,P,T

    ! formation free energy contributions
    T=Gs%T(Tind)
    Gdef=0
    do j=1,Gs%nG
        Gf = 0.
        do k=1,Gs%G(j)%nV; Gf = Gf + Gs%G(j)%fit(Tind,k) * P**(k-1.); enddo
        Gdef = Gdef - kB*T * exp(-Gf/(kB*T))
    enddo

  end subroutine


  subroutine freeEnergyDefect(Fsb,Fsd,EbulkAtoms,EdefAtoms,defAtoms,Vb,Vd,Tind,Fd)
     use constants
     implicit none
     type(Fsurf) Fsb,Fsd
     integer Tind,EbulkAtoms,EdefAtoms
     real*8, pointer :: defAtoms(:)
     integer j,k
     real*8 s,Vb,Vd,Fd,V0,BM,BMder,Ed

     ! T=0K contribution
     V0 = Fsd%E%fit(2)
     BM = Fsd%E%fit(3)*GPaTomeVAng3
     BMder = Fsd%E%fit(4)
     call EOS(Fsd%Etype,Vd,V0,BM,BMder,Ed)
     Fd = Ed

     ! finite T contributions
     do j=1,Fsd%nF
         do k=1,Fsd%F(j)%nV
           ! it is possible to supply finite T contributions for a smaller supercell
           ! than the T=0K cold curve; in such a case the 'missing' contribution
           ! is filled with the free energy of the perfect bulk (which is extensive)
           ! the s factor below takes care of this filling procedure
           ! here we do not need to take the rescaling factor s to the power of (2-k)
           ! because the parameterization is already for the larger volume that we are
           ! coming here with; we only need to rescale the free energy
           ! s F[Vlarger] // Expand
           ! a s + b s Vlarge + c s Vlarge^2 + d s Vlarge^3
           s = (EdefAtoms-defAtoms(j))/EbulkAtoms
           Fd = Fd + Fsd%F(j)%fit(Tind,k) * (Vd-s*Vb)**(k-1.) &
                 + s*Fsb%F(j)%fit(Tind,k) * Vb**(k-1.)
         enddo
     enddo
     
  end subroutine

end module free_energy_surface
