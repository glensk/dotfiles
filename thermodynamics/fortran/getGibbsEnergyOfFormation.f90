program getFormationFreeEnergy
    use constants
    use utils
    use free_energy_surface

    implicit none
    type(Fsurf) Fsb,Fsd ! Fsb: bulk, Fsd: defect
    type(GDsurf) GsDummy
    integer i,j,n,nn,nCoef,k,l,cn,nnT,nnV,nnP,o,c,LWORK
    real*8, pointer :: bulkAtoms(:),defAtoms(:)
    integer EbulkAtoms, EdefAtoms, INFO
    character*1 TRANS
    character(2) string
    character(20) dummyS
    character(10) str
    character(40) output
    real*8 dummyP,Veps2,Vmin,Vmax,Pbmin1,Pbmax1,PbminN,PbmaxN,s,Vd,Vstep,VstepOld,Vb,VbAt,vForm
    real*8 Fb1,Fb2,Fb3,Fd1,Fd2,Fd3,Pb,Pb2,Pd,Pd2,E0Form
    logical found
    real*8, allocatable :: T(:),Tlist(:),Vlist(:),Plist(:)
    real*8, allocatable :: TT(:,:,:),VV(:,:,:),PP(:,:,:),TTPV(:,:,:),VVPV(:,:,:),PPPV(:,:,:)
    real*8, pointer :: Gf(:),PV(:),Vmat(:,:),Pmat(:,:),coefsV(:,:),coefsP(:,:),GfV(:),GfP(:),WORK(:)
    real*8, pointer :: coefsPVV(:,:),coefsPVP(:,:),PVV(:),PVP(:),Vmat2(:,:),Pmat2(:,:)
    real*8 maxDev(3),fitGf,fitPV

    real*8 T1,T2,G1,G2
    real*8, allocatable :: HP(:),SP(:)

    ! read bulk and defect input
    call readInp ("_bul_",dummyS,dummyP,Fsb,GsDummy)
    call readInp ("_def_",dummyS,dummyP,Fsd,GsDummy)

    ! read additional input
    open(unit=10,file="_additional_input",status="old",err=99)
    read(10,*) nnP,nnV,nnT
    allocate(Plist(nnP),Vlist(nnV),Tlist(nnT))
    read(10,*) Plist(:); read(10,*) Vlist(:); read(10,*) Tlist(:)
    read(10,*) Vmin,Vmax,nn,nCoef
    read(10,*) EbulkAtoms,EdefAtoms
    if (Fsb%nF>0) then
      allocate(bulkAtoms(Fsb%nF),defAtoms(Fsd%nF))
      read(10,*) bulkAtoms(:)
      read(10,*) defAtoms(:)
    endif
    close(10)

    n = Fsb%nT
    allocate(T(n))
    T = Fsb%T

    ! bring each bulk contribution to same # of atoms as T=0K bulk contribution
    do j=1,Fsb%nF
      s = EbulkAtoms/bulkAtoms(j)
      do i=1,n; do k=1,Fsb%F(j)%nV
         ! the idea is as follows:
         ! we have a free energy parametrization for some small supercell:
         ! F[V_] := a + b V + c V^2 + d V^3
         ! where F and V correspond to this smaller supercell
         ! later on we are coming with a larger volume corresponding the larger supercell
         ! we therefore need to rescale V=Vlarge/s and we also need to scale the resulting free energy
         ! s F[Vlarge/s] // Expand
         ! a0 s + a1 Vlarge + (a2 Vlarge^2)/s + (a3 Vlarge^3)/s^2
         ! or:
         ! a0 s^1 + a1 Vlarge s^0 + (a2 Vlarge^2) s^(-1) + (a3 Vlarge^3) s^(-2) ...
         ! = a_(k-1) Vlarge^(k-1) s^(2-k), k=1...Vorder
         ! so the coefficients need to be taken to the power of **(2-k)
         Fsb%F(j)%fit(i,k) = Fsb%F(j)%fit(i,k) * s**(2.-k)
      enddo; enddo
    enddo

    if (Vmin==-1) Vmin = Fsb%Veq*VratioMin
    if (Vmax==-1) Vmax = Fsb%Veq*VratioMax

    Veps2 = Veps*EbulkAtoms

    ! get pressure range at lowest and highest T
    call freeEnergyPerfect(Fsb,Vmin,1,Fb1); call freeEnergyPerfect(Fsb,Vmin+Veps2,1,Fb2)
    Pbmin1 = -(Fb2-Fb1)/Veps2
    call freeEnergyPerfect(Fsb,Vmax,1,Fb1); call freeEnergyPerfect(Fsb,Vmax+Veps2,1,Fb2)
    Pbmax1 = -(Fb2-Fb1)/Veps2
    call freeEnergyPerfect(Fsb,Vmin,n,Fb1); call freeEnergyPerfect(Fsb,Vmin+Veps2,n,Fb2)
    PbminN = -(Fb2-Fb1)/Veps2
    call freeEnergyPerfect(Fsb,Vmax,n,Fb1); call freeEnergyPerfect(Fsb,Vmax+Veps2,n,Fb2)
    PbmaxN = -(Fb2-Fb1)/Veps2

    if (nn==-1) nn = nnDef
    if (nCoef==-1) nCoef = nCoefDef
    Vb=Vmin
    Vd=Fsd%Veq
    E0Form = Fsd%E0- EdefAtoms/(1.*EbulkAtoms) * Fsb%E0
    allocate(Gf(nn),GfV(nn),GfP(nn),PV(nn),PVV(nn),PVP(nn),Vmat(nn,nCoef),Pmat(nn,nCoef))
    allocate(Vmat2(nn,nCoef),Pmat2(nn,nCoef))
    allocate(coefsV(n,nCoef),coefsP(n,nCoef),coefsPVV(n,nCoef),coefsPVP(n,nCoef))
    do i=1,nnV; if (Vlist(i)==-1) Vlist(i)=Fsb%Veq; enddo
    do i=1,nnP; Plist(i)=Plist(i)*GPaTomeVAng3; enddo

    allocate(PP(nnP,n,2),VV(nnV,n,2),TT(nnT,nn,2))
    allocate(PPPV(nnP,n,2),VVPV(nnV,n,2),TTPV(nnT,nn,2))
    PP=-1;VV=-1;TT=-1
    PPPV=-1;VVPV=-1;TTPV=-1
    TRANS='N'
    LWORK=2*nCoef
    allocate(WORK(LWORK))
    maxDev(3)=0

    write (*,'(A,I0,A)') "  using ",nCoef-1,". order fit along volume"; write (*,*)
    write (*,'(A,F6.1,A,F6.1,A,F6.1,A,I4)')   "  volume range (Ang^3/at):   Vmin= ",&
    & Vmin/(1.*EbulkAtoms),"   Vmax= ",Vmax/(1.*EbulkAtoms),"   Veq= ",Fsb%Veq/(1.*EbulkAtoms),"   mesh= ",nn
    write (*,'(A,F6.1,A,F6.1,A,F4.1)') "  temperature range (K):     Tmin= ",T(1),&
    &"   Tmax= ",T(n),"   Tstep= ",T(2)-T(1)
    write (*,'(A,F6.1,A,F6.1,A,F6.1,A,I4)')   "  pressure range (GPa):      Pmin= ",&
    & Pbmin1/GPaTomeVAng3,"   Pmax= ",Pbmax1/GPaTomeVAng3,"   (T=",T(1),")"
    write (*,'(A,F6.1,A,F6.1,A,F6.1,A,I4)')   "                             Pmin= ",&
    &Pbmin1/GPaTomeVAng3,"   Pmax= ",Pbmax1/GPaTomeVAng3,"   (T=",T(n),")"
    write (*,*)
    write (6,'(A)') "  running ...  "; write (6,*) ""

    do i=1,n ! over temperature points
     do l=1,nn ! over volume points (for perfect bulk)
      Vb = Vmin + (l-1)*(Vmax-Vmin)/nn
      ! get bulk pressure
      call freeEnergyPerfect(Fsb,Vb,i,Fb1)
      call freeEnergyPerfect(Fsb,Vb+Veps2,i,Fb2)
      Pb = -(Fb2-Fb1)/Veps2
      found=.false.
      ! search along the free energy at const. T for Pb
      cn=0
      do while (.NOT.found)
        call freeEnergyDefect(Fsb,Fsd,EbulkAtoms,EdefAtoms,defAtoms,Vb,Vd,i,Fd1)
        call freeEnergyDefect(Fsb,Fsd,EbulkAtoms,EdefAtoms,defAtoms,Vb,Vd+Veps2,i,Fd2)
        Pd = -(Fd2-Fd1)/Veps2

        if (abs(Pd-Pb)<Peps) then
          found=.true.
        else
          call freeEnergyDefect(Fsb,Fsd,EbulkAtoms,EdefAtoms,defAtoms,Vb,Vd+2*Veps2,i,Fd3)
          Pd2 = -(Fd3-Fd2)/Veps2
          VstepOld = Vstep
          Vstep = (Pb-Pd)/(Pd2-Pd)*Veps2
          Vd = Vd + Vstep
        endif
        cn=cn+1
        if (cn>Maxi) then; write (*,'(A,I0,A,F5.1,A)') "cannot converge within ",Maxi," steps at T=",T(i)," K"; stop; endif
      enddo
      VbAt = Vb/(1.*EbulkAtoms)
      do c=1,nCoef; Vmat(l,c) = VbAt**(c-1.); enddo
      do c=1,nCoef; Pmat(l,c) =   Pb**(c-1.); enddo
      Vmat2=Vmat; Pmat2=Pmat
      vForm = Vd - EdefAtoms*VbAt
      Gf(l) = E0Form + (Fd1 - EdefAtoms/(1.*EbulkAtoms) * Fb1) + Pb*vForm
      GfV(l) = Gf(l); GfP(l) = Gf(l)
      PV(l) = Pb*vForm
      PVV(l) = PV(l); PVP(l) = PV(l)

      ! get Gform at constant temperature
      do o=1,nnT
        if (Tlist(o)==T(i)) then
          TT(o,l,1)=Pb/GPaTomeVAng3; TT(o,l,2)=Gf(l)
          TTPV(o,l,1)=Pb/GPaTomeVAng3; TTPV(o,l,2)=PV(l)
        endif
      enddo
     enddo ! volume points perfect bulk


     ! get V and P fit
     call DGELS(TRANS,nn,nCoef,1,Vmat,nn,GfV,nn,WORK,LWORK,INFO)
     call DGELS(TRANS,nn,nCoef,1,Pmat,nn,GfP,nn,WORK,LWORK,INFO)
     call DGELS(TRANS,nn,nCoef,1,Vmat2,nn,PVV,nn,WORK,LWORK,INFO)
     call DGELS(TRANS,nn,nCoef,1,Pmat2,nn,PVP,nn,WORK,LWORK,INFO)
     do c=1,nCoef; coefsV(i,c)=GfV(c); coefsP(i,c)=GfP(c); enddo
     do c=1,nCoef; coefsPVV(i,c)=PVV(c); coefsPVP(i,c)=PVP(c); enddo

     ! get max deviation from fit
     do l=1,nn
       fitGf = 0.
       VbAt = (Vmin + (l-1)*(Vmax-Vmin)/nn)/(1.*EbulkAtoms)
       do c=1,nCoef; fitGf=fitGf+coefsV(i,c)*VbAt**(c-1.); enddo
       if (abs(fitGf-Gf(l))>maxDev(3)) maxDev=(/ T(i),VbAt,abs(fitGf-Gf(l)) /)
     enddo

     ! get Gform and PVform at constant pressure
     do o=1,nnP
       PP(o,i,1)=T(i); PPPV(o,i,1)=T(i)
       fitGf=0.; do c=1,nCoef; fitGf=fitGf+coefsP(i,c)*Plist(o)**(c-1.); enddo
       fitPV=0.; do c=1,nCoef; fitPV=fitPV+coefsPVP(i,c)*Plist(o)**(c-1.); enddo
       PP(o,i,2)=fitGf; PPPV(o,i,2)=fitPV
     enddo
       
     ! get Gform and PVform at constant volume
     do o=1,nnV
       VV(o,i,1)=T(i); VVPV(o,i,1)=T(i)
       fitGf=0.; do c=1,nCoef; fitGf=fitGf+coefsV(i,c)*(Vlist(o)/(1.*EbulkAtoms))**(c-1.); enddo
       fitPV=0.; do c=1,nCoef; fitPV=fitPV+coefsPVV(i,c)*(Vlist(o)/(1.*EbulkAtoms))**(c-1.); enddo
       VV(o,i,2)=fitGf; VVPV(o,i,2)=fitPV
     enddo

    enddo ! temperature loop

    allocate(SP(n),HP(n))

    ! export Gform and PVform at the given pressures
    do o=1,nnP;
      if (Plist(o)/GPaTomeVAng3<1) then
        write(str,'(F6.4)') Plist(o)/GPaTomeVAng3
      else
        write(str,'(F0.4)') Plist(o)/GPaTomeVAng3
      endif
      write (output,'(A)') "output/Gform_P_"//str//"GPa"
      call removeSpaces(output)
      call expt(PP(o,:,1),PP(o,:,2),output)

      ! temperature loop to calculate S and H
      SP(1)=0; HP(1)=PP(o,1,2)
      do i=2,n
        ! o: loop over different pressures; PP: array of temperature (1) and gibbs energy of formation (2)
        T1 = PP(o,i-1,1); G1 = PP(o,i-1,2)
        T2 = PP(o,i  ,1); G2 = PP(o,i  ,2)
        SP(i) = -(G2-G1)/(T2-T1)
        HP(i) = G2 + T2*SP(i)
      enddo
      write (output,'(A)') "output/Sform_P_"//str//"GPa"
      call removeSpaces(output)
      call expt(PP(o,:,1),SP/kB,output)
      write (output,'(A)') "output/Hform_P_"//str//"GPa"
      call removeSpaces(output)
      call expt(PP(o,:,1),HP,output)

      write (output,'(A)') "output/concentration_P_"//str//"GPa"
      call removeSpaces(output)
      do i=1,n; PP(o,i,2)=exp(-PP(o,i,2)/(kB*PP(o,i,1))); enddo
      call expt(PP(o,:,1),PP(o,:,2),output)
      write (output,'(A)') "output/PVform_P_"//str//"GPa"
      call removeSpaces(output)
      call expt(PPPV(o,:,1),PPPV(o,:,2),output)
    enddo

    ! export Gform and PVform at the given volumes
    do o=1,nnV;
      write(str,'(F0.2)') Vlist(o)
      write (output,'(A)') "output/Gform_V_"//str//"Ang^3"
      call removeSpaces(output)
      call expt(VV(o,:,1),VV(o,:,2),output)
      write (output,'(A)') "output/concentration_V_"//str//"Ang^3"
      call removeSpaces(output)
      do i=1,n; VV(o,i,2)=exp(-VV(o,i,2)/(kB*VV(o,i,1))); enddo
      call expt(VV(o,:,1),VV(o,:,2),output)
      write (output,'(A)') "output/PVform_V_"//str//"Ang^3"
      call removeSpaces(output)
      call expt(VVPV(o,:,1),VVPV(o,:,2),output)
    enddo

    ! export Gform and PVform at the given temperatures
    do o=1,nnT;
      write(str,'(F0.0)') Tlist(o)
      write (output,'(A)') "output/Gform_T_"//str//"K"
      call removeSpaces(output)
      call expt(TT(o,:,1),TT(o,:,2),output)
      write (output,'(A)') "output/concentration_T_"//str//"K"
      call removeSpaces(output)
      do i=1,nn; TT(o,i,2)=exp(-TT(o,i,2)/(kB*Tlist(o))); enddo
      call expt(TT(o,:,1),TT(o,:,2),output)
      write (output,'(A)') "output/PVform_T_"//str//"K"
      call removeSpaces(output)
      call expt(TTPV(o,:,1),TTPV(o,:,2),output)
    enddo

    ! export fit
    if (EbulkAtoms==EdefAtoms+1) str='vac'
    if (EbulkAtoms==EdefAtoms+2) str='divac'
    if (EbulkAtoms==EdefAtoms-1) str='int'
    open(unit=12,file="output/Gform_"//str,status="new",err=99)
    do i=1,n
      write(12,'(F20.14)',advance='no') T(i)
      do c=1,nCoef; write(12,'(F20.14)',advance='no') coefsP(i,c); enddo
      write(12,*)
    enddo
    close(12)

    call system('rm -f _additional_input _def_input _bul_input _def_F__* _bul_F__* _def_T _bul_T')
    write (*,'(A,F5.1,A,F6.1,A,F5.2,A)') "  Done!    Max deviation of fit: ",maxDev(3),&
    &" meV at ",maxDev(1)," K and ",maxDev(2)," Ang^3/at"

    if (maxDev(3)>maxDevFit) then
      write(*,*); write(*,'(A,F0.1,A)') "  !!! WARNING:  Deviation from fit larger than ",maxDevFit," meV !!!"
    endif
    stop     
 99 continue
    write (*,*) "error opening _additional_input or output/Gform_"//str; stop
end

