
! ----------------------------------------------------------------------------------------------------------------------------------
! general idea:
! we start with F^perf(V,T) for perfect bulk (input), transform to G^perf(P,T), and add
! G^defect(P,T) (input) to obtain complete G(P,T) from which all thermodynamics at Pin is derived
!
! in detail:
! for each T we search for the supplied pressure Pin along the Helmholtz free energy of the
! perfect bulk F^perf(V), which is provided as input parametrized in polynomials (except for E0K)
! if we find Pin at V, then we calculate 4 nearby points on F^perf, the corresponding 3 pressures
! and 3 Gibbs energies G^perf, and add on top the Gibbs energy of defects which is obtained
! from the Gibbs energies of formation G^f which are provided as parametrized input (polynomials):
!
!      Fm1 = F^perf(V-Veps)            F = F^perf(V)          Fp1 = F^perf(V+Veps)          Fp2 = F^perf(V+2*Veps)
! 
! (P=-dF/dV)     Pm1 = -(F-Fm1)/Veps              P = -(Fp1-F)/Veps            Pp1 = -(Fp2-Fp1)/Veps
!
! (G=F+PV)  Gpm1 = (F+Fm1)/2 + Pm1*(V-Veps)     Gp0 = (Fp1+F)/2 + P*V      Gpp1 = (Fp2+Fp1)/2 + Pp1*(V+Veps)
! 
!           Gdm1 = -kBT*exp[-G^f(Pm1)/kBT]    Gd0 = -kBT*exp[-G^f(P)/kBT]    Gdp1 = -kBT*exp[-G^f(Pp1)/kBT]
!
!                Gm1 = Gpm1 + Gdm1                   G = Gp0 + Gd0                 Gp1 = Gpp1 + Gdp1
!
!
! relevant thermodynamic relations (Ref. "Thermodynamics of point defects and their relation with with bulk properties", Varotsos):
!  G = G^perf+Gdefs      with  G^defs = -c^eq*kB*T  
!                        with  c^eq   = exp[-G^form/kB*T]
!                        with  G^form = G^defect-N^defect*G^perf = F^defect-N^defect*F^perfect+p*v^form
!                        with  v^form = Omega^defect-N^defect*V^perfect
!
! 3 points of the Gibbs energy at 3 different pressures (and const T) are sufficient to calculate
! up to the 2. derivative which is all we need for the derived properties, in particular:
!
!                Gm1 = Gpm1 + Gdm1                   G = Gp0 + Gd0                 Gp1 = Gpp1 + Gdp1
!
!   (V=dG/dP)                    Vm1 = (G-Gm1)/(P-Pm1)        Vp1 = (Gp1-G)/(Pp1-P)
!
!
!  derived quantities at const P (i=temperature index):
!
!                     Gibbs energy:                    GP(i)    =  G
!                     equilibrium volume:              VP(i)    =  (Vm1+Vp1)/2                                  average
!                     isothermal compress.:            kT(i)    =  -2/(Vp1+Vm1)*(Vp1-Vm1)/((Pp1-Pm1)/2)         kT=-V^-1 dV/dP
!                     isothermal bulk modulus:         BT(i)    =  1/kT(i)
!                     relative expansion:              rVP(i)   =  (VP(i)-VP(1))/VP(1)*100                      r=(V(T)-V0)/V0
!                     expansion coefficient:           Vcoef(i) =  1e5 * (VP(i)-VP(i-1))/(T(i)-T(i-1))/VP(i)    c=1/V*dV/dT
!                     lattice expansion, e.g., fcc:    aP(i)    =  (4.*VP(i))**(1./3.)
!                     relative lattice expansion:      raP(i)   =  (aP(i)-aP(1))/aP(1)*100                      r=(a(T)-a0)/a0
!                     lattice expansion coefficient:   aCoef(i) =  1e5 * (aP(i)-aP(i-1))/(T(i)-T(i-1))/aP(i)    c=1/a*da/dT
!                     entropy:                         SP(i)    =  -(GP(i)-GP(i-1))/(T(i)-T(i-1))               S=-dG/dT  (at const. pressure)
!                     enthalpy:                        HP(i)    =  GP(i) + T(i)*SP(i)                           H=G+T*S
!                     isobaric heat capacity:          CP(i)    =  T(i)*(SP(i)-SP(i-1))/(T(i)-T(i-1))           C=T*dS/dT
!
! Refs for relevant thermodynamic relations are e.g. "Thermodynamics of Crystals", Wallace or
! "Thermodynamics of point defects and their relation with with bulk properties", Varotsos
!
! main equations used:  dG =  -S dT  +  V dP  +  mu dN          (here d means total derivative)
!
!                       -->  -S = dG/dT  at constant P and N    (here d means partial derivative)
!                       -->   V = dG/dP  at constant T and N    
!                       -->  mu = dG/dN  at constant T and P
! ----------------------------------------------------------------------------------------------------------------------------------

program getThermodynamics

     use utils
     use constants
     use free_energy_surface

     implicit none

     type(Fsurf) Fs
     type(GDsurf) Gs
     integer i,j,n,cn,nAtoms
     character(2) string
     character(20) PinStr, strType
     character(40) output
     real*8 V,Vm1,V0,Vp1,Pm1,P,Pp1,PinGPa,Pin,Vstep
     real*8 Fm1,F,Fp1,Fp2,Gpm1,Gp0,Gpp1,Gdm1,Gd0,Gdp1,Gm1,G,Gp1,Gd
     logical found, file_exists

     real*8, allocatable :: T(:),aP(:),raP(:),VP(:),rVP(:),aCoef(:),Vcoef(:),GP(:),SP(:),HP(:)
     real*8, allocatable :: CV(:),CP(:),kT(:),BA(:),BT(:),gamma(:),TVb(:),FvsV(:),FvsVno0K(:),PV(:)
     real*8, allocatable :: c(:),rc(:),cCoef(:),cBya(:),volExpIn(:),FV(:),GV(:),FvsT(:)

     write (6,'(A)',advance='no') "  running ...  "
     call readInp ("_",strType,PinGPa,Fs,Gs)
     Pin=PinGPa*GPaTomeVAng3
     n = Fs%nT
     allocate(T(n),aP(n),raP(n),VP(n),rVP(n),aCoef(n),Vcoef(n),GP(n),SP(n),HP(n),FvsT(n))
     allocate(CV(n),CP(n),kT(n),BA(n),BT(n),gamma(n),TVb(n),FvsV(n),FvsVno0K(n),PV(n))
     if (strType.EQ."hcp") allocate(c(n),rc(n),cCoef(n),cBya(n))
     T = Fs%T
     V=Fs%Veq; 

     ! check if _volume_expansion file exists for exporting free energy at specific volumes
     INQUIRE(FILE="_volume_expansion",EXIST=file_exists)
     if (file_exists) then
       allocate(volExpIn(n),FV(n),GV(n))
       open(unit=11,file="_volume_expansion",status="old")
       do i=1,n; read(11,*) volExpIn(i); enddo
       close(11)
     endif

     do i=1,n  ! main loop over temperature points
       ! do first free energy at specific volume if _volume_expansion file exists
       if (file_exists) then
         call freeEnergyPerfect(Fs,volExpIn(i),i,F)
         call freeEnergyPerfect(Fs,volExpIn(i)+Veps,i,Fp1)
         P = -(Fp1-F)/Veps
         FV(i) = F + Fs%E0
         GV(i) = F + Fs%E0 + P*volExpIn(i)
       endif

       found=.false.
       ! search along the Helmholtz free energy of perfect bulk at const. T for Pin
       cn=0
       do while (.NOT.found)
         call freeEnergyPerfect(Fs,V,i,F)
         call freeEnergyPerfect(Fs,V+Veps,i,Fp1)
         P = -(Fp1-F)/Veps
         if (abs(P-Pin)<Peps) then
           found=.true.
           P = Pin
         else
           ! estimate new search volume from a second pressure
           call freeEnergyPerfect(Fs,V+2.*Veps,i,Fp2)
           Pp1 = -(Fp2-Fp1)/Veps
           if (Pp1.ne.P) then
             Vstep = (Pin-P)/(Pp1-P)*Veps
           else
             if (P>Pin) then; Vstep=Veps; else; Vstep=-Veps; endif
           endif
           V = V + Vstep
         endif
         cn=cn+1
         if (cn>Maxi) then; 
           write (*,'(A,I0,A,F7.1,A)') "cannot converge within ",Maxi," steps at T=",T(i)," K";
           write(*,*);write (*,*) " goal Peps: ",Peps,"   current Peps: ",Pin-P
           stop
         endif
       enddo

       ! get remaining free energy and pressure points
       call freeEnergyPerfect(Fs,V+2.*Veps,i,Fp2)
       Pp1 = -(Fp2-Fp1)/Veps
       call freeEnergyPerfect(Fs,V-Veps,i,Fm1)
       Pm1 = -(F-Fm1)/Veps

       ! Legendre transform to Gibbs energy of perfect bulk
       Gpm1 = (F  +Fm1)/2 + Pm1*(V-Veps)
       Gp0  = (Fp1+F  )/2 + P*V
       Gpp1 = (Fp2+Fp1)/2 + Pp1*(V+Veps)

       if (Gs%nG>0) then
         ! Gibbs energy of defects
         call defectGibbsEnergy(Gs,Pm1,i,Gdm1)
         call defectGibbsEnergy(Gs,P  ,i,Gd0 )
         call defectGibbsEnergy(Gs,Pp1,i,Gdp1)
       else
         Gdm1=0; Gd0=0; Gdp1=0
       endif

       ! total Gibbs energies
       Gm1 = Gpm1 + Gdm1;    G = Gp0 + Gd0;    Gp1 = Gpp1 + Gdp1

       ! equilibrium volumes
       Vm1 = (G-Gm1)/(P-Pm1)
       Vp1 = (Gp1-G)/(Pp1-P)

       ! main quantities at single T
       GP(i) = G+Fs%E0
       VP(i) = (Vm1+Vp1)/2
       kT(i) = -2/(Vp1+Vm1)*(Vp1-Vm1)/((Pp1-Pm1)/2)
       BT(i) = 1/kT(i)

       ! PV term
       PV(i) = P*VP(i)

       ! relative volume expansion
       rVP(i)=(VP(i)-VP(1))/VP(1)
       if (i>1) then
         ! volume expansion coefficient
         Vcoef(i) = (VP(i)-VP(i-1))/(T(i)-T(i-1))/VP(i)
       else
         Vcoef(i) = 0.
       endif

       if (strType.NE."none") then
         ! lattice constant expansion
         if (strType.EQ."fcc") aP(i)=(4.*VP(i))**(1./3.);
         if (strType.EQ."bcc") aP(i)=(2.*VP(i))**(1./3.);
         ! for hcp also c and cBya expansion
         if (strType.EQ."hcp") then
           call cByaRatio(Fs,VP(i),cBya(i))
           ! hcp conversion from volume and cBya to a and c
           nAtoms=2
           aP(i) = ((nAtoms*VP(i)/sin(pi/3.))/cBya(i))**(1./3.)
           c(i) = cBya(i)*aP(i)
         endif

         ! relative lattice constant expansion
         raP(i)=(aP(i)-aP(1))/aP(1)
         if (strType.EQ."hcp") rc(i)=(c(i)-c(1))/c(1)
         ! lattice constant coefficient
         if (i>1) then
           aCoef(i) = (aP(i)-aP(i-1))/(T(i)-T(i-1))/aP(i)
           if (strType.EQ."hcp") cCoef(i) = (c(i)-c(i-1))/(T(i)-T(i-1))/c(i)
         else
           aCoef(i) = 0.
           if (strType.EQ."hcp") cCoef(i) = 0.
         endif
       endif

       if (i>1) then
         ! entropy at const. pressure
         SP(i) = -(GP(i)-GP(i-1))/(T(i)-T(i-1))
       else
         SP(i) = 0.
       endif

       ! enthalpy
       HP(i) = GP(i) + T(i)*SP(i)

       if (i>2) then
         ! heat capacity at cont. pressure
         CP(i) = T(i)*(SP(i)-SP(i-1))/(T(i)-T(i-1))
       else
         CP(i) = 0.
       endif

       ! T * V * expCoef^2 * bulkMod and isochoric heat capacity
       TVb(i) = T(i)*VP(i)*(Vcoef(i))**2*BT(i)
       CV(i) = CP(i) - TVb(i)

       if (CV(i)>0.1) then
         ! Grueneisen parameter
         gamma(i) = VP(i)*Vcoef(i)*BT(i)/CV(i)
         ! adiabatic bulk modulus
         BA(i) = BT(i)*CP(i)/CV(i)
       else
         gamma(i) = -1
         BA(i) = BT(i)
       endif
     enddo

     ! Helmholtz free energy vs. volume at Tmax including & excluding E(T=0K)
     ! NOTE: we make a small error for systems with vacancies by adding the perfect
     !       bulk free energy to the defect formation Gibbs energy, i.e., we include
     !       the P*Vform term; the reason is that we cannot separate it at this place
     !       (we would need to do this in getGibbsEnergyOfFormation.f90) but the
     !       corresponding difference should be comparatively small and FvsV is
     !       anyhow only meant for cross-check purposes
     do i=1,n
       Gd = 0
       if (Gs%nG>0) call defectGibbsEnergy(Gs,P,n,Gd)
       call freeEnergyPerfect(Fs,VP(i),-n,F)
       FvsVno0K(i) = F + Gd
       call freeEnergyPerfect(Fs,VP(i),n,F)
       FvsV(i) = F + Gd

       ! free energy at Veq at T=0K
       call freeEnergyPerfect(Fs,VP(1),i,F)
       FvsT(i) = F + Gd
     enddo


     call export(T,VP,"volume_expansion")                        ! Ang^3
     call export(T,100*rVP,"volume_expansion_relative")          ! %
     call export(T,1e5*Vcoef,"volume_expansion_coefficient")     ! 10^-5 K^-1
     call export(T,GP,"Gibbs_energy")                            ! meV
     call export(T,HP,"enthalpy")                                ! meV
     call export(T,PV,"PV_term")                                 ! meV
     call export(T,SP/kB,"entropy")                              ! kB
     call export(T,CP/kB,"heat_capacity_isobaric")               ! kB
     call export(T,CV/kB,"heat_capacity_isochoric")              ! kB
     call export(T,TVb/kB,"TVbetaSqrBT")                         ! kB
     call export(T,kT*GPaTomeVAng3,"compressibility_isothermal") ! GPa^-1
     call export(T,BT/GPaTomeVAng3,"bulk_modulus_isothermal")    ! GPa
     call export(T,BA/GPaTomeVAng3,"bulk_modulus_adiabatic")     ! GPa
     call export(T,gamma,"Grueneisen_parameter")                 ! unitless
     call export(T,FvsT,"free_energy_atVeqT0K")                  ! meV
                                                               
     call export(VP,FvsVno0K,"free_energy_vsV_atTmax_noET0K")    ! meV
     call export(VP,FvsV,"free_energy_vsV_atTmax")               ! meV

     if (strType.NE."none") then
       call export(T,aP,"aLat_expansion")                        ! Ang
       call export(T,100*raP,"aLat_expansion_relative")          ! %
       call export(T,1e5*aCoef,"aLat_expansion_coefficient")     ! 10^-5 K-1
       if (strType.EQ."hcp") then
         call export(T,c,"cLat_expansion")                       ! Ang
         call export(T,100*rc,"cLat_expansion_relative")         ! %
         call export(T,1e5*cCoef,"cLat_expansion_coefficient")   ! 10^-5 K^-1
         call export(T,cBya,"cBya")                              ! unitless
       endif
     endif

     if (file_exists) then ! if _volume_expansion file exists
       ! the appropriate file is free_energy_at_volume_expansion
       ! because it does not contain the P*V term
       call export(T,FV,"free_energy_at_volume_expansion")       ! meV
       !call export(T,GV,"Gibbs_energy_at_volume_expansion")      ! meV
     endif

     stop     
end

