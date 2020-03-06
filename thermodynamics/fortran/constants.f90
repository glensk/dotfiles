MODULE constants
   real*8 :: Peps=1E-4        ! GPa for determining equal pressures
   real*8 :: Veps=1E-4      ! Angstrom^3/atom for finite differences
   real*8 :: VratioMin=0.97   ! for determining volume range for formation energy fitting
   real*8 :: VratioMax=1.12 
   real*8 :: maxDevFit=1.0

   integer :: Maxi=50000
   integer :: nnDef=100
   integer :: nCoefDef=4

   real*8 :: GPaTomeVAng3=6.2415097 ! GPa --> meV/angstrom^3
   real*8 :: kB=0.086173422; ! Boltzmann in meV/K
   real*8 :: pi=3.1415927
end module constants
