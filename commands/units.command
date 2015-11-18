cat reg_white_dots.dat | awk '{print 1357/$1,exp(-$2/($1*0.086173422))}'  ## Temp(K) and Gform(meV) --> conz and Tm/Temp

cal / g / K  to >> J/K/mol     1cal = 4.1868Joule; mg = 24.305u = 24.305g/mol >> 1J/K/mol=1cal/g/K * 4.1868 * 24.305

cu = 63.546u >> 1J/K/mol=1cal/g/K * 4.1868 * 63.546

pt = 195.084u >> 1J/K/mol=1cal/g/K * 4.1868 * 195.084; 
You have: 4.1868 195.084 joule / mol K N_A
You want: k
    * 98.235666
    *   / 0.010179602
--> also die tabellendaten in cal/g/K * 98.235666



Silver (Ag):
    You have:
    You have: 4.1868 * 107.8682 joule / mol K N_A
    You want: k
        * 54.31765
    *   / 0.018410222
--> also die tabellendaten in cal/g/K * 54.31765

Palladium (Pd):
You have: 4.1868 106.42 joule / mol K N_A
You want: k
    * 53.588401
    *   / 0.018660755
    *

you have: k			(Bolzmann const)
you want: joule / mol K N_A	(J/mol K) 	--> * 8.314472

you have: joule / mol K N_A	(J/mol K)  k			(Bolzmann const)
you want:  k			(Bolzmann const)	--> * 0.12027222

You have: joule / mol K N_A
You want: k
	* 0.12027222
	/ 8.314472

you have: k			(Bolzmann const)
you want: meV/K							--> * 0.086173422

you have: joule / mol K N_A	(J/mol K) 
you want: meV/K			(meV/K)				-->	* 0.010364269

you have: kjoule / mol N_A 	(kJ/mol)
you want: meV			(meV)				--> 	*10.364269

you have: joule / mol N_A	(J/mol K) 
you want: meV							-->	* 0.010364269

you have: k (bolz)
you want: eV/K							--> *8.6173422e-05

You have: eV/K			(eV to Kelvin no, rather not)
You want: k							--> * 11604.506
	
You have: Hartree 2 / 93 k	(Hartree to Kelvin)
You want: K
        * 6790.853						--> * 67901.853
        / 0.00014725691



You have: 2 pi 0.125/3.75 angstrom	SPHINX --> VASP kpoints (sphinx coordinate war 0.110831, VASP war 0.125
You want: 1/bohrradius						gitterkonstante war 3.75 einer 2x2x2 sc,
        * 0.11083062						der faktor 2pi ist wegen der fourier transformation)
        / 9.0227777


1100K to eV of 32 atoms: for 32 (N=32, N*3=96, N*3-3=93, kB in eV/K = 0.0000861733)
echo "1100*0.0000861733*93/2" | bc -l >> 4.4077eV

You have: 3.942103eV*2/93 k
You want: K
        * 983.78834
        / 0.0010164788


  free  energy   TOTEN  =       -59.378999 eV
  free  energy   TOTEN  =       -59.516666 eV
                                                                                                                                                                                                                       
  free  energy   TOTEN  =       -63.526839 eV
                                                                                                                                                                                                                       
  free  energy   TOTEN  =       -59.378999 eV
  free  energy   TOTEN  =       -59.516666 eV
                                                                                                                                                                                                                       
2526 units, 72 prefixes, 56 nonlinear units


U(meV/at) im dUdL hier ist -63.526839 die erste free energy TOTEN und -59.516666 die des entsprechenden schritts
You have: (-59.516666--63.526839)eV/32
You want: meV
        * 125.31790625000005
        / 0.00797970561369796064



ein exemplarischer run:/data/glensk/v/Cu--/fcc4__/PAW_PBE/fah/2x2x2sc_kp02m02m02_230eV/V-run01/3.75Ang_1100K_EDIFF1e-2/lambda0.1_8136
allerste refer, Iteration    1(  14):   free energy    TOTEN  =      -111.70331275 eV
allerste refer, Iteration    1(  14):   energy without entropy =     -110.96301640  energy(sigma->0) =     -111.33316457
                                        energy  without entropy=     -110.963016  energy(sigma->0) =     -111.333165
% ion-electron   TOTEN  =      -111.703313  see above   
kinetic Energy EKIN   =         4.052585  (temperature 1011.35 K)
total energy   ETOTAL =      -107.650728 eV  

erster schritt, Iteration    2(  10): free energy    TOTEN  =      -107.21407167 eV
erster schritt, Iteration    2(  10): energy without entropy =     -106.62850672  energy(sigma->0) =     -106.92128919
                                      energy  without entropy=     -106.628507  energy(sigma->0) =     -106.921289

% ion-electron   TOTEN  =      -107.214072  see above
  kinetic Energy EKIN   =         3.892063  (temperature  971.30 K)
  total energy   ETOTAL =      -103.322009 eV
$units
2526 units, 72 prefixes, 56 nonlinear units

You have: (-107.21407167--111.70331275)eV/32  /31 wenn neu normale zelle (vakanz:...eV/31   /30 wenn NEU) (vak ordner mitd: /31, vak ordern ohne d:/30)
You want: meV
        * 140.28878
        / 0.0071281536
You have: (-106.62850672--110.96301640)eV/32 /31 wenn neu
You want: meV
        * 135.45343
        / 0.0073826113
You have: (-106.92128919--111.33316457)eV/32 /31 wenn neu
You want: meV
        * 137.87111
        / 0.0072531514
$head dUdL 
#  step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset
      1      10.0    971.3    971.3        140.29    137.59          2.70      2.70      0.00


$units
2526 units, 72 prefixes, 56 nonlinear units

You have: 1000 k K 3/2
You want: meV
        * 129.260134447221219
        / 0.00773633730365888301

You have: 3/2* k 1000 K   (3/2 kB T)=Ekin
You want: meV
        * 129.26013


You have: 1100 k K 3/2
You want: meV
        * 142.18615


You have: meV
You want: kjoule / mol N_A
        * 0.09648534
        / 10.364269


http://www.kbfi.ee/thz/?p=energyconverter
1 THz = 33.35640946151377 cm-1
1 THz = 4.135667426829911 meV

$units
2526 units, 72 prefixes, 56 nonlinear units

You have: angstrom^3*GPa
You want: meV
        * 6.2415096
        / 0.16021765
You have: ^C

You have: angstrom^3*0.1GPa
You want: kJ/mol N_A
        * 0.060221418
        / 16.605388

y=s0.y-(2/3)*s1.y-(1/3)*s2.y


You have: angstrom^3*0.1GPa
You want: meV
        * 0.62415096
        / 1.6021765





vacancy concentration in bulk (c=0.001 or c=0.0002; at100Kelvin or at 1360Kelvin)

$units
2526 units, 72 prefixes, 56 nonlinear units

You have: k 1000K * 0.001
You want: meV
        * 0.086173423
        / 11.604506
You have: k 1360K * 0.0002
You want: meV
        * 0.023439171
        / 42.663625




####################################################
convert slope of arrhenius plot to H

You have: -(ln(10^-3/10^-6)/(0.0012-0.0017))k*K
You want: eV
        * 1.1905298
        / 0.83996215

see for this paper of hehenkamp (Cu_vakancies_1994 paper experimental)

        1.6021765e-22 kg m^2 / s^2
You have: THz h
You want: meV
        * 4.1356673
        / 0.24179895


$units
2526 units, 72 prefixes, 56 nonlinear units

You have: hbar(hartree/(bohrradius^2 u))^(1/2)
You want: meV
        * 637.33912
        / 0.0015690234
You have: ^C

