Compiler=ifort

band_analysis.x:Band_Analysis.o auxiliary_modules.o
	$(Compiler) -o Band_Analysis.x Band_Analysis.o auxiliary_modules.o
auxiliary_modules.o:auxiliary_modules.f90
	$(Compiler) -c auxiliary_modules.f90
param.o:param.f90
	$(Compiler) -c param.f90
Band_Analysis.o:Band_Analysis.f90 auxiliary_modules.o
	$(Compiler) -c Band_Analysis.f90 auxiliary_modules.f90
clean:
	rm *.o *.mod *.x

