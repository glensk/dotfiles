Compiler=gfortran

Diff_Charge.x:Diff_Charge.o
	$(Compiler) -o Diff_Charge.x Diff_Charge.o
Diff_Charge.o:Diff_Charge.f90
	$(Compiler) -c Diff_Charge.f90
clean:
	rm *.o *.x
