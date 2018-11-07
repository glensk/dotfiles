# runner wants all the positions in bohrradius (as are the lattice vectors)
# runner wants all the energies in hartree
awk '{if (NF ==4){
printf "%s\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\n","atom",$2/0.529117,$3/0.529117,$4/0.529117,$1,0,0,$2,$3,$4
} 
else if (NF>6)
{ 
a=$3/0.529117;b=$4/0.529117;c=$5/0.529117;alpha=$6*0.017453293;beta=$7*0.017453293;gamma=$8*0.017453293; 
h00=a;h01=b*cos(gamma);h02=c*cos(beta);h11=b*sin(gamma);h12=(b*c*cos(alpha)-h02*h01)/h11;h22=sqrt(c*c-h02*h02-h12*h12);
printf("energy\t\t0.000000\n")
printf("charge\t\t0.000000\n")
printf("end\n")
printf("begin\n")
printf("c MD data from the monster under the bed \n")
printf("%s%f\t%f\t%f\n", "lattice\t\t",h00,0,0);
printf("%s%f\t%f\t%f\n", "lattice\t\t",h01,h11,0);
printf("%s%f\t%f\t%f\n", "lattice\t\t",h02,h12,h22);
}
}

' $1 | sed '1,3d'
echo "energy      0.000000"
echo "charge      0.000000"
echo "end"
