for x in $(ls testpoints.0000*); do awk 'BEGIN{avg=0;count=0;u=27.211385*1000}!/poi/{avg= avg + ($2-$3)^2*u*u; count=count + 1}END{print sqrt(avg/count)}' $x; done > /tmp/tmptest.data 
for x in $(ls trainpoints.0000*); do awk 'BEGIN{avg=0;count=0;u=27.211385*1000}!/poi/{avg= avg + ($2-$3)^2*u*u; count=count + 1}END{print sqrt(avg/count)}' $x; done > /tmp/tmptrain.data

echo "#RMSE_energy (test) [mev/atom] RMSE_energy (train) [mev/atom]"
paste  /tmp/tmptest.data /tmp/tmptrain.data | awk '!/nan/{printf "%15f %15f \n",$1,$2}' 
