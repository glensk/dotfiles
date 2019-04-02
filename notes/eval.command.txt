abc=10 x=abc
y='$'$x
echo $y               ### -> $abc
eval y='$'$x
echo $y               ### -> 10
