program dotProduct

     implicit none
     real*8 a,b,c,d,e,f,n1,n2,dot,ang,pi
     character(len=20) :: str

     pi = acos(0.D0)*2
 
     if ( iargc() .NE. 6) then
       write (0,*) "ERROR: provide two vectors as input"
       stop
     endif

     call getarg(1,str); read (str,*) a
     call getarg(2,str); read (str,*) b
     call getarg(3,str); read (str,*) c
     call getarg(4,str); read (str,*) d
     call getarg(5,str); read (str,*) e
     call getarg(6,str); read (str,*) f

     n1 = (a**2+b**2+c**2)**(0.5)
     n2 = (d**2+e**2+f**2)**(0.5)
     dot = (a*d+b*e+c*f)/n1/n2
     ang = acos(dot)*180/pi

     write (*,*) ang
end


