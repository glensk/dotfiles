program norm

     implicit none
     real*8 a,b,c,n
     character(len=20) :: str

     if ( iargc() .NE. 3) then
       write (0,*) "ERROR: provide vector as input"
       stop
     endif

     call getarg(1,str); read (str,*) a
     call getarg(2,str); read (str,*) b
     call getarg(3,str); read (str,*) c

     n = (a**2+b**2+c**2)**(0.5)

     write (*,*) n
end


