hallo wie geht es dir [ENTER]
letztes event hier [ENTER]


!$  last argument of last event         # -> hier
!^  first argument of last event        # -> event
!-1 last event                          # -> letztes event hier
!!  last event                          # -> letztes event hier
!:  last event                          # -> letztes event hier
!-2 second last event                   # -> hallo wie geht es dir 
!h  most recent event starting with h   # -> hallo wie geht es dir
!l  most recent event starting with l   # -> letztes event hier
!-1:0  == !:0                           # -> letztes
!-1:1  == !:1                           # -> event
!-1:2  == !:2                           # -> hier
!-2:0                                   # -> hallo
!-1:0-1   (range of words)              # -> letztes event
!#                                      # current event
hallo wie gehts !#                      # -> hallo wie gehts hallo wie gehts
hallo wie gehts !#:0                    # -> hallo wie geths hallo
hallo wie gehts !#:1                    # -> hallo wie geths wie
echo !:0                                # -> exannds to echo hallo
echo !:1

from http://cf.ccmr.cornell.edu/cgi-bin/w3mman2html.cgi?tcsh(1)
           0       The first (command) word
           n       The nth argument
           ^       The first argument, equivalent to `1'
           $       The last argument
           %       The word matched by an ?s? search
           x-y     A range of words
           -y      Equivalent to `0-y'
           *       Equivalent to `^-$', but  returns  nothing  if
                   the event contains only 1 word
           x*      Equivalent to `x-$'
           x-      Equivalent to `x*', but omitting the last word
                   (`$')
