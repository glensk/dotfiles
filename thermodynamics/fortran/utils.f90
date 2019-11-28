module utils
  contains
  subroutine removeSpaces (str1)
     character(40) str1,str2
     integer*4 i,ls1,ls2
     ls1=len_trim(str1)
     ls2=0
     do i = 1,ls1
       if (str1(i:i).ne.' ') then
         ls2 = ls2 + 1
         str2(ls2:ls2) = str1(i:i)
       endif
     enddo
     do i=ls2+1,40
       str2(i:i)=' '
     enddo
     str1=str2
     return
  end subroutine

end module utils
