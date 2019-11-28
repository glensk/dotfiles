SUBROUTINE INVERT(MATR,MATR2)
      IMPLICIT NONE
      real(8) :: MATR(3,3)
      real(8) :: MATR2(3,3)
      real(8) :: V1(3),V2(3),V3(3),b1(3),b2(3),b3(3),A1(3),A2(3),A3(3)

      V1 = MATR(1,:)
      V2 = MATR(2,:)
      V3 = MATR(3,:)

      b1(1) = V2(2)*V3(3)-V3(2)*V2(3)
      b1(2) = V2(3)*V3(1)-V3(3)*V2(1)
      b1(3) = V2(1)*V3(2)-V3(1)*V2(2)

      b2(1) = V3(2)*V1(3)-V1(2)*V3(3)
      b2(2) = V3(3)*V1(1)-V1(3)*V3(1)
      b2(3) = V3(1)*V1(2)-V1(1)*V3(2)

      b3(1) = V1(2)*V2(3)-V2(2)*V1(3)
      b3(2) = V1(3)*V2(1)-V2(3)*V1(1)
      b3(3) = V1(1)*V2(2)-V2(1)*V1(2)
      
      A1 = b1/DOT_PRODUCT(b1,V1)
      A2 = b2/DOT_PRODUCT(b2,V2)
      A3 = b3/DOT_PRODUCT(b3,V3)

      MATR2(:,1) = A1
      MATR2(:,2) = A2
      MATR2(:,3) = A3
END SUBROUTINE INVERT
