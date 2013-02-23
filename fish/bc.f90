

! ******************************************************************************
! 1D OUTFLOW BC'S
! ******************************************************************************
SUBROUTINE OUTFLOW_1D_X0(DATA, NX, NG, NC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NG, NC
  REAL*8, INTENT(INOUT) :: DATA(NC, NX+2*NG)
  INTEGER I
  DO I=1,NG
     DATA(:,I) = DATA(:,NG+1)
  END DO
END SUBROUTINE OUTFLOW_1D_X0

SUBROUTINE OUTFLOW_1D_X1(DATA, NX, NG, NC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NG, NC
  REAL*8, INTENT(INOUT) :: DATA(NC, NX+2*NG)
  INTEGER I
  DO I=1,NG
     DATA(:,NX+NG+I) = DATA(:,NX+NG)
  END DO
END SUBROUTINE OUTFLOW_1D_X1


! ******************************************************************************
! 2D OUTFLOW BC'S
! ******************************************************************************
SUBROUTINE OUTFLOW_2D_X0(DATA, NX, NY, NG, NC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY, NG, NC
  REAL*8, INTENT(INOUT) :: DATA(NC, NY+2*NG, NX+2*NG)
  INTEGER I
  DO I=1,NG
     DATA(:,:,I) = DATA(:,:,NG+1)
  END DO
END SUBROUTINE OUTFLOW_2D_X0

SUBROUTINE OUTFLOW_2D_X1(DATA, NX, NY, NG, NC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY, NG, NC
  REAL*8, INTENT(INOUT) :: DATA(NC, NY+2*NG, NX+2*NG)
  INTEGER I
  DO I=1,NG
     DATA(:,:,NX+NG+I) = DATA(:,:,NX+NG)
  END DO
END SUBROUTINE OUTFLOW_2D_X1

SUBROUTINE OUTFLOW_2D_Y0(DATA, NX, NY, NG, NC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY, NG, NC
  REAL*8, INTENT(INOUT) :: DATA(NC, NY+2*NG, NX+2*NG)
  INTEGER I
  DO I=1,NG
     DATA(:,:,I) = DATA(:,:,NG+1)
  END DO
END SUBROUTINE OUTFLOW_2D_Y0

SUBROUTINE OUTFLOW_2D_Y1(DATA, NX, NY, NG, NC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY, NG, NC
  REAL*8, INTENT(INOUT) :: DATA(NC, NY+2*NG, NX+2*NG)
  INTEGER I
  DO I=1,NG
     DATA(:,NY+NG+I,:) = DATA(:,NY+NG,:)
  END DO
END SUBROUTINE OUTFLOW_2D_Y1
