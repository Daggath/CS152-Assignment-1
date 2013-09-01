C     ===== LIN EQ SOLVER MAIN PROGRAM =====
C
C     A FORTRAN IV PROGRAM THAT SOLVES A SYSTEM OF LINEAR EQUATIONS. BUILDS ON
C     THE HILBERT PROGRAM SHARED IN CLASS BY DR. MAK.
C
C     AUTHORS: VENKATA MUNNANGI, SRIKANTH NARAHARI, OLIVER SEET.
C     CLASS: CS 152, TR 6:00 PM TO 7:15 PM.
C
      DOUBLE PRECISION A(4, 4)
      DOUBLE PRECISION B(4, 1)
      DOUBLE PRECISION AINV(4, 4)
      DOUBLE PRECISION X(4, 1)
      DOUBLE PRECISION AXB(4, 1)
      DOUBLE PRECISION AXAINV(4, 4)
      DOUBLE PRECISION AXX(4, 1)
C
      INTEGER N
C
      N = 4
C
C     HARD CODE THE COEFFICIENT MATRIX A
C
      A(1, 1) = 3
      A(1, 2) = 1
      A(1, 3) = -5
      A(1, 4) = 4
      A(2, 1) = 2
      A(2, 2) = -3
      A(2, 3) = 3
      A(2, 4) = -2
      A(3, 1) = 5
      A(3, 2) = -3
      A(3, 3) = 4
      A(3, 4) = 1
      A(4, 1) = -2
      A(4, 2) = 4
      A(4, 3) = -3
      A(4, 4) = -3
C
C     PRINT THE COEFFICIENT MATRIX A
C
      WRITE (6, 100)
  100 FORMAT (/'COEFFICIENT MATRIX A:'/)
      CALL PRINT (A, N, N)
C
C     HARD CODE THE RIGHT-HAND-SIDE VECTOR B
      B(1, 1) = -18
      B(2, 1) = 19
      B(3, 1) = 22
      B(4, 1) = -16
C
C     PRINT THE RIGHT-HAND-SIDE VECTOR B
C
      WRITE (6, 200)
  200 FORMAT (/'RIGHT-HAND-SIDE VECTOR B:'/)
      CALL PRINT (B, N, 1)
C
C     INVERT A
C
      CALL INVERT (N, A, AINV)
C
C     PRINT THE INVERSE OF A
C
      WRITE (6, 300)
  300 FORMAT (/'MATRIX A INVERTED:'/)
      CALL PRINT (AINV, N, N)
C
C     MULTIPLY A BY A INVERSE AND PRINT TO ENSURE WE GET THE IDENTITY MATRIX
C
      CALL MATXMAT (A, AINV, N, AXAINV)
C
      WRITE (6, 400)
  400 FORMAT (/'COMPUTED IDENTITY MATRIX:'/)
      CALL PRINT (AXAINV, N, N)
C
C     COMPUTE AND PRINT INVERSE OF A X B TO GET SOLUTION VECTOR X
C
      CALL MATXVEC (AINV, B, N, X)
C
      WRITE (6, 500)
  500 FORMAT (/'SOLUTION VECTOR X:'/)
      CALL PRINT (X, N, 1)
C
C     COMPUTE AND PRINT A X X TO VERIFY THAT WE GET BACK VECTOR B
C
      CALL MATXVEC (A, X, N, AXX)
C
      WRITE (6, 600)
  600 FORMAT (/'COMPUTED VECTOR B:'/)
      CALL PRINT (B, N, 1)
C
C
      PAUSE
      STOP
      END
C
C     ===== SUBROUTINE DECOMP =====
C
      SUBROUTINE DECOMP (N, A, LU)
      INTEGER N
      DOUBLE PRECISION A(4, 4), LU(4, 4), SCALES(4), PS(4)
C
      INTEGER I, J, K, PIVOTX, NM1, PSI, PSK, PSN, KP1
      DOUBLE PRECISION NRMROW, PIVOT, SIZE, BIGGST, MULT
      COMMON PS
C
C     INITIALIZE PS, LU, AND SCALES.
      DO 200 I = 1, N
          PS(I) = I
          NRMROW = 0.0
C
          DO 100 J = 1, N
              LU(I, J) = A(I, J)
C
C             FIND THE LARGEST ROW ELEMENT.
              IF (NRMROW .GE. DABS(LU(I, J))) GO TO 100
                  NRMROW = DABS(LU(I, J))
  100     CONTINUE
C
C         SET THE SCALING FACTOR FOR ROW EQUILIBRATION.
          IF (NRMROW .NE. 0) GO TO 110
              SCALES(I) = 0.0
              CALL SINGLR(0)
              GO TO 200
C
  110     SCALES(I) = 1.0/NRMROW
  200 CONTINUE
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING.
C     PIVOT ROW K.
      NM1 = N-1
      DO 400 K = 1, NM1
          PIVOTX = 0
          BIGGST = 0.0
C
C         GO DOWN ROWS FROM ROW K.
          DO 300 I = K, N
              PSI = PS(I)
C
C             DIVIDE BY THE BIGGEST ROW ELEMENT.
              SIZE = DABS(LU(PSI, K))*SCALES(PSI)
C
              IF (BIGGST .GE. SIZE) GO TO 300
                  BIGGST = SIZE
                  PIVOTX = I
  300     CONTINUE
C
          IF (BIGGST .NE. 0.0) GO TO 310
              CALL SINGLR(1)
              GO TO 400
C
  310     IF (PIVOTX .EQ. K) GO TO 350
C
C             EXCHANGE ROWS.
              J = PS(K)
              PS(K) = PS(PIVOTX)
              PS(PIVOTX) = J
C
C         PIVOT ELEMENT.
  350     PSK = PS(K)
          PIVOT = LU(PSK, K)
C
C         GO DOWN REST OF ROWS.
          KP1 = K+1
          DO 380 I = KP1, N
              PSI = PS(I)
              MULT = LU(PSI, K)/PIVOT
              LU(PSI, K) = MULT
C
              IF (MULT .EQ. 0.0) GO TO 380
C
C                 INNER LOOP. ONLY THE COLUMN SUBSCRIPT VARIES.
                  DO 360 J = KP1, N
                      LU(PSI, J) = LU(PSI, J) - MULT*LU(PSK, J)
  360             CONTINUE
  380     CONTINUE
  400 CONTINUE
C
C     CHECK THE BOTTOM RIGHT ELEMENT OF THE PERMUTED MATRIX.
      PSN = PS(N)
      IF (LU(PSN, N) .EQ. 0.0) CALL SINGLR(2)
c
      RETURN 
      END
C
C     ===== SUBROUTINE SOLVE =====
C
      SUBROUTINE SOLVE(N, LU, B, X)
      INTEGER N
      DOUBLE PRECISION LU(4, 4), B(4), X(4), PS(4)
C
      INTEGER I, J, IM1, PSI, NP1, IBACK, IP1
      DOUBLE PRECISION DOT
      COMMON PS
C
C     LY = B : SOLVE FOR Y.
      DO 200 I = 1, N
          IM1 = I-1
          PSI = PS(I)
C
          DOT = 0.0
          DO 100 J = 1, IM1
              DOT = DOT + LU(PSI, J)*X(J)
  100     CONTINUE
          X(I) = B(PSI) - DOT
  200 CONTINUE
C
C     UX = Y : SOLVE FOR X
      NP1 = N+1
      DO 400 IBACK = 1, N
          I = NP1-IBACK
C         I = N, N-1, N-2, ..., 2, 1
C
          IP1 = I+1
          PSI = PS(I)
C
          DOT = 0.0
          DO 300 J = IP1, N
              DOT = DOT + LU(PSI, J)*X(J)
  300     CONTINUE
          X(I) = (X(I) - DOT)/LU(PSI, I)
  400 CONTINUE
C
      RETURN
      END
C
C     ===== SUBROUTINE INVERT =====
C
      SUBROUTINE INVERT(N, A, AINV)
      INTEGER N
      DOUBLE PRECISION A(4, 4), AINV(4, 4)
C
C     COMPUTE AINV = INVERSE(A). NOTE THAT A AND AINV
C     ARE OFTEN PASSED THE SAME MATRIX.
C
      DOUBLE PRECISION LU(4, 4), B(4), X(4)
      INTEGER I, J
C
      CALL DECOMP(N, A, LU)
C
      DO 300 J = 1, N
          DO 100 I = 1, N
              B(I) = 0.0
              IF (I .EQ. J) B(I) = 1.0
  100     CONTINUE
c
          CALL SOLVE(N, LU, B, X)
C
          DO 200 I = 1, N
              AINV(I, J) = X(I)
  200     CONTINUE
  300 CONTINUE
C
      RETURN
      END
C
C     ===== SUBROUTINE SINGLR =====
C
      SUBROUTINE SINGLR(WHY)
      INTEGER WHY
C
      GO TO (1, 2, 3), WHY
C
    1 WRITE (6, 10)
   10 FORMAT ('MATRIX WITH ZERO ROW IN DECOMPOSE.')
      RETURN
C
    2 WRITE (6, 20)
   20 FORMAT ('SINGULAR MATRIX IN DECOMPOSE.')
      RETURN
C
    3 WRITE (6, 30)
   30 FORMAT ('NO CONVERGENCE IN IMPROVE.')
      RETURN
C
      END
C
C     ===== SUBROUTINE PRINT =====
C     THIS SUBROUTINE HAS BEEN MODIFIED TO ACCOMODATE THE NUMBER OF COLUMNS
C
      SUBROUTINE PRINT (A, N, M)
      INTEGER N
      INTEGER M
      DOUBLE PRECISION A(4, 4)
C
      INTEGER I, J
C
      DO 100 I = 1, N
          WRITE (6, 10) (A(I, J), J = 1, M)
   10     FORMAT(5F15.6)
  100 CONTINUE
C
      RETURN
      END
C
C     ===== SUBROUTINE MATXMAT =====
C
      SUBROUTINE MATXMAT (A, B, N, AXB)
      INTEGER N
      DOUBLE PRECISION A(4, 4)
      DOUBLE PRECISION B(4, 4)
      DOUBLE PRECISION AXB(4, 4)
C
      INTEGER I, J, K
      DOUBLE PRECISION SUM
C
      DO 100 I = 1, N
          DO 200 J = 1, N
              SUM = 0
              DO 300 K = 1, N
                  SUM = SUM + A(I, K) * B(K, J)
  300         CONTINUE
	      AXB(I, J) = SUM
  200     CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
C     ===== SUBROUTINE MATXVEC =====
C
      SUBROUTINE MATXVEC (A, V, N, AXV)
      INTEGER N
      DOUBLE PRECISION A(4, 4)
      DOUBLE PRECISION V(4, 1)
      DOUBLE PRECISION AXV(4, 1)
C
      INTEGER I, J, K
      DOUBLE PRECISION SUM
C
      DO 100 I = 1, N
          SUM = 0
          DO 300 K = 1, N
              SUM = SUM + A(I, K) * V(K, 1)
  300     CONTINUE
          AXV(I, 1) = SUM
  100 CONTINUE
C
      RETURN
      END
C
