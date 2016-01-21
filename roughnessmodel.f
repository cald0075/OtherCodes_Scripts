      PROGRAM STATICSR

C=====USED FOR CPU ALLOCATION
      IMPLICIT NONE

      INTEGER I,J
      REAL DUM
      INTEGER NK,NTHE
      PARAMETER (NK=512,NTHE=256)

      REAL Z0T
      REAL ACH,G,U10,FETCH,NUA,KAPPA,D2
      REAL SK(NK),SK2D(NK,NTHE),WVN(NK),KP,WTDP,REK(NK)
      REAL BETA(NK),BETA2D(NK,NTHE),USTAR,Z0S,UD2,PI
      REAL KMAX,KMIN,KC,KD,PEX
      REAL USTARDP,ACHDP,CDDP,Z0TDP
      REAL X,Y
      
!--------------------------
!     INITIALIZATION
!--------------------------

      ACH = 0.0144
      G = 9.8
      NUA = 1.48E-5
      KAPPA = 0.41
      D2 = 0.5
      WTDP = 100000.0
      PI = ACOS(-1.0)
      KC = 314.0

      KMAX = 314.0
      KMIN = 0.008
      DO I = 1, NK
         WVN(I) = EXP(LOG(KMIN) + (I - 1) * (LOG(KMAX) - LOG(KMIN))
     &        / (NK - 1))
      ENDDO
      
      U10 = 10
      FETCH = 500000.0
      KD = 50

      GOTO 100
      CALL JONSWAP_DIRC(U10,FETCH,G,NK,NTHE,SK2D,WVN,KP)
      GOTO 500
!      PEX = 0.000996534
!      DO I = 1, NK
!         WVN(I) = I * PEX
!      ENDDO
      GOTO 225
      WTDP = 71.875
      CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
      WTDP = 21.875
      CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
      GOTO 500
      DO I = 1, NK
         READ(39,*) WVN(I),SK(I)
         REK(I) = WVN(I)
      ENDDO
      CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &     REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)
      WRITE(36700+INT(U10),*) WTDP,USTAR**2,ACH,
     &     (KAPPA/LOG(PI/KP/Z0T))**2

 121  CONTINUE
      WTDP = 2.5
      DO I = 1, NK
         READ(39,*) WVN(I),SK(I),SK(I)
         REK(I) = WVN(I)
      ENDDO
      CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &     REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)
      WRITE(36700+INT(U10),*) WTDP,USTAR**2,ACH,
     &     (KAPPA/LOG(PI/KP/Z0T))**2
!      CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)

      GOTO 500

 225  CONTINUE
      
      U10 = 15.0
      WRITE(293,921)
      
      WTDP = 9999.9
!      WRITE(293,976) WTDP
      CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
      PRINT *,"KP = ",KP


      WTDP = 10.0
!      WRITE(293,976) WTDP
      DO I = 1, 30
         CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
         WTDP = WTDP + 5.0
      ENDDO
      GOTO 500

 122  CONTINUE
      WTDP = 10.0
      U10 = 15.0
      WRITE(293,921) 
      WRITE(293,976) WTDP
      DO I = 1, 20
         CALL JONSWAP(U10,FETCH,G,NK,SK,WVN,KP)
!         CALL SOLVER_DYN(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
!     &        D2,BETA,USTAR,Z0S,UD2,KC,KD)
         USTARDP = USTAR
         ACHDP = ACH
         CDDP = (KAPPA/LOG(PI/KP/Z0T))**2
         Z0TDP = Z0T
         WRITE(293,*) U10,USTARDP,ACHDP,CDDP,Z0TDP
         U10 = U10 + 1.0
      ENDDO
 921  FORMAT('VARIABLES=U10,')
      GOTO 500
 123  CONTINUE

      PRINT *,"START HERE!"
      WRITE(8860,*)" VARIABLES=X,Y,WTDP"
      WRITE(8861,*)" VARIABLES=X,Y,TAU"
      WRITE(8862,*)" VARIABLES=X,Y,ACH"
      WRITE(8863,*)" VARIABLES=X,Y,CD"
      WRITE(8864,*)" VARIABLES=X,Y,Z0T"

      DO I = 1, 2500
         READ (8888,*) Y,X,WTDP
         CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
         CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &        REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)
         PRINT *,"TOTAL STRESS IN SHALLOW WATER:",USTAR**2
         WRITE(8860,*) X,Y,WTDP
         WRITE(8861,*) X,Y,USTAR**2
         WRITE(8862,*) X,Y,ACH
         WRITE(8863,*) X,Y,(KAPPA/LOG(PI/KP/Z0T))**2
         WRITE(8864,*) X,Y,Z0T
      ENDDO
      PRINT *,"END HERE!"
      GOTO 500

 124  CONTINUE
      WTDP = 10
      WRITE(36700,*)" VARIABLES=WTDP,U10,ETA,ETAB" 
      WRITE(36700,900) 0.0,40,50
      DO J = 1, 50
         U10 = 4.0
         DO I = 1, 40
         CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
         CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &        REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)
         PRINT *,"TOTAL STRESS IN SHALLOW WATER:",USTAR**2
         WRITE(36700,*) WTDP,U10,USTAR**2,ACH,
     &        (KAPPA/LOG(PI/KP/Z0T))**2,Z0T
         U10 = U10 + 0.3
!         WRITE(886,*) X,Y,WTDP,USTAR**2,ACH,
!     &        (KAPPA/LOG(PI/KP/Z0T))**2,Z0T,USTARDP**2,ACHDP,CDDP,Z0TDP
         ENDDO
         WTDP = WTDP + 0.15
      ENDDO
 900  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
 931  FORMAT('VARIABLES=WTDP,TAU,ACH,CD,Z0T,TAUDP,ACHDP,CDDP,Z0TDP')

      KD = 50
!      CALL JONSWAP(U10,FETCH,G,NK,SK,WVN,KP)
!      CALL SOLVER_DYN(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,
!     &     D2,BETA,USTAR,Z0S,UD2,KC,KD)
C      CALL SOLVER(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,D2,
C     &           BETA,USTAR,Z0S,UD2,KC)
!      CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,KP,WTDP)
!      CALL JONSWAP(U10,FETCH,G,NK,SK,WVN,KP)

!      GOTO 500
!      CALL JONSWAP(U10,FETCH,G,NK,SK,WVN,KP)
!      DO I = 1, NK
!         WRITE(378,*) WVN(I),SK(I)
!      ENDDO
 
C      GOTO 500
C      DO I = 1, NK
C         READ(24,*) WVN(I),SK(I)
C      ENDDO
C      CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
C     &     KP,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)

C      CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,KP,WTDP)
C      CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
C     &     KP,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)


      GOTO 500
 100  CONTINUE
      FETCH = 200000
      U10 = 8.0
      WRITE(38,*)" VARIABLES=X,WTDP,U10,USTAR/C,ALPHACH,CDL,"
     &     //"CD10,USTARSQ"
      WRITE(38,976) U10
      WTDP = 10.0
      DO J = 1,123
!         FETCH = 500
!         WRITE(38,976) U10
         READ(88,*) X,WTDP
         DO I = 1,1

!            CALL DONELAN(U10,FETCH,G,NK,SK,WVN,KP)
!            CALL DONELAN_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
            CALL JONSWAP(U10,FETCH,G,NK,SK,WVN,KP)
!            CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
!            CALL SOLVER(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,D2,
!     &           BETA,USTAR,Z0S,UD2,KC)
            CALL SOLVER_DYN(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &           D2,BETA,USTAR,Z0S,UD2,KC,KD)
!            CALL SOLVER_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,
!     &           D2,WTDP,BETA,USTAR,Z0S,UD2,KC)

!            CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
!     &           REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)

!            CALL SOLVER_DYN(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,
!     &           D2,BETA,USTAR,Z0S,UD2,KC,KD)

C            WRITE(38,*) USTAR/SQRT(G/KP),Z0T*G/USTAR**2,
C     &           (KAPPA/LOG(PI/KP/Z0T))**2
           
!            CALL JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)
!            CALL SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
!     &           REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)

!     sep_9          
!            CALL JONSWAP_DIRC(U10,FETCH,G,NK,NTHE,SK2D,WVN,KP)
!            CALL SOLVER_2D(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,NTHE,WVN,KP,
!     &           SK2D,D2,BETA2D,USTAR,Z0S,UD2,KC)
!     sep_9

!            CALL SOLVER_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,
!     &           D2,WTDP,BETA,USTAR,Z0S,UD2,KC)
            WRITE(38,*) X,WTDP,U10,USTAR/SQRT(G/KP),Z0T*G/USTAR**2,
     &           (KAPPA/LOG(PI/KP/Z0T))**2,(KAPPA/LOG(10.0/Z0T))**2,
     &           USTAR**2
            PRINT *,I,J
!            FETCH = FETCH*EXP(LOG(1000.0)/40)
         ENDDO
!         U10 = U10 + 0.1
      ENDDO
 976  FORMAT(' ZONE T="',F11.5,' "')
 500  CONTINUE

      END PROGRAM STATICSR

C=====PROGRAM HOS_TEST END HERE

      SUBROUTINE INTERPOLATION

      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER NX,NY
      PARAMETER (NX=51,NY=51)
      REAL X(NX),Y(NY),WTDP,DX,DY
      REAL ALLDP(NX,NY),XMAX,YMAX,XMIN,YMIN
      INTEGER COUNT(NX,NY),IX(4),IY(4),REIX,REIY
      REAL W1,W2,W3,W4,REX,REY
      

      YMIN = 35.73
      XMIN = -75.45
      YMAX = 36.0
      XMAX = -75.04

      DX = (XMAX - XMIN) / (NX - 1)
      DY = (YMAX - YMIN) / (NY - 1)
      DO I = 1, NX
         X(I) = XMIN + (I - 1) * DX
      ENDDO
      DO J = 1, NY
         Y(J) = YMIN + (J - 1) * DY
      ENDDO


      END SUBROUTINE INTERPOLATION
!==========================================================================

!==========================================================================
      SUBROUTINE INTP2D(W1,W2,W3,W4,X1,Y1,X2,Y2,X,Y)
!     BY XUANTING HAO
!     THIS SUBROUTINE DETERMINES THE WEIGHTING FACTOR USING THE BILINEAR INTERPOLATION
!     02/16/2014 FIRST EDITION
!     02/22/2014 IMPROVE ROBUSTNESS

!
!                  X2,Y1 W4        X1,Y1  W1
!                       ---------------
!                       |             |
!                       |   X,Y       |
!                       |             |
!                       |             |
!                       ---------------
!               X2,Y2  W3        X1,Y2  W2
!
!   2D INTERPOLATION FORMULA:
!   F(X,Y) = W1 * F(X1,Y1) + W2 * F(X2,Y2) + W3 * F(X3,Y3) + W4 * F(X4,Y4)
!       WHERE:
!       W1 = (X - X2) * (Y - Y2) / (X2 - X1) /(Y2 - Y1)
!       W2 = (X1 - X) * (Y - Y2) / (X2 - X1) /(Y2 - Y1)
!       W3 = (X1 - X) * (Y1 - Y) / (X2 - X1) /(Y2 - Y1)
!       W4 = (X - X2) * (Y1 - Y) / (X2 - X1) /(Y2 - Y1)
      IMPLICIT NONE

      REAL X1,X2,Y1,Y2,X,Y,W1,W2,W3,W4
      IF (X1 .NE. X2 .AND. Y1 .NE. Y2) THEN
         W1 = (X - X2) * (Y - Y2) / (X2 - X1) /(Y2 - Y1)
         W2 = (X1 - X) * (Y - Y2) / (X2 - X1) /(Y2 - Y1)
         W3 = (X1 - X) * (Y1 - Y) / (X2 - X1) /(Y2 - Y1)
         W4 = (X - X2) * (Y1 - Y) / (X2 - X1) /(Y2 - Y1)
      ENDIF

      IF (X1 .EQ. X2 .AND. Y1 .NE. Y2) THEN
         W1 = (Y - Y2) / (Y1 - Y2) / 2
         W2 = (Y1 - Y) / (Y1 - Y2) / 2
         W3 = W2
         W4 = W1
      ENDIF

      IF (X1 .NE. X2 .AND. Y1 .EQ. Y2) THEN
         W1 = (X - X2) / (X1 - X2) / 2
         W2 = W1
         W3 = (X1 - X) / (X1 - X2) / 2
         W4 = W3
      ENDIF

      IF (X1 .EQ. X2 .AND. Y1 .EQ. Y2) THEN
         W1 = 0.25
         W2 = 0.25
         W3 = 0.25
         W4 = 0.25
      ENDIF

      END SUBROUTINE INTP2D

!-----------------------------------------------------------------------

!======================================================================

      SUBROUTINE JONSWAP_DIRC(U10,FETCH,G,NK,NTHE,SK2D,WVN,KP)

!------------------------------------------------------------------
!
!     BY XUANTING HAO
!
!     THIS SUBROUTINE GENERATES A DIRECTIONAL JONSWAP SPECTRUM
!
!------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER NK,NTHE
      REAL SK2D(NK,*),WVN(*),U10,FETCH,G
      REAL ALPHA_P,KP,GAMMA,SIGMA,KMAX,KMIN
      REAL OMEGA,OMEGAP,TWOPI,PI,THETA
      INTEGER I,J

      PI = ACOS(-1.0)
      TWOPI = 2 * ACOS(-1.0)
      ALPHA_P = 0.076 * (U10**2 / G / FETCH)**0.22
      KP = (G * (22**3 / U10 / FETCH)**2)**0.33333
      OMEGAP = 22.*(G**2/U10/FETCH)**(1./3.)
      GAMMA = 3.3
 976  FORMAT(' ZONE T="',F11.2,' "')
!      WRITE(888,902) FETCH,NK,NTHE
      DO J = 1, NTHE
         THETA = -PI / 2.0 + (J - 1.0) / NTHE * PI
         DO I = 1, NK
            OMEGA = SQRT(G * WVN(I))
            IF (WVN(I) .LE. KP) THEN
               SIGMA = 0.07
            ELSE
               SIGMA = 0.09
            ENDIF
            SK2D(I,J) = ALPHA_P * G**2 / OMEGA**5
     1           * EXP(-5./4.*(OMEGA/OMEGAP)**(-4.))
     1           * GAMMA**(EXP(-(OMEGA-OMEGAP)**2/2/SIGMA**2/OMEGAP**2))
     &           * (4.0 / TWOPI) * COS(THETA)**2

            SK2D(I,J) = SK2D(I,J) * G / 2. / OMEGA
!            WRITE(888,901) OMEGA,THETA,SK2D(I,J)
!            WRITE(3600+INT(U10),901) WVN(I),SK(I)
         ENDDO
      ENDDO
 901  FORMAT(25E12.4)
 902  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')

!      PRINT *,"KP = ",KP
      END SUBROUTINE JONSWAP_DIRC
!=========================================================================

!======================================================================

      SUBROUTINE JONSWAP_DIRC_FD(U10,FETCH,G,NK,NTHE,SK2D,WVN,KP,WTDP)

!------------------------------------------------------------------
!
!     BY XUANTING HAO
!
!     THIS SUBROUTINE GENERATES A DIRECTIONAL JONSWAP SPECTRUM
!     IN COASTAL WATER.
!
!------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER NK,NTHE
      REAL SK2D(NK,*),WVN(*),U10,FETCH,G,CG
      REAL ALPHA_P,KP,GAMMA,SIGMA,KMAX,KMIN
      REAL OMEGA,OMEGAP,TWOPI,PI,THETA,WTDP
      INTEGER I,J

      PI = ACOS(-1.0)
      TWOPI = 2 * ACOS(-1.0)
      ALPHA_P = 0.076 * (U10**2 / G / FETCH)**0.22
      OMEGAP = 22.*(G**2/U10/FETCH)**(1./3.)
      CALL DP2SHA(OMEGAP,KP,G,WTDP)

      GAMMA = 3.3
 976    FORMAT(' ZONE T="',F11.2,' "')
!      WRITE(888,902) FETCH,NK,NTHE
      DO J = 1, NTHE
         THETA = -PI / 2.0 + (J - 1.0) / NTHE * PI
         DO I = 1, NK
            OMEGA = SQRT(G * WVN(I))
            IF (WVN(I) .LE. KP) THEN
               SIGMA = 0.07
            ELSE
               SIGMA = 0.09
            ENDIF
            SK2D(I,J) = ALPHA_P * G**2 / OMEGA**5
     1           * EXP(-5./4.*(OMEGA/OMEGAP)**(-4.))
     1           * GAMMA**(EXP(-(OMEGA-OMEGAP)**2/2/SIGMA**2/OMEGAP**2))
     &           * (4.0 / TWOPI) * COS(THETA)**2

            SK2D(I,J) = SK2D(I,J) * G / 2. / OMEGA
!            WRITE(888,901) OMEGA,THETA,SK2D(I,J)
!            WRITE(3600+INT(U10),901) WVN(I),SK(I)
         ENDDO
      ENDDO
 901    FORMAT(25E12.4)
 902      FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')

!      PRINT *,"KP = ",KP
      END SUBROUTINE JONSWAP_DIRC_FD
!=========================================================================


!=====================================================================
      SUBROUTINE DONELAN(U10,FETCH,G,NK,SK,WVN,KP)

      IMPLICIT NONE
      INTEGER NK
      REAL SK(*),WVN(*),U10,FETCH,G
      REAL ALPHA_P,KP,GAMMA,SIGMA,KMAX,KMIN
      REAL OMEGA,OMEGAP,CP
      INTEGER I

      KP = (G * (22**3 / U10 / FETCH)**2)**0.33333
      OMEGAP = 22.*(G**2/U10/FETCH)**(1./3.)
      CP = OMEGAP / KP
      ALPHA_P = 0.006 * (U10 / CP)**0.55
      IF (U10 / CP .LT. 1) THEN
         GAMMA = 1.7
      ELSE
         GAMMA = 1.7 + 6.0 * LOG(U10 / CP) / LOG(10.0)
      ENDIF
      SIGMA = 0.08 * (1 + 4.0 / (U10 / CP)**3)
!      PRINT *,"KP=",KP
!      PRINT *,"ALPHA_P=",ALPHA_P
 976    FORMAT(' ZONE T="',F11.2,' "')
!      WRITE(3500+INT(U10),976) 99999.9
      WRITE(3600+INT(U10),976) 99999.9

      DO I = 1, NK
         OMEGA = SQRT(G * WVN(I))
         SK(I) = ALPHA_P * G**2 / OMEGA**4 / OMEGAP
     1           * EXP(-1.0*(OMEGA/OMEGAP)**(-4.))
     1           * GAMMA**(EXP(-(OMEGA-OMEGAP)**2/2/SIGMA**2/OMEGAP**2))
         IF (ISNAN(SK(I))) THEN
            PRINT *,ALPHA_P,OMEGA,WVN(I),GAMMA
            STOP
         ENDIF
!         WRITE(3500+INT(U10),*) OMEGA,SK(I)
         SK(I) = SK(I) * G / 2. / OMEGA
         WRITE(3600+INT(U10),901) WVN(I),SK(I)
      ENDDO
 901    FORMAT(25E12.4)

!      PRINT *,"KP = ",KP
      END SUBROUTINE DONELAN
!=========================================================================

!=====================================================================
      SUBROUTINE DONELAN_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)

      IMPLICIT NONE
      INTEGER NK
      REAL SK(*),WVN(*),U10,FETCH,G
      REAL ALPHA_P,KP,GAMMA,SIGMA,KMAX,KMIN,WTDP
      REAL OMEGA,OMEGAP,CP,KDP,CG,REK(*),MYTANH
      INTEGER I

      KP = (G * (22**3 / U10 / FETCH)**2)**0.33333
      OMEGAP = 22.*(G**2/U10/FETCH)**(1./3.)
      CP = OMEGAP / KP
      ALPHA_P = 0.006 * (U10 / CP)**0.55
      IF (U10 / CP .LT. 1) THEN
         GAMMA = 1.7
      ELSE
         GAMMA = 1.7 + 6.0 * LOG(U10 / CP) / LOG(10.0)
      ENDIF
      REK(1:NK)=WVN(1:NK)
      SIGMA = 0.08 * (1 + 4.0 / (U10 / CP)**3)
!      PRINT *,"KP=",KP
!      PRINT *,"ALPHA_P=",ALPHA_P
 976      FORMAT(' ZONE T="',F11.2,' "')
!      WRITE(3500+INT(U10),976) 99999.9
      WRITE(3600+INT(U10),976) 99999.9

      DO I = 1, NK
         OMEGA = ( G * WVN(I) * MYTANH(REK(I),WTDP))**0.5
         KDP = OMEGA**2 / G
         CG = (G / 2 / OMEGA) * (MYTANH(REK(I),WTDP) + REK(I) * WTDP
     &        / COSH(REK(I)*WTDP)**2)
         
         SK(I) = ALPHA_P * G**2 / OMEGA**4 / OMEGAP
     1           * EXP(-1.0*(OMEGA/OMEGAP)**(-4.))
     1           * GAMMA**(EXP(-(OMEGA-OMEGAP)**2/2/SIGMA**2/OMEGAP**2))
         IF (ISNAN(SK(I))) THEN
            PRINT *,ALPHA_P,OMEGA,WVN(I),GAMMA
            STOP
         ENDIF
!         WRITE(3500+INT(U10),*) OMEGA,SK(I)
         SK(I) = SK(I) * G / 2. / OMEGA
         WRITE(3600+INT(U10),901) WVN(I),SK(I)
      ENDDO
 901      FORMAT(25E12.4)

!      PRINT *,"KP = ",KP
      END SUBROUTINE DONELAN_FD
!=========================================================================



!=====================================================================
      SUBROUTINE JONSWAP(U10,FETCH,G,NK,SK,WVN,KP)

      IMPLICIT NONE
      INTEGER NK
      REAL SK(*),WVN(*),U10,FETCH,G
      REAL ALPHA_P,KP,GAMMA,SIGMA,KMAX,KMIN
      REAL OMEGA,OMEGAP
      INTEGER I

      ALPHA_P = 0.076 * (U10**2 / G / FETCH)**0.22
      KP = (G * (22**3 / U10 / FETCH)**2)**0.33333
      OMEGAP = 22.*(G**2/U10/FETCH)**(1./3.)
!      PRINT *,"KP=",KP
!      PRINT *,"ALPHA_P=",ALPHA_P
      GAMMA = 3.3
 976  FORMAT(' ZONE T="',F11.2,' "')
!      WRITE(3500+INT(U10),976) 99999.9
      WRITE(3600+INT(U10),976) 99999.9

      DO I = 1, NK
         OMEGA = SQRT(G * WVN(I))
         IF (WVN(I) .LE. KP) THEN
            SIGMA = 0.07
         ELSE
            SIGMA = 0.09
         ENDIF
         SK(I) = ALPHA_P * G**2 / OMEGA**5
     1           * EXP(-5./4.*(OMEGA/OMEGAP)**(-4.))
     1           * GAMMA**(EXP(-(OMEGA-OMEGAP)**2/2/SIGMA**2/OMEGAP**2))

!         WRITE(3500+INT(U10),*) OMEGA,SK(I)
         SK(I) = SK(I) * (G / 2. / OMEGA)
         WRITE(3600+INT(U10),901) WVN(I),SK(I)         
      ENDDO
 901  FORMAT(25E12.4)

!      PRINT *,"KP = ",KP
      END SUBROUTINE JONSWAP
!=========================================================================

      SUBROUTINE DP2SHA(OMEGA,K,G,WTDP)

      IMPLICIT NONE
      REAL OMEGA,K,G,WTDP,MYTANH
      
      INTEGER I,ITMAX
      REAL TMP

      TMP = K
      ITMAX = 500
      IF (K * WTDP .GT. 10) RETURN
      DO I = 1, ITMAX
         K = OMEGA**2 / G / MYTANH(K,WTDP)

         IF (ABS(TMP - K) .LT. 1.0E-8) EXIT
         TMP = K
      ENDDO

      END SUBROUTINE DP2SHA

!========================================================================
      FUNCTION MYTANH(K,WTDP)

      IMPLICIT NONE
      REAL K,WTDP
      REAL MYTANH
      IF (K * WTDP .GT. 10) THEN
         MYTANH = 1.0
      ELSE
         MYTANH = TANH(K*WTDP)
      ENDIF

      END FUNCTION MYTANH

!=====================================================================
      SUBROUTINE JONSWAP_FD(U10,FETCH,G,NK,SK,WVN,REK,KP,WTDP)

      IMPLICIT NONE
      INTEGER NK
      REAL SK(*),WVN(*),U10,FETCH,G
      REAL ALPHA_P,KP,GAMMA,SIGMA,KMAX,KMIN,WTDP,OMEGA,OMEGAP,CG
      REAL KDP
      INTEGER I
      REAL MYTANH,TMP,REK(NK)      

      ALPHA_P = 0.076 * (U10**2 / G / FETCH)**0.22
      OMEGAP = 22.*(G**2/U10/FETCH)**(1./3.)
      CALL DP2SHA(OMEGAP,KP,G,WTDP)
!      PRINT *,"KP=",KP
!      PRINT *,"ALPHA_P=",ALPHA_P
      GAMMA = 3.3

      REK(1:NK) = WVN(1:NK)
      TMP = 0
 976  FORMAT(' ZONE T="',F11.0,' "')
!      WRITE(3500+INT(U10),976) WTDP
      WRITE(3600+INT(U10),976) WTDP
      DO I = 1, NK
         OMEGA = ( G * WVN(I) * MYTANH(REK(I),WTDP))**0.5
         KDP = OMEGA**2 / G
         CG = (G / 2 / OMEGA) * (MYTANH(REK(I),WTDP) + REK(I) * WTDP  
     &          / COSH(REK(I)*WTDP)**2)

         IF ( OMEGA .LT. OMEGAP ) THEN
            SIGMA = 0.07
         ELSE
            SIGMA = 0.09
         END IF
         
         SK(I) = ALPHA_P * G**2 / OMEGA**5
     1           * EXP(-5./4.*(OMEGA/OMEGAP)**(-4.))
     1           * GAMMA**(EXP(-(OMEGA-OMEGAP)**2/2/SIGMA**2/OMEGAP**2))

!         WRITE(3500+INT(U10),*) OMEGA,SK(I)
         SK(I) = SK(I) * (G / 2. / OMEGA)! * (KDP / WVN(I))**3.0
         
!         IF (SK(I) .GT. TMP) THEN
!            TMP = SK(I) 
!            KP = WVN(I)
!         ENDIF
         WRITE(3600+INT(U10),*) REK(I),SK(I)
      ENDDO

!      PRINT *,"KP = ",KP
      END SUBROUTINE JONSWAP_FD
!=========================================================================


!=========================================================================
      SUBROUTINE SOLVER(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,D2,
     &     BETA,USTAR,Z0S,UD2,KC)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,U10,FETCH,NUA,KAPPA,SK(*),WVN(*),KP,D2,KC
      REAL,INTENT(INOUT)::Z0T,ACH,BETA(*),USTAR,Z0S,UD2

      REAL ZC(NK),ZETAC,HEAVI
      REAL T1,T2,NEW_Z0T
      REAL EPS,PI

      INTEGER I,J,IT,ITMAX

      PI = ACOS(-1.0)
      Z0T = 0.001
      ITMAX = 500
      PRINT *,"NK=",NK
      WRITE(93,976) FETCH
      DO I = 1, NK
          WRITE(93,*) WVN(I),SK(I),SQRT(G/WVN(I))
      ENDDO
!--------------------------
!     STEP (a)
!--------------------------
      DO I = 1, 100
         USTAR = KAPPA * U10 / LOG(10.0/Z0T)
         Z0T = ACH * USTAR**2 / G
      ENDDO

      EPS = 1.0
 976      FORMAT(' ZONE T="',F11.5,' "')
 901            FORMAT(25E12.4)
      DO IT = 1, ITMAX
!--------------------------
!     STEP (b)
!--------------------------
         UD2 = (USTAR / KAPPA) * LOG(D2 / Z0T)
!--------------------------
!     STEP (c)
!--------------------------
         Z0S = 0.11 * NUA / USTAR
!--------------------------
!     STEP (d),(e)
!--------------------------
         DO I = 1, NK
            ZC(I) = Z0T * EXP((KAPPA / USTAR)
     &           * SQRT(G / WVN(I)))
            ZETAC = ZC(I) * WVN(I)
            BETA(I) = PI * ZETAC * (LOG(0.281  / ZETAC))**4
     &           * HEAVI(0.281  - ZETAC) / KAPPA**2+
     &           2*LOG(0.281  / ZETAC)
            IF (I .EQ. 506) THEN
               PRINT *,"TEST HERE!"
            ENDIF
         ENDDO
!--------------------------
!     STEP (f)
!--------------------------
         T1 = 0
         DO I = 1, NK-1
            IF (WVN(I+1) .LE. KC) THEN
               T1 = T1 + 0.5 * (WVN(I+1)-WVN(I)) * (WVN(I)**2*BETA(I)
     &              *SK(I)+WVN(I+1)**2 * BETA(I+1) * SK(I+1))
            ENDIF
         ENDDO
         T1 = T1 * USTAR**2
         T2 = (KAPPA * UD2)**2 / (LOG(D2 / Z0S))**2
         NEW_Z0T = D2 * EXP(-KAPPA * UD2 / (T1 + T2)**0.5)
         EPS = ABS(NEW_Z0T - Z0T)
         Z0T = Z0T*0.9 + NEW_Z0T*0.1
!--------------------------
!     STEP (g)
!--------------------------
         USTAR = USTAR*0.5 + KAPPA * U10 / LOG(10.0/Z0T)*0.5
         IF (EPS .LT. 1.0E-8) EXIT
      END DO
      ACH = Z0T*G/USTAR**2
      RETURN
      END SUBROUTINE SOLVER
!=========================================================================

!=========================================================================
      SUBROUTINE SOLVER_2D(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,NTHE,WVN,KP,
     &     SK2D,D2,BETA2D,USTAR,Z0S,UD2,KC)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK,NTHE
      REAL,INTENT(IN)::G,U10,FETCH,NUA,KAPPA,SK2D(NK,*),WVN(*),KP,D2,KC
      REAL,INTENT(INOUT)::Z0T,ACH,BETA2D(NK,*),USTAR,Z0S,UD2

      REAL ZC(NK,NTHE),ZETAC,HEAVI,THETA
      REAL T1,T2,TS1,NEW_Z0T
      REAL EPS,PI,TWOPI

      INTEGER I,J,IT,ITMAX

      ACH = 0.0144
      PI = ACOS(-1.0)
      TWOPI = 2 * ACOS(-1.0)
      Z0T = 0.001
      ITMAX = 800
!--------------------------
!     STEP (a)
!--------------------------
      DO I = 1, 100
         USTAR = KAPPA * U10 / LOG(10.0/Z0T)
         Z0T = ACH * USTAR**2 / G
      ENDDO

      DO I = 1, NK
         DO J = 1, NTHE
            BETA2D(I,J) = 0.0            
         ENDDO
      ENDDO

      EPS = 1.0
 902  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
      WRITE(889,902) U10,NK,NTHE 
 976  FORMAT(' ZONE T="',F11.5,' "')
 901  FORMAT(25E12.4)
      DO IT = 1, ITMAX
!--------------------------
!     STEP (b)
!--------------------------
         UD2 = (USTAR / KAPPA) * LOG(D2 / Z0T)
!--------------------------
!     STEP (c)
!--------------------------
         Z0S = 0.11 * NUA / USTAR
!--------------------------
!     STEP (d),(e)
!--------------------------
         DO J = 1, NTHE
            THETA = -PI / 2.0 + (J - 1.0) / NTHE * PI
            IF (ABS(COS(THETA)) .LT. 0.002) THEN
               BETA2D(:,J) = 0.0
            ELSE
               DO I = 1, NK
                  IF ((KAPPA / USTAR / COS(THETA))
     &                 * SQRT(G / WVN(I)) .GT. 5) THEN
                     BETA2D(I,J) = 0.0
                  ELSE
                     ZC(I,J) = Z0T * EXP((KAPPA / USTAR / COS(THETA))
     &                    * SQRT(G / WVN(I)))
                     ZETAC = ZC(I,J) * WVN(I)
                     IF (0.281 .LT. ZETAC) THEN
                        BETA2D(I,J) = 2 * LOG(0.281  / ZETAC)
                     ELSE
                        BETA2D(I,J) = PI * ZETAC 
     &                       * (LOG(0.281  / ZETAC))**4
     &                       * HEAVI(0.281  - ZETAC) / KAPPA**2 +
     &                       2 * LOG(0.281  / ZETAC)
                     ENDIF
                     IF (ISNAN(BETA2D(I,J))) THEN
                        PRINT *,ZC(I,J),WVN(I),THETA,COS(THETA),
     &                       BETA2D(I,J),(KAPPA / USTAR / COS(THETA))
     &                       * SQRT(G / WVN(I))
                        STOP
!                        PRINT *,"BUGS HERE!"
                     ENDIF
                     IF (BETA2D(I,J) .LT. 0.0) THEN
                        BETA2D(I,J) = 0.0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

!         IF (IT .EQ. 1) THEN
!             DO J = 1, NTHE
!                THETA = -PI / 2.0 + (J - 1.0) / NTHE * PI
!                DO I = 1, NK
!                   WRITE(889,901) WVN(I),THETA,BETA2D(I,J)
!                ENDDO
!             ENDDO
!         ENDIF

!--------------------------
!     STEP (f)
!--------------------------
         T1 = 0
         DO J = 1, NTHE
            THETA = -PI / 2.0 + (J - 1.0) / NTHE * PI
            TS1 = 0
            DO I = 1, NK-1
               IF (WVN(I+1) .LE. KC) THEN
                  TS1 = TS1 + 0.5 * (WVN(I+1)-WVN(I)) 
     &                 * (WVN(I)**3 * BETA2D(I,J) * SK2D(I,J)  
     &                 + WVN(I+1)**3 * BETA2D(I+1,J) * SK2D(I+1,J))
               ENDIF
            ENDDO
            T1 = T1 + TS1 * (COS(THETA))**3  * (PI / NTHE)
         ENDDO
         T1 = T1 * USTAR**2
         T2 = (KAPPA * UD2)**2 / (LOG(D2 / Z0S))**2
         NEW_Z0T = D2 * EXP(-KAPPA * UD2 / (T1 + T2)**0.5)
         EPS = ABS(NEW_Z0T - Z0T)
         Z0T = Z0T*0.8 + NEW_Z0T*0.2
!--------------------------
!     STEP (g)
!--------------------------
         USTAR = USTAR*0.5 + KAPPA * U10 / LOG(10.0/Z0T)*0.5
         IF (ISNAN(USTAR)) THEN
            PRINT *,"BUGS HERE!"
            STOP
         ENDIF
         IF (EPS .LT. 1.0E-3*Z0T) EXIT
      END DO
      ACH = Z0T*G/USTAR**2
      RETURN
      END SUBROUTINE SOLVER_2D
!=========================================================================


!=========================================================================
      SUBROUTINE SOLVER_DYN_2D(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &     D2,BETA,USTAR,Z0S,UD2,KC,KD)

!-----------------------------------------------------
!
!     THIS SUBROUTINE CALCULATES THE TOTAL SURFACE ROUGHNESS WITH
!     THE DYNAMIC METHOD
!
!     07/18/2015
!
!-----------------------------------------------------

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,U10,FETCH,NUA,KAPPA,SK(*),WVN(*),D2,KC,KD
      REAL,INTENT(INOUT)::Z0T,ACH,BETA(*),USTAR,Z0S,UD2

      REAL ZC(NK),ZETAC,HEAVI,AW,REK(NK),WTDP
      REAL T1,T2,NEW_Z0T
      REAL EPS,PI

      INTEGER I,J,IT,ITMAX

      PI = ACOS(-1.0)
      Z0T = 0.001
      ITMAX = 50
      REK(1:NK) = WVN(1:NK)
      WTDP = 100000.0
!--------------------------
!     STEP (a)
!--------------------------
      DO I = 1, 100
         USTAR = KAPPA * U10 / LOG(10.0/Z0T)
         Z0T = ACH * USTAR**2 / G
      ENDDO

      EPS = 1.0
C      Z0T = 3.05E-4
C      USTAR = 0.1536
      DO IT = 1, ITMAX
!--------------------------
!     STEP (b)
!--------------------------
         UD2 = (USTAR / KAPPA) * LOG(D2 / Z0T)

!--------------------------
!     STEP (c)
!--------------------------
         Z0S = 0.11 * NUA / USTAR

!--------------------------
!     STEP (d),(e)
!--------------------------
         DO I = 1, NK
            ZC(I) = Z0T * EXP((KAPPA / USTAR)
     &           * SQRT(G / WVN(I)))
            ZETAC = ZC(I) * WVN(I)
            BETA(I) = PI * ZETAC * (LOG(0.281  / ZETAC))**4
     &           * HEAVI(0.281  - ZETAC) / KAPPA**2+
     &           2*LOG(0.281  / ZETAC)
         ENDDO

!--------------------------
!     STEP (f)
!--------------------------
         CALL FINDZ0T(NEW_Z0T,AW,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &        BETA,SK,WVN,REK,WTDP)
         EPS = ABS(NEW_Z0T - Z0T)
         Z0T = Z0T*0.9 + NEW_Z0T*0.1

!--------------------------
!     STEP (g)
!--------------------------
         USTAR = USTAR*0.5 + KAPPA * U10 / LOG(10.0/Z0T)*0.5
         IF (EPS .LT. 1.0E-8) EXIT
      END DO
      ACH = Z0T*G/USTAR**2

      RETURN
      END SUBROUTINE SOLVER_DYN_2D
!=========================================================================

!=========================================================================
      SUBROUTINE SOLVER_DYN(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &     D2,BETA,USTAR,Z0S,UD2,KC,KD)

!-----------------------------------------------------
!
!     THIS SUBROUTINE CALCULATES THE TOTAL SURFACE ROUGHNESS WITH
!     THE DYNAMIC METHOD
!
!     07/18/2015
!
!-----------------------------------------------------

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,U10,FETCH,NUA,KAPPA,SK(*),WVN(*),D2,KC,KD
      REAL,INTENT(INOUT)::Z0T,ACH,BETA(*),USTAR,Z0S,UD2

      REAL ZC(NK),ZETAC,HEAVI,AW,REK(NK),WTDP
      REAL T1,T2,NEW_Z0T
      REAL EPS,PI

      INTEGER I,J,IT,ITMAX

      PI = ACOS(-1.0)
      Z0T = 0.001
      ITMAX = 50
      REK(1:NK) = WVN(1:NK) 
      WTDP = 100000.0
!--------------------------
!     STEP (a)
!--------------------------
      DO I = 1, 100
         USTAR = KAPPA * U10 / LOG(10.0/Z0T)
         Z0T = ACH * USTAR**2 / G
      ENDDO

      EPS = 1.0
C      Z0T = 3.05E-4
C      USTAR = 0.1536
      DO IT = 1, ITMAX
!--------------------------
!     STEP (b)
!--------------------------
         UD2 = (USTAR / KAPPA) * LOG(D2 / Z0T)
         
!--------------------------
!     STEP (c)
!--------------------------
         Z0S = 0.11 * NUA / USTAR

!--------------------------
!     STEP (d),(e)
!--------------------------
         DO I = 1, NK
            ZC(I) = Z0T * EXP((KAPPA / USTAR)
     &           * SQRT(G / WVN(I)))
            ZETAC = ZC(I) * WVN(I)
            BETA(I) = PI * ZETAC * (LOG(0.281  / ZETAC))**4
     &           * HEAVI(0.281  - ZETAC) / KAPPA**2+ 
     &           2*LOG(0.281  / ZETAC)
         ENDDO

!--------------------------
!     STEP (f)
!--------------------------
         CALL FINDZ0T(NEW_Z0T,AW,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &        BETA,SK,WVN,REK,WTDP)
         EPS = ABS(NEW_Z0T - Z0T)
         Z0T = Z0T*0.9 + NEW_Z0T*0.1
         
!--------------------------
!     STEP (g)
!--------------------------
         USTAR = USTAR*0.5 + KAPPA * U10 / LOG(10.0/Z0T)*0.5
         IF (EPS .LT. 1.0E-8) EXIT
      END DO
      ACH = Z0T*G/USTAR**2

      RETURN
      END SUBROUTINE SOLVER_DYN
!=========================================================================

!=========================================================================

!=========================================================================
      SUBROUTINE SOLVER_DYN_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,
     &     REK,D2,BETA,USTAR,Z0S,UD2,KC,KD,WTDP)

!-----------------------------------------------------
!
!     THIS SUBROUTINE CALCULATES THE TOTAL SURFACE ROUGHNESS WITH
!     THE DYNAMIC METHOD
!
!     07/18/2015
!
!-----------------------------------------------------

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,U10,FETCH,NUA,KAPPA,SK(*),WVN(*),REK(*)
      REAL,INTENT(IN)::D2,KC,KD
      REAL,INTENT(INOUT)::Z0T,ACH,BETA(*),USTAR,Z0S,UD2,WTDP

      REAL ZC(NK),ZETAC,HEAVI,AW
      REAL T1,T2,NEW_Z0T
      REAL EPS,PI,C(NK)

      INTEGER I,J,IT,ITMAX
      REAL MYTANH

      PI = ACOS(-1.0)
      Z0T = 0.001
      ITMAX = 50
!      PRINT *,"NK=",NK
!      WRITE(93,976) FETCH
      DO I = 1, NK
         C(I) = SQRT(G * MYTANH(REK(I),WTDP) / REK(I))
C         WRITE(93,*) WVN(I),SK(I),C(I)
      ENDDO

!--------------------------
!     STEP (a)
!--------------------------
      DO I = 1, 100
         USTAR = KAPPA * U10 / LOG(10.0/Z0T)
         Z0T = ACH * USTAR**2 / G
      ENDDO
!      PRINT *,"Z0T=",Z0T,"USTAR=",USTAR,"ALPHA=",Z0T*G/USTAR**2
C      ENDDO

      EPS = 1.0
 976  FORMAT(' ZONE T="',F11.5,' "')
 901  FORMAT(25E12.4)
C      Z0T = 3.05E-4
C      USTAR = 0.1536
      DO IT = 1, ITMAX
C      DO WHILE (.TRUE.)  
!--------------------------
!     STEP (b)
!--------------------------
         UD2 = (USTAR / KAPPA) * LOG(D2 / Z0T)
         
!--------------------------
!     STEP (c)
!--------------------------
         Z0S = 0.11 * NUA / USTAR

!--------------------------
!     STEP (d),(e)
!--------------------------
C         WRITE(35,976) IT * 1.0
         DO I = 1, NK
            ZC(I) = Z0T * EXP((KAPPA / USTAR) * C(I))
            ZETAC = ZC(I) * REK(I)
            BETA(I) = PI * ZETAC * (LOG(0.281  / ZETAC))**4
     &           * HEAVI(0.281  - ZETAC) / KAPPA**2+ 
     &           2*LOG(0.281  / ZETAC)
C            WRITE(35,901) REK(I),BETA(I)
         ENDDO

!--------------------------
!     STEP (f)
!--------------------------

         CALL FINDZ0T(NEW_Z0T,AW,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &        BETA,SK,WVN,REK,WTDP)
         EPS = ABS(NEW_Z0T - Z0T)
         Z0T = Z0T*0.5 + NEW_Z0T*0.5
         
!--------------------------
!     STEP (g)
!--------------------------
         USTAR = USTAR*0.5 + KAPPA * U10 / LOG(10.0/Z0T)*0.5
C         PRINT *,"Z0T=",Z0T,"USTAR=",USTAR,"ALPHA=",Z0T*G/USTAR**2

         IF (EPS .LT. 1.0E-8) EXIT
      END DO
!      PRINT *,"Z0T=",Z0T,"USTAR=",USTAR,"ALPHA=",Z0T*G/USTAR**2
!      PRINT *,"EPS=",EPS
      ACH = Z0T*G/USTAR**2
!      WRITE(38,*) USTAR/SQRT(G/KP),Z0T*G/USTAR**2,
!     &     (KAPPA/LOG(PI/KP/Z0T))**2

      RETURN
      END SUBROUTINE SOLVER_DYN_FD
!=========================================================================


!=========================================================================

!=========================================================================
      SUBROUTINE SOLVER_FD(Z0T,ACH,G,U10,FETCH,NUA,KAPPA,NK,SK,WVN,KP,
     &     D2,WTDP,BETA,USTAR,Z0S,UD2,KC)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,U10,FETCH,NUA,KAPPA,SK(*),WVN(*),KP,D2,KC,WTDP
      REAL,INTENT(INOUT)::Z0T,ACH,BETA(*),USTAR,Z0S,UD2

      REAL ZC(NK),ZETAC,HEAVI,T1,T2,NEW_Z0T
      REAL EPS,PI,C(NK)

      INTEGER I,J,IT,ITMAX
      REAL MYTANH

      PI = ACOS(-1.0)
      Z0T = 0.001
      ITMAX = 500
!      PRINT *,"NK=",NK
      WRITE(93,976) FETCH
      DO I = 1, NK         
         C(I) = SQRT(G * MYTANH(WVN(I),WTDP) / WVN(I))
!         WRITE(93,*) WVN(I),SK(I),C(I)
      ENDDO
!--------------------------
!     STEP (a)
!--------------------------
      DO I = 1, 100
         USTAR = KAPPA * U10 / LOG(10.0/Z0T)
         Z0T = ACH * USTAR**2 / G
      ENDDO
C      PRINT *,"Z0T=",Z0T,"USTAR=",USTAR,"ALPHA=",Z0T*G/USTAR**2
C      ENDDO

      EPS = 1.0
 976  FORMAT(' ZONE T="',F11.5,' "')
 901  FORMAT(25E12.4)
C     Z0T = 3.05E-4
C      USTAR = 0.1536
      DO IT = 1, ITMAX
C      DO WHILE (.TRUE.)  
!--------------------------
!     STEP (b)
!--------------------------
         UD2 = (USTAR / KAPPA) * LOG(D2 / Z0T)

!--------------------------
!     STEP (c)
!--------------------------
         Z0S = 0.11 * NUA / USTAR

!--------------------------
!     STEP (d),(e)
!--------------------------
C         WRITE(35,976) IT * 1.0
         DO I = 1, NK
            ZC(I) = Z0T * EXP((KAPPA / USTAR) * C(I))
            ZETAC = ZC(I) * WVN(I)
            BETA(I) = PI * ZETAC * (LOG(0.281  / ZETAC))**4
     &           * HEAVI(0.281  - ZETAC) / KAPPA**2+ 
     &           2*LOG(0.281  / ZETAC)
C            WRITE(35,901) WVN(I),BETA(I)
         ENDDO

!--------------------------
!     STEP (f)
!--------------------------
         T1 = 0
         DO I = 1, NK-1
            IF (WVN(I+1) .LE. KC) THEN
               T1 = T1 + 0.5 * (WVN(I+1)-WVN(I)) * (WVN(I)**2*BETA(I)
     &              *SK(I)+WVN(I+1)**2 * BETA(I+1) * SK(I+1))
            ENDIF
         ENDDO
         T1 = T1 * USTAR**2        
         T2 = (KAPPA * UD2)**2 / (LOG(D2 / Z0S))**2
         NEW_Z0T = D2 * EXP(-KAPPA * UD2 / (T1 + T2)**0.5)
         EPS = ABS(NEW_Z0T - Z0T)
         Z0T = Z0T*0.9 + NEW_Z0T*0.1

!--------------------------
!     STEP (g)
!--------------------------
         USTAR = USTAR*0.5 + KAPPA * U10 / LOG(10.0/Z0T)*0.5
C         PRINT *,"Z0T=",Z0T,"USTAR=",USTAR,"ALPHA=",Z0T*G/USTAR**2

         IF (EPS .LT. 1.0E-8) EXIT
      END DO
!      PRINT *,"Z0T=",Z0T,"USTAR=",USTAR,"ALPHA=",Z0T*G/USTAR**2
!      PRINT *,"EPS=",EPS
      ACH = Z0T*G/USTAR**2
!      WRITE(38,*) USTAR/SQRT(G*TANH(WTDP*KP)/KP),Z0T*G/USTAR**2,
!     &     (KAPPA/LOG(PI/KP/Z0T))**2

      RETURN
      END SUBROUTINE SOLVER_FD
!==============================================================================


!==============================================================================
      FUNCTION HEAVI(X)
      IMPLICIT NONE
      REAL HEAVI,X

      IF (X .LT. 0) THEN
         HEAVI = 0.0
      ELSE
         HEAVI = 1.0
      ENDIF
      RETURN

      END FUNCTION HEAVI
!==============================================================================

!==============================================================================
      SUBROUTINE FINDZ0T(Z0T,AW,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &     BETA,SK,WVN,REK,WTDP)

!---------------------------------------------------------
!
!     BY XUANTING HAO, THIS SUBROUTINE CALCULATES THE COEFFICIENT AW
!     BASED ON FIVE DIFFERENT MODELS.
!
!     IDX = 1 RMS MODEL
!     IDX = 2 GEOMETRY MODEL
!     IDX = 3 STEEPNESS DEPENDENT CHARNOCK MODEL
!     IDX = 4 WAVE KINEMATICS DEPENDENT MODEL
!     IDX = 5 COMBINED KINEMATICS STEEPNESS MODEL
!
!     07/16/2015   FIRST VERSION
!
!---------------------------------------------------------
      
      IMPLICIT NONE

      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,WTDP
      REAL,INTENT(IN)::BETA(*),SK(*),WVN(*),REK(*)
      REAL,INTENT(INOUT)::Z0T,AW

      REAL SIG,EPSAW,EPS,TD,T2D,K2D
      REAL TMPTD,TMPT2D
      REAL AWA,AWB,AWC,FDA,FDB,FDC
      INTEGER I,J,IT,ITMAX

      K2D = KD / 2
      ITMAX = 20000

      AWA = 0
      AWB = 100

      EPSAW = 1E-3
      EPS = 1E-4
      CALL FDELTA(TD,AWA,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &     BETA,SK,WVN,REK,WTDP)
!      PRINT *,TD,AWA,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD
!      PRINT *,BETA(1:NK)
!      PRINT *,"SK",SK(1:NK)
!      PRINT *,"WVN",WVN(1:NK)
!      STOP
      CALL FDELTA(T2D,AWA,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,K2D,
     &     BETA,SK,WVN,REK,WTDP)
      FDA = TD - T2D
      TMPTD = TD
      TMPT2D = T2D
!      PRINT *,"TD/USTAR**2=",TD/USTAR**2
!      PRINT *,"T2D/USTAR**2=",T2D/USTAR**2

      CALL FDELTA(TD,AWB,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &     BETA,SK,WVN,REK,WTDP)
      CALL FDELTA(T2D,AWB,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,K2D,
     &     BETA,SK,WVN,REK,WTDP)
      FDB = TD - T2D
!      PRINT *,"TD/USTAR**2=",TD/USTAR**2
!      PRINT *,"T2D/USTAR**2=",T2D/USTAR**2
!      STOP

      IF (FDA * FDB .GT. 0) THEN
         PRINT *,TMPTD,TMPT2D
         PRINT *,TD,T2D
         PRINT *,"ADJUST THE INITIAL VALUE OF AWA AND AWB!"
         STOP
      ENDIF

      DO I = 1, ITMAX
         AWC = 0.5 * (AWA + AWB)
         CALL FDELTA(TD,AWC,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &        BETA,SK,WVN,REK,WTDP)
         CALL FDELTA(T2D,AWC,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,K2D,
     &        BETA,SK,WVN,REK,WTDP)
         FDC = TD - T2D
         
         IF (ABS(FDC) .LE. EPS * ABS(TD) .OR.
     &        ABS(AWB - AWA)/2.0 .LE. EPSAW) THEN
            AW = AWC
            EXIT
         ENDIF

         IF (FDA * FDC .GT. 0) THEN
            AWA = AWC
         ELSE
            AWB = AWC
         ENDIF

         CALL FDELTA(TD,AWA,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &        BETA,SK,WVN,REK,WTDP)
         CALL FDELTA(T2D,AWA,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,K2D,
     &        BETA,SK,WVN,REK,WTDP)
         FDA = TD - T2D

         CALL FDELTA(TD,AWB,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &        BETA,SK,WVN,REK,WTDP)
         CALL FDELTA(T2D,AWB,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,K2D,
     &        BETA,SK,WVN,REK,WTDP)
         FDB = TD - T2D      
      ENDDO

      
      Z0T = D2 / EXP(KAPPA * UD2 / SQRT(-TD))
!      PRINT *,"AWA,AWB",AWA,AWB
C      PRINT *,"FDC",FDC
C      PRINT *,"TD,T2D",TD,T2D
C      PRINT *,"TD/USTAR**2",TD/USTAR**2
C      PRINT *,"AW=",AW      
      RETURN
      END SUBROUTINE FINDZ0T
!==============================================================================

!==============================================================================

      SUBROUTINE FDELTA(TAOT,AW,NK,G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,
     &     BETA,SK,WVN,REK,WTDP)

!---------------------------------------------------------
!
!     BY XUANTING HAO, THIS SUBROUTINE CALCULATES THE COEFFICIENT AW
!     BASED ON FIVE DIFFERENT MODELS.
!
!     IDX = 1 RMS MODEL
!     IDX = 2 GEOMETRY MODEL
!     IDX = 3 STEEPNESS DEPENDENT CHARNOCK MODEL
!     IDX = 4 WAVE KINEMATICS DEPENDENT MODEL
!     IDX = 5 COMBINED KINEMATICS STEEPNESS MODEL
!
!     07/16/2015   FIRST VERSION
!
!---------------------------------------------------------

      IMPLICIT NONE

      REAL,INTENT(INOUT)::TAOT,AW
      INTEGER,INTENT(IN):: NK
      REAL,INTENT(IN)::G,Z0S,UD2,USTAR,KAPPA,D2,KC,KD,WTDP
      REAL,INTENT(IN)::BETA(*),SK(*),WVN(*),REK(*)

      REAL SIG
      INTEGER I,J,IT,ITMAX

      TAOT = 0

!     ADD GRID SCALE STRESS
      DO I = 1, NK - 1
         IF (REK(I) .LE. KD) THEN
            TAOT = TAOT  + 0.5 * (REK(I+1)-REK(I)) * (REK(I)**2*BETA(I)
     &              *SK(I)+REK(I+1)**2 * BETA(I+1) * SK(I+1))  
         ENDIF
      ENDDO
!      PRINT *,"TAOT GS / USTAR**2",TAOT / USTAR**2
      TAOT = TAOT * (-USTAR**2)

!     ADD SUB-GRID SCALE STRESS
      CALL EFFAMP(SIG,NK,SK,REK,KD,KC,USTAR,G,KAPPA,
     &     4,WTDP)
      TAOT = TAOT - (KAPPA * UD2)**2 
     &      / (LOG(D2 / SQRT(Z0S**2 + (AW * SIG)**2 )))**2
!      PRINT *,"TAOT SGS / USTAR**2",(KAPPA * UD2)**2
!     &      / (LOG(D2 / SQRT(Z0S**2 + (AW * SIG)**2 )))**2/USTAR**2

      END SUBROUTINE FDELTA

!=============================================================================



!=============================================================================
      SUBROUTINE EFFAMP_BENCHMARK(EETABM,SK,REK,NK,G,D2,KAPPA,Z0S,USTAR)

      IMPLICIT NONE
      REAL Z0T

      INTEGER NK
      REAL SK(*),REK(*),EETABM
      REAL ACH,G,U10,FETCH,NUA,KAPPA,D2,KP

      REAL ZC(NK),BETA(NK),USTAR,ZETAC,HEAVI
      REAL Z0S,UD2,T1,T2,NEW_Z0T
      REAL EPS,PI
      INTEGER I

      PI = ACOS(-1.0)
      DO I = 1, NK
         ZC(I) = Z0T * EXP((KAPPA / USTAR)
     &        * SQRT(G / REK(I)))
         ZETAC = ZC(I) * REK(I)
         BETA(I) = PI * ZETAC * (LOG(0.281  / ZETAC))**4
     &        * HEAVI(0.281  - ZETAC) / KAPPA**2 +
     &        2*LOG(0.281  / ZETAC)
C     WRITE(35,901) REK(I),BETA(I)
      ENDDO


      END SUBROUTINE EFFAMP_BENCHMARK
!=============================================================================

!=============================================================================
      SUBROUTINE EFFAMP(EFFETA,NK,SK,REK,KD,KC,USTAR,G,KAPPA,
     &     IDX,WTDP)

      IMPLICIT NONE

!---------------------------------------------------------
!
!     BY XUANTING HAO, THIS SUBROUTINE CALCULATES THE EFFECTIVE AMPLITUDE
!     BASED ON FIVE DIFFERENT MODELS.
!
!     IDX = 1 RMS MODEL
!     IDX = 2 GEOMETRY MODEL
!     IDX = 3 STEEPNESS DEPENDENT CHARNOCK MODEL
!     IDX = 4 WAVE KINEMATICS DEPENDENT MODEL
!     IDX = 5 COMBINED KINEMATICS STEEPNESS MODEL
!
!     07/16/2015   FIRST VERSION
!
!---------------------------------------------------------

      REAL,INTENT(INOUT)::EFFETA
      INTEGER,INTENT(IN)::NK
      REAL,INTENT(IN)::SK(*),REK(*),USTAR,G,KAPPA,KC,KD,WTDP
      INTEGER,INTENT(IN)::IDX
      
      INTEGER I,IU,IL
      REAL TWOPI,TMP1,TMP2,MYTANH

      TWOPI = 2 * ACOS(-1.0)
      EFFETA = 0
      IF (KC .LT. KD) THEN
         EFFETA = 0 
         RETURN
      ENDIF

!     COMPUTE THE UPPER LIMIT OF INTEGRATION
      IU = NK
      DO WHILE (REK(IU) .GT. KC) 
         IU = IU - 1
      ENDDO
      IF (IU .LE. 1) THEN
         PRINT *,"KC=",KC
         PRINT *,"KD=",KD
         STOP
      ENDIF

!     COMPUTE THE LOWER LIMIT OF INTEGRATION
      IL = 1
      DO WHILE (REK(IL) .LT. KD)
         IL = IL + 1
      ENDDO
      IF (IL .GE. NK) THEN
         PRINT *,"KC=",KC
         PRINT *,"KD=",KD
         STOP
      ENDIF

      SELECT CASE (IDX)
         CASE (1)
            DO I = IL, IU - 1               
               EFFETA = EFFETA + 0.5 * (REK(I+1) - REK(I)) * 
     &              (SK(I+1) + SK(I))
            ENDDO
            EFFETA = SQRT(EFFETA)
         CASE (2)
            DO I = IL, IU - 1
               EFFETA = EFFETA + 0.5 * (REK(I+1) - REK(I)) *
     &              (SK(I+1) * REK(I+1) + SK(I) * REK(I))
            ENDDO
         CASE (3)
            DO I = IL, IU - 1
               EFFETA = EFFETA + (REK(I+1) - REK(I)) *
     &              (SK(I+1) * REK(I+1)**2 + SK(I) * REK(I)**2)
            ENDDO
            EFFETA = SQRT(EFFETA) * USTAR**2 / TWOPI / G
         CASE (4) 
            DO I = IL, IU - 1
               TMP1 = -(2 * KAPPA / USTAR) * SQRT(G
     &              * MYTANH(REK(I+1) , WTDP) / REK(I+1))
               TMP2 = -(2 * KAPPA / USTAR) * SQRT(G 
     &              * MYTANH(REK(I) , WTDP) / REK(I))
               IF (TMP1 .GT. -40.0) THEN
                  EFFETA = EFFETA + 0.5 * (REK(I+1) - REK(I)) *
     &                 SK(I+1) * EXP(TMP1)
               ENDIF
               IF (TMP2 .GT. -40.0) THEN
                  EFFETA = EFFETA + 0.5 * (REK(I+1) - REK(I)) *
     &                 SK(I) * EXP(TMP2)
               ENDIF
            ENDDO
            EFFETA = SQRT(EFFETA)
         CASE (5)
            DO I = IL, IU - 1
               TMP1 = -(2 * KAPPA / USTAR) * SQRT(G 
     &              * MYTANH(REK(I+1) , WTDP) / REK(I+1))
               TMP2 = -(2 * KAPPA / USTAR) * SQRT(G 
     &              * MYTANH(REK(I) , WTDP) /REK(I))
               IF (TMP1 .GT. -40.0) THEN
                  EFFETA = EFFETA + 0.5 * (REK(I+1) - REK(I)) *
     &                 * REK(I+1)**2 * SK(I+1) * EXP(TMP1)
               ENDIF
               IF (TMP2 .GT. -40.0) THEN
                  EFFETA = EFFETA + 0.5 * (REK(I+1) - REK(I)) *
     &                 * REK(I)**2 * SK(I) * EXP(TMP2)
               ENDIF
            ENDDO
            EFFETA = SQRT(EFFETA) * USTAR**2 / TWOPI / G
         CASE DEFAULT
            PRINT *,"WRONG VALUE FOR IDX IN EFFAMP!"
            STOP
      END SELECT
      
      RETURN
      END SUBROUTINE EFFAMP

!=============================================================================

