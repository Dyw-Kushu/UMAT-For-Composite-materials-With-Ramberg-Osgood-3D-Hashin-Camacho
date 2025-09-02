      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1        RPL, DDSDDT, DRPLDE, STRAN, DSTRAN, TEMP, DTEMP,
     2        PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
     3        PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT,
     4        DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      IMPLICIT NONE
      CHARACTER*80 CMNAME
      INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS
      INTEGER NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      DOUBLE PRECISION STRESS(*), STATEV(*), DDSDDE(*)
      DOUBLE PRECISION SSE, SPD, SCD, RPL
      DOUBLE PRECISION DDSDDT(*), DRPLDE(*), STRAN(*), DSTRAN(*)
      DOUBLE PRECISION TEMP, DTEMP, PREDEF(*), DPRED(*)
      DOUBLE PRECISION PROPS(*), COORDS(*), DROT(3,3)
      DOUBLE PRECISION PNEWDT, CELENT, DFGRD0(3,3), DFGRD1(3,3)

      INTEGER i
      DOUBLE PRECISION E0(6), SIG0(6), NEXP(6)
      DOUBLE PRECISION sigma_new(6), eps_total(6), d_sigma_d_eps(6)
      DOUBLE PRECISION Xt, Xc, Yt, Yc, S12s, S13s, S23s
      DOUBLE PRECISION d_mat_tens, d_mat_comp, d_fib_tens, d_fib_comp
      DOUBLE PRECISION S11,S22,S33,S12,S13,S23,x
      DOUBLE PRECISION xi_ft, xi_fc, xi_Mt, xi_Mc
      DOUBLE PRECISION xi_hist_ft, xi_hist_fc, xi_hist_Mt, xi_hist_Mc
      DOUBLE PRECISION factor(6), Eeff(6)
      LOGICAL fib_tens_fail, fib_comp_fail, mat_tens_fail, mat_comp_fail
      DOUBLE PRECISION ONE, ZERO
      PARAMETER (ONE=1.0D0, ZERO=0.0D0)

C --- PROPS layout ---
C 1..6   : E0(1..6)
C 7..12  : SIG0(1..6)
C 13..18 : NEXP(1..6)
C 19..25 : Xt, Xc, Yt, Yc, S12, S13, S23
C 26     : d_mat_tens
C 27     : d_mat_comp
C 28     : d_fib_tens
C 29     : d_fib_comp
C STATEV(1:6)  : current stress
C STATEV(7)   : xi_hist_ft
C STATEV(8)   : xi_hist_fc
C STATEV(9)   : xi_hist_Mt
C STATEV(10)  : xi_hist_Mc
C STATEV(11)  : failure flag (0/1)
C STATEV(12:17): stiffness reduction factors for directions 1-6

      IF (NPROPS .LT. 29) THEN
         WRITE(*,*) 'UMAT ERROR: PROPS need >=29, got', NPROPS
      END IF

      DO i = 1,6
         E0(i)   = PROPS(i)
         SIG0(i) = PROPS(6 + i)
         NEXP(i) = PROPS(12 + i)
      END DO

      Xt   = PROPS(19)
      Xc   = PROPS(20)
      Yt   = PROPS(21)
      Yc   = PROPS(22)
      S12s = PROPS(23)
      S13s = PROPS(24)
      S23s = PROPS(25)

      d_mat_tens = PROPS(26)
      d_mat_comp = PROPS(27)
      d_fib_tens = PROPS(28)
      d_fib_comp = PROPS(29)

C --- strain input
      DO i = 1,6
         eps_total(i) = STRAN(i)
      END DO
      
C --- initialize STATEV
      IF (KINC .EQ. 1) THEN
          DO i=1,11
            STATEV(i) = ZERO
          END DO
          DO i=12,18
            STATEV(i) = ONE
          END DO
      END IF  
      
C --- read historical factors
      DO i = 1,6
          factor(i)= STATEV(11+i)
      END DO
      
C --- Ramberg-Osgood
      DO i = 1,6
         Eeff(i) = E0(i) * factor(i)
         x = (Eeff(i)*eps_total(i))/SIG0(i)
         sigma_new(i) = (Eeff(i)*eps_total(i)) /
     &                  ( (ONE+ABS(x)**NEXP(i))**(ONE/NEXP(i)) )
         d_sigma_d_eps(i) = Eeff(i) /
     &                  ( (ONE+ABS(x)**NEXP(i))**(ONE+ONE/NEXP(i)) )
      END DO  

C --- stress components
      S11 = sigma_new(1)
      S22 = sigma_new(2)
      S33 = sigma_new(3)
      S12 = sigma_new(4)
      S13 = sigma_new(5)
      S23 = sigma_new(6)

C --- Hashin failure
      xi_ft = ZERO
      xi_fc = ZERO
      xi_Mt = ZERO
      xi_Mc = ZERO

      IF (S11 .GE. ZERO) THEN
         xi_ft = (S11/Xt)**2 + (S12/S12s)**2 + (S13/S13s)**2
      ELSE
         xi_fc = (S11/Xc)**2
      END IF

      IF ( (S22+S33) .GE. ZERO ) THEN
         xi_Mt = ( (S22+S33)**2 )/(Yt*Yt)
     &          + (S23**2 - S22*S33)/(S23s*S23s)
     &          + (S12**2 + S13**2)/(S12s*S12s)
      ELSE
         xi_Mc = ( ( (Yc/(2.0D0*S23s))**2 - ONE ) * (S22+S33)/(Yc*Yc))
     &          + ( (S22+S33)**2 )/(4.0D0*S23s*S23s)
     &          + (S23**2 - S22*S33)/(S23s*S23s)
     &          + (S12**2 + S13**2)/(S12s*S12s)
      END IF

C --- read history xi
      xi_hist_ft = STATEV(7)
      xi_hist_fc = STATEV(8)
      xi_hist_Mt = STATEV(9)
      xi_hist_Mc = STATEV(10)

C --- update history with maximum
      xi_hist_ft = MAX(xi_hist_ft, xi_ft)
      xi_hist_fc = MAX(xi_hist_fc, xi_fc)
      xi_hist_Mt = MAX(xi_hist_Mt, xi_Mt)
      xi_hist_Mc = MAX(xi_hist_Mc, xi_Mc)

C --- failure flags from historical maxima
      fib_tens_fail = (xi_hist_ft .GE. ONE)
      fib_comp_fail = (xi_hist_fc .GE. ONE)
      mat_tens_fail = (xi_hist_Mt .GE. ONE)
      mat_comp_fail = (xi_hist_Mc .GE. ONE)
      
C --- Camacho stiffness reduction
C --- set factors according to historical flags
      IF (fib_tens_fail) THEN
         factor(1) = d_fib_tens
      ELSEIF (fib_comp_fail) THEN
         factor(1) = d_fib_comp
      END IF

      IF (mat_tens_fail) THEN
         factor(2) = d_mat_tens
         factor(3) = d_mat_tens
         factor(4) = d_mat_tens
         factor(6) = d_mat_tens
      END IF
      IF (mat_comp_fail) THEN
         factor(2) = d_mat_comp
         factor(3) = d_mat_comp
         factor(4) = d_mat_comp
         factor(6) = d_mat_comp
      END IF

C --- now recompute sigma_new and d_sigma_d_eps using Eeff = E0 * factor
      DO i = 1,6
         Eeff(i) = E0(i) * factor(i)
         x = (Eeff(i)*eps_total(i))/SIG0(i)
         sigma_new(i) = (Eeff(i)*eps_total(i)) /
     &                  ( (ONE+ABS(x)**NEXP(i))**(ONE/NEXP(i)) )
         d_sigma_d_eps(i) = Eeff(i) /
     &                  ( (ONE+ABS(x)**NEXP(i))**(ONE+ONE/NEXP(i)) )
      END DO      

C --- update STRESS
      DO i = 1,6
         STRESS(i) = sigma_new(i)
      END DO

C --- tangent DDSDDE (diagonal approximation)
      DO i = 1,36
         DDSDDE(i) = ZERO
      END DO
      DO i = 1,6
         DDSDDE(i+(i-1)*6) = d_sigma_d_eps(i)
      END DO

C --- update STATEV
      DO i = 1,6
          STATEV(i) = STRESS(i)
      END DO
      
      STATEV(7)  = xi_hist_ft
      STATEV(8)  = xi_hist_fc
      STATEV(9)  = xi_hist_Mt
      STATEV(10) = xi_hist_Mc
      STATEV(11)= ZERO
      
      IF (fib_tens_fail .OR. fib_comp_fail .OR.
     &    mat_tens_fail .OR. mat_comp_fail) THEN
         STATEV(11)=ONE
      END IF
     
      DO i = 1,6
          STATEV(11+i) = factor(i)
      END DO
      
C --- energy terms zero
      SSE=ZERO
      SPD=ZERO
      SCD=ZERO
      RPL=ZERO

      RETURN
      END
