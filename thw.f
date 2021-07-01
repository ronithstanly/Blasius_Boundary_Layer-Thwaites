!--------------------------------------------------------------------------------------------------------------------------------------------
!******************* Fortran Program for solving 2D, Laminar Boundary Layer over Flat Plate using Thwaites' Method **************************
!Date: 23 May, 2017
!Authors: Prof.Alexey Nikolaevych Kudryavtsev and Ronith Stanly
!Affiliation: Lab7, ITAM, Russia
!--------------------------------------------------------------------------------------------------------------------------------------------     
      MODULE Prec
      INTEGER, PARAMETER :: knd = 8					!Double Precision
      END MODULE Prec

      PROGRAM Thwaites
      USE Prec
      IMPLICIT REAL(knd) (A-H,O-Z)
      INTEGER, PARAMETER :: N = 100					!No. of 'x' Divisions
      REAL(knd), PARAMETER :: XL = 1.0_knd 				!Length of plate
      REAL(knd), PARAMETER :: Visc = 1.0_knd/200000.0_knd 		!Kinematic Viscosity
      REAL(knd), DIMENSION(0:N) :: X, Ue, dUe, Thet

      DO i=0,N 								!Grid loop
        X(i) = i*XL/N
        Ue(i) = 1.d0
      END DO

      dUe(0) = (Ue(1)-Ue(0))/(X(1)-X(0)) 				!One-sided Differencing for derivative of velocity at first point
      DO i=1,N-1
       dUe(i) = (Ue(i+1)-Ue(i-1))/(X(i+1)-X(i-1)) 			!Central Differencing for all in-between derivatives
      END DO
      dUe(N) = (Ue(N)-Ue(N-1))/(X(N)-X(N-1)) 				!One-sided Differencing for derivative of velocity at last point

      Thet(0) = 0.0_knd							!Momentum Thickness=0, at x=0

      S = 0.0_knd
      DO i=1,N
        S = S+0.5_knd*(Ue(i-1)**5+Ue(i)**5)*(X(i)-X(i-1))		!Numerical Integration of Eq. (4.76a) in [1] by Trapezoidal Rule
        Thet(i) = SQRT(Thet(0)**2*Ue(0)**6+0.45_knd*Visc*S)/Ue(i)**3	!Momentum Thickness:- Eq. (4.76a) in [1]
      END DO

      OPEN (10,FILE='thwaites.dat')					!Create output data file
        DO i=0,N
          Rth = Ue(i)*Thet(i)/Visc					!Momentum Thickness Reynold's Number:- Eq.(4.78c) in [1]
          CL = Thet(i)**2*dUe(i)/Visc					!Thwaites' Parameter-Lambda
          IF (CL > 0.0_knd) THEN
            H  = 2.61_knd-CL*(3.75_knd-CL*5.24_knd)			!Shape Factor:- Eq.(4.78a) in [1]
            Cf = 0.225_knd+CL*(1.61_knd-CL*(3.75_knd-CL*5.24_knd))	!Coefficient of Friction:- Eq.(4.78a) in [1]
            Cf = 2.0_knd*Cf/Rth						!Coefficient of Friction:- Eq.(4.78a) in [1]
          ELSE
            H  = 2.472_knd+0.0147_knd/(CL+0.107_knd)			!Shape Factor:- Eq.(4.78b) in [1]
            Cf = 0.225_knd+CL*(1.472_knd+0.0147_knd/(CL+0.107_knd))	!Coefficient of Friction:- Eq.(4.78b) in [1]
            Cf = 2.0_knd*Cf/Rth						!Coefficient of Friction:- Eq.(4.78b) in [1]
          END IF
          WRITE (10,*) X(i),Ue(i),Thet(i),Rth,H,Cf			!Writes output data to file
        END DO
      CLOSE (10)

      END PROGRAM Thwaites
!--------------------------------------------------------------------------------------------------------------------------------------------
!Reference:[1] Turner Cebeci and Peter Bradshaw, "Physical and Computational Aspects of Convective Heat Transfer",1984
!--------------------------------------------------------------------------------------------------------------------------------------------
