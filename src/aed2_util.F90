!###############################################################################
!#                                                                             #
!# aed2_util.F90                                                               #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created June  2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_util
!-------------------------------------------------------------------------------
! aed2_util --- shared utility functions for aed modules
!
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC find_free_lun, qsort
   PUBLIC aed2_gas_piston_velocity, aed2_oxygen_sat, aed2_n2o_sat, exp_integral
   PUBLIC aed2_bio_temp_function,fTemp_function
   PUBLIC PO4AdsorptionFraction, in_zone_set
   PUBLIC water_viscosity
!


!===============================================================================
CONTAINS



!###############################################################################
INTEGER FUNCTION find_free_lun()
!-------------------------------------------------------------------------------
! find a free logical unit number
!-------------------------------------------------------------------------------
!LOCALS
    INTEGER :: lun
    LOGICAL :: opend
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lun = 10,99
      inquire( unit = lun, opened = opend )
      IF ( .not. opend ) THEN
         find_free_lun = lun
         RETURN
      ENDIF
   ENDDO

   find_free_lun = -1
END FUNCTION find_free_lun
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_gas_piston_velocity(wshgt,wind,tem,sal,LA,schmidt_model)
!-------------------------------------------------------------------------------
! Atmospheric-surface water exchange piston velocity for O2, CO2 etc
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(IN) :: wshgt,wind
   AED_REAL,INTENT(in) :: tem,sal
   AED_REAL,INTENT(in),OPTIONAL :: LA
   INTEGER,INTENT(in),OPTIONAL  :: schmidt_model
!
!LOCALS
   ! Temporary variables
   AED_REAL :: schmidt,k_wind,k_flow,temp,salt,hgtCorrx,a,x
   INTEGER  :: schmidt_model_l
   ! Parameters
   AED_REAL,PARAMETER :: roughlength = 0.000114  ! momn roughness length (m)
!
!-------------------------------------------------------------------------------
!BEGIN

   !-----------------------------------------------
   ! Decide on Sc equation to apply
   schmidt_model_l = 2 !default
   IF (PRESENT(schmidt_model)) schmidt_model_l = schmidt_model

   ! Adjust the windspeed if the sensor height is not 10m
   hgtCorrx =  LOG(10.00 / roughLength) / LOG(wshgt / roughLength)

   !-----------------------------------------------
   ! Compute k_wind
   IF (PRESENT(LA)) THEN

      ! New option for the calculation of k_wind. Note that this has a
      ! "lake area" (LA) variable included in it.

      ! Valchon & Prairie 2013: The ecosystem size and shape dependence of gas transfer
      !                              velocity versus wind speed relationships in lakes
      ! k600 = 2.51 (±0.99) + 1.48 (±0.34) · U10 + 0.39 (±0.08) · U10 · log10 LA

      k_wind = 2.51 + 1.48*wind*hgtCorrx  +  0.39*wind*hgtCorrx*log10(LA)

   ELSE
      temp=tem
      salt=sal
      IF (temp < 0.0)       temp = 0.0; IF (temp > 38.0)      temp = 38.0
      IF (salt < 0.0)       salt = 0.0; IF (salt > 75.0)      salt = 75.0

      ! Schmidt, Sc
      ! control value : Sc = 590 at 20°C and 35 psu
      schmidt = 590.

      SELECT CASE (schmidt_model_l)
      CASE (1)
         schmidt = (0.9 + 0.1*salt/35.0)*(1953.4+temp*(-128.0+temp*(3.9918-temp*0.050091)))
      CASE (2)
         schmidt = (0.9 + salt/350.0)
         schmidt = schmidt * (2073.1 -125.62*temp +3.6276*temp*temp - 0.043219*temp*temp*temp)
      CASE (3)
         ! http://www.geo.uu.nl/Research/Geochemistry/kb/Knowledgebook/O2_transfer.pdf
         schmidt = (1.0 + 3.4e-3*salt)
         schmidt = schmidt * (1800.6 -120.1*temp +3.7818*temp*temp - 0.047608*temp*temp*temp)
      CASE (4)
         ! CH4 one from Arianto Santoso <abs11@students.waikato.ac.nz>
         schmidt = 2039.2 - (120.31*temp) + (3.4209*temp*temp) - (0.040437*temp*temp*temp)
         schmidt = schmidt / 600
       CASE (5)
         ! CH4 from Sturm et al. 2014 (ex Wanninkhof, 1992)
         schmidt = 1897.8 - (114.28*temp) + (3.2902*temp*temp) - (0.039061*temp*temp*temp)
       CASE (6)
         ! N2O from Sturm et al. 2014 (ex Wanninkhof, 1992)
         schmidt = 2055.6 - (137.11*temp) + (4.3173*temp*temp) - (0.054350*temp*temp*temp)
      END SELECT

      ! Gas transfer velocity (cm/hr)
      ! k = a u^2 (Sc/600)^-x
      ! This parameterization assumes 10m windspeed, and must be scaled by hgtCorrx
      a = 0.31
      x = 0.50
      IF( wind*hgtCorrx <3.) x = 0.66
      k_wind = a * (wind*hgtCorrx)*(wind*hgtCorrx) * (schmidt/600.0)**(-x)
   ENDIF
   ! convert to m/s
   k_wind = k_wind / 3.6e5

   !-----------------------------------------------
   ! Compute k_flow
   k_flow = zero_ !(vel**0.5)*(depth**(-0.5))

   !-----------------------------------------------
   ! piston velocity is the sum due to flow and wind
   aed2_gas_piston_velocity = k_flow + k_wind

END FUNCTION aed2_gas_piston_velocity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_oxygen_sat(salt,temp)
!-------------------------------------------------------------------------------
!  Calculated saturated oxygen concentration at salinity and temperature
! Taken from Riley and Skirrow (1974)
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: salt,temp
!
!LOCALS
   AED_REAL :: Tabs
   AED_REAL :: buf1, buf2, buf3, sol_coeff
!
!-------------------------------------------------------------------------------
!BEGIN
   buf1 = zero_ ; buf2 = zero_ ; buf3 = zero_ ; sol_coeff = zero_

   Tabs = temp + 273.15
   buf1 = -173.4292 + 249.6339 * 100.0 / Tabs + 143.3483 * LOG(Tabs/100.0)
   buf2 = 21.8492 * Tabs / 100.0
   buf3 = salt * (-0.033096 + 0.014259 * Tabs / 100.0 - 0.0017 * (Tabs / 100.0)**2.0)
   sol_coeff = buf1 - buf2 + buf3

   aed2_oxygen_sat = 1.42763 * exp(sol_coeff) !in g/m3

   !Convert to mmol/m3
   aed2_oxygen_sat = (aed2_oxygen_sat / 32.) * 1e3
END FUNCTION aed2_oxygen_sat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_n2o_sat(salt,temp)
!-------------------------------------------------------------------------------
!  gsw_N2Osol_SP_pt                            solubility of N2O in seawater
!
!  USAGE:
!   N2Osol = gsw_N2Osol_SP_pt(SP,pt)
!
!  DESCRIPTION:
!   Calculates the nitrous oxide, N2O, concentration expected at equilibrium
!   with air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
!   including saturated water vapor  This function uses the solubility
!   coefficients as listed in Hamme and Emerson (2004).
!
!   Note that this algorithm has not been approved by IOC and is not work
!   from SCOR/IAPSO Working Group 127. It is included in the GSW
!   Oceanographic Toolbox as it seems to be oceanographic best practice.
!
!  INPUT:
!   salt  =  Practical Salinity  (PSS-78)                         [ unitless ]
!   temp  =  potential temperature (ITS-90) referenced               [ deg C ]
!          to one standard atmosphere (0 dbar).
!
!  OUTPUT:
!   N2Osol = solubility of N2O                                      [ mol/L ]
!
!  AUTHOR:  Rich Pawlowicz, Paul Barker and Trevor McDougall
!                                                       [ help@teos-10.org ]
!
!  VERSION NUMBER: 3.05 (27th January 2015)
!
!  REFERENCES:
!   IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
!    seawater - 2010: Calculation and use of thermodynamic properties.
!    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
!    UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
!
!   Weiss, R.F. and B.A. Price, 1980: Nitrous oxide solubility in water and
!    seawater. Mar. Chem., 8, 347-359.
!
!   The software is available from http://www.TEOS-10.org
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: salt,temp
!
!LOCALS
   AED_REAL :: x, y, y_100, pt68, ph2odP
   AED_REAL :: a0,a1,a2,a3,b1,b2,b3,m0,m1,m2,m3
!
!-------------------------------------------------------------------------------
!BEGIN

  x = salt     !  Note that salinity argument is Practical Salinity, this is
               !  beacuse the major ionic components of seawater related to Cl
               !  are what affect the solubility of non-electrolytes in seawater

  pt68 = temp*1.00024 ! pt68 is the potential temperature in degress C on
                      ! the 1968 International Practical Temperature Scale IPTS-68.
  y = pt68 + 273.15
  y_100 = y*1e-2

  !  The coefficents below are from Table 2 of Weiss and Price (1980)
  a0 = -165.8806
  a1 =  222.8743
  a2 =  92.0792
  a3 = -1.48425
  b1 = -0.056235
  b2 =  0.031619
  b3 = -0.0048472

  m0 = 24.4543
  m1 = 67.4509
  m2 = 4.8489
  m3 = 0.000544

  ph2odP = exp(m0 - m1*100.0/y - m2*log(y_100) - m3*x) !  Moist air correction at 1 atm.

  !aed2_n2o_sat [mol/L] = (exp(a0 + a1*100.0/y + a2*log(y_100) + a3*y_100 + x*(b1 + y_100*(b2 + b3*y_100))))/(1.-ph2odP);
  aed2_n2o_sat = (exp(a0 + a1*100.0/y + a2*log(y_100) + a3*y_100*y_100 + x*(b1 + y_100*(b2 + b3*y_100))))/(1.-ph2odP);

  !Convert to mmol/m3
  aed2_n2o_sat = aed2_n2o_sat * 1e6

END FUNCTION aed2_n2o_sat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION exp_integral(inp) RESULT(E_ib)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: inp
!
!LOCALS
   AED_REAL  :: E_ib !-- Outgoing
   INTEGER   :: j
   AED_REAL  :: ff
!
!-------------------------------------------------------------------------------
!BEGIN
   ff = -1e-9
   IF(ABS(inp-10.0) < 12.0) THEN
     IF(inp==0.0) THEN
       E_ib = inp
     ELSE
       j  = 10+2*IABS(INT(inp))
       ff = 1.0/(REAL(j+1)**2.0)
       DO WHILE(j/=0)
         ff = (ff*REAL(j)*inp+1.0)/REAL(j*j)
         j  = j-1
       ENDDO
       ff   = ff*inp+LOG(1.781072418*ABS(inp))
       E_ib = ff
     ENDIF
   ELSE
     j = 5 + 20 / IABS(INT(inp))
     ff = inp
     DO WHILE(j/=0)
       ff = (1.0/(1.0/ff-1.0/REAL(j)))+inp
       j = j-1
     ENDDO
     ff  = EXP(inp)/ff
     E_ib = ff
   ENDIF

END FUNCTION exp_integral
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_bio_temp_function(numg, theta, T_std, T_opt, T_max, aTn, bTn, kTn, name)
!-------------------------------------------------------------------------------
! Numerically solver for continuos temperature function
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)       :: numg        ! Number of groups
   AED_REAL,INTENT(in)      :: theta(:)
   AED_REAL,INTENT(inout)   :: T_std(:), T_opt(:), T_max(:)
   AED_REAL,INTENT(out)     :: aTn(:), bTn(:), kTn(:)
   CHARACTER(64),INTENT(in) :: name(:)
!
!LOCALS
   AED_REAL :: Ts     ! Min. temperature where fT(Ts)=I (usually 1)
   AED_REAL :: To     ! Optimum temperature where d(fT(To))/dT=0
   AED_REAL :: Tm     ! Maximum temperature where fT(Tm)=0
   AED_REAL :: in     ! Constant for fT(Ts)=in
   AED_REAL :: v      ! Constant v
   AED_REAL :: k,a,b  ! Model constants
   AED_REAL :: G      ! Function fT()
   AED_REAL :: devG       ! Derivative of fT()
   AED_REAL :: a0,a1,a2   ! Dummies
   AED_REAL :: tol        ! Tolerance
   INTEGER :: group       ! Group counter
   INTEGER :: i           ! Counters
   AED_REAL,PARAMETER :: t20=20.0
   LOGICAL,PARAMETER :: curvef=.true. ! T : f(T)=v**(T-20) at T=Tsta
                                      ! F : f(T) = 1 at T=Tsta

   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: value
!
!-------------------------------------------------------------------------------
!BEGIN
    write(*,"('Estimating temperature functions for phytoplankton - ')")
    write(*,"(' Temperature function of the form :',/)")
    write(*,"('    fT = v^(T-20)-v^(k(T-a))+b',/)")

    tol   = 0.05

    DO group=1,numg

      ! Set the constants for the correct group
      v = theta(group)

      IF(v < 1.01) THEN
        print "(/,2X,'WARNING: theta_growth for group ',I2,' < 1.01',/)",group
      ENDIF

      Tm = T_max(group)
      Ts = T_std(group)
      To = T_opt(group)


      IF (Ts<0.0 .AND. To<0.0 .AND. Tm<0.0) THEN
        ! The user inputs the values of kTn, aTn and bTn directly
        kTn(group) = -Ts
        bTn(group) = -Tm
        aTn(group) = -To

        ALLOCATE(value(401,1))
        ! Calculate the temperature function using 0.1 deg C intervals

        DO i = 0,400
          b = REAL(i)/10.0
          value(i+1,1) = v**(b-20) - v**(kTn(group) * (b - aTn(group))) + bTn(group)
        ENDDO

        ! Find the values of Tsta, T_opt and T_max from the temp function
        a=0.0
        DO i=1,SIZE(value,1)
          b=REAL(i-1)/10.0
          IF(value(i,1)>0.0) THEN
            T_max(group) = b
          ENDIF
          IF(value(i,1)>a) THEN
            T_opt(group) = b
            a=value(i,1)
          ENDIF
          IF(value(i,1)>v**(b-20)-tol .and. value(i,1)<v**(b-20)+tol) THEN
            T_std(group) = b
          ENDIF
        ENDDO
        DEALLOCATE(value)


      ELSE
        in = 1.0
        a0 = v**(Ts-t20)
        a1 = v**(To-t20)
        a2 = v**(Tm-t20)

        ! Perform the iteration to find the constants.
        ! First approximation of k.
        k = 6.0
        i = 0
        G = tol + 1.0
        ! Do the iterations until -tol < G < tol
        DO WHILE((G <= -tol) .OR. (G >= tol))
          i=i+1
          IF(i==100) THEN  ! Increases the tolerance if more than 100
            i=0            ! iterations are performed.
            tol=tol+0.01
          ENDIF
          IF(curvef) THEN
            ! Use the condition f(T)=v**(T-20) at T=Tsta
            G = k * v**(k * To) * a2 - a1 * (v**(k * Tm) - v**(k * Ts))
            devG = v**(k * To) * a2 * (in + k * To * log(v)) - a1 * log(v) &
              * (Tm * v**(k * Tm) - Ts * v**(k * Ts))
          ELSE
            ! Use the condition f(T)=1 at T=Tsta
            G = k * v**(k * To) * (a0 - a2 - in) - a1 * (v**(k * Ts) &
              - v**(k * Tm))
            devG = (a0 - a2 - in) * v**(k * To) * (in + k * To * log(v)) - a1 &
              * log(v) * (Ts * v**(k * Ts) - Tm * v**(k * Tm))
          ENDIF
          ! Find the next iteration of k
          k = k - G / devG
        ENDDO

        ! Get the remaining model constants
        IF(k/=0.0) THEN
          a=-log(a1/(k*v**(k*To)))/(k*log(v))
          IF(curvef) THEN
            b=v**(k*(Ts-a))
          ELSE
            b=in+v**(k*(Ts-a))-a0
          ENDIF
        ELSE
          a=0.0
          b=0.0
        ENDIF

        ! Set the model constants to the calculated values
        kTn(group) = k
        aTn(group) = a
        bTn(group) = b
      ENDIF

      IF (kTn(group) < 0.1 .AND. bTn(group) > 100.0) THEN
            PRINT *,'Cannot solve for fT for: ', name(group)
            STOP
      ENDIF

    ENDDO

END SUBROUTINE aed2_bio_temp_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION fTemp_function(method,T_max,T_std,theta,aTn,bTn,kTn,temp) RESULT(fT)
!-------------------------------------------------------------------------------
! Generic temperature function for phytoplankton and zooplankton taking into
! account a decrease in production above T_opt.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: method
   AED_REAL,INTENT(in) :: T_max, T_std,theta,aTn,bTn,kTn
   AED_REAL,INTENT(in) :: temp  ! Temperature
!
!LOCALS
   AED_REAL  :: fT        !-- Value of the temperature function
   AED_REAL,PARAMETER  :: tp = 20.0
!
!-------------------------------------------------------------------------------
!BEGIN
   fT = one_

   IF ( method /= 1 ) RETURN

   IF (temp > T_max) THEN
       fT = zero_
   ELSEIF ( temp < T_std ) THEN
       IF (ABS(temp-tp) > 1+MINEXPONENT(temp)/2) THEN
         fT = theta**(temp-tp)
       ENDIF
   ELSE
      IF (ABS(temp-tp) > 1 + MINEXPONENT(temp)/2 .AND. &
          ABS((kTn*(temp-aTn)) + bTn) > 1 + MINEXPONENT(temp)/2) THEN
        fT = theta**(temp-tp) - theta**(kTn*(temp - aTn)) + bTn
      ENDIF
   ENDIF
END FUNCTION fTemp_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# A fortran implementation of the quicksort algorithm.                        #
!#                                                                             #
!###############################################################################
RECURSIVE SUBROUTINE qsort(RA,IA,start,end)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: RA(:)
   INTEGER,INTENT(inout) :: IA(:)
   INTEGER,INTENT(in) :: start,end
!
!LOCALS
  INTEGER :: p, l, r
!
!-------------------------------------------------------------------------------
!BEGIN
  IF ( start .LT. end ) THEN
     l=start+1
     r=end
     p = IA(start);

     DO WHILE(l<r)
        IF (cmp(IA(l), p) .LE. 0) THEN
           l=l+1;
        ELSEIF (cmp(IA(r), p) .GE. 0) THEN
           r=r-1
        ELSE
           CALL swap(IA(l), IA(r))
        ENDIF
     ENDDO
     IF (cmp(IA(l), p) .LT. 0 ) THEN
        CALL swap(IA(l), IA(start))
        l=l-1
     ELSE
        l=l-1
        CALL swap(IA(l), IA(start))
     ENDIF

     CALL qsort(RA,IA,start,l)
     CALL qsort(RA,IA,r,end)
  ENDIF

CONTAINS

   !############################################################################
   SUBROUTINE swap(a, b)
   !----------------------------------------------------------------------------
     INTEGER,intent(inout) :: a, b
     INTEGER t
   !----------------------------------------------------------------------------
   !BEGIN
     t = a
     a = b
     b = t
   END SUBROUTINE swap

   !############################################################################
   INTEGER FUNCTION cmp(l,r)
   !----------------------------------------------------------------------------
      INTEGER,INTENT(in)::l,r
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( RA(l) .LT. RA(r) ) THEN
         cmp = -1
      ELSEIF ( RA(l) .EQ. RA(r) ) THEN
         cmp = 0
      ELSE
         cmp = 1
      ENDIF
   END FUNCTION cmp
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE qsort
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE PO4AdsorptionFraction(PO4AdsorptionModel, &
                                 PO4tot_,            &
                                 ParticleConc_,      &
                                 Kpo4p,K,Qm,         &
                                 PO4dis,PO4par,      &
                                 thepH)
!-------------------------------------------------------------------------------
! Routine to compute fraction of PO4 adsorped to sediment/particulate concentration
!-------------------------------------------------------------------------------
    INTEGER,  INTENT(IN)  :: PO4AdsorptionModel
    AED_REAL, INTENT(IN)  :: PO4tot_, ParticleConc_
    AED_REAL, INTENT(IN)  :: Kpo4p,K,Qm
    AED_REAL, INTENT(OUT) :: PO4dis, PO4par
    AED_REAL, INTENT(IN), OPTIONAL :: thepH

!-------------------------------------------------------------------------------
!LOCALS
    AED_REAL :: buffer, f_pH, pH
    AED_REAL :: PO4tot, ParticleConc
    AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10

!
!-------------------------------------------------------------------------------
!BEGIN
   PO4dis   = zero_
   PO4par   = zero_
   buffer   = zero_
   f_pH     = one_

   ! calculate the total possible PO4 for sorption, and solids
   PO4tot        = MAX(one_e_neg_ten, PO4tot_ )       ! Co in Chao (mg)
   ParticleConc  = MAX(one_e_neg_ten, ParticleConc_ ) ! s in Chao  (mg = mol/L * g/mol * mg/g)


   IF(PO4AdsorptionModel == 1) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Ji 2008:
     !
     ! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.

     PO4par = (Kpo4p*ParticleConc) / (one_+Kpo4p*ParticleConc) * PO4tot
     PO4dis = one_ / (one_+Kpo4p*ParticleConc) * PO4tot


   ELSEIF(PO4AdsorptionModel == 2) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Chao et al. 2010:
     !
     ! Chao, X. et al. 2010. Three-dimensional numerical simulation of
     !   water quality and sediment associated processes with application
     !   to a Mississippi delta lake. J. Environ. Manage. 91 p1456-1466.
     IF(PRESENT(thepH)) THEN
       pH = thepH
       IF(thepH > 11.) pH = 11.0
       IF(thepH < 3.)  pH = 3.0

       ! -0.0094x2 + 0.0428x + 0.9574
       ! (ursula.salmon@uwa.edu.au: fPH for PO4 sorption to Fe in Mine Lakes)
       f_pH = -0.0094*pH*pH + 0.0428*pH + 0.9574
     ELSE
       f_pH = one_
     END IF

     ! calculate particulate fraction based on quadratic solution

     ! Chao Eq 16
     buffer = SQRT(((PO4tot+(1./K)-(ParticleConc*Qm*f_pH)))**2. + (4.*f_pH*ParticleConc*Qm/K))
     PO4par  = 0.5 * ((PO4tot+(1./K)+(ParticleConc*Qm*f_pH))  - buffer  )

     ! Check for stupid solutions
     IF(PO4par > PO4tot) PO4par = PO4tot
     IF(PO4par < zero_) PO4par = zero_

     ! Now set dissolved portion
     PO4dis = PO4tot - PO4par

   ELSE
     !-----------------------------------------------------
     ! No model is selected

     PO4dis = PO4tot
     PO4par = zero_

   END IF

END SUBROUTINE PO4AdsorptionFraction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION in_zone_set(matz, active_zones)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: matz
   AED_REAL,INTENT(in) :: active_zones(:)
!
!LOCALS
   INTEGER :: i, l
   LOGICAL :: res
!BEGIN
!-------------------------------------------------------------------------------
   res = .FALSE.
   l = size(active_zones)
   DO i=1,l
      IF ( active_zones(i) == matz ) THEN
         res = .TRUE.
         EXIT
      ENDIF
   ENDDO

   in_zone_set = res
END FUNCTION in_zone_set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION water_viscosity(temperature) RESULT(mu)
!-------------------------------------------------------------------------------
! Calculates the molecular viscosity of water for a given temperature
!
! From Table A.1b, FLUID_MECHANICS With Engineering Applications
! by Robert L. Daugherty and Joseph B. Franzini,
! however, note these values are common in most fluid mechanics texts.
! NOTE: N s / m^2  = kg / m / s
!
!  Temp (C)     Viscosity (N s / m^2) x 10^3
!  --------     ---------
!      0          1.781
!      5          1.518
!     10          1.307
!     15          1.139
!     20          1.002
!     25          0.890
!     30          0.798
!     40          0.653
!     50          0.547
!     60          0.466
!     70          0.404
!     80          0.354
!     90          0.315
!    100          0.282
!
!-------------------------------------------------------------------------------
!ARGUMENTS
  AED_REAL,INTENT(inout)  :: temperature
  AED_REAL :: mu
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   !-- Check for non-sensical temperatures
   IF( temperature<zero_ ) temperature = 0.0
   IF( temperature>100.0 ) temperature = 100.0

   IF( temperature<=20.0 ) THEN
     ! 0C to 20C
     ! y = 0.0008 * x^2 - 0.0556 * x + 1.7789
     ! r^2 = 0.9999
     mu = 0.0008 * temperature**2. - 0.0556 * temperature + 1.7789

   ELSEIF(temperature <= 60) THEN
     ! 20C to 60C
     ! y = 0.0002 * x^2 - 0.0323 * x + 1.5471
     ! r^2 = 0.9997
     mu = 0.0002 * temperature**2. - 0.0323 * temperature + 1.5471
   ELSE
     ! 60C to 100C
     ! y = 0.00006 * x^2 - 0.0141 * x + 1.1026
     ! r^2 = 0.9995
     mu = 0.00006 * temperature**2. - 0.0141 * temperature + 1.1026
   ENDIF

   ! Now convert to units of: N s / m^2
   mu = mu / 1e3

END FUNCTION water_viscosity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed2_util
