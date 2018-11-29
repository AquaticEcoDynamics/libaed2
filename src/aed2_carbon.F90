!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ______   | || |      __      | || |  _______     | |         !
!         | |   .' ___  |  | || |     /  \     | || | |_   __ \    | |         !
!         | |  / .'   \_|  | || |    / /\ \    | || |   | |__) |   | |         !
!         | |  | |         | || |   / ____ \   | || |   |  __ /    | |         !
!         | |  \ `.___.'\  | || | _/ /    \ \_ | || |  _| |  \ \_  | |         !
!         | |   `._____.'  | || ||____|  |____|| || | |____| |___| | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################
!#                                                                             #
!# aed2_carbon.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2018 -  The University of Western Australia               #
!#                                                                             #
!#   GLM is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   GLM is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_carbon
!-------------------------------------------------------------------------------
! aed2_carbon --- carbon biogeochemical model
!
! The AED module carbon contains equations that describe exchange of
! dissolved inorganic carbon across the air/water interface and sediment flux,
! and simulation of methane.
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_util,  ONLY: aed2_gas_piston_velocity

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_carbon_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_carbon_data_t
      !# Variable identifiers
      INTEGER  :: id_dic, id_pH, id_ch4, id_oxy, id_talk
      INTEGER  :: id_Fsed_dic, id_Fsed_ch4
      INTEGER  :: id_temp, id_salt
      INTEGER  :: id_wind, id_vel, id_depth
      INTEGER  :: id_ch4ox, id_pco2
      INTEGER  :: id_sed_dic
      INTEGER  :: id_atm_co2, id_atm_ch4
      INTEGER  :: id_par, id_extc, id_dz

      !# Model parameters
      AED_REAL :: Fsed_dic, Ksed_dic, theta_sed_dic
      AED_REAL :: Fsed_ch4, Ksed_ch4, theta_sed_ch4
      AED_REAL :: Rch4ox, Kch4ox, vTch4ox, atm_co2, atm_ch4, ionic
      AED_REAL :: maxMPBProdn, IkMPB

      LOGICAL  :: use_oxy, use_sed_model_dic, use_sed_model_ch4
      LOGICAL  :: simDIC, simCH4

      INTEGER  :: alk_mode, co2_model, co2_piston_model, ch4_piston_model

     CONTAINS
         PROCEDURE :: define            => aed2_define_carbon
         PROCEDURE :: calculate_surface => aed2_calculate_surface_carbon
         PROCEDURE :: calculate         => aed2_calculate_carbon
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_carbon
         PROCEDURE :: equilibrate       => aed2_equilibrate_carbon
!        PROCEDURE :: mobility          => aed2_mobility_carbon
!        PROCEDURE :: light_extinction  => aed2_light_extinction_carbon
!        PROCEDURE :: delete            => aed2_delete_carbon
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_carbon(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS

   INTEGER  :: status

   INTEGER           :: co2_model        = 1
   INTEGER           :: alk_mode         = 1
   INTEGER           :: co2_piston_model = 1
   INTEGER           :: ch4_piston_model = 1
   AED_REAL          :: pH_initial       = 7.5
   AED_REAL          :: ionic            = 0.0
   AED_REAL          :: dic_initial      = 1000.0
   AED_REAL          :: Fsed_dic         = 0.0
   AED_REAL          :: Ksed_dic         = 30.0
   AED_REAL          :: theta_sed_dic    = 1.0
   CHARACTER(len=64) :: Fsed_dic_variable=''
   AED_REAL          :: ch4_initial      = 4.5
   AED_REAL          :: Fsed_ch4         = 0.0
   AED_REAL          :: Ksed_ch4         = 30.0
   AED_REAL          :: theta_sed_ch4    = 1.0
   CHARACTER(len=64) :: Fsed_ch4_variable=''
   AED_REAL          :: Rch4ox           = 0.01
   AED_REAL          :: Kch4ox           = 0.01
   AED_REAL          :: vTch4ox          = 1.05
   AED_REAL          :: atm_co2          = 367e-6
   AED_REAL          :: atm_ch4          = 1.76e-6
   CHARACTER(len=64) :: methane_reactant_variable=''

   AED_REAL :: maxMPBProdn =  40.0   ! mmolC/m2/day
   AED_REAL :: IkMPB       = 180.0   ! Light sensitivity of MPB

   NAMELIST /aed2_carbon/ dic_initial,pH_initial,ch4_initial,ionic,         &
                         Fsed_dic,Ksed_dic,theta_sed_dic,Fsed_dic_variable, &
                         Fsed_ch4,Ksed_ch4,theta_sed_ch4,Fsed_ch4_variable, &
                         atm_co2,atm_ch4,Rch4ox,Kch4ox,vTch4ox,             &
                         methane_reactant_variable, &
                         maxMPBProdn, IkMPB, &
                         co2_model, alk_mode, co2_piston_model, ch4_piston_model

!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_carbon initialization"

   !# Read the namelist
   read(namlst,nml=aed2_carbon,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist aed2_carbon'
      STOP
   ENDIF

   !# Store parameter values in our own derived type
   !  NB: all rates must be provided in values per day,
   !  and are converted here to values per second.
   data%Fsed_dic      = Fsed_dic/secs_per_day
   data%Ksed_dic      = Ksed_dic
   data%theta_sed_dic = theta_sed_dic
   data%ionic         = ionic
   data%Fsed_ch4      = Fsed_ch4/secs_per_day
   data%Ksed_ch4      = Ksed_ch4
   data%theta_sed_ch4 = theta_sed_ch4
   data%Rch4ox        = Rch4ox/secs_per_day
   data%Kch4ox        = Kch4ox
   data%vTch4ox       = vTch4ox
   data%atm_co2        = atm_co2
   data%atm_ch4        = atm_ch4
   data%simDIC        = .false.
   data%simCH4        = .false.
   data%maxMPBProdn   = maxMPBProdn
   data%IkMPB         = IkMPB
   data%co2_model     = co2_model
   data%alk_mode      = alk_mode
   data%co2_piston_model = co2_piston_model
   data%ch4_piston_model = ch4_piston_model

   !# Register state variables
   IF (dic_initial>MISVAL) THEN
      data%id_dic = aed2_define_variable('dic','mmol/m**3','dissolved inorganic carbon',     &
                                       dic_initial,minimum=zero_)
      data%simDIC = .true.
      data%id_pH = aed2_define_variable('pH','-','pH',     &
                                       pH_initial,minimum=zero_)
   ENDIF

   IF (ch4_initial>MISVAL) THEN
      data%id_ch4 = aed2_define_variable('ch4','mmol/m**3','methane',    &
                                     ch4_initial,minimum=zero_)
      data%simCH4 = .true.
   ENDIF

   !# Register external state variable dependencies
   data%use_oxy = methane_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) THEN
      data%id_oxy = aed2_locate_variable(methane_reactant_variable)
   ENDIF

   data%use_sed_model_dic = Fsed_dic_variable .NE. ''
   IF (data%use_sed_model_dic) &
      data%id_Fsed_dic = aed2_locate_global_sheet(Fsed_dic_variable)

   data%use_sed_model_ch4 = Fsed_ch4_variable .NE. ''
   IF (data%use_sed_model_ch4) &
      data%id_Fsed_ch4 = aed2_locate_global_sheet(Fsed_ch4_variable)

   !# Register diagnostic variables
   data%id_pco2 = aed2_define_diag_variable('pCO2','atm', 'pCO2')

   data%id_ch4ox = aed2_define_diag_variable('ch4ox','mmol/m**3/d', 'methane oxidation rate')
   data%id_sed_dic = aed2_define_sheet_diag_variable('sed_dic','mmol/m**2/d',        &
                                                      'CO2 exchange across sed/water interface')

   data%id_atm_co2 = aed2_define_sheet_diag_variable('atm_co2_flux',            &
                             'mmol/m**2/d', 'CO2 exchange across atm/water interface')
   data%id_atm_ch4 = aed2_define_sheet_diag_variable('atm_ch4_flux',            &
                             'mmol/m**2/d', 'CH4 exchange across atm/water interface')

   !# Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')
   data%id_wind = aed2_locate_global_sheet('wind_speed')
   data%id_extc = aed2_locate_global('extc_coef')
   data%id_par  = aed2_locate_global('par')
   data%id_dz   = aed2_locate_global('layer_ht')
   data%id_vel  = aed2_locate_global('cell_vel') ! needed for k600
   data%id_depth= aed2_locate_global('layer_ht')

END SUBROUTINE aed2_define_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: dic,ch4,oxy,temp
   AED_REAL :: ch4oxidation
!
!-------------------------------------------------------------------------------
!BEGIN

   IF(data%simDIC .AND. data%simCH4) THEN
      ! Retrieve current (local) state variable values.
      dic = _STATE_VAR_(data%id_dic)! carbon
      ch4 = _STATE_VAR_(data%id_ch4)! carbon

      !# Retrieve current dependent state variable values.
      IF (data%use_oxy) THEN ! & use_oxy
         oxy = _STATE_VAR_(data%id_oxy)! oxygen
      ELSE
         oxy = 0.0
      ENDIF

      !# Retrieve current environmental conditions.
      temp = _STATE_VAR_(data%id_temp) ! temperature

      !# Define some intermediate quantities units mmol C/m3/day
      ch4oxidation = aed2_carbon_fch4ox(data%use_oxy,data%Rch4ox,data%Kch4ox,data%vTch4ox,oxy,temp)

      !# Set temporal derivatives
      _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (ch4*ch4oxidation)
      _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (-ch4*ch4oxidation)

      !# If an externally maintained oxygen pool is present, take nitrification from it
      IF (data%use_oxy) then ! & use_oxy
         _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-(32./12.)*ch4*ch4oxidation)
      ENDIF

      !# Export diagnostic variables
      _DIAG_VAR_(data%id_ch4ox) =  ch4*ch4oxidation*secs_per_day
   ENDIF

END SUBROUTINE aed2_calculate_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, S, T, depth, vel

   ! State
   AED_REAL :: dic,ph,ch4,talk = 0.,TCO2 = 0.

   ! Temporary variables
   AED_REAL :: pCO2 = 0.,FCO2,FCH4,henry
   AED_REAL :: Ko,kCH4,KCO2, CH4solub
   AED_REAL :: Tabs,windHt,atm
   AED_REAL :: A1,A2,A3,A4,B1,B2,B3,logC
   AED_REAL :: a,b,c,dcf
   AED_REAL :: ca, bc, cb, carba, bicarb, carb, om_cal, om_arg
   AED_REAL :: p00,p10,p01,p20,p11,p02

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simDIC .AND. .NOT.data%simCH4) RETURN


   !----------------------------------------------------------------------------
   !# Get dependent state variables from physical driver
   windHt = 10.
   wind   = _STATE_VAR_S_(data%id_wind) ! Wind speed at 10 m above surface (m/s)
   temp   = _STATE_VAR_(data%id_temp)   ! Temperature (degrees Celsius)
   salt   = _STATE_VAR_(data%id_salt)   ! Salinity (psu)
   depth  = MAX( _STATE_VAR_(data%id_depth), one_ )
   IF (data%id_vel > 0 ) THEN
     vel = _STATE_VAR_(data%id_vel)
   ELSE
    ! vel  = SQRT(_STATE_VAR_(data%id_E_tau)/_STATE_VAR_(data%id_E_dens))
    ! vel = vel/0.41 * log(depth/0.01)
     vel = 0.0001
   ENDIF

   Tabs = temp + 273.15

   !# Solubility, Ko (mol/L/atm)
   Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
                                 + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
   Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
   Ko = exp(Ko)


   !----------------------------------------------------------------------------
   !# CO2 concentration in the surface layer, depends on DIC, TA
   IF(data%simDIC) THEN

      !# Retrieve current (local) state variable values.
     dic = _STATE_VAR_(data%id_dic)! Concentration of carbon in surface layer
     ph  = _STATE_VAR_(data%id_pH)  ! Concentration of carbon in surface layer


     IF( data%co2_model == 1 ) THEN
       !# Use the Haltafall CO2 code for computing pCO2 & pH

       S=salt; T=temp

       IF( data%alk_mode == 1 ) THEN
         ! talk = 520.1 + 51.24*S  ! Atlantic (Millero 1998) from fabm, not suitable for estuaries
         ! talk = 1136.1 + 1.2*S*S + 2.8*S !Chesapeake Bay (George et al., 2013)
         talk =  1627.4 + 22.176*S   !regression from Naomi's data on Caboolture
         a    =  8.24493d-1 - 4.0899d-3*T + 7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
         b    = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2
         c    =  4.8314d-4
         dcf  = (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
                - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3

         ! next, use the CO2 module (same to CDIAC module) to calc pCO2
         talk = talk / 1.0D6      ! change unit to mol/kgSW
         TCO2 = dic / (1.0D6*dcf) ! change unit to mol/kgSW

       ELSEIF( data%alk_mode == 2 ) THEN
         p00  =       1063
         p10  =      1.751
         p01  =   -0.05369
         p20  =     0.2266
         p11  =  -0.001252
         p02  =  0.0002546
         a    =  8.24493d-1 - 4.0899d-3*T + 7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
         b    = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2
         c    =  4.8314d-4
         dcf  = (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
                      - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3

         talk = p00 + p10*S + p01*dic + p20*S**2 + p11*dic*S + p02*dic**2
         TCO2 = dic / (1.0D6*dcf) ! change unit to mol/kgSW

       ENDIF

       CALL CO2DYN ( TCO2, talk, T, S, pCO2, pH, HENRY, ca, bc, cb)

         ! Adjust outputs back to units used in the parent model code (e.g. mmol/m3) if appropriate
         ! note the output pCO2 is in unit of ATM
!       pCO2 = pCO2*1.0D6   ! partial pressure of co2 in water
!       _STATE_VAR_(data%id_talk) = talk*(1.0D6)           ! total alkalinity (umol/kg)

     ELSEIF ( data%co2_model == 2 ) THEN
       !# Use the Butler CO2 code for computing pCO2 & pH
       pCO2 = aed2_carbon_co2(data%ionic,temp,dic,ph)*1e-6 / Ko  !(=atm), use Yanti's script for pCO2

     ELSEIF ( data%co2_model == 0 ) THEN
       !# Use the aed2_geochem module for computing pCO2 & pH
       pCO2 = _DIAG_VAR_(data%id_pco2)  ! this diagnostic is getting set in aed2_geocehmistry

     ENDIF

     _DIAG_VAR_(data%id_pco2) = pCO2

     !# Now compute piston velocity, k
     kCO2 = aed2_gas_piston_velocity(windHt,wind,temp,salt,  &
         vel=vel,depth=depth,schmidt_model=2,piston_model=data%co2_piston_model)


     !# Now compute the CO2 flux
     ! FCO2 = kCO2 * Ko * (pCO2 - PCO2a)
     ! pCO2a = 367e-6 atm (Keeling & Wharf, 1999)
     ! mmol/m2/s = m/s * mmol/m3/atm * atm
     FCO2 = kCO2 * (1e6*Ko) * (pCO2 - data%atm_co2)
     ! FCO2 = - kCO2 * Ko*1e6 * ((pCO2 * 1e-6) - data%atm_co2) ! dCO2/dt

     !--------------------------------------------------------------------------
     !# Transfer surface exchange value to AED2 (mmmol/m2/s) converted by driver
     _FLUX_VAR_T_(data%id_dic) = -FCO2

     !# Also store co2 flux across the atm/water interface as a
     !  diagnostic variable (mmmol/m2/d)
     _DIAG_VAR_S_(data%id_atm_co2) = FCO2*secs_per_day

   END IF


   !----------------------------------------------------------------------------
   !# CH4 flux
   IF(data%simCH4) THEN
     ! Algorithm from Arianto Santoso <abs11@students.waikato.ac.nz>

     ! Concentration of methane in surface layer
     ch4 = _STATE_VAR_(data%id_ch4)

     ! Piston velocity for CH4
     kCH4 = aed2_gas_piston_velocity(windHt,wind,temp,salt, &
         vel=vel,depth=depth,schmidt_model=4,piston_model=data%ch4_piston_model)

     ! Solubility, Ko (mol/L/atm)
     atm = data%atm_ch4    ! 1.76 * 1e-6 !## current atmospheric CH4 data from NOAA (in ppm? atm)

     A1 = -415.2807
     A2 =  596.8104
     A3 =  379.2599
     A4 =  -62.0757
     B1 =   -0.05916
     B2 =    0.032174
     B3 =   -0.0048198

     logC = (log(atm)) + A1                                                   &
          + (A2 * (100./Tabs)) + (A3 * log (Tabs/100.)) + (A4 * (Tabs/100.))  &
          + salt * (B1 + (B2  * (Tabs/100.)) + (B3 * (Tabs/100.)*(Tabs/100.)))

     CH4solub = exp(logC) * 1e-3

     !# Now compute methane atm flux  (mmol/m2/s = m/s * mmol/m3)
     FCH4 = kCH4 *  (ch4 - CH4solub)

     !----------------------------------------------------------------------------
     !# Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
     _FLUX_VAR_T_(data%id_ch4) = -FCH4

     !# Also store CH4 flux across the atm/water interface as
     !  diagnostic variable (mmmol/m2/d)
     _DIAG_VAR_S_(data%id_atm_ch4) = FCH4*secs_per_day

   END IF
   !----------------------------------------------------------------------------

END SUBROUTINE aed2_calculate_surface_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED carbon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, par, extc, dz

   ! State
   AED_REAL :: dic, oxy, mpb, ph

   ! Temporary variables
   AED_REAL :: dic_flux, ch4_flux, Fsed_dic, Fsed_ch4
   !AED_REAL, PARAMETER :: maxMPBProdn = 40.     ! mmolC/m2/day                     !
   !AED_REAL, PARAMETER :: IkMPB       = 180.0   ! Light sensitivity of MPB  !

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simDIC) RETURN


   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   par  = _STATE_VAR_(data%id_par)  ! local par
   dz   = _STATE_VAR_(data%id_dz)   ! local layer depth
   extc = _STATE_VAR_(data%id_extc) ! local extinction

    ! Retrieve current (local) state variable values.
   dic  = _STATE_VAR_(data%id_dic)! carbon
   pH   = _STATE_VAR_(data%id_pH)! pH

   IF ( data%use_sed_model_dic ) THEN
      Fsed_dic = _STATE_VAR_S_(data%id_Fsed_dic)
   ELSE
       Fsed_dic = data%Fsed_dic
   ENDIF
   IF ( data%use_sed_model_ch4 ) THEN
      Fsed_ch4 = _STATE_VAR_S_(data%id_Fsed_ch4)
   ELSE
       Fsed_ch4 = data%Fsed_ch4
   ENDIF

   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      dic_flux = Fsed_dic * oxy/(data%Ksed_dic+oxy) * (data%theta_sed_dic**(temp-20.0))
      ch4_flux = Fsed_ch4 * data%Ksed_ch4/(data%Ksed_ch4+oxy) * (data%theta_sed_ch4**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      dic_flux = Fsed_dic * (data%theta_sed_dic**(temp-20.0))
      ch4_flux = Fsed_ch4 * (data%theta_sed_ch4**(temp-20.0))
   ENDIF


  !! Allow photosynthetic production of CO2 in the benthos due to MPB if light and suitable pH
  !par = par * (exp(-extc*dz))
  ! IF( par > 50. .AND. pH > 5.5 .AND. pH < 9.6 ) THEN
  !   mpb = (data%maxMPBProdn/secs_per_day)*(1.0-exp(-par/data%IkMPB)) * (data%theta_sed_dic**(temp-20.0))
  !   dic_flux = Fsed_dic - mpb
  !   IF (data%use_oxy) THEN
  !     _FLUX_VAR_(data%id_oxy) =  _FLUX_VAR_(data%id_oxy) + mpb
  !   ENDIF
  ! ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   !_SET_BOTTOM_FLUX_(data%id_dic,dic_flux/secs_per_day)
   !_SET_SED_FLUX_(data%id_dic,dic_flux)
   _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (dic_flux)
   _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (ch4_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_dic) = _FLUX_VAR_B_(data%id_ben_dic) + (-dic_flux/secs_per_day)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_dic) = dic_flux

END SUBROUTINE aed2_calculate_benthic_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_equilibrate_carbon(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update pH after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_carbon_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! State
   AED_REAL :: dic, pH, pCO2, temp, salt
   AED_REAL :: S,T,a,b,c,dcf,talk = 0.,TCO2 = 0.,ca,bc,cb,HENRY

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT.data%simDIC) RETURN

    pCO2 = zero_
    pH   = zero_

    !# Retrieve current (local) state variable values.
    dic  = _STATE_VAR_(data%id_dic) ! Concentration of DIC in the cell
    salt = _STATE_VAR_(data%id_salt) ! Concentration of DIC in the cell
    temp = _STATE_VAR_(data%id_temp) ! Concentration of DIC in the cell


    IF( data%co2_model == 1 ) THEN
      !# Use the Haltafall CO2 code for computing pCO2 & pH

      S=salt; T=temp

      IF( data%alk_mode == 1 ) THEN

        S=salt; T=temp

        ! talk = 520.1 + 51.24*S  ! Atlantic (Millero 1998) from fabm, not suitable for estuaries
        ! talk = 1136.1 + 1.2*S*S + 2.8*S !Chesapeake Bay (George et al., 2013)
        talk =  1627.4 + 22.176*S   !regression from Naomi's data on Caboolture
        a    =  8.24493d-1 - 4.0899d-3*T + 7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
        b    = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2
        c    =  4.8314d-4
        dcf  = (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
                 - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3

        ! next, use the CO2 module (same to CDIAC module) to calc pCO2
        talk  = talk / 1.0D6      ! change unit to mol/kgSW
        TCO2  = dic / (1.0D6*dcf) ! change unit to mol/kgSW
      ENDIF

      CALL CO2DYN ( TCO2, talk, T, S, pCO2, pH, HENRY, ca, bc, cb)

    ENDIF

    !# SET PCO2 & pH as returned
    _DIAG_VAR_(data%id_pco2) = pCO2
    _STATE_VAR_(data%id_pH)  =  pH

END SUBROUTINE aed2_equilibrate_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_carbon_fch4ox(use_oxy,Rch4ox,Kch4ox,vTch4ox,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for methane oxidation
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rch4ox,Kch4ox,vTch4ox,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      aed2_carbon_fch4ox = Rch4ox * oxy/(Kch4ox+oxy) * (vTch4ox**(temp-20.0))
   ELSE
      aed2_carbon_fch4ox = Rch4ox * (vTch4ox**(temp-20.0))
   ENDIF

END FUNCTION aed2_carbon_fch4ox
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed2_carbon_co2(ionic, temp, dic, pH)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T, from Butler
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL, INTENT(IN) :: ionic, dic, temp, pH
!
!LOCALS
   ! Temporary variables
   AED_REAL :: K_h, Kw, Ka1, Ka2, i_f
   AED_REAL :: H, CO2, HCO3, CO3, TA
!-------------------------------------------------------------------------------
!BEGIN

   ! Acidity constants temperature dependence

   ! pKh  =  -0.000075324675x2 + 0.016279653680x + 1.110424242424
   ! pKa1 = 0.000142121212x2 - 0.012648181818x + 6.577539393939
   ! pKa2 =  0.000113679654x2 - 0.014687186147x + 10.625769696970
   ! pKw  =   0.000201991342x2 - 0.043419653680x + 14.949709090909

   K_h = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
   Ka1 =  0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
   Ka2 =  0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
   Kw  =  0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909


   ! Ionic strength dependence

   ! 1st calculate function f
   i_f = (((SQRT(ionic)) / (1+SQRT(ionic))) -0.20*ionic) * &
                       (298.0/(temp+273.))**0.666667

   ! pKh = pKh(0) + bI
   ! b = 0.105 (Butler, 1982)
   K_h = K_h + 0.105*ionic

   ! pKw = pKw(0) - f
   Kw = Kw - i_f

   ! pKa1 = pKa1(0) - f - bI
   Ka1 = Ka1 - i_f - 0.105*ionic

   !pKa2 = pKa2(0) - 2f
   Ka2 = Ka2 + 2.0*i_f

   ! Convert from pK etc to Kh, Kw, Ka1, Ka2
   K_h  = 10.**(-K_h)
   Ka1 = 10.**(-Ka1)
   Ka2 = 10.**(-Ka2)
   Kw  = 10.**(-Kw)


   ! Calculate the speciation to know the molar mass of DIC                                                             !
   H    = 10.**(-pH)
   CO3  = (Ka1*Ka2)/(H*H + Ka1*H + Ka1*Ka2)
   HCO3 = (Ka1*H)/(H*H + Ka1*H + Ka1*Ka2)
   CO2  = (H*H)/(H*H + Ka1*H + Ka1*Ka2)


   ! and update speciation (mol C/L)
   CO3  = dic*CO3
   HCO3 = dic*HCO3
   CO2  = dic*CO2

   ! calculate TA for the previous timestep
   TA = dic * (Ka1*H + 2.0*Ka1*Ka2) / (H*H + Ka1*H + Ka1*Ka2)
   TA = TA + (Kw/H) - H

   aed2_carbon_co2 = CO2

END FUNCTION aed2_carbon_co2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE CO2DYN ( TCO2, TA, T, S, pCO2, pH, HENRY, ca, bc, cb)
!-------------------------------------------------------------------------------
!# This subroutine acts as an interface to the Haltafall iteration, by setting
!  options etc. In this case TCO2 and TA are set and pH and pCO2 are returned
IMPLICIT NONE
AED_REAL :: PRSS, pH, AKVAL, CONCS,                                            &
                  TCO2, TA, T, S, pCO2,                                        &
                  SOLBTY, CCO2,                                                &
                  A1, A2, A3, B1, B2, B3, TK, TK1, SOL1, SOL2,                 &
                  HENRY, ca, bc, cb
INTEGER :: MCONC, MKVAL, ICONST, ICALC

! INPUT PARAMETERS:
PARAMETER ( MCONC = 9, MKVAL = 4 )

DIMENSION AKVAL(MKVAL), CONCS(MCONC)
!-------------------------------------------------------------------------------

  ICONST   = 6
  PRSS     = 1.0d0
  CONCS(1) = TCO2
  CONCS(2) = TA
  ICALC    = 1

  CALL POLYCO(PRSS,T,S,CONCS,MCONC,AKVAL,MKVAL,ICALC,ICONST)

  pCO2     = CONCS(3)
  pH       = CONCS(4)
  ca       = CONCS(5)
  bc       = CONCS(6)
  cb       = CONCS(7)
  HENRY    = AKVAL(1)

  RETURN
END SUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE POLYCO(PD,TD,SD,CONCS,NCONC,AKVAL,NKVAL,ICALC,ICONST)
!     -----------------------------------------------------------------
! MASTER SUBROUTINE FOR CALCULATION OF THE CO2 SYSTEM THERMODYNAMICS
!

! EXPLANATION OF POLYCO PARAMETERS
!       P - PRESSURE IN ATMOSPHERES (P<>1 NOT YET CODED)
!       T - TEMPERATURE IN DEG.C
!       S - SALINITY IN PPT
!      CONCS(1) - TOTAL C (MOL/KG)
!      CONCS(2) - TOTAL ALKALINITY (MOL/KG)
!      CONCS(3) - PCO2 (ATM)
!      CONCS(4) - PH
!      CONCS(5) - {H2CO3} (MOL/KG)
!      CONCS(6) - {HCO3} (MOL/KG)
!      CONCS(7) - {CO3} (MOL/KG)
!        CONCS(8) - CARBONATE ALKALINITY  ) FOR ICONST = 4,5,6
!        CONCS(9) - BORATE ALKALINITY     )       ONLY
!         NCONC - SIZE OF CONCS ARRAY (7 FOR ICONST=1,2,3; 9 FOR ICONST
!     AKVAL(1) - KP (HENRY'S LAW CONSTANT) (MOL/KG/ATM)
!     AKVAL(2) - K1C (H2CO3 DISSOCIATION) (MOL/KG)
!       AKVAL(3) - K2C (HCO3 DISSOCIATION) (MOL/KG)
!       AKVAL(4) - KB (B(OH)3 DISSOCIATION) (MOL/KG)  FOR ICONST=4,5,6
!        NKVAL - SIZE OF AKVAL ARRAY (3 FOR ICONST=1,2,3; 4 FOR ICONST=
!        ICALC - SELECTION OF THE TWO INPUT PARAMETERS:
!       ICALC = 1  TOTAL C AND ALKALINITY
!       ICALC = 2  TOTAL C AND PCO2
!       ICALC = 3  TOTAL C AND PH
!       ICALC = 4  ALKALINITY AND PCO2
!       ICALC = 5  ALKALINITY AND PH
!       ICALC = 6  PCO2 AND PH
!       ICALC = 7  CALCULATE CONSTANTS AKVAL ONLY
!       ICONST - SELECTION OF PH SCALE AND COMPONENTS:
!       ICONST = 1  NBS PH SCALE
!       ICONST = 2  HANSSON'S SCALE (SWS WITHOUT FLUORIDE)
!       ICONST = 3  SWS PH SCALE
!       ICONST = 4  AS 1 BUT INCLUDING BORATE IN THE CALCULATION
!       ICONST = 5  AS 2 BUT INCLUDING BORATE IN THE CALCULATION
!       ICONST = 6  AS 3 BUT INCLUDING BORATE IN THE CALCULATION

!  NOTE: FOR ICONST=1,2,3 CONCS(2) REPRESENTS CARBONATE ALKALINITY SINC
!        BORATE IS NOT INCLUDED IN THE CALCULATION. FOR ICONST=4,5,6 CO
!        REPRESENTS TOTAL ALKALINITY (CARBONATE + BORATE), THE COMPONEN
!        WHICH ARE GIVEN IN CONCS(8) AND CONCS(9)
!     -----------------------------------------------------------------


IMPLICIT NONE

AED_REAL :: PMIN, PMAX, SMIN, SMAX, TMIN, TMAX, CONCS,  &
&      AKVAL, PD, TD, SD, P, T, S, BTOT
INTEGER MINJC, MAXJC, MINJK, MAXJK, MINCAL, MAXCAL, MINCON,  &
&      MAXCON, NCONC, NKVAL, ICALC, ICONST, IC
LOGICAL BORON

PARAMETER(MINJC=7,MAXJC=9,MINJK=3,MAXJK=4)
PARAMETER(MINCAL=1,MAXCAL=7,MINCON=1,MAXCON=6)
PARAMETER(PMIN=0.99999D0,PMAX=1.00001D0,SMIN=0.0D0,  &
&  SMAX=45.0D0,TMIN=-4.0D0,TMAX=40.0D0)
DIMENSION CONCS(NCONC),AKVAL(NKVAL)
!     -----------------------------------------------------------------

P = PD
S = SD
T = TD

!     IF(T.LT.TMIN.OR.T.GT.TMAX)WRITE (*,*) P, S, T, TMIN, TMAX
IF(P.LT.PMIN.OR.P.GT.PMAX) write (*,*) 'fatal error'
IF(S.LT.SMIN.OR.S.GT.SMAX) write (*,*) 'fatal error'
IF(T.LT.TMIN.OR.T.GT.TMAX) write (*,*) 'fatal error'
IF(ICALC.LT.MINCAL.OR.ICALC.GT.MAXCAL)  &
&  write (*,*) 'fatal error'
IF(ICONST.LT.MINCON.OR.ICONST.GT.MAXCON)  &
&  write (*,*) 'fatal error'
BORON=(ICONST.GT.3)
IF(BORON) THEN
  IC=ICONST-3
  BTOT=0.0004128D0*S/35.0D0
  IF(NCONC.NE.MAXJC) write (*,*) 'fatal error'
  IF(NKVAL.NE.MAXJK) write (*,*) 'fatal error'
ELSE
  IC=ICONST
  IF(NCONC.NE.MINJC) write (*,*) 'fatal error'
  IF(NKVAL.NE.MINJK) write (*,*) 'fatal error'
ENDIF

CALL CO2SET(P,T,S,AKVAL,NKVAL,IC)
IF(ICALC.LT.MAXCAL)  &
& CALL CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT)

!if (concs(4).eq.100.) then
!   write (*,*) 'S,T,P',S,T,P
!   write (*,*) 'CONCS',CONCS
!   write (*,*) 'fatal error'
!end if

RETURN
END SUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE CO2SET(P,T,S,AKVAL,NKVAL,IC)
!     -----------------

! Routine to calculate CO2 system constants under the conditions set by
! P,S,T     (NOTE: PRESSURE <> 1ATM IS NOT YET CODED)

! I. Calculate constants at P=1 and S=0 using

!      ln K0  =  A + B/TK + C ln TK
!                                   (where TK is in Kelvin)

! II. Calculate constants at P=1 and salinity S using

!    ln K  =  ln K0 + (a0 + a1/TK + a2 ln TK) S**1/2
!                 + (b0 + b1TK + b2TK**2) S

! The sources of the coefficients are as follows:

!  IC=                  1                    2                  3
!               (NBS pH scale)        (SWS pH scale       (SWS pH scale
!                                       with no F)

!  KP            WEISS (1974)           WEISS(1974)          WEISS(1974

!  K1C )      MEHRBACH ACC. TO       HANSSON ACC. TO    HANSSON AND MEH
!  K2C )       MILLERO (1979)         MILLERO (1979)      ACC. TO DICKS
!   KB )                                                 AND MILLERO (1
!                                                         (K1C AND K2C
!                                                           HANSSON ACC
!                                                            MILLERO (1
!                                                              (KB ONLY

! ***
!      IMPLICIT real(rk) (A-H,O-Z)

!     Modified by jcb 17/02/10 to use OCMIP calculations of K1, K2, Kb.
!     Differences are subtle rather than significant
!use fabm_types

IMPLICIT NONE

INTEGER MAXK, MAXCON, NKVAL, ICON, IC, IK
! ***
PARAMETER(MAXK=4,MAXCON=3)
AED_REAL,DIMENSION(:) :: A(MAXK),B(MAXK),C(MAXK)
AED_REAL,DIMENSION(:) :: A0(MAXK,MAXCON),A1(MAXK,MAXCON),A2(MAXK,MAXCON)
AED_REAL,DIMENSION(:) :: B0(MAXK,MAXCON),B1(MAXK,MAXCON),B2(MAXK,MAXCON)
AED_REAL,DIMENSION(:) :: AKVAL(NKVAL)
AED_REAL              :: P,T,S,VAL,TK
AED_REAL              :: dlogTK, S2, sqrtS, S15, k1, k2, kb
DATA A/-167.8108D0, 290.9097D0, 207.6548D0, 148.0248D0/
DATA B/9345.17D0, -14554.21D0, -11843.79D0, -8966.9D0/
DATA C/23.3585D0, -45.0575D0, -33.6485D0, -24.4344D0/
DATA (A0(1,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (A0(2,ICON),ICON=1,MAXCON) /0.0221D0, 0.5709D0, -45.8076D0/
DATA (A0(3,ICON),ICON=1,MAXCON) /0.9805D0, 1.4853D0, -39.5492D0/
DATA (A0(4,ICON),ICON=1,MAXCON) /0.0473D0, 0.5998D0, 0.5998D0/
DATA (A1(1,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (A1(2,ICON),ICON=1,MAXCON) /34.02D0, -84.25D0, 1935.07D0/
DATA (A1(3,ICON),ICON=1,MAXCON) /-92.65D0, -192.69D0, 1590.14D0/
DATA (A1(4,ICON),ICON=1,MAXCON) /49.10D0, -75.25D0, -75.25D0/
DATA (A2(1,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (A2(2,ICON),ICON=1,MAXCON) /2*0.0D0,6.9513D0/
DATA (A2(3,ICON),ICON=1,MAXCON) /2*0.0D0,6.1523D0/
DATA (A2(4,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (B0(1,ICON),ICON=1,MAXCON) /3*0.023517D0/
DATA (B0(2,ICON),ICON=1,MAXCON) /0.0D0,-0.01632D0,-0.01566D0/
DATA (B0(3,ICON),ICON=1,MAXCON) /-0.03294D0,-0.05058D0,-0.04997D0/
DATA (B0(4,ICON),ICON=1,MAXCON) /0.0D0, -0.01767D0, -0.01767D0/
DATA (B1(1,ICON),ICON=1,MAXCON) /3*-2.3656D-4/
DATA (B1(2,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (B1(3,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (B1(4,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (B2(1,ICON),ICON=1,MAXCON) /3*4.7036D-7/
DATA (B2(2,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (B2(3,ICON),ICON=1,MAXCON) /3*0.0D0/
DATA (B2(4,ICON),ICON=1,MAXCON) /3*0.0D0/

TK=T+273.15D0
DO 100 IK=1,NKVAL
  VAL=A(IK) + B(IK)/TK + C(IK)*LOG(TK)
  VAL=VAL + (A0(IK,IC) + A1(IK,IC)/TK + A2(IK,IC)*LOG(TK))*SQRT(S)
  VAL=VAL + (B0(IK,IC) + B1(IK,IC)*TK + B2(IK,IC)*TK*TK)*S
  AKVAL(IK)=EXP(VAL)
100    CONTINUE

IF (IC .EQ. 3) THEN
!  Calculation of constants as used in the OCMIP process for ICONST = 3 or 6
!  see http://www.ipsl.jussieu.fr/OCMIP/
!  added jcb 17/02/10

!  Derive simple terms used more than once
dlogTK = log(TK)
S2 = S*S
sqrtS = sqrt(S)
S15 = S**1.5
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
k1=10**(-1*(3670.7/TK - 62.008 + 9.7944*dlogTK - &
&      0.0118 * S + 0.000116*S2))
k2=10**(-1*(1394.7/TK + 4.777 - &
&      0.0184*S + 0.000118*S2))
! kb = [H][BO2]/[HBO2]
! Millero p.669 (1995) using data from Dickson (1990)
kb=exp((-8966.90 - 2890.53*sqrtS - 77.942*S + &
&      1.728*S15 - 0.0996*S2)/TK + &
&      (148.0248 + 137.1942*sqrtS + 1.62142*S) + &
&      (-24.4344 - 25.085*sqrtS - 0.2474*S) * &
&      dlogTK + 0.053105*sqrtS*TK)
! replace haltafall calculations with OCMIP calculations
AKVAL(2) = k1
AKVAL(3) = k2
AKVAL(4) = kb
END IF ! section implimenting OCMIP coefficients

RETURN
END SUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT)
!     -----------------

! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975)
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE.


!      IMPLICIT real(rk) (A-H,O-Z)
!use fabm_types

IMPLICIT NONE

INTEGER NCONC, NKVAL, ICALC, II, KARL, LQ
AED_REAL              :: CTOT,ALK,PCO2,PH,H2CO3,HCO3,CO3,ALKC
AED_REAL              :: ALKB,AKP,AK1C,AK2C,AKB,BTOT
AED_REAL              :: AKR,AHPLUS
AED_REAL              :: PROD,tol1,tol2,tol3,tol4,steg,fak
AED_REAL              :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z
AED_REAL,DIMENSION(:) :: CONCS(NCONC),AKVAL(NKVAL),CONCS2(9),AKVAL2(4)
EQUIVALENCE (CTOT  , CONCS2(1)), (ALK   , CONCS2(2)),  &
&            (PCO2  , CONCS2(3)), (PH    , CONCS2(4)),  &
&            (H2CO3 , CONCS2(5)), (HCO3  , CONCS2(6)),  &
&            (CO3   , CONCS2(7)), (ALKC  , CONCS2(8)),  &
&            (ALKB  , CONCS2(9)),                       &
&            (AKP   , AKVAL2(1)), (AK1C  , AKVAL2(2)),  &
&            (AK2C  , AKVAL2(3)), (AKB   , AKVAL2(4))
LOGICAL BORON,DONE
integer :: it
integer,parameter :: maxiter = 1000

btot = 0.0 ; x1 = 0.0 ; x2 = 0.0 ; y1 = 0.0 ; y2 = 0.0  !## CAB [-Wmaybe-uninitialized]

DO 100 II=1,NCONC
  CONCS2(II)=CONCS(II)
100    CONTINUE
DO 110 II=1,NKVAL
  AKVAL2(II)=AKVAL(II)
110    CONTINUE
AKR = AK1C/AK2C
AHPLUS=10.0D0**(-PH)
PROD=AKR*AKP*PCO2

IF(BORON) THEN

  IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN
!         *** ALK, BTOT AND CTOT OR PCO2 FIXED ***
!         *** ITERATIVE CALCULATION NECESSARY HERE

!         SET INITIAL GUESSES AND TOLERANCE
    H2CO3=PCO2*AKP
    CO3=ALK/10.0D0
    AHPLUS=1.0D-8
    ALKB=BTOT
    TOL1=ALK/1.0D5
    TOL2=H2CO3/1.0D5
    TOL3=CTOT/1.0D5
    TOL4=BTOT/1.0D5

!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS
    KARL=1
    STEG=2.0D0
    FAK=1.0D0
    STEGBY=0.4D0
    do it=1,maxiter
       DONE=.TRUE.
       IF(ICALC.EQ.4) THEN
!         *** PCO2 IS FIXED ***
         Y=AHPLUS*AHPLUS*CO3/(AK1C*AK2C)
         IF(ABS(Y-H2CO3).GT.TOL2) THEN
           CO3=CO3*H2CO3/Y
           DONE=.FALSE.
         ENDIF
       ELSEIF(ICALC.EQ.1) THEN
!           *** CTOT IS FIXED ***
         Y=CO3*(1.0D0+AHPLUS/AK2C+AHPLUS*AHPLUS/(AK1C*AK2C))
         IF(ABS(Y-CTOT).GT.TOL3) THEN
           CO3=CO3*CTOT/Y
           DONE=.FALSE.
         ENDIF
       ENDIF
       Y=ALKB*(1.0D0+AHPLUS/AKB)
       IF(ABS(Y-BTOT).GT.TOL4) THEN
         ALKB=ALKB*BTOT/Y
         DONE=.FALSE.
       ENDIF

! Alkalinity is equivalent to -(total H+), so the sign of W is opposite
! to that normally used

       Y=CO3*(2.0D0+AHPLUS/AK2C)+ALKB
       IF(ABS(Y-ALK).GT.TOL1) THEN
         DONE=.FALSE.
         X=LOG(AHPLUS)
         W=SIGN(1.0D0,Y-ALK)
         IF(W.GE.0.0D0) THEN
           X1=X
           Y1=Y
         ELSE
           X2=X
           Y2=Y
         ENDIF
         LQ=KARL
         IF(LQ.EQ.1) THEN
           KARL=2*NINT(W)
         ELSEIF(IABS(LQ).EQ.2.AND.(LQ*W).LT.0.) THEN
           FAK=0.5D0
           KARL=3
         ENDIF
         IF(KARL.EQ.3.AND.STEG.LT.STEGBY) THEN
           W=(X2-X1)/(Y2-Y1)
           X=X1+W*(ALK-Y1)
         ELSE
           STEG=STEG*FAK
           X=X+STEG*W
         ENDIF
         AHPLUS=EXP(X)
       ENDIF
       IF(DONE) exit
    end do
    if (.not. DONE) then
      PH = 100.
    else
       HCO3=CO3*AHPLUS/AK2C
       IF(ICALC.EQ.4) THEN
         CTOT=H2CO3+HCO3+CO3
       ELSEIF(ICALC.EQ.1) THEN
         H2CO3=HCO3*AHPLUS/AK1C
         PCO2=H2CO3/AKP
       ENDIF
       PH=-LOG10(AHPLUS)
       ALKC=ALK-ALKB
    end if
  ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***
    Y=SQRT(PROD*(PROD-4.0D0*AKP*PCO2+4.0D0*CTOT))
    H2CO3=PCO2*AKP
    HCO3=(Y-PROD)/2.0D0
    CO3=CTOT-H2CO3-HCO3
    ALKC=HCO3+2.0D0*CO3
    AHPLUS=AK1C*H2CO3/HCO3
    PH=-LOG10(AHPLUS)
    ALKB=BTOT/(1.0D0+AHPLUS/AKB)
    ALK=ALKC+ALKB
  ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT, PH AND BTOT FIXED ***
    FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
    CO3=FACTOR*AK1C*AK2C
    HCO3=FACTOR*AK1C*AHPLUS
    H2CO3=FACTOR*AHPLUS*AHPLUS
    PCO2=H2CO3/AKP
    ALKC=HCO3+2.0D0*CO3
    ALKB=BTOT/(1.0D0+AHPLUS/AKB)
    ALK=ALKC+ALKB
  ELSEIF(ICALC.EQ.5) THEN
!         *** ALK, PH AND BTOT FIXED ***
    ALKB=BTOT/(1.0D0+AHPLUS/AKB)
    ALKC=ALK-ALKB
    HCO3=ALKC/(1.0D0+2.0D0*AK2C/AHPLUS)
    CO3=HCO3*AK2C/AHPLUS
    H2CO3=HCO3*AHPLUS/AK1C
    PCO2=H2CO3/AKP
    CTOT=H2CO3+HCO3+CO3
  ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2, PH AND BTOT FIXED ***
    ALKB=BTOT/(1.0D0+AHPLUS/AKB)
    H2CO3=PCO2*AKP
    HCO3=H2CO3*AK1C/AHPLUS
    CO3=HCO3*AK2C/AHPLUS
    CTOT=H2CO3+HCO3+CO3
    ALKC=HCO3+2.0D0*CO3
    ALK=ALKC+ALKB
  ENDIF
ELSE
  IF(ICALC.EQ.1) THEN
!         *** CTOT AND ALK FIXED ***
    TERM=4.0D0*ALK+CTOT*AKR-ALK*AKR
    Z=SQRT(TERM*TERM+4.0D0*(AKR-4.0D0)*ALK*ALK)
    CO3=(ALK*AKR-CTOT*AKR-4.0D0*ALK+Z)/(2.0D0*(AKR-4.0D0))
    HCO3=(CTOT*AKR-Z)/(AKR-4.0D0)
    H2CO3=CTOT-ALK+CO3
    PCO2=H2CO3/AKP
    PH=-LOG10(AK1C*H2CO3/HCO3)
  ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT AND PCO2 FIXED ***
    Y=SQRT(PROD*(PROD-4.0D0*AKP*PCO2+4.0D0*CTOT))
    H2CO3=PCO2*AKP
    HCO3=(Y-PROD)/2.0D0
    CO3=CTOT-H2CO3-HCO3
    ALK=HCO3+2.0D0*CO3
    PH=-LOG10(AK1C*H2CO3/HCO3)
  ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT AND PH FIXED ***
    FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
    CO3=FACTOR*AK1C*AK2C
    HCO3=FACTOR*AK1C*AHPLUS
    H2CO3=FACTOR*AHPLUS*AHPLUS
    PCO2=H2CO3/AKP
    ALK=HCO3+2.0D0*CO3
  ELSEIF(ICALC.EQ.4) THEN
!         *** ALK AND PCO2 FIXED ***
    TERM=SQRT((8.0D0*ALK+PROD)*PROD)
    CO3=ALK/2.0D0+PROD/8.0D0-TERM/8.0D0
    HCO3=-PROD/4.0D0+TERM/4.0D0
    H2CO3=PCO2*AKP
    CTOT=CO3+HCO3+H2CO3
    PH=-LOG10(AK1C*H2CO3/HCO3)
  ELSEIF(ICALC.EQ.5) THEN
!         *** ALK AND PH FIXED ***
    HCO3=ALK/(1.0D0+2.0D0*AK2C/AHPLUS)
    CO3=HCO3*AK2C/AHPLUS
    H2CO3=HCO3*AHPLUS/AK1C
    PCO2=H2CO3/AKP
    CTOT=H2CO3+HCO3+CO3
  ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2 AND PH FIXED ***
    H2CO3=PCO2*AKP
    HCO3=H2CO3*AK1C/AHPLUS
    CO3=HCO3*AK2C/AHPLUS
    CTOT=H2CO3+HCO3+CO3
    ALK=HCO3+2.0D0*CO3
  ENDIF
ENDIF


DO 120 II=1,NCONC
  CONCS(II)=CONCS2(II)
120    CONTINUE
RETURN
END SUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE CaCO3_Saturation (Tc, S, D, CO3, Om_cal, Om_arg)
!-----------------------------------------------------------------------

! Routine to calculate the saturation state of calcite and aragonite
! Inputs:
!               Tc      Temperature (C)
!               S       Salinity
!               D       Depth (m)
!               CO3     Carbonate ion concentration (mol.kg-1 ie /1D6)
!
! Outputs
!               Om_cal  Calite saturation
!               Om_arg  Aragonite saturation
!
! Intermediates
!               K_cal   Stoichiometric solubility product for calcite
!               K_arg   Stoichiometric solubility product for aragonite
!               Ca      Calcium 2+ concentration (mol.kg-1)
!               P       Pressure (bars)
!
! Source
!       Zeebe & Wolf-Gladrow 2001 following Mucci (1983)
!       with pressure corrections from Millero (1995)
!       Code tested against reference values given in Z & W-G
!       Built Jerry Blackford, 2008
!
!  use fabm_types

  IMPLICIT None
  AED_REAL :: Tc, Tk, Kelvin, S, D, Ca, CO3
  AED_REAL :: logKspc, Kspc, Om_cal
  AED_REAL :: logKspa, Kspa, Om_arg
  AED_REAL :: tmp1, tmp2, tmp3
  AED_REAL :: dV, dK, P, R

! setup
  Kelvin = 273.15
  Tk = Tc + Kelvin
  Ca = 0.01028    ! Currently oceanic mean value at S=25, needs refining)
  R = 83.131      !(cm3.bar.mol-1.K-1)
  P = D / 10.0    !pressure in bars

! calculate K for calcite
  tmp1 = -171.9065 - (0.077993*Tk) + (2839.319/Tk) + 71.595*log10(Tk)
  tmp2 = + (-0.77712 + (0.0028426*Tk) + (178.34/Tk))*SQRT(S)
  tmp3 = - (0.07711*S) + (0.0041249*(S**1.5))
  logKspc = tmp1 + tmp2 + tmp3
  Kspc = 10.0**logKspc

! correction for pressure for calcite
  IF ( D .GT. 0) THEN
    dV = -48.76 + 0.5304*Tc
    dK = -11.76/1.0D3 + (0.3692/1.0D3) * Tc
    tmp1 = -(dV/(R*Tk))*P + (0.5*dK/(R*Tk))*P*P
    Kspc = Kspc*exp(tmp1)
    logKspc = log10(Kspc)
  END IF

! calculate K for aragonite
  tmp1 = -171.945 - 0.077993*Tk + 2903.293 / Tk + 71.595* log10(Tk)
  tmp2 = + (-0.068393 + 0.0017276*Tk + 88.135/Tk)*SQRT(S)
  tmp3 = - 0.10018*S + 0.0059415*S**1.5
  logKspa = tmp1 + tmp2 + tmp3
  Kspa = 10.0**logKspa

! correction for pressure for aragonite
  IF ( D .GT. 0) THEN
    dV = -46.00 + 0.5304*Tc
    dK = -11.76/1.0D3 + (0.3692/1.0D3) * Tc
    tmp1 = -(dV/(R*Tk))*P + (0.5*dK/(R*Tk))*P*P
    Kspa = Kspa*exp(tmp1)
    logKspa = log10(Kspa)
  END IF

! calculate saturation states
  Om_cal = (CO3 * Ca) / Kspc
  Om_arg = (CO3 * Ca) / Kspa

RETURN
END SUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_carbon
