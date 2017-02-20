!###############################################################################
!#                                                                             #
!# aed2_nitrogen.F90                                                           #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created 9 May 2011                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_nitrogen
!-------------------------------------------------------------------------------
! aed2_nitrogen --- nitrogen biogeochemical model
!
! Nitrogen module contains equations for nitrification and deitrification
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,  ONLY: aed2_gas_piston_velocity, aed2_n2o_sat

   IMPLICIT NONE

   PRIVATE    ! By default make everything private
!
   PUBLIC aed2_nitrogen_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_nitrogen_data_t
      !# Variable identifiers
      INTEGER  :: id_nit, id_amm, id_n2o
      INTEGER  :: id_oxy, id_denit_product
      INTEGER  :: id_temp, id_salt, id_wind, id_E_depth
      INTEGER  :: id_Fsed_amm, id_Fsed_nit
      INTEGER  :: id_nitrif, id_denit, id_n2op
      INTEGER  :: id_sed_amm, id_sed_nit
      INTEGER  :: id_atm_n2o

      !# Model parameters
      AED_REAL :: Rnitrif,Rdenit,Fsed_amm,Fsed_nit,Knitrif,Kdenit,Ksed_amm,Ksed_nit, &
                  theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit,Cn2o_atm, Rn2o
      LOGICAL  :: use_oxy,use_no2,use_sed_model_amm, use_sed_model_nit
      LOGICAL  :: simN2O, simNitrfpH
      INTEGER  :: oxy_lim

     CONTAINS
         PROCEDURE :: define            => aed2_define_nitrogen
         PROCEDURE :: calculate         => aed2_calculate_nitrogen
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_nitrogen
!        PROCEDURE :: mobility          => aed2_mobility_nitrogen
!        PROCEDURE :: light_extinction  => aed2_light_extinction_nitrogen
!        PROCEDURE :: delete            => aed2_delete_nitrogen

   END TYPE

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_nitrogen(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_nitrogen_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER           :: status
   INTEGER           :: oxy_lim = 1
   LOGICAL           :: simN2O = .false.
   AED_REAL          :: nit_initial=4.5
   AED_REAL          :: nit_min=zero_
   AED_REAL          :: nit_max=1e6
   AED_REAL          :: amm_initial=4.5
   AED_REAL          :: amm_min=zero_
   AED_REAL          :: amm_max=1e6
   AED_REAL          :: n2o_initial=0.5
   AED_REAL          :: n2o_min=zero_
   AED_REAL          :: n2o_max=1e6
   AED_REAL          :: Rnitrif = 0.01
   AED_REAL          :: Rdenit = 0.01
   AED_REAL          :: Rn2o = 0.0015
   AED_REAL          :: Fsed_amm = 3.5
   AED_REAL          :: Fsed_nit = 3.5
   AED_REAL          :: Knitrif = 150.0
   AED_REAL          :: Kdenit = 150.0
   AED_REAL          :: Ksed_amm = 30.0
   AED_REAL          :: Ksed_nit = 30.0
   AED_REAL          :: theta_nitrif = 1.0
   AED_REAL          :: theta_denit = 1.0
   AED_REAL          :: theta_sed_amm = 1.0
   AED_REAL          :: theta_sed_nit = 1.0
   AED_REAL          :: Cn2o_atm = 1.0
   CHARACTER(len=64) :: nitrif_reactant_variable=''
   CHARACTER(len=64) :: denit_product_variable=''
   CHARACTER(len=64) :: Fsed_amm_variable=''
   CHARACTER(len=64) :: Fsed_nit_variable=''


   NAMELIST /aed2_nitrogen/ nit_initial,nit_min, nit_max,                 &
                    amm_initial, amm_min, amm_max,                        &
                    n2o_initial, n2o_min, n2o_max,                        &
                    Rnitrif,Rdenit,Fsed_amm,Fsed_nit,                     &
                    Knitrif,Kdenit,Ksed_amm,Ksed_nit,                     &
                    theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit, &
                    nitrif_reactant_variable,denit_product_variable,      &
                    Fsed_amm_variable, Fsed_nit_variable,                 &
                    simN2O, Cn2o_atm, oxy_lim, Rn2o
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_nitrogen initialization"

   ! Read the namelist
   read(namlst,nml=aed2_nitrogen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_nitrogen'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   data%simNitrfpH = .false.
   data%simN2O = simN2O
   data%oxy_lim = oxy_lim

   data%Rnitrif  = Rnitrif/secs_per_day
   data%Rdenit   = Rdenit/secs_per_day
   data%Rn2o     = Rn2o/secs_per_day
   data%Fsed_amm = Fsed_amm/secs_per_day
   data%Fsed_nit = Fsed_nit/secs_per_day
   data%Knitrif  = Knitrif
   data%Kdenit   = Kdenit
   data%Ksed_amm  = Ksed_amm
   data%Ksed_nit  = Ksed_nit
   data%theta_nitrif = theta_nitrif
   data%theta_denit  = theta_denit
   data%theta_sed_amm = theta_sed_amm
   data%theta_sed_nit = theta_sed_nit
   data%Cn2o_atm = Cn2o_atm


   ! Register state variables
   data%id_amm = aed2_define_variable('amm','mmol/m**3','ammonium',            &
                                    amm_initial,minimum=amm_min, maximum=amm_max)
   data%id_nit = aed2_define_variable('nit','mmol/m**3','nitrate',             &
                                    nit_initial,minimum=nit_min, maximum=nit_max)


   IF( simN2O ) THEN
     data%id_n2o = aed2_define_variable('n2o','mmol/m**3','nitrous oxide',     &
                                    n2o_initial,minimum=n2o_min, maximum=n2o_max)
   ENDIF


   ! Register external state variable dependencies
   data%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) THEN
     data%id_oxy = aed2_locate_variable(nitrif_reactant_variable)
   ENDIF

   !data%use_no2 = denit_product_variable .NE. '' !This means n2 module switched on
   !IF (data%use_no2) data%id_denit_product = aed2_locate_variable(denit_product_variable)

   data%use_sed_model_amm = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_amm) &
     data%id_Fsed_amm = aed2_locate_global_sheet(Fsed_amm_variable)
   data%use_sed_model_nit = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_nit) &
     data%id_Fsed_nit = aed2_locate_global_sheet(Fsed_nit_variable)


   ! Register diagnostic variables
   data%id_nitrif = aed2_define_diag_variable('nitrif','mmol/m**3/d', &
                                                         'nitrification rate')
   data%id_denit = aed2_define_diag_variable('denit','mmol/m**3/d', &
                                                         'de-nitrification rate')
   data%id_n2op = aed2_define_diag_variable('n2oprod','mmol/m**3/d', &
                                                       'n2o prod rate')
   data%id_sed_amm = aed2_define_sheet_diag_variable('sed_amm','mmol/m**2/d', &
                                                         'ammonium sediment flux')
   data%id_sed_nit = aed2_define_sheet_diag_variable('sed_nit','mmol/m**2/d', &
                                                         'nitrate sediment flux')
   data%id_atm_n2o = aed2_define_sheet_diag_variable('atm_n2o_flux','mmol/m**2/d', &
                                                     'N2O atmospheric flux')

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')
   data%id_wind = aed2_locate_global_sheet('wind_speed') ! Wind speed at 10 m above surface (m/s)
   data%id_E_depth = aed2_locate_global('layer_ht')
   ! check here to see if oxy is simulated if simN2O is also on


END SUBROUTINE aed2_define_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_nitrogen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_nitrogen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: amm,nit,n2o,oxy,temp,pH !State variables
   AED_REAL           :: nitrification,denitrification
   AED_REAL           :: denit_n2o_prod,denit_n2o_cons,nit_n2o_prod
   AED_REAL,PARAMETER :: Xon = 3. !ratio of oxygen to nitrogen utilised during nitrification
   AED_REAL,PARAMETER :: Knev = 3. !Nevison nitrification O2 threshold
   AED_REAL,PARAMETER :: aa = 0.26 !Nevison nitrification parameter
   AED_REAL,PARAMETER :: bb = -0.0006 !Nevison nitrification parameter
   AED_REAL,PARAMETER :: Xnc = 16./106. !OM stoichiomtery
   AED_REAL,PARAMETER :: Kn2oc = 0.3 !N2O consumption O2 poisoning
   AED_REAL,PARAMETER :: Kno3 = 5.0 !Denit NO3 half-sat
   AED_REAL,PARAMETER :: remin = 0.0 !reminerlaisation of OM (get from OM)
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
   amm = _STATE_VAR_(data%id_amm)! ammonium
   nit = _STATE_VAR_(data%id_nit)! nitrate
   IF (data%use_oxy) THEN ! & use_oxy
      oxy = _STATE_VAR_(data%id_oxy)! oxygen
   ELSE
      oxy = zero_
   ENDIF
   pH = 7.

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp) ! temperature

   ! Define process rates in units mmol N/m3/s
   nitrification = fnitrif(data%use_oxy,data%Rnitrif,data%Knitrif,data%theta_nitrif,oxy,temp)
   IF( data%simNitrfpH ) nitrification = nitrification* NitrfpHFunction(pH)

   denitrification = fdenit(data%Rdenit, &
                            data%use_oxy,data%oxy_lim, &
                            data%Kdenit,data%theta_denit,Kno3, &
                            oxy,temp,nit)

   IF( data%simN2O ) THEN
     denit_n2o_prod = 0.5 * nit * denitrification !* Xnc
     denit_n2o_cons = data%Rn2o * n2o * exp(-oxy/Kn2oc)
     nit_n2o_prod = zero_
     IF(oxy>Knev) nit_n2o_prod = ((aa/oxy)+bb)*nitrification*amm
     !IF(oxy>Knev) nit_n2o_prod = ((aa/oxy)+bb)*remin*Xnc
   ENDIF

   ! Set temporal derivatives
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (-amm*nitrification)
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + (amm*nitrification - nit*denitrification)

   IF( data%simN2O ) &
    _FLUX_VAR_(data%id_n2o) = _FLUX_VAR_(data%id_n2o) + (denit_n2o_prod - denit_n2o_cons + nit_n2o_prod)

   ! If an externally maintained oxygen pool is linked, take nitrification from it
   IF (data%use_oxy) &
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-Xon*amm*nitrification)

   ! Export diagnostic variables
   _DIAG_VAR_(data%id_nitrif) = amm*nitrification*secs_per_day
   _DIAG_VAR_(data%id_denit) = nit*denitrification*secs_per_day
   _DIAG_VAR_(data%id_n2op) = (denit_n2o_prod+nit_n2o_prod)*secs_per_day


END SUBROUTINE aed2_calculate_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface_nitrogen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-water exchange for the aed nitrogen model (N2O)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind

   ! State
   AED_REAL :: n2o

   ! Temporary variables
   AED_REAL :: n2o_atm_flux = zero_
   AED_REAL :: Cn2o_air = zero_    !N2O in the air phase
   AED_REAL :: kn2o_trans = zero_
   AED_REAL :: windHt, vel, depth
   AED_REAL :: f_pres  = 1.0      ! Pressure correction function only applicable at high altitudes
!
!-------------------------------------------------------------------------------
!BEGIN

   IF( .NOT.data%simN2O ) RETURN

   !Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10 m above surface (m/s)
   windHt = 10.
   vel  = 0.
   depth = _STATE_VAR_(data%id_E_depth)

    ! Retrieve current (local) state variable values.
   n2o = _STATE_VAR_(data%id_n2o)! Concentration of N2O in surface layer

   !kn2o_trans = aed2_gas_piston_velocity(windHt,wind,temp,salt)
   kn2o_trans = 0.77*( (vel**0.5)*(depth**(-0.5)) + 0.266*wind**2. )  ! transfer velocity k of Ho et al. 2016

   ! First get the oxygen concentration in the air phase at interface
   ! Taken from Riley and Skirrow (1974)
   f_pres = 1.0
   !Cn2o_air = data%Cn2o_atm
   Cn2o_air = aed2_n2o_sat(salt,temp)

   ! Get the oxygen flux
   n2o_atm_flux = kn2o_trans * (Cn2o_air - n2o)

   ! Transfer surface exchange value to AED2 (mmmol/m2) converted by driver.
   _FLUX_VAR_T_(data%id_n2o) = n2o_atm_flux

   ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).
   _DIAG_VAR_S_(data%id_atm_n2o) = n2o_atm_flux

END SUBROUTINE aed2_calculate_surface_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_nitrogen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_nitrogen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: amm,nit,oxy

   ! Temporary variables
   AED_REAL :: amm_flux,nit_flux
   AED_REAL :: Fsed_amm, Fsed_nit
!
!-------------------------------------------------------------------------------
!BEGIN

  ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   amm = _STATE_VAR_(data%id_amm) ! ammonium
   nit = _STATE_VAR_(data%id_nit) ! nitrate

   IF (data%use_sed_model_amm) THEN
      Fsed_amm = _STATE_VAR_S_(data%id_Fsed_amm)
   ELSE
      Fsed_amm = data%Fsed_amm
   ENDIF
   IF (data%use_sed_model_nit) THEN
      Fsed_nit = _STATE_VAR_S_(data%id_Fsed_nit)
   ELSE
      Fsed_nit = data%Fsed_nit
   ENDIF

   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      amm_flux = Fsed_amm * data%Ksed_amm/(data%Ksed_amm+oxy) * (data%theta_sed_amm**(temp-20.0))
      nit_flux = Fsed_nit * oxy/(data%Ksed_nit+oxy) * (data%theta_sed_nit**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      oxy = 0.
      amm_flux = Fsed_amm * (data%theta_sed_amm**(temp-20.0))
      nit_flux = Fsed_nit * (data%theta_sed_nit**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   !_SET_BOTTOM_FLUX_(data%id_amm,amm_flux/secs_per_day)
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (amm_flux)
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + (nit_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_amm) = _FLUX_VAR_B_(data%id_ben_amm) + (-amm_flux/secs_per_day)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_amm) = amm_flux*secs_per_day
   _DIAG_VAR_S_(data%id_sed_nit) = nit_flux*secs_per_day
END SUBROUTINE aed2_calculate_benthic_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fnitrif(use_oxy,Rnitrif,Knitrif,theta_nitrif,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rnitrif,Knitrif,theta_nitrif,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fnitrif = Rnitrif * oxy/(Knitrif+oxy) * (theta_nitrif**(temp-20.0))
   ELSE
      fnitrif = Rnitrif * (theta_nitrif**(temp-20.0))
   ENDIF
END FUNCTION fnitrif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdenit(Rdenit,use_oxy,oxy_lim,Kdenit,theta_denit,Kdenitnit,oxy,temp,nit)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL, INTENT(in) :: use_oxy
   INTEGER, INTENT(in) :: oxy_lim
   AED_REAL,INTENT(in) :: Rdenit,Kdenit,theta_denit,oxy,temp,nit,Kdenitnit
   AED_REAL :: fT, fDO, fNO3
!
!-------------------------------------------------------------------------------
!BEGIN

   fT = (theta_denit**(temp-20.0))

   fDO = one_
   IF (use_oxy) THEN
     IF (oxy_lim == 1) fDO = Kdenit/(Kdenit+oxy)
     IF (oxy_lim == 2) fDO = exp(-oxy/Kdenit)
   ENDIF

   IF(Kdenitnit==zero_)THEN
     fNO3 = one_
   ELSE
     fNO3 = nit/(Kdenitnit+nit)
   ENDIF

   fdenit = Rdenit * fDO * fT * fNO3

END FUNCTION fdenit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
 FUNCTION NitrfpHFunction(pH) RESULT(limitation)
   !------------------------------------------------------------------------------!
   ! Dependence of nitrification rate on pH                                       !
   !------------------------------------------------------------------------------!
   !-- Incoming
   AED_REAL, INTENT(IN) :: pH              ! pH in the water column
   !-- Returns the salinity function
   AED_REAL :: limitation
   !-- Local
   AED_REAL, PARAMETER :: NITpHOptMin = 7.1        ! min pH of optimum range
   AED_REAL, PARAMETER :: NITpHOptMax = 7.9        ! max pH of optimum range
   AED_REAL, PARAMETER :: NITpHTolMax = 9.0        ! upper pH tolerance
   AED_REAL, PARAMETER :: NITpHTolMin = 5.5        ! lower pH tolerance
   AED_REAL :: tmp1,tmp2

   !____________________________________________________________________!
   !                                                                    !
   !  =1                    .!---------------!.                         !
   !                      !                     !                       !
   !                    !                         !                     !
   !                   !                           !                    !
   !                  !                             !                   !
   !                 !                               !                  !
   !                !                                 !                 !
   !_______________!___________________________________!________________!
   !               !         !               !         !                !
   !                    NITpHOptMin      NITpHOptMax                    !
   !           NITpHTolMax                         NITpHTolMax          !
   !____________________________________________________________________!


   !---------------------------------------------------------------------------!
   ! pH is within the tolerance; no limitation.                                !
   IF(pH >= NITpHOptMin .AND. pH <= NITpHOptMax) THEN
     limitation = one_
   ENDIF

   ! pH is greater than the upper bound of the optimum region
   IF(pH > NITpHOptMax) THEN
     limitation = (-pH*pH+2.0*NITpHOptMax*pH -                      &
                  2.0*NITpHOptMax*NITpHTolMax+NITpHTolMax*NITpHTolMax)/  &
                                         ((NITpHOptMax)*(NITpHOptMax))
   ENDIF

   ! pH is less than the lower bound of optimum region (NITpHOptMin)
   ! but greater than minimum tolerance NITpHTolMin
   tmp1 = zero_
   tmp2 = zero_
   IF(pH < NITpHOptMin .AND. pH > NITpHTolMin) THEN
     tmp1 = pH-NITpHTolMin
     tmp2 = NITpHOptMin-NITpHTolMin
     limitation = (2.0*tmp1/NITpHOptMin-(tmp1*tmp1/(NITpHOptMin*NITpHOptMin))) &
                / (2.0*tmp2/NITpHOptMin-(tmp2*tmp2/(NITpHOptMin*NITpHOptMin)))
   ENDIF

   ! pH is less than the NITpHTolMin
   IF(pH <= NITpHTolMin) THEN
     limitation = zero_
   ENDIF

   ! ensure we don't go negative
   IF(limitation <= zero_) THEN
     limitation = zero_
   ENDIF

 END FUNCTION NitrfpHFunction
!------------------------------------------------------------------------------!


END MODULE aed2_nitrogen
