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
      INTEGER  :: id_temp, id_salt, id_wind, id_E_depth, id_E_tau, id_E_dens
      INTEGER  :: id_Fsed_amm, id_Fsed_nit, id_Fsed_n2o
      INTEGER  :: id_nitrif, id_denit, id_n2op, id_anammox, id_dnra
      INTEGER  :: id_sed_amm, id_sed_nit, id_sed_n2o
      INTEGER  :: id_atm_n2o

      !# Model parameters
      AED_REAL :: Rnitrif,Rdenit,Ranammox,Rn2o,Rdnra, &
                  Knitrif,Kdenit,Kanmx_nit,Kanmx_amm,Kdnra_oxy, &
                  theta_nitrif,theta_denit, &
                  Fsed_amm,Fsed_nit,Fsed_n2o,Ksed_amm,Ksed_nit,Ksed_n2o, &
                  theta_sed_amm,theta_sed_nit, &
                  Cn2o_atm
      LOGICAL  :: use_oxy, use_sed_model_amm, use_sed_model_nit
      LOGICAL  :: simN2O, simNitrfpH, simNitrfLight
      INTEGER  :: oxy_lim

     CONTAINS
         PROCEDURE :: define            => aed2_define_nitrogen
         PROCEDURE :: calculate         => aed2_calculate_nitrogen
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_nitrogen
         PROCEDURE :: calculate_surface => aed2_calculate_surface_nitrogen
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
   LOGICAL           :: simNitrfpH = .false.
   LOGICAL           :: simNitrfLight = .false.
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
   AED_REAL          :: Ranammox = 0.0
   AED_REAL          :: Rdnra = 0.0
   AED_REAL          :: Knitrif = 150.0
   AED_REAL          :: Kdenit = 150.0
   AED_REAL          :: Kanmx_nit = 150.0
   AED_REAL          :: Kanmx_amm = 150.0
   AED_REAL          :: Kdnra_oxy = 150.0
   AED_REAL          :: theta_nitrif = 1.0
   AED_REAL          :: theta_denit = 1.0
   AED_REAL          :: Fsed_amm = 3.5
   AED_REAL          :: Fsed_nit = 3.5
   AED_REAL          :: Fsed_n2o = 3.5
   AED_REAL          :: Ksed_amm = 30.0
   AED_REAL          :: Ksed_nit = 30.0
   AED_REAL          :: Ksed_n2o = 30.0
   AED_REAL          :: theta_sed_amm = 1.0
   AED_REAL          :: theta_sed_nit = 1.0
   AED_REAL          :: Cn2o_atm = 0.32 * 1e-6 !## current atmospheric N2O data (in ppm)

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
                    simN2O, Cn2o_atm, oxy_lim, Rn2o, Fsed_n2o, Ksed_n2o,  &
                    Ranammox, Rdnra, Kanmx_nit, Kanmx_amm, Kdnra_oxy,     &
                    simNitrfpH, simNitrfLight
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_nitrogen initialization"

   !-----------------------------------------------
   ! Read the namelist
   read(namlst,nml=aed2_nitrogen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_nitrogen'

   !-----------------------------------------------
   ! Store config options and parameter values in module's own derived type
   ! Note: all rates must be provided in values per day,
   ! and are converted for internal use as values per second.
   data%simNitrfpH = simNitrfpH
   data%simNitrfLight = simNitrfLight
   data%simN2O = simN2O
   data%oxy_lim = oxy_lim

   data%Rnitrif  = Rnitrif/secs_per_day
   data%Rdenit   = Rdenit/secs_per_day
   data%Rn2o     = Rn2o/secs_per_day
   data%Ranammox = Ranammox/secs_per_day
   data%Rdnra    = Rdnra/secs_per_day
   data%Knitrif  = Knitrif
   data%Kdenit   = Kdenit
   data%Kanmx_nit  = Kanmx_nit
   data%Kanmx_amm  = Kanmx_amm
   data%Kdnra_oxy  = Kdnra_oxy
   data%theta_nitrif = theta_nitrif
   data%theta_denit  = theta_denit
   data%Cn2o_atm   = Cn2o_atm

   data%Fsed_amm = Fsed_amm/secs_per_day
   data%Fsed_nit = Fsed_nit/secs_per_day
   data%Fsed_n2o = Fsed_n2o/secs_per_day
   data%Ksed_amm  = Ksed_amm
   data%Ksed_nit  = Ksed_nit
   data%Ksed_n2o  = Ksed_n2o
   data%theta_sed_amm = theta_sed_amm
   data%theta_sed_nit = theta_sed_nit

   !-----------------------------------------------
   ! Register state variables
   data%id_amm = aed2_define_variable('amm','mmol/m**3','ammonium',            &
                                    amm_initial,minimum=amm_min, maximum=amm_max)
   data%id_nit = aed2_define_variable('nit','mmol/m**3','nitrate',             &
                                    nit_initial,minimum=nit_min, maximum=nit_max)

   IF( simN2O ) THEN
     ! TODO check here to see if oxy is simulated if simN2O is also on
     IF (nitrif_reactant_variable .NE. '') THEN
       data%id_n2o = aed2_define_variable('n2o','mmol/m**3','nitrous oxide',     &
                                    n2o_initial,minimum=n2o_min, maximum=n2o_max)
      ELSE
        print *,'simN2O is true, however oxygen not provided so variable disabled'
      ENDIF
   ENDIF

   !-----------------------------------------------
   ! Register external state variable dependencies
   data%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) THEN
     data%id_oxy = aed2_locate_variable(nitrif_reactant_variable)
   ENDIF

   data%use_sed_model_amm = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_amm) &
     data%id_Fsed_amm = aed2_locate_global_sheet(Fsed_amm_variable)
   data%use_sed_model_nit = Fsed_amm_variable .NE. ''
   IF (data%use_sed_model_nit) &
     data%id_Fsed_nit = aed2_locate_global_sheet(Fsed_nit_variable)

     !-----------------------------------------------
   ! Register diagnostic variables
   data%id_nitrif = aed2_define_diag_variable('nitrif','mmol/m**3/d','nitrification rate')
   data%id_denit = aed2_define_diag_variable('denit','mmol/m**3/d','de-nitrification rate')
   data%id_anammox = aed2_define_diag_variable('anammox','mmol/m**3/d','anammox rate')
   data%id_dnra = aed2_define_diag_variable('dnra','mmol/m**3/d','dnra rate')
   data%id_sed_amm = aed2_define_sheet_diag_variable('sed_amm','mmol/m**2/d','ammonium sediment flux')
   data%id_sed_nit = aed2_define_sheet_diag_variable('sed_nit','mmol/m**2/d','nitrate sediment flux')
   IF( simN2O ) THEN
    data%id_n2op = aed2_define_diag_variable('n2oprod','mmol/m**3/d','n2o prod rate')
    data%id_atm_n2o = aed2_define_sheet_diag_variable('atm_n2o_flux','mmol/m**2/d','n2o atmospheric flux')
    data%id_sed_n2o = aed2_define_sheet_diag_variable('sed_n2o','mmol/m**2/d','n2o sediment flux')
   ENDIF

   !-----------------------------------------------
   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')
   data%id_wind = aed2_locate_global_sheet('wind_speed') ! @ 10 m above surface
   data%id_E_depth = aed2_locate_global('layer_ht')
   !data%id_E_vel = aed2_locate_global('velocity') ! needed for k600
   data%id_E_tau = aed2_locate_global('taub') ! tau to be converted to velocity
   data%id_E_dens = aed2_locate_global('density') ! density needed for tau-vel

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
   AED_REAL           :: nitrification,denitrification,anammox,dnra
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

   !-----------------------------------------------
   ! Set current (local) state variable values.
   amm = _STATE_VAR_(data%id_amm)    ! ammonium
   nit = _STATE_VAR_(data%id_nit)    ! nitrate
   IF( data%use_oxy ) THEN
      oxy = _STATE_VAR_(data%id_oxy) ! oxygen
   ELSE
      oxy = 300.0
   ENDIF
   pH = 7.

   !-----------------------------------------------
   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp) ! temperature

   !-----------------------------------------------
   ! Define process rates in units mmol N/m3/s

   !## nitrification
   nitrification = amm * fnitrif(data%use_oxy,data%Rnitrif,data%Knitrif,data%theta_nitrif,oxy,temp)
   IF( data%simNitrfpH ) nitrification = nitrification* NitrfpHFunction(pH)
   !IF( data%simNitrfLight ) nitrification = nitrification* NitrfLightFunction(I) ! Capone Pg 238

   !## de-nitrification
   denitrification = nit * fdenit(data%Rdenit, &
                            data%use_oxy,data%oxy_lim, &
                            data%Kdenit,data%theta_denit,Kno3, &
                            oxy,temp,nit)

   IF( data%simN2O ) THEN
     ! Babbin style model to capture intermediate N2O pool
     denit_n2o_prod = 0.5 * denitrification !* Xnc
     denit_n2o_cons = data%Rn2o * n2o * exp(-oxy/Kn2oc)
     nit_n2o_prod = zero_
     IF(oxy>Knev) nit_n2o_prod = ((aa/oxy)+bb)*nitrification
     !IF(oxy>Knev) nit_n2o_prod = ((aa/oxy)+bb)*remin*Xnc
   ENDIF

   !## anammox (NO2 + NH4 +CO2 -> N2); assuming nit ~ NO2
   anammox = zero_
   IF( data%use_oxy .AND. oxy < 1e-1*(1e3/32.) ) THEN
     anammox = data%Ranammox * nit/(data%Kanmx_nit+nit) * amm/(data%Kanmx_amm+amm)
   ENDIF

   !## dissasimilatory nitrate reduction to ammonia (DNRA)
   dnra = zero_
   IF( data%use_oxy ) THEN
     dnra = data%Rdnra * data%Kdnra_oxy/(data%Kdnra_oxy+oxy) * nit
   ENDIF

   !-----------------------------------------------
   ! Set temporal derivatives
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) -nitrification - anammox + dnra
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + nitrification - denitrification - anammox - dnra
   IF( data%simN2O ) &
    _FLUX_VAR_(data%id_n2o) = _FLUX_VAR_(data%id_n2o) + (denit_n2o_prod - denit_n2o_cons + nit_n2o_prod)

   ! If an externally maintained oxygen pool is linked, take nitrification from it
   IF (data%use_oxy) &
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-Xon*nitrification)

   !! add annamox use of dic & dnra creation + concsumption of DOC

   !-----------------------------------------------
   ! Export diagnostic variables
   _DIAG_VAR_(data%id_nitrif) = nitrification*secs_per_day
   _DIAG_VAR_(data%id_anammox) = anammox*secs_per_day
   _DIAG_VAR_(data%id_denit) = denitrification*secs_per_day
   _DIAG_VAR_(data%id_dnra) = dnra*secs_per_day
   IF( data%simN2O ) _DIAG_VAR_(data%id_n2op) = (denit_n2o_prod+nit_n2o_prod)*secs_per_day


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
   AED_REAL :: temp, salt, wind, vel, depth
   ! State
   AED_REAL :: n2o
   ! Temporary variables
   AED_REAL :: n2o_atm_flux = zero_
   AED_REAL :: Cn2o_air = zero_    !N2O in the air phase
   AED_REAL :: kn2o_trans = zero_
   AED_REAL :: windHt
   AED_REAL :: f_pres  = 1.0      ! Pressure correction function (unsued as only applicable at high altitudes)
!
!-------------------------------------------------------------------------------
!BEGIN

   IF( .NOT.data%simN2O ) RETURN

   !-----------------------------------------------
   !Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10 m above surface (m/s)
   windHt = 10.
   depth = MAX( _STATE_VAR_(data%id_E_depth), one_ )
!   vel  = SQRT(_STATE_VAR_(data%id_E_tau)/_STATE_VAR_(data%id_E_dens))
!   vel = vel/0.41 * log(depth/0.01)
   vel = 0.0001
   !-----------------------------------------------
   ! Retrieve current (local) state variable values.
   n2o = _STATE_VAR_(data%id_n2o)! Concentration of N2O in surface layer

   !-----------------------------------------------
   ! Get the surface piston velocity  (THIS IS BASED ON 2D FLOWS)
   !kn2o_trans = aed2_gas_piston_velocity(windHt,wind,temp,salt,schmidt_model=?)
   kn2o_trans = 0.77*( (vel**0.5)*(depth**(-0.5)) + 0.266*wind**2. ) / 3.6e5  ! transfer velocity k of Ho et al. 2016

   !-----------------------------------------------
   ! First get the N2O concentration in the air-phase at interface
   ! C_N2O = F x' P ; Capone (2008) pg 56
   F_pres = 1.0
   Cn2o_air = data%Cn2o_atm
   Cn2o_air = Cn2o_air * aed2_n2o_sat(salt,temp)

   !-----------------------------------------------
   ! Get the N2O flux:    [ mmol/m2/s = m/s * mmol/m3 ]
   n2o_atm_flux = kn2o_trans * (Cn2o_air - n2o) * F_pres

   !-----------------------------------------------
   ! Transfer surface exchange value to AED2 (mmmol/m2), converted by driver.
   _FLUX_VAR_T_(data%id_n2o) = n2o_atm_flux

   !-----------------------------------------------
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
   AED_REAL :: amm,nit,n2o,oxy
   ! Temporary variables
   AED_REAL :: amm_flux, nit_flux, n2o_flux
   AED_REAL :: Fsed_amm, Fsed_nit, Fsed_n2o
   AED_REAL :: fTa, fTo
!
!-------------------------------------------------------------------------------
!BEGIN

   !-----------------------------------------------
   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   !! Retrieve current (local) state variable values.
   !amm = _STATE_VAR_(data%id_amm) ! ammonium
   !nit = _STATE_VAR_(data%id_nit) ! nitrate

   !-----------------------------------------------
   ! Set the maximum flux (@20C) to use in this cell
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
   Fsed_n2o = data%Fsed_n2o

   !-----------------------------------------------
   ! Compute temperature scaling
   fTa = data%theta_sed_amm**(temp-20.0)
   fTo = data%theta_sed_nit**(temp-20.0)

   !-----------------------------------------------
   ! Compute actual flux based on oxygen and temperature
   n2o_flux = zero_
   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      amm_flux = Fsed_amm * data%Ksed_amm/(data%Ksed_amm+oxy) * fTa
      nit_flux = Fsed_nit *           oxy/(data%Ksed_nit+oxy) * fTo
      IF( data%simN2O ) n2o_flux = Fsed_n2o * data%Ksed_n2o/(data%Ksed_n2o+oxy) * fTa
   ELSE
      ! Sediment flux dependent on temperature only.
      amm_flux = Fsed_amm * fTa
      nit_flux = Fsed_nit * fTo
   ENDIF

   !-----------------------------------------------
   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + amm_flux
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + nit_flux
   IF( data%simN2O ) _FLUX_VAR_(data%id_n2o) = _FLUX_VAR_(data%id_n2o) + n2o_flux

   !-----------------------------------------------
   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_amm) = _FLUX_VAR_B_(data%id_ben_amm) + (-amm_flux/secs_per_day)

   !-----------------------------------------------
   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_amm) = amm_flux*secs_per_day
   _DIAG_VAR_S_(data%id_sed_nit) = nit_flux*secs_per_day
   IF( data%simN2O ) _DIAG_VAR_S_(data%id_sed_n2o) = n2o_flux*secs_per_day

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




!###############################################################################
 FUNCTION NitrfLightFunction(light) RESULT(limitation)
   !----------------------------------------------------------------------------
   ! Dependence of nitrification rate on light
   !----------------------------------------------------------------------------
   !-- Incoming
   AED_REAL, INTENT(IN) :: light              ! pH in the water column
   !-- Returns the salinity function
   AED_REAL :: limitation
   !-- Local

     limitation = one_   ! TO BE COMPLETED - SEE CAPONE 2008

 END FUNCTION NitrfLightFunction
!------------------------------------------------------------------------------!


END MODULE aed2_nitrogen
