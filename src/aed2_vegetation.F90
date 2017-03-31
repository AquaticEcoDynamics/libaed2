!###############################################################################
!#                                                                             #
!# aed2_vegetation.F90                                                         #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!#                                                                             #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created January 2015                                                        #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#define _PHYLEN_ 3
#define _PHYMOD_ 'PHY'
#define _OGMPOC_ 'OGM'

MODULE aed2_vegetation
!-------------------------------------------------------------------------------
!  aed2_vegetation --- multi-group vegetation model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util !,ONLY : find_free_lun,aed2_bio_temp_function,fTemp_function,qsort
   USE aed2_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed2_vegetation_data_t



   TYPE macrophyte_data
     INTEGER       :: growthForm
     CHARACTER(64) :: m_name
     AED_REAL      :: m0
     AED_REAL      :: R_growth
     INTEGER       :: fT_Method
     AED_REAL      :: theta_growth
     AED_REAL      :: T_std
     AED_REAL      :: T_opt
     AED_REAL      :: T_max
     INTEGER       :: lightModel
     AED_REAL      :: I_K
     AED_REAL      :: I_S
     AED_REAL      :: KeMAC
     AED_REAL      :: f_pr
     AED_REAL      :: R_resp
     AED_REAL      :: theta_resp
     INTEGER       :: salTol
     AED_REAL      :: S_bep
     AED_REAL      :: S_maxsp
     AED_REAL      :: S_opt
     AED_REAL      :: K_CD
     AED_REAL      :: f_bg
     AED_REAL      :: k_omega
     AED_REAL      :: Xcc
     AED_REAL      :: K_N
     AED_REAL      :: X_ncon
     AED_REAL      :: K_P
     AED_REAL      :: X_pcon
     !
     !AED_REAL      :: root ...
   END TYPE

   TYPE,extends(aed2_model_data_t) :: aed2_vegetation_data_t
      !# Variable identifiers
      INTEGER  :: id_veg(MAX_VEG_TYPES)
      INTEGER  :: id_vegfrac(MAX_VEG_TYPES)
      INTEGER  :: id_l_litter_land, id_l_litter_water

      !# Environmental variables
      INTEGER :: id_E_rain, id_E_area, id_E_matz, id_E_bath, id_E_salt, id_E_nearlevel, id_E_airtemp, id_E_I0
      AED_REAL,ALLOCATABLE :: active_zones(:)

      !# Diagnostic variables
      INTEGER :: id_l_lai, id_gpp, id_lfall, id_resp, id_tveg
      INTEGER :: id_atm_co2

      !# Dependant variable IDs
      INTEGER :: id_l_depth, id_l_Sb, id_l_St, id_l_Ssat, id_l_theta, id_l_Stop, id_l_capz
      INTEGER :: id_l_phreatic, id_l_qss, id_l_qse, id_l_qcap, id_l_qper, id_l_wt
      INTEGER :: id_l_soiltemp

      !# Model parameters
      INTEGER  :: num_veg, n_zones
      LOGICAL  :: simLitterLand, simLitterWater
      LOGICAL  :: simVegFeedback, simStaticBiomass, simVegFrac
      TYPE(macrophyte_data),DIMENSION(:),ALLOCATABLE :: vegdata
      AED_REAL :: int_max, lai_max

     CONTAINS
         PROCEDURE :: define             => aed2_define_vegetation
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_vegetation
!        PROCEDURE :: calculate_dry      => aed2_calculate_riparian_vegetation
!        PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_vegetation
!        PROCEDURE :: mobility           => aed2_mobility_vegetation
!        PROCEDURE :: light_extinction   => aed2_light_extinction_vegetation
         PROCEDURE :: rain_loss          => aed2_rain_loss_vegetation
         PROCEDURE :: light_shading      => aed2_light_shading_vegetation
!        PROCEDURE :: delete             => aed2_delete_vegetation

   END TYPE

   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed2_vegetation_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count   !Number of vegetation groups
   INTEGER,INTENT(in)          :: list(*) !List of vegetation groups to simulate
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   AED_REAL :: Pveg_prey(MAX_ZOOP_PREY)

   TYPE(macrophyte_data)  :: veg(MAX_ZOOP_TYPES)
   NAMELIST /vegetation_data/ veg
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD',iostat=status)
    IF (status /= 0) STOP 'Error opening vegdata_params namelist file'
    read(tfil,nml=vegetation_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist vegdata_params'

    data%num_veg = count
    allocate(data%vegdata(count))

    !# loop through each group, read parameters & define
    DO i=1,count
      data%vegdata(i)%growthForm   = 1
      data%vegdata(i)%m_name       = veg(list(i))%m_name
      data%vegdata(i)%m0           = veg(list(i))%m0
      data%vegdata(i)%R_growth     = veg(list(i))%R_growth/secs_per_day
      data%vegdata(i)%fT_Method    = veg(list(i))%fT_Method
      data%vegdata(i)%theta_growth = veg(list(i))%theta_growth
      data%vegdata(i)%T_std        = veg(list(i))%T_std
      data%vegdata(i)%T_opt        = veg(list(i))%T_opt
      data%vegdata(i)%T_max        = veg(list(i))%T_max
      data%vegdata(i)%lightModel   = veg(list(i))%lightModel
      data%vegdata(i)%I_K          = veg(list(i))%I_K
      data%vegdata(i)%I_S          = veg(list(i))%I_S
      data%vegdata(i)%KeMAC        = veg(list(i))%KeMAC
      data%vegdata(i)%f_pr         = veg(list(i))%f_pr
      data%vegdata(i)%R_resp       = veg(list(i))%R_resp/secs_per_day
      data%vegdata(i)%theta_resp   = veg(list(i))%theta_resp
      data%vegdata(i)%salTol       = veg(list(i))%salTol
      data%vegdata(i)%S_bep        = veg(list(i))%S_bep
      data%vegdata(i)%S_maxsp      = veg(list(i))%S_maxsp
      data%vegdata(i)%S_opt        = veg(list(i))%S_opt
      data%vegdata(i)%K_CD         = veg(list(i))%K_CD
      data%vegdata(i)%f_bg         = veg(list(i))%f_bg
      data%vegdata(i)%k_omega      = veg(list(i))%k_omega
      data%vegdata(i)%Xcc          = veg(list(i))%Xcc
      data%vegdata(i)%K_N          = veg(list(i))%K_N
      data%vegdata(i)%X_ncon       = veg(list(i))%X_ncon
      data%vegdata(i)%K_P          = veg(list(i))%K_P
      data%vegdata(i)%X_pcon       = veg(list(i))%X_pcon

      ! Register group as a state variable
      data%id_veg(i) = aed2_define_sheet_variable(                             &
                              veg(list(i))%m_name,                             &
                              'mmolC/m**2', 'vegetation biomass',              &
                              veg(list(i))%m0,                                 &
                              minimum=zero_)

      data%id_vegfrac(i) = aed2_define_sheet_diag_variable(                    &
                              TRIM(veg(list(i))%m_name)//'_cover' ,            &
                              '%', 'vegetation cover')
    ENDDO
!
END SUBROUTINE aed2_vegetation_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_define_vegetation(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the vegetation biogeochemical model
!
!  Here, the aed2_vegetation namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: veg_i
   INTEGER  :: status
   INTEGER  :: num_veg
   INTEGER  :: the_veg(MAX_VEG_TYPES)
   LOGICAL  :: simVegFeedback, simStaticBiomass, simVegFrac
   INTEGER  :: n_zones = 0, active_zones(MAX_ZONES), i
   CHARACTER(len=64)  :: lai_link_variable = 'LND_lai'
   CHARACTER(len=64)  :: litter_target_variable_land  = ''
   CHARACTER(len=64)  :: litter_target_variable_water = ''
   CHARACTER(len=128) :: dbase='aed2_vegetation_pars.nml'
   AED_REAL :: int_max, lai_max

   NAMELIST /aed2_vegetation/ num_veg, the_veg, dbase, n_zones, active_zones,  &
                              litter_target_variable_land, litter_target_variable_water, &
                              simVegFeedback, simStaticBiomass, simVegFrac, &
                              int_max, lai_max, lai_link_variable

!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed2_vegetation initialization"

   ! Read the namelist
   read(namlst,nml=aed2_vegetation,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_vegetation'

   data%num_veg = 0
   data%int_max = int_max; data%lai_max = lai_max

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
       ALLOCATE(data%active_zones(n_zones))
       DO i=1,n_zones
          data%active_zones(i) = active_zones(i)
       ENDDO
   ENDIF

   ! Store species specific parameter values in modules own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed2_vegetation_load_params to values per second.
   CALL aed2_vegetation_load_params(data, dbase, num_veg, the_veg)

   ! Register link to nutrient pools, if variable names are provided in namelist.
   data%simLitterLand = litter_target_variable_land .NE. ''
   IF (data%simLitterLand) THEN
     data%id_l_litter_land = aed2_locate_variable(litter_target_variable_land)
   ENDIF
   data%simLitterWater = litter_target_variable_water .NE. ''
   IF (data%simLitterWater) THEN
     data%id_l_litter_water = aed2_locate_variable(litter_target_variable_water)
   ENDIF


   ! Register diagnostic variables
   IF ( lai_link_variable .EQ. '' ) &
   data%id_l_lai = aed2_define_sheet_diag_variable('lai','mmolC/m**2/d',  'vegetation LAI')
   data%id_gpp   = aed2_define_sheet_diag_variable('gpp','mmolC/m**2/d',  'net vegetation productivity')
   data%id_resp  = aed2_define_sheet_diag_variable('resp','mmolC/m**2/d','net vegetation respiration')
   data%id_lfall = aed2_define_sheet_diag_variable('lfall','mmolC/m**2/d','net vegetation mortality')
   data%id_atm_co2 = aed2_define_sheet_diag_variable('atm_co2','mmolC/m**2/d',  'co2 exchange to the atmosphere')
   data%id_tveg = aed2_define_sheet_diag_variable('veg','mmolC/m**2',  'total veg biomass')

   ! Register module dependencies
   IF ( .NOT. lai_link_variable .EQ. '' ) &
   data%id_l_lai      = aed2_locate_global_sheet(TRIM(lai_link_variable))
   data%id_l_depth    = aed2_locate_global_sheet('LND_depth')         !,'m','soil depth (to datum)')
   data%id_l_phreatic = aed2_locate_global_sheet('LND_phreatic')      !,'m','depth of phreatic surface below surface')
   data%id_l_wt       = aed2_locate_global_sheet('LND_wt')            !,'m','depth of phreatic surface below surface')
   data%id_l_Sb       = aed2_locate_global_sheet('LND_Sb')            !,'mm','total capacity for water storage')
   data%id_l_St       = aed2_locate_global_sheet('LND_St')            ! 'mm', 'total soil water storage at time t')
   data%id_l_Ssat     = aed2_locate_global_sheet('LND_Ssat')          !,'mm','saturated zone soil water storage')
   data%id_l_Stop     = aed2_locate_global_sheet('LND_Stop')          !,'mm','unsaturated storage capacity')
   data%id_l_capz     = aed2_locate_global_sheet('LND_capz')          !,'mm','unsaturated storage capacity')
   data%id_l_theta    = aed2_locate_global_sheet('LND_theta')         !,'-','unsaturated moisture content')
   data%id_l_qss      = aed2_locate_global_sheet('LND_qss')           !,'mm','sat zone seepage')
   data%id_l_qse      = aed2_locate_global_sheet('LND_qs')            !,'mm','surface runoff')
   data%id_l_qcap     = aed2_locate_global_sheet('LND_qcap')          !,'mm','capillarity')
   data%id_l_qper     = aed2_locate_global_sheet('LND_qper')          !,'mm','recharge')
   data%id_l_soiltemp = aed2_locate_global_sheet('LND_soiltemp_10cm') !,'C','soiltemp_10cm')

   !# Register environmental dependencies
   data%id_E_rain = aed2_locate_global_sheet('rain')                  ! daily rainfall
   data%id_E_area = aed2_locate_global_sheet('layer_area')            ! cell area
   data%id_E_matz = aed2_locate_global_sheet('material')              ! material index
   data%id_E_bath = aed2_locate_global_sheet('bathy')                 ! cell bathy
   data%id_E_salt = aed2_locate_global('salinity')                    ! salinity of overlying water
   data%id_E_nearlevel= aed2_locate_global_sheet('nearest_depth')
   data%id_E_airtemp = aed2_locate_global_sheet('air_temp')
   data%id_E_I0 = aed2_locate_global_sheet('par_sf')

!   data%id_tem = aed2_locate_global('temperature')
!   data%id_dz = aed2_locate_global('layer_ht')

   !# NOTE: Initialisation occurs in first call of calculate_riparian

END SUBROUTINE aed2_define_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_initialize_vegetation(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to set initial state of ASS variables                                !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: veg_i
!-------------------------------------------------------------------------------
!BEGIN
   !---------------------------------------------------------------------------
!   ! Prime local cell hydrology object with data from global arrays
!   CALL SetSoilHydrology(data, column, SoilCol)

   !---------------------------------------------------------------------------

   DO veg_i=1,data%num_veg
     ! Retrieve current (local) sta te variable values
     _STATE_VAR_S_(data%id_veg(veg_i)) = 1e6 * _DIAG_VAR_S_(data%id_l_lai) / data%num_veg
   END DO


END SUBROUTINE aed2_initialize_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian_vegetation(data,column,layer_idx, pc_wet)
!-------------------------------------------------------------------------------
! Right hand sides of vegetation biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in)    :: data
   TYPE (aed2_column_t),          INTENT(inout) :: column(:)
   INTEGER,                       INTENT(in)    :: layer_idx
   AED_REAL,                      INTENT(in)    :: pc_wet
!
!LOCALS
   INTEGER  :: veg_i
   AED_REAL :: veg, pc_cover, lai !State variables
   AED_REAL :: veg_flux, photsynthesis, respiration, mortality, litterfall
   AED_REAL :: tveg, atemp, salinity, Io, theta, phreatic_depth, matz
   AED_REAL :: fI, fT, fSal, fW
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_E_matz)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Retrieve current environmental conditions
   salinity = _STATE_VAR_(data%id_E_salt)   ! local salinity
   atemp = _STATE_VAR_S_(data%id_E_airtemp)      ! local air temperature
   Io = _STATE_VAR_S_(data%id_E_I0)         ! surface short wave radiation
   theta = _STATE_VAR_S_(data%id_l_theta)+0.01   ! soil moisture in the unsaturated zone
   phreatic_depth = _STATE_VAR_S_(data%id_l_phreatic)   ! depth to water table

   ! Initialise cumulative biomass diagnostics
   _DIAG_VAR_S_(data%id_l_lai) = zero_
   _DIAG_VAR_S_(data%id_atm_co2) = zero_
   _DIAG_VAR_S_(data%id_gpp) = zero_
   _DIAG_VAR_S_(data%id_resp) = zero_
   _DIAG_VAR_S_(data%id_lfall) = zero_
   tveg = zero_

   DO veg_i=1,data%num_veg

      ! Retrieve current (local) state variable values
      veg = _STATE_VAR_S_(data%id_veg(veg_i)) ! vegetation group i
      pc_cover = _STATE_VAR_S_(data%id_vegfrac(veg_i)) ! vegetation cover i

      !# Photosynthesis
      !extc = _STATE_VAR_(data%id_extc)
      !dz   = _STATE_VAR_(data%id_dz)     ! dz = 0.5
      !fI   = photosynthesis_irradiance(data%mphydata(mphy_i)%lightModel, &
                    ! data%mphydata(mphy_i)%I_K, data%mphydata(mphy_i)%I_S, par, extc, Io, dz)
      ! Light
      fI   = 1.0 * (Io/2000.)
      IF( Io <100 ) fI = 0.0

      ! Temperature
      fT   = 1.03 ** (atemp-20.)

      ! Water availability
      IF(pc_wet<0.1) THEN
        ! Photosynthesis depends on vadose zone soil moisture availability
        fW   = MAX(zero_, (-4.3717 * theta * theta + 5.1364 * theta - 0.5455))
      ELSE
        ! Depending on species, innundation limites ability for photosynthesis
        fW   = 0.
        IF (data%vegdata(veg_i)%growthForm == 3 ) THEN !!?
          fW   = 1.
        ENDIF
      ENDIF

      photsynthesis = data%vegdata(veg_i)%R_growth * fI * fT * fW

      ! Salinity stress effect on respiration
      fSal = 1.0

      ! Respiration and general metabolic loss
      respiration = bio_respiration(data%vegdata(veg_i)%R_resp, data%vegdata(veg_i)%theta_resp, atemp) * fSal

      litterfall = bio_respiration(0.003, 1.10, atemp)

      IF( .NOT.data%simStaticBiomass ) THEN
        veg_flux = (photsynthesis - respiration - litterfall) * veg * pc_cover

        ! Set bottom fluxes for the pelagic (change per surface area per second)
        _FLUX_VAR_B_(data%id_veg(veg_i)) = _FLUX_VAR_B_(data%id_veg(veg_i)) + veg_flux
      ENDIF
      IF( data%simVegFeedback ) THEN
        !    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + mphy_flux
        !    _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - mphy_flux
        IF(pc_wet<0.1) THEN
          ! Drop literfall in soil litter
          _FLUX_VAR_B_(data%id_l_litter_land) = _FLUX_VAR_B_(data%id_l_litter_land) + litterfall*veg*pc_cover
        ELSE
          ! Drop literfall into water organic matter
          _FLUX_VAR_T_(data%id_l_litter_water) = _FLUX_VAR_T_(data%id_l_litter_water) + litterfall*veg*pc_cover
        ENDIF
      ENDIF

      ! Export diagnostic variables
      _DIAG_VAR_S_(data%id_atm_co2) = _DIAG_VAR_S_(data%id_atm_co2) + (respiration-photsynthesis)*veg*pc_cover
      _DIAG_VAR_S_(data%id_gpp) = _DIAG_VAR_S_(data%id_gpp) + photsynthesis*veg*pc_cover
      _DIAG_VAR_S_(data%id_resp) = _DIAG_VAR_S_(data%id_resp) + respiration*veg*pc_cover
      _DIAG_VAR_S_(data%id_lfall) = _DIAG_VAR_S_(data%id_lfall) + litterfall*veg*pc_cover
      tveg = tveg + veg
   ENDDO

   _DIAG_VAR_S_(data%id_tveg) = tveg
   _DIAG_VAR_S_(data%id_l_lai) = tveg


END SUBROUTINE aed2_calculate_riparian_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_light_shading_vegetation(data,column,layer_idx,shade_frac)
!-------------------------------------------------------------------------------
! Get the light shading fraction due to trees and/or emergent vegetation
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: shade_frac
!
!LOCALS
   AED_REAL :: veg, lai
   INTEGER  :: veg_i
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Retrieve current (local) state variable values.
   !DO veg_i=1,data%num_veg
   !   veg = veg + _STATE_VAR_S_(data%id_veg(veg_i))
   !END DO

   ! needs updating for overstorey vs understory

   lai = _DIAG_VAR_S_(data%id_l_lai)

   ! Self-shading with explicit contribution from background OM concentration.
   shade_frac = shade_frac + (1.-exp(-0.005*lai))

END SUBROUTINE aed2_light_shading_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE aed2_rain_loss_vegetation(data,column,layer_idx,infil)
!-------------------------------------------------------------------------------
! Get the soil moisture deficit, so host model can infiltrate rain accordingly
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!LOCALS
   AED_REAL :: smd

!-------------------------------------------------------------------------------
!BEGIN
   infil = data%int_max * ( _DIAG_VAR_S_(data%id_l_lai)/data%lai_max )

END SUBROUTINE aed2_rain_loss_vegetation




!###############################################################################
PURE AED_REAL FUNCTION vegetation_photosynthesis(data,rsa2,temp,kR2,yearday)
  !function x=p_a(rsa2,temp,kR2,yearday)
!-------------------------------------------------------------------------------
!%%%######          Running and Coughlan, 1998             #####
!% %%%######         Lohammar et al., 1980                  #####
!% %%%%%%%%%%%%%%%%    ASSIMILATION      %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!ARGUMENTS
CLASS (aed2_vegetation_data_t),INTENT(in) :: data  ! Module data, with params
AED_REAL, INTENT(IN) :: rsa2,temp,kR2,yearday
!
!LOCALS
AED_REAL , PARAMETER :: ext=0.5
AED_REAL , PARAMETER :: lai_proj=2.2   !2.2 is the coef that changes lai to projected lai in running 1988
AED_REAL , PARAMETER :: Q0=432.         !running 1988 kJ/m2d
AED_REAL , PARAMETER :: Q05=9730.       !running 1988 kJ/m2d
AED_REAL , PARAMETER :: tmax=37.        !running 1988 cent.
AED_REAL , PARAMETER :: tmin=0.         !running 1988
AED_REAL , PARAMETER :: CCmax=0.0016   !running 1988 m/s max canopy conductance of CO2
AED_REAL , PARAMETER :: Cair=7.10e-4   !Lohamar 1980 kg/m3  6.98
AED_REAL , PARAMETER :: Cleaf=0.        !Lohamar 1980 kg/m3
AED_REAL , PARAMETER :: lat=32.
AED_REAL :: drad, CMq, CMt, CM, CMmax, CC, dayl, psn

!relation betweem biomass and litterfall
! kl3u=10*10^-5;
! krt2u=10^-6;
! krt3u=2*10^-6;
!
!-------------------------------------------------------------------------------
!BEGIN
  drad=rsa2    !drad=canopy average daily radiation MJ/m2d

  !instead of using running attenuation, i'm using the canopy radiation
  !from Feikema et al, 2010:how much is absorved by canopy 2=overstory;
  !3=understory; but the rest is like running
  ! I'm multipling by 10^3 because the parameters
  !Q0 and Q0.5 are in kJ/m2d and rsa2 is in MJ/m2d!

  CMq=(drad*1e3-Q0)/(drad*1e3+Q05)

  if (CMq>1.) CMq=1.0
  if (CMq<0.) CMq=0.0

  CMt=((tmax-temp)*(temp-tmin))/(tmax*tmax)
  if (CMt>1.) CMt=1.0
  if (CMt<0.) CMt=0.0

  CM=CMq*CMt*CMmax
  if (CM>1.) CM=1.

  CC=CCmax*kR2
  if (CC>1.) CC=1.0

  dayl=((exp(7.42+0.045*lat))*(sin(yearday-79)*0.01721)+43200);

  PSN=((CC*CM)/(CC+CM))*dayl !potential Canopy photosynthesis [m]
  if ( (CC+CM)==0. ) PSN=0.
  !% if PSN<0;PSN=0;end;
  !% if PSN>0.02;PSN=0.02;end;
  vegetation_photosynthesis=PSN

END FUNCTION vegetation_photosynthesis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





END MODULE aed2_vegetation
