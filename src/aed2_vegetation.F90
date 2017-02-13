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
   USE aed2_util,ONLY : find_free_lun,aed2_bio_temp_function,fTemp_function,qsort

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
   END TYPE

   TYPE,extends(aed2_model_data_t) :: aed2_vegetation_data_t
      !# Variable identifiers
      INTEGER  :: id_veg(MAX_VEG_TYPES)
      INTEGER  :: id_vegfrac(MAX_VEG_TYPES)
      INTEGER  :: id_l_litter_land,id_l_litter_water

      !# Environmental variables
      INTEGER :: id_E_rain, id_E_area, id_E_matz, id_E_bath, id_E_salt, id_E_nearlevel

      !# Diagnostic variables
      INTEGER :: id_lai, id_gpp, id_lfall, id_resp

      !# Dependant variable IDs
      INTEGER :: id_l_depth, id_l_Sb, id_l_St, id_l_Ssat, id_l_theta, id_l_Stop, id_l_capz
      INTEGER :: id_l_phreatic, id_l_qss, id_l_qse, id_l_qcap, id_l_qper, id_l_wt
      INTEGER :: id_l_soiltemp


      !# Model parameters
      INTEGER  :: num_veg
      LOGICAL  :: simLitterLand, simLitterWater
      LOGICAL  :: simVegFeedback, simStaticBiomass, simVegFrac
      TYPE(macrophyte_data),DIMENSION(:),ALLOCATABLE :: vegdata

     CONTAINS
         PROCEDURE :: define             => aed2_define_vegetation
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_vegetation
!        PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_vegetation
!        PROCEDURE :: mobility           => aed2_mobility_vegetation
!        PROCEDURE :: light_extinction   => aed2_light_extinction_vegetation
!        PROCEDURE :: light_shading      => aed2_light_shading_vegetation
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
   CHARACTER(len=64)  :: litter_target_variable_land  = ''
   CHARACTER(len=64)  :: litter_target_variable_water = ''
   CHARACTER(len=128) :: dbase='aed2_vegetation_pars.nml'

   NAMELIST /aed2_vegetation/ num_veg, the_veg, dbase, &
                              litter_target_variable_land, litter_target_variable_water, &
                              simVegFeedback, simStaticBiomass, simVegFrac

!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed2_vegetation initialization"
   print *,"  WARNING! aed2_vegetation model is currently under development"

   ! Read the namelist
   read(namlst,nml=aed2_vegetation,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_vegetation'

    data%num_veg = 0

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
   data%id_lai   = aed2_define_sheet_diag_variable('lai','mmolC/m**2/d',  'vegetation LAI')
   data%id_gpp   = aed2_define_sheet_diag_variable('gpp','mmolC/m**2/d',  'net vegetation productivity')
   data%id_resp  = aed2_define_sheet_diag_variable('resp','mmolC/m**2/d','net vegetation respiration')
   data%id_lfall = aed2_define_sheet_diag_variable('lfall','mmolC/m**2/d','net vegetation mortality')

   ! Register module dependencies
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

   !# NOTE: Initialisation occurs in first call of calculate_riparian

END SUBROUTINE aed2_define_vegetation
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
   AED_REAL           :: veg,temp,salinity !State variables
   AED_REAL           :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for veg groups
   AED_REAL           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL           :: pon, pop, poc             !Mortaility and literfall state variables
   AED_REAL           :: FGrowth_Limitation, f_Temp, f_Salinity
   AED_REAL           :: photsynthesis, respiration, mortality !Growth & decay functions
   AED_REAL           :: pon_excr, pop_excr, poc_excr !POM excretion rates
   AED_REAL           :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   INTEGER            :: veg_i

   !CAB added
   AED_REAL :: f_dens, W, Imax, psuedofaeces, ingestion, excretion, egestion, iteg, R20, oxy

!
!-------------------------------------------------------------------------------
!BEGIN

RETURN
   ! Retrieve current environmental conditions.
   salinity = _STATE_VAR_(data%id_E_salt)   ! local salinity


END SUBROUTINE aed2_calculate_riparian_vegetation
!END SUBROUTINE aed2_calculate_benthic_vegetation
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
   AED_REAL :: veg
   INTEGER  :: veg_i
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Retrieve current (local) state variable values.
   DO veg_i=1,data%num_veg
     veg = veg + _STATE_VAR_S_(data%id_veg(veg_i))
   END DO

   ! Self-shading with explicit contribution from background OM concentration.
   shade_frac = shade_frac + veg*0.0001 !coefficient

END SUBROUTINE aed2_light_shading_vegetation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
FUNCTION aed2_vegetation_respiration(data,veg_i,iteg,temp,sal) RESULT(resp)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER  :: veg_i ! Invertebrate group
   AED_REAL, INTENT(IN)                :: temp  ! Temp value being used
   AED_REAL, INTENT(IN)                :: sal   ! Salinity value being used
   AED_REAL, INTENT(IN)                :: iteg  ! Ingestion-Egestion
!
!LOCALS
   AED_REAL :: W, TmaxR, maxTR, VV,WW,YY,XX,fT
   AED_REAL :: resp, Q, R20
!
!-------------------------------------------------------------------------------
!BEGIN

     !Make this better
     resp = data%vegdata(veg_i)%R_resp

!   ! Get the salinity limitation.
!   resp = resp * fSalinity_Limitation(data%vegdata,veg_i,sal)


END FUNCTION aed2_vegetation_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





!###############################################################################
PURE AED_REAL FUNCTION fTemp_function_veg(data,veg,temp)
!-------------------------------------------------------------------------------
! Temperature growth multiplier for vegdata
!
!ARGUMENTS
   CLASS (aed2_vegetation_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                 :: veg   ! Veg functional group
   AED_REAL, INTENT(IN)                :: temp  ! Temp value being used
!
!LOCALS
   AED_REAL , PARAMETER :: a = 1.00
   AED_REAL             :: MINt,Tmin,Tmax,MAXt

!
!-------------------------------------------------------------------------------
!BEGIN
       fTemp_function_veg = zero_

 END FUNCTION fTemp_function_veg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++










END MODULE aed2_vegetation
