!###############################################################################
!#                                                                             #
!# aed2_soilbgc.F90                                                                #
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
!# Created Feb 2017                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_soilbgc
!------------------------------------------------------------------------------+
! AED2 module for Riparian Soil Biogeochemistry                                |
!------------------------------------------------------------------------------+
   USE aed2_core
   USE aed2_util
   USE aed2_riptypes

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_soilbgc_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_soilbgc_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_pom0(:), id_pom(:)
      INTEGER,ALLOCATABLE :: id_eh(:), id_ch4(:)
      INTEGER :: id_uzdom, id_szdom, id_toc, id_pomt, id_litter
      INTEGER :: id_wettime, id_drytime
      INTEGER :: id_soilbgctracer

      !# Environmental variables
      INTEGER :: id_E_rain, id_E_area, id_E_matz, id_E_bath, id_E_salt, id_E_nearlevel, id_E_airtemp

      !# Diagnostic variables
      INTEGER :: id_rwet, id_rchg, id_bflw, id_omox, id_so4r, id_pml, id_sflw
      INTEGER :: id_atm_co2, id_atm_ch4, id_atm_n2o

      !# Dependant variable IDs
      INTEGER :: id_o_doc, id_c_dic
      INTEGER :: id_l_depth, id_l_Sb, id_l_St, id_l_Ssat, id_l_theta, id_l_Stop, id_l_capz
      INTEGER :: id_l_phreatic, id_l_qss, id_l_qse, id_l_qcap, id_l_qper, id_l_wt
      INTEGER :: id_l_soiltemp

      LOGICAL :: simProfiles

      !# ASS params
      INTEGER  :: n_zones, nlay
      AED_REAL,ALLOCATABLE :: active_zones(:)
      AED_REAL :: pom_0(MAX_ASS_PARAMS)
      AED_REAL :: Rom(MAX_ASS_PARAMS)
      AED_REAL :: flux_bf(MAX_ASS_PARAMS)
      AED_REAL :: flux_rn(MAX_ASS_PARAMS)
      AED_REAL :: flux_rn_max(MAX_ASS_PARAMS)
      AED_REAL :: flux_rw_a(MAX_ASS_PARAMS)
      AED_REAL :: flux_rw_d(MAX_ASS_PARAMS)
      AED_REAL :: zOM(MAX_ASS_PARAMS)
      AED_REAL :: Porosity(MAX_ASS_PARAMS),Density(MAX_ASS_PARAMS)
      AED_REAL :: X_hso4

     CONTAINS
         PROCEDURE :: define             => aed2_define_soilbgc
         PROCEDURE :: initialize         => aed2_initialize_soilbgc
!        PROCEDURE :: calculate_surface  => aed2_calculate_surface_soilbgc
!        PROCEDURE :: calculate          => aed2_calculate_soilbgc
!        PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_soilbgc
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_soilbgc
!        PROCEDURE :: calculate_dry      => aed2_calculate_dry_soilbgc
!        PROCEDURE :: equilibrate        => aed2_equilibrate_soilbgc
!        PROCEDURE :: mobility           => aed2_mobility_soilbgc
!        PROCEDURE :: light_extinction   => aed2_light_extinction_soilbgc
!        PROCEDURE :: delete             => aed2_delete_soilbgc
   END TYPE


!-------------------------------------------------------------------------------
!MODULE VARIABLES

   TYPE(SoilUnit) :: SoilCol
   INTEGER, PARAMETER :: OMCLAYL = 4
   INTEGER, PARAMETER :: OMCLAYH = 5
   INTEGER, PARAMETER :: OMSANDL = 6
   INTEGER, PARAMETER :: OMSANDH = 7
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_soilbgc(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_soilbgc_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER :: status, i

   !# ASS params
   INTEGER  :: n_zones, nlay
   INTEGER  :: active_zones(MAX_ZONES)
   AED_REAL :: pom_0(MAX_ASS_PARAMS)
   AED_REAL :: Rom(MAX_ASS_PARAMS)
   AED_REAL :: flux_bf(MAX_ASS_PARAMS)
   AED_REAL :: flux_rn(MAX_ASS_PARAMS)
   AED_REAL :: flux_rn_max(MAX_ASS_PARAMS)
   AED_REAL :: flux_rw_a(MAX_ASS_PARAMS)
   AED_REAL :: flux_rw_d(MAX_ASS_PARAMS)
   AED_REAL :: zOM(MAX_ASS_PARAMS)
   AED_REAL :: Porosity(MAX_ASS_PARAMS),Density(MAX_ASS_PARAMS)
   AED_REAL :: X_hso4

   LOGICAL  :: simProfiles

   CHARACTER(len=64) :: dom_link = ''
   CHARACTER(len=64) :: dic_link = ''
   CHARACTER(4) :: trac_name

   NAMELIST /aed2_soilbgc/ n_zones, active_zones, nlay, dom_link, simProfiles, &
                       pom_0, Rom, flux_bf, flux_rn, flux_rn_max, flux_rw_a, flux_rw_d, zOM, &
                       Porosity, Density, dic_link, X_hso4
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_soilbgc initialization"

   simProfiles = .FALSE.

   ! Read the namelist
   read(namlst,nml=aed2_soilbgc,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_soilbgc'

   data%simProfiles = simProfiles
   data%n_zones = n_zones         ! THIS NEEDS TO BE 4 AT THE MOMENT!
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   data%X_hso4 = X_hso4
   DO i=1,n_zones
      data%pom_0(i)  = pom_0(i)
      data%Rom(i)  = Rom(i)
      data%flux_bf(i)  = flux_bf(i)
      data%flux_rn(i) = flux_rn(i)
      data%flux_rn_max(i) = flux_rn_max(i)
      data%flux_rw_a(i) = flux_rw_a(i)
      data%flux_rw_d(i)  = flux_rw_d(i)
      data%zOM(i)  = zOM(i)
      data%nlay  = nlay
      data%Porosity(i) = Porosity(i)
      data%Density(i) = Density(i)
   ENDDO

   ALLOCATE(data%id_pom0(nlay));  ALLOCATE(data%id_eh(nlay))
   ALLOCATE(data%id_pom(nlay));  ALLOCATE(data%id_ch4(nlay))

   !# Register state variables
   data%id_litter = aed2_define_sheet_variable('litter','mmolC/m**2',          &
                                               'litter carbon density',        &
                                               zero_,                          &
                                               minimum=zero_)

   data%id_soilbgctracer = aed2_define_variable('tracer', 'mmol/m**3', 'dom tracer', &
                                                                zero_,minimum=zero_)
   IF ( dom_link .EQ. '' ) &
   data%id_o_doc = aed2_define_variable('dom_leachate', 'mmol/m**3', 'dom_leachate', &
                                                                zero_,minimum=zero_)
     ! Register diagnostic variables
   data%id_pomt  = aed2_define_sheet_diag_variable('pom_total','mol C/kg','soil total potential reactive carbon (depth averaged)')
   data%id_uzdom = aed2_define_sheet_diag_variable('uzdom','mol C/kg','unsat zone DOM')
   data%id_szdom = aed2_define_sheet_diag_variable('szdom','mol C/kg','sat zone DOM')
   data%id_toc   = aed2_define_sheet_diag_variable('toc','molH+/L','total actual acidity (porewater)')
   trac_name = 'lay0'
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_pom0(i) = aed2_define_sheet_diag_variable('pom0_'//TRIM(trac_name),'molH+/kg','starting particualte organic matter')
   ENDDO
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_pom(i)  = aed2_define_sheet_diag_variable('pom_'//TRIM(trac_name),'molH+/m2','particulate organic matter')
   ENDDO
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_eh(i) = aed2_define_sheet_diag_variable('anc0_'//TRIM(trac_name),'molH+/kg','starting acid neutralising capacity')
   ENDDO
   DO i=1,nlay
     trac_name(4:4) = CHAR(ICHAR('0') + i)
     data%id_ch4(i)  = aed2_define_sheet_diag_variable('anc_'//TRIM(trac_name),'molH+/m2','acid neutralising capacity')
   ENDDO
   data%id_rwet  = aed2_define_sheet_diag_variable('rwet','molH+/m2/day','rewetting flux of dom')
   data%id_rchg  = aed2_define_sheet_diag_variable('rchg','molH+/m2/day','percolation of dom')
   data%id_bflw  = aed2_define_sheet_diag_variable('bflw','molH+/m2/day','baseflow of dom')
   data%id_sflw  = aed2_define_sheet_diag_variable('sflw','molH+/m2/day','surface flow of dom')
   data%id_omox = aed2_define_sheet_diag_variable('omox','/day','OM oxidation rate')
   data%id_so4r  = aed2_define_sheet_diag_variable('so4r','molH+/m2/day','acid neutralising capacity')
   data%id_pml   = aed2_define_sheet_diag_variable('pml','m','past maximum groundwater level')
   data%id_wettime = aed2_define_sheet_diag_variable('wettime','day','time cell has been innundated')
   data%id_drytime = aed2_define_sheet_diag_variable('drytime','day','time cell has been exposed')
   data%id_atm_co2 = aed2_define_sheet_diag_variable('atm_co2_flux','mmolC/m**2/d',  'co2 exchange to the atmosphere')
   data%id_atm_ch4 = aed2_define_sheet_diag_variable('atm_ch4_flux','mmolC/m**2/d',  'ch2 exchange to the atmosphere')
   data%id_atm_n2o = aed2_define_sheet_diag_variable('atm_n2o_flux','mmolC/m**2/d',  'n2o exchange to the atmosphere')

   ! Register module dependencies
   IF ( .NOT. dom_link .EQ. '' ) &
   data%id_o_doc  = aed2_locate_global(TRIM(dom_link))
   IF ( .NOT. dic_link .EQ. '' ) &
   data%id_c_dic      = aed2_locate_global(TRIM(dic_link))
   data%id_l_depth    = aed2_locate_global_sheet('LND_depth')   !,'m','soil depth (to datum)')
   data%id_l_phreatic = aed2_locate_global_sheet('LND_phreatic')!,'m','depth of phreatic surface below surface')
   data%id_l_wt       = aed2_locate_global_sheet('LND_wt')      !,'m','depth of phreatic surface below surface')
   data%id_l_Sb       = aed2_locate_global_sheet('LND_Sb')      !,'mm','total capacity for water storage')
   data%id_l_St       = aed2_locate_global_sheet('LND_St')      ! 'mm', 'total soil water storage at time t')
   data%id_l_Ssat     = aed2_locate_global_sheet('LND_Ssat')    !,'mm','saturated zone soil water storage')
   data%id_l_Stop     = aed2_locate_global_sheet('LND_Stop')    !,'mm','unsaturated storage capacity')
   data%id_l_capz     = aed2_locate_global_sheet('LND_capz')    !,'mm','unsaturated storage capacity')
   data%id_l_theta    = aed2_locate_global_sheet('LND_theta')   !,'-','unsaturated moisture content')
   data%id_l_qss      = aed2_locate_global_sheet('LND_qss')     !,'mm','sat zone seepage')
   data%id_l_qse      = aed2_locate_global_sheet('LND_qs')     !,'mm','surface runoff')
   data%id_l_qcap     = aed2_locate_global_sheet('LND_qcap')    !,'mm','capillarity')
   data%id_l_qper     = aed2_locate_global_sheet('LND_qper')    !,'mm','recharge')
   data%id_l_soiltemp = aed2_locate_global_sheet('LND_soiltemp_10cm') !,'C','soiltemp_10cm')

   !# Register environmental dependencies
   data%id_E_rain = aed2_locate_global_sheet('rain')        ! daily rainfall
   data%id_E_area = aed2_locate_global_sheet('layer_area')  ! cell area
   data%id_E_matz = aed2_locate_global_sheet('material')    ! material index
   data%id_E_bath = aed2_locate_global_sheet('bathy')       ! cell bathy
   data%id_E_salt = aed2_locate_global('salinity')          ! salinity of overlying water
   data%id_E_nearlevel= aed2_locate_global_sheet('nearest_depth')
   data%id_E_airtemp =  aed2_locate_global_sheet('air_temp')

   ! Initialisation occurs in first call of calculate_riparian

END SUBROUTINE aed2_define_soilbgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_initialize_soilbgc(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to set initial state of soil BGC variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_soilbgc_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: POMt, ANCt
!-------------------------------------------------------------------------------
!BEGIN
   !---------------------------------------------------------------------------
   ! Prime local cell hydrology object with data from global arrays
   CALL SetSoilHydrology(data, column, SoilCol)

   !---------------------------------------------------------------------------
   ! Initialise vertical profiles of ANC and POM
   _DIAG_VAR_S_(data%id_pml) = SoilCol%PhreaticHgt  !Reset
   ANCt = 600. !_DIAG_VAR_S_(data%id_eh(1))
   POMt = _DIAG_VAR_S_(data%id_pom0(1))
   CALL SetSoilPOMProfile(data, column, SoilCol, POMt, ANCt)
   _DIAG_VAR_S_(data%id_uzdom)= zero_


END SUBROUTINE aed2_initialize_soilbgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian_soilbgc(data, column, layer_idx, pc_wet)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of soil biogeochemsitry and determine
! the flux to the water column from exposed or re-wetted soils
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_soilbgc_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   AED_REAL :: rain, salt, bathy, area, matz, soiltemp, atem
   INTEGER  :: var, i, sub, omz

   AED_REAL :: newDOM, domFlux, DOMmnlzn, avgomox, pom_vol
   AED_REAL :: Kass, sedDensity, pom_depth, moist, maxqse, OxdnRate,  NeutRate
   AED_REAL :: avgLevel, oldlevel, newlevel, SZDepth, UZDepth
   AED_REAL :: decomposition, litter

   AED_REAL :: dom_flux0, dom_flux1, dom_flux2, dom_flux3, dom_flux4, dom_flux5, dom_flux6, dom_flux7
   AED_REAL :: dom_gen, acid_potl, dom_minl, dom_anc

   AED_REAL :: DOMRWET, DOMRCHG, DOMBFLW, DOMSFLW, OMOX, SO4REDN
   AED_REAL :: POMt, TOC, SZDOM, UZDOM, SO4, DIC
!
!-------------------------------------------------------------------------------
!BEGIN

   !---------------------------------------------------------------------------
   ! Set local cell properties
   matz = _STATE_VAR_S_(data%id_E_matz)
   IF(.NOT.in_zone_set(matz,data%active_zones)) RETURN
   omz = 1
   DO i =1,data%n_zones
     IF( INT(matz) == data%active_zones(i) )THEN
       omz = i
     ENDIF
   ENDDO
   avgLevel = _STATE_VAR_S_(data%id_E_nearlevel)
   bathy = _STATE_VAR_S_(data%id_E_bath)
   IF( ABS(bathy) > 1e9 ) RETURN
   rain = _STATE_VAR_S_(data%id_E_rain)
   atem = _STATE_VAR_S_(data%id_E_airtemp)
   area = _STATE_VAR_S_(data%id_E_area)
   soiltemp = _STATE_VAR_S_(data%id_l_soiltemp)


   !---------------------------------------------------------------------------
   ! Zero dummy (working) variables
   dom_anc   = zero_
   dom_gen   = zero_
   acid_potl  = zero_
   dom_flux0 = zero_
   dom_flux1 = zero_
   dom_flux2 = zero_
   dom_flux3 = zero_
   dom_flux4 = zero_
   dom_flux5 = zero_
   dom_flux6 = zero_
   dom_flux7 = zero_
   dom_minl  = zero_
   pom_vol   = zero_
   avgomox   = zero_
   domFlux   = zero_
   newDOM    = zero_
   maxqse    = zero_

   ! Zero diagnostic variables
   OMOX    = zero_
   DOMRWET = zero_
   DOMRCHG = zero_
   DOMBFLW = zero_
   SO4REDN = zero_

   !---------------------------------------------------------------------------
   ! Get local OM properties for this cell
   TOC   = _DIAG_VAR_S_(data%id_toc)
   POMt  = _DIAG_VAR_S_(data%id_pom0(1))
   UZDOM = _DIAG_VAR_S_(data%id_uzdom)
   SZDOM = _DIAG_VAR_S_(data%id_szdom)

   !---------------------------------------------------------------------------
   ! Prime local cell hydrology object with data from global arrays
   CALL SetSoilHydrology(data, column, SoilCol)
   IF(SoilCol%qse>maxqse) THEN
      maxqse = SoilCol%qse
   END IF

   !---------------------------------------------------------------------------
   ! Get POM depth & soil volume relevant to the oxidisable layer
   pom_depth = MAX(SoilCol%Bathy - SoilCol%pastMaxLevel, zero_)
   pom_vol   = pom_depth * SoilCol%Area

   SZDepth    = SoilCol%Depth - SoilCol%PhreaticDepth
   UZDepth    = SoilCol%PhreaticDepth


   !---------------------------------------------------------------------------
   ! Dry sediment : Check for exposure of new POM to oxygen and do oxidation, etc
   IF(pc_wet<0.1) THEN

       ! Get WT level from SoilHydrology
       newLevel = SoilCol%PhreaticHgt

       !-----------------------------------------------------------------------!
       !-- Update surface litter
       litter = _STATE_VAR_S_(data%id_litter)
       decomposition = 0.003/86400 * 1.08**(atem-20.) * litter

       _FLUX_VAR_B_(data%id_litter) = _FLUX_VAR_B_(data%id_litter) - decomposition

       ! Partition decomposition into CO2 flux, DOM creation and addtion to soil POM
       _DIAG_VAR_S_(data%id_atm_co2) = 0.69*decomposition  ! assume mainly aerobic
       _DIAG_VAR_S_(data%id_atm_n2o) = 0.01*decomposition  ! assume mainly aerobic
       _DIAG_VAR_S_(data%id_pom(1)) = _DIAG_VAR_S_(data%id_pom(1)) + 0.25*decomposition
       UZDOM =  UZDOM + 0.05*decomposition

       !-----------------------------------------------------------------------!
       !-- Update redox profile
       CALL UpdateRedoxProfile(SoilCol,column)

       !-----------------------------------------------------------------------!
       !-- Change in POM profile and partitioning between anaerobic and aerobic
       !   due to change water table and reactions
       CALL UpdatePOMProfile( SoilCol,           &
                               column,           &
                               newDOM,           &
                               avgomox,          &
                               data%Rom(omz),    &
                               temp=soiltemp     )

       POMt = zero_
       DO i=1,data%nlay
         POMt = POMt + _DIAG_VAR_S_(data%id_pom(i)) + _DIAG_VAR_S_(data%id_pom0(i))
       ENDDO
       IF( pom_depth > 0.05 ) THEN
         _DIAG_VAR_S_(data%id_pomt) = POMt / (pom_depth * SoilCol%Density * 1e-3)
       ENDIF

       ! Increase available acidity
       UZDOM =  UZDOM + newDOM
       dom_gen = dom_gen + newDOM
       OMOX = avgomox

       newDOM = zero_ !consumption/mineralisation here
       !CALL UpdateRedoxProfile( SoilCol, &
       !                         column, &
       !                      newAcidity )
       !
       !ANCt = zero_
       !DO i=1,data%nlay
       !   ANCt = ANCt + _DIAG_VAR_S_(data%id_ch4(i))
       !ENDDO

       ! Consume available acidity
       UZDOM = UZDOM - newDOM
       dom_anc = dom_anc - newDOM

       ! Reset pastMaxLevel if necessary
       oldLevel = newLevel
       IF(SoilCol%PhreaticHgt < SoilCol%pastMaxLevel) THEN
           ! Newly exposed sediment is yet to contibute
           SoilCol%pastMaxLevel = SoilCol%PhreaticHgt
       END IF

       !-----------------------------------------------------------------------
       !-- Decrease in available DOM in SZ due to (anaerobic) mineralisation

       !newANC = data%flux_rn_max(assz) * MIN(SZDepth,0.5) * area * DDT * 1e-3
       DOMmnlzn = data%flux_rn_max(omz) * MIN(SZDepth,0.5) * DDT * 1e-3
       IF(DOMmnlzn > SZDOM) DOMmnlzn=SZDOM

       ! Consume available acidty
       SZDOM = SZDOM - DOMmnlzn

       dom_minl = dom_minl - DOMmnlzn

       !-----------------------------------------------------------------------
       !-- Decrease in available DOM in UZ due to recharge/discharge/etc

       ! R-RCG : recharge/percolation   [UNITS=> mol/L/timestep = mol/L * ((m/ts)/m) ]
       IF( SoilCol%recharge>zero_ .AND. UZDepth>0.1 ) THEN

         domFlux = UZDOM * 0.5 * ( MIN(SoilCol%recharge/UZDepth,1.0) **data%flux_rn(omz) )

         UZDOM = UZDOM - domFlux
         SZDOM = SZDOM + domFlux

         dom_flux0 = dom_flux0 + domFlux
         DOMRCHG = domFlux/DDT  ! molH+/L/day
       END IF

       ! R-BFLW : baseflow  [ UNITS>= mol/timestep = mol/L * m/day /m *day/ts ]
       IF( SoilCol%qss>zero_ .AND. SZDepth>0.01 ) THEN

        IF(SZDOM>zero_)THEN
           domFlux = data%flux_bf(omz) * SZDOM * ( SoilCol%qss / MIN(SZDepth,0.5) )
        ELSE
           domFlux = zero_
        ENDIF

         IF ( domFlux > SZDOM ) THEN
           domFlux = SZDOM
         END IF

         SZDOM = SZDOM - domFlux
         dom_flux1 = dom_flux1 + domFlux

        ! set diagnostic flux and update seepage tracer into the water: mol/m2/day
         DOMBFLW = domFlux/DDT
        _FLUX_VAR_R_(data%id_soilbgctracer) = _FLUX_VAR_R_(data%id_soilbgctracer) + DOMBFLW/secs_per_day


        ! _FLUX_VAR_R_(data%id_g_ubalchg) = _FLUX_VAR_R_(data%id_g_ubalchg) + 1e-3*ASSBFLW/secs_per_day
        _FLUX_VAR_R_(data%id_o_doc) = _FLUX_VAR_R_(data%id_o_doc) + DOMBFLW/secs_per_day
        IF ( data%id_c_dic>0 ) THEN
         _FLUX_VAR_R_(data%id_c_dic) = _FLUX_VAR_R_(data%id_c_dic) + (DOMBFLW/secs_per_day) !* data%X_hso4
        ENDIF

       END IF


       ! R-DIS : If water table moves up, move some UZDOM -> SZDOM
       IF(newLevel > oldLevel .AND. UZDepth > 0.01 ) THEN

         domFlux = MAX((( newLevel - oldLevel ) / UZDepth),1.0) * UZDOM

         UZDOM = UZDOM - domFlux
         SZDOM = SZDOM + domFlux

         dom_flux6 = dom_flux6 + domFlux
       END IF


       ! R-POND : ponding and sat excess  [ UNITS>= mol/timestep = mol/L * m/day /m *day/ts ]
       IF( SoilCol%qse>zero_ .AND. UZDOM>zero_ .AND. UZDepth>0.1 ) THEN

         !  [UNITS=> mol/L/timestep = mol/L *m/ts *day/tstep *m2 * ??]
         domFlux = UZDOM * ( MIN(SoilCol%qse/UZDepth,1.0) **data%flux_rn(omz) )

         UZDOM = UZDOM - domFlux
         dom_flux5 = dom_flux5 + domFlux

        ! set diagnostic flux and update seepage tracer into the water: mol/m2/day
         DOMSFLW = domFlux/DDT

        _FLUX_VAR_R_(data%id_soilbgctracer) = _FLUX_VAR_R_(data%id_soilbgctracer) + DOMSFLW/secs_per_day

        ! ! Store acidity flux in DissFluxRates  [ UNITS = mol /m2 /day ]
        _FLUX_VAR_R_(data%id_o_doc) = _FLUX_VAR_R_(data%id_o_doc) + DOMSFLW/secs_per_day
        ! _FLUX_VAR_R_(data%id_g_ubalchg) = _FLUX_VAR_R_(data%id_g_ubalchg) + 1e-3*ASSSFLW/secs_per_day
        IF ( data%id_c_dic>0 ) THEN
         ! SO4 flux associated with pyrite oxidation (SHOULD THIS BE MULTIPLIED BY 1000 as SO4 in mmol/m3?)
         _FLUX_VAR_R_(data%id_c_dic) = _FLUX_VAR_R_(data%id_c_dic) + (DOMSFLW/secs_per_day) !* data%X_hso4
        ENDIF

       END IF




       DOMRWET = 0.0
       _DIAG_VAR_S_(data%id_drytime) = _DIAG_VAR_S_(data%id_drytime) + DDT
       _DIAG_VAR_S_(data%id_wettime) = zero_

       TOC = zero_
       TOC = UZDOM + SZDOM + POMt !? * 1e3 / ( MAX(UZDepth,0.01)  * SoilCol%Density )

     !-------------------------------------------------------------------------
     ! R-WET : Newly re-wetted sediment
     ! Flux rate is initially high
     ELSE IF(pc_wet>0.1 .AND. _DIAG_VAR_S_(data%id_wettime)<=0.25) THEN

       salt = _STATE_VAR_(data%id_E_salt)
       domFlux = GetDOMFluxRate(data%flux_rw_a(omz),salt,omz,1)

       DOMRWET = domFlux

       _FLUX_VAR_B_(data%id_o_doc) = _FLUX_VAR_B_(data%id_o_doc) + DOMRWET/secs_per_day
       IF ( data%id_c_dic>0 ) THEN
         _FLUX_VAR_B_(data%id_c_dic) = _FLUX_VAR_B_(data%id_c_dic) + (DOMRWET/secs_per_day) !* data%X_hso4
       ENDIF

       dom_flux3 = dom_flux3 + domFlux*DDT

       _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT


       SZDOM = SZDOM + UZDOM
       UZDOM = zero_


     !-------------------------------------------------------------------------
     ! R-WET : Old re-wetted sediment - keep fluxing remainder (via diffusion)
     ELSE IF(pc_wet>0.1 .AND. _DIAG_VAR_S_(data%id_wettime)>0.25 ) THEN

       DIC = _STATE_VAR_(data%id_c_dic)
       salt = _STATE_VAR_(data%id_E_salt)
       IF(salt < 1.)  salt = 1.0
       IF(salt > 40.) salt = 40.0

       IF( _DIAG_VAR_S_(data%id_wettime)>90. .OR.  _DIAG_VAR_S_(data%id_drytime)<5.0) THEN
         ! Old innundation .or. soil not dry for very long before innundation

         ! SO4 reduction of 5mmol H+/day from Koschorreck M, Tittel J.
         domFlux = -0.005

         ! Scale SO4 reduction down if SO4 is limiting (>280mg/L SO4 is non-limiting)
         domFlux = domFlux * ( SO4 ) / ( SO4+(153.*(1e3/96.)) )

         SO4REDN = domFlux
         dom_flux7 = dom_flux7 + domFlux*DDT

         ! Past dryness and PASS production is reset by now
         _DIAG_VAR_S_(data%id_drytime) = zero_

       ELSEIF(_DIAG_VAR_S_(data%id_wettime) < 1.) THEN
         ! Innundated within last 1 days
         domFlux = GetDOMFluxRate(data%flux_rw_a(omz),salt,omz,1)

         dom_flux3 = dom_flux3 + domFlux*DDT

       ELSE
         ! Innundated within last 1-90 days
         domFlux = GetDOMFluxRate(data%flux_rw_d(omz),salt,omz,2)

         dom_flux4 = dom_flux4 + domFlux*DDT
       END IF

       ! Store acidity flux in ASSRWET for update within routine
       DOMRWET = domFlux * 1e-3

       _FLUX_VAR_B_(data%id_o_doc) = _FLUX_VAR_B_(data%id_o_doc) + DOMRWET/secs_per_day
       IF ( data%id_c_dic>0 ) THEN
         _FLUX_VAR_B_(data%id_c_dic) = _FLUX_VAR_B_(data%id_c_dic) + (DOMRWET/secs_per_day) * data%X_hso4
       ENDIF

       ! Update the time-counter for inundation
       _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT

     END IF
     !------------------------------


   ! Update the main arrays
   _DIAG_VAR_S_(data%id_uzdom) = UZDOM
   _DIAG_VAR_S_(data%id_szdom) = SZDOM
   _DIAG_VAR_S_(data%id_toc)    = TOC
   _DIAG_VAR_S_(data%id_rwet)   = DOMRWET
   _DIAG_VAR_S_(data%id_rchg)   = DOMRCHG
   _DIAG_VAR_S_(data%id_bflw)   = DOMBFLW
   _DIAG_VAR_S_(data%id_sflw)   = DOMSFLW
   _DIAG_VAR_S_(data%id_omox)   = OMOX
   _DIAG_VAR_S_(data%id_so4r)   = SO4REDN
   _DIAG_VAR_S_(data%id_pml)    = SoilCol%pastMaxLevel  ! _DIAG_VAR_S_(data%id_l_wt)



CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 SUBROUTINE UpdatePOMProfile(theSoil, column, newDOM, avgomox, Rom, temp)
   !-- Incoming
   TYPE(SoilUnit)                     :: theSoil
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   AED_REAL                           :: newDOM, avgomox, Rom
   AED_REAL, OPTIONAL                 :: temp
   !-- Local
   INTEGER  :: lay
   AED_REAL :: dep, depm1, middep, newAcidity, OxdnRate

   !---------------------------------------------------------------------------!

   newDOM = zero_
   avgomox = zero_

   DO lay = 1,theSoil%nlay

       newAcidity = zero_

!      dep   = theSoil%Bathy  + ( lay    * ( theSoil%Depth / theSoil%nlay ) )
!      depm1 = theSoil%Bathy  + ((lay-1) * ( theSoil%Depth / theSoil%nlay ) )
       dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
       depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )

       middep = (dep+depm1)/2.

       !-----------------------------------------------------------------------!
       ! Check for newly exposed layers and update POM_ox based on POM_0
       IF(middep > theSoil%pastMaxLevel .AND. _DIAG_VAR_S_(data%id_pom0(lay)) > zero_ ) THEN

         _DIAG_VAR_S_(data%id_pom(lay)) = (dep - depm1)  &
                      * _DIAG_VAR_S_(data%id_pom0(lay))  * theSoil%Density * 1e-3
         ! Once a layer is added to the POM_ox array it cannot be repeated.
         _DIAG_VAR_S_(data%id_pom0(lay)) = zero_

         IF( middep > 0.5 )  _DIAG_VAR_S_(data%id_pom(lay)) = zero_
       END IF


       !-----------------------------------------------------------------------!
       !-- Increase in actual acidity due to oxidation process
       OxdnRate = GetOMOxidnRate( Rom, theSoil%Moisture(lay),  theSoil%Substrate )
       IF( PRESENT(temp) )  OxdnRate = OxdnRate * 1.05**(temp-20.)

       newAcidity = _DIAG_VAR_S_(data%id_pom(lay)) * OxdnRate * DDT

       ! Reduce POM0 based on used amount
       _DIAG_VAR_S_(data%id_pom(lay)) = _DIAG_VAR_S_(data%id_pom(lay)) - newAcidity

       ! Increase available acidity
    !   newAASS = newAASS + newAcidity

       avgomox  =  avgomox + OxdnRate
   END DO

   ! Profile averaged OM oxidation rate (for plots)
   avgomox =  avgomox/theSoil%nlay

 END SUBROUTINE UpdatePOMProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 SUBROUTINE UpdateRedoxProfile(theSoil,column)
   !-- Incoming
   TYPE(SoilUnit)       :: theSoil
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   !-- Local
   INTEGER :: lay
   AED_REAL :: dep, depm1, middep, dEh, wfps
   AED_REAL, PARAMETER :: CR = 100.    !100mv/day (Zhang et al)
   AED_REAL, PARAMETER :: aeren = 0.1

   !---------------------------------------------------------------------------!

   DO lay = 1,theSoil%nlay

     dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
     depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )
     middep = (dep+depm1)/2.

     !IF(middep > theSoil%pastMaxLevel)
     IF(middep > theSoil%PhreaticHgt)THEN !~?
         ! decrease Eh for soil below the water table
         dEh = CR * ( aeren-1. )
     ELSE
         wfps = theSoil%Moisture(lay)/theSoil%Porosity
         dEh = CR * ( aeren+1.-wfps )
     ENDIF

     _DIAG_VAR_S_(data%id_eh(lay)) = _DIAG_VAR_S_(data%id_eh(lay)) + dEh*DDT
     IF(_DIAG_VAR_S_(data%id_eh(lay)) < -300) _DIAG_VAR_S_(data%id_eh(lay)) = -300
     IF(_DIAG_VAR_S_(data%id_eh(lay)) > 600) _DIAG_VAR_S_(data%id_eh(lay)) = 600

   END DO

 END SUBROUTINE UpdateRedoxProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 SUBROUTINE UpdateMethaneProfile(theSoil,column,temp)
   !-- Incoming
   TYPE(SoilUnit)       :: theSoil
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   AED_REAL :: temp
   !-- Local
   INTEGER :: lay
   AED_REAL :: dM,Mprd,Moxd,Mdfs,Mebl,Mplt, ch4
   AED_REAL :: dep, depm1, middep, fEhMP, fEhMO, fT
   AED_REAL, PARAMETER :: CM = 1. !???
   AED_REAL, PARAMETER :: Kch4 = 1. !???

   !---------------------------------------------------------------------------!

   DO lay = 1,theSoil%nlay

     dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
     depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )
     middep = (dep+depm1)/2.

     Mprd = zero_
     Moxd = zero_
     Mdfs = zero_
     Mebl = zero_
     Mplt = zero_

     ch4 = _DIAG_VAR_S_(data%id_ch4(lay))

     ! Production
     fEhMP = zero_
     IF( _DIAG_VAR_S_(data%id_eh(lay)) < -100. ) THEN
       fEhMP = (_DIAG_VAR_S_(data%id_eh(lay))- -100.)*-0.01
       IF( fEhMP>one_ ) fEhMP=one_
     ENDIF
     fT = 1.05**(temp-20.)
     Mprd = CM * fT * fEhMP

     ! Oxidation
     fEhMO = zero_
     IF( _DIAG_VAR_S_(data%id_eh(lay)) > -200. ) THEN
       fEhMO = (_DIAG_VAR_S_(data%id_eh(lay))+ 200.)*0.0025
       IF( fEhMO>one_ ) fEhMO=one_
     ENDIF
     fT = 1.05**(temp-20.)
     Moxd = (ch4/(Kch4+ch4)) * fT * fEhMO

     ! Diffusion


     ! Ebullition


     dM = Mprd-Moxd-Mdfs-Mebl-Mplt

     _DIAG_VAR_S_(data%id_ch4(lay)) = _DIAG_VAR_S_(data%id_ch4(lay)) + dM*DDT
     IF(_DIAG_VAR_S_(data%id_ch4(lay)) < zero_) _DIAG_VAR_S_(data%id_ch4(lay)) = zero_
     IF(_DIAG_VAR_S_(data%id_ch4(lay)) > 1e6) _DIAG_VAR_S_(data%id_ch4(lay)) = 1e6

   END DO

 END SUBROUTINE UpdateMethaneProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 SUBROUTINE UpdateANCProfile(theSoil,column,conAASS)
   !-- Incoming
   TYPE(SoilUnit)       :: theSoil
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   !AED_REAL, DIMENSION(20)   :: theANC
   !AED_REAL, DIMENSION(20)   :: ANC0
   AED_REAL                  :: conAASS
   !-- Local
   INTEGER :: lay
   AED_REAL :: dep, depm1, middep, newAcidity, NeutRate, theTAA

   !---------------------------------------------------------------------------!

   theTAA = conAASS
   conAASS = zero_

   DO lay = 1,theSoil%nlay

       newAcidity = zero_

      !dep   = theSoil%Bathy + ( lay * ( theSoil%Depth/theSoil%nlay ))
      !depm1 = theSoil%Bathy + ((lay-1) * ( theSoil%Depth/theSoil%nlay ))
       dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
       depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )

       middep = (dep+depm1)/2.

       !-----------------------------------------------------------------------!
       ! Check for newly exposed layers and update ANCcontent
       IF(middep > theSoil%pastMaxLevel .AND. _DIAG_VAR_S_(data%id_eh(lay)) > zero_ ) THEN

         _DIAG_VAR_S_(data%id_ch4(lay)) = (dep - depm1) & !* theSoil%Area &
                      * _DIAG_VAR_S_(data%id_eh(lay))  * theSoil%Density * 1e-3
         ! Once a layer is added to the ANCcontent array it cannot be repeated.
         _DIAG_VAR_S_(data%id_eh(lay)) = zero_

       END IF

       !-----------------------------------------------------------------------!
       !-- Decrease in actual acidity due to nuetralisation process
       NeutRate = 0.00027  * 2.0   !20% /year

       IF(theTAA*1e3 > 1.0) THEN
         ! [ UNITS = mol/L * /day * day/timestep = mol /L /timestep ]
         newAcidity = _DIAG_VAR_S_(data%id_ch4(lay)) * NeutRate * DDT
       ELSE
         newAcidity = zero_
       END IF

       ! Reduce ANC based on used amount
        _DIAG_VAR_S_(data%id_ch4(lay)) =  _DIAG_VAR_S_(data%id_ch4(lay)) - newAcidity

       ! Increase consumed acidty
       conAASS = conAASS + newAcidity


   END DO


 END SUBROUTINE UpdateANCProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END SUBROUTINE aed2_calculate_riparian_soilbgc


!###############################################################################
SUBROUTINE SetSoilHydrology(data, column, theSoil)
!-------------------------------------------------------------------------------
! Primes the local cell Bucket model "object" with relevant parameters based on!
! its material zone id number                                                  !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_soilbgc_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   TYPE(SoilUnit),INTENT(inout) :: theSoil
   !INTEGER,INTENT(in) :: MatZoneID
!LOCALS
   INTEGER :: lay, zindex, i
   AED_REAL :: dep, depm1, middep, botmoist, topmoist

!-------------------------------------------------------------------------------
!BEGIN

   theSoil%UnitNum = -1
   theSoil%Substrate = INT(_STATE_VAR_S_(data%id_E_matz))
   zindex = 1
   DO i =1,data%n_zones
     IF( theSoil%Substrate == INT(data%active_zones(i)) )THEN
       zindex = i
     ENDIF
   ENDDO

   ! Physical properties
   theSoil%Depth = _DIAG_VAR_S_(data%id_l_depth)
   theSoil%Area = _STATE_VAR_S_(data%id_E_area)
   theSoil%Bathy = _STATE_VAR_S_(data%id_E_bath)
   theSoil%Sb = _STATE_VAR_S_(data%id_l_Sb)
   theSoil%Porosity = data%Porosity(zindex)
   theSoil%Density = data%Density(zindex)

   ! Hydrologic dynamics (these have been previously sent by aed2_land)
   theSoil%St = _DIAG_VAR_S_(data%id_l_St)
   theSoil%Ssat = _DIAG_VAR_S_(data%id_l_Ssat)
   theSoil%Sus = MAX( theSoil%St-theSoil%Ssat , zero_ )
   theSoil%S_top = MAX( _DIAG_VAR_S_(data%id_l_Stop) , zero_ )
   theSoil%theta = _DIAG_VAR_S_(data%id_l_theta)
   !theSoil%Ucap = _DIAG_VAR_S_(data%id_l_Ucap)

   theSoil%PhreaticDepth =  _DIAG_VAR_S_(data%id_l_phreatic)
   theSoil%PhreaticHgt = _DIAG_VAR_S_(data%id_l_wt) !theSoil%Bathy - theSoil%PhreaticDepth
   theSoil%CapHgt = _DIAG_VAR_S_(data%id_l_capz)
   theSoil%Dcap = 0.5*  (theSoil%CapHgt - theSoil%PhreaticHgt)
   theSoil%Dtrn = 0.5*  (theSoil%CapHgt - theSoil%PhreaticHgt)
   theSoil%pastMaxLevel = _DIAG_VAR_S_(data%id_pml)

   theSoil%qss = _DIAG_VAR_S_(data%id_l_qss)
   theSoil%qse = _DIAG_VAR_S_(data%id_l_qse)
   theSoil%recharge = _DIAG_VAR_S_(data%id_l_qper)
   theSoil%qsuc = _DIAG_VAR_S_(data%id_l_qcap)

   ! Moisture disaggregation model
   theSoil%nlay = data%nlay

   IF(.NOT.ALLOCATED(theSoil%Moisture)) ALLOCATE( theSoil%Moisture(theSoil%nlay) )

   theSoil%Moisture(:) = theSoil%theta


   ! Update Layer Moisture Fractions (wt% water)
   botmoist = 0.9
   topmoist = (theSoil%S_top) /  MAX( (theSoil%Bathy-theSoil%CapHgt) * theSoil%Porosity,0.001)

   DO lay = 1,theSoil%nlay

         dep   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
         depm1 = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )

         middep = (dep+depm1)/2.

         theSoil%Moisture(lay) = 0.0

         IF( middep  <= theSoil%PhreaticHgt ) THEN
            ! Saturated Zone: % VWC/theta
            theSoil%Moisture(lay) = 1.000
         ELSE IF(middep > theSoil%PhreaticHgt .AND. middep <=  theSoil%PhreaticHgt+theSoil%Dcap) THEN
            ! Near saturated region above the water table: % VWC/theta
            theSoil%Moisture(lay) = botmoist
         ELSE IF(middep >  theSoil%PhreaticHgt+theSoil%Dcap .AND. middep <= theSoil%CapHgt) THEN
            ! Transition region: % VWC/theta
            !print *,'middep',topmoist,theSoil%Ztrn, (botmoist-topmoist) , ((middep)-theSoil%Bathy),(theSoil%Dcap - middep)
            theSoil%Moisture(lay) = botmoist - (botmoist-topmoist) *  (middep-(theSoil%PhreaticHgt+theSoil%Dcap)) / theSoil%Dtrn  ! linear gradient
         ELSE
            ! Top of unsaturated zone: % VWC/theta
            theSoil%Moisture(lay) = topmoist !&
            ! MIN(theSoil%S_top / MAX(theSoil%Sb - theSoil%Ssat - theSoil%S_trn - theSoil%S_cap,1e-2),1.0)
         ENDIF

         !print *,'moist',lay,dep,middep,theSoil%Moisture(lay)

         ! convert to % wtWC
         theSoil%Moisture(lay) = theSoil%Moisture(lay)*theSoil%Porosity * 1e3 /        &
            ( theSoil%Moisture(lay)*theSoil%Porosity*1e3 + (1.-theSoil%Porosity)*theSoil%Density )

         !theSoil%Moisture(lay) =theSoil%theta  !TEMPP
   ENDDO

END SUBROUTINE SetSoilHydrology
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 FUNCTION GetOMOxidnRate(MaxRate, Moisture, SOIL) RESULT (OxdnRate)
   !-- Incoming
   AED_REAL :: MaxRate,Moisture   ! Max rate of FeS2 oxidation for this soil
   AED_REAL :: OxdnRate           ! The actual oxidation rate calculated
   INTEGER  :: SOIL               ! Integer soil type
   !-- Local
   AED_REAL,PARAMETER   :: X=0.0333 ! Parameter - Jeff Turner ASS study
   AED_REAL             :: MC

   !---------------------------------------------------------------------------!

   MC = Moisture

   OxdnRate = 0.0

   ! SANDY ASS
   IF(SOIL == OMSANDL .OR. SOIL == OMSANDH) THEN

     OxdnRate = -9.7011*MC**3. + 2.1949*MC**2. + 0.0025*MC + 0.0006

     OxdnRate = MaxRate* OxdnRate

   ! CLAYEY ASS
   ELSE IF(SOIL == OMCLAYL .OR. SOIL == OMCLAYH) THEN

     IF(MC > 0.225 .AND. MC < 0.48) THEN
       OxdnRate = -0.0142*MC + 0.0068
     ELSE IF(MC <= 0.225) THEN
       OxdnRate = 0.0142*MC
     ELSE
       OxdnRate = 0.0
     END IF

     OxdnRate = MaxRate* OxdnRate

   ! UNSURE
   ELSE

     OxdnRate = ( MaxRate - X*(Moisture*0.5) )/2.0

   END IF

   ! CHECK
   IF(OxdnRate<0.0) OxdnRate = 0.0
   IF(OxdnRate>0.5) OxdnRate = 0.5

 END FUNCTION GetOMOxidnRate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
 FUNCTION GetDOMFluxRate(FWRate,Salinity,SOIL,theTime) RESULT (FluxRate)
   !-- Incoming
   AED_REAL   :: FWRate,Salinity    ! Rate of FW release of H for this soil
   AED_REAL   :: FluxRate           ! The actual flux rate calculated
   INTEGER    :: SOIL, theTime      ! Integer soil type & Time
   !-- Local
   AED_REAL,PARAMETER   :: BC1=0.436 ! Parameter - CSIRO ASS study CLAY
   AED_REAL,PARAMETER   :: BS1=0.152 ! Parameter - CSIRO ASS study SAND
   AED_REAL,PARAMETER   :: BCt=0.033 ! Parameter - CSIRO ASS study CLAY
   AED_REAL,PARAMETER   :: BSt=0.006 ! Parameter - CSIRO ASS study SAND
   AED_REAL             :: kk

   !---------------------------------------------------------------------------!
   IF(theTime == 1) THEN
     IF(SOIL == OMCLAYL .OR. SOIL == OMCLAYH) THEN
        kk = BC1
     ELSE
        kk = BS1
     END IF
   ELSE
     IF(SOIL == OMCLAYL .OR. SOIL == OMCLAYH) THEN
        kk = BCt
     ELSE
        kk = BSt
     END IF
   END IF

   IF(FWRate < 1e-10) THEN
     FluxRate = 0.00
   ELSE
     FluxRate =  kk * Salinity/35. + FWRate
   END IF

 END FUNCTION GetDOMFluxRate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE SetSoilPOMProfile(data, column, theSoil, PASSt, ANCt)
!-------------------------------------------------------------------------------
! Primes the local cell Bucket model "object" with relevant parameters based on!
! its material zone id number                                                  !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_soilbgc_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   TYPE(SoilUnit),INTENT(inout) :: theSoil ! passing SoilUnit into theSoil
   AED_REAL,INTENT(in) :: PASSt, ANCt
!LOCALS
   INTEGER :: lay, matz, zindex, i
   AED_REAL :: dep(theSoil%nlay), depm1(theSoil%nlay), middep(theSoil%nlay)
!-------------------------------------------------------------------------------
!BEGIN

   matz = INT(_STATE_VAR_S_(data%id_E_matz))
   IF (.NOT.in_zone_set(REAL(matz),data%active_zones) ) RETURN
   zindex = 1
   DO i =1,data%n_zones
     IF( matz == INT(data%active_zones(i)) )THEN
       zindex = i
     ENDIF
   ENDDO


   ! Set the grid depths of layers below surface
   DO lay = 1,theSoil%nlay
     dep(lay)   =  ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
     depm1(lay) =  (  lay *    ( theSoil%Depth / theSoil%nlay ) )
!     dep(lay)   = theSoil%Bathy  - ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
!     depm1(lay) = theSoil%Bathy  - (  lay *    ( theSoil%Depth / theSoil%nlay ) )
!     dep(lay)    = theSoil%Bathy  + (  lay *    ( theSoil%Depth / theSoil%nlay ) )
!     depm1(lay)  = theSoil%Bathy  + ( (lay-1) * ( theSoil%Depth / theSoil%nlay ) )
     middep(lay) = ( dep(lay)+depm1(lay) )/2.
   END DO

   ! Set active arrays pass and anc to 0 (no oxidation as yet) & pass0 and anc0 to top value
   DO lay = 1,theSoil%nlay
     _DIAG_VAR_S_(data%id_pom(lay)) = zero_
     _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt
     _DIAG_VAR_S_(data%id_ch4(lay)) = zero_
     _DIAG_VAR_S_(data%id_eh(lay)) = ANCt
   ENDDO
   _DIAG_VAR_S_(data%id_pomt) = zero_



   IF (.NOT. data%simProfiles)  RETURN


   ! Lake Alex == 1
   IF( data%zOM(zindex)>0.5 .AND. data%zOM(zindex)<1.5) THEN

       IF(theSoil%Substrate == OMCLAYL .OR. theSoil%Substrate == OMCLAYH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt       ! Surface concs
              _DIAG_VAR_S_(data%id_eh(lay))  = ANCt        ! Surface concs
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 1.009705089
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.187140007
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.83412953
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.487488851
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.667303624
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.761085667
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.500477718
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.034682483
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.333651812
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.308279299
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.166825906
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.581876115
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.215260078
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.855472932
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.205900944
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.01769224
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           END IF
         END DO

       ELSE IF(theSoil%Substrate == OMSANDL .OR. theSoil%Substrate == OMSANDH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt       ! Surface concs
              _DIAG_VAR_S_(data%id_eh(lay))  = ANCt        ! Surface concs
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.025127245
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.128028612
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.162224986
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.040779625
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.109189895
              _DIAG_VAR_S_(data%id_eh(lay))  = 2.741532621
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.056154803
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.44407992
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.087351916
              _DIAG_VAR_S_(data%id_eh(lay))  = 4.750826297
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 4.057572674
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.092343454
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.874113298
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.690556051
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_eh(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_eh(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.945188106
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.791124722
           END IF
         END DO

       END IF

     ! Currency == 2
     ELSEIF( data%zOM(zindex)>1.5 .AND. data%zOM(zindex) <2.5) THEN

       IF(theSoil%Substrate == OMCLAYL .OR. theSoil%Substrate == OMCLAYH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt      ! Surface concs
              _DIAG_VAR_S_(data%id_eh(lay))  = ANCt        ! Surface concs
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.582346104
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.094914577
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.808629162
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.092250308
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.579018412
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.108235921
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.290133148
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.275751823
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.288573293
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.274895451
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.106070183
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.517534219
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.172624024
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.716688314
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.205900944
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.927165551
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.287013437
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.119892097
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.287013437
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.119892097
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.089919073
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) =  0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  =  0.089919073
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) =  0.155985564
              _DIAG_VAR_S_(data%id_eh(lay))  =  0.089919073
           END IF
         END DO

       ELSE IF(theSoil%Substrate == OMSANDL .OR. theSoil%Substrate == OMSANDH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt       ! Surface concs
              _DIAG_VAR_S_(data%id_eh(lay))  = ANCt        ! Surface concs
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.026740382
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.188401867
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.03327692
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.254371066
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.043675958
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.365670896
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.047419611
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.463183135
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.054074995
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.744663137
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.147666334
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.731341792
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.2214995
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.006094515
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.2214995
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.006094515
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.215260078
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.670729677
           END IF
         END DO

       END IF


     ! Lake Albert == 3
!CAB was KASS
     ELSEIF( data%zOM(zindex)>2.5 .AND. data%zOM(zindex) <3.5) THEN

       IF(theSoil%Substrate == OMCLAYL .OR. theSoil%Substrate == OMCLAYH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt       ! Surface concs
              _DIAG_VAR_S_(data%id_eh(lay))  = ANCt        ! Surface concs
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 1.009705089
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.187140007
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 1.000955436
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.213892035
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.83412953
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.487488851
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.667303624
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.761085667
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.500477718
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.034682483
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.333651812
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.308279299
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.166825906
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.581876115
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.215260078
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.855472932
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.205900944
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.01769224
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.190302388
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.965818518
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.187182676
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.342130858
           END IF
         END DO

       ELSE IF(theSoil%Substrate == OMSANDL .OR. theSoil%Substrate == OMSANDH) THEN

         DO lay = 1,theSoil%nlay
           IF(middep(lay)-dep(1) >0.0 .AND. middep(lay)-dep(1) <=0.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = PASSt       ! Surface concs
              _DIAG_VAR_S_(data%id_eh(lay))  = ANCt        ! Surface concs
           ELSE IF(middep(lay)-dep(1) >0.2 .AND. middep(lay)-dep(1) <=0.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.025127245
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.128028612
           ELSE IF(middep(lay)-dep(1) >0.3 .AND. middep(lay)-dep(1) <=0.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.4 .AND. middep(lay)-dep(1) <=0.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.039191336
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.151751814
           ELSE IF(middep(lay)-dep(1) >0.5 .AND. middep(lay)-dep(1) <=0.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.6 .AND. middep(lay)-dep(1) <=0.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.7 .AND. middep(lay)-dep(1) <=0.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.093591338
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.272543826
           ELSE IF(middep(lay)-dep(1) >0.8 .AND. middep(lay)-dep(1) <=0.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.162224986
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.040779625
           ELSE IF(middep(lay)-dep(1) >0.9 .AND. middep(lay)-dep(1) <=1.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.109189895
              _DIAG_VAR_S_(data%id_eh(lay))  = 2.741532621
           ELSE IF(middep(lay)-dep(1) >1.0 .AND. middep(lay)-dep(1) <=1.1) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.056154803
              _DIAG_VAR_S_(data%id_eh(lay))  = 5.44407992
           ELSE IF(middep(lay)-dep(1) >1.1 .AND. middep(lay)-dep(1) <=1.2) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.087351916
              _DIAG_VAR_S_(data%id_eh(lay))  = 4.750826297
           ELSE IF(middep(lay)-dep(1) >1.2 .AND. middep(lay)-dep(1) <=1.3) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.118549028
              _DIAG_VAR_S_(data%id_eh(lay))  = 4.057572674
           ELSE IF(middep(lay)-dep(1) >1.3 .AND. middep(lay)-dep(1) <=1.4) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.092343454
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.874113298
           ELSE IF(middep(lay)-dep(1) >1.4 .AND. middep(lay)-dep(1) <=1.5) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_eh(lay))  = 3.690556051
           ELSE IF(middep(lay)-dep(1) >1.5 .AND. middep(lay)-dep(1) <=1.6) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_eh(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.6 .AND. middep(lay)-dep(1) <=1.7) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.047835573
              _DIAG_VAR_S_(data%id_eh(lay))  = 2.724078941
           ELSE IF(middep(lay)-dep(1) >1.7 .AND. middep(lay)-dep(1) <=1.8) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.065513937
              _DIAG_VAR_S_(data%id_eh(lay))  = 1.945188106
           ELSE IF(middep(lay)-dep(1) >1.8 .AND. middep(lay)-dep(1) <=1.9) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >1.9 .AND. middep(lay)-dep(1) <=2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.791124722
           ELSE IF(middep(lay)-dep(1) >2.0) THEN
              _DIAG_VAR_S_(data%id_pom0(lay)) = 0.012478845
              _DIAG_VAR_S_(data%id_eh(lay))  = 0.791124722
           END IF
         END DO

       END IF

     ! Constant
     ELSE

     END IF


END SUBROUTINE SetSoilPOMProfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




END MODULE aed2_soilbgc
