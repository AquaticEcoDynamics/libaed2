!###############################################################################
!#                                                                             #
!# aed2_habitat.F90                                                            #
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
!# Created March 2016                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_habitat
!-------------------------------------------------------------------------------
! aed2_habitat --- habitat model
!
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_habitat_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_habitat_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_bird, id_mtox
      INTEGER :: id_chsi, id_chpl, id_chfl, id_chsd
      INTEGER :: id_rhsi, id_rhpl, id_rhfl, id_rhsd, id_rhtr, id_rhsp
      INTEGER :: id_d_rupfs,id_d_rupft,id_d_rupfl,id_d_rupfa,id_d_rupfd

      !# Dependencies
      INTEGER :: id_l_ph, id_l_hab, id_l_aass, id_l_rveg, id_l_bveg
      INTEGER :: id_l_salg, id_l_falg
      INTEGER :: id_l_otrc, id_l_oxy, id_l_sav
      INTEGER, ALLOCATABLE :: id_l_mtox(:)
      !# Environment variables
      INTEGER :: id_E_temp, id_E_salt, id_E_bathy, id_E_matz, id_E_depth, id_E_nearlevel, id_E_extc, id_E_Io

      !# Model parameters
      LOGICAL :: simBirdForaging,simBenthicProd,simFishTolerance,simMetalTox,simCyanoRisk !,simMBORisk, simRuppiaGerm
      LOGICAL :: simCrabHabitat,simRuppiaHabitat
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      INTEGER :: num_mtox

     CONTAINS
         PROCEDURE :: define             => aed2_define_habitat
         PROCEDURE :: calculate          => aed2_calculate_habitat
         PROCEDURE :: calculate_benthic  => aed2_calculate_benthic_habitat
         PROCEDURE :: calculate_riparian => aed2_calculate_riparian_habitat
!        PROCEDURE :: mobility           => aed2_mobility_habitat
!        PROCEDURE :: light_extinction   => aed2_light_extinction_habitat
!        PROCEDURE :: delete             => aed2_delete_habitat

   END TYPE

   LOGiCAL :: extra_diag

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_habitat(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_habitat_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, status, num_habitats, num_mtox
   LOGICAL :: simBirdForaging, simBenthicProd, simFishTolerance, simMetalTox, simCyanoRisk, simCrabHabitat, simRuppiaHabitat
   AED_REAL          :: mtox_lims(10)
   CHARACTER(len=64) :: bird_acid_link, bird_habs_link, bird_aass_link, bird_rveg_link, bird_bveg_link
   CHARACTER(len=64) :: rhsi_salg_link, rhsi_falg_link, chsi_otrc_link, chsi_oxy_link, chsi_veg_link
   CHARACTER(len=40) :: mtox_acid_link, mtox_aass_link, mtox_vars(10)


   NAMELIST /aed2_habitat/ simBirdForaging, simBenthicProd, simFishTolerance, &
                           simMetalTox, mtox_vars, mtox_lims,   &
                           simCyanoRisk, simCrabHabitat, &
                           simRuppiaHabitat, extra_diag !, rhsi_falg_link,
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_habitat initialization"
   print *,"  WARNING! aed2_habitat model is currently under development"

   ! Default
   simBirdForaging = .false.
   simBenthicProd = .false.
   simFishTolerance = .false.
   simMetalTox = .false.
   simCyanoRisk = .false.
   simCrabHabitat = .false.
   simRuppiaHabitat = .false.

   extra_diag = .false.

   ! Read the namelist
   read(namlst,nml=aed2_habitat,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_habitat'

   ! Store parameter values in our own derived type
!  data%num_habitats = num_habitats
   data%simFishTolerance = simFishTolerance
   data%simBirdForaging = simBirdForaging
   data%simMetalTox = simMetalTox
   data%simBenthicProd = simBenthicProd
   data%simCyanoRisk = simCyanoRisk
   data%simCrabHabitat = simCrabHabitat
   data%simRuppiaHabitat = simRuppiaHabitat


   ! Define variables and dependencies
   IF( simBirdForaging ) THEN
     data%id_bird =  aed2_define_sheet_diag_variable('bird_wading','%', 'Suitability')

     bird_acid_link = 'CAR_pH'
     bird_habs_link = 'PHY_grn'
     bird_aass_link = 'ASS_uzaass'
     bird_rveg_link = '' ! 'VEG_xxx'
     bird_bveg_link = '' ! 'MAC_xxx'

     data%id_l_ph    = aed2_locate_global(TRIM(bird_acid_link))
     data%id_l_hab   = aed2_locate_global(TRIM(bird_habs_link))
     data%id_l_aass  = aed2_locate_global_sheet(TRIM(bird_aass_link))
     data%id_l_rveg  = 0 !aed2_locate_global_sheet(TRIM(bird_rveg_link))
     data%id_l_bveg  = 0 !aed2_locate_global_sheet(TRIM(bird_rveg_link))
   ENDIF
   IF( simMetalTox ) THEN
     data%id_mtox =  aed2_define_sheet_diag_variable('toxicity','-', 'Suitability')

     mtox_acid_link = 'CAR_pH'
     mtox_aass_link = 'ASS_uzaass'

     mtox_vars = '' ;  mtox_lims = 1.0
     DO i=1,10 ; IF (mtox_vars(i)  .EQ. '' ) THEN ; num_mtox = i-1 ; EXIT ; ENDIF ; ENDDO
     ALLOCATE(data%id_l_mtox(num_mtox)); ALLOCATE(data%mtox_lims(num_mtox))
     data%num_mtox = num_mtox
     DO i=1,data%num_mtox
       data%id_l_mtox(i) =  aed2_locate_variable(mtox_vars(i))
       data%mtox_lims(i) =  mtox_lims(i)
       !print*,'Tox : ', TRIM(tfe_vars(i)), ' * ', data%tfe_varscale(i)
     ENDDO
     IF ( .NOT.simBirdForaging ) data%id_l_ph   = aed2_locate_global(TRIM(mtox_acid_link))
     IF ( .NOT.simBirdForaging ) data%id_l_aass = aed2_locate_global_sheet(TRIM(mtox_aass_link))
   ENDIF
   IF( simCrabHabitat ) THEN
     data%id_chsi =  aed2_define_sheet_diag_variable('crab_hsi','-', 'Crab Habitat Suitability Index')
     data%id_chpl =  aed2_define_sheet_diag_variable('crab_hsi_larvae','-', 'Crab Habitat Suitability - larval connectivity')
     data%id_chfl =  aed2_define_sheet_diag_variable('crab_hsi_wq','-', 'Crab Habitat Suitability - water quality')
     data%id_chsd =  aed2_define_sheet_diag_variable('crab_hsi_ben','-', 'Crab Habitat Suitability - benthic habitat')

     chsi_otrc_link = 'TRC_ss1'   ! ocean larvae tracer
     chsi_oxy_link = 'OXY_oxy'    ! oxygen
     chsi_veg_link = 'MAC_sav'    ! submerged aquatic vegetation

     data%id_l_otrc  = aed2_locate_global(TRIM(chsi_otrc_link))
     data%id_l_oxy  = aed2_locate_global(TRIM(chsi_oxy_link))
     data%id_l_sav  = aed2_locate_global_sheet(TRIM(chsi_veg_link))
   ENDIF
   IF( simRuppiaHabitat ) THEN
     data%id_rhsi =  aed2_define_sheet_diag_variable('ruppia_hsi','-', 'Ruppia Habitat Suitability Index')
     data%id_rhpl =  aed2_define_sheet_diag_variable('ruppia_hsi_plant','-', 'Ruppia Habitat Suitability - plant')
     data%id_rhfl =  aed2_define_sheet_diag_variable('ruppia_hsi_flower','-', 'Ruppia Habitat Suitability - flowering')
     data%id_rhsd =  aed2_define_sheet_diag_variable('ruppia_hsi_seed','-', 'Ruppia Habitat Suitability - seed germination')
     data%id_rhtr =  aed2_define_sheet_diag_variable('ruppia_hsi_turion','-', 'Ruppia Habitat Suitability - turion')
     data%id_rhsp =  aed2_define_sheet_diag_variable('ruppia_hsi_sprout','-', 'Ruppia Habitat Suitability - sprout')

     rhsi_falg_link = 'MAG_ulva_ben'
     rhsi_salg_link = 'MAG_ulva'

     data%id_l_salg  = aed2_locate_global(TRIM(rhsi_salg_link))
     data%id_l_falg  = aed2_locate_global_sheet(TRIM(rhsi_falg_link))

     if( extra_diag )then
       data%id_d_rupfs = aed2_define_sheet_diag_variable('ruppia_hsi_fsal','-', 'Ruppia Habitat Suitability - fSal')
       data%id_d_rupft = aed2_define_sheet_diag_variable('ruppia_hsi_ftem','-', 'Ruppia Habitat Suitability - fTem')
       data%id_d_rupfl = aed2_define_sheet_diag_variable('ruppia_hsi_flgt','-', 'Ruppia Habitat Suitability - fLgt')
       data%id_d_rupfa = aed2_define_sheet_diag_variable('ruppia_hsi_falg','-', 'Ruppia Habitat Suitability - fAlg')
       data%id_d_rupfd = aed2_define_sheet_diag_variable('ruppia_hsi_fdep','-', 'Ruppia Habitat Suitability - fDep')
     endif

   ENDIF


   ! Register environmental dependencies
   data%id_E_salt = aed2_locate_global('salinity')
   data%id_E_extc = aed2_locate_global('extc_coef')
   data%id_E_temp = aed2_locate_global('temperature')
   data%id_E_matz = aed2_locate_global_sheet('material')
   data%id_E_bathy = aed2_locate_global_sheet('bathy')
   data%id_E_depth = aed2_locate_global('layer_ht')
   IF(simBirdForaging) data%id_E_nearlevel = aed2_locate_global_sheet('nearest_depth')


END SUBROUTINE aed2_define_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_habitat(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: trc
!
!-------------------------------------------------------------------------------
!BEGIN
!  DO i=1,data%num_habitats
!     trc = _STATE_VAR_(data%id_ss(i))
!     _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + data%decay(i)*trc
!  ENDDO

!  IF (data%id_retain < 1) RETURN
!  _FLUX_VAR_(data%id_retain) = _FLUX_VAR_(data%id_retain) + 1.0
END SUBROUTINE aed2_calculate_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_habitat(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED habitat.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss, bottom_stress

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0, resus_flux = 0., dummy_tau
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN


END SUBROUTINE aed2_calculate_benthic_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian_habitat(data,column,layer_idx, pc_wet)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED habitat.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wlevel, extc, bathy, matz

   ! State
   AED_REAL :: depth, ph, hab, sdepth, uzaass, aass, conc, mtox

   ! Temporary variables
   INTEGER  :: i
   AED_REAL :: bird_acid, bird_soil, bird_rveg, bird_bveg, bird_salt, bird_dept, bird_habs
   AED_REAL :: euphotic

   ! Parameters
   AED_REAL, PARAMETER :: crit_soil_acidity = 100.000
   AED_REAL, PARAMETER :: crit_water_ph = 6.0
   AED_REAL, PARAMETER :: crit_soil_type = 5.5 ! (matz 6&7 is sand)
   AED_REAL, PARAMETER :: crit_rveg_depth = 0.6
   AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
   AED_REAL, PARAMETER :: crit_salinity = 50.0
   AED_REAL, PARAMETER :: crit_leg_depth = 0.12
   AED_REAL, PARAMETER :: crit_hab_conc = 500.

   AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp,falg,Io,vel
   AED_REAL :: limitation(5,6)

!-------------------------------------------------------------------------------
!BEGIN
   !-- HABITAT TEMPLATE 1: Bird Foraging
   IF( data%simBirdForaging ) THEN

     bird_acid=one_; bird_soil=one_; bird_rveg=one_;
     bird_bveg=one_; bird_salt=one_; bird_dept=one_; bird_habs=one_

     IF( data%id_l_ph>0 .AND. pc_wet>0.1 ) THEN
       ph = _STATE_VAR_(data%id_l_pH)
     ELSE
       ph = 8.0
     ENDIF
     IF( data%id_E_depth>0 ) THEN
       depth = _STATE_VAR_(data%id_E_depth)
     ELSE
       depth = zero_
     ENDIF
     sdepth = zero_
     IF( pc_wet < 0.1 ) THEN
       bathy = _STATE_VAR_S_(data%id_E_bathy)
       wlevel = _STATE_VAR_S_(data%id_E_nearlevel)
       sdepth = bathy - wlevel
       depth = zero_
       aass = zero_
     ELSE
       extc = _STATE_VAR_(data%id_E_extc)
       euphotic = log(0.9) / (-1.*extc)
       aass = zero_
     ENDIF
     IF( data%id_l_hab>0 ) THEN
       hab = _STATE_VAR_(data%id_l_hab)
     ELSE
       hab = zero_
     ENDIF

     ! Effect of acid
     IF( ph<crit_water_ph .OR. aass>crit_soil_acidity) bird_acid = zero_

     ! Requirement for mud-flat
     IF( matz > crit_soil_type ) bird_soil = zero_

     ! Presence of riparian/emergent vegetation (veg is located >0.6mAHD)
     IF( bathy < crit_rveg_depth ) bird_rveg = zero_

     ! Presence or liklihood of benthic vegetation
     IF( (depth-sav_height) > euphotic .OR. pc_wet < 0.1 ) bird_bveg = zero_

     ! Check for tolerable salinity
     IF( salt > crit_salinity ) bird_salt = zero_

     ! Check for relevant leg/beak depth
     IF( depth > crit_leg_depth .OR. sdepth>crit_leg_depth ) bird_dept = zero_

     ! Limit suitability of Harmful Algal Bloom (HAB)
     IF( hab > crit_hab_conc ) bird_habs = zero_

     _DIAG_VAR_S_(data%id_bird) =  (bird_acid+bird_soil+bird_rveg+bird_bveg+bird_habs)/5. &
                                   * bird_salt * bird_dept

   ENDIF


   !-- HABITAT TEMPLATE 2: Metal Toxicity
   IF( data%simMetalTox ) THEN

     mtox = zero_
     DO  i=1,data%num_mtox
        conc = _STATE_VAR_( data%id_l_mtox(i) )
        IF ( conc > data%mtox_lims(i) ) mtox = one_
     ENDDO
     IF( pc_wet < 0.1 ) THEN
        aass = _STATE_VAR_S_( data%id_l_aass )
        IF ( aass < crit_soil_acidity ) mtox = one_
     ELSE
        aass = _STATE_VAR_( data%id_l_ph )
        IF ( aass < crit_water_ph ) mtox = one_
     END IF

     _DIAG_VAR_S_(data%id_mtox) = mtox

   ENDIF


   !-- HABITAT TEMPLATE 3: Fish Tolerance
   IF( data%simFishTolerance ) THEN

   !  _DIAG_VAR_S_(data%id_fish) = mtox

   ENDIF


   !-- HABITAT TEMPLATE 4: Crab Habitat
   IF( data%simCrabHabitat ) THEN

   !  _DIAG_VAR_S_(data%id_fish) = mtox

   ENDIF

   !-- HABITAT TEMPLATE 4: Ruppia Habitat
   IF( data%simRuppiaHabitat ) THEN

     depth = _STATE_VAR_(data%id_E_depth)  ! metres
     salt  = _STATE_VAR_(data%id_E_salt)   ! salinity g/L
     temp  = _STATE_VAR_(data%id_E_temp)   ! degC
     extc  = _STATE_VAR_(data%id_E_salt)   ! /m
     Io    = _STATE_VAR_S_(data%id_E_Io)   ! W/m2
     falg  =(_STATE_VAR_S_(data%id_l_falg) + _STATE_VAR_(data%id_l_salg)*depth) * 12. * 1e-3 / 0.5  ! convert mmolC/m2 to gDW/m2
     vel   = 0. !

     CALL ruppia_habitat_suitability(data, &
                                     rhpl,rhfl,rhsd,rhtr,rhsp,&
                                     depth,salt,temp,extc,falg,Io,vel,pc_wet,&
                                     limitation)

     _DIAG_VAR_S_(data%id_rhpl) = rhpl
     _DIAG_VAR_S_(data%id_rhfl) = rhfl
     _DIAG_VAR_S_(data%id_rhsd) = rhsd
     _DIAG_VAR_S_(data%id_rhtr) = rhtr
     _DIAG_VAR_S_(data%id_rhsp) = rhsp

     ! Overall HSI : Habitat Suitability Index (Issue here re time integration)
     _DIAG_VAR_S_(data%id_rhsi) = (rhpl+rhfl+rhsd+rhtr+rhsp)/5.

     iF( extra_diag ) THEN
       _DIAG_VAR_S_(data%id_d_rupfs) = limitation(1,1)
       _DIAG_VAR_S_(data%id_d_rupft) = limitation(1,2)
       _DIAG_VAR_S_(data%id_d_rupfl) = limitation(1,3)
       _DIAG_VAR_S_(data%id_d_rupfa) = limitation(1,4)
       _DIAG_VAR_S_(data%id_d_rupfd) = limitation(1,5)
      !_DIAG_VAR_S_(data%id_d_rupfm) = limitation(1,6)
     ENDiF
   ENDIF


END SUBROUTINE aed2_calculate_riparian_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_habitat(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: ss
   INTEGER  :: ss_i
!
!-----------------------------------------------------------------------
!BEGIN
!  DO ss_i=1,ubound(data%id_ss,1)
!     ! Retrieve current (local) state variable values.
!     ss = _STATE_VAR_(data%id_ss(ss_i))

!     ! Self-shading with explicit contribution from background habitat concentration.
!     extinction = extinction + (data%Ke_ss(ss_i)*ss)
!  ENDDO
END SUBROUTINE aed2_light_extinction_habitat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE ruppia_habitat_suitability(data,rhpl,rhfl,rhsd,rhtr,rhsp,depth,salt,temp,extc,fa,Io,vel,pc_wet,limitation)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_habitat_data_t),INTENT(in) :: data
   AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp
   AED_REAL :: depth,salt,temp,extc,fa,Io,vel,pc_wet
   AED_REAL :: limitation(:,:)
!
!LOCALS
   AED_REAL :: rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess
   AED_REAL :: light
!
!-----------------------------------------------------------------------
!BEGIN

     rupp_salt=one_; rupp_temp=one_; rupp_lght=one_;
     rupp_falg=one_; rupp_matz=one_; rupp_dess=one_;

     IF( depth<0.1 ) THEN
       light = 100.
     ELSE
       light = 100. * exp(-extc*(depth-0.08))
     ENDIF

     !-- First do ADULT tolerance
     IF( pc_wet < 0.1 ) THEN
       ! Dry cell - set dessication factor
       rupp_dess = zero_    ! maybe need a time counter here.

     ELSE
       ! Wet cell
       rupp_salt = ruppia_salinity(salt, "adult")
       rupp_temp = ruppia_temp(temp, "adult")
       rupp_lght = ruppia_light(light, "adult")
       rupp_falg = ruppia_filalgae(fa, "adult")
       rupp_dess = ruppia_depth(depth, "adult")

     ENDIF
     ! Adult plant habitat suitability
     rhpl = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)

     limitation(1,1) = rupp_salt
     limitation(1,2) = rupp_temp
     limitation(1,3) = rupp_lght
     limitation(1,4) = rupp_falg
     limitation(1,5) = rupp_dess
     limitation(1,6) = rupp_matz


     !-- Second do FLOWER tolerance
     IF( pc_wet < 0.1 ) THEN
       ! Dry cell - set dessication factor
       rupp_dess = zero_    ! maybe need a time counter here.

     ELSE
       ! Wet cell
       rupp_salt = ruppia_salinity(salt, "flower")
       rupp_temp = ruppia_temp(temp, "flower")
       rupp_lght = ruppia_light(light, "flower")
       rupp_falg = ruppia_filalgae(fa, "flower")
       rupp_dess = ruppia_depth(depth, "flower")

     ENDIF
     ! Habitat suitability for flowering
     rhfl = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)

     !-- Third do seed germination
     IF( pc_wet < 0.1 ) THEN
       ! Dry cell - set dessication factor

       rupp_dess = 0.5    ! maybe need a time counter here.

     ELSE
       ! Wet cell
       rupp_salt = ruppia_salinity(salt, "seed")
       rupp_temp = ruppia_temp(temp, "seed")
       rupp_lght = one_
       rupp_falg = one_
       rupp_dess = ruppia_depth(depth, "seed")

     ENDIF
     ! Habitat suitability for seed germination
     rhsd = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)


     !-- Fourth do sediment suitability for turions/turion sprouting
     IF( pc_wet < 0.1 ) THEN
       ! Dry cell - set dessication factor

       rupp_dess = zero_    ! maybe need a time counter here.

     ELSE
       ! Wet cell
       rupp_salt = ruppia_salinity(salt, "turion")
       rupp_temp = ruppia_temp(temp, "turion")
       rupp_lght = one_
       rupp_falg = one_
       rupp_dess = ruppia_depth(depth, "turion")

     ENDIF
     !IF( .NOT. _STATE_VAR_S_(data%id_E_matz)==in_zone_set ) rupp_matz = zero_
     ! Habitat suitability for seed germination
    rhtr = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)

    !-- Fifth do sprout tolerance
    IF( pc_wet < 0.1 ) THEN
      ! Dry cell - set dessication factor

      rupp_dess = zero_    ! maybe need a time counter here.

    ELSE
      ! Wet cell
      rupp_salt = ruppia_salinity(salt, "sprout")
      rupp_temp = ruppia_temp(temp, "sprout")
      rupp_lght = ruppia_light(light, "sprout")
      rupp_falg = ruppia_filalgae(fa, "sprout")
      rupp_dess = ruppia_depth(depth, "sprout")

    ENDIF
    ! Habitat suitability for flowering
    rhsp = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)

  !---------------------------------------------------------------------
  CONTAINS

  !#############################################################################
  AED_REAL FUNCTION ruppia_salinity(salt,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: salt
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

     ruppia_salinity = one_

     IF( TRIM(stage)=="seed" ) THEN

       !<19 unsuitable
       !19-40 suboptimal
       !40-60 optimal
       !60-85 suboptimal
       !>85 unsuitable
       IF( salt<=19. ) THEN
         ruppia_salinity = zero_
       ELSE IF ( salt>19. .AND. salt<=40.  ) THEN
         ruppia_salinity = 0. + ( (salt-19.)/(40.-19.) )
       ELSE IF ( salt>40. .AND. salt<=60. ) THEN
         ruppia_salinity = one_
       ELSE IF ( salt>60. .AND. salt<=85. ) THEN
         ruppia_salinity = 1. - ( (salt-60.)/(85.-60.) )
       ELSE IF ( salt>85. ) THEN
         ruppia_salinity = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="sprout" ) THEN

       !<20 suboptimal
       !20-75 optimal
       !75 - 130 suboptimal
       !>130 unsuitable
       IF( salt<=0. ) THEN
         ruppia_salinity = zero_
       ELSE IF ( salt>0. .AND. salt<=20.  ) THEN
         ruppia_salinity = 0. + ( (salt-0.)/(20.-0.) )
       ELSE IF ( salt>20. .AND. salt<=75. ) THEN
         ruppia_salinity = one_
       ELSE IF ( salt>75. .AND. salt<=130. ) THEN
         ruppia_salinity = 1. - ( (salt-75.)/(130.-75.) )
       ELSE IF ( salt>130. ) THEN
         ruppia_salinity = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="adult" ) THEN

       !0 - 10 unsuitable
       !10 - 71 suboptimal
       !72 - 123 optimal
       !123 - 230 suboptimal
       !>230 unsuitable
       IF( salt<=10. ) THEN
         ruppia_salinity = zero_
       ELSE IF ( salt>10. .AND. salt<=72.  ) THEN
         ruppia_salinity = 0. + ( (salt-10.)/(72.-10.) )
       ELSE IF ( salt>72. .AND. salt<=123. ) THEN
         ruppia_salinity = one_
       ELSE IF ( salt>123. .AND. salt<=230. ) THEN
         ruppia_salinity = 1. - ( (salt-123.)/(230.-123.) )
       ELSE IF ( salt>230. ) THEN
         ruppia_salinity = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="flower" ) THEN

       !<47 unsuitable
       !47-50 suboptimal
       !50 - 62 optimal
       !62 - 70 suboptimal
       !>70 unsuitable
       IF( salt<=47. ) THEN
         ruppia_salinity = zero_
       ELSE IF ( salt>47. .AND. salt<=50.  ) THEN
         ruppia_salinity = 0. + ( (salt-47.)/(50.-47.) )
       ELSE IF ( salt>50. .AND. salt<=62. ) THEN
         ruppia_salinity = one_
       ELSE IF ( salt>62. .AND. salt<=70. ) THEN
         ruppia_salinity = 1. - ( (salt-62.)/(70.-62.) )
       ELSE IF ( salt>70. ) THEN
         ruppia_salinity = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="turion" ) THEN

       !<70 unsuitable
       !70 - 124 suboptimal
       !124 - 160 optimal
       !160 - 230 suboptimal
       !>230 unsuitable
       IF( salt<=70. ) THEN
         ruppia_salinity = zero_
       ELSE IF ( salt>70. .AND. salt<=124.  ) THEN
         ruppia_salinity = 0. + ( (salt-70.)/(124.-70.) )
       ELSE IF ( salt>124. .AND. salt<=160. ) THEN
         ruppia_salinity = one_
       ELSE IF ( salt>160. .AND. salt<=230. ) THEN
         ruppia_salinity = 1. - ( (salt-160.)/(230.-160.) )
       ELSE IF ( salt>230. ) THEN
         ruppia_salinity = zero_
       ENDIF

     ENDIF

  END FUNCTION ruppia_salinity
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !#############################################################################
  AED_REAL FUNCTION ruppia_temp(temp,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: temp
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

     ruppia_temp = one_

     IF( TRIM(stage)=="seed"   .OR. &
         TRIM(stage)=="sprout" .OR. &
         TRIM(stage)=="adult"  .OR. &
         TRIM(stage)=="flower" .OR. &
         TRIM(stage)=="turion"  ) THEN

       !<4 unsuitable
       !4 - 10 suboptimal
       !10 - 20 optimal
       !20-30 suboptimal
       !>30 unsuitable
       IF( temp<=4. ) THEN
         ruppia_temp = zero_
       ELSE IF ( temp>4. .AND. temp<=10.  ) THEN
         ruppia_temp = 0. + ( (temp-4.)/(10.-4.) )
       ELSE IF ( temp>10. .AND. temp<=20. ) THEN
         ruppia_temp = one_
       ELSE IF ( temp>20. .AND. temp<=30. ) THEN
         ruppia_temp = 1. - ( (temp-20.)/(30.-20.) )
       ELSE IF ( temp>30. ) THEN
         ruppia_temp = zero_
       ENDIF

     ENDIF

  END FUNCTION ruppia_temp
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !#############################################################################
  AED_REAL FUNCTION ruppia_light(light,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: light
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

     ruppia_light = one_

     IF( TRIM(stage)=="sprout" .OR. &
         TRIM(stage)=="adult"  ) THEN

       !0 - 7.5 unsuitable
       !7.5 - 24 suboptimal
       !>24 optimal
       IF( light<=7.5 ) THEN
         ruppia_light = zero_
       ELSE IF ( light>7.5 .AND. light<=24.  ) THEN
         ruppia_light = 0. + ( (light-7.5)/(24.-7.5) )
       ELSE IF ( light>24. ) THEN
         ruppia_light = one_
       ENDIF

     ENDIF

  END FUNCTION ruppia_light
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !#############################################################################
  AED_REAL FUNCTION ruppia_filalgae(fa,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: fa
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

     ruppia_filalgae = one_

     IF( TRIM(stage)=="sprout" .OR. &
         TRIM(stage)=="adult"  .OR. &
         TRIM(stage)=="flower" ) THEN

       IF ( fa>=0. .AND. fa<=50. ) THEN
         ruppia_filalgae = one_
       ELSE IF ( fa>50. .AND. fa<=100. ) THEN
         ruppia_filalgae = 1. - ( (fa-50.)/(100.-50.) )
       ELSE IF ( fa>100. ) THEN
         ruppia_filalgae = zero_
       ENDIF

     ENDIF

  END FUNCTION ruppia_filalgae
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !#############################################################################
  AED_REAL FUNCTION ruppia_depth(depth,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: depth
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

     ruppia_depth = one_

     IF( TRIM(stage)=="adult" ) THEN

       IF ( depth>=0. .AND. depth<=0.1 ) THEN
         ruppia_depth = zero_
       ELSE IF ( depth>0.1 .AND. depth<=0.3 ) THEN
         ruppia_depth = 1. - ( (depth-0.1)/(0.3-0.1) )
       ELSE IF ( depth>0.3 ) THEN
         ruppia_depth = one_
       ENDIF

     ELSEIF( TRIM(stage)=="sprout" ) THEN
 
       IF ( depth<=0.01 ) THEN
         ruppia_depth = zero_
       ELSE IF ( depth>0.01 .AND. depth<=0.2 ) THEN
         ruppia_depth = 1. - ( (depth-0.01)/(0.2-0.01) )
       ELSE IF ( depth>0.2 ) THEN
         ruppia_depth = one_
       ENDIF
     
     ELSEIF( TRIM(stage)=="flower" ) THEN
 
       IF ( depth<=0.01 ) THEN
         ruppia_depth = zero_
       ELSE IF ( depth>0.01 .AND. depth<=0.1 ) THEN
         ruppia_depth = 1. - ( (depth-0.01)/(0.1-0.01) )
       ELSE IF ( depth>0.1 ) THEN
         ruppia_depth = one_
       ENDIF
     
     ENDIF

  END FUNCTION ruppia_depth
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE ruppia_habitat_suitability
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed2_habitat
