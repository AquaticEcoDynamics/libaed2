!###############################################################################
!#                                                                             #
!# aed2_organic_matter.F90                                                     #
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
!# Created 6 June 2011                                                         #
!# Updated 11 May 2015 : Add refractory groups, based on Swan estuary model    #
!# Updated 18 Oct 2016 : Add resuspension of POM and cleaned up input          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_organic_matter
!-------------------------------------------------------------------------------
! aed2_organic_matter --- organic matter biogeochemical model
!
! The Organic Matter (OM) module contains equations for mineralisation
! and breakdown of particulate and dissolved organic matter pools, as well as
! settling, resuspension and other processes
!
! Users can optionally configure refractory pools
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_util,ONLY : water_viscosity

   IMPLICIT NONE

   PRIVATE  ! By default make everything private
!
   PUBLIC aed2_organic_matter_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_organic_matter_data_t

      !# Variable identifiers
      INTEGER  :: id_poc,id_doc           ! particulate & dissolved organic carbon
      INTEGER  :: id_pon,id_don           ! particulate & dissolved organic nitrogen
      INTEGER  :: id_pop,id_dop           ! particulate & dissolved organic phosphorus
      INTEGER  :: id_docr,id_donr,id_dopr ! dissolved refractory OM pools
      INTEGER  :: id_cpom                 ! coarse particulate OM (bulk) pool
      INTEGER  :: id_photolysis           ! photolysis

      INTEGER  :: id_oxy,id_amm,id_frp,id_dic
      INTEGER  :: id_Fsed_pon, id_Fsed_don ! sed. rate organic nitrogen
      INTEGER  :: id_Fsed_pop, id_Fsed_dop ! sed. rate organic phosphorus
      INTEGER  :: id_Fsed_poc, id_Fsed_doc ! sed. rate organic carbon
      INTEGER  :: id_Psed_poc, id_Psed_pon, id_Psed_pop, id_Psed_cpom ! sedimentation rates
      INTEGER  :: id_temp, id_salt, id_vis, id_uva, id_uvb, id_rho
      INTEGER  :: id_pon_miner, id_don_miner
      INTEGER  :: id_pop_miner, id_dop_miner
      INTEGER  :: id_poc_miner, id_doc_miner
      INTEGER  :: id_sed_pon, id_sed_don
      INTEGER  :: id_sed_pop, id_sed_dop
      INTEGER  :: id_sed_poc, id_sed_doc
      INTEGER  :: id_bod, id_cdom
      INTEGER  :: id_l_resus

      !# Model parameters
      AED_REAL :: Rpoc_hydrol, Rdoc_minerl, Rpon_hydrol, &
                          Rdon_minerl, Rpop_hydrol, Rdop_minerl, &
                          theta_hydrol, theta_minerl, Kpom_hydrol, Kdom_minerl
      AED_REAL :: Rdocr_miner, Rdonr_miner, Rdopr_miner, Rcpom_bdown, &
                          X_cpom_n, X_cpom_p
      AED_REAL :: KeDOM, KePOM, KeDOMR, KeCPOM, photo_fmin
      AED_REAL :: w_pom, d_pom, rho_pom, w_cpom, d_cpom, rho_cpom
      AED_REAL :: sedimentOMfrac, Xsc, Xsn, Xsp
      AED_REAL :: Ksed_dom, theta_sed_dom, Fsed_doc, Fsed_don, Fsed_dop

      !# Model options
      INTEGER  :: resuspension, settling
      LOGICAL  :: simRPools,extra_diag,simphotolysis
      LOGICAL  :: simSedimentOM
      LOGICAL  :: use_oxy,use_amm,use_frp,use_dic,use_Fsed_model,use_Psed_model

     CONTAINS
         PROCEDURE :: define            => aed2_define_organic_matter
         PROCEDURE :: calculate         => aed2_calculate_organic_matter
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_organic_matter
         PROCEDURE :: mobility          => aed2_mobility_organic_matter
         PROCEDURE :: light_extinction  => aed2_light_extinction_organic_matter
!        PROCEDURE :: delete            => aed2_delete_organic_matter

   END TYPE

   AED_REAL :: c

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_define_organic_matter(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                              :: namlst
   CLASS (aed2_organic_matter_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status
   !-- Default initial values and limits
   AED_REAL                  :: poc_initial  = 100.0
   AED_REAL                  :: doc_initial  = 100.0
   AED_REAL                  :: poc_min      = zero_
   AED_REAL                  :: poc_max      = nan_
   AED_REAL                  :: doc_min      = zero_
   AED_REAL                  :: doc_max      = nan_
   AED_REAL                  :: pon_initial  = 16.0
   AED_REAL                  :: don_initial  = 16.0
   AED_REAL                  :: pon_min      = zero_
   AED_REAL                  :: pon_max      = nan_
   AED_REAL                  :: don_min      = zero_
   AED_REAL                  :: don_max      = nan_
   AED_REAL                  :: pop_initial  = 1.0
   AED_REAL                  :: dop_initial  = 1.0
   AED_REAL                  :: pop_min      = zero_
   AED_REAL                  :: pop_max      = nan_
   AED_REAL                  :: dop_min      = zero_
   AED_REAL                  :: dop_max      = nan_
   AED_REAL                  :: docr_initial = zero_
   AED_REAL                  :: donr_initial = zero_
   AED_REAL                  :: dopr_initial = zero_
   AED_REAL                  :: cpom_initial = zero_
   !-- Breakdown and mineralisation (basic pool)
   AED_REAL                  :: Rpoc_hydrol
   AED_REAL                  :: Rdoc_minerl
   AED_REAL                  :: Rpon_hydrol
   AED_REAL                  :: Rdon_minerl
   AED_REAL                  :: Rpop_hydrol
   AED_REAL                  :: Rdop_minerl
   AED_REAL                  :: theta_hydrol
   AED_REAL                  :: theta_minerl
   AED_REAL                  :: Kpom_hydrol
   AED_REAL                  :: Kdom_minerl
   CHARACTER(len=64)         :: doc_miner_product_variable=''
   CHARACTER(len=64)         :: doc_miner_reactant_variable=''
   CHARACTER(len=64)         :: don_miner_product_variable=''
   CHARACTER(len=64)         :: dop_miner_product_variable=''
   !-- Refractory organic matter (optional)
   LOGICAL                   :: simRPools = .false.
   AED_REAL                  :: Rdocr_miner
   AED_REAL                  :: Rdonr_miner
   AED_REAL                  :: Rdopr_miner
   AED_REAL                  :: Rcpom_bdown
   AED_REAL                  :: X_cpom_n
   AED_REAL                  :: X_cpom_p
   CHARACTER(len=64)         :: docr_miner_product_variable='OGM_doc'
   CHARACTER(len=64)         :: donr_miner_product_variable='OGM_don'
   CHARACTER(len=64)         :: dopr_miner_product_variable='OGM_dop'
   !-- Light related parameters
   AED_REAL                  :: KeDOM  = 0.0001
   AED_REAL                  :: KePOM  = 0.0001
   AED_REAL                  :: KeDOMR = 0.0001
   AED_REAL                  :: KeCPOM = 0.0001
   LOGICAL                   :: simPhotolysis = .FALSE.
   AED_REAL                  :: photo_fmin = 0.9, photo_c = 7.52
   !-- Particle settling parameters
   INTEGER                   :: settling = 0
   AED_REAL                  :: w_pom
   AED_REAL                  :: d_pom
   AED_REAL                  :: rho_pom
   AED_REAL                  :: w_cpom
   AED_REAL                  :: d_cpom
   AED_REAL                  :: rho_cpom
   CHARACTER(len=64)         :: Psed_poc_variable=''
   CHARACTER(len=64)         :: Psed_pon_variable=''
   CHARACTER(len=64)         :: Psed_pop_variable=''
   CHARACTER(len=64)         :: Psed_cpom_variable=''
   !-- Sediment interaction parameters (basic model)
   INTEGER                   :: resuspension = 0
   AED_REAL                  :: sedimentOMfrac
   AED_REAL                  :: Xsc
   AED_REAL                  :: Xsn
   AED_REAL                  :: Xsp
   AED_REAL                  :: Fsed_doc
   AED_REAL                  :: Fsed_don
   AED_REAL                  :: Fsed_dop
   AED_REAL                  :: Ksed_dom
   AED_REAL                  :: theta_sed_dom
   CHARACTER(len=64)         :: Fsed_poc_variable=''
   CHARACTER(len=64)         :: Fsed_doc_variable=''
   CHARACTER(len=64)         :: Fsed_pon_variable=''
   CHARACTER(len=64)         :: Fsed_don_variable=''
   CHARACTER(len=64)         :: Fsed_pop_variable=''
   CHARACTER(len=64)         :: Fsed_dop_variable=''
   CHARACTER(len=64)         :: resus_link=''
   !-- Misc options
   LOGICAL                   :: extra_diag = .FALSE.

   !AED_REAL                  :: Rpoc_miner  = 0.01
   !AED_REAL                  :: Rdoc_miner  = 0.01
   !AED_REAL                  :: Fsed_doc    = 30.0
   !AED_REAL                  :: Kpoc_miner  = 30.0
   !AED_REAL                  :: Kdoc_miner  = 30.0
   !AED_REAL                  :: Ksed_doc    = 4.5
   !AED_REAL                  :: theta_poc_miner = 1.0
   !AED_REAL                  :: theta_doc_miner = 1.0
   !AED_REAL                  :: theta_sed_doc   = 1.0
   !AED_REAL                  :: w_pon       = 0.0
   !AED_REAL                  :: Rpon_miner  = 0.01
   !AED_REAL                  :: Rdon_miner  = 0.01
   !AED_REAL                  :: Fsed_pon    =  0.0
   !AED_REAL                  :: Fsed_don    = 30.0
   !AED_REAL                  :: Kpon_miner  = 30.0
   !AED_REAL                  :: Kdon_miner  = 30.0
   !AED_REAL                  :: Ksed_don    = 4.5
   !AED_REAL                  :: theta_pon_miner = 1.0
   !AED_REAL                  :: theta_don_miner = 1.0
   !AED_REAL                  :: theta_sed_don   = 1.0
   !AED_REAL                  :: w_pop       = 0.0
   !AED_REAL                  :: Rpop_miner  = 0.01
   !AED_REAL                  :: Rdop_miner  = 0.01
   !AED_REAL                  :: Fsed_pop    =  0.0
   !AED_REAL                  :: Fsed_dop    = 30.0
   !AED_REAL                  :: Kpop_miner  = 30.0
   !AED_REAL                  :: Kdop_miner  = 30.0
   !AED_REAL                  :: Ksed_dop    = 4.5
   !AED_REAL                  :: theta_pop_miner = 1.0
   !AED_REAL                  :: theta_dop_miner = 1.0
   !AED_REAL                  :: theta_sed_dop   = 1.0
   !AED_REAL                  :: Rdocr_miner  = 0.01
   !AED_REAL                  :: Rdonr_miner  = 0.01
   !AED_REAL                  :: Rdopr_miner  = 0.01
   !AED_REAL                  :: Rcpom_bdown  = 0.01
   !AED_REAL                  :: X_cpom_n     = 0.0
   !AED_REAL                  :: X_cpom_p     = 0.0

   NAMELIST /aed2_organic_matter/ &
             poc_min, poc_max, doc_min, doc_max, &
             pon_min, pon_max, don_min, don_max, &
             pop_min, pop_max, dop_min, dop_max, &
             poc_initial,doc_initial,pon_initial,&
             don_initial,pop_initial,dop_initial,&
             docr_initial, donr_initial, dopr_initial, cpom_initial, &
             Rpoc_hydrol,Rdoc_minerl,Rpon_hydrol,&
             Rdon_minerl,Rpop_hydrol,Rdop_minerl,&
             theta_hydrol,theta_minerl,Kpom_hydrol,Kdom_minerl, &
             doc_miner_reactant_variable, doc_miner_product_variable, &
             don_miner_product_variable,dop_miner_product_variable,&
             simRPools,Rdocr_miner,Rdonr_miner,Rdopr_miner,&
             Rcpom_bdown,X_cpom_n,X_cpom_p,&
             ! docr_miner_product_variable, donr_miner_product_variable,dopr_miner_product_variable, &
             KeDOM, KePOM, KeDOMR, KeCPOM, &
             simphotolysis, photo_fmin, photo_c, &
             settling,w_pom,d_pom,rho_pom,w_cpom,d_cpom,rho_cpom,&
             Psed_poc_variable, Psed_pon_variable, Psed_pop_variable,Psed_cpom_variable,&
             resuspension, resus_link, sedimentOMfrac, Xsc, Xsn, Xsp,&
             Fsed_doc,Fsed_don,Fsed_dop,Ksed_dom,theta_sed_dom,&
             Fsed_poc_variable, Fsed_doc_variable,&
             Fsed_pop_variable, Fsed_dop_variable,&
             Fsed_pon_variable, Fsed_don_variable,&
             extra_diag

!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_organic_matter initialization"

   ! Set default parameters here


   ! Read the namelist
   read(namlst,nml=aed2_organic_matter,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_organic_matter'

   ! Store parameter values in the modules dat structure
   !   NB: all rates are provided in values per day,
   !       and are converted here to values per second.

   !-- Breakdown and mineralisation (basic pool)
   data%Rpoc_hydrol     = Rpoc_hydrol/secs_per_day
   data%Rdoc_minerl     = Rdoc_minerl/secs_per_day
   data%Rpon_hydrol     = Rpon_hydrol/secs_per_day
   data%Rdon_minerl     = Rdon_minerl/secs_per_day
   data%Rpop_hydrol     = Rpop_hydrol/secs_per_day
   data%Rdop_minerl     = Rdop_minerl/secs_per_day
   data%theta_hydrol    = theta_hydrol
   data%theta_minerl    = theta_minerl
   data%Kpom_hydrol     = Kpom_hydrol
   data%Kdom_minerl     = Kdom_minerl
   !-- Refractory organic matter (optional)
   data%simRPools       = simRPools
   data%Rdocr_miner     = Rdocr_miner/secs_per_day
   data%Rdonr_miner     = Rdonr_miner/secs_per_day
   data%Rdopr_miner     = Rdopr_miner/secs_per_day
   data%Rcpom_bdown     = Rcpom_bdown/secs_per_day
   data%X_cpom_n        = X_cpom_n
   data%X_cpom_p        = X_cpom_p
   !-- Light related parameters
   data%KeDOM           = KeDOM
   data%KePOM           = KePOM
   data%KeDOMR          = KeDOMR
   data%KeCPOM          = KeCPOM
   data%simphotolysis   = simphotolysis
   data%photo_fmin      = photo_fmin
   c = photo_c
   !-- Particle settling parameters
   data%settling        = settling
   data%w_pom           = w_pom/secs_per_day
   data%d_pom           = d_pom
   data%rho_pom         = rho_pom
   data%w_cpom          = w_cpom/secs_per_day
   data%d_cpom          = d_cpom
   data%rho_cpom        = rho_cpom
   IF( settling==0 ) THEN
     w_pom=zero_; w_cpom=zero_
   ENDIF
   !-- Sediment interaction parameters (basic model)
   data%resuspension    = resuspension
   data%sedimentOMfrac  = sedimentOMfrac
   data%Xsc             = Xsc
   data%Xsn             = Xsn
   data%Xsp             = Xsp
   data%Fsed_doc        = Fsed_doc/secs_per_day
   data%Fsed_don        = Fsed_don/secs_per_day
   data%Fsed_dop        = Fsed_dop/secs_per_day
   data%Ksed_dom        = Ksed_dom
   data%theta_sed_dom   = theta_sed_dom
   !-- Misc options
   data%extra_diag      = extra_diag

   ! Register state variables
   data%id_doc = aed2_define_variable('doc','mmol/m**3','dissolved organic carbon', &
                                    doc_initial,minimum=doc_min,maximum=doc_max)
   data%id_poc = aed2_define_variable('poc','mmol/m**3','particulate organic carbon', &
                                    poc_initial,minimum=poc_min,maximum=poc_max,mobility=data%w_pom)
   data%id_don = aed2_define_variable('don','mmol/m**3','dissolved organic nitrogen', &
                                    don_initial,minimum=don_min,maximum=don_max)
   data%id_pon = aed2_define_variable('pon','mmol/m**3','particulate organic nitrogen', &
                                    pon_initial,minimum=pon_min,maximum=pon_max,mobility=data%w_pom)
   data%id_dop = aed2_define_variable('dop','mmol/m**3','dissolved organic phosphorus', &
                                    dop_initial,minimum=dop_min,maximum=dop_max)
   data%id_pop = aed2_define_variable('pop','mmol/m**3','particulate organic phosphorus', &
                                    pop_initial,minimum=pop_min,maximum=pop_max,mobility=data%w_pom)
   IF (simRPools) THEN
     data%id_docr = aed2_define_variable('docr','mmol/m**3','refractory dissolved organic carbon', &
                                    docr_initial,minimum=zero_)
     data%id_donr = aed2_define_variable('donr','mmol/m**3','refractory dissolved organic nitrogen', &
                                    donr_initial,minimum=zero_)
     data%id_dopr = aed2_define_variable('dopr','mmol/m**3','refractory dissolved organic phosphorus', &
                                    dopr_initial,minimum=zero_)
     data%id_cpom = aed2_define_variable('cpom','mmol/m**3','coarse particulate matter', &
                                    cpom_initial,minimum=zero_,mobility=data%w_cpom)
   ENDIF


   ! Register external state variable dependencies
   !-- carbon
   data%use_oxy = doc_miner_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%use_oxy) &
     data%id_oxy = aed2_locate_variable(doc_miner_reactant_variable)
   data%use_dic = doc_miner_product_variable .NE. '' !This means carbon module switched on
   IF (data%use_dic) &
     data%id_dic = aed2_locate_variable(doc_miner_product_variable)
   !-- nitrogen
   data%use_amm = don_miner_product_variable .NE. '' !This means nitrogen module switched on
   IF (data%use_amm) &
     data%id_amm = aed2_locate_variable(don_miner_product_variable)
   !-- phosphorus
   data%use_frp = dop_miner_product_variable .NE. '' !This means phosphorus module switched on
   IF (data%use_frp) &
     data%id_frp = aed2_locate_variable(dop_miner_product_variable)

   !-- sediment flux link variables
   data%id_Fsed_pon = -1 ; data%id_Fsed_pop = -1 ; data%id_Fsed_poc = -1
   data%id_Fsed_don = -1 ; data%id_Fsed_dop = -1 ; data%id_Fsed_doc = -1
   data%use_Fsed_model = Fsed_pon_variable .NE. ''
   IF (data%use_Fsed_model) THEN
     data%id_Fsed_pon = aed2_locate_global_sheet(Fsed_pon_variable)
     IF (Fsed_don_variable .NE. '') &
        data%id_Fsed_don = aed2_locate_global_sheet(Fsed_don_variable)
     IF (Fsed_pop_variable .NE. '') &
        data%id_Fsed_pop = aed2_locate_global_sheet(Fsed_pop_variable)
     IF (Fsed_dop_variable .NE. '') &
        data%id_Fsed_dop = aed2_locate_global_sheet(Fsed_dop_variable)
     IF (Fsed_poc_variable .NE. '') &
        data%id_Fsed_poc = aed2_locate_global_sheet(Fsed_poc_variable)
     IF (Fsed_doc_variable .NE. '') &
        data%id_Fsed_doc = aed2_locate_global_sheet(Fsed_doc_variable)
   ENDIF

   !-- sedimentation link variables
   data%id_Psed_poc = -1 ; data%id_Psed_pon = -1 ; data%id_Psed_pop = -1
   data%use_Psed_model = Psed_poc_variable .NE. ''
   IF (data%use_Psed_model) THEN
     ! locate link variables
     data%id_Psed_poc = aed2_locate_global_sheet(Psed_poc_variable)
     IF (Psed_pon_variable .NE. '') &
        data%id_Psed_pon = aed2_locate_global_sheet(Psed_pon_variable)
     IF (Psed_pop_variable .NE. '') &
        data%id_Psed_pop = aed2_locate_global_sheet(Psed_pop_variable)
   ELSE
     ! make some diagnostics for internal use
     data%id_Psed_poc = aed2_define_sheet_diag_variable('Psed_poc','mmol/m**2/s',  'POC sedimentation')
     data%id_Psed_pon = aed2_define_sheet_diag_variable('Psed_pon','mmol/m**2/s',  'POC sedimentation')
     data%id_Psed_pop = aed2_define_sheet_diag_variable('Psed_pop','mmol/m**2/s',  'POC sedimentation')
     IF (simRPools) &
       data%id_Psed_cpom = aed2_define_sheet_diag_variable('Psed_cpom','mmol/m**2/s',  'CPOM sedimentation')
   ENDIF

   !-- resuspension link variable
   IF ( data%resuspension>0 .AND. .NOT.resus_link .EQ. '' ) THEN
      data%id_l_resus  = aed2_locate_global(TRIM(resus_link)) ! ('TRC_resus')
   ELSE
      data%resuspension = 0
   ENDIF

   !-- light
   IF (simRPools) THEN
     data%id_vis= aed2_locate_global('par')
     data%id_uva= aed2_locate_global('uva')
     data%id_uvb= aed2_locate_global('uvb')
   ENDIF

   ! Register diagnostic variables
   data%id_cdom = aed2_define_diag_variable('CDOM','/m',  'Chromophoric DOM (CDOM)')
   IF (extra_diag) THEN
     data%id_poc_miner = aed2_define_diag_variable('poc_miner','mmol/m**3/d',  'POC mineralisation')
     data%id_doc_miner = aed2_define_diag_variable('doc_miner','mmol/m**3/d',  'DOC mineralisation')
     data%id_sed_poc = aed2_define_sheet_diag_variable('sed_poc','mmol/m**2/d',  'Net POC sediment flux')
     data%id_sed_doc = aed2_define_sheet_diag_variable('sed_doc','mmol/m**2/d',  'DOC sediment flux')

     data%id_pon_miner = aed2_define_diag_variable('pon_miner','mmol/m**3/d',  'PON mineralisation')
     data%id_don_miner = aed2_define_diag_variable('don_miner','mmol/m**3/d',  'DON mineralisation')
     data%id_sed_pon = aed2_define_sheet_diag_variable('sed_pon','mmol/m**2/d',  'Net PON sediment flux')
     data%id_sed_don = aed2_define_sheet_diag_variable('sed_don','mmol/m**2/d',  'DON sediment flux')

     data%id_pop_miner = aed2_define_diag_variable('pop_miner','mmol/m**3/d',  'POP mineralisation')
     data%id_dop_miner = aed2_define_diag_variable('dop_miner','mmol/m**3/d',  'DOP mineralisation')
     data%id_sed_pop = aed2_define_sheet_diag_variable('sed_pop','mmol/m**2/d',  'Net POP sediment flux')
     data%id_sed_dop = aed2_define_sheet_diag_variable('sed_dop','mmol/m**2/d',  'DOP sediment flux')

     data%id_bod = aed2_define_diag_variable('BOD5','mg O2/L',  'Biochemical Oxygen Demand (BOD)')
     IF ( simphotolysis ) data%id_photolysis = &
       aed2_define_diag_variable('photolysis','mmol C/m3/d',  'photolysis rate of breakdown of DOC')
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')
   data%id_rho = aed2_locate_global('density')

END SUBROUTINE aed2_define_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION photo(Q, cdom, band)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: Q          ! [W/m2] incoming light intensity
   AED_REAL,INTENT(in) :: cdom       ! [/m] current cdom absorbance (based on doc)
   INTEGER, INTENT(in) :: band       ! [-] bandwidth identifier (VIS, UVA, UVB)
!LOCALS
!  AED_REAL,PARAMETER :: c = 7.52    ! in NML as photo_c and module global c
   AED_REAL,PARAMETER :: d = 0.0122  ! ref
   AED_REAL,PARAMETER :: S = 0.0188  ! ref
   AED_REAL,PARAMETER :: x = 440.    ! [nm] wavelength at which cdom is measured
   AED_REAL,PARAMETER :: lamda(3) = (/440., 358., 298./) ! [nm] wavelengths
   ! vis = 390-700 nanometers
   ! uva = 315-400 nanometers
   ! uvb = 280-315 nanometers
!
   AED_REAL :: phi, alpha, QQ
!
!-------------------------------------------------------------------------------
!BEGIN
   ! convert incoming light to quantum units: [mol photons /m2 /hr]
   QQ = Q * 4.53 * 1e-6 * 3.6e3

   ! apparent quantum yield [mol produced /(mol absorbed photons) /nm-1]
   phi = c * 10.**(-(d*lamda(band)))

   ! absorption coefficient (/m)
   alpha = cdom*exp(S*(x-lamda(band)))

   ! return the rate at which DOC is broken down [mmol C /m3 /s]
   photo = phi * QQ * alpha * 1e3 / 3.6e3

END FUNCTION photo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_organic_matter(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_organic_matter model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: pon,don,amm,oxy,temp !State variables
   AED_REAL :: pon_hydrolysis, don_mineralisation
   AED_REAL :: pop,dop,frp !State variables
   AED_REAL :: pop_hydrolysis, dop_mineralisation
   AED_REAL :: poc,doc,dic !State variables
   AED_REAL :: poc_hydrolysis, doc_mineralisation
   AED_REAL :: docr,donr,dopr,cpom,cdom
   AED_REAL :: docr_mineralisation, donr_mineralisation
   AED_REAL :: dopr_mineralisation, cpom_breakdown
   AED_REAL :: photolysis
   AED_REAL :: vis, uva, uvb, photo_fmin, cdoc

   AED_REAL, PARAMETER :: Yoxy_doc_miner = 32./12. !ratio of oxygen to carbon utilised during doc mineralisation

!-----------------------------------------------------------------------
!BEGIN
   !CALL log_message('model aed2_organic_matter enter do loop successfully.')

   photolysis = zero_
   photo_fmin = data%photo_fmin

   ! Retrieve current (local) state variable values.
   pon = _STATE_VAR_(data%id_pon)! particulate organic nitrogen
   don = _STATE_VAR_(data%id_don)! dissolved organic nitrogen
   pop = _STATE_VAR_(data%id_pop)! particulate organic phosphorus
   dop = _STATE_VAR_(data%id_dop)! dissolved organic phosphorus
   poc = _STATE_VAR_(data%id_poc)! particulate organic carbon
   doc = _STATE_VAR_(data%id_doc)! dissolved organic carbon

   docr = zero_
   IF(data%simRPools) THEN
      docr = _STATE_VAR_(data%id_docr)
      donr = _STATE_VAR_(data%id_donr)
      dopr = _STATE_VAR_(data%id_dopr)
      cpom = _STATE_VAR_(data%id_cpom)
      vis = _STATE_VAR_(data%id_vis)
      uva = _STATE_VAR_(data%id_uva)
      uvb = _STATE_VAR_(data%id_uvb)
      cdom = _DIAG_VAR_(data%id_cdom)
   ENDIF


   IF (data%use_oxy) THEN ! & use_oxy
      oxy = _STATE_VAR_(data%id_oxy)! oxygen
   ELSE
      oxy = 300.0
   ENDIF
   IF (data%use_dic) THEN ! & use_dic
      dic = _STATE_VAR_(data%id_dic)! disolved inorganic carbon
   ELSE
      dic = 100.0
   ENDIF
   IF (data%use_amm) THEN ! & use_amm
      amm = _STATE_VAR_(data%id_amm)! ammonium
   ELSE
      amm = 0.0
   ENDIF
   IF (data%use_frp) THEN ! & use_frp
      frp = _STATE_VAR_(data%id_frp)! phosphate
   ELSE
      frp = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_temp) ! temperature

   ! Define some intermediate quantities (units mmol N/m3/day)
   pon_hydrolysis = fpon_miner(data%use_oxy,data%Rpon_hydrol,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)
   don_mineralisation = fdon_miner(data%use_oxy,data%Rdon_minerl,data%Kdom_minerl,data%theta_minerl,oxy,temp)
   pop_hydrolysis = fpop_miner(data%use_oxy,data%Rpop_hydrol,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)
   dop_mineralisation = fdop_miner(data%use_oxy,data%Rdop_minerl,data%Kdom_minerl,data%theta_minerl,oxy,temp)
   poc_hydrolysis = fpoc_miner(data%use_oxy,data%Rpoc_hydrol,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)
   doc_mineralisation = fdoc_miner(data%use_oxy,data%Rdoc_minerl,data%Kdom_minerl,data%theta_minerl,oxy,temp)
   IF(data%simRPools) THEN
      docr_mineralisation = fdoc_miner(data%use_oxy,data%Rdocr_miner,data%Kdom_minerl,data%theta_minerl,oxy,temp)
      donr_mineralisation = fdoc_miner(data%use_oxy,data%Rdonr_miner,data%Kdom_minerl,data%theta_minerl,oxy,temp)
      dopr_mineralisation = fdoc_miner(data%use_oxy,data%Rdopr_miner,data%Kdom_minerl,data%theta_minerl,oxy,temp)
      cpom_breakdown = fpoc_miner(data%use_oxy,data%Rcpom_bdown,data%Kpom_hydrol,data%theta_hydrol,oxy,temp)
      IF ( data%simphotolysis ) THEN
         photolysis = photo(vis,cdom,1) + photo(uva,cdom,2) + photo(uvb,cdom,3)
         !# Limit photolysis to 90% of doc pool within 1 hour
         IF(photolysis > 0.9*docr/3.6e3) photolysis = 0.9*docr/3.6e3
      ENDIF
   ENDIF

   ! Set temporal derivatives
   _FLUX_VAR_(data%id_pon)  = _FLUX_VAR_(data%id_pon) + (-pon*pon_hydrolysis)
   _FLUX_VAR_(data%id_don)  = _FLUX_VAR_(data%id_don) + (pon*pon_hydrolysis-don*don_mineralisation)
   _FLUX_VAR_(data%id_pop)  = _FLUX_VAR_(data%id_pop) + (-pop*pop_hydrolysis)
   _FLUX_VAR_(data%id_dop)  = _FLUX_VAR_(data%id_dop) + (pop*pop_hydrolysis-dop*dop_mineralisation)
   _FLUX_VAR_(data%id_poc)  = _FLUX_VAR_(data%id_poc) + (-poc*poc_hydrolysis)
   _FLUX_VAR_(data%id_doc)  = _FLUX_VAR_(data%id_doc) + (poc*poc_hydrolysis-doc*doc_mineralisation)
   IF(data%simRPools) THEN
      docr = MAX(1.e-3,docr)
      _FLUX_VAR_(data%id_docr) = _FLUX_VAR_(data%id_docr) - docr*docr_mineralisation - photolysis
      _FLUX_VAR_(data%id_doc)  = _FLUX_VAR_(data%id_doc)  + docr*docr_mineralisation + photolysis*photo_fmin
      _FLUX_VAR_(data%id_donr) = _FLUX_VAR_(data%id_donr) - donr*donr_mineralisation - photolysis*(donr/docr)
      _FLUX_VAR_(data%id_don)  = _FLUX_VAR_(data%id_don)  + donr*donr_mineralisation + photolysis*photo_fmin*(donr/docr)
      _FLUX_VAR_(data%id_dopr) = _FLUX_VAR_(data%id_dopr) - dopr*dopr_mineralisation - photolysis*(dopr/docr)
      _FLUX_VAR_(data%id_dop)  = _FLUX_VAR_(data%id_dop)  + dopr*dopr_mineralisation + photolysis*photo_fmin*(dopr/docr)
      _FLUX_VAR_(data%id_cpom) = _FLUX_VAR_(data%id_cpom) - cpom*cpom_breakdown
      _FLUX_VAR_(data%id_poc)  = _FLUX_VAR_(data%id_poc) + (cpom*cpom_breakdown)
      _FLUX_VAR_(data%id_pon)  = _FLUX_VAR_(data%id_pon) + (cpom*cpom_breakdown)*data%X_cpom_n
      _FLUX_VAR_(data%id_pop)  = _FLUX_VAR_(data%id_pop) + (cpom*cpom_breakdown)*data%X_cpom_p
   ENDIF


   ! If an externally maintained oxygen pool is present, take mineralisation from it
   IF (data%use_oxy) THEN
      !_FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-Yoxy_doc_miner*doc*doc_mineralisation)
      _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (-doc*doc_mineralisation)
   ENDIF
   if (data%use_dic) THEN
      _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (doc*doc_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + photolysis*(1.-photo_fmin)
   ENDIF
   IF (data%use_amm) THEN
      _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + (don*don_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + photolysis*(1.-photo_fmin)*(donr/docr)
    ! IF ( data%simRPools ) _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + photolysis*(1.-photo_fmin)*MAX((donr/docr),1e-3)
   ENDIF
   IF (data%use_frp) THEN
      _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + (dop*dop_mineralisation)
      IF ( data%simRPools ) _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + photolysis*(1.-photo_fmin)*(dopr/docr)
   ENDIF

   !# Export diagnostic variables
   !-- CDOM computed as a function of DOC amount, as empirically defined by
   !   Kostoglidis et al 2005 for Swan-Canning
   IF ( data%simRPools ) THEN
     cdoc = MIN(doc+docr,1.e4)
   ELSE
     cdoc = MIN(doc,1.e4)
   ENDIF
  _DIAG_VAR_(data%id_cdom) = 0.35*exp(0.1922*(cdoc)*(12./1e3))
   !-- Extra diagnostics
   IF (data%extra_diag) THEN
     _DIAG_VAR_(data%id_poc_miner) = poc_hydrolysis*poc*secs_per_day
     _DIAG_VAR_(data%id_pon_miner) = pon_hydrolysis*pon*secs_per_day
     _DIAG_VAR_(data%id_pop_miner) = pop_hydrolysis*pop*secs_per_day
     _DIAG_VAR_(data%id_doc_miner) = doc_mineralisation*doc*secs_per_day
     _DIAG_VAR_(data%id_don_miner) = don_mineralisation*don*secs_per_day
     _DIAG_VAR_(data%id_dop_miner) = dop_mineralisation*dop*secs_per_day
     IF ( data%simphotolysis ) &
        _DIAG_VAR_(data%id_photolysis) = photolysis*secs_per_day
     ! BOD5 is computed as the amount of oxygen consumed over a 5-day period (mgO2/L)
     _DIAG_VAR_(data%id_bod) = Yoxy_doc_miner*(doc*doc_mineralisation*secs_per_day*12./1e3)*5.0
   ENDIF

END SUBROUTINE aed2_calculate_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_organic_matter(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, oxy

   ! State
   AED_REAL :: pon,don
   AED_REAL :: pop,dop
   AED_REAL :: poc,doc

   ! Temporary variables
   AED_REAL :: poc_flux,doc_flux
   AED_REAL :: pon_flux,don_flux
   AED_REAL :: pop_flux,dop_flux
   AED_REAL :: Fsed_poc,Fsed_doc
   AED_REAL :: Fsed_pon,Fsed_don
   AED_REAL :: Fsed_pop,Fsed_dop
   AED_REAL :: Psed_poc, Psed_pon, Psed_pop
   AED_REAL :: fT, fDO, fDOM

!-------------------------------------------------------------------------------
!BEGIN


   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp)                   ! local temperature
   oxy = 300. ; IF (data%use_oxy) oxy = _STATE_VAR_(data%id_oxy)

   ! Retrieve current (local) state variable values.
   poc = _STATE_VAR_(data%id_poc) ! particulate organic carbon
   doc = _STATE_VAR_(data%id_doc) ! dissolved organic carbon
   pon = _STATE_VAR_(data%id_pon) ! particulate organic nitrogen
   don = _STATE_VAR_(data%id_don) ! dissolved organic nitrogen
   pop = _STATE_VAR_(data%id_pop) ! particulate organic phosphorus
   dop = _STATE_VAR_(data%id_dop) ! dissolved organic phosphorus

   ! Set flux rate of organic matter pools
   IF (data%use_Fsed_model) THEN
      ! Get the flux value from the linked variable in aed_sedflux
      Fsed_pon = _STATE_VAR_S_(data%id_Fsed_pon)
      IF ( data%id_Fsed_don .NE. -1 ) Fsed_don = _STATE_VAR_S_(data%id_Fsed_don)
      IF ( data%id_Fsed_pop .NE. -1 ) Fsed_pop = _STATE_VAR_S_(data%id_Fsed_pop)
      IF ( data%id_Fsed_dop .NE. -1 ) Fsed_dop = _STATE_VAR_S_(data%id_Fsed_dop)
      IF ( data%id_Fsed_poc .NE. -1 ) Fsed_poc = _STATE_VAR_S_(data%id_Fsed_poc)
      IF ( data%id_Fsed_doc .NE. -1 ) Fsed_doc = _STATE_VAR_S_(data%id_Fsed_doc)
   ELSE
      ! Compute directly
      !-- DOM fluxes
      fT = data%theta_sed_dom**(temp-20.0)
      fDO = data%Ksed_dom/(data%Ksed_dom+oxy)
      fDOM = 1.
      Fsed_doc = data%Fsed_doc * fDO * fT * fDOM
      Fsed_don = data%Fsed_don * fDO * fT * fDOM
      Fsed_dop = data%Fsed_dop * fDO * fT * fDOM
      !-- POM fluxes
      Fsed_poc = zero_ !data%Fsed_pom * sedimentOMfrac * data%Xsc *(bottom_stress - data%tau_0) / data%tau_r
      Fsed_pon = zero_ !data%Fsed_pom * sedimentOMfrac * data%Xsn *(bottom_stress - data%tau_0) / data%tau_r
      Fsed_pop = zero_ !data%Fsed_pom * sedimentOMfrac * data%Xsp *(bottom_stress - data%tau_0) / data%tau_r
      IF( data%resuspension>0 )THEN
        Fsed_poc = _DIAG_VAR_S_(data%id_l_resus) * data%sedimentOMfrac * data%Xsc
        Fsed_pon = _DIAG_VAR_S_(data%id_l_resus) * data%sedimentOMfrac * data%Xsn
        Fsed_pop = _DIAG_VAR_S_(data%id_l_resus) * data%sedimentOMfrac * data%Xsp
      ENDIF
   ENDIF

   ! Set bottom fluxes for the pelagic variables, note this doesnt consdier
   ! sedimentation of particles (change per surface area per second)
   _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + Fsed_doc
   _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + Fsed_don
   _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) + Fsed_dop
   _FLUX_VAR_(data%id_poc) = _FLUX_VAR_(data%id_poc) + Fsed_poc
   _FLUX_VAR_(data%id_pon) = _FLUX_VAR_(data%id_pon) + Fsed_pon
   _FLUX_VAR_(data%id_pop) = _FLUX_VAR_(data%id_pop) + Fsed_pop

   ! Get sedimentation flux (mmmol/m2/s) loss into the benthos
   IF (data%use_Psed_model) THEN
       !-- linked variable requires populating
       Psed_poc = data%w_pom * max(zero_,poc)
       IF ( data%id_Psed_pon .NE. -1 ) Psed_pon = data%w_pom * max(zero_,pon)
       IF ( data%id_Psed_pop .NE. -1 ) Psed_pop = data%w_pom * max(zero_,pop)
       !?
   ELSE
       !-- diagnostic was set in mobility
       Psed_poc = _DIAG_VAR_S_(data%id_Psed_poc)
       Psed_pon = _DIAG_VAR_S_(data%id_Psed_pon)
       Psed_pop = _DIAG_VAR_S_(data%id_Psed_pop)
   ENDIF

   ! Set sink & source terms for the sediment pools (change per surface area per second)
   ! Note that this must include fluxes to and from the pelagic.
   IF( data%simSedimentOM ) THEN
    !_FLUX_VAR_B_(data%id_sed_poc) = _FLUX_VAR_B_(data%id_sed_poc) - poc_flux
    !_FLUX_VAR_B_(data%id_sed_pon) = _FLUX_VAR_B_(data%id_sed_pon) - pon_flux
    !_FLUX_VAR_B_(data%id_sed_pop) = _FLUX_VAR_B_(data%id_sed_pop) - pop_flux
    !_FLUX_VAR_B_(data%id_sed_doc) = _FLUX_VAR_B_(data%id_sed_doc) - doc_flux
    !_FLUX_VAR_B_(data%id_sed_don) = _FLUX_VAR_B_(data%id_sed_don) - don_flux
    !_FLUX_VAR_B_(data%id_sed_dop) = _FLUX_VAR_B_(data%id_sed_dop) - dop_flux
  ENDIF

   ! Also store net sediment fluxes as diagnostic variable.
   IF (data%extra_diag) THEN
      _DIAG_VAR_S_(data%id_sed_poc) = Fsed_poc + Psed_poc ! resus & settling
      _DIAG_VAR_S_(data%id_sed_doc) = Fsed_doc            ! dissolved flux
      _DIAG_VAR_S_(data%id_sed_pon) = Fsed_poc + Psed_poc
      _DIAG_VAR_S_(data%id_sed_don) = Fsed_don
      _DIAG_VAR_S_(data%id_sed_pop) = Fsed_poc + Psed_poc
      _DIAG_VAR_S_(data%id_sed_dop) = Fsed_dop
   ENDIF
END SUBROUTINE aed2_calculate_benthic_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_organic_matter(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: doc, poc, cpom, cdom
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
   doc = _STATE_VAR_(data%id_doc)
   poc = _STATE_VAR_(data%id_poc)

   ! Self-shading with explicit contribution from background OM concentration.
   extinction = extinction + (data%KeDOM*doc +data%KePOM*poc)

   IF (data%simRPools) THEN
     cdom = _DIAG_VAR_(data%id_cdom)    ! CDOM is in "/m" units here
     cpom = _STATE_VAR_(data%id_cpom)   ! CPOM is in "mmol/m3" units here

     extinction = extinction + (data%KeDOMR*cdom)
     extinction = extinction + (data%KeCPOM*cpom)
   ENDIF
END SUBROUTINE aed2_light_extinction_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_mobility_organic_matter(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_organic_matter_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: vvel, vvel_cpom
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: rho_pom, temp
!
!-------------------------------------------------------------------------------
!BEGIN
     ! settling = 0 : no settling
     ! settling = 1 : constant settling @ w_pom
     ! settling = 2 : constant settling @ w_pom, corrected for variable density
     ! settling = 3 : settling based on Stoke's Law (calculated below)

     SELECT CASE (data%settling)

        CASE ( _MOB_OFF_ )
          ! disable settling by settign vertical velocity to 0
          vvel = zero_
          vvel_cpom = zero_

        CASE ( _MOB_CONST_ )
          ! constant settling velocity using user provided value
          vvel = data%w_pom
          vvel_cpom = data%w_cpom

        CASE ( _MOB_TEMP_ )
          ! constant settling velocity @20C corrected for density changes
          pw = _STATE_VAR_(data%id_rho)
          temp = _STATE_VAR_(data%id_temp)
          mu = water_viscosity(temp)
          mu20 = 0.001002  ! N s/m2
          pw20 = 998.2000  ! kg/m3 (assuming freshwater)
          vvel = data%w_pom*mu20*pw / ( mu*pw20 )
          vvel_cpom = data%w_cpom*mu20*pw / ( mu*pw20 )

        CASE ( _MOB_STOKES_ )
          ! settling velocity based on Stokes Law calculation and cell density
          pw = _STATE_VAR_(data%id_rho)              ! water density
          temp = _STATE_VAR_(data%id_temp)
          mu = water_viscosity(temp)                 ! water dynamic viscosity
          rho_pom = data%rho_pom
          vvel = -9.807*(data%d_pom**2.)*( rho_pom-pw ) / ( 18.*mu )
          IF(data%simRPools) &
          vvel_cpom = -9.807*(data%d_cpom**2.)*( data%rho_cpom-pw ) / ( 18.*mu )

        CASE DEFAULT
          ! unknown settling/migration option selection
          vvel = data%w_pom
          vvel_cpom = data%w_cpom

      END SELECT
      ! set global mobility array (m/s), later used to compute settling
      mobility(data%id_poc) = vvel
      mobility(data%id_pon) = vvel
      mobility(data%id_pop) = vvel
      IF(data%simRPools) mobility(data%id_cpom) = vvel_cpom
      ! set sedimentation flux (mmmol/m2) for later use/reporting
      _DIAG_VAR_S_(data%id_Psed_poc) = mobility(data%id_poc)*_STATE_VAR_(data%id_poc)
      _DIAG_VAR_S_(data%id_Psed_pon) = mobility(data%id_pon)*_STATE_VAR_(data%id_pon)
      _DIAG_VAR_S_(data%id_Psed_pop) = mobility(data%id_pop)*_STATE_VAR_(data%id_pop)
      IF(data%simRPools) _DIAG_VAR_S_(data%id_Psed_cpom) = mobility(data%id_cpom)*_STATE_VAR_(data%id_cpom)

END SUBROUTINE aed2_mobility_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fpon_miner(use_oxy,Rpon_miner,Kpon_miner,theta_pon_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Nitrogen
!
! Michaelis-Menten formulation for mineralisation
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpon_miner,Kpon_miner,theta_pon_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpon_miner = Rpon_miner * oxy/(Kpon_miner+oxy) * (theta_pon_miner**(temp-20.0))
   ELSE
      fpon_miner = Rpon_miner * (theta_pon_miner**(temp-20.0))
   ENDIF

END FUNCTION fpon_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdon_miner(use_oxy,Rdon_miner,Kdon_miner,theta_don_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdon_miner,Kdon_miner,theta_don_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdon_miner = Rdon_miner * oxy/(Kdon_miner+oxy) * (theta_don_miner**(temp-20.0))
   ELSE
      fdon_miner = Rdon_miner * (theta_don_miner**(temp-20.0))
   ENDIF

END FUNCTION fdon_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fpop_miner(use_oxy,Rpop_miner,Kpop_miner,theta_pop_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Phosphorus
!
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpop_miner,Kpop_miner,theta_pop_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpop_miner = Rpop_miner * oxy/(Kpop_miner+oxy) * (theta_pop_miner**(temp-20.0))
   ELSE
      fpop_miner = Rpop_miner * (theta_pop_miner**(temp-20.0))
   ENDIF

END FUNCTION fpop_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdop_miner(use_oxy,Rdop_miner,Kdop_miner,theta_dop_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdop_miner,Kdop_miner,theta_dop_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdop_miner = Rdop_miner * oxy/(Kdop_miner+oxy) * (theta_dop_miner**(temp-20.0))
   ELSE
      fdop_miner = Rdop_miner * (theta_dop_miner**(temp-20.0))
   ENDIF

END FUNCTION fdop_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fpoc_miner(use_oxy,Rpoc_miner,Kpoc_miner,theta_poc_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Carbon
!
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpoc_miner,Kpoc_miner,theta_poc_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpoc_miner = Rpoc_miner * oxy/(Kpoc_miner+oxy) * (theta_poc_miner**(temp-20.0))
   ELSE
      fpoc_miner = Rpoc_miner * (theta_poc_miner**(temp-20.0))
   ENDIF

END FUNCTION fpoc_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdoc_miner(use_oxy,Rdoc_miner,Kdoc_miner,theta_doc_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdoc_miner,Kdoc_miner,theta_doc_miner,oxy,temp
!
!-----------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdoc_miner = Rdoc_miner * oxy/(Kdoc_miner+oxy) * (theta_doc_miner**(temp-20.0))
   ELSE
      fdoc_miner = Rdoc_miner * (theta_doc_miner**(temp-20.0))
   ENDIF

END FUNCTION fdoc_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_organic_matter
