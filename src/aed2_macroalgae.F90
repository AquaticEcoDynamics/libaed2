!###############################################################################
!#                                                                             #
!# aed2_macroalgae.F90                                                      #
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
!# Created August 2011                                                         #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_macroalgae
!-------------------------------------------------------------------------------
!  aed2_macroalgae --- macroalgae biogeochemical model
!-------------------------------------------------------------------------------
   USE aed2_core
   USE aed2_util,ONLY : find_free_lun, &
                        exp_integral, &
                        aed2_bio_temp_function, &
                        fTemp_function,fSal_function, &
                        water_viscosity, in_zone_set
   USE aed2_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private

   PUBLIC aed2_macroalgae_data_t


   TYPE,extends(aed2_model_data_t) :: aed2_macroalgae_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_p(:), id_pben(:)
      INTEGER,ALLOCATABLE :: id_in(:), id_inben(:)
      INTEGER,ALLOCATABLE :: id_ip(:), id_ipben(:)
      INTEGER,ALLOCATABLE :: id_rho(:)
      INTEGER,ALLOCATABLE :: id_NtoP(:)
      INTEGER,ALLOCATABLE :: id_vvel(:)
      INTEGER,ALLOCATABLE :: id_fT(:), id_fI(:), id_fNit(:), &
                             id_fPho(:), id_fSil(:), id_fSal(:)
      INTEGER,ALLOCATABLE :: id_fT_ben(:), id_fI_ben(:), id_fNit_ben(:), &
                             id_fPho_ben(:), id_fSal_ben(:)
      INTEGER :: id_Pexctarget,id_Pmorttarget,id_Pupttarget(1:2)
      INTEGER :: id_Nexctarget,id_Nmorttarget,id_Nupttarget(1:4)
      INTEGER :: id_Cexctarget,id_Cmorttarget,id_Cupttarget
      INTEGER :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      INTEGER :: id_DOupttarget
      INTEGER :: id_par, id_I_0, id_extc, id_taub, id_sedzone
      INTEGER :: id_tem, id_sal, id_dz, id_dens
      INTEGER :: id_GPP, id_NCP, id_PPR, id_NPR, id_dPAR
      INTEGER :: id_TMALG, id_TCHLA, id_TIN, id_TIP, id_MPB, id_d_MPB, id_d_BPP
      INTEGER :: id_NUP, id_PUP, id_CUP
      INTEGER :: id_mhsi
      INTEGER :: n_zones
      AED_REAL, ALLOCATABLE :: active_zones(:)

      !# Model parameters
      INTEGER  :: num_malgae
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: phytos
      ! LOGICAL  :: do_exc,do_mort,do_upt, do_N2uptake
      LOGICAL  :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL  :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL  :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL  :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER  :: nnup, npup
      AED_REAL :: dic_per_n
      AED_REAL :: min_rho,max_rho
      AED_REAL,ALLOCATABLE :: resuspension(:)
      AED_REAL :: slough_stress
      INTEGER :: simMalgHSI

     CONTAINS
         PROCEDURE :: define            => aed2_define_macroalgae
         PROCEDURE :: calculate         => aed2_calculate_macroalgae
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_macroalgae
         PROCEDURE :: mobility          => aed2_mobility_macroalgae
         PROCEDURE :: light_extinction  => aed2_light_extinction_macroalgae
!        PROCEDURE :: delete            => aed2_delete_macroalgae

   END TYPE

   AED_REAL :: dtlim = 0.9 * 3600
   LOGICAL  :: extra_diag = .false.

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed2_macroalgae_load_params(data, dbase, count, list, settling, resuspension)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macroalgae_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: settling(*)
   AED_REAL,INTENT(in)         :: resuspension(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(phyto_nml_data) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /malgae_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD', iostat=status)
    IF (status /= 0) STOP 'Cannot open phyto_data namelist file'
    read(tfil,nml=malgae_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist phyto_data'

    data%num_malgae = count
    ALLOCATE(data%phytos(count))
    ALLOCATE(data%id_p(count)) ; data%id_p(:) = 0
    ALLOCATE(data%id_in(count)) ; data%id_in(:) = 0
    ALLOCATE(data%id_ip(count)) ; data%id_ip(:) = 0
    ALLOCATE(data%id_rho(count)) ; data%id_rho(:) = 0
    ALLOCATE(data%id_NtoP(count)) ; data%id_NtoP(:) = 0
    ALLOCATE(data%id_pben(count)) ; data%id_p(:) = 0
    ALLOCATE(data%id_inben(count)) ; data%id_inben(:) = 0
    ALLOCATE(data%id_ipben(count)) ; data%id_ipben(:) = 0
    IF (extra_diag) THEN
       ALLOCATE(data%id_fT(count)) ; data%id_fT(:) = 0
       ALLOCATE(data%id_fI(count)) ; data%id_fI(:) = 0
       ALLOCATE(data%id_fNit(count)) ; data%id_fNit(:) = 0
       ALLOCATE(data%id_fPho(count)) ; data%id_fPho(:) = 0
       ALLOCATE(data%id_fSil(count)) ; data%id_fSil(:) = 0
       ALLOCATE(data%id_fSal(count)) ; data%id_fSal(:) = 0
       ALLOCATE(data%id_vvel(count)) ; data%id_vvel(:) = 0
       ALLOCATE(data%id_fT_ben(count)) ; data%id_fT_ben(:) = 0
       ALLOCATE(data%id_fI_ben(count)) ; data%id_fI_ben(:) = 0
       ALLOCATE(data%id_fNit_ben(count)) ; data%id_fNit_ben(:) = 0
       ALLOCATE(data%id_fPho_ben(count)) ; data%id_fPho_ben(:) = 0
       ALLOCATE(data%id_fSal_ben(count)) ; data%id_fSal_ben(:) = 0
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%phytos(i)%p_name       = pd(list(i))%p_name
       data%phytos(i)%p0           = pd(list(i))%p0
       data%phytos(i)%w_p          = pd(list(i))%w_p/secs_per_day
       data%phytos(i)%settling      = settling(i)
       data%phytos(i)%resuspension = resuspension(i)
       data%phytos(i)%Xcc          = pd(list(i))%Xcc
       data%phytos(i)%R_growth     = pd(list(i))%R_growth/secs_per_day
       data%phytos(i)%fT_Method    = pd(list(i))%fT_Method
       data%phytos(i)%theta_growth = pd(list(i))%theta_growth
       data%phytos(i)%T_std        = pd(list(i))%T_std
       data%phytos(i)%T_opt        = pd(list(i))%T_opt
       data%phytos(i)%T_max        = pd(list(i))%T_max
       data%phytos(i)%lightModel   = pd(list(i))%lightModel
       data%phytos(i)%I_K          = pd(list(i))%I_K
       data%phytos(i)%I_S          = pd(list(i))%I_S
       data%phytos(i)%KePHY        = pd(list(i))%KePHY
       data%phytos(i)%f_pr         = pd(list(i))%f_pr
       data%phytos(i)%R_resp       = pd(list(i))%R_resp/secs_per_day
       data%phytos(i)%theta_resp   = pd(list(i))%theta_resp
       data%phytos(i)%k_fres       = pd(list(i))%k_fres
       data%phytos(i)%k_fdom       = pd(list(i))%k_fdom
       data%phytos(i)%salTol       = pd(list(i))%salTol
       data%phytos(i)%S_bep        = pd(list(i))%S_bep
       data%phytos(i)%S_maxsp      = pd(list(i))%S_maxsp
       data%phytos(i)%S_opt        = pd(list(i))%S_opt
       data%phytos(i)%simDINUptake = pd(list(i))%simDINUptake
       data%phytos(i)%simDONUptake = pd(list(i))%simDONUptake
       data%phytos(i)%simNFixation = pd(list(i))%simNFixation
       data%phytos(i)%simINDynamics= pd(list(i))%simINDynamics
       data%phytos(i)%N_o          = pd(list(i))%N_o
       data%phytos(i)%K_N          = pd(list(i))%K_N
       data%phytos(i)%X_ncon       = pd(list(i))%X_ncon
       data%phytos(i)%X_nmin       = pd(list(i))%X_nmin
       data%phytos(i)%X_nmax       = pd(list(i))%X_nmax
       data%phytos(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_per_day
       data%phytos(i)%k_nfix       = pd(list(i))%k_nfix
       data%phytos(i)%R_nfix       = pd(list(i))%R_nfix/secs_per_day
       data%phytos(i)%simDIPUptake = pd(list(i))%simDIPUptake
       data%phytos(i)%simIPDynamics= pd(list(i))%simIPDynamics
       data%phytos(i)%P_0          = pd(list(i))%P_0
       data%phytos(i)%K_P          = pd(list(i))%K_P
       data%phytos(i)%X_pcon       = pd(list(i))%X_pcon
       data%phytos(i)%X_pmin       = pd(list(i))%X_pmin
       data%phytos(i)%X_pmax       = pd(list(i))%X_pmax
       data%phytos(i)%R_puptake    = pd(list(i))%R_puptake/secs_per_day
       data%phytos(i)%simSiUptake  = pd(list(i))%simSiUptake
       data%phytos(i)%Si_0         = pd(list(i))%Si_0
       data%phytos(i)%K_Si         = pd(list(i))%K_Si
       data%phytos(i)%X_sicon      = pd(list(i))%X_sicon

       data%phytos(i)%c1           = 0.0124/60.   ! From Chung et al (2014)
       data%phytos(i)%c3           = 0.0230/60.   !  "
       data%phytos(i)%f1           = 0.675        ! Ross and Sharples (2007)
       data%phytos(i)%f2           = 0.750        !  "
       data%phytos(i)%d_phy        = 1e-5

       ! Register group as a state variable
       data%id_p(i) = aed2_define_variable(                                    &
                              TRIM(data%phytos(i)%p_name),                     &
                              'mmol/m**3',                                     &
                              'macroalgae '//TRIM(data%phytos(i)%p_name),   &
                              pd(list(i))%p_initial,                           &
                              minimum=pd(list(i))%p0,                          &
                              mobility = data%phytos(i)%w_p)

       ! Register rho (density) group as a state variable, if required
         IF (data%phytos(i)%settling == _MOB_STOKES_) THEN
           data%id_rho(i) = aed2_define_variable(                               &
                              TRIM(data%phytos(i)%p_name)//'_rho',             &
                              'kg/m**3',                                       &
                           'macroalgae '//TRIM(data%phytos(i)%p_name)//'_rho', &
                              (data%min_rho+data%max_rho)/2.,                  &
                              minimum=data%min_rho,                            &
                              mobility = data%phytos(i)%w_p)
         ENDIF

       ! Register internal nitrogen group as a state variable, if required
       IF (data%phytos(i)%simINDynamics /= 0) THEN
          IF(data%phytos(i)%simINDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_ncon
          ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_nmin
          ENDIF
          ! Register IN group as a state variable
          data%id_in(i) = aed2_define_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IN',              &
                              'mmol/m**3',                                     &
                         'macroalgae '//TRIM(data%phytos(i)%p_name)//'_IN', &
                              pd(list(i))%p_initial*data%phytos(i)%X_ncon,     &
                              minimum=minNut,                                  &
                              mobility = data%phytos(i)%w_p)

       ENDIF

       ! Register internal phosphorus group as a state variable, if required
       IF (data%phytos(i)%simIPDynamics /= 0) THEN
          IF(data%phytos(i)%simIPDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_pcon
          ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_pmin
          ENDIF
          ! Register IP group as a state variable
          data%id_ip(i) = aed2_define_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IP',              &
                              'mmol/m**3',                                     &
                         'macroalgae '//TRIM(data%phytos(i)%p_name)//'_IP', &
                              pd(list(i))%p_initial*data%phytos(i)%X_pcon,     &
                              minimum=minNut,                                  &
                              mobility = data%phytos(i)%w_p)
       ENDIF

       ! Group specific diagnostic variables
       data%id_NtoP(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_NtoP','-', 'internal n:p ratio')
       IF (extra_diag) THEN
          data%id_fI(i)   = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fI', '-', 'fI (0-1)')
          data%id_fNit(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fNit', '-', 'fNit (0-1)')
          data%id_fPho(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fPho', '-', 'fPho (0-1)')
          data%id_fSil(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSil', '-', 'fSil (0-1)')
          data%id_fT(i)   = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fT', '-', 'fT (>0)')
          data%id_fSal(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSal', '-', 'fSal (>1)')
          ! Register vertical velocity diagnostic, where relevant
          IF (data%phytos(i)%settling == _MOB_STOKES_ ) THEN
            data%id_vvel(i) = aed2_define_diag_variable( TRIM(data%phytos(i)%p_name)//'_vvel', 'm/s', 'vertical velocity')

          ENDIF
       ENDIF

       ! Allocate benthic variables if the selcted group has attached/benthic growth form
       IF (data%phytos(i)%settling == _MOB_ATTACHED_) THEN
         data%id_pben(i) = aed2_define_sheet_variable(                       &
                                TRIM(data%phytos(i)%p_name)//'_ben',                 &
                                'mmolC/m**2',                                &
                                'macroalgae '//TRIM(data%phytos(i)%p_name),  &
                                pd(list(i))%p_initial,                       &
                                minimum=pd(list(i))%p0,                   &
                                maximum=8e3 )

         IF (data%phytos(i)%simIPDynamics /= 0) THEN
            IF(data%phytos(i)%simIPDynamics == 1)THEN
              minNut = data%phytos(i)%p0*data%phytos(i)%X_pcon
            ELSE
              minNut = data%phytos(i)%p0*data%phytos(i)%X_pmin
            ENDIF
            ! Register IP group as a state variable
            data%id_ipben(i) = aed2_define_sheet_variable(                     &
                              TRIM(data%phytos(i)%p_name)//'_IP_ben',          &
                              'mmol/m**2',                                     &
                              'macroalgae '//TRIM(data%phytos(i)%p_name)//'_IP_ben', &
                              pd(list(i))%p_initial*data%phytos(i)%X_pcon,     &
                              minimum=minNut)
         ENDIF
         IF (data%phytos(i)%simINDynamics /= 0) THEN
           IF(data%phytos(i)%simINDynamics == 1)THEN
            minNut = data%phytos(i)%p0*data%phytos(i)%X_ncon
           ELSE
            minNut = data%phytos(i)%p0*data%phytos(i)%X_nmin
           ENDIF
           ! Register IN group as a state variable
           data%id_inben(i) = aed2_define_sheet_variable(                     &
                            TRIM(data%phytos(i)%p_name)//'_IN_ben',          &
                            'mmol/m**2',                                     &
                            'macroalgae '//TRIM(data%phytos(i)%p_name)//'_IN_ben', &
                            pd(list(i))%p_initial*data%phytos(i)%X_ncon,     &
                            minimum=minNut)
         ENDIF
         ! Group specific diagnostic variables
         IF (extra_diag) THEN
          data%id_fI_ben(i)   = aed2_define_sheet_diag_variable( TRIM(data%phytos(i)%p_name)//'_fI_ben', '-', 'fI (0-1)')
          data%id_fNit_ben(i) = aed2_define_sheet_diag_variable( TRIM(data%phytos(i)%p_name)//'_fNit_ben', '-', 'fNit (0-1)')
          data%id_fPho_ben(i) = aed2_define_sheet_diag_variable( TRIM(data%phytos(i)%p_name)//'_fPho_ben', '-', 'fPho (0-1)')
          data%id_fT_ben(i)   = aed2_define_sheet_diag_variable( TRIM(data%phytos(i)%p_name)//'_fT_ben', '-', 'fT (>0)')
          data%id_fSal_ben(i) = aed2_define_sheet_diag_variable( TRIM(data%phytos(i)%p_name)//'_fSal_ben', '-', 'fSal (-)')
         ENDIF
       ENDIF



    ENDDO
END SUBROUTINE aed2_macroalgae_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_macroalgae(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the macroalgae biogeochemical model
!
!  Here, the aed2_p_m namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macroalgae_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: num_malgae
   INTEGER  :: the_malgae(MAX_PHYTO_TYPES)
   INTEGER  :: settling(MAX_PHYTO_TYPES)
   AED_REAL :: resuspension(MAX_PHYTO_TYPES)
   CHARACTER(len=64)  :: p_excretion_target_variable=''
   CHARACTER(len=64)  :: p_mortality_target_variable=''
   CHARACTER(len=64)  :: p1_uptake_target_variable=''
   CHARACTER(len=64)  :: p2_uptake_target_variable=''
   CHARACTER(len=64)  :: n_excretion_target_variable=''
   CHARACTER(len=64)  :: n_mortality_target_variable=''
   CHARACTER(len=64)  :: n1_uptake_target_variable=''
   CHARACTER(len=64)  :: n2_uptake_target_variable=''
   CHARACTER(len=64)  :: n3_uptake_target_variable=''
   CHARACTER(len=64)  :: n4_uptake_target_variable=''
   CHARACTER(len=64)  :: c_excretion_target_variable=''
   CHARACTER(len=64)  :: c_mortality_target_variable=''
   CHARACTER(len=64)  :: c_uptake_target_variable=''
   CHARACTER(len=64)  :: do_uptake_target_variable=''
   CHARACTER(len=64)  :: si_excretion_target_variable=''
   CHARACTER(len=64)  :: si_mortality_target_variable=''
   CHARACTER(len=64)  :: si_uptake_target_variable=''
   CHARACTER(len=128) :: dbase='aed2_malgae_pars.nml'
   AED_REAL           :: zerolimitfudgefactor = 0.9 * 3600
   AED_REAL           :: min_rho, max_rho
   LOGICAL            :: extra_debug = .false.
   AED_REAL           :: slough_stress = one_
   INTEGER            :: simMalgHSI = 0
   INTEGER            :: i, n_zones
   INTEGER            :: active_zones(1000)

   NAMELIST /aed2_macroalgae/ num_malgae, the_malgae, settling,resuspension,&
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    si_excretion_target_variable,si_mortality_target_variable, &
                      si_uptake_target_variable,                               &
                    dbase, zerolimitfudgefactor, extra_debug, extra_diag,      &
                    slough_stress, simMalgHSI, n_zones, active_zones
!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed2_macroalgae initialization"

   ! Default settings
   settling = _MOB_CONST_
   resuspension = zero_

   ! Read the namelist, and set module parameters
   read(namlst,nml=aed2_macroalgae,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_macroalgae'

   dtlim = zerolimitfudgefactor
   IF( extra_debug ) extra_diag = .true.       ! legacy use of extra_debug
   data%min_rho = min_rho ; data%max_rho = max_rho
   data%slough_stress = slough_stress
   data%simMalgHSI = simMalgHSI
   data%n_zones = n_zones
   IF( n_zones>0 ) THEN
     ALLOCATE(data%active_zones(n_zones))
     DO i=1,n_zones
       data%active_zones(i) = active_zones(i)
     ENDDO
   ENDIF

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   !     and are converted in here to values per second.
   CALL aed2_macroalgae_load_params(data,dbase,num_malgae,the_malgae,settling,resuspension)

   CALL aed2_bio_temp_function(data%num_malgae,              &
                               data%phytos%theta_growth,     &
                               data%phytos%T_std,            &
                               data%phytos%T_opt,            &
                               data%phytos%T_max,            &
                               data%phytos%aTn,              &
                               data%phytos%bTn,              &
                               data%phytos%kTn,              &
                               data%phytos%p_name)

   ! Register link to nutrient pools, if variable names are provided in namelist.
   data%do_Pexc = p_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Pexctarget = aed2_locate_variable( p_excretion_target_variable)
   ENDIF
   data%do_Nexc = n_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Nexctarget = aed2_locate_variable( n_excretion_target_variable)
   ENDIF
   data%do_Cexc = c_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Cexctarget = aed2_locate_variable( c_excretion_target_variable)
   ENDIF
   data%do_Siexc = si_excretion_target_variable .NE. ''
   IF (data%do_Siexc) THEN
     data%id_Siexctarget = aed2_locate_variable( si_excretion_target_variable)
   ENDIF

   data%do_Pmort = p_mortality_target_variable .NE. ''
   IF (data%do_Pmort) THEN
     data%id_Pmorttarget = aed2_locate_variable( p_mortality_target_variable)
   ENDIF
   data%do_Nmort = n_mortality_target_variable .NE. ''
   IF (data%do_Nmort) THEN
     data%id_Nmorttarget = aed2_locate_variable( n_mortality_target_variable)
   ENDIF
   data%do_Cmort = c_mortality_target_variable .NE. ''
   IF (data%do_Cmort) THEN
     data%id_Cmorttarget = aed2_locate_variable( c_mortality_target_variable)
   ENDIF
   data%do_Simort = si_mortality_target_variable .NE. ''
   IF (data%do_Simort) THEN
     data%id_Simorttarget = aed2_locate_variable( si_mortality_target_variable)
   ENDIF

   data%npup = 0
   IF (p1_uptake_target_variable .NE. '') data%npup = 1
   IF (p2_uptake_target_variable .NE. '') data%npup = 2
   data%do_Puptake = .FALSE.
   IF (data%npup>0) data%do_Puptake=.TRUE.
   IF (data%do_Puptake) THEN
     IF (data%npup>0) data%id_Pupttarget(1) = aed2_locate_variable( p1_uptake_target_variable); ifrp=1
     IF (data%npup>1) data%id_Pupttarget(2) = aed2_locate_variable( p2_uptake_target_variable); idop=2
   ENDIF
   data%nnup = 0
   IF (n1_uptake_target_variable .NE. '') data%nnup = 1
   IF (n2_uptake_target_variable .NE. '') data%nnup = 2
   IF (n3_uptake_target_variable .NE. '') data%nnup = 3
   IF (n4_uptake_target_variable .NE. '') data%nnup = 4
   data%do_Nuptake = .false.
   IF (data%nnup>0) data%do_Nuptake=.true.
   IF (data%do_Nuptake) THEN
     IF (data%nnup>0) data%id_Nupttarget(1) = aed2_locate_variable( n1_uptake_target_variable); ino3=1
     IF (data%nnup>1) data%id_Nupttarget(2) = aed2_locate_variable( n2_uptake_target_variable); inh4=2
     IF (data%nnup>2) data%id_Nupttarget(3) = aed2_locate_variable( n3_uptake_target_variable); idon=3
     IF (data%nnup>3) data%id_Nupttarget(4) = aed2_locate_variable( n4_uptake_target_variable); in2 =4
   ENDIF
   data%do_Cuptake = c_uptake_target_variable .NE. ''
   IF (data%do_Cuptake) THEN
     data%id_Cupttarget = aed2_locate_variable( c_uptake_target_variable)
   ENDIF
   data%do_DOuptake = do_uptake_target_variable .NE. ''
   IF (data%do_DOuptake) THEN
     data%id_DOupttarget = aed2_locate_variable( do_uptake_target_variable)
   ENDIF
   data%do_Siuptake = si_uptake_target_variable .NE. ''
   IF (data%do_Siuptake) THEN
     data%id_Siupttarget = aed2_locate_variable( si_uptake_target_variable)
   ENDIF

   ! Register diagnostic variables
   data%id_GPP = aed2_define_diag_variable('GPP','mmol/m**3/d',  'gross primary production')
  !data%id_NCP = aed2_define_diag_variable('NCP','mmol/m**3/d',  'net community production')
  !data%id_PPR = aed2_define_diag_variable('PPR','-','macroalgae p/r ratio (gross)')
  !data%id_NPR = aed2_define_diag_variable('NPR','-','macroalgae p/r ratio (net)')

  !data%id_NUP = aed2_define_diag_variable('NUP','mmol/m**3/d','nitrogen uptake')
  !data%id_PUP = aed2_define_diag_variable('PUP','mmol/m**3/d','phosphorous uptake')
  !data%id_CUP = aed2_define_diag_variable('CUP','mmol/m**3/d','carbon uptake')

   data%id_dPAR = aed2_define_diag_variable('PAR','W/m**2', 'photosynthetically active radiation')
   data%id_TCHLA = aed2_define_diag_variable('EXTC','/m', 'extinction')
   data%id_TMALG = aed2_define_diag_variable('TMALG','g DW/m**2', 'total macroalgae')
   data%id_TIN = aed2_define_diag_variable('IN','mmol/m**3', 'total phy nitrogen')
   data%id_TIP = aed2_define_diag_variable('IP','mmol/m**3', 'total phy phosphorus')

   IF ( simMalgHSI>0 ) &
     data%id_mhsi = aed2_define_sheet_diag_variable('HSI','-', 'macroalgae habitat suitability')

   ! Register environmental dependencies
   data%id_tem = aed2_locate_global('temperature')
   data%id_sal = aed2_locate_global('salinity')
   data%id_par = aed2_locate_global('par')
   data%id_I_0 = aed2_locate_global_sheet('par_sf')
   data%id_dz = aed2_locate_global('layer_ht')
   data%id_extc = aed2_locate_global('extc_coef')
   data%id_dens = aed2_locate_global('density')
   data%id_taub = aed2_locate_global_sheet('taub')
   data%id_sedzone = aed2_locate_global_sheet('sed_zone')

END SUBROUTINE aed2_define_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_macroalgae(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of macroalgae biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: phy, tphy, tin, tip, tchla
   AED_REAL :: INi, IPi
   AED_REAL :: pup
   AED_REAL :: no3up,nh4up
   AED_REAL :: cup, rsiup
   AED_REAL :: temp, par, Io, salinity, extc, dz, salt
   AED_REAL :: primprod(data%num_malgae), exudation(data%num_malgae), &
               a_nfix(data%num_malgae), respiration(data%num_malgae)
   AED_REAL :: cuptake(data%num_malgae), cexcretion(data%num_malgae), cmortality(data%num_malgae)
   AED_REAL :: nuptake(data%num_malgae,1:4), nexcretion(data%num_malgae), nmortality(data%num_malgae)
   AED_REAL :: puptake(data%num_malgae,1:2), pexcretion(data%num_malgae), pmortality(data%num_malgae)
   AED_REAL :: siuptake(data%num_malgae), siexcretion(data%num_malgae), simortality(data%num_malgae)
   AED_REAL :: fT, fNit, fPho, fSil, fI, fXl, fSal, PNf
   AED_REAL :: upTot,net_cuptake

   INTEGER  :: phy_i,c
   AED_REAL :: flux, available

!-------------------------------------------------------------------------------
!BEGIN


!RETURN

   ! Retrieve current environmental conditions.
   temp = _STATE_VAR_(data%id_tem)    ! local temperature
   salinity = _STATE_VAR_(data%id_sal)! local salinity
   par = _STATE_VAR_(data%id_par)     ! local photosynthetically active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation

   pup = 0.
   ! Retrieve current (local) state variable values.
   IF (data%do_Puptake)  pup = _STATE_VAR_(data%id_Pupttarget(1))

   no3up = 0.
   nh4up = 0.
   IF (data%do_Nuptake) THEN
       no3up = _STATE_VAR_(data%id_Nupttarget(1))
       nh4up = _STATE_VAR_(data%id_Nupttarget(2))
   ENDIF
   cup = 0.
   IF (data%do_Cuptake)  cup = _STATE_VAR_(data%id_Cupttarget)
   rsiup = 0.
   IF (data%do_Siuptake)  rsiup = _STATE_VAR_(data%id_Siupttarget)

   tphy = 0.0
   tchla = 0.0
   tin  = 0.0
   tip  = 0.0

   INi = 0.
   IPi = 0.

   DO phy_i=1,data%num_malgae

      primprod(phy_i)    = zero_
      exudation(phy_i)   = zero_
      a_nfix(phy_i)      = zero_
      respiration(phy_i) = zero_

      cuptake(phy_i)     = zero_
      cexcretion(phy_i)  = zero_
      cmortality(phy_i)  = zero_
      nuptake(phy_i,:)   = zero_
      nexcretion(phy_i)  = zero_
      nmortality(phy_i)  = zero_
      puptake(phy_i,:)   = zero_
      pexcretion(phy_i)  = zero_
      pmortality(phy_i)  = zero_

      ! Retrieve this macroalgae group
      phy = _STATE_VAR_(data%id_p(phy_i))

      ! Get the temperature limitation function
      fT = fTemp_function(data%phytos(phy_i)%fT_Method,    &
                          data%phytos(phy_i)%T_max,        &
                          data%phytos(phy_i)%T_std,        &
                          data%phytos(phy_i)%theta_growth, &
                          data%phytos(phy_i)%aTn,          &
                          data%phytos(phy_i)%bTn,          &
                          data%phytos(phy_i)%kTn,temp)

      !fSal = fTemp_function(salinity, minS, Smin, Smax, maxS )
      fSal = one_ ! fSal_function(salinity, 25., 30., 45., 80. )
      salt = salinity
      IF( salt<=5. ) THEN
         fSal = zero_
      ELSE IF ( salt>5. .AND. salt<=18.  ) THEN
         fSal = 0. + ( (salt-5.)/(18.-5.) )
      ELSE IF ( salt>18. .AND. salt<=40. ) THEN
         fSal = one_
      ELSE IF ( salt>40. .AND. salt<=65. ) THEN
         fSal = 1. - ( (salt-40.)/(65.-40.) )
      ELSE IF ( salt>65. ) THEN
         fSal = zero_
      ENDIF

      ! Get the light and nutrient limitation.
      ! NITROGEN.
      fNit = 0.0
      IF(data%phytos(phy_i)%simINDynamics /= 0) THEN
         ! IN variable available
         INi = _STATE_VAR_(data%id_in(phy_i))
      ELSE
         ! Assumed constant IN:
         INi = phy*data%phytos(phy_i)%X_ncon
      END IF

      ! Estimate fN limitation from IN or ext N value
      IF(data%phytos(phy_i)%simINDynamics > 1) THEN
         IF (phy > data%phytos(phy_i)%p0) THEN
            fNit = INi / phy
            fNit = phyto_fN(data%phytos,phy_i,IN=fNit)
         ENDIF
         IF (phy > zero_ .AND. phy <= data%phytos(phy_i)%p0) THEN
            fNit = phyto_fN(data%phytos,phy_i,din=no3up+nh4up)
         ENDIF
      ELSE
         fNit = phyto_fN(data%phytos,phy_i,din=no3up+nh4up)
      ENDIF
      IF (data%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixer: apply no N limitation. N Fixation ability
         ! depends on DIN concentration
         a_nfix = (one_ - fNit)
         fNit = one_
      ENDIF


      ! PHOSPHOROUS.
      fPho = zero_
      IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! IP variable available
         IPi = _STATE_VAR_(data%id_ip(phy_i))
      ELSE
         ! Assumed constant IP:
         IPi = phy*data%phytos(phy_i)%X_pcon
      END IF

      ! Estimate fP limitation from IP or ext P value
      IF (data%phytos(phy_i)%simIPDynamics > 1) THEN
         IF (phy > data%phytos(phy_i)%p0) THEN
            fPho = IPi / phy
            fPho = phyto_fP(data%phytos,phy_i,IP=fPho)
         ENDIF
         IF (phy > zero_ .AND. phy <= data%phytos(phy_i)%p0) THEN
            fPho = phyto_fP(data%phytos,phy_i,frp=pup)
         ENDIF
      ELSE
         fPho = phyto_fP(data%phytos,phy_i,frp=pup)
      ENDIF

      ! SILICA.
      fSil = one_ !phyto_fSi(data%phytos,phy_i,rsiup)


      ! LIGHT
      extc = _STATE_VAR_(data%id_extc)
      ! dz = 0.5     !MH: to fix
      dz = _STATE_VAR_(data%id_dz)
      fI = photosynthesis_irradiance(data%phytos(phy_i)%lightModel, &
               data%phytos(phy_i)%I_K, data%phytos(phy_i)%I_S, par, extc, Io, dz)
      ! fI = 0.1


      ! METAL AND TOXIC EFFECTS
      fXl = 1.0

      ! Primary production rate
      primprod(phy_i) = data%phytos(phy_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl * fSal

      ! Adjust primary production rate for nitrogen fixers
      IF (data%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixing species, and the growth rate to  must be reduced
         ! to compensate for the increased metabolic cost of this process
         primprod(phy_i) = primprod(phy_i) * (data%phytos(phy_i)%k_nfix + &
                           (1.0-a_nfix(phy_i))*(1.0-data%phytos(phy_i)%k_nfix))
      ENDIF


      ! Respiration and general metabolic loss
      respiration(phy_i) = bio_respiration(data%phytos(phy_i)%R_resp,data%phytos(phy_i)%theta_resp,temp)

      ! Salinity stress effect on respiration
      !fSal =  phyto_salinity(data%phytos,phy_i,salinity)
      !respiration(phy_i) = respiration(phy_i) * fSal

      ! photo-exudation
      exudation(phy_i) = primprod(phy_i)*data%phytos(phy_i)%f_pr

      ! Limit respiration if at the min biomass to prevent
      ! leak in the C mass balance
      IF (phy <= data%phytos(phy_i)%p0) THEN
         respiration(phy_i) = zero_
         exudation(phy_i) = zero_
      ENDIF

      ! write(*,"(4X,'limitations (fT,fI,fN,fP,fSi,Io, par, mu): ',9F9.2)")fT,fI,fNit,fPho,fSil,Io,par,primprod*secs_per_day


      ! Carbon uptake and excretion
      cuptake(phy_i)    = -primprod(phy_i) * phy
      cexcretion(phy_i) = (data%phytos(phy_i)%k_fdom*(1.0-data%phytos(phy_i)%k_fres)*respiration(phy_i)+exudation(phy_i)) * phy
      cmortality(phy_i) = ((1.0-data%phytos(phy_i)%k_fdom)*(1.0-data%phytos(phy_i)%k_fres)*respiration(phy_i)) * phy

      ! Nitrogen uptake and excretion
      CALL phyto_internal_nitrogen(data%phytos,phy_i,data%do_N2uptake,phy,INi,primprod(phy_i),&
                             fT,no3up,nh4up,a_nfix(phy_i),respiration(phy_i),exudation(phy_i),PNf,&
                                   nuptake(phy_i,:),nexcretion(phy_i),nmortality(phy_i))

      ! Phosphorus uptake and excretion
      CALL phyto_internal_phosphorus(data%phytos,phy_i,data%npup,phy,IPi,primprod(phy_i),&
                                 fT,pup,respiration(phy_i),exudation(phy_i),&
                                         puptake(phy_i,:),pexcretion(phy_i),pmortality(phy_i))

      ! Silica uptake and excretion
      IF (data%phytos(phy_i)%simSiUptake > 0) THEN
         siuptake(phy_i)    =-data%phytos(phy_i)%X_sicon * primprod(phy_i) * phy
         siexcretion(phy_i) = data%phytos(phy_i)%X_sicon * (data%phytos(phy_i)%k_fdom*respiration(phy_i)+exudation(phy_i)) * phy
         simortality(phy_i) = data%phytos(phy_i)%X_sicon * ((1.0-data%phytos(phy_i)%k_fdom)*respiration(phy_i)) * phy
      ELSE
         siuptake(phy_i)    = zero_
         siexcretion(phy_i) = zero_
         simortality(phy_i) = zero_
      ENDIF

      ! Diagnostic info
      _DIAG_VAR_(data%id_NtoP(phy_i)) =  INi/IPi

      IF (extra_diag) THEN
         _DIAG_VAR_(data%id_fT(phy_i)) =  fT
         _DIAG_VAR_(data%id_fI(phy_i)) =  fI
         _DIAG_VAR_(data%id_fNit(phy_i)) =  fNit
         _DIAG_VAR_(data%id_fPho(phy_i)) =  fPho
         _DIAG_VAR_(data%id_fSil(phy_i)) =  fSil
         _DIAG_VAR_(data%id_fSal(phy_i)) =  fSal
      ENDIF
   END DO


   !-----------------------------------------------------------------
   ! Check uptake values for availability to prevent -ve numbers

   ! pup   - p available
   ! no3up - no3 available
   ! nh4up - nh4 available
   ! cup   - c available
   ! rsiup - Si available

   IF (data%do_Puptake) THEN
      upTot = sum(puptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= pup ) THEN
         DO phy_i=1,data%num_malgae
            puptake(phy_i,1) = (pup*0.99/dtlim) * (puptake(phy_i,1)/upTot)
         ENDDO
      ENDIF
   ENDIF

   IF (data%do_Nuptake) THEN
      upTot = sum(nuptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= no3up ) THEN
         DO phy_i=1,data%num_malgae
            nuptake(phy_i,1) = (no3up*0.99/dtlim) * (nuptake(phy_i,1)/upTot)
         ENDDO
      ENDIF

      upTot = sum(nuptake(:,2))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= nh4up ) THEN
         DO phy_i=1,data%num_malgae
            nuptake(phy_i,2) = (nh4up*0.99/dtlim) * (nuptake(phy_i,2)/upTot)
         ENDDO
      ENDIF
   ENDIF
   IF (data%do_Cuptake) THEN
      upTot = sum(cuptake)*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= cup ) THEN
         DO phy_i=1,data%num_malgae
            cuptake(phy_i) = (cup*0.99/dtlim) * (cuptake(phy_i)/upTot)
         ENDDO
      ENDIF
   ENDIF
!  IF (data%do_DOuptake) THEN
!     !
!  ENDIF
   IF (data%do_Siuptake) THEN
      upTot = sum(siuptake)*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= rsiup ) THEN
         DO phy_i=1,data%num_malgae
            siuptake(phy_i) = (rsiup*0.99/dtlim) * (siuptake(phy_i)/upTot)
         ENDDO
      ENDIF
   ENDIF

   !-----------------------------------------------------------------
   ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER
   net_cuptake = zero_
   DO phy_i=1,data%num_malgae

      !# macroalgae PRODUCTION & RESPIRATION
      phy = _STATE_VAR_(data%id_p(phy_i))
      flux = (primprod(phy_i) - respiration(phy_i) - exudation(phy_i)) * phy
      available = MAX(zero_, phy - data%phytos(phy_i)%p0)
      IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
      _FLUX_VAR_(data%id_p(phy_i)) = _FLUX_VAR_(data%id_p(phy_i)) + ( flux)

      !# macroalgae INTERNAL NITROGEN
      IF (data%phytos(phy_i)%simINDynamics /= 0) THEN
         ! _FLUX_VAR_(data%id_in(phy_i)) = _FLUX_VAR_(data%id_in(phy_i)) + ( (-sum(nuptake) - nexcretion(phy_i) - nmortality(phy_i) )*INi )
         INi = _STATE_VAR_(data%id_in(phy_i))
         flux = (-sum(nuptake(phy_i,:)) - nexcretion(phy_i) - nmortality(phy_i) )
         available = MAX(zero_, INi - data%phytos(phy_i)%X_nmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_in(phy_i)) = _FLUX_VAR_(data%id_in(phy_i)) + ( flux)
      ENDIF

      !# macroalgae INTERNAL PHOSPHORUS
      IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! _FLUX_VAR_(data%id_ip(phy_i)) = _FLUX_VAR_(data%id_ip(phy_i)) + ( (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) ) )
         IPi = _STATE_VAR_(data%id_ip(phy_i))
         flux = (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) )
         available = MAX(zero_, IPi - data%phytos(phy_i)%X_pmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_ip(phy_i)) = _FLUX_VAR_(data%id_ip(phy_i)) + ( flux)
      ENDIF

      !# macroalgae CELL DENSITY
      IF ( data%id_rho(phy_i)>0 ) THEN
         ! density increases during carbohydrate creation (daytime)
         flux = zero_
         IF( par>zero_ ) THEN
           flux = data%phytos(phy_i)%c1 * &
              (one_ - EXP(-par/data%phytos(phy_i)%I_K) ) - data%phytos(phy_i)%c3
         ELSE
           ! darkness
           flux = -data%phytos(phy_i)%c3
         ENDIF
        _FLUX_VAR_(data%id_rho(phy_i)) = _FLUX_VAR_(data%id_rho(phy_i)) + flux
         ! check maximum/minimum density are not exceeded
         IF( _STATE_VAR_(data%id_rho(phy_i))>data%max_rho ) THEN
            _FLUX_VAR_(data%id_rho(phy_i)) =zero_
            _STATE_VAR_(data%id_rho(phy_i))=data%max_rho
         ENDIF
         IF( _STATE_VAR_(data%id_rho(phy_i))<data%min_rho ) THEN
             _FLUX_VAR_(data%id_rho(phy_i)) =zero_
             _STATE_VAR_(data%id_rho(phy_i))=data%min_rho
         ENDIF
      ENDIF

      ! BIOGEOCHEMICAL FEEDBACKS
      ! Now manage uptake of nutrients, CO2 and DO - these cumulative fluxes already limited above loop
      IF (data%do_Puptake) THEN
         DO c = 1,data%npup
            _FLUX_VAR_(data%id_Pupttarget(c)) = _FLUX_VAR_(data%id_Pupttarget(c)) + ( puptake(phy_i,c))
         ENDDO
      ENDIF
      IF (data%do_Nuptake) THEN
         DO c = 1,data%nnup
            _FLUX_VAR_(data%id_Nupttarget(c)) = _FLUX_VAR_(data%id_Nupttarget(c)) + ( nuptake(phy_i,c))
         ENDDO
      ENDIF
      IF (data%do_Cuptake) THEN
         _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) + (  cuptake(phy_i) - respiration(phy_i)*data%phytos(phy_i)%k_fres*phy )
         net_cuptake = net_cuptake + (cuptake(phy_i) - respiration(phy_i)*data%phytos(phy_i)%k_fres*phy)
      ENDIF
      IF (data%do_DOuptake) THEN
         _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) + ( -cuptake(phy_i) + respiration(phy_i)*data%phytos(phy_i)%k_fres*phy )
      ENDIF
      IF (data%do_Siuptake) THEN
         _FLUX_VAR_(data%id_Siupttarget) = _FLUX_VAR_(data%id_Siupttarget) + ( siuptake(phy_i))
      ENDIF
      ! Now manage mortality contributions to POM
      IF (data%do_Pmort) THEN
         _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + (pmortality(phy_i))
      ENDIF
      IF (data%do_Nmort) THEN
         _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + (nmortality(phy_i))
      ENDIF
      IF (data%do_Cmort) THEN
         _FLUX_VAR_(data%id_Cmorttarget) = _FLUX_VAR_(data%id_Cmorttarget) + (cmortality(phy_i))
      ENDIF
      IF (data%do_Simort) THEN
         _FLUX_VAR_(data%id_Simorttarget) = _FLUX_VAR_(data%id_Simorttarget) + (simortality(phy_i))
      ENDIF
      ! Now manage excretion/exudation contributions to DOM
      IF (data%do_Pexc) THEN
         _FLUX_VAR_(data%id_Pexctarget) = _FLUX_VAR_(data%id_Pexctarget) + (pexcretion(phy_i))
      ENDIF
      IF (data%do_Nexc) THEN
         _FLUX_VAR_(data%id_Nexctarget) = _FLUX_VAR_(data%id_Nexctarget) + (nexcretion(phy_i))
      ENDIF
      IF (data%do_Cexc) THEN
         _FLUX_VAR_(data%id_Cexctarget) = _FLUX_VAR_(data%id_Cexctarget) + (cexcretion(phy_i))
      ENDIF
      IF (data%do_Siexc) THEN
         _FLUX_VAR_(data%id_Siexctarget) = _FLUX_VAR_(data%id_Siexctarget) + (siexcretion(phy_i))
      ENDIF

      !-----------------------------------------------------------------
      ! UPDATE DIAGNOSTIC VARIABLES

      ! total macroalgae carbon
      tphy = tphy + phy
      ! total internal nutrients
      tin = tin + INi
      tip = tip + IPi
   ENDDO

   ! Set diagnostic arrays for combined assemblage properties
   _DIAG_VAR_(data%id_GPP) =  sum(cuptake)*secs_per_day
  !_DIAG_VAR_(data%id_NCP) =  net_cuptake*secs_per_day
  !_DIAG_VAR_(data%id_PPR) =  -999. !sum(cuptake) / ( sum(cuptake) - net_cuptake)
  !_DIAG_VAR_(data%id_NPR) =  -999. !net_cuptake / ( sum(cuptake) - net_cuptake)
  !_DIAG_VAR_(data%id_NUP) =  sum(nuptake)*secs_per_day
  !_DIAG_VAR_(data%id_PUP) =  sum(puptake)*secs_per_day
  !_DIAG_VAR_(data%id_CUP) =  sum(cuptake)*secs_per_day

  !_DIAG_VAR_(data%id_dPAR) =  par
  !_DIAG_VAR_(data%id_TCHLA)=  tchla
!  _DIAG_VAR_(data%id_TMALG) =  tphy
   _DIAG_VAR_(data%id_TIN)   =  tin
   _DIAG_VAR_(data%id_TIP)   =  tip

END SUBROUTINE aed2_calculate_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_benthic_macroalgae(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of macroalgae.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: malg_i,phy_i,c
   AED_REAL :: malg,temp,extc,par,dz,Io,fI,bottom_stress,depth,light,matz        ! State
   AED_REAL :: malg_flux,malg_prod,malg_resp                    ! Fluxes

   AED_REAL :: tphy, tin, tip, tchla
   AED_REAL :: INi, IPi
   AED_REAL :: pup
   AED_REAL :: no3up,nh4up
   AED_REAL :: cup, rsiup
   AED_REAL :: salinity, salt
   AED_REAL :: primprod(data%num_malgae), exudation(data%num_malgae), &
               a_nfix(data%num_malgae), respiration(data%num_malgae)
   AED_REAL :: cuptake(data%num_malgae), cexcretion(data%num_malgae), cmortality(data%num_malgae)
   AED_REAL :: nuptake(data%num_malgae,1:4), nexcretion(data%num_malgae), nmortality(data%num_malgae)
   AED_REAL :: puptake(data%num_malgae,1:2), pexcretion(data%num_malgae), pmortality(data%num_malgae)
   AED_REAL :: siuptake(data%num_malgae), siexcretion(data%num_malgae), simortality(data%num_malgae)
   AED_REAL :: fT, fNit, fPho, fSil, fXl, fSal, PNf
   AED_REAL :: upTot,net_cuptake,available,flux
!
!-------------------------------------------------------------------------------
!BEGIN

  ! Benthic light fraction and extinction, for diagnostics
  extc = _STATE_VAR_(data%id_extc) ! cell extinction
  depth = _STATE_VAR_(data%id_dz)  ! cell depth
  IF( depth<0.1 ) THEN
    light = 100.
  ELSE
    light = 100. * exp(-extc*(depth-0.08))
  ENDIF
  _DIAG_VAR_(data%id_dPAR) =  light
  _DIAG_VAR_(data%id_TCHLA) = extc

  matz = _STATE_VAR_S_(data%id_sedzone)

  ! Loop through selected groups, and determine if benthic/attahced is active
  DO malg_i=1,data%num_malgae

     fSal = zero_; fI = zero_; fNit = zero_; fPho = zero_
     ! Process benthic / attached macroalgae
     IF ( data%phytos(malg_i)%settling == _MOB_ATTACHED_ .AND. &
                    in_zone_set(matz,data%active_zones) ) THEN

       ! Get local conditions
       salinity = _STATE_VAR_(data%id_sal)  ! local salinity 
       temp = _STATE_VAR_(data%id_tem)  ! local temperature
       dz   = _STATE_VAR_(data%id_dz)     ! cell depth
       par  = _STATE_VAR_(data%id_par)   ! local photosynthetically active radiation
       Io   = _STATE_VAR_S_(data%id_I_0)  ! surface short wave radiation
       bottom_stress = _STATE_VAR_S_(data%id_taub)

       phy_i = malg_i

       !---- FROM ABOVE ----
       pup = 0.
       ! Retrieve current (local) state variable values.
       IF (data%do_Puptake)  pup = _STATE_VAR_(data%id_Pupttarget(1))

       no3up = 0.
       nh4up = 0.
       IF (data%do_Nuptake) THEN
         no3up = _STATE_VAR_(data%id_Nupttarget(1))
         nh4up = _STATE_VAR_(data%id_Nupttarget(2))
       ENDIF
       cup = 0.
       IF (data%do_Cuptake)  cup = _STATE_VAR_(data%id_Cupttarget)

       tphy = 0.0
       tchla = 0.0
       tin  = 0.0
       tip  = 0.0

       INi = 0.
       IPi = 0.

       primprod(phy_i)    = zero_
       exudation(phy_i)   = zero_
       a_nfix(phy_i)      = zero_
       respiration(phy_i) = zero_

       cuptake(phy_i)     = zero_
       cexcretion(phy_i)  = zero_
       cmortality(phy_i)  = zero_
       nuptake(phy_i,:)   = zero_
       nexcretion(phy_i)  = zero_
       nmortality(phy_i)  = zero_
       puptake(phy_i,:)   = zero_
       pexcretion(phy_i)  = zero_
       pmortality(phy_i)  = zero_

       ! Retrieve this macroalgae group
       malg = _STATE_VAR_S_(data%id_pben(malg_i)) ! local malg density

       ! Get the temperature limitation function
       fT = fTemp_function(data%phytos(phy_i)%fT_Method,    &
                           data%phytos(phy_i)%T_max,        &
                           data%phytos(phy_i)%T_std,        &
                           data%phytos(phy_i)%theta_growth, &
                           data%phytos(phy_i)%aTn,          &
                           data%phytos(phy_i)%bTn,          &
                           data%phytos(phy_i)%kTn,temp)

       !fSal = fSal_function(salinity, minS, Smin, Smax, maxS )
       fSal = one_ !fSal_function(salinity, 5., 18., 35., 45. )

       salt = salinity
       IF( salt<=5. ) THEN
         fSal = zero_
       ELSE IF ( salt>5. .AND. salt<=18.  ) THEN
         fSal = 0. + ( (salt-5.)/(18.-5.) )
       ELSE IF ( salt>18. .AND. salt<=40. ) THEN
         fSal = one_
       ELSE IF ( salt>40. .AND. salt<=85. ) THEN
         fSal = 1. - ( (salt-40.)/(85.-40.) )
       ELSE IF ( salt>85. ) THEN
         fSal = zero_
       ENDIF
       
       ! Get the light and nutrient limitation.
       ! NITROGEN.
       fNit = 0.0
       IF(data%phytos(phy_i)%simINDynamics /= 0) THEN
          ! IN variable available
          INi = _STATE_VAR_S_(data%id_inben(phy_i))
       ELSE
          ! Assumed constant IN:
          INi = malg*data%phytos(phy_i)%X_ncon
       END IF

       ! Estimate fN limitation from IN or ext N value
       IF(data%phytos(phy_i)%simINDynamics > 1) THEN
          IF (malg > data%phytos(phy_i)%p0) THEN
             fNit = INi / malg
             fNit = phyto_fN(data%phytos,phy_i,IN=fNit)
          ENDIF
          IF (malg > zero_ .AND. malg <= data%phytos(phy_i)%p0) THEN
             fNit = phyto_fN(data%phytos,phy_i,din=no3up+nh4up)
          ENDIF
       ELSE
          fNit = phyto_fN(data%phytos,phy_i,din=no3up+nh4up)
       ENDIF
       IF (data%phytos(phy_i)%simNFixation /= 0) THEN
          ! Nitrogen fixer: apply no N limitation. N Fixation ability
          ! depends on DIN concentration
          a_nfix = (one_ - fNit)
          fNit = one_
       ENDIF

       ! PHOSPHOROUS.
       fPho = zero_
       IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
          ! IP variable available
          IPi = _STATE_VAR_S_(data%id_ipben(phy_i))
       ELSE
          ! Assumed constant IP:
          IPi = malg*data%phytos(phy_i)%X_pcon
       END IF

       ! Estimate fP limitation from IP or ext P value
       IF (data%phytos(phy_i)%simIPDynamics > 1) THEN
          IF (malg > data%phytos(phy_i)%p0) THEN
             fPho = IPi / malg
             fPho = phyto_fP(data%phytos,phy_i,IP=fPho)
          ENDIF
          IF (malg > zero_ .AND. malg <= data%phytos(phy_i)%p0) THEN
             fPho = phyto_fP(data%phytos,phy_i,frp=pup)
          ENDIF
       ELSE
          fPho = phyto_fP(data%phytos,phy_i,frp=pup)
       ENDIF
       fPho = 1.0

       ! SILICA.
       fSil = 1.0

       ! METAL AND TOXIC EFFECTS
       fXl = 1.0
       IF(bottom_stress>0.09) fXl =0.0

       ! Compute P-I
       fI = photosynthesis_irradiance(0, &  ! 0 is vertical integral
            data%phytos(malg_i)%I_K, data%phytos(malg_i)%I_S, par, extc, Io, dz)

       IF (extra_diag) THEN
         _DIAG_VAR_S_(data%id_fT_ben(phy_i)) =  fT
         _DIAG_VAR_S_(data%id_fI_ben(phy_i)) =  fI
         _DIAG_VAR_S_(data%id_fNit_ben(phy_i)) =  fNit
         _DIAG_VAR_S_(data%id_fPho_ben(phy_i)) =  fPho
         _DIAG_VAR_S_(data%id_fSal_ben(phy_i)) =  fSal
       ENDIF

      ! Primary production rate
      primprod(phy_i) = data%phytos(phy_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl * fSal

      ! Respiration and general metabolic loss
      respiration(phy_i) = bio_respiration(data%phytos(phy_i)%R_resp,data%phytos(phy_i)%theta_resp,temp)

      ! Photo-exudation
      exudation(phy_i) = primprod(phy_i)*data%phytos(phy_i)%f_pr

      ! Limit respiration if at the min biomass to prevent leak in the C mass balance
      IF (malg <= data%phytos(phy_i)%p0) THEN
         respiration(phy_i) = zero_
         exudation(phy_i) = zero_
      ENDIF

      ! Carbon uptake and excretion
      cuptake(phy_i)    = -primprod(phy_i) * malg
      cexcretion(phy_i) = (data%phytos(phy_i)%k_fdom*(1.0-data%phytos(phy_i)%k_fres)*respiration(phy_i)+exudation(phy_i)) * malg
      cmortality(phy_i) = ((1.0-data%phytos(phy_i)%k_fdom)*(1.0-data%phytos(phy_i)%k_fres)*respiration(phy_i)) * malg

      ! Nitrogen uptake and excretion
      CALL phyto_internal_nitrogen(data%phytos,phy_i,data%do_N2uptake,malg,INi,primprod(phy_i),&
                             fT,no3up,nh4up,a_nfix(phy_i),respiration(phy_i),exudation(phy_i),PNf,&
                                   nuptake(phy_i,:),nexcretion(phy_i),nmortality(phy_i))

      ! Phosphorus uptake and excretion
      CALL phyto_internal_phosphorus(data%phytos,phy_i,data%npup,malg,IPi,primprod(phy_i),&
                                 fT,pup,respiration(phy_i),exudation(phy_i),&
                                         puptake(phy_i,:),pexcretion(phy_i),pmortality(phy_i))

      malg_flux = (primprod(phy_i)-respiration(phy_i)-exudation(phy_i))*malg


       ! Update the attached biomass, and water O2/CO2/Nuts
       _FLUX_VAR_B_(data%id_pben(malg_i)) = _FLUX_VAR_B_(data%id_pben(malg_i)) + malg_flux
       !# macroalgae INTERNAL NITROGEN
       IF (data%phytos(phy_i)%simINDynamics /= 0) THEN
          INi = _STATE_VAR_S_(data%id_inben(phy_i))
          flux = (-sum(nuptake(phy_i,:)) - nexcretion(phy_i) - nmortality(phy_i) )
          available = MAX(zero_, INi - data%phytos(phy_i)%X_nmin*malg)
          IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
          _FLUX_VAR_B_(data%id_inben(phy_i)) = _FLUX_VAR_B_(data%id_inben(phy_i)) + ( flux)
       ENDIF
       !# macroalgae INTERNAL PHOSPHORUS
       IF (data%phytos(phy_i)%simIPDynamics /= 0) THEN
          IPi = _STATE_VAR_S_(data%id_ipben(phy_i))
          flux = (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) )
          available = MAX(zero_, IPi - data%phytos(phy_i)%X_pmin*malg)
          IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
          _FLUX_VAR_B_(data%id_ipben(phy_i)) = _FLUX_VAR_B_(data%id_ipben(phy_i)) + ( flux)
       ENDIF
       IF (data%do_DOuptake) THEN
         _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) + malg_flux
       ENDIF
       IF (data%do_Cuptake) THEN
         _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) - malg_flux
       ENDIF
       IF (data%do_Puptake) THEN
          DO c = 1,data%npup
             _FLUX_VAR_(data%id_Pupttarget(c)) = _FLUX_VAR_(data%id_Pupttarget(c)) + ( puptake(phy_i,c))
          ENDDO
       ENDIF
       IF (data%do_Nuptake) THEN
          DO c = 1,data%nnup
             _FLUX_VAR_(data%id_Nupttarget(c)) = _FLUX_VAR_(data%id_Nupttarget(c)) + ( nuptake(phy_i,c))
          ENDDO
       ENDIF
       ! OM exretion here


       ! Redistribute into the water column if sloughing occurs.
      IF( bottom_stress>data%slough_stress*2. ) THEN
         malg_flux = 0.3*malg
         _FLUX_VAR_(data%id_p(malg_i)) = _FLUX_VAR_(data%id_p(malg_i)) + malg_flux
         _FLUX_VAR_B_(data%id_pben(malg_i)) = _FLUX_VAR_B_(data%id_pben(malg_i)) - malg_flux

       ELSEIF( bottom_stress>data%slough_stress .AND. salinity>60.) THEN
         malg_flux = 0.67*malg
         _FLUX_VAR_(data%id_p(malg_i)) = _FLUX_VAR_(data%id_p(malg_i)) + malg_flux
         _FLUX_VAR_B_(data%id_pben(malg_i)) = _FLUX_VAR_B_(data%id_pben(malg_i)) - malg_flux

         !% in & ip here
       ENDIF

       ! Update the diagnostic variables
       _DIAG_VAR_(data%id_TMALG) =  (malg+_STATE_VAR_(data%id_p(malg_i))*depth )*(12.*1e-3/0.5)
     ENDIF

     IF( data%simMalgHSI == malg_i) THEN
       _DIAG_VAR_S_(data%id_mhsi) = min( fSal, fI, fNit, fPho ) * fT
     ENDIF

   ENDDO
   !_DIAG_VAR_(data%id_TMALG) =  _DIAG_VAR_(data%id_TMALG)*dz *(12.*1e-3/0.5)

END SUBROUTINE aed2_calculate_benthic_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_mobility_macroalgae(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values for macroalgae cells
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: temp, par, rho_p, Io
   AED_REAL :: vvel
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: IN, IC, Q, Qmax
   INTEGER  :: phy_i
!
!-------------------------------------------------------------------------------
!BEGIN

   DO phy_i=1,data%num_malgae
      SELECT CASE (data%phytos(phy_i)%settling)

         CASE ( _MOB_OFF_ )
            ! disable settling by setting vertical velocity to 0
            vvel = zero_

         CASE ( _MOB_CONST_ )
            ! constant settling velocity using user provided value
            vvel = data%phytos(phy_i)%w_p

         CASE ( _MOB_TEMP_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_dens)
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%phytos(phy_i)%w_p*mu20*pw / ( mu*pw20 )

         CASE ( _MOB_STOKES_ )
            ! settling velocity based on Stokes Law calculation and cell density
            pw = _STATE_VAR_(data%id_dens)       ! water density
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)                 ! water dynamic viscosity
            IF( data%id_rho(phy_i)>0 ) THEN
              rho_p = _STATE_VAR_(data%id_rho(phy_i))  ! cell density
            ELSE
              rho_p = data%phytos(phy_i)%rho_phy
            ENDIF
            vvel = -9.807*(data%phytos(phy_i)%d_phy**2.)*( rho_p-pw ) / ( 18.*mu )

          CASE ( _MOB_ATTACHED_ )
            ! ! settling velocity based on Stokes Law calculation and cell density
            ! pw = _STATE_VAR_(data%id_dens)       ! water density
            ! temp = _STATE_VAR_(data%id_tem)
            ! mu = water_viscosity(temp)                 ! water dynamic viscosity
            ! IF( data%id_rho(phy_i)>0 ) THEN
            !   rho_p = _STATE_VAR_(data%id_rho(phy_i))  ! cell density
            ! ELSE
            !   rho_p = data%phytos(phy_i)%rho_phy
            ! ENDIF
            ! vvel = -9.807*(data%phytos(phy_i)%d_phy**2.)*( rho_p-pw ) / ( 18.*mu )

            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_dens)
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%phytos(phy_i)%w_p*mu20*pw / ( mu*pw20 )

         CASE DEFAULT
            ! unknown settling/migration option selection
            vvel =  zero_

      END SELECT
      ! set global mobility array
      mobility(data%id_p(phy_i)) = vvel
      IF(extra_diag .AND. data%id_vvel(phy_i)>0) _DIAG_VAR_(data%id_vvel(phy_i)) = vvel
    ENDDO
END SUBROUTINE aed2_mobility_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_macroalgae(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: phy,malg,dz
   INTEGER  :: phy_i
!
!-----------------------------------------------------------------------
!BEGIN
   DO phy_i=1,data%num_malgae
      ! Retrieve current (local) state variable values.
      phy = _STATE_VAR_(data%id_p(phy_i))! macroalgae

      ! Self-shading with explicit contribution from background macroalgae concentration.
      extinction = extinction + (data%phytos(phy_i)%KePHY*phy)

      ! Process benthic / attached macroalgae
      IF ( data%phytos(phy_i)%settling == _MOB_ATTACHED_ ) THEN
        dz   = MAX(_STATE_VAR_(data%id_dz),0.1)     ! cell depth
        malg = _STATE_VAR_S_(data%id_pben(phy_i))/dz ! attached macroalgae
        extinction = extinction + (data%phytos(phy_i)%KePHY*malg)
      ENDIF
   ENDDO

END SUBROUTINE aed2_light_extinction_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_macroalgae
