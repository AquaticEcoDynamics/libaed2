!###############################################################################
!#                                                                             #
!# aed2_tracer.F90                                                             #
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
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed2.h"

!
MODULE aed2_tracer
!-------------------------------------------------------------------------------
! aed2_tracer --- tracer biogeochemical model
!
! The AED2 module tracer contains equations that describe a
! soluble or particle tracer, including decay, sediment interaction, and &
! resupension and settling
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE

   PUBLIC aed2_tracer_data_t

   TYPE,extends(aed2_model_data_t) :: aed2_tracer_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_ss(:), id_sfss(:)
      INTEGER :: id_retain
      INTEGER :: id_l_bot, id_tau_0, id_epsilon
      INTEGER :: id_temp, id_taub, id_salt, id_rho
      INTEGER :: id_d_taub, id_resus

      !# Module configuration
      INTEGER :: num_tracers
      INTEGER :: resuspension, settling

      !# Model parameters
      AED_REAL,ALLOCATABLE :: decay(:), Fsed(:), Ke_ss(:)
      AED_REAL,ALLOCATABLE :: w_ss(:), rho_ss(:), d_ss(:)
      AED_REAL,ALLOCATABLE :: fs(:), tau_0(:)
      AED_REAL             :: epsilon, tau_0_min, kTau_0, tau_r

     CONTAINS
         PROCEDURE :: define            => aed2_define_tracer
         PROCEDURE :: calculate         => aed2_calculate_tracer
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_tracer
         PROCEDURE :: mobility          => aed2_mobility_tracer
         PROCEDURE :: light_extinction  => aed2_light_extinction_tracer
        !PROCEDURE :: delete            => aed2_delete_tracer

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_tracer(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read in and the variables simulated
!  by the model are registered with AED2 core.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_tracer_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status,i
   LOGICAL  :: retention_time = .FALSE.
   INTEGER  :: resuspension = 0
   INTEGER  :: settling = 0
   INTEGER  :: num_tracers
   AED_REAL :: trace_initial = zero_
   AED_REAL :: decay(100)
   AED_REAL :: Fsed(100)
   AED_REAL :: Ke_ss(100)
   AED_REAL :: w_ss(100), rho_ss(100), d_ss(100)
   AED_REAL :: epsilon, tau_r, kTau_0
   AED_REAL :: tau_0(100)
   AED_REAL :: fs(100)
   CHARACTER(len=64) :: macrophyte_link_var = ''
   CHARACTER(4) :: trac_name

   NAMELIST /aed2_tracer/ num_tracers, decay, Fsed, Ke_ss, &
                          settling, w_ss, rho_ss, d_ss, &
                          resuspension, epsilon, tau_0, tau_r, Ktau_0, &
                          macrophyte_link_var, fs, &
                          trace_initial, retention_time
!
!-------------------------------------------------------------------------------
!BEGIN
   ! set default parameter values
   decay = zero_
   Fsed = zero_
   Ke_ss = 0.02
   w_ss = zero_
   d_ss = 1e-6
   rho_ss = 1.6e3
   epsilon = 0.02
   tau_r = 1.0
   tau_0 = 0.04
   kTau_0 = 1.0
   fs = 1.0

   ! Read the namelist
   read(namlst,nml=aed2_tracer,iostat=status)
   IF (status /= 0) STOP 'ERROR reading namelist aed2_tracer'

   ! Store parameter values in our own derived type
   data%num_tracers = num_tracers
   data%resuspension = resuspension
   data%settling = settling

   data%epsilon = epsilon
   data%tau_r = tau_r
   data%kTau_0 = kTau_0

   ! Setup tracers
   IF ( num_tracers > 0 ) THEN
      ALLOCATE(data%id_ss(num_tracers))
      ALLOCATE(data%decay(num_tracers)) ; data%decay(1:num_tracers) = decay(1:num_tracers)
      ALLOCATE(data%Fsed(num_tracers))  ; data%Fsed(1:num_tracers)  = Fsed(1:num_tracers)
      ALLOCATE(data%Ke_ss(num_tracers)) ; data%Ke_ss(1:num_tracers) = Ke_ss(1:num_tracers)

      ALLOCATE(data%w_ss(num_tracers))  ; data%w_ss(1:num_tracers)  = w_ss(1:num_tracers)
      ALLOCATE(data%d_ss(num_tracers))  ; data%d_ss(1:num_tracers)  = d_ss(1:num_tracers)
      ALLOCATE(data%rho_ss(num_tracers)); data%rho_ss(1:num_tracers)= rho_ss(1:num_tracers)

      ALLOCATE(data%tau_0(num_tracers)) ; data%tau_0(1:num_tracers) = tau_0(1:num_tracers)
      ALLOCATE(data%fs(num_tracers))    ; data%fs(1:num_tracers)    = fs(1:num_tracers)

      trac_name = 'ss0'
      ! Register state variables
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
                                             ! divide settling by secs_per_day to convert m/d to m/s
         data%id_ss(i) = aed2_define_variable(TRIM(trac_name),'mmol/m**3','tracer', &
                                                  trace_initial,minimum=zero_,mobility=(w_ss(i)/secs_per_day))
      ENDDO
   ENDIF

   ! Setup bottom arrays if spatially variable resuspension
   IF ( resuspension == 2 ) THEN
      data%id_tau_0 =  aed2_define_sheet_diag_variable('tau_0','N/m**2', 'dynamic bottom drag')
      data%id_epsilon =  aed2_define_sheet_diag_variable('epsilon','g/m**2/s', 'max resuspension rate')

      ALLOCATE(data%id_sfss(num_tracers))

      trac_name = 'fs0'
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_sfss(i) =  aed2_define_sheet_diag_variable(TRIM(trac_name),'-', 'sediment fraction of sed size')
      ENDDO

      IF ( macrophyte_link_var .NE. '' ) THEN
         data%id_l_bot = aed2_locate_sheet_variable(macrophyte_link_var)
         IF ( data%id_l_bot .LE. 0 ) THEN
            print *, "Macrophyte Link Variable ", TRIM(macrophyte_link_var), " is not defined."
            STOP
         ENDIF
      ENDIF
   ENDIF


   ! Retention time
   IF (retention_time) THEN
      data%id_retain = aed2_define_variable("ret",'sec','tracer',trace_initial,minimum=zero_)
   ELSE
      data%id_retain = -1
   ENDIF

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature')
   data%id_salt = aed2_locate_global('salinity')
   IF ( settling > 1 ) THEN
      data%id_rho = aed2_locate_global('density')
   ENDIF
   IF ( resuspension > 0 ) THEN
      data%id_taub = aed2_locate_global_sheet('taub')
      data%id_d_taub = aed2_define_sheet_diag_variable('d_taub','N/m**2',  'taub diagnostic')
      data%id_resus =  aed2_define_sheet_diag_variable('resus','g/m**2/s', 'resuspension rate')
    ENDIF

END SUBROUTINE aed2_define_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_tracer(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: trc
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i=1,data%num_tracers
      trc = _STATE_VAR_(data%id_ss(i))
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + data%decay(i)*trc
   ENDDO

   IF (data%id_retain < 1) RETURN
   _FLUX_VAR_(data%id_retain) = _FLUX_VAR_(data%id_retain) + 1.0
END SUBROUTINE aed2_calculate_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_tracer(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED tracer.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss, bottom_stress

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0, resus_flux = 0.
   AED_REAL :: dummy_eps, dummy_tau
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. ALLOCATED(data%id_ss) ) RETURN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   IF ( data%resuspension  > 0) THEN
      bottom_stress = _STATE_VAR_S_(data%id_taub)
      bottom_stress = MIN(bottom_stress, 100.)
      _DIAG_VAR_S_(data%id_d_taub) = bottom_stress
      _DIAG_VAR_S_(data%id_resus) = zero_
   ENDIF

   IF ( data%resuspension == 2 .AND. data%id_l_bot > 0 ) &
      _DIAG_VAR_S_(data%id_tau_0) = data%tau_0(1) + data%kTau_0 * _STATE_VAR_S_(data%id_l_bot)


   DO i=1,ubound(data%id_ss,1)

      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(i))

      ! Resuspension
      IF ( data%resuspension > 0 ) THEN

         !IF ( data%resuspension == 2 .AND. i == 1 ) THEN
         !   dummy_tau = _DIAG_VAR_S_(data%id_tau_0)
         !ELSE
         !   dummy_tau =  data%tau_0(i)
         !ENDIF

         IF ( data%resuspension == 2 ) THEN
            dummy_tau = data%tau_0(i) + data%kTau_0 * _STATE_VAR_S_(data%id_l_bot)
            dummy_eps = data%epsilon * _DIAG_VAR_S_(data%id_sfss(i))
         ELSE
            dummy_tau = data%tau_0(i)
            dummy_eps = data%epsilon * data%fs(i)
         ENDIF

         IF ( bottom_stress > dummy_tau ) THEN
            resus_flux = dummy_eps * (bottom_stress - dummy_tau) / data%tau_r
         ELSE
            resus_flux = zero_
         ENDIF
         _DIAG_VAR_S_(data%id_resus) = _DIAG_VAR_S_(data%id_resus) + resus_flux
      ENDIF

      ! Sediment "flux" (not sedimentation) dependent on temperature only.
      ss_flux = data%Fsed(i) * (theta_sed_ss**(temp-20.0))

      ! Transfer sediment flux value to model.
      _FLUX_VAR_(data%id_ss(i)) = _FLUX_VAR_(data%id_ss(i)) + ss_flux + resus_flux

   ENDDO

END SUBROUTINE aed2_calculate_benthic_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction_tracer(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
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
   DO ss_i=1,ubound(data%id_ss,1)
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(ss_i))

      ! Self-shading with contribution from background tracer concentration.
      extinction = extinction + (data%Ke_ss(ss_i)*ss)
   ENDDO
END SUBROUTINE aed2_light_extinction_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_mobility_tracer(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement velocities (+ve up; -ve down)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_tracer_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: vvel
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: temp
!
!-------------------------------------------------------------------------------
!BEGIN

   DO i=1,data%num_tracers
      SELECT CASE (data%settling)

         CASE ( _MOB_CONST_ )
            ! constant settling velocity using user provided value
            vvel = data%w_ss

         CASE ( _MOB_TEMP_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_rho)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%w_ss*mu20*pw / ( mu*pw20 )

         CASE ( _MOB_STOKES_ )
            ! settling velocity based on Stokes Law calculation and cell density
            pw = _STATE_VAR_(data%id_rho)              ! water density
            mu = water_viscosity(temp)                 ! water dynamic viscosity
            IF( data%id_rho(phy_i)>0 ) THEN
              rho_p = _STATE_VAR_(data%id_rho(phy_i))  ! cell density
            ELSE
              rho_p = data%phytos(phy_i)%rho_phy
            ENDIF
            vvel = -9.807*(data%phytos(phy_i)%d_phy**2.)*( rho_p-pw ) / ( 18.*mu )

         CASE DEFAULT
            ! unknown settling/migration option selection
            vvel = data%w_ss

      END SELECT
      ! set global mobility array
      mobility(data%id_ss(i)) = vvel
   ENDDO

END SUBROUTINE aed2_mobility_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_tracer
