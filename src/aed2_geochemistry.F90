!###############################################################################
!#                                                                             #
!# aed2_geochemistry.F90                                                       #
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


MODULE aed2_geochemistry
!-------------------------------------------------------------------------------
! aed2_geochemistry --- geochemistry model
!
! The AED module geochemistry contains equations that describe exchange of
! soluable reactive geochemistry across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_gclib,ONLY : AED_GC_Input,printTables
   USE aed2_gcsolver,ONLY : ConfigEquilibriumSolver, &
                            GetListOfGeochemDiagnostics, &
                            InitialiseGCProperties, &
                            UpdateEquilibration, &
                            returnGCDerivedVector, &
                            simManganRedox, &
                            simArsenicRedox, &
                            simSulfurRedox, &
                            simCarbonRedox, &
                            simIronRedox


   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_geochemistry_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_geochemistry_data_t
      !# Variable identifiers
      INTEGER  :: id_comp(MAX_GC_COMPONENTS), id_mins(MAX_GC_MINERALS)
      INTEGER  :: id_cdep(MAX_GC_COMPONENTS), id_mdep(MAX_GC_MINERALS)
      INTEGER  :: id_ubalchg ,id_pH, id_c_pco2, id_o_oxy
      INTEGER  :: id_temp,id_sal
      INTEGER  :: id_sed_dic
      INTEGER  :: id_gcdiag(MAX_GC_COMPONENTS), id_noncon
      INTEGER  :: id_totC

      !# Model parameters
      INTEGER  :: num_comp, num_mins

      LOGICAL :: component_linked(MAX_GC_COMPONENTS),mineral_linked(MAX_GC_MINERALS)
      LOGICAL :: simEq
      CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: listDissTransVars,listPartTransVars
      AED_REAL,          DIMENSION(:), ALLOCATABLE :: DissComp,PartComp
      AED_REAL :: speciation_dt
      AED_REAL :: Fsed_gch(MAX_GC_COMPONENTS)
      AED_REAL :: w_gch(MAX_GC_MINERALS)

     CONTAINS
         PROCEDURE :: define            => aed2_define_geochemistry
         PROCEDURE :: initialize        => aed2_initialize_geochemistry
         PROCEDURE :: calculate         => aed2_calculate_geochemistry
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_geochemistry
         PROCEDURE :: equilibrate       => aed2_equilibrate_geochemistry
!        PROCEDURE :: mobility          => aed2_mobility_geochemistry
!        PROCEDURE :: light_extinction  => aed2_light_extinction_geochemistry
!        PROCEDURE :: delete            => aed2_delete_geochemistry

   END TYPE

!===============================================================================
CONTAINS







!###############################################################################
SUBROUTINE aed2_define_geochemistry(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status, i
   INTEGER  :: speciation_dt
   INTEGER  :: num_components,num_minerals
   INTEGER  :: nDissTransportables, nPartTransportables
   LOGICAL  :: simEq = .TRUE.
   AED_REAL          :: min
   AED_REAL          :: dis_initial(MAX_GC_COMPONENTS)
   AED_REAL          :: Fsed_gch(MAX_GC_COMPONENTS)
   AED_REAL          :: min_initial(MAX_GC_MINERALS)
   AED_REAL          :: w_gch(MAX_GC_MINERALS)
   AED_REAL          :: pH_initial, pe_initial
   CHARACTER(len=64) :: geochem_file
   CHARACTER(len=64) :: dis_components(MAX_GC_COMPONENTS)
   CHARACTER(len=64) :: component_link(MAX_GC_COMPONENTS)
   CHARACTER(len=64) :: speciesOutput(10)
   CHARACTER(len=64) :: the_minerals(MAX_GC_MINERALS)
   CHARACTER(len=64) :: mineral_link(MAX_GC_MINERALS)
   CHARACTER(len=64) :: carbon_ph_link = 'CAR_pH'
   CHARACTER(len=64) :: carbon_pco2_link = 'CAR_pCO2'
   CHARACTER(len=64) :: oxy_link = 'OXY_oxy'
   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: diagnosticList

   NAMELIST /aed2_geochemistry/ speciation_dt, geochem_file,                             &
                    num_components, dis_components, component_link, Fsed_gch,dis_initial,&
                    num_minerals, the_minerals, mineral_link, w_gch, min_initial,        &
                    carbon_ph_link, pH_initial, speciesOutput, simEq

!-------------------------------------------------------------------------------
!&aed2_geochemistry
!  speciation_dt  = 10
!  geochem_file   = 'geochem_data.dat'
!  num_components = 5,
!  dis_components = 'DIC','Ca','PO4','FeIII','FeII'
!  component_link = 'aed2_carbon_dic','','aed2_phosphorus_frp','',''
!  Fsed_gch       = 0.,0.,0.,0.,0.
!  num_minerals   = 2
!  the_minerals   = 'Calcite', 'FeOH3A'
!  mineral_link   = '',''
!  w_gch          = 0.,0.
!/
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_geochemistry initialization"

   print *,"  WARNING! aed2_geochemistry model is currently under development"
   ! MH:JOBS
   ! remove pe
   ! species outputs
   ! sediment flux
   ! redox

   !----------------------------------------------------------------------------
   ! Initialise variables
   component_link(:) = ''
   data%component_linked(:) = .FALSE.
   data%mineral_linked(:) = .FALSE.
   data%simEq = simEq

   dis_initial = 1.0  ! default, overwritten by namelist
   min_initial = 0.1  ! default, overwritten by namelist
   pH_initial  = 7.5  ! default, overwritten by namelist
   pe_initial  = 8.0

   !----------------------------------------------------------------------------
   ! Read the namelist
   read(namlst,nml=aed2_geochemistry,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_geochemistry'


   data%speciation_dt = speciation_dt  ! Note this is now managed in FV_AED2 or GLM_AED2
   speciesOutput = ''
   speciesOutput(1) = 'NONCON'
   !speciesOutput(1) = 'HCO3-'

   !----------------------------------------------------------------------------
   ! Now load the geochem database

   !CALL aed2_geochem_load_params(data, geochem_file, modelinfo)
   CALL AED_GC_Input(geochem_file)
   !CALL printTables()

   !----------------------------------------------------------------------------
   ! Configure the geochemical solver
   CALL ConfigEquilibriumSolver( num_components,  num_minerals,                &
                                 dis_components(1:num_components),             &
                                 the_minerals(1:num_minerals),                 &
                                 nDissTransportables, nPartTransportables,     &
                                 data%listDissTransVars,data%listPartTransVars)

   data%num_comp = nDissTransportables
   data%num_mins = nPartTransportables

   !----------------------------------------------------------------------------
   ! Initialise the module level geochemical values ready for the registration

   ALLOCATE(data%DissComp(nDissTransportables))
   ALLOCATE(data%PartComp(nPartTransportables))

   data%DissComp = zero_
   DO i=1,nDissTransportables
     data%DissComp(i) = dis_initial(i)
     data%Fsed_gch(i) = Fsed_gch(i) / secs_per_day
   END DO
   component_link(num_components+1) = carbon_ph_link  ! Special pH var
   data%DissComp(num_components+1) = pH_initial
   data%DissComp(num_components+2) = pe_initial

   data%PartComp = zero_
   DO i=1,num_minerals
     data%PartComp(i) = min_initial(i)
     data%w_gch(i) = w_gch(i)
   END DO

   CALL InitialiseGCProperties(data%DissComp, data%PartComp, 2)

   print *,'data%DissComp',data%DissComp
   print *,'data%PartComp',data%PartComp

   !----------------------------------------------------------------------------
   ! Process dis components adding as state vars or dependancies as appropriate
   DO i=1,nDissTransportables
      !print *,'i,',i,component_link(i)
      IF ( component_link(i) .EQ. '' ) THEN
         min = zero_
         IF ( TRIM(data%listDissTransVars(i)) == 'ubalchg' ) min=nan_
         ! Register state variables
         data%id_comp(i) = aed2_define_variable(                               &
          !                          TRIM(dis_components(i)),                  &
                                    TRIM(data%listDissTransVars(i)),           &
                                    'mmol/m**3','geochemistry',                &
                                    data%DissComp(i),                          &
                                    minimum=min)
      ELSE
         ! Register external state variable dependencies
         data%id_cdep(i) = aed2_locate_variable( TRIM(component_link(i)) )
         data%component_linked(i) = .true.
      ENDIF
   ENDDO

   !----------------------------------------------------------------------------
   ! Process minerals adding as state vars or dependancies as appropriate
   DO i=1,num_minerals
      IF ( mineral_link(i) .EQ. '' ) THEN
         ! Register state variables
         data%id_mins(i) = aed2_define_variable(                               &
                                    TRIM(the_minerals(i)),                     &
                                    'mmol/m**3','geochemistry',                &
                                    data%PartComp(i),                          &
                                    minimum=zero_,                             &
                                    mobility=data%w_gch(i))
      ELSE
         ! Register external state variable dependencies
         data%id_mdep(i) = aed2_locate_variable( mineral_link(i))
         data%mineral_linked(i) = .true.
      ENDIF
   ENDDO


   !----------------------------------------------------------------------------
   ! Register diagnostic variables

   CALL GetListOfGeochemDiagnostics(speciesOutput,diagnosticList)

   DO i=1,SIZE(diagnosticList)
     data%id_gcdiag(i) = aed2_define_diag_variable( diagnosticList(i), &
                         '?mmol/m**3?', 'Geochemistry Diagnostic')
   END DO

   !MH solution to get aed_carbon's pCO2 updated for atm exchange ...
   data%id_c_pco2 = aed2_locate_global(carbon_pco2_link)

   data%id_noncon = aed2_define_diag_variable( 'noncon_mh', &
                         'niter', 'non-convergence status')

   data%id_o_oxy = aed2_locate_global(oxy_link)

   data%id_pH = aed2_locate_global(carbon_ph_link)
   !----------------------------------------------------------------------------

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global( 'temperature' )
   data%id_sal = aed2_locate_global( 'salinity' )


   !----------------------------------------------------------------------------

END SUBROUTINE aed2_define_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed2_initialize_geochemistry(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: temp, tss

   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
!-------------------------------------------------------------------------------
!BEGIN

   !-- Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   !-- Retrieve current (local) state variable values into array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = _STATE_VAR_(data%id_comp(i))
      ELSE
          dissConcs(i) = _STATE_VAR_(data%id_cdep(i))
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = _STATE_VAR_(data%id_mins(i))
      ELSE
          partConcs(i) = _STATE_VAR_(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Redo geochemical equilibration, now spatial initialisation is done
   CALL InitialiseGCProperties(dissConcs, partConcs, 2, inTemp=REAL(temp))

   !-- Copy back into main AED2 arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         _STATE_VAR_(data%id_comp(i)) =  dissConcs(i)
      ELSE
         _STATE_VAR_(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         _STATE_VAR_(data%id_mins(i)) =  partConcs(i)
      ELSE
         _STATE_VAR_(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO


END SUBROUTINE aed2_initialize_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_geochemistry model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: reduction,oxidation

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current (local) state variable values.
!  dic = _STATE_VAR_(data%id_dic)! geochemistry

   reduction = zero_
   oxidation = zero_

   !-- 1. Iron  ---------------------------------------------------------------!
   IF (simIronRedox) THEN
     !TBC
   END IF

!  _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) + (diff_dic)


END SUBROUTINE aed2_calculate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_calculate_benthic_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED geochemistry.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
     AED_REAL :: temp,ph
   ! State
     AED_REAL :: dic,oxy
   ! Parameters
     AED_REAL, PARAMETER :: KpHS = 2.0   ! about half of maximum flux at pH5.
     AED_REAL, PARAMETER :: KDOs = 2.0*1e3/16.
   ! Temporary variables
     AED_REAL :: gch_flux,oxyEffect,pHEffect
     INTEGER :: i

!-------------------------------------------------------------------------------
!BEGIN
   RETURN
   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   ! Retrieve current (local) state variable values.
   ph = _STATE_VAR_(data%id_ph)     ! local pH
   oxy = _STATE_VAR_(data%id_o_oxy)
   !IF (data%use_oxy) oxy = _STATE_VAR_(data%id_oxy)

   oxyEffect = one_
   pHEffect = one_

   DO i=1,data%num_comp

      ! Sediment flux dependent on oxygen
      ! Special hard-coded oxy dependent release FeII
      IF( i==2 ) THEN
        oxyEffect = 1e-8
        oxyEffect = ( KdoS / (KdoS + oxy) )
      ENDIF

      ! Special hard-coded pH dependent release for FeIII and Al
      IF( i==3 .OR. i==4 ) THEN
        pHEffect = 1e-8
        IF( pH<6.0 ) pHEffect = ( abs(ph-7.0) / (KpHS + abs(ph-7.0)) )
      ENDIF

      gch_flux = data%Fsed_gch(i) * 1.05**(temp-20.0) * oxyEffect * pHEffect

      !gch_flux = gch_flux * data%Fsed_gch(i) / (data%Fsed_gch(i) + oxy)

      IF (.NOT.data%component_linked(i)) THEN
        ! geochem module variables mmol/m2/day
         _FLUX_VAR_(data%id_comp(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ELSE
        ! other module variables
        !_FLUX_VAR_(data%id_cdep(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ENDIF
   ENDDO


END SUBROUTINE aed2_calculate_benthic_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed2_equilibrate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, sal
   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
   REAL :: pco2,nc

!-------------------------------------------------------------------------------
!BEGIN

   !-- Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   sal = _STATE_VAR_(data%id_sal) ! local salinity

   !-- Retrieve current gch state variable values into work array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = _STATE_VAR_(data%id_comp(i))
      ELSE
          dissConcs(i) = _STATE_VAR_(data%id_cdep(i))
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = _STATE_VAR_(data%id_mins(i))
      ELSE
          partConcs(i) = _STATE_VAR_(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Do geochemical equilibration
   IF (data%simEq) &
      CALL UpdateEquilibration(dissConcs, partConcs, concMode=2, &
                               inTemp=REAL(temp), inSalt=REAL(sal), &
                               stoEq=.true., upDerv=.true.)


   !-- Copy back into main AED2 arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         _STATE_VAR_(data%id_comp(i)) =  dissConcs(i)
      ELSE
         _STATE_VAR_(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         _STATE_VAR_(data%id_mins(i)) =  partConcs(i)
      ELSE
         _STATE_VAR_(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO

   !-- Update diagnostic arrays
   IF( returnGCDerivedVector("pCO2",pco2) > 0) THEN
      !print *,'pco2: ',pco2,data%id_c_pco2
     _DIAG_VAR_(data%id_c_pco2) = pco2
!     _DIAG_VAR_(data%id_gcdiag(2)) = pco2
   ENDIF
!   IF( returnGCDerivedVector("NONCON",nc) > 0) THEN
!     _DIAG_VAR_(data%id_noncon) = nc
!     _DIAG_VAR_(data%id_gcdiag(3)) = nc
!   ENDIF


END SUBROUTINE aed2_equilibrate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed2_geochemistry
