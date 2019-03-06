!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ____     | || |  ____  ____  | || |  ____  ____  | |         !
!         | |   .'    `.   | || | |_  _||_  _| | || | |_  _||_  _| | |         !
!         | |  /  .--.  \  | || |   \ \  / /   | || |   \ \  / /   | |         !
!         | |  | |    | |  | || |    > `' <    | || |    \ \/ /    | |         !
!         | |  \  `--'  /  | || |  _/ /'`\ \_  | || |    _|  |_    | |         !
!         | |   `.____.'   | || | |____||____| | || |   |______|   | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################
!#                                                                             #
!# aed2_oxygen.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2018 -  The University of Western Australia               #
!#                                                                             #
!#   GLM is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   GLM is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created May 2011                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"

MODULE aed2_oxygen
!-------------------------------------------------------------------------------
! aed2_oxygen --- oxygen biogeochemical model
!
! The AED module oxygen contains equations that describe exchange of
! oxygen across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_util,  ONLY: aed2_gas_piston_velocity, aed2_oxygen_sat

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_oxygen_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_oxygen_data_t
      !# Variable identifiers
      INTEGER  :: id_oxy
      INTEGER  :: id_temp, id_salt
      INTEGER  :: id_wind
      INTEGER  :: id_Fsed_oxy
      INTEGER  :: id_oxy_sat !, id_atm_oxy_exch3d
      INTEGER  :: id_atm_oxy_exch
      INTEGER  :: id_sed_oxy
      INTEGER  :: id_larea, id_lht, id_cell_vel
      INTEGER  :: oxy_piston_model

      !# Model parameters
      AED_REAL :: Fsed_oxy,Ksed_oxy,theta_sed_oxy
      LOGICAL  :: use_sed_model

     CONTAINS
         PROCEDURE :: define            => aed2_define_oxygen
         PROCEDURE :: calculate_surface => aed2_calculate_surface_oxygen
         PROCEDURE :: calculate         => aed2_calculate_oxygen
         PROCEDURE :: calculate_benthic => aed2_calculate_benthic_oxygen
!        PROCEDURE :: mobility          => aed2_mobility_oxygen
!        PROCEDURE :: light_extinction  => aed2_light_extinction_oxygen
         PROCEDURE :: particle_bgc      => aed2_particle_bgc_oxygen
!        PROCEDURE :: delete            => aed2_delete_oxygen

   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_oxygen(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the aed2_oxygen model
!
!  Here, the oxygen namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_oxygen_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: oxy_piston_model =1
   AED_REAL :: oxy_initial=300.
   AED_REAL :: oxy_min=0.
   AED_REAL :: oxy_max=nan_
   AED_REAL :: Fsed_oxy = -20.0
   AED_REAL :: Ksed_oxy = 30.0
   AED_REAL :: theta_sed_oxy = 1.0
   CHARACTER(len=64) :: Fsed_oxy_variable=''

   NAMELIST /aed2_oxygen/ oxy_initial,oxy_min,oxy_max,Fsed_oxy,Ksed_oxy,theta_sed_oxy,  &
                         Fsed_oxy_variable,oxy_piston_model
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_oxygen initialization"

   ! Read the namelist
   read(namlst,nml=aed2_oxygen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_oxygen'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   data%Ksed_oxy = Ksed_oxy
   data%Fsed_oxy = Fsed_oxy/secs_per_day
   data%theta_sed_oxy = theta_sed_oxy
   data%use_sed_model = Fsed_oxy_variable .NE. ''
   data%oxy_piston_model = oxy_piston_model

   ! Register state variables
   data%id_oxy = aed2_define_variable('oxy','mmol/m**3','oxygen',   &
                                    oxy_initial,minimum=oxy_min,maximum=oxy_max)

   ! Register link to external pools

   IF (data%use_sed_model) data%id_Fsed_oxy = aed2_locate_global_sheet(Fsed_oxy_variable)

   ! Register diagnostic variables
   data%id_sed_oxy = aed2_define_sheet_diag_variable(        &
                     'sed_oxy', 'mmol/m**2/d', 'O2 exchange across sed/water interface')

   data%id_atm_oxy_exch = aed2_define_sheet_diag_variable(   &
                     'atm_oxy_flux', 'mmol/m**2/d', 'O2 exchange across atm/water interface')

!  data%id_atm_oxy_exch3d = aed2_define_sheet_diag_variable( &
!                    'atm_oxy_exch3d', 'mmol/m**2/d', 'Oxygen exchange across atm/water interface')

   data%id_oxy_sat = aed2_define_diag_variable(              &
                     'sat', '%', 'oxygen saturation')

   ! Register environmental dependencies
   data%id_temp = aed2_locate_global('temperature') ! Temperature (degrees Celsius)
   data%id_salt = aed2_locate_global('salinity') ! Salinity (psu)
!  data%id_pres = aed2_locate_global_sheet('pressure') ! Pressure (dbar = 10 kPa)
   data%id_wind = aed2_locate_global_sheet('wind_speed') ! Wind speed at 10 m above surface (m/s)
   data%id_larea = aed2_locate_global_sheet('layer_area')
   data%id_lht = aed2_locate_global('layer_ht')
   IF( oxy_piston_model>3 )data%id_cell_vel= aed2_locate_global('cell_vel')! needed for k600

END SUBROUTINE aed2_define_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface_oxygen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-water exchange for the aed oxygen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_oxygen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, vel, depth

   ! State
   AED_REAL :: oxy

   ! Temporary variables
   AED_REAL :: oxy_atm_flux = zero_  ! Surface atm flux of O2
   AED_REAL :: Coxy_air = zero_      ! Dissolved oxygen in the air phase
   AED_REAL :: koxy_trans = zero_    ! k600 for O2
   AED_REAL :: windHt  = 10.0        ! Height of U10 sensor
   AED_REAL :: f_pres  = 1.0         ! Pressure correction, only applicable at high altitudes
!
!-------------------------------------------------------------------------------
!BEGIN
   !Get dependent state variables from physical driver
   temp  = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt  = _STATE_VAR_(data%id_salt)    ! Salinity (psu)
   wind  = _STATE_VAR_S_(data%id_wind)  ! Wind speed at 10 m above surface (m/s)
   windHt= 10.
   depth = MAX( _STATE_VAR_(data%id_lht), one_ )
   IF (data%id_cell_vel > 0 ) THEN
     vel = _STATE_VAR_(data%id_cell_vel)
   ELSE
     vel = 0.0001
   ENDIF

    ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! Concentration of oxygen in surface layer

  !koxy_trans = aed2_gas_piston_velocity(windHt,wind,temp,salt)
   koxy_trans = aed2_gas_piston_velocity(windHt,wind,temp,salt,               &
                                         vel=vel,                             &
                                         depth=depth,                         &
                                         schmidt_model=2,                     &
                                         piston_model=data%oxy_piston_model)


   ! First get the oxygen concentration in the air phase at the interface
   ! (taken from Riley and Skirrow, 1974)
   f_pres = 1.0
   Coxy_air = f_pres * aed2_oxygen_sat(salt,temp)

   ! Get the oxygen flux
   oxy_atm_flux = koxy_trans * (Coxy_air - oxy)

   ! Transfer surface exchange value to AED2 (mmmol/m2/s) converted by driver
   _FLUX_VAR_T_(data%id_oxy) = oxy_atm_flux

   ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2/day)
   _DIAG_VAR_S_(data%id_atm_oxy_exch) = oxy_atm_flux * secs_per_day
   _DIAG_VAR_(data%id_oxy_sat) =  Coxy_air
   
END SUBROUTINE aed2_calculate_surface_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_oxygen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_oxygen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS(aed2_oxygen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: oxy, temp, salt
   AED_REAL :: diff_oxy, f_pres, coxy_sat

!-------------------------------------------------------------------------------
!BEGIN
   ! Get dependent state variables from physical driver
   temp = _STATE_VAR_(data%id_temp)    ! Temperature (degrees Celsius)
   salt = _STATE_VAR_(data%id_salt)    ! Salinity (psu)

   ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen

   ! Set temporal derivatives
   diff_oxy = 0.

   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (diff_oxy)


   ! Compute the oxygen saturation for diagnostic output
   f_pres = 1.0
   coxy_sat = f_pres * aed2_oxygen_sat(salt,temp)

   ! Export diagnostic variables
   _DIAG_VAR_(data%id_oxy_sat) =  (oxy/coxy_sat)*100.

   ! If an externally maintained pool is present, change the pool according
END SUBROUTINE aed2_calculate_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic_oxygen(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED oxygen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_oxygen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: oxy

   ! Temporary variables
   AED_REAL :: oxy_flux, Fsed_oxy
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   oxy = _STATE_VAR_(data%id_oxy)! oxygen

   IF (data%use_sed_model) THEN
       !Fsed_oxy = _STATE_VAR_S_(data%id_Fsed_oxy)
       Fsed_oxy = _DIAG_VAR_S_(data%id_Fsed_oxy)
   ELSE
       Fsed_oxy = data%Fsed_oxy
   ENDIF

    ! Sediment flux dependent on oxygen and temperature
   oxy_flux = Fsed_oxy * MIN(3.,oxy/(data%Ksed_oxy+oxy) * (data%theta_sed_oxy**(temp-20.0)))

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to AED2.
   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + (oxy_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_oxy) = _FLUX_VAR_B_(data%id_ben_oxy) + (-oxy_flux)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_oxy) = oxy_flux * secs_per_day

END SUBROUTINE aed2_calculate_benthic_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_particle_bgc_oxygen(data,column,layer_idx,ppid,partcl)
!ARGUMENTS
   CLASS (aed2_oxygen_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   AED_REAL,DIMENSION(:),INTENT(inout) :: partcl
!
!LOCALS
   AED_REAL k, mass, pod
!
!-------------------------------------------------------------------------------
!BEGIN
!print*,"aed2_particle_bgc_oxygen ", trim(data%aed2_model_name)
   k = 1   !* for now

   mass = partcl(15) ! I think this is mass [CAB]
   pod = mass * k

!  _STATE_VAR_(data%id_oxy) = pod / ( _STATE_VAR_(data%id_lht) * _STATE_VAR_S_(data%id_larea) )

!  _DIAG_VAR_(data%id_ptm_n) = -pod * (???)

!  ! work out oxygen consumption by an individual particle, if its active
!  IF(ptm%group(g)%particle(p)%istat ==1)
!     oxy_decay = ptm%group(g)%particle(p)%mass * data%coef_something

!  ! add to generic flux doo-daddy, for returning
!  flux_part (id_oxy) = flux_part (id_oxy) + oxy_decay

END SUBROUTINE aed2_particle_bgc_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_oxygen
