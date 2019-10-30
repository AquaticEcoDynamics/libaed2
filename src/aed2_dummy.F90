!###############################################################################
!#                                                                             #
!# aed2_dummy.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2019 -  The University of Western Australia               #
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
!# Created July 2019                                                           #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |  ________    | || | ____    ____ | || | ____    ____ | |         !
!         | | |_   ___ `.  | || ||_  |    |  _|| || ||_   \  /   _|| |         !
!         | |   | |   `. \ | || |  | |    | |  | || |  |   \/   |  | |         !
!         | |   | |    | | | || |  | |    , |  | || |  | |\  /| |  | |         !
!         | |  _| |___.' / | || |  | `.__.  ,  | || | _| |_\/_| |_ | |         !
!         | | |________.'  | || |  |_______.   | || ||_____||_____|| |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed2.h"

!
MODULE aed2_dummy
!-------------------------------------------------------------------------------
! aed2_dummy --- dummy model
!
! The AED module "dummy" contains only variables to provide vars required in
! other modules but we dont provide (usually for debugging purposes
!-------------------------------------------------------------------------------
   USE aed2_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed2_dummy_data_t
!
   TYPE,extends(aed2_model_data_t) :: aed2_dummy_data_t
      !# Variable identifiers
      INTEGER  :: num_v, num_dv, num_sv, num_dsv
      INTEGER,ALLOCATABLE :: id_dummy_v(:), id_dummy_dv(:),           &
                             id_dummy_sv(:), id_dummy_dsv(:)

     CONTAINS
         PROCEDURE :: define            => aed2_define_dummy
         PROCEDURE :: calculate         => aed2_calculate_dummy
   END TYPE

! MODULE GLOBALS
   INTEGER :: diag_level = 10

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed2_define_dummy(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed2_dummy_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i, num_v, num_dv, num_sv, num_dsv

   CHARACTER(len=4),POINTER :: prefix => null()

!  %% NAMELIST
   CHARACTER(len=40) :: dm_vars(100)
   AED_REAL          :: dm_max(100)
   AED_REAL          :: dm_min(100)
   AED_REAL          :: dm_init(100)
   CHARACTER(len=40) :: dm_dvars(100)
   CHARACTER(len=40) :: dm_svars(100)
   AED_REAL          :: dm_smax(100)
   AED_REAL          :: dm_smin(100)
   AED_REAL          :: dm_sinit(100)
   CHARACTER(len=40) :: dm_dsvars(100)
!  %% END NAMELIST

   NAMELIST /aed2_dummy/ dm_vars, dm_max, dm_min, dm_init,             &
                         dm_dvars,                                     &
                         dm_svars, dm_smax, dm_smin, dm_sinit,         &
                         dm_dsvars
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed2_dummy initialization"

   dm_vars = ''    ; dm_dvars = ''
   dm_max = NaN_   ; dm_min = NaN_   ; dm_init = 0
   dm_svars = ''   ; dm_dsvars = ''
   dm_smax = NaN_  ; dm_smin = NaN_  ; dm_sinit = 0

   num_v = 0 ; num_dv = 0 ; num_sv = 0 ; num_dsv = 0

   ! Read the namelist
   read(namlst,nml=aed2_dummy,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed2_dummy'

   DO i=1,100 ; IF (dm_vars(i)  .EQ. '' ) THEN  ; num_v  = i-1  ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_dvars(i) .EQ. '' ) THEN  ; num_dv = i-1  ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_svars(i)  .EQ. '' ) THEN ; num_sv  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_dsvars(i) .EQ. '' ) THEN ; num_dsv = i-1 ; EXIT ; ENDIF ; ENDDO

   ALLOCATE(data%id_dummy_v(num_v))
   ALLOCATE(data%id_dummy_dv(num_dv))
   ALLOCATE(data%id_dummy_sv(num_sv))
   ALLOCATE(data%id_dummy_dsv(num_dsv))

   data%num_v   = num_v
   data%num_dv  = num_dv
   data%num_sv  = num_sv
   data%num_dsv = num_dsv

   CALL aed2_set_prefix(prefix)

   ! Register state variables
   DO i=1,data%num_v
      data%id_dummy_v(i) = aed2_define_variable(dm_vars(i), '', '', dm_init(i), dm_min(i), dm_max(i), 0.)
   ENDDO

   DO i=1,data%num_sv
      data%id_dummy_sv(i) = aed2_define_sheet_variable(dm_svars(i), '', '', dm_sinit(i), dm_smin(i), dm_smax(i), .FALSE.)
   ENDDO

   DO i=1,data%num_dv
      data%id_dummy_dv(i) = aed2_define_diag_variable(dm_dvars(i), '', '')
   ENDDO

   DO i=1,data%num_dsv
      data%id_dummy_dsv(i) = aed2_define_sheet_diag_variable(dm_dsvars(i), '', '', .FALSE.)
   ENDDO

END SUBROUTINE aed2_define_dummy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_dummy(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed2_dummy model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_dummy_data_t),INTENT(in) :: data
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: i, count
   AED_REAL :: val, tot

!-------------------------------------------------------------------------------
!BEGIN
END SUBROUTINE aed2_calculate_dummy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed2_dummy
