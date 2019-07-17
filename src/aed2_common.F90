!###############################################################################
!#                                                                             #
!# aed2_common.F90                                                             #
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
!# Created Aug 2013                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"


!###############################################################################
MODULE aed2_common
!-------------------------------------------------------------------------------
   USE aed2_core

   USE aed2_sedflux
   USE aed2_chlorophylla
   USE aed2_oxygen
   USE ufz_oxygen
   USE aed2_silica
   USE aed2_carbon
   USE aed2_nitrogen
   USE aed2_phosphorus
   USE aed2_organic_matter
   USE aed2_phytoplankton
   USE aed2_zooplankton
   USE aed2_tracer
   USE aed2_noncohesive
   USE aed2_totals
   USE aed2_dummy

   IMPLICIT NONE

   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   PUBLIC aed2_new_model, aed2_define_model, aed2_build_model, aed2_print_version

   !#---------------------------------------------------------------------------

   PUBLIC aed2_calculate, aed2_calculate_surface, aed2_calculate_benthic
   PUBLIC aed2_light_extinction, aed2_delete, aed2_equilibrate
   PUBLIC aed2_initialize, aed2_calculate_riparian, aed2_calculate_dry
   PUBLIC aed2_mobility, aed2_rain_loss, aed2_light_shading
   PUBLIC aed2_bio_drag, aed2_particle_bgc

   !# Re-export these from aed2_core.
   PUBLIC aed2_model_data_t, aed2_variable_t, aed2_column_t
   PUBLIC aed2_init_core, aed2_get_var, aed2_core_status

   PUBLIC zero_, one_, nan_, secs_per_day, misval_

   !#---------------------------------------------------------------------------

   CLASS(aed2_model_data_t), POINTER :: model_list => null()
   CLASS(aed2_model_data_t), POINTER :: last_model => null()

CONTAINS
!===============================================================================


!###############################################################################
FUNCTION aed2_new_model(modelname) RESULT(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
   CHARACTER(len=4) :: prefix
   LOGICAL :: is_plus = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)

   SELECT CASE (modelname)
      CASE ('aed2_sedflux');        prefix = 'SDF'; ALLOCATE(aed2_sedflux_data_t::model)
      CASE ('aed2_chlorophylla');   prefix = 'CHL'; ALLOCATE(aed2_chla_data_t::model)
      CASE ('aed2_oxygen');         prefix = 'OXY'; ALLOCATE(aed2_oxygen_data_t::model)
      CASE ('ufz_oxygen');          prefix = 'UOX'; ALLOCATE(ufz_oxygen_data_t::model)
      CASE ('aed2_silica');         prefix = 'SIL'; ALLOCATE(aed2_silica_data_t::model)
      CASE ('aed2_carbon');         prefix = 'CAR'; ALLOCATE(aed2_carbon_data_t::model)
      CASE ('aed2_nitrogen');       prefix = 'NIT'; ALLOCATE(aed2_nitrogen_data_t::model)
      CASE ('aed2_phosphorus');     prefix = 'PHS'; ALLOCATE(aed2_phosphorus_data_t::model)
      CASE ('aed2_organic_matter'); prefix = 'OGM'; ALLOCATE(aed2_organic_matter_data_t::model)
      CASE ('aed2_phytoplankton');  prefix = 'PHY'; ALLOCATE(aed2_phytoplankton_data_t::model)
      CASE ('aed2_zooplankton');    prefix = 'ZOO'; ALLOCATE(aed2_zooplankton_data_t::model)
      CASE ('aed2_tracer');         prefix = 'TRC'; ALLOCATE(aed2_tracer_data_t::model)
      CASE ('aed2_noncohesive');    prefix = 'NCS'; ALLOCATE(aed2_noncohesive_data_t::model)
      CASE ('aed2_totals');         prefix = 'TOT'; ALLOCATE(aed2_totals_data_t::model)
      CASE ('aed2_dummy');          prefix = 'DUM'; ALLOCATE(aed2_dummy_data_t::model)
      CASE ('aed2_land');           is_plus = .TRUE.
      CASE ('aed2_ass');            is_plus = .TRUE.
      CASE ('aed2_soilbgc');        is_plus = .TRUE.
!     CASE ('aed2_cladophora');     is_plus = .TRUE.
      CASE ('aed2_macrophyte');     is_plus = .TRUE.
      CASE ('aed2_macroalgae');     is_plus = .TRUE.
      CASE ('aed2_iron');           is_plus = .TRUE.
      CASE ('aed2_isotope');        is_plus = .TRUE.
      CASE ('aed2_isotope_c');      is_plus = .TRUE.
      CASE ('aed2_radon');          is_plus = .TRUE.
      CASE ('aed2_sulfur');         is_plus = .TRUE.
      CASE ('aed2_geochemistry');   is_plus = .TRUE.
      CASE ('aed2_seddiagenesis');  is_plus = .TRUE.
      CASE ('aed2_vegetation');     is_plus = .TRUE.
      CASE ('aed2_habitat');        is_plus = .TRUE.
      CASE ('csiro_optical');       is_plus = .TRUE.
      CASE ('aed2_bivalve');        is_plus = .TRUE.
      CASE ('aed2_pathogens');      is_plus = .TRUE.
      CASE ('aed2_test');           is_plus = .TRUE.
      CASE DEFAULT;                 print *,'*** Unknown module ', TRIM(modelname)
   END SELECT

   IF ( is_plus ) &
      print*,"To use ",TRIM(modelname)," you will need aed2+"

   IF (ASSOCIATED(model)) THEN
      model%aed2_model_name = modelname
      model%aed2_model_prefix = prefix

      IF ( .NOT. ASSOCIATED(model_list) ) model_list => model
      IF ( ASSOCIATED(last_model) ) last_model%next => model
      last_model => model
   ENDIF
END FUNCTION aed2_new_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_print_version
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed2 version ", TRIM(AED2_VERSION)
   print*,"    libaed2+ not found "
#ifdef __INTEL_COMPILER
   print*,"    libaed2 built using intel fortran version ", __INTEL_COMPILER
#else
# ifdef __PGI
   print*,"    libaed2 built using pgfortran version ", __PGIC__, '.', __PGIC_MINOR__, '.', __PGIC_PATCHLEVEL__
# else
   print*,"    libaed2 built using gfortran version ", __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
# endif
#endif
END SUBROUTINE aed2_print_version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_build_model(model, namlst, do_prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed2_model_data_t),POINTER :: model
   INTEGER,INTENT(in)      :: namlst
   LOGICAL,INTENT(in)      :: do_prefix
!
!LOCALS
   CHARACTER(len=4),POINTER :: prefix_p
!
!-------------------------------------------------------------------------------
!BEGIN
   cur_model_name = model%aed2_model_name
   IF ( do_prefix ) THEN
      prefix_p => model%aed2_model_prefix
      CALL aed2_set_prefix(prefix_p)
   ENDIF
   CALL model%define(namlst)
   IF ( do_prefix ) THEN
      prefix_p => null()
      CALL aed2_set_prefix(prefix_p)
   ENDIF
   cur_model_name = ''
END SUBROUTINE aed2_build_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_define_model(modelname, namlst)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
   INTEGER,INTENT(in)      :: namlst
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)
   model => aed2_new_model(modelname)
   IF ( ASSOCIATED(model) ) CALL aed2_build_model(model, namlst, .TRUE.)
END SUBROUTINE aed2_define_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# These are wrappers for the individual models.                               #
!#                                                                             #
!###############################################################################


!###############################################################################
SUBROUTINE aed2_initialize(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%initialize(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_initialize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_surface(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_surface(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_surface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_benthic(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_benthic(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_riparian(column, layer_idx, pc_wet)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_riparian(column, layer_idx, pc_wet)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_calculate_dry(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_dry(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_calculate_dry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_equilibrate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%equilibrate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed2_equilibrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed2_validate(column, layer_idx) RESULT(valid)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
   LOGICAL :: valid
!-------------------------------------------------------------------------------
   valid = .TRUE.
   model => model_list
   DO WHILE (ASSOCIATED(model))
      IF (.NOT. model%validate(column, layer_idx)) valid = .FALSE.
      model => model%next
   ENDDO
END FUNCTION aed2_validate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_extinction(column, layer_idx, extinction)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   extinction = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_extinction(column, layer_idx, extinction)
      model => model%next
   ENDDO
END SUBROUTINE aed2_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_mobility(column, layer_idx, mobility)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   !mobility = zero_ !MH leave this as is in case default settling vals provided
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%mobility(column, layer_idx, mobility)
      model => model%next
   ENDDO
END SUBROUTINE aed2_mobility
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_rain_loss(column, layer_idx, infil)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   infil = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%rain_loss(column, layer_idx, infil)
      model => model%next
   ENDDO
END SUBROUTINE aed2_rain_loss
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_light_shading(column, layer_idx, shade_frac)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: shade_frac
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   shade_frac = one_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_shading(column, layer_idx, shade_frac)
      model => model%next
   ENDDO
END SUBROUTINE aed2_light_shading
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_bio_drag(column, layer_idx, drag)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   drag = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%bio_drag(column, layer_idx, drag)
      model => model%next
   ENDDO
END SUBROUTINE aed2_bio_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_particle_bgc(column, layer_idx, ppid, partcl)
!-------------------------------------------------------------------------------
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   AED_REAL,DIMENSION(:),INTENT(inout) :: partcl
!
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   !ppid = 0
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%particle_bgc(column, layer_idx, ppid, partcl)
      model => model%next
   ENDDO
END SUBROUTINE aed2_particle_bgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_delete
!-------------------------------------------------------------------------------
!LOCALS
   CLASS (aed2_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%delete
      model => model%next
   ENDDO
END SUBROUTINE aed2_delete
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed2_common
