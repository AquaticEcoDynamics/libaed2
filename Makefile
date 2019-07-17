###############################################################################
#                                                                             #
# Makefile to build the aed2 water quality library                            #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013 - 2018 -  The University of Western Australia               #
#                                                                             #
#   GLM is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GLM is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

srcdir=src
incdir=include

ifeq ($(SINGLE),true)
  TARGET=lib/libaed2_s.a
  objdir=obj_s
  moddir=mod_s
else
  TARGET=lib/libaed2.a
  objdir=obj
  moddir=mod
endif

INCLUDES=-I${incdir}

#ifneq ("$(wildcard ../libaed2-plus/Makefile)","")
#  HAVEPLUS=-DHAVE_PLUS
#endif

ifeq ($(F90),ifort)
  INCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback
  OPT_FFLAGS=-O3
  FFLAGS=-fPIC -warn all -module ${moddir} -i-static -mp1 -stand f08 -warn nounused $(DEFINES) $(INCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check
  endif
  ifeq ($(SINGLE),true)
    FFLAGS+=-real-size 32
  else
    FFLAGS+=-real-size 64
  endif
else ifeq ($(F90),pgfortran)
  DEBUG_FFLAGS=-g
  OPT_FFLAGS=-O3
  FFLAGS=-fPIC -module ${moddir} $(DEFINES) $(INCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-r8
else
  DEBUG_FFLAGS=-g -fbacktrace
  OPT_FFLAGS=-O3
  FFLAGS=-fPIC -Wall -J ${moddir} -ffree-line-length-none -std=f2008 $(DEFINES) $(INCLUDES) -fall-intrinsics -Wno-unused -Wno-unused-dummy-argument -fno-range-check -Wno-integer-division
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-fcheck=all
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
endif

ifeq ($(DEBUG),true)
  DEBUG_CFLAGS=-g
  OPT_CFLAGS=
  OPT_FFLAGS=
else
  DEBUG_FFLAGS=
  DEBUG_CFLAGS=
  # OPT_CFLAGS=-O4 -Ofast -frounding-math
  OPT_CFLAGS=-O3
  # OPT_CFLAGS=
  # OPT_FFLAGS=
endif

ifeq ($(SINGLE),true)
  FFLAGS += -DSINGLE=1
endif


FFLAGS+=$(DEBUG_FFLAGS) $(OPT_FFLAGS) $(HAVEPLUS)


OBJS = \
${objdir}/aed2_core.o \
${objdir}/aed2_util.o \
${objdir}/aed2_csv_reader.o \
${objdir}/aed2_sedflux.o \
${objdir}/aed2_chlorophylla.o \
${objdir}/aed2_oxygen.o \
${objdir}/ufz_oxygen.o \
${objdir}/aed2_silica.o \
${objdir}/aed2_carbon.o \
${objdir}/aed2_nitrogen.o \
${objdir}/aed2_phosphorus.o \
${objdir}/aed2_organic_matter.o \
${objdir}/aed2_bio_utils.o \
${objdir}/aed2_phytoplankton.o \
${objdir}/aed2_zoop_utils.o \
${objdir}/aed2_zooplankton.o \
${objdir}/aed2_tracer.o \
${objdir}/aed2_noncohesive.o \
${objdir}/aed2_totals.o \
${objdir}/aed2_dummy.o \
${objdir}/aed2_common.o


all: $(TARGET)


lib:
	@mkdir lib

${moddir}:
	@mkdir ${moddir}

${objdir}:
	@mkdir ${objdir}

${TARGET}: ${objdir} ${moddir} ${OBJS} lib
	ar rv $@ ${OBJS}
	ranlib $@

clean: ${objdir}
	@touch ${objdir}/1.o 1.i90
	@/bin/rm ${objdir}/*.o *.i90

distclean: clean
	@touch lib mod_s mod
	@/bin/rm -rf lib
	@/bin/rm -rf obj obj_s
	@/bin/rm -rf mod mod_s

${objdir}/%.o: ${srcdir}/%.F90 ${srcdir}/aed2_core.F90 ${incdir}/aed2.h
	$(F90) $(FFLAGS) -g -c $< -o $@
