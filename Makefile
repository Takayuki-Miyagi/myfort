#--------------------------------------------------
# Make file for MyLibrary
#FC=ifort
FC=gfortran
INSTLDIR=$(HOME)
MODDIR = mod
FDEP=makedepf90
OS = Linux
DEBUG=on
ifneq (,$(findstring arwin,$(shell uname)))
  OS = OSX
endif

ifeq ($(FC),gfortran)
  ifeq ($(OS), OSX)
    LFLAGS+= -I/usr/local/include -L/usr/local/lib
  endif
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS+= -O3 -fopenmp
  ifeq ($(DEBUG), on)
    FFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
  endif
  MODOUT=-J$(MODDIR)
endif

ifeq ($(FC),ifort)
  LFLAGS+= -mkl -lgsl -lz
  FFLAGS+= -O3 -heap-arrays -qopenmp
  #FFLAGS+= -O3 -heap-arrays -fopenmp
  MODOUT=-module $(MODDIR)
endif

#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
#SRCDIR_LinAlg = submodules/LinAlgf90/src
#SRCDIR_spline = submodules/NdSpline/src
#SRCDIR_2Body = src/TwoBody
#SRCDIR_3Body = src/ThreeBody
#SRCDIR_ABody = src/ABody
#DEPDIR = .
OBJDIR = obj

#MODDIR = $(SRCDIR)
#MODDIR = $(SRCDIR_2Body)
#MODDIR = $(SRCDIR_3Body)
#MODDIR = $(SRCDIR_ABody)

SRCC:=$(wildcard $(SRCDIR)/*.c)
SRCF77:=$(wildcard $(SRCDIR)/*.f)
SRCF90:=$(wildcard $(SRCDIR)/*.f90)
SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJC:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC))))
OBJF77:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77))))
OBJF90:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90))))
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCC) $(SRCF77) $(SRCF90) $(SRCF95)
OBJS= $(OBJC) $(OBJF77) $(OBJF90) $(OBJF95)

#SRCC_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.c)
#SRCF77_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f)
#SRCF90_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f90)
#SRCF95_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.F90)
#OBJC_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_LinAlg))))
#OBJF77_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_LinAlg))))
#OBJF90_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_LinAlg))))
#OBJF95_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_LinAlg))))
#SRCS_LinAlg= $(SRCC_LinAlg) $(SRCF77_LinAlg) $(SRCF90_LinAlg) $(SRCF95_LinAlg)
#OBJS_LinAlg= $(OBJC_LinAlg) $(OBJF77_LinAlg) $(OBJF90_LinAlg) $(OBJF95_LinAlg)
#
#SRCF95_spline:=$(wildcard $(SRCDIR_spline)/*.F90)
#OBJF95_spline:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_spline))))
#SRCS_spline=  $(SRCF95_spline)
#OBJS_spline=  $(OBJF95_spline)
#
#SRCC_2Body:=$(wildcard $(SRCDIR_2Body)/*.c)
#SRCF77_2Body:=$(wildcard $(SRCDIR_2Body)/*.f)
#SRCF90_2Body:=$(wildcard $(SRCDIR_2Body)/*.f90)
#SRCF95_2Body:=$(wildcard $(SRCDIR_2Body)/*.F90)
#OBJC_2Body:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_2Body))))
#OBJF77_2Body:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_2Body))))
#OBJF90_2Body:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_2Body))))
#OBJF95_2Body:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_2Body))))
#SRCS_2Body= $(SRCC_2Body) $(SRCF77_2Body) $(SRCF90_2Body) $(SRCF95_2Body)
#OBJS_2Body= $(OBJC_2Body) $(OBJF77_2Body) $(OBJF90_2Body) $(OBJF95_2Body)
#
#SRCC_3Body:=$(wildcard $(SRCDIR_3Body)/*.c)
#SRCF77_3Body:=$(wildcard $(SRCDIR_3Body)/*.f)
#SRCF90_3Body:=$(wildcard $(SRCDIR_3Body)/*.f90)
#SRCF95_3Body:=$(wildcard $(SRCDIR_3Body)/*.F90)
#OBJC_3Body:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_3Body))))
#OBJF77_3Body:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_3Body))))
#OBJF90_3Body:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_3Body))))
#OBJF95_3Body:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_3Body))))
#SRCS_3Body= $(SRCC_3Body) $(SRCF77_3Body) $(SRCF90_3Body) $(SRCF95_3Body)
#OBJS_3Body= $(OBJC_3Body) $(OBJF77_3Body) $(OBJF90_3Body) $(OBJF95_3Body)
#
#SRCC_ABody:=$(wildcard $(SRCDIR_ABody)/*.c)
#SRCF77_ABody:=$(wildcard $(SRCDIR_ABody)/*.f)
#SRCF90_ABody:=$(wildcard $(SRCDIR_ABody)/*.f90)
#SRCF95_ABody:=$(wildcard $(SRCDIR_ABody)/*.F90)
#OBJC_ABody:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_ABody))))
#OBJF77_ABody:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_ABody))))
#OBJF90_ABody:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_ABody))))
#OBJF95_ABody:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_ABody))))
#SRCS_ABody= $(SRCC_ABody) $(SRCF77_ABody) $(SRCF90_ABody) $(SRCF95_ABody)
#OBJS_ABody= $(OBJC_ABody) $(OBJF77_ABody) $(OBJF90_ABody) $(OBJF95_ABody)


SRCS_ALL = $(SRCS) $(SRCS_LinAlg) $(SRCS_spline) $(SRCS_2Body) $(SRCS_3Body) $(SRCS_ABody)
OBJS_ALL = $(OBJS) $(OBJS_LinAlg) $(OBJS_spline) $(OBJS_2Body) $(OBJS_3Body) $(OBJS_ABody)

#$(info $(SRCS_ALL))
#$(info )
#$(info $(OBJS_ALL))

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs libmyfort
libmyfort: $(OBJS_ALL)
	$(FC) $(FFLAGS) $(DFLAGS) -shared -fpic -o libmyfort.so $^ $(LFLAGS)
	@echo "#####################################################################################"
	@echo "To complete the installation, do 'make install'."
	@echo "#####################################################################################"

$(OBJDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CFLAGS)  -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(MODOUT)  -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<

test:
	$(FC) test.f90 -o test.exe -I$(INSTLDIR)/include -lmyfort

dep:
	$(FDEP) $(SRCS_ALL) -b $(OBJDIR)/ > makefile.d

dirs:
	if test -d $(OBJDIR); then \
		: ; \
	else \
		mkdir $(OBJDIR); \
	fi
	if test -d $(MODDIR); then \
		: ; \
	else \
		mkdir $(MODDIR); \
	fi

install:
	@if [ ! -d $(INSTLDIR)/lib ] ; then \
	  mkdir $(INSTLDIR)/lib; \
	fi
	@if [ ! -d $(INSTLDIR)/include ] ; then \
	  mkdir $(INSTLDIR)/include; \
	fi
	ln -sf $(PWD)/libmyfort.so $(INSTLDIR)/lib/libmyfort.so
	@for f in $(MODDIR)/*.mod; do \
	 g=$$(basename $$f); \
	 ln -sf $(PWD)/$(MODDIR)/$$g $(INSTLDIR)/include/$$g; \
	done
	@printf "####################################################################\n"
	@printf " make sure that '$(INSTLDIR)' is included in LIBRARY_PATH and LD_LIBRARY_PATH \n"
	@printf "####################################################################\n"

clean:
	rm -r obj
	rm -r mod
	rm -f libmyfort.so
	rm -f $(PWD)/*.exe
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS_ALL)

#--------------------------------------------------
-include $(wildcard *.d)
