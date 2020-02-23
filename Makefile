#--------------------------------------------------
# Make file for MyLibrary
#FC=ifort
FC=gfortran
INSTLDIR=$(HOME)
MODDIR = mod
FDEP=makedepf90
OS = Linux
DEBUG=off

ifneq (,$(findstring arwin,$(shell uname)))
  OS = OSX
endif

ifeq ($(FC),gfortran)
  ifeq ($(OS), OSX)
    LFLAGS+= -I/usr/local/include -L/usr/local/lib
  endif
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS+= -O3 -fopenmp -fPIC
  ifeq ($(DEBUG), on)
    FFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
  endif
  MODOUT=-J$(MODDIR)
  LINT= -fdefault-integer-8
endif

ifeq ($(FC),ifort)
  LFLAGS+= -mkl -lgsl -lz
  FFLAGS+= -O3 -heap-arrays -qopenmp
  #FFLAGS+= -O3 -heap-arrays -fopenmp
  MODOUT=-module $(MODDIR)
  LINT= -i8
  ifeq ($(DEBUG), on)
    FFLAGS+=-check all
  endif
endif

#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
OBJDIR = obj

SRCF90:=$(wildcard $(SRCDIR)/*.f90)
SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJF90:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90))))
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCF90) $(SRCF95)
OBJS= $(OBJF90) $(OBJF95)
OBJS+= $(OBJDIR)/dvode_f90_m.o

SRCS_ALL = $(SRCS)
OBJS_ALL = $(OBJS)

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs libmyfort test
libmyfort: $(OBJS_ALL)
	$(FC) $(FFLAGS) $(DFLAGS) -shared -o libmyfort.so $^ -shared -fPIC $(LFLAGS)
	@echo "#####################################################################################"
	@echo "To complete the installation, do 'make install'."
	@echo "#####################################################################################"

$(OBJDIR)/dvode_f90_m.o:$(SRCDIR)/dvode/dvode_f90_m.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(LINT) $(MODOUT)  -o $@ -c $(SRCDIR)/dvode/dvode_f90_m.f90
$(OBJDIR)/renormalization.o:$(SRCDIR)/renormalization.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(LINT) $(MODOUT)  -o $@ -c $(SRCDIR)/renormalization.f90
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<


test:
	$(FC) test.f90 -o test.exe -I$(INSTLDIR)/include/myfort -lmyfort

dep:
	$(FDEP) $(SRCS_ALL) -b $(OBJDIR)/ > makefile.d
	@echo "obj/dvode_f90_m.o : src/dvode/dvode_f90_m.f90" >> makefile.d
	@echo obj/renormalization.o : src/renormalization.f90 obj/dvode_f90_m.o obj/linear_algebra.o obj/profiler.o >> makefile.d

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
	 ln -sf $(PWD)/$(MODDIR)/$$g $(INSTLDIR)/include/myfort/$$g; \
	done
	@printf "####################################################################\n"
	@printf " make sure that '$(INSTLDIR)' is included in LIBRARY_PATH and LD_LIBRARY_PATH \n"
	@printf "####################################################################\n"

clean:
	@if [ -d $(OBJDIR) ] ; then rm -r $(OBJDIR); fi
	@if [ -d $(MODDIR) ] ; then rm -r $(MODDIR); fi
	rm -f libmyfort.so
#--------------------------------------------------
-include $(wildcard *.d)
