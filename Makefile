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
OBJDIR = obj

SRCF90:=$(wildcard $(SRCDIR)/*.f90)
SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJF90:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90))))
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCF90) $(SRCF95)
OBJS= $(OBJF90) $(OBJF95)

SRCS_ALL = $(SRCS)
OBJS_ALL = $(OBJS)

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
