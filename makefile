# ======================================== #

# Directory for SSHT
SSHTDIR	= ${SSHT}
# Directory for FFTW
FFTWDIR	= ${FFTW}
# Directory for MATLAB
MLAB	=  /Applications/MATLAB_R2011b.app
# Directory for DOXYGEN
DOXYGEN_PATH=/Applications/Doxygen.app/Contents/Resources/doxygen

# Compiler and options
CC	= gcc
OPT	= -Wall -O3 -g -DS2LET_VERSION=\"1.0\" -DS2LET_BUILD=\"`svnversion -n .`\"
UNAME := $(shell uname)

# ======================================== #

# === MATLAB ===
ifeq ($(UNAME), Linux)
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib
  MEXEXT	= mexa64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif
ifeq ($(UNAME), Darwin)
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib
  MEXEXT	= mexmaci64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif

# === S2LET ===
S2LETDIR = .
S2LETLIB = $(S2LETDIR)/lib/c
S2LETINC = $(S2LETDIR)/include/c
S2LETBIN = $(S2LETDIR)/bin/c
S2LETLIBN= s2let
S2LETSRC = $(S2LETDIR)/src/c
S2LETOBJ = $(S2LETSRC)

# === SSHT ===
SSHTLIB	= $(SSHTDIR)/lib/c
SSHTINC	= $(SSHTDIR)/include/c
SSHTLIBN= ssht

# === FFTW ===
FFTWINC	    = $(FFTWDIR)/include
FFTWLIB     = $(FFTWDIR)/lib
FFTWLIBNM   = fftw3

# ======================================== #

S2LETSRCMAT	= $(S2LETDIR)/src/matlab
S2LETOBJMAT  = $(S2LETSRCMAT)
S2LETOBJMEX  = $(S2LETSRCMAT)

vpath %.c $(S2LETSRC)
vpath %.h $(S2LETSRC)
vpath %_mex.c $(S2LETSRCMAT)

LDFLAGS = -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(S2LETLIB) -l$(S2LETLIBN) -lm

LDFLAGSMEX = -I/usr/local/include -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(S2LETLIB) -l$(S2LETLIBN)

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC) -I$(S2LETINC)

S2LETOBJS= $(S2LETOBJ)/s2let_tilling.o	\
	  $(S2LETOBJ)/s2let_axisym.o 		\
	  $(S2LETOBJ)/s2let_math.o

S2LETOBJSMAT = $(S2LETOBJMAT)/s2let_axisym_tilling_mex.o	\
	  $(S2LETOBJMAT)/s2let_axisym_analysis_mex.o		\
	  $(S2LETOBJMAT)/s2let_axisym_synthesis_mex.o	

S2LETOBJSMEX = $(S2LETOBJMEX)/s2let_axisym_tilling_mex.$(MEXEXT)	\
	  $(S2LETOBJMEX)/s2let_axisym_analysis_mex.$(MEXEXT)	\
	  $(S2LETOBJMEX)/s2let_axisym_synthesis_mex.$(MEXEXT)

$(S2LETOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(S2LETOBJMAT)/%_mex.o: %_mex.c $(S2LETLIB)/lib$(S2LETLIBN).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC} 

$(S2LETOBJMEX)/%_mex.$(MEXEXT): $(S2LETOBJMAT)/%_mex.o $(S2LETLIB)/lib$(S2LETLIBN).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

# ======================================== #

.PHONY: default
default: lib test about tidy

.PHONY: matlab
matlab: $(S2LETOBJSMEX)

.PHONY: all
all: lib matlab doc test about tidy

.PHONY: lib
lib: $(S2LETLIB)/lib$(S2LETLIBN).a
$(S2LETLIB)/lib$(S2LETLIBN).a: $(S2LETOBJS)
	ar -r $(S2LETLIB)/lib$(S2LETLIBN).a $(S2LETOBJS)

.PHONY: lib test
test: $(S2LETBIN)/s2let_test
$(S2LETBIN)/s2let_test: $(S2LETOBJ)/s2let_test.o $(S2LETLIB)/lib$(S2LETLIBN).a
	$(CC) $(OPT) $< -o $(S2LETBIN)/s2let_test $(LDFLAGS)
	$(S2LETBIN)/s2let_test

.PHONY: about
about: $(S2LETBIN)/s2let_about
$(S2LETBIN)/s2let_about: $(S2LETOBJ)/s2let_about.o 
	$(CC) $(OPT) $< -o $(S2LETBIN)/s2let_about
	$(S2LETBIN)/s2let_about

.PHONY: doc
doc:
	$(DOXYGEN_PATH) $(S2LETDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(S2LETDIR)/doc/html/*

.PHONY: clean
clean:	tidy cleandoc
	rm -f $(S2LETLIB)/lib$(S2LETLIBN).a
	rm -f $(S2LETOBJMEX)/*_mex.$(MEXEXT)
	rm -f $(S2LETBIN)/s2let_test
	rm -f $(S2LETBIN)/s2let_about

.PHONY: tidy
tidy:
	rm -f $(S2LETOBJ)/*.o
	rm -f $(S2LETOBJMEX)/*.o
	rm -f *~ 

# ======================================== #