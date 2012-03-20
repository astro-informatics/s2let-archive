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
OPT	= -Wall -O3
UNAME := $(shell uname)

# ======================================== #

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

vpath %.c $(S2LETSRC)
vpath %.h $(S2LETSRC)

LDFLAGS = -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(S2LETLIB) -l$(S2LETLIBN) -lm

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC) -I$(S2LETINC)

S2LETOBJS= $(S2LETOBJ)/s2let_tilling.o		\
	  $(S2LETOBJ)/s2let_axisym.o 				\
	  $(S2LETOBJ)/s2let_math.o

$(S2LETOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

# ======================================== #

.PHONY: default
default: lib test tidy

.PHONY: all
all: lib doc test tidy

.PHONY: lib
lib: $(S2LETLIB)/lib$(S2LETLIBN).a
$(S2LETLIB)/lib$(S2LETLIBN).a: $(S2LETOBJS)
	ar -r $(S2LETLIB)/lib$(S2LETLIBN).a $(S2LETOBJS)

.PHONY: test
lib: $(S2LETBIN)/s2let_test
$(S2LETBIN)/s2let_test: $(S2LETOBJ)/s2let_test.o $(S2LETLIB)/lib$(S2LETLIBN).a
	$(CC) $(OPT) $< -o $(S2LETBIN)/s2let_test $(LDFLAGS)
	$(S2LETBIN)/s2let_test


.PHONY: doc
doc:
	$(DOXYGEN_PATH) $(S2LETDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(S2LETDIR)/doc/html/*

.PHONY: clean
clean:	tidy cleandoc
	rm -f $(S2LETLIB)/lib$(S2LETLIBN).a
	rm -f $(S2LETBIN)/s2let_test

.PHONY: tidy
tidy:
	rm -f $(S2LETOBJ)/*.o
	rm -f *~ 

# ======================================== #