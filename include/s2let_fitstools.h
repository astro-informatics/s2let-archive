// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_FITSTOOLS
#define S2LET_FITSTOOLS

int read_mw_bandlimit(char* filename);
void read_mw_map(double* f, char* file, int L);
void write_mw_map(char* file, double* f, int L);
void printerror(int status);

#endif
