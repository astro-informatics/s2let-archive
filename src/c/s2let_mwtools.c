// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"

/*!
 * Read MW resolution / band-limit parameter from a FITS file.
 *
 * \param[in]  file Filename.
 * \retval int resolution parameter
 */
int read_mw_bandlimit(char* filename)
{
  long     naxes, *naxis, Lread;
  int      status, hdutype, nfound;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) ) 
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );
  
  if (hdutype != BINARY_TBL) 
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) ) 
    printerror( status );
  
  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status) 
       || nfound != naxes ) 
           printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) ) 
    printerror(status);

  if ( fits_close_file(fptr, &status) ) 
    printerror( status );

  return Lread;

}

/*!
 * Read MW map from a FITS file.
 *
 * \param[out]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void read_mw_map(double* f, char* filename, int L)
{
  long     naxes, *naxis, npix, npercol, irow, Lread;
  int      status, hdutype, nfound, anynul;
  float    nulval;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) ) 
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) 
    printerror( status );

  if (hdutype != BINARY_TBL) 
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) ) 
    printerror( status );
  
  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status) 
       || nfound != naxes ) 
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) ) 
    printerror(status);

  if( Lread != L )
  	printf("Attention : read L = %li but you specified L = %i\n",Lread,L);

  npix = L * (2*L-1);

  if ( (npix%naxis[1]) != 0 ) 
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);

  npercol = npix/naxis[1];
  nulval = -1.6375e30;
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval, 
		       &(f[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }

  if ( fits_close_file(fptr, &status) ) 
    printerror( status );

}

/*!
 * Write MW map from a FITS file.
 *
 * \param[in]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void write_mw_map(char* filename, double* f, int L)
{
  fitsfile *fptr;  
  int status, hdutype;
  long firstrow, firstelem;
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  int tfields   = 1;
  char extname[] = "BINTABLE";  
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1E" };
  char *tunit[] = { " " };
 
  long npix = L * (2*L-1);
  status = 0;
  
  if (fits_create_file(&fptr, filename, &status)) 
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
	    __FILE__, __LINE__);
  
  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n", 
	    __FILE__, __LINE__);
 
  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) ) 
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n", 
	    __FILE__, __LINE__);
  
  if ( fits_create_tbl( fptr, BINARY_TBL, npix, tfields, ttype, tform,
			tunit, extname, &status) )
    fprintf(stderr, "%s (%d): Could not create new binary table.\n", 
	    __FILE__, __LINE__);
  
  if (fits_write_key(fptr, TINT, "L", &L,
		     "Resolution parameter", &status))  
    fprintf(stderr, "%s (%d): Could not write L keyword.\n", 
	    __FILE__, __LINE__);
  
  firstrow  = 1; 
  firstelem = 1; 
  
  if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, npix, f, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);
  
  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n", 
	    __FILE__, __LINE__);
  
}


void printerror(int status)
{
 if (status){
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}