! ----------------------------------------------------------- !

SUBROUTINE healpix_forward_real( alm, map, nside, L )

  use healpix_types
  use healpix_modules
  implicit none

  integer(i4b) 	:: nside, nlmax, L, el, em
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: map
  complex(dp), dimension(0:L*L-1) 	  :: alm
  real(dp), allocatable, dimension(:,:) :: dw8 
  real(dp), dimension(2) :: z
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp

  nlmax = L - 1
  npix = nside2npix(nside)

  allocate(dw8(1:2*nside, 1:1))
  allocate(alm_temp(1:1, 0:nlmax, 0:nlmax))
  dw8 = 1.0_dp
  z = 0.0_dp

  call map2alm(nside, nlmax, nlmax, map(0:npix-1), &
       alm_temp(1:1,0:nlmax, 0:nlmax), z , dw8 )

  do el = 0, nlmax
     do em = 0, el
        alm( el * el + el + em ) = alm_temp(1, el, em)
        alm( el * el + el - em ) =  (-1.0)**real(-em) * CONJG( alm_temp(1, el, em) )
     enddo
  enddo

  deallocate(dw8, alm_temp)

END SUBROUTINE healpix_forward_real

! ----------------------------------------------------------- !

SUBROUTINE healpix_inverse_real( map, alm, nside, L )

  use healpix_types
  use healpix_modules
  implicit none

  integer(i4b) 	:: nside, nlmax, L, el, em
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: map
  complex(dp), dimension(0:L*L-1) 	  :: alm
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp

  nlmax = L - 1
  npix = nside2npix(nside)

  allocate(alm_temp(1:1, 0:nlmax, 0:nlmax))

  do el = 0, nlmax
     do em = 0, el
        alm_temp(1, el, em) = alm( el * el + el + em )
     enddo
  enddo

  call alm2map(nside, nlmax, nlmax, alm_temp, map)

  deallocate(alm_temp)

END SUBROUTINE healpix_inverse_real

! ----------------------------------------------------------- !

SUBROUTINE read_healpix_map(map, file, nside)

  use healpix_types
  use healpix_modules

  character(len=filenamelen) :: file
  integer :: nside, npix
  real(dp), dimension(0:12*nside*nside-1,1:1) :: map
  real(dp) :: nullval
  logical :: anynull

  npix = nside2npix(nside)
  call read_bintab(file, map, npix, 1, nullval, anynull)

END SUBROUTINE read_healpix_map

! ----------------------------------------------------------- !

SUBROUTINE write_healpix_map(file, map, nside)

  use healpix_types
  use healpix_modules

  character(len=filenamelen) :: file
  CHARACTER(len=80), DIMENSION(1:120):: header
  integer :: nside, npix
  real(dp), dimension(0:12*nside*nside-1,1:1) :: map

  npix = nside2npix(nside)
  header = ''
  call write_minimal_header(header, 'MAP', nside=nside, ordering='ring')
  call write_bintab(map, npix, 1, header, 120, file)

END SUBROUTINE write_healpix_map

! ----------------------------------------------------------- !
