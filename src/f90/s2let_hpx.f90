! ----------------------------------------------------------- !

SUBROUTINE s2let_hpx_map2alm( map, alm, nside, L )
  use healpix_types
  use healpix_modules
  implicit none
  integer(i4b) 	:: nside, nlmax, L, npix, el, em
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

END SUBROUTINE s2let_hpx_map2alm

! ----------------------------------------------------------- !

SUBROUTINE s2let_hpx_alm2map( alm, map, nside, L )
  use healpix_types
  use healpix_modules
  implicit none
  integer(i4b) 	:: nside, nlmax, L, npix, el, em
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

END SUBROUTINE s2let_hpx_alm2map

! ----------------------------------------------------------- !

SUBROUTINE  s2let_hpx_read_map(map, file, nside)
  use healpix_types
  use healpix_modules
  character(len=filenamelen) :: file
  integer :: nside, npix
  real(dp), dimension(0:12*nside*nside-1,1:1) :: map
  real(dp) :: nullval
  logical :: anynull

  npix = nside2npix(nside)
  call read_bintab(file, map, npix, 1, nullval, anynull)

END SUBROUTINE s2let_hpx_read_map

! ----------------------------------------------------------- !

SUBROUTINE s2let_hpx_write_map(file, map, nside)
  use healpix_types
  use healpix_modules
  character(len=filenamelen) :: file
  CHARACTER(len=80), DIMENSION(1:120):: header
  integer :: nside, npix
  real(dp), dimension(0:12*nside*nside-1,1:1) :: map

  npix = nside2npix(nside)
  header = ''
  call write_bintab(map, npix, 1, header, 120, file)

END SUBROUTINE s2let_hpx_write_map

! ----------------------------------------------------------- !
