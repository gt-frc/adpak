      program plotit

c/    This program uses the PLPLOT library to plot the results of
c/    the comparison between theory and interpolation from the
c/    ADPAK package. It was written as a standalone program, so it
c/    can operate with data from different computer platforms that
c/    may not have the PLPLOT libraries.
c/    The code reads the file outplt.dat which must be present in
c/    the directory from which we run the code.

      use plplot

      integer nin, npnts, maxpnts
      parameter (maxpnts = 100)
      integer, parameter :: dp = kind(1.0d0)

      real(kind=dp) tekev(maxpnts), Lz_tab(maxpnts), Lz_th(maxpnts)

      data nin /10/

      character*10 dev
      character*2 spec
      character*80 title

c/    Read data:
      open (nin, file = 'outplt.dat', status = 'old')

      npnts = 0
      read (nin, '(37x, i4)') nz
      do i = 1, maxpnts
	 read (nin, *, end = 100) tekev(i), Lz_th(i), Lz_tab(i)
	 npnts = npnts + 1
      enddo

 100  continue

c/    Determine what species we have:

      spec = '??'

      if (nz.eq.2) spec='He'
      if (nz.eq.3) spec='Li'
      if (nz.eq.4) spec='Be'
      if (nz.eq.5) spec='B '
      if (nz.eq.6) spec='C '
      if (nz.eq.7) spec='N '
      if (nz.eq.8) spec='O '
      if (nz.eq.9) spec='F '
      if (nz.eq.10) spec='Ne'
      if (nz.eq.11) spec='Na'
      if (nz.eq.12) spec='Mg'
      if (nz.eq.13) spec='Al'
      if (nz.eq.14) spec='Si'
      if (nz.eq.16) spec='S '
      if (nz.eq.18) spec='Ar'
      if (nz.eq.19) spec='K '
      if (nz.eq.20) spec='Ca'
      if (nz.eq.21) spec='Sc'
      if (nz.eq.22) spec='Ti'
      if (nz.eq.23) spec='V '
      if (nz.eq.24) spec='Cr'
      if (nz.eq.26) spec='Fe'
      if (nz.eq.28) spec='Ni'
      if (nz.eq.29) spec='Cu'
      if (nz.eq.30) spec='Zn'
      if (nz.eq.33) spec='As'
      if (nz.eq.36) spec='Kr'
      if (nz.eq.37) spec='Rb'
      if (nz.eq.40) spec='Zr'
      if (nz.eq.41) spec='Nb'
      if (nz.eq.42) spec='Mo'
      if (nz.eq.45) spec='Rh'
      if (nz.eq.47) spec='Ag'
      if (nz.eq.50) spec='Sn'
      if (nz.eq.54) spec='Xe'
      if (nz.eq.55) spec='Cs'
      if (nz.eq.56) spec='Ba'
      if (nz.eq.64) spec='Gd'
      if (nz.eq.73) spec='Ta'
      if (nz.eq.74) spec='W '
      if (nz.eq.77) spec='Ir'
      if (nz.eq.79) spec='Au'
      if (nz.eq.80) spec='Hg'
      if (nz.eq.83) spec='Bi'
      if (nz.eq.86) spec='Rn'
      if (nz.eq.90) spec='Th'
      if (nz.eq.92) spec='U '

c/    Write title:

      write (title, 5000) spec

c/    Determine scales (maximum always 100 keV):

      xmin = alog10(tekev(1))
      xmax = 2.0
      ymax = -1.0e+30
      ymin = 1.0e+30

      do i = 1, npnts
	 ymax = amax1(ymax, Lz_th(i), Lz_tab(i))
	 ymin = amin1(ymin, Lz_th(i), Lz_tab(i))
      enddo

      ymax = alog10(ymax)
      ymin = alog10(ymin)

      ymax = aint(ymax)
      ymin = aint(ymin) - 1

c/    Convert to logarithms:
      do i = 1, npnts
	 tekev(i) = alog10(tekev(i))
	 Lz_th(i) = alog10(Lz_th(i))
	 Lz_tab(i) = alog10(Lz_tab(i))
      enddo

c/    Start plotting program:
c     call plsdev(dev)

      call plinit()

      call plenv (xmin, xmax, ymin, ymax, 0, 30)
      call pllab ('T#de#u (keV)', 'Cooling Rate (ergs cm#u3#d/sec)',
     .   title)
      call plline(npnts, tekev, Lz_th)
      call plpoin(npnts, tekev, Lz_tab, 9)

      call plend

 5000 format ('Cooling Rates for ', a2)

      stop
      end
