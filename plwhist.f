c/    This routine plots results from the WHIST runs using
c/    the PLPLOT package on the SUN Workstation.
c/    Written by John Mandrekas, GIT, 1-19-95

      integer mxp
      parameter (mxp = 201, mxr = 51, maxt = 5)
      real tim(mxp), paux(mxp), prad(mxp), palf(mxp), pout(mxp),
     .   wstored(mxp), tau(mxp), Q(mxp)
      real te0(mxp), ti0(mxp), teav(mxp), tiav(mxp), ne0(mxp),
     .   neav(mxp), niav(mxp), ali(mxp)
      real rho(mxr), ctot(mxr,maxt), jtot(mxr,maxt), qmhd(mxr,maxt),
     .   jbs(mxr,maxt)
      real pfuse(mxr), pconde(mxr), pconve(mxr), pbe(mxr), prfe(mxr),
     .   pie(mxr), pdive(mxr), prade(mxr), pfusi(mxr), pcondi(mxr),
     .   pconvi(mxr), pbi(mxr), prfi(mxr), pdivi(mxr), r(mxr)
      real pfusei(mxr), pcondei(mxr), pconvei(mxr), pbei(mxr),
     .   prfei(mxr), pohmi(mxr), pdivei(mxr), pradei(mxr), pfusii(mxr),
     .   pcondii(mxr), pconvii(mxr), pbii(mxr), prfii(mxr), pdivii(mxr)
      real tekev(mxr), tikev(mxr), n_e(mxr), n_i(mxr), zeff(mxr),
     .   chi_e(mxr), chi_i(mxr), dfc(mxr)
      real rmhd(mxr), elong(mxr), triang(mxr), vol(mxr), vpr(mxr),
     .   gth(mxr), gph(mxr), grho(mxr), F(mxr), shift(mxr)
      integer time(3), date(3)
      real qmax, jmax, cmax, bsmax, zefmax, chimax

      namelist /dbg1/ denmax, tkevmax, zefmax, chimax

      character*20  dummy
      character*10 dev
      character*40 datetime
      character*25 title(5)
      character*80 bigTitle

      data dev/'tek'/, bigTitle/'This is the Main Title'/
      data nin0/3/, ntty/6/, ndbg /10/
      data nin1/11/, nin2/12/, nin3/13/, nin4/14/, nin5/15/, nin6/16/,
     .   nin7/17/, nin8/18/, nin9/19/, nin10/20/, nin11/21/, nin12/22/



      nprof = 5

c/    First, read data from input file:

      open (nin0, file = 'reslts.dat', status = 'unknown')
      open (nin1, file = 'pwrs2d.txt', status = 'old')
      open (nin2, file = 'prof2d.txt', status = 'old')
      open (nin3, file = 'ctot.plt', status = 'old')
      open (nin4, file = 'jtot.plt', status = 'old')
      open (nin5, file = 'qMHD.plt', status = 'old')
      open (nin6, file = 'jbs.plt', status = 'old')
      open (nin7, file = 'elecBal.plt', status = 'old')
      open (nin8, file = 'elecIntBal.plt', status = 'old')
      open (nin9, file = 'ionBal.plt', status = 'old')
      open (nin10, file = 'ionIntBal.plt', status = 'old')
      open (nin11, file = 'plasProfs.plt', status = 'old')
      open (nin12, file = 'metrics.plt', status = 'old')
      open (ndbg, file = 'debug.plots',  status = 'unknown')

c/    Initializations:

      npnts = 0
      nrho = 0
      nr = 0
      nrmhd = 0

c-----Get system data and time:           ! SunOS version
      call itime (time)
      call idate (date)
c-----Write date and time:
      write (datetime, 750) date(2), date(1), date(3), time(1),
     .  time(2), time(3)

c/    Read time-average results from file reslts.txt (if it exists!)

      idat = 0
      read (nin0, '(A80)', err=50) bigTitle
      read (nin0, '(A10)' ,err=50) dev
      do i = 1, 5
        read (nin0, 500,err=50) title(i)
      enddo
      idat = 1
  50  continue

c/    READ DATA FROM INPUT FILES:

      read (nin1,'(1x,a5)') dummy
      read (nin2,'(1x,a5)') dummy
      do i = 1, mxp
	 read (nin1, *, end = 100) tim(i), paux(i), prad(i), palf(i),
     .      pout(i), wstored(i), tau(i), Q(i)
	 npnts = npnts + 1
         read (nin2, *) tdum, te0(i), teav(i), ti0(i), tiav(i), ne0(i),
     .      neav(i), niav(i), ali(i)
      enddo

      if (npnts.gt.mxp) then
	 write (ntty, 1000) npnts, mxp
	 stop
      endif

 100  continue

c/    Read Current profile files:

      do i = 1, mxr
	 read (nin3, *, end = 101) rho(i), (ctot(i,j), j = 1, nprof)
	 read (nin4, *) rdum, (jtot(i,j), j = 1, nprof)
	 read (nin5, *) rdum, (qmhd(i,j), j = 1, nprof)
	 read (nin6, *) rdum, (jbs(i,j), j = 1, nprof)
	 nrho = nrho + 1
      enddo

      if (nrho.gt.mxr) then
	 write (ntty, 2000) nrho, mxr
	 stop
      endif

 101  continue

c/    Read electron and ion energy balance files:

      read (nin7,'(1x,a5)') dummy
      read (nin8,'(1x,a5)') dummy
      read (nin9,'(1x,a5)') dummy
      read (nin10,'(1x,a5)') dummy

      do i = 1, mxr
	 read (nin7, *, end = 102) r(i), pfuse(i), pconde(i), pconve(i),
     .      pbe(i), prfe(i), pie(i), pdive(i), prade(i)
	 read (nin8, *) rdum, pfusei(i), pcondei(i), pconvei(i),
     .      pbei(i), prfei(i), pdivei(i), pradei(i), pohmi(i)
	 read (nin9,*) rdum, pfusi(i), pcondi(i), pconvi(i),
     .      pbi(i), prfi(i), pdivi(i)
	 read (nin10, *) rdum, pfusii(i), pcondii(i), pconvii(i),
     .      pbii(i), prfii(i), pdivii(i)
	 nr = nr + 1
      enddo

      if (nr.gt.mxr) then
	 write (ntty, 3000) nr, mxr
	 stop
      endif

 102  continue

c/    Read file plasProfs.plt with T, n, Zeff and chi profiles:

      read (nin11,'(1x,a5)') dummy

      do i = 1, nr
	 read (nin11,*) rdum, tekev(i), tikev(i), n_e(i),
     .      n_i(i), zeff(i), chi_e(i), chi_i(i), dfc(i)
      enddo
      chi_e(nr) = chi_e(nr-1)
      chi_i(nr) = chi_i(nr-1)
      dfc(nr) = dfc(nr-1)

c/    Read MHD metric quantities:

      read (nin12,'(1x,a5)') dummy

      do i = 1, mxr
	 read (nin12,*, end = 103) rmhd(i), elong(i), triang(i), vol(i),
     .     vpr(i), gth(i), gph(i), grho(i), F(i), shift(i)
	 nrmhd = nrmhd + 1
      enddo

 103  continue

c/    START OF PLOTTING STATEMENTS *************************************

c/    Determine scales for plotting:

      xmin = tim(1)
      xmax = tim(npnts)

      palfmax = -1.0e30
      pradmax = -1.0e30
      denmax = -1.0e30
      tkevmax = -1.0e30

      do i = 1, npnts
	 palfmax = amax1(palf(i), palfmax)
	 pradmax = amax1(prad(i), pradmax)
	 denmax = amax1(ne0(i), denmax)
	 tkevmax = amax1(tkevmax, te0(i), ti0(i))
      enddo

      ymin = 0.0
      pwmax = amax1(palfmax, pradmax)
      pwmax = getmax(pwmax)
      denmax = getmax(denmax)
      tkevmax = getmax(tkevmax)

c/    Initialize plotting program:

c/    Get output device (default is Tektronix terminal)

      call plsdev(dev)

c/    PLOTS OF TIME HISTORY OF GLOBAL QUANTITIES:

      call plinit()

c/    Write title page if idat = 1:
      if (idat.eq.1) then
	dx = 0.0
        call pladv(0)
        call plvsta
        call plwind(0.0, 1.0, 0.0, 1.0)
        call plptex(0.5,0.9,0.1,0.0,0.5,bigTitle)
        call plptex(0.5,0.8,0.1,0.0,0.5,datetime)
        call plptex(0.5,0.7,0.1,0.0,0.5,
     .    'Time Averages of Global Power Quantities')
	do i = 1, 5
          call plptex(0.5,0.5-dx,0.1,0.0,0.5,title(i))
	  dx = dx + 0.1
        enddo
      endif

      call plenv (xmin, xmax, ymin, pwmax, 0, 0)

      call pllab ('Time (s)', 'Powers (MW)', 'Power Balance')

      call pllsty(1)
      call plline(npnts,tim,pout)
      call pllsty(2)
      call plline (npnts, tim, palf)
      call pllsty(3)
      call plline (npnts, tim, prad)
      call pllsty(4)
      call plline (npnts, tim, paux)

c/    Next page:

      call pllsty(1)
      call plenv (xmin, xmax, ymin, tkevmax, 0, 0)
      call pllab ('Time (s)', 'Temperatures (keV)',
     .   'Peak and Average Electron and Ion Temperatures')
      call plline (npnts, tim, te0)
      call plline (npnts, tim, teav)
      call pllsty(2)
      call plline (npnts, tim, ti0)
      call plline (npnts, tim, tiav)

      call pllsty(1)
      call plenv (xmin, xmax, ymin, denmax, 0, 0)
      call pllab ('Time (s)', 'Densities (10#u20#d m#u-3#d)',
     .   'Peak and Average Densities')
      call plline (npnts, tim, ne0)
      call plline (npnts, tim, neav)
      call pllsty(2)
      call plline (npnts, tim, niav)


      call plend

c/    PLOTS OF MHD-RELATED QUANTITIES:

c/    Change units (Currents in MA and densities in MA/m2)
      do i = 1, nrho
	 do j = 1, nprof
	   ctot(i,j) = 1.0e-6*ctot(i,j)
	   jtot(i,j) = 1.0e-2*jtot(i,j)
	   jbs(i,j) = 1.0e-2*jbs(i,j)
         enddo
      enddo

c/    Determine scales for plotting:
      rmin = 0.0
      rmax = 1.0

      qmax = -1.0e30
      cmax = -1.0e30
      jmax = -1.0e30
      bsmax = -1.0e30

      do i = 1, nrho
	 do j = 1, nprof
	   qmax = amax1(qmax, qmhd(i,j))
	   cmax = amax1(cmax, ctot(i,j))
	   jmax = amax1(jmax, jtot(i,j))
	   bsmax = amax1(bsmax, jbs(i,j))
         enddo
      enddo

      qmax = getmax(qmax)
      cmax = getmax(cmax)
      jmax = getmax(jmax)
      bsmax = getmax(bsmax)

      call plssub(2,2)
      call plinit()

c/    Total Currents :

      call pllsty(1)
      call plenv (rmin, rmax, ymin, cmax, 0, 0)

      call pllab ('Normalized Radius #gr#fn', 'Currents (MA)',
     .  'Currents at different times')

      lj = 0
      do j = 1, nprof
	 jj = nprof - j + 1
	 lj = lj + 1
	 if (lj.gt.8) lj = 8
	 call pllsty(lj)
         call plline(nrho,rho,ctot(1,jj))
      enddo

c/    Total current densities:
c/    (Solid line --> more recent profile)

      call pllsty(1)
      call plenv (rmin, rmax, ymin, jmax, 0, 0)

      call pllab ('Normalized Radius #gr#fn',
     .  'Current Density (MA/m#u2#d)',
     .  'Current Densities at different times')

      lj = 0
      do j = 1, nprof
	 jj = nprof - j + 1
	 lj = lj + 1
	 if (lj.gt.8) lj = 8
	 call pllsty(lj)
         call plline(nrho,rho,jtot(1,jj))
      enddo

c/    Bootstrap Current:

      call pllsty(1)
      call plenv (rmin, rmax, ymin, bsmax, 0, 0)

      call pllab ('Normalized Radius #gr#fn',
     .  'Bootstrap (MA/m#u2#d)',
     .  'Bootstrap Densities at different times')

      lj = 0
      do j = 1, nprof
	 jj = nprof - j + 1
	 lj = lj + 1
	 if (lj.gt.8) lj = 8
	 call pllsty(lj)
         call plline(nrho,rho,jbs(1,jj))
      enddo

c/    Safety Factor:

      call pllsty(1)
      call plenv (rmin, rmax, ymin, qmax, 0, 0)

      call pllab ('Normalized Radius #gr#fn',
     .  'q', 'MHD safety factor')

      lj = 0
      do j = 1, nprof
	 jj = nprof - j + 1
	 lj = lj + 1
	 if (lj.gt.8) lj = 8
	 call pllsty(lj)
         call plline(nrho,rho,qmhd(1,jj))
      enddo

c/    PLOTS FOR ELECTRON AND ION BALANCE PARAMETERS (LOCAL AND
C/    INTEGRATED)

c/    First, determine scales for plotting:

      xmin = 0.0
      xmax = 1.10

      pemax = -1.0e30
      pemin = 1.0e30
      pemaxi = -1.0e30
      pemini = 1.0e30
      pimax = -1.0e30
      pimin = 1.0e30
      pimaxi = -1.0e30
      pimini = 1.0e30

      do i = 1, nr
	pemax = amax1(pemax, pfuse(i), pconde(i), pconve(i), pbe(i),
     .     prfe(i))
	pemin = amin1(pemin, prade(i), pconde(i), pconve(i))
	pemaxi = amax1(pemaxi, pfusei(i), pcondei(i), pconvei(i),
     .     pbei(i), prfei(i))
	pemini = amin1(pemini, pradei(i), pcondei(i), pconvei(i))
	pimax = amax1(pimax, pfusi(i), pcondi(i), pconvi(i), pbi(i),
     .     prfi(i))
	pimin = amin1(pimin, pcondi(i), pconvi(i))
	pimaxi = amax1(pimaxi, pfusii(i), pcondii(i), pconvii(i),
     .     pbii(i), prfii(i))
	pimini = amin1(pimini, pcondii(i), pconvii(i))
      enddo

      pemax = getmax(pemax)
      pemaxi = getmax(pemaxi)
      pimax = getmax(pimax)
      pimaxi = getmax(pimaxi)

      if (pemin.lt.0) then
	 pemin = - getmax(abs(pemin))
      else
	 pemin = 0.0
      endif
      if (pemini.lt.0) then
	 pemini = - getmax(abs(pemini))
      else
	 pemini = 0.0
      endif
      if (pimin.lt.0) then
	 pimin = - getmax(abs(pimin))
      else
	 pimin = 0.0
      endif
      if (pimini.lt.0) then
	 pimini = - getmax(abs(pimini))
      else
	 pimini = 0.0
      endif


c/    START PLOTTING!

c/    Electron Quantities:

      call pllsty(1)
      call plenv (xmin, xmax, pemin, pemax, 0, 0)

      call pllab ('Normalized Radius #gr#fn', 'Powers (eV/s)',
     .  'Local Electron Power Balance')

      ls = 0
      call plline(nr,r,pfuse)
      call pllsty(2)
      call plline(nr,r,pconde)
      call pllsty(3)
      call plline(nr,r,pconve)
      if (pbei(nr).gt.0) then
         call pllsty(4)
	 call plline(nr,r,pbe)
	 ls = 1
      endif
      if (prfei(nr).gt.0) then
         call pllsty(4+ls)
	 call plline(nr,r,prfe)
	 ls = ls + 1
      endif
      call pllsty(4+ls)
      call plline(nr,r,prade)

c/    Integrated Power Balance for electrons:

      call pllsty(1)
      call plenv (xmin, xmax, pemini, pemaxi, 0, 0)

      call pllab ('Normalized Radius #gr#fn', 'Powers (MW)',
     .  'Integrated Electron Power Balance')

      ls = 0
      call plline(nr,r,pfusei)
      call pllsty(2)
      call plline(nr,r,pcondei)
      call pllsty(3)
      call plline(nr,r,pconvei)
      if (pbei(nr).gt.0) then
         call pllsty(4)
	 call plline(nr,r,pbei)
	 ls = 1
      endif
      if (prfei(nr).gt.0) then
         call pllsty(4+ls)
	 call plline(nr,r,prfei)
	 ls = ls + 1
      endif
      call pllsty(4+ls)
      call plline(nr,r,pradei)

c/    Ion Quantities:

      call pllsty(1)
      call plenv (xmin, xmax, pimin, pimax, 0, 0)

      call pllab ('Normalized Radius #gr#fn', 'Powers (eV/s)',
     .  'Local Ion Power Balance')

      ls = 0
      call plline(nr,r,pfusi)
      call pllsty(2)
      call plline(nr,r,pcondi)
      call pllsty(3)
      call plline(nr,r,pconvi)
      if (pbii(nr).gt.0) then
         call pllsty(4)
	 call plline(nr,r,pbi)
	 ls = 1
      endif
      if (prfii(nr).gt.0) then
         call pllsty(4+ls)
	 call plline(nr,r,prfi)
	 ls = ls + 1
      endif

c/    Integrated Power Balance for ions:

      call pllsty(1)
      call plenv (xmin, xmax, pimini, pimaxi, 0, 0)

      call pllab ('Normalized Radius #gr#fn', 'Powers (MW)',
     .  'Integrated Ion Power Balance')

      ls = 0
      call plline(nr,r,pfusii)
      call pllsty(2)
      call plline(nr,r,pcondii)
      call pllsty(3)
      call plline(nr,r,pconvii)
      if (pbii(nr).gt.0) then
         call pllsty(4)
	 call plline(nr,r,pbii)
	 ls = 1
      endif
      if (prfii(nr).gt.0) then
         call pllsty(4+ls)
	 call plline(nr,r,prfii)
	 ls = ls + 1
      endif

c/    PLOT PROFILES FOR T_e, T_i, n_e, n_i, Zeff, and chi's:

c/    Determine plotting scales:

      ymin = 0.0
      denmax = -1.0e30
      tkevmax = -1.0e30
      zefmax = -1.0e30
      chimax = -1.0e30

      do i = 1, nr
	 denmax = amax1(denmax, n_e(i))
	 tkevmax = amax1(tkevmax, tekev(i), tikev(i))
         zefmax = amax1(zefmax, zeff(i))
	 chimax = amax1(chimax, chi_e(i), chi_i(i), dfc(i))
      enddo

      denmax = 1.0e20 * getmax(1.0e-20 * denmax)
      tkevmax = getmax(tkevmax)
      zefmax = getmax(zefmax)
      chimax = getmax(chimax)

c/    Write debugging information:
      write (ndbg, dbg1)

c/    Temperatures:

      call pllsty(1)
      call plenv (xmin, xmax, ymin, tkevmax, 0, 0)

      call pllab ('Normalized Radius #gr#fn', 'T#de#u, T#di#u (kev)',
     .  'Temperature Profiles')

      call plline(nr,r,tekev)
      call pllsty(2)
      call plline(nr,r,tikev)

c/    Densities:

      call pllsty(1)
      call plenv(xmin, xmax, ymin, denmax, 0, 0)
      call pllab ('Normalized Radius #gr#fn',
     .    'n#de#u, n#di#u (m#u-3#d)', 'Density Profiles')

      call plline(nr,r,n_e)
      call pllsty(2)
      call plline(nr,r,n_i)

c/    Zeff profile:

      call pllsty(1)
      call plenv(xmin, xmax, ymin, zefmax, 0, 0)
      call pllab ('Normalized Radius #gr#fn',
     .    'Z#deff#u', 'Plasma Z#deff#u Profile')

      call plline(nr,r,zeff)

c/    Chi profiles:

      call plenv(xmin, xmax, ymin, chimax, 0, 0)
      call pllab ('Normalized Radius #gr#fn',
     .    '#gx#fn#de#u, #gx#fn#di#u, D (m#u2#d/s)',
     .    'Transport Coefficients')

      call pllsty(1)
      call plline(nr,r,chi_e)
      call pllsty(2)
      call plline(nr,r,chi_i)
      call pllsty(3)
      call plline(nr,r,dfc)

c/    Plot of MHD metric quantities:

c/    Convert units to MKS and redefine shift first:
      do i = 1, nrmhd
	 vpr(i) = 1.0e-4 * vpr(i)
	 vol(i) = 1.0e-6 * vol(i)
	 shift(i) = 1.0e-2 * (shift(nrmhd) - shift(i))
      enddo

c/    Plotting scales...


      elmax = -1.e30
      delmax = -1.0e30
      vpmax  = -1.0e30
      gmax = -1.0e30
      smax = -1.0e30
      fmax = -1.0e30

      do i = 1, nrmhd
	elmax = amax1(elmax, elong(i))
	delmax = amax1(delmax, triang(i))
	vpmax = amax1(vpmax, vpr(i))
	gmax = amax1(gmax, grho(i), gth(i), gph(i))
	smax = amax1(smax, shift(i))
	fmax = amax1(fmax, F(i))
      enddo

      elmax = getmax(elmax)
      vpmax = getmax(vpmax)
      delmax = getmax(delmax)
      gmax = getmax(gmax)
      smax = getmax(smax)
      fmax = getmax(fmax)

c/    Plot elongation and triangularity:

      call pladv(0)
      call plvsta
      call plwind(0.0, 1.0, 0.0, elmax)
      call pllsty(1)
      call plbox('bcnst',0.0, 0, 'bnstv', 0.0, 0)
      call plline(nrmhd,rmhd,elong)

      call plwind(0.0, 1.0, 0.0, delmax)
      call plbox(' ',0.0, 0, 'cmstv', 0.0, 0)
      call pllsty(2)
      call plline(nrmhd,rmhd,triang)

      call pllsty(1)
      call plmtex('b', 3.5, 0.5, 0.5, '#gr#fn')
      call plmtex('l', 4.5, 0.5, 0.5, 'Elongation #gk#fn')
      call plmtex('r', 4.0, 0.5, 0.5, 'Triangularity #gd#fn')

c/    Plot V' and F:

      call pladv(0)
      call plvsta
      call plwind(0.0, 1.0, 0.0, vpmax)
      call plbox('bcnst',0.0, 0, 'bnstv', 0.0, 0)
      call plline(nrmhd,rmhd,vpr)

      call plwind(0.0, 1.0, 0.0, fmax)
      call plbox(' ',0.0, 0, 'cmstv', 0.0, 0)
      call pllsty(2)
      call plline(nrmhd,rmhd,F)

      call pllsty(1)
      call plmtex('b', 3.5, 0.5, 0.5, '#gr#fn')
      call plmtex('l', 5.5, 0.5, 0.5, 'V#dpr#u (m#u2#d)')
      call plmtex('r', 4.0, 0.5, 0.5, 'F (T-m)')

c/    Plot G's:

      call pllsty(1)
      call plenv (0.0, 1.0, 0.0, gmax, 0, 0)

      call pllab ('#gr#fn', 'G#d#gy#u#fn, G#d#gf#u#fn, G#d#gr#u#fn ',
     .  'Metric Tensor Coefficients')

      call plline(nrmhd,rmhd,gth)
      call pllsty(2)
      call plline(nrmhd,rmhd,gph)
      call pllsty(3)
      call plline(nrmhd,rmhd,grho)

c/    Plot Shafranov shift:

      call pllsty(1)
      call plenv (0.0, 1.0, 0.0, smax, 0, 0)

      call pllab ('#gr#fn', '#gD#fn (m)', 'Shift')

      call plline(nrmhd,rmhd,shift)

      call plend

  500 format (a25)
  750 format (1x, i2.2, '/', i2.2, '/', i4, ', ',
     .        i2.2, ':', i2.2, ':', i2.2)
 1000 format (1x, 'Number of time data points = ', i3, ' while mxp = ',
     .  i3, /, 1x, 'Increase value of mxp and recompile')
 2000 format (1x, 'Number of radial data points = ', i3,
     .' while mxr = ', i3, /, 1x, 'Increase value of mxr and recompile')
 3000 format (1x, 'Number of radial data points = ', i3,
     .' while mxr = ', i3, /, 1x, 'Increase value of mxr and recompile')

      end

c/**********************************************************************

      real function getmax(p)

c/    This function returns an appropriate scale for the quantity p

      real p, s, power, pmod, pincr
      integer is, indx

      s = alog10 (p)
      is = aint (s)
      power = 10 ** is
      pmod = amod(p,power)
      indx = aint (p / power)
      if (pmod.lt.(0.5*power)) then
	 pincr = 0.5 * power
      else
         pincr = power
      endif
      getmax = power * indx + pincr

      return
      end
