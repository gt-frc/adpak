      program adpak

c/    Main Driver for the ADPAK code based on the original ADLHAX routine
c/    by R. Hulse. This new version introduces Namelists, defines I/O
c/    units, has more comments, and also creates the blockdata file
c/    needed by the impurity transport routine.
c/    01/17/95, John Mandrekas, GIT
c/    08/08/97, jm, included output file for fitting programs:
c/    08/26/98, jm, replaced zne with ane in table output.

c/    INPUT VARIABLES
c/    ---------------
c/    imode      : Determines what kind of run we have:
c/                 0 --> calculate atomic coefficients for zte, zne
c/                 1 --> create blockdata file for transport code
c/    inucz      : Atomic number, Z, of desired element
c/    zte        : Electron temperature in keV
c/    zne        : Electron density in cm-3
c/    laden      : Flag for excitation energy and ion. potential calc.
c/                 If 0 --> use Mayer formalism
c/                 If 1 --> use More formalism
c/    ladtip     : If 1 --> use tabulated ionization potentials
c/    leci       : Flag for calculation of ionization rates:
c/                 1 --> XSNQ
c/                 2 --> Belfast group (for H through O)
c/                 3 --> Younger (Scandium and Fe)
c/    ldrmlt     : Dielectronic multiplier
c/                 0 --> Use CDNN, CDNM arrays as given
c/                 1 --> Set CDNN, CDNM equal to 1
c/                 2 --> Use Y. Hahn factors (input) to set up CDNN, CDNM
c/    yhnna  |
c/    yhnnb  |
c/    yhnnc   \
c/             } : Hahn input coefficients (if ldrmlt = 2)
c/    yhnma   /
c/    yhnmb  |
c/    yhnmc  |
c/    ncxb       : Number of NB components (if 0 --> No NB's_)
c/    ncxopt     : Selects cross sections to be used
c/                 1 --> OSAS
c/                 2 --> GJ
c/                 3 --> OSCT
c/    ivunit     : Units for NB energy
c/                 1 --> cm/s
c/                 2 --> keV / amu
c/    anneut     : Array of neutral densities (cm-3)
c/    vneut      : Array of neutral energies / velocities
c/
c/    The following input variables are needed to create the blockdata
c/    file outblk.dat needed in the impurity transport calculations:

c/    nte        : Number of electron temperatures in table
c/    nne        : Number of electron densities in table
c/    tei        : Array of nte Te values in keV
c/    anei       : Array of nne ne values in /cm^3
c/    nmin       : Exponent of minimum temperature for tables (10^nmin)
c/    nmax       : Exponent of maximum temperature for tables (10^nmax)
c/
c/    CALCULATED OUTPUT VARIABLES (selection):
c/    ----------------------------------------
c/    nvalnc(jq)   : Principal quantum number of highest occupied shell
c/    apn(jq,jn)   : Number of electrons in shell JN for species with
c/                   charge JQ-1
c/    ein(jq,jn)   : Ionization potential in keV for electron in shell JN
c/    rclion(jq)   : collisional ionization rates for charge state JQ (sec-1)
c/    rrarec(jq)   : radiative recombination rates for state jq (sec-1)
c/    rdirec(jq)   : dielectronic recomb. rate for state jq (sec-1)
c/    rcxrec(jq)   : charge exchange recomb. rate for state jq (sec-1)
c/    radrrc(jq)   : radiation rate for radiative recomb. process (W)
c/    raddrc(jq)   :          "        dielectronic recomb. process (W)
c/    radclx(jq)   :          "        collisional excitation process (W)
c/    radbrm(jq)   :          "        bremsstrahlung process (W)
c/    radbcx(jq)   :          "        beam CX process        (W)
c/    cerad        : Total radiation rate from all states (W)
c/    ceavgz       : <Z>
c/
c/    Note that ADPAK outputs all rates per ion (i.e., they are
c/    multiplied by the electron density, n_e, already). The older
c/    version of the code from INEL used the same variable names for
c/    the "pure" rates, i.e the <sigmaV>. To go from the new rates
c/    to the old, we must divide by n_e, which is done before writing
c/    things to the output.
C
C     REV: 3/30/83 TO ADD BELFAST IONIZATION RATES
C     REV: 10/10/82 CREATED FROM OLD BEQLST FOR AD PACKAGE
C
      common /adsdat/ nucz, nspc, apn(100,10), aqn(100,10),
     .    nvalnc(100), sigma(10,10), qn(100,10), en(100,10),
     .    ein(100,5), enm(100,5,10), fnm(100,5,10), fnn(100), enn(100)
C
      common /adelec/ rclion(100), rrarec(100), rdirec(100),
     .    cizlos(100), radrrc(100), raddrc(100), radclx(100),
     .    radbrm(100)
C
      common /adneut/ rcxrec(100), radbcx(100)
C
      common /adparm/ laden, ladtip, ltipok, leci, leciok
C
      common /addiel/ cdnn(100), cdnm(100), ldrmlt, yhnna, yhnnb, yhnnc,
     .    yhnma, yhnmb, yhnmc
C
      common /adceq/ cefrac(100), ceavgz, cerad
C
      dimension ihmomy(2), ihtip(2), ihvu(2), ihxs(5), iheci(3)
C
      dimension anneut(5), vneut(5)

      real tei(100), anei(100), te_tabl(100), radre(100), radrt(100),
     .     zavg(100)
      real telog, lzlog
      CHARACTER(len=100) arg

      double precision ainzrjq, ainzejq, radrjq, recrjq
      integer nin, nout, nplt1, nplt2, iout

c/    Include table files:
      include 'params.inc'
      include 'table.inc'
      include 'dzfrcn.inc'

      data nin /10/, iout /12/, nout/14/, nplt1/16/, nplt2/17/
      data cvtrm1 /1.e-6/
      data cvtrm2 /1.e+6/

      namelist /inp/ inucz, zte, zne, laden, ladtip, leci, ldrmlt,
     .    yhnna, yhnnb, yhnnc, yhnma, yhnmb, yhnmc, ncxb, ncxopt,
     .    ivunit, anneut, vneut, nte, nne, tei, anei, imode, nmin,
     .    nmax
C
c/    Do file I/O:
      CALL getarg(1, arg)
      open (nin,file =TRIM(arg),status ='old')
      open (iout,file='outadpk.txt',status='unknown')
      open (nout,file='outblk.dat ',status='unknown')
      open (nplt1,file='outplt.txt ',status='unknown')
      open (nplt2,file='adfits.txt',status='unknown')

c/    Define some parameters used in the output:

      ihmomy(1) = 4hMAYR
      ihmomy(2) = 4hMORE
      ihtip(1) = 3hOUT
      ihtip(2) = 3h
      ihvu(1) = 4hCM/S
      ihvu(2) = 4hKV/A
      ihxs(1) = 4hOSAS
      ihxs(2) = 4hGJ
      ihxs(3) = 4hOSCT
      iheci(1) = 4hXSNQ
      iheci(2) = 4hBFST
      iheci(3) = 3hYGR
C
c/    Initialize...

      inucz = 26
      zte = 2.5
      zne = 2.e13
      anneut(1) = 2.0e08
      vneut(1) = 2.0e08
      ivunit = 1
      ncxb = 0
      ncxopt = 3
C
      ldrmlt = 1
      yhnna = 0.2
      yhnnb = 0.7
      yhnnc = 1.4
      yhnma = 0.1
      yhnmb = 0.55
      yhnmc = 1.0
C
      laden = 0
      ladtip = 0
C
      leci = 1
C
c/  READ INPUT PARAMETERS:

      read (nin,inp)

      write (*,*) "inucz = ",inucz
c/  Stop execution if INCUCZ = 0 :
      if (inucz.eq.0) stop

C*******************************************************************

c/    Start calculations:

c/    Setup basic atomic structure data for desired element:

      call adset(inucz)
      call adbcxr(anneut,vneut,ncxb,ncxopt,ncxerr,ivunit)
      aneplt = zne

      if (imode.eq.1) go to 1111

      call aderc(zte,zne)

      if (ncxerr.ne.0) then
         write (iout,8000) ncxerr
         stop
      endif
C
c/    Perform coronal equilibrium calculations:

      call adce
C
C*******************************************************************
C
C     PRINTOUT RATES TO DESIRED UNIT
C
      write (iout,8100) nucz, zte, zne, ncxopt, ihxs(ncxopt)
      if (ncxb.le.0) write (iout,8200)
      if (ncxb.gt.0) write (iout,8300) (anneut(j),j = 1,ncxb)
      if (ncxb.gt.0) write (iout,8400) ihvu(ivunit), (vneut(j),j = 1,
     .    ncxb)
C
      write (iout,8500) ihmomy(laden+1), ihtip(ltipok+1)
C
      write (iout,8600) iheci(leciok)
C
      if (ldrmlt.eq.1) write (iout,8700)
      if (ldrmlt.eq.2) write (iout,8800) yhnna, yhnnb, yhnnc, yhnma,
     .    yhnmb, yhnmc
C
      write (iout,8900)
C
C
      do 10 jq = 1, nspc
C
         ivalnc = nvalnc(jq)
         zclion = rclion(jq) / zne
         zrarec = rrarec(jq) / zne
         zdirec = rdirec(jq) / zne
         zcxrec = rcxrec(jq) / zne
C
         totrec = rrarec(jq) + rdirec(jq) + rcxrec(jq)
         zrtot = totrec / zne
C
         write (iout,9000) jq, ivalnc, apn(jq,ivalnc), ein(jq,ivalnc),
     .       rclion(jq), rrarec(jq), rdirec(jq), rcxrec(jq), totrec
C
         write (iout,9100) zclion, zrarec, zdirec, zcxrec, zrtot
C
   10 continue
C
      write (iout,9200)
C
C
C     PRINTOUT NORMALIZED SPECIES FRACTIONS AND RADIATION RATES
C
      write (iout,9300)
C
C
      do 20 jq = 1, nspc
C
         zrdrrc = radrrc(jq) / zne
         zrddrc = raddrc(jq) / zne
         zrdclx = radclx(jq) / zne
         zrdbrm = radbrm(jq) / zne
         zrdbcx = radbcx(jq) / zne
C
         totrad = radrrc(jq) + raddrc(jq) + radclx(jq) + radbrm(jq) +
     .       radbcx(jq)
         ztrad = totrad / zne
C
         write (iout,9400) jq, cefrac(jq), radrrc(jq), raddrc(jq),
     .       radclx(jq), radbrm(jq), radbcx(jq), totrad
C
         write (iout,9500) zrdrrc, zrddrc, zrdclx, zrdbrm, zrdbcx, ztrad
C
   20 continue
C
C
      zcerad = cerad / zne
C
      write (iout,9600) ceavgz, zcerad, cerad
      stop

c/    Control is transfered to the statement below if imode = 1:

 1111 continue

c/    This part of the code creates the blockdata file for the transport
c/    code:

      do ite = 1, nte
        do ine = 1, nne

          zte = tei(ite)
	  zne = anei(ine)

          altei(ite) = alogxx(zte)
	  alnei(ine) = alogxx(zne*cvtrm2)

          call aderc(zte,zne)

          do jq = 1, nspc

            alinzr(jq,ite) = alogxx(rclion(jq)*cvtrm1/zne)
            alrecr(jq,ite,ine) =
     .            alogxx((rrarec(jq)+rdirec(jq)+rcxrec(jq))/zne*cvtrm1)
            alradr(jq,ite,ine) = alogxx((radrrc(jq)+raddrc(jq)+
     .             radclx(jq)+radbrm(jq)+radbcx(jq))/zne*cvtrm1)
	    write (*,*) "zne = ",zne,"alradr = ",alradr(jq,ite,ine)
          enddo

        enddo

      enddo

c/    WRITE THE BLOCKDATA FILE:

         write (nout,5700)
         write (nout,6400) nte, nne, nspc
         write (nout,5800)
         write (nout,5900) nte
	 write (*,*) (altei(i),i = 1,nte)
         write (nout,6500) (altei(i),i = 1,nte)
         nleft = nte - 5 * (nte/5)
         if (nleft.ge.0) then
            backspace nout
            if (nleft.eq.0) then
               nleft = 5
               backspace nout
            endif
            go to (50,60,70,80,90), nleft
   50       write (nout,6600) (altei(ite),ite = nte-nleft+1,nte)
            go to 100
   60       write (nout,6700) (altei(ite),ite = nte-nleft+1,nte)
            go to 100
   70       write (nout,6800) (altei(ite),ite = nte-nleft+1,nte)
            go to 100
   80       write (nout,6900) (altei(ite),ite = nte-nleft+1,nte)
            go to 100
   90       write (nout,7000) (altei(ite),ite = nte-nleft+1,nte)
         endif
  100    continue
         write (nout,5800)
         write (nout,6000) nne
         write (nout,6500) (alnei(i),i = 1,nne)
         nleft = nne - 5 * (nne/5)
         if (nleft.ge.0) then
            backspace nout
            if (nleft.eq.0) then
               nleft = 5
               backspace nout
            endif
            go to (110,120,130,140,150), nleft
  110       write (nout,6600) (alnei(ine),ine = nne-nleft+1,nne)
            go to 160
  120       write (nout,6700) (alnei(ine),ine = nne-nleft+1,nne)
            go to 160
  130       write (nout,6800) (alnei(ine),ine = nne-nleft+1,nne)
            go to 160
  140       write (nout,6900) (alnei(ine),ine = nne-nleft+1,nne)
            go to 160
  150       write (nout,7000) (alnei(ine),ine = nne-nleft+1,nne)
         endif
  160    continue
         write (nout,5800)
c
         do 230 jq = 1, nspc
            write (nout,6100) jq, nte
            write (nout,6500) (alinzr(jq,ite),ite = 1,nte)
            nleft = nte - 5 * (nte/5)
            if (nleft.ge.0) then
               backspace nout
               if (nleft.eq.0) then
                  nleft = 5
                  backspace nout
               endif
               go to (170,180,190,200,210), nleft
  170          write (nout,6600) (alinzr(jq,ite),ite = nte-nleft+1,nte)
               go to 220
  180          write (nout,6700) (alinzr(jq,ite),ite = nte-nleft+1,nte)
               go to 220
  190          write (nout,6800) (alinzr(jq,ite),ite = nte-nleft+1,nte)
               go to 220
  200          write (nout,6900) (alinzr(jq,ite),ite = nte-nleft+1,nte)
               go to 220
  210          write (nout,7000) (alinzr(jq,ite),ite = nte-nleft+1,nte)
            endif
  220       continue
            write (nout,5800)
  230    continue
         do 300 jq = 1, nspc
            write (nout,6200) jq, nte, nne
	    write (nout,6500) ((alradr(jq,ite,ine),ite = 1,nte),
     .  	 ine = 1, nne)
            nleft = nte * nne - 5 * (nte*nne/5)
            if (nleft.ge.0) then
               backspace nout
               if (nleft.eq.0) then
                  nleft = 5
                  backspace nout
               endif
               go to (240,250,260,270,280), nleft
  240          write (nout,6600) (alradr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 290
  250          write (nout,6700) (alradr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 290
  260          write (nout,6800) (alradr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 290
  270          write (nout,6900) (alradr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 290
  280          write (nout,7000) (alradr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
            endif
  290       continue
            write (nout,5800)
  300    continue
         do 370 jq = 1, nspc
            write (nout,6300) jq, nte, nne
            write (nout,6500) ((alrecr(jq,ite,ine),ite = 1,nte),ine = 1,
     .          nne)
            nleft = nte * nne - 5 * (nte*nne/5)
            if (nleft.ge.0) then
               backspace nout
               if (nleft.eq.0) then
                  nleft = 5
                  backspace nout
               endif
               go to (310,320,330,340,350), nleft
  310          write (nout,6600) (alrecr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 360
  320          write (nout,6700) (alrecr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 360
  330          write (nout,6800) (alrecr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 360
  340          write (nout,6900) (alrecr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
               go to 360
  350          write (nout,7000) (alrecr(jq,ite,nne),ite = nte-nleft+1,
     .             nte)
            endif
  360       continue
            write (nout,5800)
  370    continue

c/    CREATE PLOT FILE FOR COMPARISON:

c/    Create data for comparison plotting between calculated and
c/    interpolated values of the cooling rate:
c/    In order to be able to compare with the graphs in Post's paper,
c/    the different quantities are converted into cgs units, and the
c/    temperatures are in keV.

      ane = aneplt
      indx = 0

      nPeriods = nmax - nmin
      do j = 1, nPeriods
         zte = 0.0
	 nexp = nmin + j - 1
	 if (nexp.ge.0) delte = 10 ** nexp
	 if (nexp.lt.0) then
	    nexp = -nexp
	    delte = 1.0 / (10 ** nexp)
	 endif
	 imax = 9
	 if (j.eq.nPeriods) imax = 10
         do  i = 1, imax
           zte = zte + delte
           indx = indx + 1
           te_tabl(indx) = zte
           call aderc(zte,ane)
           call adce
           radre(indx) = 1.e7 * cerad / ane

c * * compare with tables
c
           anemks = ane * cvtrm2

           do  jq = 1, nspc
             call lookup(jq,zte,anemks,ainzrjq,ainzejq,radrjq,recrjq)
             tradr(jq) =  radrjq
	     ainz(jq) = anemks * ainzrjq
	     rec(jq) = anemks * recrjq
           enddo

c/         calculate coronal equilibrium fractions:

	   call cesolve(nspc)

	   zsum = 0.0
           sumrt = 0.0
           do  jq = 1, nspc
              sumrt = sumrt + fracz(jq) * tradr(jq)
	      zsum = zsum + (jq-1) * fracz(jq)
           enddo

           radrt(indx) = 1.e13 * sumrt
	   zavg(indx) = zsum

         enddo
      enddo


c/    Write results to plot files:
      write (nplt1,7100) inucz
      write (nplt2, 4000)
      do i = 1, indx
	 te_ev = 1000.0 * te_tabl(i)
	 telog = alog(te_ev)
	 lzlog = alog(radrt(i))
         write (nplt1,7200) te_tabl(i), radre(i), radrt(i)
	 if (1.0e-3*te_ev.GE.tei(1))
     .      write (nplt2,4500) telog, lzlog, zavg(i)
      enddo

 5700 format (6X,'data nte(4),nne(4),nchrgsr(4) /')
 5800 format (5X,'.     /')
 5900 format (6X,'data (altei(4,i),i=1,',I5,') /')
 6000 format (6X,'data (alnei(4,i),i=1,',I5,') /')
 6100 format (6X,'data (alinzr(4,',I3,',ite),ite=1,',I5,') /')
 6200 format (6X,'data ((alradr(4,',I3,',ite,ine),ite=1,',I5,'),ine=1,',
     .    I5,') /')
 6300 format (6X,'data ((alrecr(4,',I3,',ite,ine),ite=1,',I5,'),ine=1,',
     .    I5,') /')
 6400 format (5X,'.',I5,',',1X,I5,',',1X,I5)
 6500 format (5X,'.',5(1X,1PE11.4,','),200(/5X,'.',5(1X,1PE11.4,',')))
 6600 format (5X,'.',1X,1PE11.4)
 6700 format (5X,'.',1X,1PE11.4,',',1X,1PE11.4)
 6800 format (5X,'.',1X,1PE11.4,',',1X,1PE11.4,',',1X,1PE11.4)
 6900 format (5X,'.',1X,1PE11.4,',',1X,1PE11.4,',',1X,1PE11.4,',',1X,1
     .    PE11.4)
 7000 format (5X,'.',1X,1PE11.4,',',1X,1PE11.4,',',1X,1PE11.4,',',1X,1
     .    PE11.4,',',1X,1PE11.4)
 7100 format (5X,'Te_keV',6X,'Theory',8X,'Table', 1x, i4)
 7200 format (3(1X,E12.5))

 8000 format (1X,'NCXERR =',I3)
 8100 format (1X,10('-------'),'-'//12X,'Z =',I3,'  TE =',F7.3,'  NE =',
     .    1PE9.2,'  NCXOPT =',I2,' (',A4,')'/)
 8200 format (25X,'***** NO BEAMS *****')
 8300 format (1X,'BEAM DENSITY(CM-3):  ',1P5E11.2)
 8400 format (1X,'BEAM VELOCITY(',A4,'):',1P5E11.2)
 8500 format (/1X,5X,A4,' ENERGY LEVELS WITH',A3,
     .    ' TABULATED IONIZATION POTENTIALS ADDED')
 8600 format (/1X,'IONIZATION RATES: ',A4)
 8700 format (/1X,25X,'DIELECTRONIC FACTORS = 1.0')
 8800 format (/1X,'HAHN DIELECTRONIC FACTORS: NN =',3F7.3,', NM =',3F7.3
     .    )
 8900 format (/,1X,'SPC',2X,'VS',2X,'VSP',3X,'ION POT',4X,'RCLION',5X,
     .    'RRAREC',5X,'RDIREC',5X,'RCXREC',5X,'TOTREC'/)
 9000 format (1X,I2,I5,F5.0,F9.3,1X,1P5E11.2)
 9100 format (1X,23X,5(' (',1PE8.2')'))
 9200 format (/)
 9300 format (/1X,'SPC',2X,'CEQ FRAC',2X,'RADRRC',5X,'RADDRC',5X,
     .    'RADCLX',5X,'RADBRM',5X,'RADBCX',6X,'TOTAL'/)
 9400 format (1X,I2,1PE10.2,1PE10.2,1P5E11.2)
 9500 format (1X,12X,6(' (',1PE8.2,')'))
 9600 format (/6X,'<Z> =',F6.2,6X,'RAD COEFF =',1PE9.2,' W-CM3',6X,
     .    'RAD =',1PE9.2,' W/ION'///)
 4000 format (1x, 'log_TeV', 5x, 'log_Lz', 5x, 'Zavg')
 4500 format (1x, 3(e12.5,1x))
      stop
      end
