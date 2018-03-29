      subroutine lookup(jq,te,ane,ainzrjq,ainzejq,radrjq,recrjq)

      include 'params.inc'
      include 'table.inc'
      double precision ainzrjq, ainzejq, radrjq, recrjq

      alte = alogxx(te)
      if (alte.lt.altei(1)) alte = altei(1)
      if (alte.gt.altei(nte)) alte = altei(nte)
      alne = alogxx(ane)
      if (alne.lt.alnei(1)) alne = alnei(1)
      if (alne.gt.alnei(nne)) alne = alnei(nne)

      do 10 i = 1, nte
         itep = i
         item = i - 1
         if (altei(i).gt.alte) go to 20
   10 continue
c
   20 f = (alte-altei(item)) / (altei(itep)-altei(item))
c
      do 30 i = 1, nne
         inep = i
         inem = i - 1
         if (alnei(i).gt.alne) go to 40
   30 continue
c
   40 g = (alne-alnei(inem)) / (alnei(inep)-alnei(inem))
c
      ainzrjq = alinzr(jq,item) + f * (alinzr(jq,itep)-alinzr(jq,item))
      ainzejq = alinze(jq,item) + f * (alinze(jq,itep)-alinze(jq,item))

      radrjq = (1-f) * (1-g) * alradr(jq,item,inem) + (1-f) * g * 
     .    alradr(jq,item,inep) + f * (1-g) * alradr(jq,itep,inem) + f *
     .    g * alradr(jq,itep,inep)
      recrjq = (1-f) * (1-g) * alrecr(jq,item,inem) + (1-f) * g * 
     .    alrecr(jq,item,inep) + f * (1-g) * alrecr(jq,itep,inem) + f *
     .    g * alrecr(jq,itep,inep)

      recrjq = 1.0d+01 ** recrjq
      radrjq = 1.0d+01 ** radrjq
      ainzejq = 1.0d+01 ** ainzejq
      ainzrjq = 1.0d+01 ** ainzrjq
      return
      end
      
     
      function alogxx(a)

      if (a.lt.1.0e-37) a = 1.0e-37

      alogxx = alog10(a)
 
      return
      end
      

      subroutine cesolve(n)

c/    calculates steady state (CE) ionic species fractions fracz(j)

      include 'params.inc'
      include 'dzfrcn.inc'
      double precision sum, prod, ratio

      sum = 0.0d+00
      nspc = n
      nspc1 = nspc - 1
      prod = 1.0d+00
      do i = 1, nspc1
	 j = nspc - i + 1
	 ratio = ainz(j-1) / rec(j)
	 sum = 1.0 + ratio * prod
	 prod = sum
      enddo

      fracz(1) = 1.0 / sum
      do i = 1, nspc1
	 fracz(i+1) = fracz(i) * ainz(i) / rec(i+1)
      enddo

      return
      end
