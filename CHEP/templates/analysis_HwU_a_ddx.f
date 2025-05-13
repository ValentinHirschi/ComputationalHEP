c
c Example analysis for "p p > e+ ve [QCD]" process.
c Example analysis for "p p > e- ve~ [QCD]" process.
c Example analysis for "p p > mu+ vm [QCD]" process.
c Example analysis for "p p > mu- vm~ [QCD]" process.
c Example analysis for "p p > ta+ vt [QCD]" process.
c Example analysis for "p p > ta- vt~ [QCD]" process.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nwgt
      character*(*) weights_info(*)
      integer i,kk,l
      character*6 cc(2)
      data cc/'|T@NLO','|T@LO '/
      call HwU_inithist(nwgt,weights_info)
      do i=1,2
        l=(i-1)*8
        call HwU_book(l+1,'total rate '//cc(i), 5,0.5d0,5.5d0)
        call HwU_book(l+2,'cos(ThetaQG) '//cc(i), 100,-1d0,1d0)
        call HwU_book(l+3,'cos(ThetaQXG) '//cc(i), 100,-1d0,1d0)
        call HwU_book(l+4,'Eg '//cc(i), 100,0d0,100d0)
      enddo
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(dummy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision dummy
      call HwU_write_file
      return                
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'nexternal.inc'
      integer istatus(nexternal)
      integer iPDG(nexternal)
      double precision p(0:4,nexternal)
      double precision wgts(*)
      integer ibody
      double precision wgt,var
      integer i,kk,l
      double precision CosThetaQG,CosThetaQXG,EG
      double precision cosThetaAnalysis
      external cosThetaAnalysis
      if (nexternal.ne.4) then
         write (*,*) 'error #1 in analysis_fill: '/
     &        /'only for process "z > d dx QCD=0 [QCD]"'
         stop 1
      endif
      if (p(0,4).eq.0.0d0) then
          CosThetaQG = 0.0d0
          CosThetaQXG = 0.0d0
          EG = 0.0d0
      else
          CosThetaQG = cosThetaAnalysis(p(0,4),p(0,2))
          CosThetaQXG = cosThetaAnalysis(p(0,4),p(0,2))
          EG = p(0,4)
      endif
      var    = 1.d0
      do i=1,2
         l=(i-1)*8
         if (ibody.ne.3 .and.i.eq.2) cycle
         call HwU_fill(l+1,var,wgts)
         call HwU_fill(l+2,CosThetaQG,wgts)
         call HwU_fill(l+3,CosThetaQXG,wgts)
         call HwU_fill(l+4,EG,wgts)
      enddo
 999  return      
      end
      
      function cosThetaAnalysis(p1, p2)
          implicit none
          double precision cosThetaAnalysis, p1(0:3), p2(0:3), p1_norm, p2_norm
          p1_norm = dsqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
          p2_norm = dsqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))
          cosThetaAnalysis = (p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))/(p1_norm*p2_norm)
          return
      end
