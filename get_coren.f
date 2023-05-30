      program main
      integer  nsq(5,5),  ix,iy, iev,icnf,ishx,ishy,istep, ncnf,
     *             nsqall(5,5,100),  i,j,k,   l,m,mm, kev(100)
      real*8 xx(5), yy(5), x0,y0,  xa,ya, delr, d1,d2
      common/expdat/ xx,yy, x0,y0,  xa,ya, delr, nsq
      
      real xc,yc,h,   xcmi,xcma,ycmi,ycma, xaxis(100),yaxis(100),
     *      stot100,stot300,atot100,atot300
      real err,fl,  rmax,  chtot(10),chtotsq(10),est
      data chtot/494417.0,559874.0,543855.0,514927.0,463112.0,
     *            475309.0,496181.0,535286.0,691261.0,515000.0/,
     *    chtotsq/407639.0,467228.0,446556.0,424948.0,382137.0,
     *            426369.0,405786.0,447263.0,593508.0,422671.0/
      double precision    fff,  xx0(5),yy0(5),fmin,ff(5),x00,y00
      data xx0/5.0d0,15.0d0,-15.0d0,15.0d0,-15.0d0/
      data yy0/5.0d0,-15.0d0,15.0d0,15.0d0,-15.0d0/
*      data xx0/5.0d0,25.0d0,-25.0d0,25.0d0,-25.0d0/
*      data yy0/5.0d0,-25.0d0,25.0d0,25.0d0,-25.0d0/
      common/func/ fff
      real*8  a0,a1,a2,a3,eest
      real*8  co0,c10,c20,c30,sco0,sc10,sc20,sc30
      common/energy/ a0,a1,a2,a3,eest
      real     ee(5)

      external ch
      external ldffit
      external futil


*        write(*,*)  ''
*      OPEN(UNIT=10,FILE='get_coren.in',STATUS='OLD',
*     *                FORM='FORMATTED')
*      read(10,*) ncnf
*      write(*,*) '   >> ncnf=',ncnf
*      read(10,*) istep
*      write(*,*) '   >> istep=',istep
*      close(10)
        write(*,*)  ''
      open(unit=10,file='get_coren.in',status='OLD',form='FORMATTED')
      read(10,*) co0
      write(*,*) ' >> co0=',co0
      read(10,*) c10
      write(*,*) ' >> c10=',c10
      read(10,*) c20
      write(*,*) ' >> c20=',c20
      read(10,*) c30
      write(*,*) ' >> c30=',c30
      read(10,*) sco0
      write(*,*) ' >> sco0=',sco0
      read(10,*) sc10
      write(*,*) ' >> sc10=',sc10
      read(10,*) sc20
      write(*,*) ' >> sc20=',sc20
      read(10,*) sc30
      write(*,*) ' >> sc30=',sc30
      close(10)
        write(*,*)  ''


      open(unit=20,file='cnfgs',status='OLD',form='FORMATTED')
      read(20,*) ncnf,istep
      write(*,*) '    ncnf=',ncnf,'  istep=',istep
      h = 2.0*istep

      do k=1,ncnf
*        read(20,301)  icnf,iev,ic,jc,xaxis(k),yaxis(k)
        read(20,301)  icnf,kev(k),ic,jc,xaxis(k),yaxis(k)
	read(20,400)  ((nsqall(i,j,k),i=1,5),j=1,5)
      enddo
      close(20)

      do i=1,5
        xx(i) = h*(i - 3)
        yy(i) = h*(i - 3)
      enddo
      write(*,*)  ' xx:' 
      write(*,*) (xx(i),i=1,5)
      write(*,*)  ' yy:' 
      write(*,*) (yy(i),i=1,5)
      write(*,*)  ''
      write(*,*)  ''

      write(*,*)  ''
      open(unit=30,file='delr',status='NEW',form='FORMATTED')
      open(unit=31,file='deln',status='NEW',form='FORMATTED')

      atot100 = 0.0
      atot300 = 0.0
      stot100 = 0.0
      stot300 = 0.0
      d1 = 0.0d0
      d2 = 0.0d0
      do k=1,ncnf     !  loop in configurations begins
*      do k=1,1
        do i=1,5
	  do j=1,5
	    nsq(i,j) = nsqall(i,j,k)
	  enddo
	enddo
	xa = xaxis(k)
	ya = yaxis(k)
        write(*,*)  ''
        write(*,*)  ''
        write(*,*)  '         ***   xa=',xa,' ya=',ya

      l = 0
  77  l = l + 1
      if(l.ge.6)  then
        fmin = 1.0e5
        do m=1,5
	  if(fmin.gt.ff(m))  then
	    fmin = ff(m)
	    mm = m
            write(*,*)  '           ******  mm=',mm,'  fmin=',fmin
	  endif
	enddo
	x00 = xx0(mm)
        y00 = yy0(mm)
      else
        x00 = xx0(l)
        y00 = yy0(l)
      endif

 88   call mninit(5,6,7)

      call mnparm(1,'x0', x00  , 1.1d0, -100.0d0, 100.0d0, ierr)
      call mnparm(2,'y0', y00  , 1.1d0, -100.0d0, 100.0d0, ierr)
      call mnparm(3,'co', co0,sco0,100.0d0, 5.0d6, ierr)
      call mnparm(4,'c1', c10,sc10,-1.5d0, 10.0d0, ierr)
      call mnparm(5,'c2', c20,sc20, -0.2d0, 0.2d0, ierr)
      call mnparm(6,'c3', c30,sc30, -0.02d0, 0.02d0, ierr)

**    call mncomd(frfit,'set pri 2',ierr,0)
      call mncomd(ldffit,'set pri 1',ierr,0)
      if(ierr.gt.0) write(6,*) '  set print: ierr=',ierr
      call mncomd(ldffit,'set str 2',ierr,0)
      if(ierr.gt.0) write(6,*) '  set strat: ierr=',ierr
      call mncomd(ldffit,'show eps',ierr,0)
      if(ierr.gt.0) write(6,*) '  show eps: ierr=',ierr

      call mncomd(ldffit,'migrad',ierr,0)
      if(ierr.gt.0) write(6,*) '  migrad: ierr=',ierr
*            keycnt = 1                !  register 1st core location error
      call mncomd(ldffit,'exit',ierr,0)
      if(ierr.gt.0) write(6,*) '   exit: ierr=',ierr
*            keycnt = 0


      if(l.ge.6) then
        go to 79
      endif
      if(fff.gt.100.0d0) then
        ff(l) = fff
	go to 77
      endif

   79   d1 = d1 + delr
        d2 = d2 + delr*delr

        write(30,100) delr
      write(*,*)  ''
      write(*,*)  '         delr=',delr
      write(*,*)  ''

        do i=1,5

          rmax = 100.0*i

	  call qaout(ch,0.0,rmax,0.0,1.0e-5,est,
     *                err,nf,fl)
          eest = 0.25d0*est       !  integration result
	  ee(i) = eest

	  if(i.eq.1) then
	    atot100 = atot100 + chtotsq(kev(k))
	    stot100 = stot100 + (eest - chtotsq(kev(k)))**2
	  endif
	  if(i.eq.3) then
	    atot300 = atot300 + chtot(kev(k))
	    stot300 = stot300 + (eest - chtot(kev(k)))**2
	  endif

*      write(*,*)  ''
          write(*,*)  '   rmax=',rmax,' integral=',eest

        enddo
      write(*,*)  ''
        write(31,101) ee

      enddo            !  loop in configurations ends

      close(30)
      close(31)

      write(*,*)  ''
      write(*,*)  ''
      d1 = d1/ncnf
      d2 = sqrt(d2/ncnf)
      write(*,*) '   d1=',d1,'  d2=',d2
      write(*,*)  ''

      atot100 = atot100/100
      atot300 = atot300/100
      stot100 = sqrt(stot100/100)
      stot300 = sqrt(stot300/100)
      write(*,*)  '    atot100=',atot100,'  atot300=',atot300
      write(*,*)  '    stot100=',stot100,'  stot300=',stot300
      write(*,*)  '    dtot100=',stot100/atot100,'  dtot300=',
     *                stot300/atot300
      write(*,*)  ''

      stop
*  100 format(1pd12.5)
  100 format(f10.5)
  101 format(5(1x,f10.1))
  301 format(1x,i4,1x,i2,1x,i2,1x,i2,1x,f5.1,1x,f5.1)
  400 format(5(1x,i8))
      end

***************************************************************************

      subroutine ldffit(npar,g,f,z,ifl,futil)
*      implicit double precision (a-h,o-z)
      integer ifl,npar, ix,iy, key
      common/print/ key
      double precision g,f,z,      fff
      dimension z(npar),g(npar)
      common/func/ fff

      integer nsq(5,5)
      real*8 xx(5), yy(5), x0,y0,  xa,ya, delr
      common/expdat/ xx,yy, x0,y0,  xa,ya, delr, nsq

*      real*8 co,c1,c2,c3
*      common/norma/ co,c1,c2,c3
      real*8  co,c1,c2,c3,eest
      common/energy/ co,c1,c2,c3,eest

      real*8  dx,dy,rho,add,chldf

      external futil

      save

      if(ifl.eq.3) then
        x0 = z(1)
        y0 = z(2)
        co = z(3)
        c1 = z(4)
        c2 = z(5)
        c3 = z(6)
	delr = sqrt((x0-xa)**2 + (y0-ya)**2)
*        write(*,*)  '   x0=',x0,' y0=',y0,' delr=',delr,
*     *               '   xa=',xa,' ya=',ya,'  fff=',fff
        write(*,100)  x0,y0,delr,xa,ya,fff
      endif


      f = 0.0d0
      x0 = z(1)
      y0 = z(2)
      co = z(3)
      c1 = z(4)
      c2 = z(5)
      c3 = z(6)

      do ix=1,5
	dx = x0 - xx(ix)
        do iy=1,5
	  dy = y0 - yy(iy)
	  rho = sqrt(dx**2 + dy**2)
          chldf = 1.0d0 * chld(sngl(rho))
	  add = (nsq(ix,iy) - chldf)**2/dabs(chldf)
**      write(*,*)  ' =>  rho=',rho,' chldf=',chldf,' add=',add
	  f = f + add
	enddo
      enddo
      fff = f
**      write(*,*)  ' x0=',x0,' y0=',y0,' co=',co,
**     *            ' c1=',c1,' c2=',c2,' c3=',c3,' f=',f

      return
  100 format('   x0=',1pd11.4,' y0=',1pd11.4,' delr=',1pd11.4,
     *       '   xa=',1pd11.4,' ya=',1pd11.4,'  fff=',1pd11.4)
      end

***************************************************************************

      real function chld(rperp)
      real*8  a0,a1,a2,a3,eest,        ch
      common/energy/ a0,a1,a2,a3,eest
*      real*8 co,c1,c2,c3
*      common/norma/ co,c1,c2,c3

      real rperp

      save

      ch = sngl(1.0 + rperp*(a1 + rperp*(a2 + rperp*a3)))
**      write(*,*)  ' ==>> ch=',ch
*      ch = rperp*sngl(a0)/ch   ! 1PeV proton, vert
      chld = sngl(a0)/ch   ! 1PeV proton, vert
**      write(*,*)  '     ==>> chld=',chld

      return
      end

***************************************************************************

      real function ch(rperp)
      real*8  a0,a1,a2,a3,eest
      common/energy/ a0,a1,a2,a3,eest
*      real*8 co,c1,c2,c3
*      common/norma/ co,c1,c2,c3

      real rperp

      save

      ch = sngl(1.0 + rperp*(a1 + rperp*(a2 + rperp*a3)))
*      ch =  2.0*3.141592654 * rperp * rperp * sngl(a0)/ch   ! 1PeV proton, vert
      ch =  2.0*3.141592654 * rperp * sngl(a0)/ch   ! 1PeV proton, vert

      return
      end

************************************************************************************************************

       logical function intrac()

c       logical intrac
       intrac = .false.

       return
       end

C     -------------------------------------------------------------------------------------------------------

      SUBROUTINE QAOUT(FUN,A,B,ABSERR,RELERR,RESULT,ERREST,NOFUN,FLAG)                                                             
      DIMENSION QRIGHT(31),F(8),X(8),FSAVE(4,30),XSAVE(4,30)          
      LEVMIN=1                                                        
      LEVMAX=30                                                       
      LEVOUT=6                                                        
      NOMAX=5000                                                      
      NOFIN=NOMAX-8*(LEVMAX-LEVOUT+2**(LEVOUT+1))                     
      W0=14./45.                                                      
      W1=64./45.                                                      
      W2=24./45                                                       
C     W3=41984./14175.                                                
C     W4=-18160./14175.                                               
      FLAG=0.                                                         
      RESULT=0.                                                       
      COR11=0.                                                        
      ERREST=0.                                                       
      AREA=0.                                                         
      NOFUN=0                                                         
      IF(A.EQ.B) RETURN                                               
      LEV=0                                                           
      NIM=1                                                           
      X0=A                                                            
      X(8)=B                                                          
      QPREV=0.                                                        
      F0=FUN(X0)                                                   
      STONE=(B-A)*.125                                                
      X(4)=(X0+X(8))*.5                                               
      X(2)=(X0+X(4))*.5                                               
      X(6)=(X(4)+X(8))*.5                                             
      X(1)=(X0+X(2))*.5                                               
      X(3)=(X(2)+X(4))*.5                                             
      X(5)=(X(4)+X(6))*.5                                             
      X(7)=(X(6)+X(8))*.5                                             
      DO 25 J=2,8,2                                                   
         F(J)=FUN(X(J))                                            
   25 CONTINUE                                                        
      NOFUN=9                                                         
   30 X(1)=(X0+X(2))*.5                                               
      F(1)=FUN(X(1))                                               
      DO 35 J=3,7,2                                                   
         X(J)=(X(J-1)+X(J+1))*.5                                     
         F(J)=FUN(X(J))                                           
   35 CONTINUE                                                       
      NOFUN=NOFUN+4                                                  
      STEP=(X(8)-X0)*.125                                            
      QLEFT=(W0*(F0+F(4))+W1*(F(1)+F(3))+W2*F(2))*STEP               
C    1W4*F(4))*STEP                                                  
      QRIGHT(LEV+1)=(W0*(F(4)+F(8))+W1*(F(5)+F(7))+W2*F(6))*STEP     
C    1W3*(F(11)+F(13))+W4*F(12))*STEP                                
      QNOW=QLEFT+QRIGHT(LEV+1)                                       
      QDIFF=QNOW-QPREV                                               
      AREA=AREA+QDIFF                                                
      ESTERR=ABS(QDIFF)/64.                                          
      TOLERR=AMAX1(ABSERR,RELERR*ABS(AREA))*(STEP/STONE)             
      IF(LEV.LT.LEVMIN) GO TO 50                                     
      IF(LEV.GE.LEVMAX) GO TO 62                                     
      IF(NOFUN.GT.NOFIN) GO TO 60                                    
      IF(ESTERR.LE.TOLERR) GO TO 70                                  
   50 NIM=2*NIM                                                      
      LEV=LEV+1                                                      
      DO 52 I=1,4                                                    
            FSAVE(I,LEV)=F(I+4)                                      
            XSAVE(I,LEV)=X(I+4)                                      
   52 CONTINUE                                                       
      QPREV=QLEFT                                                    
      DO 55 I=1,4                                                    
         J=-I                                                        
         F(2*J+10)=F(J+5)                                            
         X(2*J+10)=X(J+5)                                            
   55 CONTINUE                                                       
      GO TO 30                                                       
   60 NOFIN=2*NOFIN                                                  
      LEVMAX=LEVOUT                                                    
      FLAG=FLAG+(B-X0)/(B-A)                                           
      GO TO 70                                                         
   62 FLAG=FLAG+1.                                                     
   70 RESULT=RESULT+QNOW                                               
      ERREST=ERREST+ESTERR                                             
      COR11=COR11+QDIFF/63.                                            
   72 IF(NIM.EQ.2*(NIM/2)) GO TO 75                                    
      NIM=NIM/2                                                        
      LEV=LEV-1                                                        
      GO TO 72                                                         
   75 NIM=NIM+1                                                        
      IF(LEV.LE.0) GO TO 80                                            
      QPREV=QRIGHT(LEV)                                                
      X0=X(8)                                                          
      F0=F(8)                                                          
      DO 78 I=1,4                                                      
         F(2*I)=FSAVE(I,LEV)                                           
         X(2*I)=XSAVE(I,LEV)                                           
   78 CONTINUE                                                         
      GO TO 30                                                         
   80 RESULT=RESULT+COR11                                              
      IF(FLAG.LE.1.) GO TO 81                                          
      PRINT 90,FLAG,NOFUN                                              
C     WRITE (12,90) FLAG,NOFUN 
   90 FORMAT('     ** WARNING:FLAG=',E12.6,' NF=',I6,'    IN QAOUT **')
   81 RETURN                                                           
      END                                                              

*******************************************************************

      double precision FUNCTION FUTIL()
C
**      IMPLICIT NONE
C
C     Description:-
C     =============
C
C     This is essentially a dummy routine. The name of this
C     routine is passed in MINUIT in the fitting process as
C     and external (by FCNHV1 and FCN1D). This is called the
C     "utilitary routine". I am not sure what it is used for but
C     this version does not do anything.
C
C     Author:- unknown
C     ========
C
C     Creation Date: unknown
C     ======================
C
C     Revisions:-
C     ===========
C         Date    Name  Description
C         ----    ----  ----------------------------------------
C
C     Arguments:-
C     ===========
C
**      DOUBLE PRECISION FUTIL
C
C     Implicit inputs, outputs, side effects:-
C     ========================================
C
C        Called by: could be called by a number of routines in
C                   the MINUIT package. The name of this routine
C                   is passed to MINUIT as an external.
C        Calls    : (none)
C
C     Global Specifications:-
C     =======================
C
C     External Specifications:-
C     =========================
C
C     Local Specifications:-
C     ======================
C
C     Executable Statements:-
C     =======================
C
      FUTIL = 0.
C
      RETURN
      END                                   ! FUTIL.

**********************************************


