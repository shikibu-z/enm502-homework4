      program fem
c
      implicit real*8 (a-h,o-z)
      common/mat1/a(1800,600),b(1800),t(1800),x(1800),y(1800),det,nm(4)
      common/mat2/dpdx(4),dpdy(4),bfn(4,9),dbfc(4,9),dbfe(4,9),xx,yy,
     $     ajac
      common/mat3/gapt(3),wo(3),tf(2,3),dtf(2,3),nint1,nlc,nuc,iband
      common/mat4/axis,bxis,xfact,yfact,c1,c2,c3
      common/mat5/nnx,nny,nxel,nyel,nntol,nodtol,nrhs
c     
c *** open input output files
c     
      open(1,file='input')
      open(2,file='output')
      open(61,file='xc.m')
      open(62,file='yc.m')
      open(63,file='sol.m')
c     
c *** read input data
c
      read(1,*)nxel,nyel
      read(1,*)niter,error
      read(1,*)c1,c2,c3
      read(1,*)axis,bxis,xfact,yfact
c       
c *** calculate various parameters
c
      nrhs=1
      nint1=3
      
      nnx=nxel+1
      nny=nyel+1
      nuc=nnx+1
      nlc=nnx+1
      iband=2*nlc+1
      nodtol=nnx*nny
      nntol=nxel*nyel
c
c *** Set Gauss Quadrature properties
c
      gapt(1)=0.77459666241483d0
      gapt(2)=0.d0
      gapt(3)=-0.77459666241483d0
         
      wo(1)=0.555555555555556d0
      wo(2)=0.888888888888889d0
      wo(3)=0.555555555555556d0
c
c *** Choose initial guess
c
      do 15 i=1,nodtol
         t(i)=0.5d0
 15   continue
c        
c *** print input data and initial guess - OUTPUT
c
      write(2,1)nntol,nodtol,nxel,nyel,nnx,nny
      write(2,2)niter,error,c1,c2,c3
      write(2,3)axis,bxis,xfact,yfact,nuc,nlc,iband,nint1
      write(2,*)(gapt(i),wo(i),i=1,nint1)
c
c *** initialize bilinear basis functions - call the FIRST PART of
c *** the subroutine BASIS which loads in the values of the basis functions
c *** at the gauss point positions in the UNIT ELEMENT
c
      l=0
      do 100 i=1,nint1
         do 100 j=1,nint1
            l=l+1
            call basis(1,l,gapt(i),gapt(j))
 100  continue
c
c *** initialize linear basis functions - use for boundary integrals
c *** in subroutine BOUND for flux conditions
c
      do 200 i=1,nint1
         tf(1,i)=0.5d0*(1.d0-gapt(i))
         tf(2,i)=0.5d0*(1.d0+gapt(i))
         dtf(1,i)=-0.5d0
         dtf(2,i)=0.5d0
 200  continue
c
c *** calculate the mesh points - see routine for comments
c
      call mesh

c**************************************************************************
c *** newton iteration loop - MAIN CODE LOOP (for steady state)
c**************************************************************************

      iter=0
 300  continue
      iter=iter+1
      write(2,6)iter
c
c *** form the JACOBIAN matrix a and RESIDUAL vector b
c
      call domi
c     
c *** insert boundary conditions
c
      call bound
c
c *** solve the system of linear equations - DO NOT CHANGE
c
      call band
c
c *** update the SOLUTION vector t and check convergence
c *** Compute the L2 and Linf norms of the update vector to determine
c *** whether convergence has been reached - criterion is ERROR
c *** If convergence is not achieved and # of NR iterations > ITER, also exit
c
      sum=0.d0
      error1=0.d0
      do 400 i=1,nodtol
         t(i)=t(i)-b(i)
         sum=sum+b(i)*b(i)
         error1=max(dabs(b(i)),error1)
 400  continue
 
      write(2,7)sum,error1,det
      if(sum.le.error)go to 500
      if(iter.lt.niter)go to 300
      write(2,8)
      go to 600
 500  continue
      write(2,9)
 600  continue
c
c *** Once the code gets here - WE ARE DONE (converged or not)
c *** print results to output files
c
      write(2,10)
      write(2,5)(t(i),i=1,nodtol)
c
c *** print results to output files - MATLAB FORMAT
c *** The output is: (i)   X vector
c ***                (ii)  Y vector
c ***                (iii) Solution vector
c

      do 1000 i=1,nodtol
         write(61,*) x(i)
         write(62,*) y(i)
         write(63,*) t(i)
 1000 continue
c
c *** Formatting Statements for Output Files
c

 1    format(/,1x,'____________________________________________',//,1x,
     1     '    input data ',//,1x,
     2     '____________________________________________',//,1x,
     3     'nntol =',i5,5x,'nodtol =',i5,//,1x,'nxel =',i5,5x,
     4     'nyel =',i5,//,1x,'nnx =',i5,5x,'nny =',i5,/)
 2    format(/,1x,'total number of iterations (niter) =',i5,//,1x,
     1     'error bounds (error) =',e15.8,//,1x,'c1 =',f12.6,5x,'c2 =',
     2     f12.6,//,1x,'c3 =',f12.6,/)
 3    format(/,1x,'axis =',f12.6,5x,'bxis =',f12.6,//,1x,
     1     'xfact =',f12.6,5x,'yfact =',f12.6,//,1x,
     2     'nuc =',i5,5x,'nlc ',i5,5x,'iband =',i5,5x,//,1x,
     3     'number of integration points =',i5,/)
 5    format(/,5(1x,f12.6,2x))
 6    format(/,1x,'newtons iteration number (iter) =',i5,/)
 7    format(/,1x,'sum =',e15.8,2x,'error1 =',e15.8,2x,'det =',
     1     e15.8,/)
 8    format(/,1x,'iterations has not converged',/)
 9    format(/,1x,'---- iterations has converged ----',/)
 10   format(/,1x,'--------------------------------------------',//,10x,
     1     'results',//,1x,
     2    '____________________________________________',/)
 444  format(12e15.7)
      end

c****************************************************************************
c******                      END OF MAIN CODE                           *****
c****************************************************************************

c******************************************
c *** Subroutine - Evaluate BASIS FUNCTIONS
c******************************************

      subroutine basis(itype,ig,c,e)
      implicit real*8 (a-h,o-z)
      common/mat1/a(1800,600),b(1800),t(1800),x(1800),y(1800),det,nm(4)
      common/mat2/dpdx(4),dpdy(4),bfn(4,9),dbfc(4,9),dbfe(4,9),xx,yy,
     $     ajac
      go to(1,2),itype
c
c *** if ITYPE = 1 then goto line 1 (just below), otherwise go to line 2
c
 1    continue
c
c *** calculate basis functions and partial dervatives
c *** at a given c,e in transformed coordinates
c
      bfn(1,ig)=(1.d0-c)*(1.d0-e)/4.d0
      bfn(2,ig)=(1.d0+c)*(1.d0-e)/4.d0
      bfn(3,ig)=(1.d0-c)*(1.d0+e)/4.d0
      bfn(4,ig)=(1.d0+c)*(1.d0+e)/4.d0
      dbfc(1,ig)=-(1.d0-e)/4.d0
      dbfc(2,ig)=(1.d0-e)/4.d0
      dbfc(3,ig)=-(1.d0+e)/4.d0
      dbfc(4,ig)=(1.d0+e)/4.d0
      dbfe(1,ig)=-(1.d0-c)/4.d0
      dbfe(2,ig)=-(1.d0+c)/4.d0
      dbfe(3,ig)=(1.d0-c)/4.d0
      dbfe(4,ig)=(1.d0+c)/4.d0
      return
c
c *** When ITYPE = 2 jump straight to this point
c
 2    continue
c
c *** calculate partial dervatives of transformation
c *** at a point c,e in nelem
c
      xx=0.d0
      yy=0.d0
      dxdc=0.d0
      dxde=0.d0
      dydc=0.d0
      dyde=0.d0
      do 3 i=1,4
         xx=xx+bfn(i,ig)*x(nm(i))
         yy=yy+bfn(i,ig)*y(nm(i))
         dxdc=dxdc+dbfc(i,ig)*x(nm(i))
         dxde=dxde+dbfe(i,ig)*x(nm(i))
         dydc=dydc+dbfc(i,ig)*y(nm(i))
         dyde=dyde+dbfe(i,ig)*y(nm(i))
 3    continue
c     
c *** calculate jacobian of transformation from UNIT to REAL element
c
      cjac=dxdc*dyde-dxde*dydc
      ajac=abs(cjac)
c
c *** calculate dervatives of basis functions in x,y coordinates
c
      do 4 i=1,4
         dpdx(i)=(dbfc(i,ig)*dyde-dbfe(i,ig)*dydc)/cjac
         dpdy(i)=(dbfe(i,ig)*dxdc-dbfc(i,ig)*dxde)/cjac
 4    continue
      return
      end

c**********************************************
c *** Subroutine - Evaluate MESH POINTS ON GRID
c**********************************************

      subroutine mesh
      implicit real*8 (a-h,o-z)
      common/mat1/a(1800,600),b(1800),t(1800),x(1800),y(1800),det,nm(4)
      common/mat4/axis,bxis,xfact,yfact,c1,c2,c3
      common/mat5/nnx,nny,nxel,nyel,nntol,nodtol,nrhs
c
c *** evaluate the mesh points at each ij node
c *** Origin is at (0,0)
c *** AXIS, BXIS = are the X and Y lengths
c *** XFACT, YFACT are refinement parameters - set to 1 for uniform grid
c
      do 20 i=1,nnx
         do 20 j=1,nny
            x(nnum(i,j))=axis*(dfloat(i-1)/dfloat(nxel))**xfact
            y(nnum(i,j))=bxis*(dfloat(j-1)/dfloat(nyel))**yfact
 20   continue

      return
      end
     
c**********************************************************
c *** FUNCTION - GIVEN (i,j Index, find single global index
c**********************************************************

      integer function nnum(i,j)
      implicit real*8 (a-h,o-z)
      common/mat5/nnx,nny,nxel,nyel,nntol,nodtol,nrhs
c
c *** calculate the node number of node i,j in element nelem
c
      nnum=(j-1)*(nxel+1)+i
      return
      end

c*************************************************
c *** Subroutine - ASSEMBLE JACOBIAN AND RESIDUALS
c*************************************************

      subroutine domi

      implicit real*8 (a-h,o-z)

      common/mat1/a(1800,600),b(1800),t(1800),x(1800),y(1800),det,nm(4)
      common/mat2/dpdx(4),dpdy(4),bfn(4,9),dbfc(4,9),dbfe(4,9),xx,yy,
     $     ajac
      common/mat3/gapt(3),wo(3),tf(2,3),dtf(2,3),nint1,nlc,nuc,iband
      common/mat5/nnx,nny,nxel,nyel,nntol,nodtol,nrhs

      dimension temp(4,4),temp1(4)
c
c *** initialize the jacobian matrix a and residuals vector b
c
      do 20 i=1,nodtol
         b(i)=0.d0
         do 10 j=1,iband
            a(i,j)=0.d0
 10      continue
 20   continue
c
c *** iterate over elements in domain - LEVEL 1
c
      jtally=0
      ifact=0
      do 100 nelem=1,nntol
c
c *** find the node numbers of element # NELEM
c
         jtally=jtally+1
         if(jtally .lt. nxel+1) go to 30
         jtally=1
         ifact=ifact+1
 30      continue
         nm(1)=ifact*(nxel+1)+jtally
         nm(2)=nm(1)+1
         nm(3)=(ifact+1)*(nxel+1)+jtally
         nm(4)=nm(3)+1
c
c *** initialize working matrices TEMP (elemental Jacobian) and 
c *** TEMP1 (elemental Residual) for element integration
c
         do 40 iw=1,4
            temp1(iw)=0.d0
            do 40 jw=1,4
               temp(iw,jw)=0.d0
 40      continue
c
c *** iterate over each gauss point in element - LEVEL 2
c
         kk=0
         do 60 i=1,nint1
            do 60 j=1,nint1
               kk=kk+1
c
c *** calculate dervatives of basis functions and transformation
c *** jacobian at the gaussian points in x,y coordinates
c
               call basis(2,kk,c,e)
               wet=wo(i)*wo(j)*ajac
c              
c *** calculate dependent variable and partial dervatives at
c *** the gaussian integration points
c     
               dtdx=0.d0
               dtdy=0.d0
               do 51 ii=1,4
                  dtdx=dtdx+t(nm(ii))*dpdx(ii)
                  dtdy=dtdy+t(nm(ii))*dpdy(ii)
 51            continue
c
c *** iterate over weighting functions - LEVEL 3
c
               do 60 iw=1,4
                  bifn=bfn(iw,kk)
                  dpix=dpdx(iw)
                  dpiy=dpdy(iw)
                  ti=dpix*dtdx+dpiy*dtdy
c                 
c *** form the working residual vector in nelem
c
                  temp1(iw)=temp1(iw)+ti*wet
c     
c *** iterate over each node - LEVEL 4
c
                  do 60 jw=1,4
                     dpjy=dpdy(jw)
                     dpjx=dpdx(jw)
                     ti=dpix*dpjx+dpiy*dpjy
c
c *** form the working jacobian matrix in element
c
                     temp(iw,jw)=temp(iw,jw)+ti*wet
 60      continue
c
c *** store the ELEMENTAL integration matrix and vector in the global 
c *** matrix a and vector b
c
         do 90 i=1,4
            irow=nm(i)
            b(irow)=b(irow)+temp1(i)
            do 90 j=1,4
               icol=nm(j)
               iloc=icol-irow+nlc+1
               a(irow,iloc)=a(irow,iloc)+temp(i,j)
 90      continue
c
c *** Finish loop over elements (LEVEL 1)
c
 100  continue

      return
      end

c*************************************************
c *** Subroutine - IMPOSE BOUNDARY CONDITIONS
c*************************************************

      subroutine bound
      implicit real*8 (a-h,o-z)
      common/mat1/a(1800,600),b(1800),t(1800),x(1800),y(1800),det,nm(4)
      common/mat3/gapt(3),wo(3),tf(2,3),dtf(2,3),nint1,nlc,nuc,iband
      common/mat4/axis,bxis,xfact,yfact,c1,c2,c3
      common/mat5/nnx,nny,nxel,nyel,nntol,nodtol,nrhs

      dimension p(2,2),r(2),n(2)

c***********************************************************************
c *** Insert NATURAL BC at the boundary x=1 (RIGHT)
c***********************************************************************

c$$$      do 40 ii=1,nyel
c$$$c
c$$$c *** find the node numbers of element ii *****
c$$$c
c$$$         n(1)=nnum(nnx,ii)
c$$$         n(2)=nnum(nnx,ii+1)
c$$$         ajacy=0.5d0*(y(n(2))-y(n(1)))
c$$$c     
c$$$c *** intialize working areas p and r *****
c$$$c
c$$$         do 10 i=1,2
c$$$            r(i)=0.d0
c$$$            do 10 j=1,2
c$$$               p(i,j)=0.d0
c$$$ 10      continue
c$$$c
c$$$c *** iterate over guassian points *****
c$$$c
c$$$         do 20 ig=1,nint1
c$$$            funt=t(n(1))*tf(1,ig)+t(n(2))*tf(2,ig)
c$$$c           
c$$$c *** iterate over each node *****
c$$$c
c$$$            do 20 i=1,2
c$$$               r(i)=r(i)+wo(ig)*tf(i,ig)*ajacy*(funt-c3)
c$$$c
c$$$c *** iterate over each weighting function *****
c$$$c
c$$$               do 20 j=1,2
c$$$                  p(i,j)=p(i,j)+wo(ig)*tf(i,ig)*tf(j,ig)*ajacy
c$$$ 20      continue
c$$$c
c$$$c *** store boundary conditions in global matrix *****
c$$$c
c$$$         do 30 i=1,2
c$$$            irow=n(i)
c$$$            b(irow)=b(irow)+r(i)
c$$$            do 30 j=1,2
c$$$               iloc=n(j)-irow+nlc+1
c$$$               a(irow,iloc)=a(irow,iloc)+p(i,j)
c$$$ 30      continue
c$$$ 40   continue

c*********************************************************************
c *** Insert NATURAL BC at the boundary Y=1 (TOP)
c*********************************************************************

c$$$      do 400 ii=1,nxel
c$$$c
c$$$c *** find the node numbers of element ii *****
c$$$c
c$$$         n(1)=nnum(ii,nny)
c$$$         n(2)=nnum(ii+1,nny)
c$$$         ajacy=0.5d0*(x(n(2))-x(n(1)))
c$$$c     
c$$$c *** intialize working areas p and r *****
c$$$c
c$$$         do 100 i=1,2
c$$$            r(i)=0.d0
c$$$            do 100 j=1,2
c$$$               p(i,j)=0.d0
c$$$ 100     continue
c$$$c
c$$$c *** iterate over guassian points *****
c$$$c
c$$$         do 200 ig=1,nint1
c$$$            funt=t(n(1))*tf(1,ig)+t(n(2))*tf(2,ig)
c$$$            xx=x(n(1))*tf(1,ig)+x(n(2))*tf(2,ig)
c$$$c
c$$$c *** iterate over each node *****
c$$$c
c$$$            do 200 i=1,2
c$$$               r(i)=r(i)+wo(ig)*tf(i,ig)*ajacy*(funt-c3)
c$$$c
c$$$c *** iterate over each weighting function *****
c$$$c
c$$$               do 200 j=1,2
c$$$                  p(i,j)=p(i,j)+wo(ig)*tf(i,ig)*tf(j,ig)*ajacy
c$$$ 200     continue
c$$$c
c$$$c *** store boundary conditions in global matrix *****
c$$$c
c$$$         do 300 i=1,2
c$$$            irow=n(i)
c$$$            b(irow)=b(irow)+r(i)
c$$$            do 300 j=1,2
c$$$               iloc=n(j)-irow+nlc+1
c$$$               a(irow,iloc)=a(irow,iloc)+p(i,j)
c$$$ 300     continue
c$$$ 400  continue

c**********************************************************************
c *** Insert ESSENTIAL  BC at the boundary Y=0 (BOTTOM)
c**********************************************************************

      irow=nnum(1,1)
      b(irow)=t(irow)-c2
      do 120 icol=1,iband
         a(irow,icol)=0.d0
 120  continue
      a(irow,nlc+1)=1.d0
     
      do 140 ii=1,nxel
         irow=nnum(ii+1,1)
         b(irow)=t(irow)-c2
         do 130 icol=1,iband
            a(irow,icol)=0.d0
 130     continue
         a(irow,nlc+1)=1.d0
 140  continue

c**********************************************************************
c *** Insert ESSENTIAL BC (=0) at the boundary Y=1 (TOP)
c**********************************************************************

      irow=nnum(1,nny)
      b(irow)=t(irow)-0.d0
      do 220 icol=1,iband
         a(irow,icol)=0.d0
 220  continue
      a(irow,nlc+1)=1.d0
     
      do 240 ii=1,nxel
         irow=nnum(ii+1,nny)
         b(irow)=t(irow)-0.d0
         do 230 icol=1,iband
            a(irow,icol)=0.d0
 230     continue
         a(irow,nlc+1)=1.d0
 240  continue

      return
      end

c**********************************************************************
c *** SUBROUTINE - BAND SOLVER (DO NOT TOUCH THIS!!!!!!!!!)
c**********************************************************************

      subroutine band
      implicit real*8 (a-h,o-z)
c
c------------------------------------------------------------------------
c
      common/mat1/c(1800,600),v(1800,1),t1(1800),t2(1800),t3(1800)
     &     ,det,n4(4)
      common/mat3/gapt(3),wo(3),tf(2,3),dtf(2,3),nint1,nlc,nuc,iband
      common/mat5/nnx,nny,nxel,nyel,nntol,nod,nvarg
c
c------------------------------------------------------------------------
c
      mc=nod
      nc=iband
      nv=nvarg
c
c     ***** initialize det  *****
c
      eps=1.0d-15
      icount=0
      det=1.d0
c
c	*** prepare the matrix c for processing by shifting undefined   ***
c	*** elements out of the upper-left-hand corner and inserting    ***
c	*** zeros in the lower right-hand-corner                        ***
c
      nc1=nc+1
      lr=nc1/2
      mr=lr-1
      ii=mc+1
      do 60 ir=1,mr
         ii=ii-1
         nr=lr-ir
         jj=nc1
         do 40 jr=1,nr
            do 20 jc=2,nc
               p=c(ir,jc)
               k=jc-1
               c(ir,k)=p
 20         continue
            c(ir,nc)=0.d0
            jj=jj-1
            c(ii,jj)=0.d0
 40      continue
 60   continue
c
c	*** use row operations to eleminate the lower triangular part of a ***
c	*** apply these same operations to the right-hand-sides            ***
c
      do 400 ic=1,mc
         ipiv=ic
         piv=c(ic,1)
         pivmax=dabs(piv)
         kr=ic+1
         if(kr.gt.lr)go to 140
c     
c	***** find the largest possible pivot in the current column.
c
         do 120	ir=kr,lr
            pivmag=dabs(c(ir,1))
            if(pivmax.ge.pivmag)go to 120
            ipiv=ir
            pivmax=pivmag
 120     continue
         piv=c(ipiv,1)
c
c	***** check the pivot magnitude *******
c
 140     continue
         if(pivmax.ge.eps)go to 150
c     
c	**** a nonzero pivot smaller than eps has been found ****
c
c        print*,' zero pivot ', piv
         return
c
c	**** if necessary , swap row ipiv with row ic ****
c
 150     continue
         if(ipiv.eq.ic)go to 200
         det=-det
         do 160 jc=1,nc
            t=c(ic,jc)
            c(ic,jc)=c(ipiv,jc)
            c(ipiv,jc)=t
 160     continue
         do 180 jv=1,nv
            t=v(ic,jv)
            v(ic,jv)=v(ipiv,jv)
            v(ipiv,jv)=t
 180     continue
c
c	***** update the determinant value *****
c
 200     continue
         if(dabs(det).ge.1.0d+25)icount=icount+1
         if(dabs(det).ge.1.0d+25)det=det/dabs(det)
         det=det*piv
c
c	***** normalize the pivot row *****
c
         theta=1./piv
         do 220 jc=2,nc
            c(ic,jc)=c(ic,jc)*theta
 220     continue
         do 240 jv=1,nv
            v(ic,jv)=v(ic,jv)*theta
 240     continue
c
c	*** eliminate the lower triangular elements in the current column ***
c
         if(kr.gt.lr)go to 400
         do 380 ir=kr,lr
            t=c(ir,1)
            do 230 jc=2,nc
               k=jc-1
               c(ir,k)=c(ir,jc)-t*c(ic,jc)
 230        continue
            c(ir,nc)=0.d0
            do 360 jv=1,nv
               v(ir,jv)=v(ir,jv)-t*v(ic,jv)
 360        continue
 380     continue
         if(lr.eq.mc)go to 400
         lr=lr+1
 400  continue
c
c	***** triangularization is complete *****
c
c	
c	***** back-substitute to compute the solution vector(s) ******
c
      kr=2
      lc=mc-1
      do 480 ic=1,lc
         iv=mc-ic
         do 460 jv=1,nv
            ii=iv
            do 440 jc=2,kr
               ii=ii+1
               v(iv,jv)=v(iv,jv)-c(iv,jc)*v(ii,jv)
 440        continue
 460     continue
         if(kr.eq.nc) go to 480
         kr=kr+1
 480  continue
c
c	****** the matrix equation is solved ******
c
      return
      end
