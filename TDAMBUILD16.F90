	  Program BUILDDRIVER
	  USE DAMINITIAL_T
	  USE DAMBUILD_T
	  implicit none
	  integer(kint)::iden,ierr,i,j,iw,jw,dum
	  real :: tarray(2), tiempo, dtime
	  logical,external:: accept
!      character, external:: ncase
	  namelist / options / lgradient, lderiv2, largo, lexact, &
				& lmaxi,filename,newfile,iswindows
                         
	  tiempo = dtime(tarray)
	  lgradient  = .true.
	  lderiv2    = .true.
	  largo      = .false.		! IF .TRUE. LONG-RANGE POTENTIAL
	  lexact     = .false.		! IF .TRUE. "EXACT" POTENTIAL IS TABULATED
	  iswindows  = .false.
	  lmaxi      = 15
	  filename   = ""
	  newfile    = ""
       
	  read(5,options)     ! Read from standard input
	  read(5,*) projectname

! Decides if the user has provided projectname with full path.
	  if (iswindows) then
		dirsep = "\\"
		i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
		if (i .eq. 0) then	! This is intended for MinGW, whose directory separator in windows is also /
			dirsep = "/"
			i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
		endif
	  else
		dirsep = "/"
		i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
	  end if  
	  if (len_trim(filename).eq.0) then
		filename	 = projectname
	  else
		filename = projectname(1:i)//trim(filename)
	  endif
	  if (iswindows) then
		filenamewin = filename
		do j = 1, len_trim(filenamewin)
			i = index(filenamewin,'/') ; if (i == 0) exit
			filenamewin = filenamewin(:i-1) // '\' // filenamewin(i+1:)
		end do
	  endif
      call predampot
      call system ("rm -f initialgeo.xyz finalgeo.xyz")
!Read the new file and convert the coordinates in au
      open(unit=301, file = trim(adjustl(newfile)),iostat=ierr)
	  if (ierr.ne.0) call error(1,"error reading guest file geometry")
      read(301,*,err=302,end=302) nb
      read(301,*,err=302,end=302)
      allocate (atwtn(nb),bs(nb),bx(nb),by(nb),bz(nb),bq(nb),atvdwn(nb),stat=ierr)
      allocate (bxl(nb),byl(nb),bzl(nb),stat=ierr)
      allocate (dist(ncen),stat=ierr)
      do i = 1,nb
            read(301,*,err=302,end=302) bs(i),bx(i),by(i),bz(i), bq(i)
			bx(i)=bx(i)*angau
			by(i)=by(i)*angau
			bz(i)=bz(i)*angau
      end do
302   close(301)
        
      do i=1,nb
       call ncase(bs(i))   ! Case correction of Atomic symbol
      end do 

!_______Printing initial geometry______________________________________________________
      open(8,file="initialgeo.xyz",status="new",iostat=ierr)
	  if (ierr.ne.0) call error(1,"error final file of optimized geometry")
        write(8,*)nb+ncen
        write(8,*)  "Energy"             
        do i = 1,ncen
            write(8,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
        end do
        do i=1,nb
            write(8,"(a2,3f15.10)") bs(i),bx(i)*auang,by(i)*auang,bz(i)*auang
        end do
        close(8)
!--------------------------------------------------------------------------------------

! Allocate an array for atomic mass and van der Waals radii of atoms. 
! atmnms are atomic names defined in global module of size 103. If the sym of the 
! atom in the new file matches with the declared atomic symbol, the 
! local variable for at wt atwtn and at vdw atvdwn are updated.  
        do i=1,nb
        	do j =1,103
            	if (trim(adjustl(bs(i))).ne.trim(adjustl(atmnms(j)))) cycle
            	if (trim(adjustl(bs(i))).eq.trim(adjustl(atmnms(j)))) then
                	atwtn(i)  = atmwts(j)
                	atvdwn(i) = atmvdw(j)*0.8d0
                	exit
            	end if
        	end do
        end do
189   continue
      open(11,file="visualize.xyz")           
      Enj=0.0d0 ! Initialize Electrostatic interaction energy 
      flag=0
      do 
          dum=1
          do i=1,nb
              bxl(i)=bx(i)
              byl(i)=by(i)
              bzl(i)=bz(i)
          end do 
          if (accept(dum)) call ROTATION
           write(11,*) nb
           write(11,*) Enj  
           do i=1,nb           
      	    write(11,"(a2,3f15.10)") bs(i),bx(i)*auang,by(i)*auang,bz(i)*auang
          end do
          if (accept(dum)) call TRANSLATION
           write(11,*) nb
           write(11,*) Enj  
           do i=1,nb           
      	  write(11,"(a2,3f15.10)") bs(i),bx(i)*auang,by(i)*auang,bz(i)*auang
          end do
          if (accept(dum).eqv..false.) exit    
	  if (flag==-1) exit
      end do   
      call ROTATION
      open(18,file="finalgeo.xyz",status="new",iostat=ierr)
      write(18,*) nb+ncen
      write(18,*) Enj             
      do iw = 1,ncen
        write(18,"(a2,3f15.10)") atmnms( int(zn(iw))), (rcen(jw,iw)*0.529177249d0,jw=1,3)
      end do
      do iw=1,nb
      	write(18,"(a2,3f15.10)") bs(iw),bx(iw)*auang,by(iw)*auang,bz(iw)*auang
      end do
      close(18)
      close(11)
      tiempo = dtime(tarray)
      write(iout,*)""
	  write(iout,"(1x,'Timing in seconds of processor ',' (user, system, total):',/5x,'(', &
            &e12.5,',',e12.5,',',e12.5')')") &
            &tarray(1), tarray(2), tarray(1)+tarray(2)
	  write(iout,*)""
      close(iden)
      close(iout)
      stop
      end

!******************************************
      SUBROUTINE ROTATION  
!******************************************      
      USE DAMBUILD_T
      IMPLICIT NONE 
      REAL(KREAL), ALLOCATABLE :: enjr(:,:,:)
      REAL(KREAL), ALLOCATABLE :: btx(:),bty(:),btz(:)
      REAL(KREAL), ALLOCATABLE :: rcom(:),rcox(:),rcoy(:),rcoz(:)
      REAL(KREAL), DIMENSION(3,3) :: rt
      REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
      REAL(KREAL) :: theta,phi,ang,ua,va,wa,eval
      INTEGER(KINT), ALLOCATABLE :: ex(:)
      INTEGER(KINT) :: i,j,k,l,a1,a2,a3,i1,j1,k1,as,dm,is,js,ierr,dum
      LOGICAL,external:: accept
!--------------------------------------------------
	  call COMASS
      allocate (rcox(nb),rcoy(nb),rcoz(nb),rcom(nb),stat=ierr)
	  do i=1,nb
	    rcox(i)=bx(i)-cox
	    rcoy(i)=by(i)-coy
	    rcoz(i)=bz(i)-coz
	    rcom(i)=dsqrt(rcox(i)**2+rcoy(i)**2+rcoz(i)**2)
	  end do
	
	  as=20
	  a1=180/as+1   !Number of theta's 
	  a2=360/as     !Number of phi's
	  a3=360/as     !Number of rotation angle's
	  allocate (btx(nb),bty(nb),btz(nb),stat=ierr)
 	  allocate (enjr(a1,a2,a3),stat=ierr)
	  enjr = 0.0
!	  Choosing a rotation matrix of general form.
!	  Rotating the axis of rotation using theta and phi
!	  Rotating molecule around every axis using ang
	  do i=0,180,as
	    i1=i/as+1
	    theta = real(i)*p_i/180.0d0
	    do j=0,345,as
	       j1=j/as+1
	       if (i.eq.0.and.j.gt.0) exit	
	       if (i.eq.180.and.j.gt.0) exit	
	
	       phi = real(j)*p_i/180.0d0
	       ua=sin(theta)*cos(phi)
	       va=sin(theta)*sin(phi)
	       wa=cos(theta)

	       do k = 0,345,as
		    k1=k/as+1
		    ang = real(k)*p_i/180.0d0
		    rt(1,1) = ua**2+(1-ua**2)*cos(ang)
		    rt(1,2) = ua*va*(1-cos(ang)) - wa*sin(ang)
		    rt(1,3) = ua*wa*(1-cos(ang)) + va*sin(ang)
		    rt(2,1) = ua*va*(1-cos(ang)) + wa*sin(ang)
		    rt(2,2) = va**2+(1-va**2)*cos(ang)
		    rt(2,3) = va*wa*(1-cos(ang)) - ua*sin(ang)
		    rt(3,1) = ua*wa*(1-cos(ang)) - va*sin(ang)
		    rt(3,2) = va*wa*(1-cos(ang)) + ua*sin(ang)
		    rt(3,3) = wa**2+(1-wa**2)*cos(ang)
		
!            write(11,*) nb
!	        write(11,*)	
		    do l=1,nb
	   	      btx(l) = cox + (rcox(l)*rt(1,1) + rcoy(l)*rt(1,2) + rcoz(l)*rt(1,3))
	   	      bty(l) = coy + (rcox(l)*rt(2,1) + rcoy(l)*rt(2,2) + rcoz(l)*rt(2,3))
	   	      btz(l) = coz + (rcox(l)*rt(3,1) + rcoy(l)*rt(3,2) + rcoz(l)*rt(3,3))
!              write(11,*) bs(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
        	end do	
            do l=1,nb
                bxl(l)=btx(l) 
                byl(l)=bty(l) 
                bzl(l)=btz(l)
            end do 
            dum=2
             if (accept(dum).eqv..false.) then
	   	        enjr(i1,j1,k1)= +10000.0d0
                cycle
             end if    
		    do l=1,nb
    	      CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
		      & dyytot,dyztot,dzztot,btx(l),bty(l),btz(l))
		      enjr(i1,j1,k1)=enjr(i1,j1,k1)+vtot*bq(l)
	!		print*,"vtot,bq(l)"
	!		print*,vtot,bq(l)
!              write(6,*)"One move of rotation produces energy = ",enjr(i1,j1,k1)
		    end do

	      end do
	   end do
	end do

	dm=size(shape(enjr))
	allocate (ex(dm),stat=ierr)
	ex=minloc(enjr)
	eval=minval(enjr)
     !print*, "OK before checking energy in rotation"
    if ((Enj-eval).gt.1.0e-8) then  ! If previous energy is smaller than new, then update
        enj=eval
!	print*,"ok after checking energy in rotation"
!	    write(6,*)ex(1),ex(2),ex(3)
	    write(6,*) "Best energy of all rotations      = ",eval
	    i= (ex(1)-1)*as
	    theta = real(i)*p_i/180.0d0
	    j=(ex(2)-1)*as
	    phi = real(j)*p_i/180.0d0
	    ua=sin(theta)*cos(phi)
	    va=sin(theta)*sin(phi)
	    wa=cos(theta)
	    k=(ex(3)-1)*as
	    ang = real(k)*p_i/180.0d0

	    rt(1,1) = ua**2+(1-ua**2)*cos(ang)
	    rt(1,2) = ua*va*(1-cos(ang)) - wa*sin(ang)
	    rt(1,3) = ua*wa*(1-cos(ang)) + va*sin(ang)
	    rt(2,1) = ua*va*(1-cos(ang)) + wa*sin(ang)
	    rt(2,2) = va**2+(1-va**2)*cos(ang)
	    rt(2,3) = va*wa*(1-cos(ang)) - ua*sin(ang)
	    rt(3,1) = ua*wa*(1-cos(ang)) - va*sin(ang)
	    rt(3,2) = va*wa*(1-cos(ang)) + ua*sin(ang)
	    rt(3,3) = wa**2+(1-wa**2)*cos(ang)

	    do l=1,nb
	      bx(l) = cox + (rcox(l)*rt(1,1) + rcoy(l)*rt(1,2) + rcoz(l)*rt(1,3))
	      by(l) = coy + (rcox(l)*rt(2,1) + rcoy(l)*rt(2,2) + rcoz(l)*rt(2,3))
	      bz(l) = coz + (rcox(l)*rt(3,1) + rcoy(l)*rt(3,2) + rcoz(l)*rt(3,3))
        end do
    else 
	flag=-1 
!	return
    end if 
	deallocate (btx,bty,btz,ex,enjr)
    deallocate (rcox,rcoy,rcoz,rcom)
	return
	end subroutine 	

!******************************************
!> @details Creates guess points for MED and MESP topography.
      SUBROUTINE TRANSLATION 
!******************************************      
      USE DAMBUILD_T
      IMPLICIT NONE 
      REAL(KREAL), ALLOCATABLE:: enjt(:,:,:)
      REAL(KREAL), ALLOCATABLE:: btx(:),bty(:),btz(:)
      REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
      REAL(KREAL) :: theta,phi,ang,ua,va,wa,eval
      INTEGER(KINT) :: i,j,k,l,n,ierr,as,np,nt,dm,is,js,dum
      INTEGER(KINT) :: a1,a2,i1,j1
      INTEGER(KINT), ALLOCATABLE :: ex(:)
      LOGICAL,EXTERNAL:: accept
		
	  if (.not.allocated(btx)) allocate (btx(nb),stat=ierr)
	  if (.not.allocated(bty)) allocate (bty(nb),stat=ierr)
	  if (.not.allocated(btz)) allocate (btz(nb),stat=ierr)
          
	   allocate (enjt(3,3,3),stat=ierr)
	   enjt = 0.0d0
        do i= -1,1
          do j=-1,1
            do k=-1,1
!             if (i.eq.0.and.j.eq.0.and.k.eq.0) cycle
!             write(11,*) nb
!             write(11,*)
             do l = 1,nb
                btx(l)= bx(l)+i*0.5 
                bty(l)= by(l)+j*0.5
                btz(l)= bz(l)+k*0.5
!                write(11,*) bs(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
	         end do

             do l=1,nb
                bxl(l)=btx(l)
                byl(l)=bty(l)
                bzl(l)=btz(l)
             end do 
             
             dum=3
             if (accept(dum).eqv..false.) then
	   	        enjt(i+2,j+2,k+2)= +10000.0d0
                cycle
             end if    
		     do l=1,nb
    	   	    CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
	   	        & dyytot,dyztot,dzztot,btx(l),bty(l),btz(l))
	   	        enjt(i+2,j+2,k+2)=enjt(i+2,j+2,k+2)+vtot*bq(l)
	!		print*,"vtot,bq(l)"
	!		print*,vtot,bq(l)
!                write(6,*)"One valid move shows energy= ",enjt(i1,j1)
		     end do
           end do   
         end do
       end do
       dm=size(shape(enjt))
	   allocate (ex(dm),stat=ierr)
	   ex=minloc(enjt)
	   eval=minval(enjt)
       if ((Enj-Eval).gt.1.0e-8) then
     ! print*, "Ok in translation"
        Enj=eval
!	    write(6,*) ex(1),ex(2),ex(3)
	    write(6,*) "Best energy of all translatory moves", eval
        i=ex(1)-2
        j=ex(2)-2
        k=ex(3)-2
 
        do l = 1,nb
          bx(l)= bx(l)+i*0.5
          by(l)= by(l)+j*0.5
          bz(l)= bz(l)+k*0.5
        end do
       else 
	flag=-1
	end if
	  deallocate (btx,bty,btz,ex,enjt)
      return 
      end subroutine

!=================================
!    Center of Mass
!=================================
    SUBROUTINE COMASS
	USE DAMBUILD_T
	IMPLICIT NONE
	INTEGER:: i
	
        cox=0.0d0; coy=0.0d0; coz=0.0d0; tmass=0.0d0
      	do i=1,nb
           tmass=tmass+atwtn(i)
           cox=cox+atwtn(i)*bx(i)
           coy=coy+atwtn(i)*by(i)
           coz=coz+atwtn(i)*bz(i)
        end do
        cox=cox/tmass; coy=coy/tmass; coz=coz/tmass
	return
	end subroutine	
    
!==================================
! Case correction of Atomic Symbol
!==================================
    SUBROUTINE NCASE(CHSYM)
    IMPLICIT NONE
    CHARACTER (LEN=*) , INTENT(IN OUT) :: CHSYM
    INTEGER                            :: I,IC,NLEN

    NLEN = LEN(TRIM(ADJUSTL(CHSYM)))
    IF (NLEN.EQ.2) THEN
        IC = ICHAR(CHSYM(1:1))
        IF (IC .GE. IACHAR("a") .AND. IC .LE. IACHAR("z")) CHSYM(1:1) = CHAR(IC-32)
        IC = ICHAR(CHSYM(2:2))
        IF (IC .GE. IACHAR("A") .AND. IC .LE. IACHAR("Z")) CHSYM(2:2) = CHAR(IC+32)
    ELSE IF (NLEN.EQ.1) then
        IC = ICHAR(CHSYM(1:1))
        IF (IC .GE. IACHAR("a") .AND. IC .LE. IACHAR("z")) CHSYM(1:1) = CHAR(IC-32)
    END IF
    RETURN    
    END 

!====================================
! Van der Waals overlap check
!====================================
	LOGICAL FUNCTION ACCEPT(DUM)
	USE DAMBUILD_T
    IMPLICIT NONE
    INTEGER(KINT):: ju,ku,ierr,dum

    ACCEPT=.TRUE.
      DO JU=1,NCEN
         DO KU=1,NB
      	    DIST(JU) = SQRT((RCEN(1,JU)-BXL(KU))**2+(RCEN(2,JU)-BYL(KU))**2+(RCEN(3,JU)-BZL(KU))**2)
            IF ((ATVDW(JU)+ATVDWN(KU)).GT.DIST(JU)) THEN
                ACCEPT=.FALSE.
                GO TO 98 
            END IF   
         END DO
98    END DO       
    RETURN
    END 


!******************************************
!> @details Creates guess points for MED and MESP topography.
!      SUBROUTINE TUNE 
!******************************************      
!      USE DAMINITIAL_T
!      USE DAMBUILD_T
!      IMPLICIT NONE 
!      REAL(KREAL), ALLOCATABLE:: enjt(:),rcom(:),rcox(:),rcoy(:),rcoz(:)
!      REAL(KREAL), ALLOCATABLE:: btx(:),bty(:),btz(:)
!      REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
!      REAL(KREAL):: eval
!      REAL(KREAL), DIMENSION(2):: frvcom,frvxcom,frvycom,frvzcom
!      INTEGER(KINT):: i,j,k,l,n,ierr,as,np,nt,dm,is,js,nopt
!      CHARACTER::dummy*10
!      INTEGER(KINT):: ierr,ig,newpt
!!      REAL(KREAL):: xmin,ymin,zmin,xmax,ymax,zmax
!      CHARACTER(2), ALLOCATABLE:: cnsym(:)
!      REAL(KREAL), ALLOCATABLE:: Vt(:), drvxt(:), drvyt(:), drvzt(:)
!	  REAL(KREAL), ALLOCATABLE:: xnpt(:),ynpt(:),znpt(:), atwtn(:),atvdwn(:),qnpt(:),dist(:)
!!      REAL(KREAL) :: X1,Y1,Z1,den,vtot,drvxtot,drvytot,drvztot,
!      REAL(KREAL) :: dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot,enerji
!!      REAL(KREAL) :: Sr,bx,by,bz,cdx,cdy,cdz
!      REAL(KREAL):: cnox,cnoy,cnoz,tnmass
!	  REAL(KREAL),ALLOCATABLE::	frvxcom(:),frvycom(:),frvzcom(:)
!	  INTEGER(KINT):: i,j,k
!		allocate (vt(newpt),drvxt(newpt),drvyt(newpt),drvzt(newpt),stat=ierr) 
!		
!	  frvxcom(1)=1000.0
!	  frvycom(1)=1000.0
!	  frvzcom(1)=1000.0
!	  frvxcom(2)=0.0E+0
!	  frvycom(2)=0.0E+0
!	  frvzcom(2)=0.0E+0
!
!        allocate (frvxcom(10000),frvycom(10000),frvzcom(10000),stat=ierr)
!        frvxcom=0.0E+0;frvycom=0.0E+0;frvzcom=0.0E+0;enerji=0.0E+0
!        do i =1,10000
!            do j=1,ncen
!                do k=1,newpt
!                    dist(j) = sqrt((rcen(1,j)-xnpt(k))**2+(rcen(2,j)-ynpt(k))**2+(rcen(3,j)-znpt(k))**2)
!                    if ((atvdw(j)+atvdwn(k)).ge.dist(j)) go to 168
!                end do
!            end do
!		    do ig = 1,newpt
!			    CALL DAMPOT(vt(ig),drvxt(ig),drvyt(ig),drvzt(ig),&
!			    &			dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xnpt(ig),ynpt(ig),znpt(ig))
!		    end do
!
!		    do ig = 1,newpt
!			    frvxcom(i)=frvxcom(i)+drvxt(ig)*qnpt(ig)
!			    frvycom(i)=frvycom(i)+drvyt(ig)*qnpt(ig)
!			    frvzcom(i)=frvzcom(i)+drvzt(ig)*qnpt(ig)
!                enerji=enerji+vt(ig)*qnpt(ig)
!		    end do
!            write(7,*)frvxcom(i),frvycom(i),frvzcom(i),enerji
!            if (i.eq.1) then
!                do ig = 1,newpt
!                    xnpt(ig)=(xnpt(ig)+frvxcom(i)*100)
!                    ynpt(ig)=(ynpt(ig)+frvycom(i)*100)
!                    znpt(ig)=(znpt(ig)+frvzcom(i)*100)
!                end do
!            end if
!            if (abs(frvxcom(i)).le.abs(frvxcom(i-1)).and.i.gt.1.and.abs(frvxcom(i)).gt.10e-9) then
!                do ig = 1,newpt
!                    xnpt(ig)=(xnpt(ig)+frvxcom(i)*100)
!                end do
!            end if    
!            if (abs(frvycom(i)).le.abs(frvycom(i-1)).and.i.gt.1.and.abs(frvycom(i)).gt.10e-9) then
!                do ig = 1,newpt
!                    ynpt(ig)=(ynpt(ig)+frvycom(i)*100)
!                end do
!            end if    
!            if (abs(frvzcom(i)).le.abs(frvzcom(i-1)).and.i.gt.1.and.abs(frvzcom(i)).gt.10e-9) then
!                do ig = 1,newpt
!                    znpt(ig)=(znpt(ig)+frvzcom(i)*100)
!                end do
!            end if    
!        end do
!168     continue     
!
!	do while (abs(frvxcom(2)).le.abs(frvxcom(1)).or. &
!	       &  abs(frvycom(2)).le.abs(frvycom(1)).or. &
!	       &  abs(frvzcom(2)).le.abs(frvzcom(1)))
!	   do l=1,nb
!    	   	CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
!	   	    & dyytot,dyztot,dzztot,bx(l),by(l),bz(l))
!	   	    frvxcom(2)=frvxcom(2)+drvx*bq(l)
!	   	    frvycom(2)=frvycom(2)+drvy*bq(l)
!	   	    frvzcom(2)=frvzcom(2)+drvz*bq(l)
!	   	    enj=enj+vtot*bq(l)
!	   end do
!!	frvcom(2)=sqrt(frvxcom(2)**2+frvycom(2)**2+frvzcom(2)**2)
!
!	   if (.not.allocated(btx)) allocate (btx(nb),stat=ierr)
!	   if (.not.allocated(bty)) allocate (bty(nb),stat=ierr)
!	   if (.not.allocated(btz)) allocate (btz(nb),stat=ierr)
!           if (abs(frvxcom(2)).le.abs(frvxcom(1)).and.abs(frvxcom(2)).gt.10e-9) then
!            	do i = 1,nb
!                   btx(i)=(bx(i)+frvxcom(2)*100)
!                end do
!           end if    
!           if (abs(frvycom(2)).le.abs(frvycom(1)).and.abs(frvycom(2)).gt.10e-9) then
!               	do i = 1,nb
!                   bty(i)=(by(i)+frvycom(2)*100)
!                end do
!           end if    
!           if (abs(frvzcom(2)).le.abs(frvzcom(1)).and.abs(frvzcom(2)).gt.10e-9) then
!            	do i = 1,nb
!                   btz(i)=(bz(i)+frvzcom(2)*100)
!                end do
!           end if    
!	   write(11,*) frvxcom(2),frvycom(2),frvzcom(2),frvcom(2),enj
!
!	   frvxcom(1)=frvxcom(2)
!	   frvycom(1)=frvycom(2)
!	   frvzcom(1)=frvzcom(2)
!
!           do j=1,ncen
!              do k=1,nb
!              	dist(j) = sqrt((rcen(1,j)-btx(k))**2+(rcen(2,j)-bty(k))**2+(rcen(3,j)-btz(k))**2)
!              	if ((atvdwn(j)+atvdwn(k)).ge.dist(j)) return
!              end do
!           end do
!
!	   do i=1,nb
!	      bx(i)=btx(i)
!	      by(i)=bty(i)
!          bz(i)=btz(i)
!	   end do
!
!	  end do
!      return 
!      end subroutine
!!
