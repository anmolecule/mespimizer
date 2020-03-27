        Program BUILDDRIVER
        USE DAMINITIAL_T
        USE DAMBUILD_T
        implicit none
        integer(kint)::iden,ierr,i,j,k,iw,jw,dum,nbatms,nb,flag
        real :: tarray(2), tiempo, dtime
        integer:: nmols
        REAL(KREAL) :: Enj,tssize 
!       character, external:: ncase
        namelist / options / lgradient, lderiv2, largo, lexact, &
                & lmaxi,filename,newfile,iswindows,tssize
                         
        tiempo = dtime(tarray)
        lgradient  = .true.
        lderiv2    = .true.
        largo      = .false.        ! IF .TRUE. LONG-RANGE POTENTIAL
        lexact     = .false.        ! IF .TRUE. "EXACT" POTENTIAL IS TABULATED
        iswindows  = .false.
        lmaxi      = 15
        filename   = ""
        newfile    = ""
        tssize     = 0.5
       
        read(5,options)     ! Read from standard input
        read(5,*) projectname

! Decides if the user has provided projectname with full path.
        if (iswindows) then
            dirsep = "\\"
            i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
            if (i .eq. 0) then  ! This is intended for MinGW, whose directory separator in windows is also /
                dirsep = "/"
                i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
            endif
        else
            dirsep = "/"
            i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        end if  
        if (len_trim(filename).eq.0) then
            filename = projectname
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
        !print *,"Number of centers",ncen

        !Read the new file and convert the coordinates in au
        open(unit=301, file = trim(adjustl(newfile)),iostat=ierr)
        if (ierr.ne.0) call error(1,"error reading guest file geometry")
        nmols=0
        do
            read(301,*,err=302,end=302) nb
            nmols=nmols+1
            molecules(nmols)%natoms = nb
            read(301,*,err=302,end=302)
            allocate (molecules(nmols)%bawt(nb),molecules(nmols)%bs(nb),molecules(nmols)%bx(nb),&
            molecules(nmols)%by(nb),molecules(nmols)%bz(nb),molecules(nmols)%bq(nb),molecules(nmols)%bvdw(nb),stat=ierr)
            do i = 1,nb
                read(301,*,err=302,end=302) molecules(nmols)%bs(i),molecules(nmols)%bx(i),molecules(nmols)%by(i),&
                molecules(nmols)%bz(i), molecules(nmols)%bq(i)
                molecules(nmols)%bx(i)=molecules(nmols)%bx(i)*angau
                molecules(nmols)%by(i)=molecules(nmols)%by(i)*angau
                molecules(nmols)%bz(i)=molecules(nmols)%bz(i)*angau
            end do
        enddo
302     close(301)
        print *,"Number of molecules detected in input file ",nmols

        ! Allocate an array for atomic mass and van der Waals radii of atoms. 
        ! atmnms are atomic names defined in global module of size 103. If the sym of the 
        ! atom in the new file matches with the declared atomic symbol, the 
        ! local variable for at wt atwtn and at vdw atvdwn are updated.  
        do i=1,nmols
            do j = 1,molecules(i)%natoms
                call ncase(molecules(i)%bs(j))   ! Case correction of Atomic symbol
                do k =1,103
                    if (trim(adjustl(molecules(i)%bs(j))).ne.trim(adjustl(atmnms(k)))) cycle
                    if (trim(adjustl(molecules(i)%bs(j))).eq.trim(adjustl(atmnms(k)))) then
                        molecules(i)%bawt(j) = atmwts(k)
                        molecules(i)%bvdw(j) = atmvdw(k)*0.8d0
                        exit
                    end if
                end do
            enddo
        enddo

        call optimize(nmols,tssize)

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
     
        SUBROUTINE OPTIMIZE(nmols,tssize) 
        USE DAMINITIAL_T
        USE DAMBUILD_T
        IMPLICIT NONE
        CHARACTER(2), ALLOCATABLE:: bsl(:)
        REAL(KREAL), ALLOCATABLE:: bxl(:),byl(:),bzl(:),bql(:),bawtl(:),bvdwl(:)
        integer:: nmols
        REAL(KREAL) :: Enj,tssize 
        integer(kint)::ierr,i,j,k,nbatms,flag,icurr
        logical,external:: accept
        
        call system ("rm -f initialgeo.xyz finalgeo.xyz")
        Enj=0.0d0 ! Initialize Electrostatic interaction energy 
        call writeclustergeom("initialgeo.xyz",nmols,Enj)

        open(11,file="visualize.xyz")           
        do i = 1,nmols
            icurr=i
			nbatms= molecules(i)%natoms
            allocate (bsl(nbatms),bxl(nbatms),byl(nbatms),bzl(nbatms),bql(nbatms),bawtl(nbatms),bvdwl(nbatms),stat=ierr)
            flag=0
            do j=1,nbatms
				bsl(j)=molecules(i)%bs(j)
                bxl(j)=molecules(i)%bx(j)
                byl(j)=molecules(i)%by(j)
                bzl(j)=molecules(i)%bz(j)
                bql(j)=molecules(i)%bq(j)
                bawtl(j)=molecules(i)%bawt(j)
                bvdwl(j)=molecules(i)%bvdw(j)
            end do 
            if (accept(icurr,nbatms,bxl,byl,bzl,bvdwl)) &
                call ROTATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,bawtl,bvdwl,flag,enj)
            do j=1,nbatms
				molecules(i)%bs(j)=bsl(j)
                molecules(i)%bx(j)=bxl(j)
                molecules(i)%by(j)=byl(j)
                molecules(i)%bz(j)=bzl(j)
                molecules(i)%bq(j)=bql(j)
                molecules(i)%bawt(j)=bawtl(j)
                molecules(i)%bvdw(j)=bvdwl(j)
            end do 
            if (flag/=-1) then
                write(11,*) nbatms
                write(11,*) Enj  
                do j=1,nbatms           
                    write(11,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
				    molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                end do
            endif

            flag = 0
            do while (flag/=-1) 
                do j=1,nbatms
					bsl(j)=molecules(i)%bs(j)
                    bxl(j)=molecules(i)%bx(j)
                    byl(j)=molecules(i)%by(j)
                    bzl(j)=molecules(i)%bz(j)
                    bql(j)=molecules(i)%bq(j)
                    bawtl(j)=molecules(i)%bawt(j)
                    bvdwl(j)=molecules(i)%bvdw(j)
                end do 
                if (accept(icurr,nbatms,bxl,byl,bzl,bvdwl)) &
                    call TRANSLATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,bawtl,bvdwl,flag,tssize,enj)
                do j=1,nbatms
					molecules(i)%bs(j)=bsl(j)
                    molecules(i)%bx(j)=bxl(j)
                    molecules(i)%by(j)=byl(j)
                    molecules(i)%bz(j)=bzl(j)
                    molecules(i)%bq(j)=bql(j)
                    molecules(i)%bawt(j)=bawtl(j)
                    molecules(i)%bvdw(j)=bvdwl(j)
                end do 

                if (.not. (accept(icurr,nbatms,bxl,byl,bzl,bvdwl))) flag = -1
                if (flag/=-1) then
                    write(11,*) nbatms
                    write(11,*) Enj  
                    do j=1,nbatms           
                        write(11,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
					    molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                    end do
                endif
            end do

            flag=0
            do j=1,nbatms
				bsl(j)=molecules(i)%bs(j)
                bxl(j)=molecules(i)%bx(j)
                byl(j)=molecules(i)%by(j)
                bzl(j)=molecules(i)%bz(j)
                bql(j)=molecules(i)%bq(j)
                bawtl(j)=molecules(i)%bawt(j)
                bvdwl(j)=molecules(i)%bvdw(j)
            end do 
            call ROTATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,bawtl,bvdwl,flag,enj)
            do j=1,nbatms
				molecules(i)%bs(j)=bsl(j)
                molecules(i)%bx(j)=bxl(j)
                molecules(i)%by(j)=byl(j)
                molecules(i)%bz(j)=bzl(j)
                molecules(i)%bq(j)=bql(j)
                molecules(i)%bawt(j)=bawtl(j)
                molecules(i)%bvdw(j)=bvdwl(j)
            end do
            if (flag/=-1) then
                write(11,*) nbatms
                write(11,*) Enj  
                do j=1,nbatms           
                    write(11,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
				    molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                end do
            endif

            deallocate(bsl,bxl,byl,bzl,bql,bawtl,bvdwl)
        enddo
!
        close(11)
        call writeclustergeom("finalgeo.xyz",nmols,Enj)

        end subroutine

!!******************************************
        SUBROUTINE ROTATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,batwtl,batvdwl,flag,enj)
!******************************************      
        USE DAMBUILD_T
        IMPLICIT NONE 
        INTEGER::nbatms,flag,icurr
        REAL(KREAL):: enj
        REAL(KREAL), ALLOCATABLE :: enjr(:,:,:)
        REAL(KREAL), DIMENSION(nbatms):: btx,bty,btz
        REAL(KREAL), DIMENSION(nbatms):: rcom,rcox,rcoy,rcoz
        CHARACTER(2), DIMENSION(nbatms):: bsl
        REAL(KREAL), DIMENSION(nbatms):: bxl,byl,bzl,bql,batwtl,batvdwl
        REAL(KREAL), DIMENSION(3,3) :: rt
        REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
        REAL(KREAL) :: theta,phi,ang,ua,va,wa,eval
        REAL(KREAL) :: coxl,coyl,cozl
        INTEGER(KINT), ALLOCATABLE :: ex(:)
        INTEGER(KINT) :: i,j,k,l,a1,a2,a3,i1,j1,k1,as,dm,is,js,ierr,dum
        LOGICAL,external:: accept
!--------------------------------------------------
        call COMASS(nbatms,bxl,byl,bzl,batwtl,coxl,coyl,cozl)
        do i=1,nbatms
            rcox(i)=bxl(i)-coxl
            rcoy(i)=byl(i)-coyl
            rcoz(i)=bzl(i)-cozl
            rcom(i)=dsqrt(rcox(i)**2+rcoy(i)**2+rcoz(i)**2)
        end do
   
        as=20
        a1=180/as+1   !Number of theta's 
        a2=360/as     !Number of phi's
        a3=360/as     !Number of rotation angle's
!       allocate (btx(nb),bty(nb),btz(nb),stat=ierr)
        allocate (enjr(a1,a2,a3),stat=ierr)
        enjr = 0.0
!       Choosing a rotation matrix of general form.
!       Rotating the axis of rotation using theta and phi
!       Rotating molecule around every axis using ang
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
            
!                     write(11,*) nb
!                   write(11,*) 
                    !print *, nbatms
					!print *, ""
                    do l=1,nbatms
                        btx(l) = coxl + (rcox(l)*rt(1,1) + rcoy(l)*rt(1,2) + rcoz(l)*rt(1,3))
                        bty(l) = coyl + (rcox(l)*rt(2,1) + rcoy(l)*rt(2,2) + rcoz(l)*rt(2,3))
                        btz(l) = cozl + (rcox(l)*rt(3,1) + rcoy(l)*rt(3,2) + rcoz(l)*rt(3,3))
                     !   print *, bsl(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
                    end do  
                    do l=1,nbatms
                        bxl(l)=btx(l) 
                        byl(l)=bty(l) 
                        bzl(l)=btz(l)
                    end do 
                    dum=2
                    if (accept(icurr,nbatms,bxl,byl,bzl,batvdwl).eqv..false.) then
                        enjr(i1,j1,k1)= +10000.0d0
                        cycle
                    end if    
                    do l=1,nbatms
                        CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
                        & dyytot,dyztot,dzztot,btx(l),bty(l),btz(l))
                        enjr(i1,j1,k1)=enjr(i1,j1,k1)+vtot*bql(l)
!                       print*,"vtot,bq(l)"
!                       print*,vtot,bq(l)
                        !write(6,*)"One move of rotation produces energy = ",enjr(i1,j1,k1)
                    end do
                end do
            end do
        end do

        dm=size(shape(enjr))
        allocate (ex(dm),stat=ierr)
        ex=minloc(enjr)
        eval=minval(enjr)
        if (eval<enj .and. (enj-eval).gt.1e-6) then
            enj=eval
!           write(6,*)ex(1),ex(2),ex(3)
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

            do l=1,nbatms
                bxl(l) = coxl + (rcox(l)*rt(1,1) + rcoy(l)*rt(1,2) + rcoz(l)*rt(1,3))
                byl(l) = coyl + (rcox(l)*rt(2,1) + rcoy(l)*rt(2,2) + rcoz(l)*rt(2,3))
                bzl(l) = cozl + (rcox(l)*rt(3,1) + rcoy(l)*rt(3,2) + rcoz(l)*rt(3,3))
            end do
        elseif ((enj-eval).le.1e-6) then
            flag=-1 
        end if 
        deallocate(ex)
        return
        end subroutine  

!!******************************************
!!> @details Creates guess points for MED and MESP topography.
        SUBROUTINE TRANSLATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,batwtl,batvdwl,flag,tssize,enj)
!!******************************************      
        USE DAMBUILD_T
        IMPLICIT NONE 
        INTEGER::nbatms,flag,icurr
        REAL(KREAL):: enj
        REAL(KREAL), DIMENSION(3,3,3) :: enjt
        REAL(KREAL), DIMENSION(nbatms):: btx,bty,btz
        REAL(KREAL), DIMENSION(nbatms):: rcom,rcox,rcoy,rcoz
        CHARACTER(2), DIMENSION(nbatms):: bsl
        REAL(KREAL), DIMENSION(nbatms):: bxl,byl,bzl,bql,batwtl,batvdwl
        REAL(KREAL), DIMENSION(3,3) :: rt
        REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
        REAL(KREAL) :: theta,phi,ang,ua,va,wa,eval,tssize
        REAL(KREAL) :: coxl,coyl,cozl
        INTEGER(KINT), ALLOCATABLE :: ex(:)
        INTEGER(KINT) :: i,j,k,l,a1,a2,a3,i1,j1,k1,as,dm,is,js,ierr,dum,np,nt
        LOGICAL,external:: accept
!       
             
        enjt = 0.0d0
        do i= -1,1
            do j=-1,1
                do k=-1,1
!                   if (i.eq.0.and.j.eq.0.and.k.eq.0) cycle
!                   write(11,*) nb
!                   write(11,*)
                    do l = 1,nbatms
                        btx(l)= bxl(l)+i*tssize 
                        bty(l)= byl(l)+j*tssize
                        btz(l)= bzl(l)+k*tssize
!                       write(11,*) bs(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
                    end do

                  !  do l=1,nb
                  !     bxl(l)=btx(l)
                  !     byl(l)=bty(l)
                  !     bzl(l)=btz(l)
                  !  end do 
                    if (accept(icurr,nbatms,bxl,byl,bzl,batvdwl).eqv..false.) then
                        enjt(i+2,j+2,k+2)= +10000.0d0
                        cycle
                    end if    
                    do l=1,nbatms
                        CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
                        & dyytot,dyztot,dzztot,btx(l),bty(l),btz(l))
                        enjt(i+2,j+2,k+2)=enjt(i+2,j+2,k+2)+vtot*bql(l)
                    end do
                end do   
             end do
        end do
        dm=size(shape(enjt))
        allocate (ex(dm),stat=ierr)
        ex=minloc(enjt)
        eval=minval(enjt)
        if (eval<enj .and. (enj-eval).gt.1e-6) then
            ! print*, "Ok in translation"
            Enj=eval
!           write(6,*) ex(1),ex(2),ex(3)
            write(6,*) "Best energy of all translatory moves", eval
            i=ex(1)-2
            j=ex(2)-2
            k=ex(3)-2
 
            do l = 1,nbatms
                bxl(l)= bxl(l)+i*tssize
                byl(l)= byl(l)+j*tssize
                bzl(l)= bzl(l)+k*tssize
            end do
        elseif ((enj-eval).le.1e-6) then
            flag=-1
        end if
        deallocate(ex)
        return 
        end subroutine

!!=================================
!!    Center of Mass
!!=================================
        SUBROUTINE COMASS(ncatms,cx,cy,cz,catwt,coxl,coyl,cozl)
        USE DAMBUILD_T
        IMPLICIT NONE
        INTEGER:: i,ncatms
		REAL(KREAL):: coxl,coyl,cozl,tmass
        REAL(KREAL), DIMENSION(ncatms):: cx,cy,cz,catwt
       
        coxl=0.0d0; coyl=0.0d0; cozl=0.0d0; tmass=0.0d0
        do i=1,ncatms
            tmass=tmass+catwt(i)
            coxl=coxl+catwt(i)*cx(i)
            coyl=coyl+catwt(i)*cy(i)
            cozl=cozl+catwt(i)*cz(i)
        end do
        coxl=coxl/tmass; coyl=coyl/tmass; cozl=cozl/tmass
        return
        END SUBROUTINE  
!    
!!==================================
!! Case correction of Atomic Symbol
!!==================================
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
!
!!====================================
!! Van der Waals overlap check
!!====================================
    logical function accept(icurr,nbatms,bxl,byl,bzl,bvdwl)
    use dambuild_t
    implicit none
    integer(kint):: iu,ju,ku,ierr,dum,nbatms,icurr
    real(kreal), dimension(nbatms):: bxl,byl,bzl,bvdwl
    real(kreal), dimension(nbatms):: rcom,rcox,rcoy,rcoz
    real(kreal), dimension(ncen):: DIST
    
    accept=.true.
    do ju=1,ncen
        do ku=1,nbatms
            dist(ju) = sqrt((rcen(1,ju)-bxl(ku))**2+(rcen(2,ju)-byl(ku))**2+(rcen(3,ju)-bzl(ku))**2)
            if ((atvdw(ju)+bvdwl(ku)).gt.dist(ju)) then
                accept=.false.
                goto 98 
            end if   
        end do
98  end do       
    do iu = 1, icurr-1
        do ju=1,molecules(iu)%natoms
            do ku=1,nbatms
                dist(ju) = sqrt((molecules(iu)%bx(ju)-bxl(ku))**2+&
                                (molecules(iu)%by(ju)-byl(ku))**2+&
                                (molecules(iu)%bz(ju)-bzl(ku))**2)
                if ((molecules(iu)%bvdw(ju)+bvdwl(ku)).gt.dist(ju)) then
                    accept=.false.
                    goto 99 
                end if   
            end do
99      end do
    end do

    return
    end 
!
    SUBROUTINE WRITECLUSTERGEOM(dfilename,nmols,ener)
    USE DAMINITIAL_T
    USE DAMBUILD_T
    integer:: nmols,totalatoms,ierr
    real(kreal):: ener
    character(len=*)::dfilename

    open(unit=8, file = trim(adjustl(dfilename)),iostat=ierr)

    if (ierr.ne.0) call error(1,"Error in creating new file with &
    initial geometry of the cluster.")
    totalatoms=0
    do i =1,nmols
        totalatoms=totalatoms+ molecules(i)%natoms
    enddo   
    totalatoms=totalatoms+ncen
    write(8,*) totalatoms
    write(8,*) "Energy= ",ener             
    do i = 1,ncen
        write(8,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
    end do
    do i=1,nmols
        do j = 1,molecules(i)%natoms
            write(8,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
            molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
        enddo
    enddo
    close(8)
    END SUBROUTINE    
!
!!******************************************
!!> @details Creates guess points for MED and MESP topography.
!!      SUBROUTINE TUNE 
!!******************************************      
!!      USE DAMINITIAL_T
!!      USE DAMBUILD_T
!!      IMPLICIT NONE 
!!      REAL(KREAL), ALLOCATABLE:: enjt(:),rcom(:),rcox(:),rcoy(:),rcoz(:)
!!      REAL(KREAL), ALLOCATABLE:: btx(:),bty(:),btz(:)
!!      REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
!!      REAL(KREAL):: eval
!!      REAL(KREAL), DIMENSION(2):: frvcom,frvxcom,frvycom,frvzcom
!!      INTEGER(KINT):: i,j,k,l,n,ierr,as,np,nt,dm,is,js,nopt
!!      CHARACTER::dummy*10
!!      INTEGER(KINT):: ierr,ig,newpt
!!!      REAL(KREAL):: xmin,ymin,zmin,xmax,ymax,zmax
!!      CHARACTER(2), ALLOCATABLE:: cnsym(:)
!!      REAL(KREAL), ALLOCATABLE:: Vt(:), drvxt(:), drvyt(:), drvzt(:)
!!    REAL(KREAL), ALLOCATABLE:: xnpt(:),ynpt(:),znpt(:), atwtn(:),atvdwn(:),qnpt(:),dist(:)
!!!      REAL(KREAL) :: X1,Y1,Z1,den,vtot,drvxtot,drvytot,drvztot,
!!      REAL(KREAL) :: dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot,enerji
!!!      REAL(KREAL) :: Sr,bx,by,bz,cdx,cdy,cdz
!!      REAL(KREAL):: cnox,cnoy,cnoz,tnmass
!!    REAL(KREAL),ALLOCATABLE:: frvxcom(:),frvycom(:),frvzcom(:)
!!    INTEGER(KINT):: i,j,k
!!      allocate (vt(newpt),drvxt(newpt),drvyt(newpt),drvzt(newpt),stat=ierr) 
!!      
!!    frvxcom(1)=1000.0
!!    frvycom(1)=1000.0
!!    frvzcom(1)=1000.0
!!    frvxcom(2)=0.0E+0
!!    frvycom(2)=0.0E+0
!!    frvzcom(2)=0.0E+0
!!
!!        allocate (frvxcom(10000),frvycom(10000),frvzcom(10000),stat=ierr)
!!        frvxcom=0.0E+0;frvycom=0.0E+0;frvzcom=0.0E+0;enerji=0.0E+0
!!        do i =1,10000
!!            do j=1,ncen
!!                do k=1,newpt
!!                    dist(j) = sqrt((rcen(1,j)-xnpt(k))**2+(rcen(2,j)-ynpt(k))**2+(rcen(3,j)-znpt(k))**2)
!!                    if ((atvdw(j)+atvdwn(k)).ge.dist(j)) go to 168
!!                end do
!!            end do
!!          do ig = 1,newpt
!!              CALL DAMPOT(vt(ig),drvxt(ig),drvyt(ig),drvzt(ig),&
!!              &           dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xnpt(ig),ynpt(ig),znpt(ig))
!!          end do
!!
!!          do ig = 1,newpt
!!              frvxcom(i)=frvxcom(i)+drvxt(ig)*qnpt(ig)
!!              frvycom(i)=frvycom(i)+drvyt(ig)*qnpt(ig)
!!              frvzcom(i)=frvzcom(i)+drvzt(ig)*qnpt(ig)
!!                enerji=enerji+vt(ig)*qnpt(ig)
!!          end do
!!            write(7,*)frvxcom(i),frvycom(i),frvzcom(i),enerji
!!            if (i.eq.1) then
!!                do ig = 1,newpt
!!                    xnpt(ig)=(xnpt(ig)+frvxcom(i)*100)
!!                    ynpt(ig)=(ynpt(ig)+frvycom(i)*100)
!!                    znpt(ig)=(znpt(ig)+frvzcom(i)*100)
!!                end do
!!            end if
!!            if (abs(frvxcom(i)).le.abs(frvxcom(i-1)).and.i.gt.1.and.abs(frvxcom(i)).gt.10e-9) then
!!                do ig = 1,newpt
!!                    xnpt(ig)=(xnpt(ig)+frvxcom(i)*100)
!!                end do
!!            end if    
!!            if (abs(frvycom(i)).le.abs(frvycom(i-1)).and.i.gt.1.and.abs(frvycom(i)).gt.10e-9) then
!!                do ig = 1,newpt
!!                    ynpt(ig)=(ynpt(ig)+frvycom(i)*100)
!!                end do
!!            end if    
!!            if (abs(frvzcom(i)).le.abs(frvzcom(i-1)).and.i.gt.1.and.abs(frvzcom(i)).gt.10e-9) then
!!                do ig = 1,newpt
!!                    znpt(ig)=(znpt(ig)+frvzcom(i)*100)
!!                end do
!!            end if    
!!        end do
!!168     continue     
!!
!!  do while (abs(frvxcom(2)).le.abs(frvxcom(1)).or. &
!!         &  abs(frvycom(2)).le.abs(frvycom(1)).or. &
!!         &  abs(frvzcom(2)).le.abs(frvzcom(1)))
!!     do l=1,nb
!!          CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
!!          & dyytot,dyztot,dzztot,bx(l),by(l),bz(l))
!!          frvxcom(2)=frvxcom(2)+drvx*bq(l)
!!          frvycom(2)=frvycom(2)+drvy*bq(l)
!!          frvzcom(2)=frvzcom(2)+drvz*bq(l)
!!          enj=enj+vtot*bq(l)
!!     end do
!!! frvcom(2)=sqrt(frvxcom(2)**2+frvycom(2)**2+frvzcom(2)**2)
!!
!!     if (.not.allocated(btx)) allocate (btx(nb),stat=ierr)
!!     if (.not.allocated(bty)) allocate (bty(nb),stat=ierr)
!!     if (.not.allocated(btz)) allocate (btz(nb),stat=ierr)
!!           if (abs(frvxcom(2)).le.abs(frvxcom(1)).and.abs(frvxcom(2)).gt.10e-9) then
!!              do i = 1,nb
!!                   btx(i)=(bx(i)+frvxcom(2)*100)
!!                end do
!!           end if    
!!           if (abs(frvycom(2)).le.abs(frvycom(1)).and.abs(frvycom(2)).gt.10e-9) then
!!                  do i = 1,nb
!!                   bty(i)=(by(i)+frvycom(2)*100)
!!                end do
!!           end if    
!!           if (abs(frvzcom(2)).le.abs(frvzcom(1)).and.abs(frvzcom(2)).gt.10e-9) then
!!              do i = 1,nb
!!                   btz(i)=(bz(i)+frvzcom(2)*100)
!!                end do
!!           end if    
!!     write(11,*) frvxcom(2),frvycom(2),frvzcom(2),frvcom(2),enj
!!
!!     frvxcom(1)=frvxcom(2)
!!     frvycom(1)=frvycom(2)
!!     frvzcom(1)=frvzcom(2)
!!
!!           do j=1,ncen
!!              do k=1,nb
!!                  dist(j) = sqrt((rcen(1,j)-btx(k))**2+(rcen(2,j)-bty(k))**2+(rcen(3,j)-btz(k))**2)
!!                  if ((atvdwn(j)+atvdwn(k)).ge.dist(j)) return
!!              end do
!!           end do
!!
!!     do i=1,nb
!!        bx(i)=btx(i)
!!        by(i)=bty(i)
!!          bz(i)=btz(i)
!!     end do
!!
!!    end do
!!      return 
!!      end subroutine
!!!
