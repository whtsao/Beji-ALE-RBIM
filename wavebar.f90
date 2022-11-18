!===============================================================================
!	SOLVE 2D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY ALE-RBIM WITH RK4
!	WEN-HAUI TSAO 2022 LSU
!===============================================================================
      PROGRAM WATER_WAVE
      IMPLICIT NONE
      INTEGER SMTYP,NNRL,MDEG,IR,IL
	  INTEGER I,J,NT,QRULE,NGA,NARC,NTIM,OUTSTEP,NFIELD,NPL,WTYP,NNODE,ICON,NITER,NWG,E1LOC,OUTYP,ALETYP,NSTYP
	  INTEGER,ALLOCATABLE::ND(:),NS(:),BUNTYP(:),IPIV(:)
	  REAL*8 DEP,CLEN,D_OUT,GRAV,MU,ARDZONE,WIDTH,THO,DELTTIME,WGX(20),WGY(20)
	  REAL*8 TIME,DIS,VEL,ACC,P_ATM,FOR,DDIS,DTEMP,WC,E1,E2,ETOL,EB
      REAL*8 WAVE_PERIOD,WAVE_HEIGHT,PSI,AMP,OMEGA,WAVE_K,WAVE_C,temp
	  REAL*8,ALLOCATABLE::COOR(:,:),SIDE_L(:)
	  REAL*8,ALLOCATABLE::NODE(:,:),WNODE(:),NORM(:,:),JCB(:),PHI(:),PPHI(:),PHIT(:),PPHIT(:)
	  REAL*8,ALLOCATABLE::KER1(:,:),KER2(:,:),DP(:,:),DPDT(:),D2P(:,:),D2PDT(:),ACCMO(:),PR(:)
      REAL*8,ALLOCATABLE::DPDS(:),DPDSS(:),DPNDS(:),DPTDS(:)
	  REAL*8,ALLOCATABLE::PHIT_TEMP(:)
	  REAL*8,ALLOCATABLE::AWT(:),ART(:) !WT(:),RT(:),SHA1(:),SHA2(:),SH(:,:),
      REAL*8,ALLOCATABLE::G1(:,:),H1(:,:),EYE(:)
	  REAL*8,ALLOCATABLE::DI(:),VE(:),AC(:)
      REAL*8,ALLOCATABLE::K1(:,:),K2(:,:),K3(:,:),K4(:,:)

      OPEN(UNIT=1,FILE='1.ipt',STATUS='OLD')
      OPEN(UNIT=2,FILE='2.ipt',STATUS='OLD')
      OPEN(UNIT=3,FILE='3.ipt',STATUS='OLD')
      OPEN(UNIT=5,FILE='IO.DAT')
      OPEN(UNIT=6,FILE='S.DAT')
      OPEN(UNIT=7,FILE='WG.DAT')
!      OPEN(UNIT=8,FILE='P.DAT')
!      OPEN(UNIT=9,FILE='F.DAT')
!      OPEN(UNIT=10,FILE='E.DAT')
!      OPEN(UNIT=11,FILE='DOMAIN.DAT')
	  OPEN(UNIT=21,FILE='ERR.DAT')
	  OPEN(UNIT=22,FILE='ABORT.TXT')
	  OPEN(UNIT=23,FILE='CFL.DAT')
      OPEN(UNIT=99,FILE='TEST.TXT')

!---SMOOTHING PARAMETERS
    CALL INPUT_3(SMTYP,NNRL,MDEG,NARC,ALETYP,NSTYP)
    ALLOCATE(AWT(NARC),ART(NARC))
    AWT=0.D0
    ART=0.D0
    CALL GAUSS(AWT,ART,NARC) ! this is the gaussian for calculating arc length of cubic spline of an element
    
!---TOPOGRAGHY AND WAVE TYPE
	CALL INPUT_2(NPL,WTYP,OUTYP,NWG,WGX,WAVE_PERIOD,WAVE_HEIGHT,PSI)
	ALLOCATE(ND(NPL),NS(NPL),BUNTYP(NPL),COOR(NPL,2),SIDE_L(NPL)) !NELEM(NPL),ME(NPL),
    ! NELEM and ME are useless in RBIM because no element is needed
	!NELEM=0
    ND=0
	NS=0
	BUNTYP=0
	COOR=0.D0
	SIDE_L=0.D0

!---INPUT ALL KINDS OF PARAMETERS
	CALL INPUT_1(QRULE,NPL,COOR,NFIELD,NNODE,ND,NS,BUNTYP,DEP,CLEN,NGA,GRAV,MU,ARDZONE,WIDTH,THO,&
			  &NTIM,DELTTIME,OUTSTEP,ICON,NITER,ETOL,EB)
    ! no more needed
    !ALLOCATE(LN(NELM,2),LENG(NELM))
	!ALLOCATE(WT(NGA),RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA))
    
    ! size has to change, NELN -> NNODE
    ALLOCATE(NORM(NNODE,2),JCB(NNODE))
    
    ! stay here
	ALLOCATE(NODE(NNODE,2),WNODE(NNODE),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),PHIT_TEMP(NNODE))
	ALLOCATE(KER1(NNODE,NNODE),KER2(NNODE,NNODE),DP(NNODE,2))
	ALLOCATE(DPDS(NNODE),DPDSS(NNODE),DPNDS(NNODE),DPTDS(NNODE),DPDT(NNODE),PR(NNODE),ACCMO(NNODE),D2P(NNODE,2),D2PDT(NNODE))
	ALLOCATE(DI(NTIM),VE(NTIM),AC(NTIM))
    ALLOCATE(G1(NNODE,NNODE),H1(NNODE,NNODE),EYE(NNODE),IPIV(NNODE))
    ALLOCATE(K1(NNODE,3),K2(NNODE,3),K3(NNODE,3),K4(NNODE,3))
!	LN=0
	NODE=0.D0
    WNODE=0.D0
	NORM=0.D0
	JCB=0.D0
!	LENG=0.D0
	PHI=0.D0
	PPHI=0.D0
	PHIT=0.D0
	PPHIT=0.D0
	KER1=0.D0
	KER2=0.D0
	DP=0.D0
	DPDS=0.D0
    DPDSS=0.D0
    DPNDS=0.D0
    DPTDS=0.D0
	DPDT=0.D0
	PR=0.D0
	ACCMO=0.D0
    D2P=0.D0
    D2PDT=0.D0
	PHIT_TEMP=0.D0
!	WT=0.D0
!	RT=0.D0
!	SHA1=0.D0
!	SHA2=0.D0
!	SH=0.D0
	DI=0.D0
	VE=0.D0
	AC=0.D0
    G1=0.D0
    H1=0.D0
    EYE=1.D0
    IPIV=0
    K1=0.D0
    K2=0.D0
    K3=0.D0
    K4=0.D0
    
!---PREPARE IF OUTLET IS A WALL (IF A WALL, NO ITERATION NEEDED)
    IF (OUTYP==0)THEN
        NITER=1
    END IF
        
!---SHAPE FUNCTION AND MESH
!    CALL GAUSS(WT,RT,NGA)
!    CALL SHAP(SHA1,SHA2,SH,NGA,RT)
	CALL LENGTH(NPL,COOR,SIDE_L)
	CALL MESH(QRULE,NPL,NNODE,ND,COOR,SIDE_L,NODE,WNODE)
    WRITE(*,*) 'PASS MESH'

!---GENERATE WAVES
	SELECT CASE (WTYP)
    CASE(1)
        CALL PERIODIC(DEP,WAVE_PERIOD,WAVE_HEIGHT,AMP)
        OMEGA=2.D0*DACOS(-1.D0)/WAVE_PERIOD
		DO I=1,NTIM
		 TIME=(I-1)*DELTTIME
		 DI(I)=AMP*DCOS(OMEGA*TIME+PSI)-AMP !AMP*SIN(OMEGA*TIME+PSI)
		 VE(I)=-AMP*OMEGA*DSIN(OMEGA*TIME+PSI) !AMP*OMEGA*COS(OMEGA*TIME+PSI)
		 AC(I)=-AMP*OMEGA**2*DCOS(OMEGA*TIME+PSI) !-AMP*OMEGA**2*SIN(OMEGA*TIME+PSI)
		 END DO
	CASE(2)
		NT=OMEGA/DELTTIME
		DO I=1,NT+1
		TIME=(I-1)*DELTTIME
		DI(I)=0.5D0*AMP-0.5D0*AMP*DCOS(DACOS(-1.D0)/OMEGA*TIME) !STR/DUR*TIME
		VE(I)=DACOS(-1.D0)/OMEGA*0.5D0*AMP*DSIN(DACOS(-1.D0)/OMEGA*TIME) !STR/DUR
		AC(I)=(DACOS(-1.D0)/OMEGA)**2*0.5D0*AMP*DCOS(DACOS(-1.D0)/OMEGA*TIME) !0.D0
		END DO
		DI(NT+2:NTIM)=DI(NT+1)
    CASE(3)
        NT=10.D0/DELTTIME
        CALL SOLITARY(DEP,WAVE_HEIGHT,GRAV,1,1,NT,DELTTIME,DI,VE,AC)
        DI(NT+1:NTIM)=DI(NT)
        VE(NT+1:NTIM)=VE(NT)
        AC(NT+1:NTIM)=AC(NT)
    END SELECT      
    
DO NT=1,NTIM-1
	WRITE(*,*) NT,'TH'
     TIME=(NT-1)*DELTTIME
!---WAVEMAKER DIS, VEL, ACC, PHASE LAG (IN RADIUS)
	 DIS=DI(NT)
     VEL=VE(NT)
     ACC=AC(NT)
	 WRITE(5,"(7(E15.8,1X))") TIME,DIS,VEL,ACC

!---CALCULATE WAVE SPEED FOR RADIATION CONDITION
	D_OUT=NODE(NS(1),2)-NODE(NS(2),2) !COOR(NPL-2,2)
	CALL WAVE_SPD(GRAV,OMEGA,D_OUT,WC)

!**************IMPLICIT SOLVER FOR PHIT AND PPHIT ON THE ABSORPTION SIDE**************
DO J=1,NITER

    PHIT_TEMP=PHIT
    !---SOLVE PHIT
    !CALL SOLVE_PHIT(NARC,NPL,NNODE,ND,NS,NSTYP,BUNTYP,OUTYP,EB,WC,VEL,ACC,DEP,GRAV,MU,CLEN,ARDZONE,&
    !               &EYE,AWT,ART,NODE,WNODE,PHI,PHIT)

    !---SOLVE BY RK4 METHOD
    CALL RK4_KI(DELTTIME,EB,WC,VEL,VE(NT+1),DEP,GRAV,MU,CLEN,ARDZONE,&
               &NARC,NPL,NNODE,ND,NS,NSTYP,BUNTYP,OUTYP,QRULE,EYE,NODE,PHI,WNODE,AWT,ART,K1,K2,K3,K4)
    WRITE(*,*) 'PASS RK4 KI'

    CALL RK4_SUM(NPL,NNODE,ND,NS,QRULE,DELTTIME,K1,K2,K3,K4,NODE,PHI)
    WRITE(*,*) 'PASS RK4 SUM'

    !---CHECK IF RADIATION CONDITION CONVERGES, IF INCLUDED
    IF (OUTYP==1)THEN
      !---SOLVE FOR PHIT TO SEE IF IT IS CONVERGED
      !CALL SOLVE_PHIT()
      CALL CONVERGE(ND(2),PHIT_TEMP(NS(1)+1:NS(2)),PHIT(NS(1)+1:NS(2)),E1LOC,E1,E2)
        WRITE(*,*) J,E1,E2
	    IF(E1<=ETOL.AND.E2<ETOL)THEN
		    WRITE(*,*) 'PASS CONVERGED',J
		    WRITE(21,*) TIME,J,E1,E2
		    GOTO 205
	    ELSE IF(J>=NITER)THEN
		    WRITE(22,*) 'CG FAIL',TIME,E1LOC,E1,E2
		    WRITE(*,*) 'CONVERGE FAIL'
		    STOP
	    END IF
    END IF
END DO
!**************************************************************************************

205 CONTINUE

!---COMPUTE PRESSURE AND OTHERS
    !++++NEED TO CALCULATE THE PHIT ON THE WALL, CAN BE ADDED IN SUBROUTINE PRESSURE
    !CALL PRESSURE(ICON,TIME,THO,GRAV,DEP,NPL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
    
!---SMOOTHING FOR FREE-SURFACE NODE AND PHI
      IF (SMTYP==1)THEN
          CALL FS_SMOOTH(NNRL,MDEG,NS(1),NODE(1:NS(1),2))
          CALL FS_SMOOTH(NNRL,MDEG,NS(1),PHI(1:NS(1)))
          WRITE(*,*) 'PASS SMOOTH'
      END IF

!-----WRITE BOUNDARY-----
      IF(MOD(NT,OUTSTEP)==0)THEN
          WRITE(6,"(5000(E15.8,1X))") NODE(:,1) 
          WRITE(6,"(5000(E15.8,1X))") NODE(:,2)
      END IF      
      
!-----WAVE GAUGES (IN CM)-----
      DO I=1,NWG
          CALL BWLOC(-WGX(I),NS(1),-NODE(1:NS(1),1),0,IR,IL)
          TEMP=(WGX(I)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))
          WGY(I)=NODE(IL,2)+TEMP*(NODE(IR,2)-NODE(IL,2))-DEP
      END DO
      WRITE(7,"(11(E15.8,1X))") TIME,WGY(1:NWG)*100.D0      
      
      
END DO



STOP
END

!**********************************************************************
      SUBROUTINE RK4_KI(DT,EB,WC,VEL0,VELN,DEP,GRAV,MU,CLEN,ARDZONE,&
                       &NARC,NPL,NNODE,ND,NS,NSTYP,BUNTYP,OUTYP,QRULE,EYE,NODE1,PHI1,WNODE1,AWT,ART,&
                       &K1,K2,K3,K4)
!**********************************************************************
      IMPLICIT NONE
      !---INPUT
      INTEGER NARC,NPL,NNODE,ND(NPL),NS(NPL),NSTYP,BUNTYP(NPL),OUTYP,QRULE
      REAL*8 DT,DT_MID,EB,WC,VEL0,VELN,DEP,GRAV,MU,CLEN,ARDZONE
      REAL*8 EYE(NNODE),AWT(NARC),ART(NARC)
      REAL*8 NODE1(NNODE,2),WNODE1(NNODE),PHI1(NNODE),PPHI1(NNODE)
      !---OUTPUT
      REAL*8 K1(NNODE,3),K2(NNODE,3),K3(NNODE,3),K4(NNODE,3)
      !---LOCAL
      INTEGER I
      REAL*8 VEL_MID
      REAL*8 JCB(NNODE),NORM(NNODE,2),DP(NNODE,2),DPDS(NNODE),DPDSS(NNODE),DPDT(NNODE),PHIT(NNODE),CV(NNODE,2)
      REAL*8 KER1(NNODE,NNODE),KER2(NNODE,NNODE)
      REAL*8 NODE2(NNODE,2),WNODE2(NNODE),PHI2(NNODE),PPHI2(NNODE)
      REAL*8 NODE3(NNODE,2),WNODE3(NNODE),PHI3(NNODE),PPHI3(NNODE)
      REAL*8 NODE4(NNODE,2),WNODE4(NNODE),PHI4(NNODE),PPHI4(NNODE)

      VEL_MID=0.5D0*(VEL0+VELN)
      DT_MID=0.5D0*DT
      NODE2=NODE1
      NODE3=NODE1
      NODE4=NODE1
      WNODE2=WNODE1
      WNODE3=WNODE1
      WNODE4=WNODE1

      !---GET K1
      CALL KERNEL_RBIM(KER1,KER2,NODE1,NORM,WNODE1,JCB,NPL,NNODE,ND,NS,EYE,EB,NSTYP) ! USE 1
      CALL BOUND(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHI1,PPHI1,PHIT,VEL0,WC) ! USE 1
      CALL SOLVE_LAPACK(NPL,PHI1,PPHI1,KER1,KER2,NNODE,ND,BUNTYP) ! USE 1
      CALL TANGENT_DERIV(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE1,PHI1,DPDS,DPDSS) ! USE 1
      CALL FIRST_ORDER_TERMS(NARC,NPL,ND,NS,NNODE,AWT,ART,DEP,GRAV,MU,CLEN,ARDZONE,NODE1,NORM,&
                            &PHI1,PPHI1,DPDS,PHIT,DP,DPDT) ! USE 1
      CALL ALE_VEL(NPL,NNODE,NS,NORM,DP,DPDT,K1(:,1:2),K1(:,3),CV) ! USE 1
      CALL ALE_BFD(DT_MID,NPL,NNODE,NS,NODE1,PHI1,K1(:,1:2),K1(:,3),NODE2,PHI2) ! GET 2
      CALL REMESH(QRULE,NPL,ND,NS,NNODE,NODE2,WNODE2) ! GET 2
    
      !---GET K2
      CALL KERNEL_RBIM(KER1,KER2,NODE2,NORM,WNODE2,JCB,NPL,NNODE,ND,NS,EYE,EB,NSTYP)
      CALL BOUND(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHI2,PPHI2,PHIT,VEL_MID,WC)
      CALL SOLVE_LAPACK(NPL,PHI2,PPHI2,KER1,KER2,NNODE,ND,BUNTYP)
      CALL TANGENT_DERIV(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE2,PHI2,DPDS,DPDSS)
      CALL FIRST_ORDER_TERMS(NARC,NPL,ND,NS,NNODE,AWT,ART,DEP,GRAV,MU,CLEN,ARDZONE,NODE2,NORM,&
                            &PHI2,PPHI2,DPDS,PHIT,DP,DPDT)
      CALL ALE_VEL(NPL,NNODE,NS,NORM,DP,DPDT,K2(:,1:2),K2(:,3),CV)
      CALL ALE_BFD(DT_MID,NPL,NNODE,NS,NODE1,PHI1,K2(:,1:2),K2(:,3),NODE3,PHI3)
      CALL REMESH(QRULE,NPL,ND,NS,NNODE,NODE3,WNODE3)
      
      !---GET K3
      CALL KERNEL_RBIM(KER1,KER2,NODE3,NORM,WNODE3,JCB,NPL,NNODE,ND,NS,EYE,EB,NSTYP)
      CALL BOUND(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHI3,PPHI3,PHIT,VEL_MID,WC)
      CALL SOLVE_LAPACK(NPL,PHI3,PPHI3,KER1,KER2,NNODE,ND,BUNTYP)
      CALL TANGENT_DERIV(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE3,PHI3,DPDS,DPDSS)
      CALL FIRST_ORDER_TERMS(NARC,NPL,ND,NS,NNODE,AWT,ART,DEP,GRAV,MU,CLEN,ARDZONE,NODE3,NORM,&
                            &PHI3,PPHI3,DPDS,PHIT,DP,DPDT)
      CALL ALE_VEL(NPL,NNODE,NS,NORM,DP,DPDT,K3(:,1:2),K3(:,3),CV)
      CALL ALE_BFD(DT,NPL,NNODE,NS,NODE1,PHI1,K3(:,1:2),K3(:,3),NODE4,PHI4)
      CALL REMESH(QRULE,NPL,ND,NS,NNODE,NODE4,WNODE4)
      
      !---GET K4
      CALL KERNEL_RBIM(KER1,KER2,NODE4,NORM,WNODE4,JCB,NPL,NNODE,ND,NS,EYE,EB,NSTYP)
      CALL BOUND(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHI4,PPHI4,PHIT,VELN,WC)
      CALL SOLVE_LAPACK(NPL,PHI4,PPHI4,KER1,KER2,NNODE,ND,BUNTYP)
      CALL TANGENT_DERIV(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE4,PHI4,DPDS,DPDSS)
      CALL FIRST_ORDER_TERMS(NARC,NPL,ND,NS,NNODE,AWT,ART,DEP,GRAV,MU,CLEN,ARDZONE,NODE4,NORM,&
                            &PHI4,PPHI4,DPDS,PHIT,DP,DPDT)
      CALL ALE_VEL(NPL,NNODE,NS,NORM,DP,DPDT,K4(:,1:2),K4(:,3),CV)
      
      RETURN
    END
!**********************************************************************
      SUBROUTINE RK4_SUM(NPL,NNODE,ND,NS,QRULE,DELTTIME,K1,K2,K3,K4,NODE,PHI)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,NPL,NNODE,ND(NPL),NS(NPL),QRULE
      REAL*8 DELTTIME
      REAL*8 K1(NNODE,3),K2(NNODE,3),K3(NNODE,3),K4(NNODE,3)
      REAL*8 NODE(NNODE,2),PHI(NNODE),WNODE(NNODE)
      
      !---UPDATE THE FREE SURFACE WITH RK4
      DO I=1,NS(1)
          NODE(I,1)=NODE(I,1)+DELTTIME*(K1(I,1)+2.D0*K2(I,1)+2.D0*K3(I,1)+K4(I,1))/6.D0
          NODE(I,2)=NODE(I,2)+DELTTIME*(K1(I,2)+2.D0*K2(I,2)+2.D0*K3(I,2)+K4(I,2))/6.D0
          PHI(I)=PHI(I)+DELTTIME*(K1(I,3)+2.D0*K2(I,3)+2.D0*K3(I,3)+K4(I,3))/6.D0
      END DO
      
      !---REMESH FOR THE NEXT TIME STEP
      CALL REMESH(QRULE,NPL,ND,NS,NNODE,NODE,WNODE) 
      
      RETURN
    END
!**********************************************************************
    SUBROUTINE SOLVE_PHIT(NARC,NPL,NNODE,ND,NS,NSTYP,BUNTYP,OUTYP,EB,WC,VEL,ACC,DEP,GRAV,MU,CLEN,ARDZONE,&
                        &EYE,AWT,ART,NODE,WNODE,PHI,PHIT)
!**********************************************************************
    IMPLICIT NONE
    INTEGER NARC,NPL,NNODE,ND(NPL),NS(NPL),NSTYP,BUNTYP(NPL),OUTYP
    REAL*8 EB,WC,VEL,ACC,DEP,GRAV,MU,CLEN,ARDZONE
    REAL*8 EYE(NNODE),AWT(NARC),ART(NARC)
    REAL*8 NODE(NNODE,2),WNODE(NNODE),JCB(NNODE),NORM(NNODE,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
    REAL*8 DP(NNODE,2),DPDT(NNODE),DPDS(NNODE),DPDSS(NNODE),DPNDS(NNODE),DPNDSS(NNODE),ACCMO(NNODE)
    REAL*8 KER1(NNODE,NNODE),KER2(NNODE,NNODE)

    !---COMPUTE VELOCITY---
    CALL KERNEL_RBIM(KER1,KER2,NODE,NORM,WNODE,JCB,NPL,NNODE,ND,NS,EYE,EB,NSTYP)
    CALL BOUND(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHIT,PPHIT,PHIT,VEL,WC)
    CALL SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,ND,BUNTYP)
    CALL TANGENT_DERIV(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE,PHI,DPDS,DPDSS)
    CALL TANGENT_DERIV_N(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE,PPHI,DPNDS,DPNDSS)
    CALL FIRST_ORDER_TERMS(NARC,NPL,ND,NS,NNODE,AWT,ART,DEP,GRAV,MU,CLEN,ARDZONE,NODE,NORM,PHI,PPHI,DPDS,PHIT,DP,DPDT)
    
    !---COMPUTE PHIT
    CALL ACCBC(NNODE,NORM,DP,DPDSS,DPNDS,ACCMO)
    CALL BOUNDT(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHIT,PPHIT,DPDS,DPDSS,WC,ACC,ACCMO)
    CALL SOLVE_LAPACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,ND,BUNTYP)
      
      RETURN
    END
!**********************************************************************
      SUBROUTINE GET_NORM_JCB(NPL,NNODE,ND,NS,NODE,NORM,JCB)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,NPL,NNODE,ND(NPL),NS(NPL)
      REAL*8 NODE(NNODE,2),NORM(NNODE,2),JCB(NNODE)
      REAL*8 B(ND(1)),C(ND(1)),D(ND(1)),Y0,Y1,Y2,DX,DY,DS

!---CALCULATE THE SLOPE OF FREE SURFACE
      CALL SPLINE(ND(1),NODE(1:ND(1),1),NODE(1:ND(1),2),B,C,D)
      DO I=1,ND(1)
          CALL SEVAL(ND(1),NODE(I,1),NODE(1:ND(1),1),NODE(1:ND(1),2),B,C,D,Y0,Y1,Y2)
          JCB(I)=DSQRT(1.D0+Y1**2)
          NORM(I,1)=-Y1/JCB(I)
          NORM(I,2)=1.D0/JCB(I)
      END DO

!---CALCULATE THE SLOPE OF FIXED BOUNDARY
      DO I=2,NPL
          DX=NODE(NS(I),1)-NODE(NS(I-1)+1,1)
          DY=NODE(NS(I),2)-NODE(NS(I-1)+1,2)
          DS=DSQRT(DX*DX+DY*DY)
          DO J=NS(I-1)+1,NS(I)
              JCB(J)=1.D0
              NORM(J,1)=-DY/DS
              NORM(J,2)=DX/DS              
          END DO
      END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE HEADLINE(ID,IREAD)
!**********************************************************************
      CHARACTER*2 ID
      ID =  '*'
      DO WHILE (ID .EQ. '*')
      READ(IREAD,'(A1)') ID
      END DO
      RETURN
      END
!**********************************************************************
SUBROUTINE INPUT_1(QRULE,NPL,COOR,NFIELD,NNODE,ND,NS,BUNTYP,DEP,CLEN,NGA,GRAV,MU,ARDZONE,WIDTH,THO,&
                  &NTIM,DELTTIME,OUTSTEP,ICON,NITER,ETOL,EB)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER I,J,QRULE,NPL,NFIELD,NNODE,NTIM,ICON,NITER,NGA,OUTSTEP !,NELM
	  INTEGER ND(NPL),NS(NPL),BUNTYP(NPL) !NELEM(NPL),ME(NPL),
	  REAL*8 COOR(NPL,2),DEP,CLEN,ENDTIME,DELTTIME,GRAV,MU,ARDZONE,WIDTH,THO,ETOL,EB
      CHARACTER*2 ID
         ID = '*'
!---NODAL COORDINATES
         CALL HEADLINE(ID,1)
         READ(1,*)  ((COOR(I,J),J=1,2),I=1,NPL)
!---ELEMENT NODE NUMBER
        CALL HEADLINE(ID,1)
        READ(1,*) (ND(I),I=1,NPL)
        NNODE = 0
        !NELM  = 0
		!ME(1) = NELEM(1)
        DO I=1,NPL
         !NELM = NELM+NELEM(I)
         NNODE = NNODE+ND(I)
        END DO

		!DO I=2,NPL
		! ME(I)=ME(I-1)+NELEM(I)
		!END DO
        NS(1)=ND(1)
        DO I=2,NPL
            NS(I)=NS(I-1)+ND(I)
        END DO

	  NFIELD=(ND(1)-2)*(ND(NPL)-2)

!---QUADRATURE RULE
         CALL HEADLINE(ID,1)
        READ(1,*) QRULE
      
!---BOUNDARY TYPE
         CALL HEADLINE(ID,1)
        READ(1,*) (BUNTYP(I),I=1,NPL)
!---READ THE GAUSSIAN INTEGRATION POINT NO.
         CALL HEADLINE(ID,1)
         READ(1,*)  NGA
!---READ GRAV ACC, MAX OF WAVE DAMPING COEF, ZONE OF WAVE DAMPING, TANK WIDTH, FLUID DENSITY
         CALL HEADLINE(ID,1)
         READ(1,*)  GRAV,MU,ARDZONE,WIDTH,THO
!---READ TIME,TIME INTERVAL,OUTPUT STEP,INDEX OF BERNOULLI CONSTANT, NUMBER IF ITERATION,TOLERANCE ERROR, DISTANCE OF NEAR SINGULARITY
         CALL HEADLINE(ID,1)
         READ(1,*)  ENDTIME,DELTTIME,OUTSTEP,ICON,NITER,ETOL,EB
		NTIM=ENDTIME/DELTTIME+1
		DEP=MAXVAL(COOR(:,2))
        CLEN=MAXVAL(COOR(:,1))

      RETURN
      END
!**********************************************************************
      SUBROUTINE INPUT_2(NPL,WTYP,OUTYP,NWG,WGX,WAVE_PERIOD,WAVE_HEIGHT,PSI)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER NPL,WTYP,OUTYP,NWG
	  REAL*8 WAVE_PERIOD,WAVE_HEIGHT,PSI,WGX(10)
      CHARACTER*2 ID
         ID = '*'
!---NUMBER OF PLANES
         CALL HEADLINE(ID,2)
         READ(2,*) NPL
!---WAVE GENERATION: 1=PERIODIC; 2=SOLITARY
         CALL HEADLINE(ID,2)
         READ(2,*) WTYP
!---READ WAVE_PERIOD,WAVE_HEIGHT
         CALL HEADLINE(ID,2)
         READ(2,*)  WAVE_PERIOD,WAVE_HEIGHT,PSI
!---WAVE GENERATION: 0=WALL; 1=RADIATION
         CALL HEADLINE(ID,2)
         READ(2,*) OUTYP
!---WAVE GAUGE NUMBER
         CALL HEADLINE(ID,2)
         READ(2,*) NWG
		 IF (NWG>10)THEN
		 WRITE(*,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 WRITE(22,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 STOP
		 END IF
!---WAVE GAUGE LOCATION
         CALL HEADLINE(ID,2)
         READ(2,*) WGX(1:NWG)

      RETURN
    END
!**********************************************************************
      SUBROUTINE INPUT_3(SMTYP,NNRL,MDEG,NARC,ALETYP,NSTYP)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER SMTYP,NNRL,MDEG,NARC,ALETYP,NSTYP
      CHARACTER*2 ID
         ID = '*'
!---Do you need free-surface smoothing: 1 = yes; 0 = no
         CALL HEADLINE(ID,3)
         READ(3,*) SMTYP
!---number of neighboring node and degree of polynomial for SG filter
!---note that nl = nr and nl + nr < m
         CALL HEADLINE(ID,3)
         READ(3,*) NNRL,MDEG

!---number of gaussian quadrature for calculating arc length
         CALL HEADLINE(ID,3)
         READ(3,*) NARC
         
!---do you use ALE approach? 0 = NO; 1 = YES
         CALL HEADLINE(ID,3)
         READ(3,*) ALETYP
         
!---do you remove near singularity? 0 = NO; 1 = YES
         CALL HEADLINE(ID,3)
         READ(3,*) NSTYP
      RETURN
    END
!**********************************************************************
      SUBROUTINE LENGTH(NPL,COOR1,SIDE_L1)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
	  INTEGER NPL
      REAL*8  COOR1(NPL,2),SIDE_L1(NPL)
      DO I=1,NPL-1
        SIDE_L1(I)=DSQRT((COOR1(I+1,1)-COOR1(I,1))**2+(COOR1(I+1,2)-COOR1(I,2))**2)
      END DO
        SIDE_L1(NPL)=DSQRT((COOR1(NPL,1)-COOR1(1,1))**2+(COOR1(NPL,2)-COOR1(1,2))**2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(QRULE,NPL,NNODE,ND,COOR,LENG,NODE,WNODE)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER QRULE,NPL,ND(NPL),NNODE !,NELM,NELEM(NPL),LN(NELM,2)
      REAL*8 SX,SY,NORM,DELT,LENG(NPL),COOR(NPL,2),NODE(NNODE,2),WNODE(NNODE),S(5000),W(5000)

S=0.D0
W=0.D0
K=0
DO I=1,NPL-1
	J=NPL-I
    !NEL(I) = NELEM(I)+1
    !DELT=LENG(J)/NELEM(I)
	SX=COOR(J,1)-COOR(J+1,1)
	SY=COOR(J,2)-COOR(J+1,2)
	NORM=DSQRT(SX**2+SY**2)
    IF (QRULE==0)THEN
        CALL TRAP(ND(I),S,W,0.D0,NORM)
    ELSE IF (QRULE==1)THEN
        CALL SIMP(ND(I),S,W,0.D0,NORM)
    ELSE IF (QRULE==2)THEN
        CALL LOBATTO(ND(I),S,W,0.D0,NORM)
    END IF
	!SX=SX/NORM
	!SY=SY/NORM
    DO L=1,ND(I)
       NODE(L+K,1)=COOR(J+1,1)+S(L)*SX/NORM
       NODE(L+K,2)=COOR(J+1,2)+S(L)*SY/NORM
       WNODE(L+K)=W(L)
    END DO
    K=K+ND(I)
END DO

    !NEL(NPL) = NELEM(NPL)+1
    !DELT=LENG(NPL)/NELEM(NPL)
	SX=COOR(NPL,1)-COOR(1,1)
	SY=COOR(NPL,2)-COOR(1,2)
	NORM=DSQRT(SX**2+SY**2)
    IF (QRULE==0)THEN
        CALL TRAP(ND(NPL),S,W,0.D0,NORM)
    ELSE IF (QRULE==1)THEN
        CALL SIMP(ND(NPL),S,W,0.D0,NORM)
    ELSE IF (QRULE==2)THEN
        CALL LOBATTO(ND(NPL),S,W,0.D0,NORM)
    END IF
    !SX=SX/NORM
	!SY=SY/NORM
    DO I=1,ND(NPL)
      NODE(I+K,1)=COOR(1,1)+S(I)*SX/NORM
      NODE(I+K,2)=COOR(1,2)+S(I)*SY/NORM
      WNODE(I+K)=W(I)
    END DO

!----TO GIVE THE LOCAL ELEMENT NODE NUMBER (USELESS IN RBIM)
!      L=1
!	  N=1
!      DO I=1,NPL
!       DO J=1,NELEM(I)
!        LN(N,1)=L
!        LN(N,2)=L+1
!        L=L+1
!		N=N+1
!       END DO
!       L=L+1
!      END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE SHAP(SHA1,SHA2,SH,NGA,RT)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NGA
      REAL*8 RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
      DO M=1,NGA
        SHA1(M)=0.5D0*(1-RT(M))
        SHA2(M)=0.5D0*(1+RT(M))

        SH(1,M)=SHA1(M)
        SH(2,M)=SHA2(M)
      END DO
	RETURN
	END 
!**********************************************************************
      SUBROUTINE REMESH(QRULE,NPL,ND,NS,NNODE,NODE,WNODE)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,K,L,QRULE,NPL,NNODE,ND(NPL),NS(NPL)
      REAL*8 SX,SY,DS,TEMP
      REAL*8 NODE(NNODE,2),WNODE(NNODE),S(5000),W(5000)
      REAL*8 B(ND(1)),C(ND(1)),D(ND(1)),Y0,Y1,Y2

!----- REMESH FOR THE FREE SURFACE
    CALL SPLINE(ND(1),NODE(1:ND(1),1),NODE(1:ND(1),2),B,C,D)
    IF (QRULE==0)THEN
        CALL TRAP(ND(1),S,W,NODE(1,1),NODE(ND(1),1))
    ELSE IF (QRULE==1)THEN
        CALL SIMP(ND(1),S,W,NODE(1,1),NODE(ND(1),1))
    ELSE IF (QRULE==2)THEN
        CALL LOBATTO(ND(1),S,W,NODE(1,1),NODE(ND(1),1))
    END IF
    
    DO I=1,ND(1)
        CALL SEVAL(ND(1),S(I),NODE(1:ND(1),1),NODE(1:ND(1),2),B,C,D,Y0,Y1,Y2)
        NODE(I,1)=S(I)
        NODE(I,2)=Y0
        WNODE(I)=W(I)
    END DO

!------ENSURE DUPLICATE POINT ON END NODE OF FREE SURFACE-----
    NODE(NNODE,:)=NODE(1,:)
	NODE(NS(1)+1,:)=NODE(NS(1),:)

!------BOTTOM END NODE GOES WITH FREE SURFACE REMAIN VERTICAL WALL-----
	NODE(NS(NPL-1)+1,1)=NODE(1,1)
	NODE(NS(2),1)=NODE(NS(1),1)
    
!------ENSURE DUPLICATE POINT ON END NODE OF VERTICAL WALLS...seems unnecessary
!    NODE(NS(NPL-1),:)=NODE(NS(NPL-1)+1,:)
!	NODE(NS(2)+1,:)=NODE(NS(2),:)
    
!----- REMESH FOR ALL WETTED BOUNDARIES
K=NS(1)
DO I=2,NPL
	SX=NODE(NS(I),1)-NODE(NS(I-1),1)
	SY=NODE(NS(I),2)-NODE(NS(I-1),2)
	DS=DSQRT(SX**2+SY**2)
    
    IF (QRULE==0)THEN
        CALL TRAP(ND(I),S,W,0.D0,DS)
    ELSE IF (QRULE==1)THEN
        CALL SIMP(ND(I),S,W,0.D0,DS)
    ELSE IF (QRULE==2)THEN
        CALL LOBATTO(ND(I),S,W,0.D0,DS)
    END IF
    DO L=1,ND(I)
       NODE(L+K,1)=NODE(NS(I-1),1)+S(L)*SX/DS
       NODE(L+K,2)=NODE(NS(I-1),2)+S(L)*SY/DS
       WNODE(L+K)=W(L)
    END DO
    K=K+ND(I)
END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE BOUND(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHI,PPHI,PHIT,VEL,WC)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8 (A-H,O-Z)
       INTEGER NPL,NNODE,ND(NPL),NS(NPL),BUNTYP(NPL),OUTYP
	   REAL*8 R,VEL,WC
       REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			   PHI(J)=PHI(J)
			   PPHI(J)=0.D0
            ELSE
			  IF (I==2)THEN
                  IF(OUTYP==0)THEN
                      PPHI(J)=0.D0
                      PHI(J)=0.D0
                  ELSE
                      PPHI(J)=-PHIT(J)/WC
                      PHI(J)=0.D0
                  END IF
			  ELSE IF (I==NPL)THEN
			   PPHI(J)=-VEL
			   PHI(J)=0.D0
			  ELSE
			   PPHI(J)=0.D0
			   PHI(J)=0.D0
			  END IF
            END IF
          END DO
          N=N+ND(I)
       END DO

      RETURN
    END
!**********************************************************************
SUBROUTINE KERNEL_RBIM(KER1,KER2,NODE,NORM,WNODE,JCB,NPL,NNODE,ND,NS,EYE,EB,NSTYP)
!**********************************************************************
IMPLICIT NONE
INTEGER I,J,K,NPL,NNODE,ND(NPL),NS(NPL),NSTYP
REAL*8 L(NPL),XP,EB,DX,DY,RD,SIGMA(NNODE),EYE(NNODE)
REAL*8 KER1(NNODE,NNODE),KER2(NNODE,NNODE),KER3(NNODE),KER4(NNODE)
REAL*8 NODE(NNODE,2),WNODE(NNODE),JCB(NNODE),NORM(NNODE,2)

!***CALCULATE THE NORMAL VECTOR AND JACOBIAN
CALL GET_NORM_JCB(NPL,NNODE,ND,NS,NODE,NORM,JCB)

KER1=0.D0
KER2=0.D0
KER3=0.D0
KER4=0.D0

!*****SUBTRACTION AND ADDITION TECHNIQUE FOR DE-SINGULARIZATION*****
L(1)=NODE(NS(1),1)-NODE(1,1)
DO I=1,ND(1)
    !---subtracion term for free surface
	DO J=1,ND(1)
	    IF(I/=J)THEN
	    KER3(I)=KER3(I)+DLOG(ABS(NODE(I,1)-NODE(J,1)))*WNODE(J)
	    END IF
    END DO
    !---addition term for free surface
    IF(I.EQ.1.OR.I.EQ.ND(1))THEN
        KER4(I)=L(1)*DLOG(L(1))-L(1)
    ELSE
        XP=NODE(I,1)-NODE(1,1)        
        KER4(I)=(L(1)-XP)*DLOG(L(1)-XP)+XP*LOG(XP)-L(1)
    END IF
END DO

DO K=2,NPL
    L(K)=DSQRT((NODE(NS(K-1)+1,1)-NODE(NS(K),1))**2+(NODE(NS(K-1)+1,2)-NODE(NS(K),2))**2)
    DO I=NS(K-1)+1,NS(K)
        !---subtracion term for all wetted walls (all straight lines)
        DO J=NS(K-1)+1,NS(K)
	        IF(I/=J)THEN
                RD=DSQRT((NODE(J,1)-NODE(I,1))**2+(NODE(J,2)-NODE(I,2))**2)
                KER3(I)=KER3(I)+DLOG(RD)*WNODE(J)
            END IF
        END DO
        !---addition term for all wetted walls (all straight lines)
        IF(I.EQ.NS(K-1)+1.OR.I.EQ.NS(K))THEN
            KER4(I)=L(K)*DLOG(L(K))-L(K)
        ELSE
            XP=DSQRT((NODE(I,1)-NODE(NS(K-1)+1,1))**2+(NODE(I,2)-NODE(NS(K-1)+1,2))**2)            
            KER4(I)=(L(K)-XP)*DLOG(L(K)-XP)+XP*LOG(XP)-L(K)
        END IF
    END DO
END DO

!*****THE SURFACE KERNELS*****
DO I=1,NNODE
	DO J=1,NNODE
        DX=NODE(J,1)-NODE(I,1)
        DY=NODE(J,2)-NODE(I,2)
        RD=DSQRT(DX*DX+DY*DY)
		IF (RD<=0.000001D0) THEN
		    KER1(I,J)=0.D0
		    KER2(I,J)=DLOG(JCB(J))*JCB(J)*WNODE(J)-KER3(J)*JCB(J)+KER4(J)*JCB(J)
		ELSE
		    KER1(I,J)=(DX*NORM(J,1)+DY*NORM(J,2))*JCB(J)*WNODE(J)/RD**2
		    KER2(I,J)=DLOG(RD)*JCB(J)*WNODE(J)
		END IF
	END DO
END DO

!*****SUBTRACTION AND ADDITION TECHNIQUE FOR DE-NEAR-SINGULARIZATION*****
    IF (NSTYP==1)THEN
        CALL NEAR(EB,NPL,NNODE,ND,NS,NODE,WNODE,NORM,L,KER1,KER2)
    END IF

!*****DIAGONAL TERMS OF KER1*****
    DO I=1,NNODE
        KER1(I,I)=0.D0
    END DO
    CALL DGEMM('N','N',NNODE,1,NNODE,1.D0,KER1,NNODE,EYE,NNODE,0.D0,SIGMA,NNODE)
    DO I=1,NNODE
        KER1(I,I)=-SIGMA(I)
    END DO  

RETURN
END
!**********************************************************************
SUBROUTINE NEAR(EB,NPL,NNODE,ND,NS,NODE,WNODE,NORM,LENG,KER1,KER2)
!**********************************************************************
IMPLICIT NONE
INTEGER I,J,C,K,L,D,M,NPL,NNODE,ND(NPL),NS(NPL),ID
REAL*8 EB,LENG(NPL),SGN
REAL*8 R1(2),RD1,R2(2),RD2,DX,DY,RD3,DST,L1,L2,EPS,TAU,BETA1,BETA2,SUB1,SUB2,ADD1,ADD2
REAL*8 NODE(NNODE,2),WNODE(NNODE),NORM(NNODE,2),KER1(NNODE,NNODE),KER2(NNODE,NNODE)

C=0
DO I=1,NPL
    DO J=C+1,NS(I) ! J is the index of P on surface I
        D=0
        DO K=1,NPL

            IF (K.NE.I)THEN ! if q'=D+1 or NS(K) is on surface K =! surface I
                ! search for the closest node to p, i.e., q', and its distance
                CALL search_closest(NODE(J,:),ND(K),NODE(D+1:NS(K),:),D+1,ID,DST) ! get q'=ID
                IF (DST>0.00001D0.AND.DST<=EB)THEN! now we have to deal with near singularity                
                    R1(1)=NODE(D+1,1)-NODE(J,1)
                    R1(2)=NODE(D+1,2)-NODE(J,2)
                    RD1=DSQRT(R1(1)*R1(1)+R1(2)*R1(2))
                    R2(1)=NODE(NS(K),1)-NODE(J,1)
                    R2(2)=NODE(NS(K),2)-NODE(J,2)
                    RD2=DSQRT(R2(1)*R2(1)+R2(2)*R2(2))
                    DX=NODE(NS(K),1)-NODE(D+1,1)
                    DY=NODE(NS(K),2)-NODE(D+1,2)
                    RD3=DSQRT(DX*DX+DY*DY)                
                    DX=NODE(ID,1)-NODE(J,1)
                    DY=NODE(ID,2)-NODE(J,2)
                    SGN=DX*NORM(ID,1)+DY*NORM(ID,2)
                    EPS=DABS(SGN)
                    
                    ! calculate the subtend angle by surface K at p
                    IF(EPS>=0.00001D0)THEN
                        ADD1=SIGN(DACOS(0.5D0*(RD1**2+RD2**2-RD3**2)/RD1/RD2),SGN)
                    ELSE
                        ADD1=0.D0
                    END IF
                    
                    L1=R1(1)*NORM(ID,2)-R1(2)*NORM(ID,1)
                    L2=R2(1)*NORM(ID,2)-R2(2)*NORM(ID,1)
                    SGN=SIGN(1.D0,L1*L2)
                    L1=DABS(L1)
                    L2=DABS(L2)
                    BETA1=DATAN(L1/EPS)
                    BETA2=DATAN(L2/EPS)
                    
                    ! calculate the analytical line integral of a source
                    IF (SGN>0.D0.AND.L2>L1)THEN
                        ADD2=0.5D0*L2*DLOG(L2**2+EPS**2)-0.5D0*L1*DLOG(L1**2+EPS**2)-L2+L1+EPS*(BETA2-BETA1)
                    ELSE IF (SGN>0.D0.AND.L1>L2)THEN
                        ADD2=0.5D0*L1*DLOG(L1**2+EPS**2)-0.5D0*L2*DLOG(L2**2+EPS**2)-L1+L2+EPS*(BETA1-BETA2)                    
                    ELSE IF (SGN<0.D0)THEN
                        ADD2=0.5D0*L1*DLOG(L1**2+EPS**2)+0.5D0*L2*DLOG(L2**2+EPS**2)-(L1+L2)+EPS*(BETA1+BETA2)
                    END IF
                    
                    SUB1=0.D0
                    SUB2=0.D0
                    DO L=D+1,NS(K)
                        SUB1=SUB1+KER1(J,L)
                        DX=NODE(L,1)-NODE(J,1)
                        DY=NODE(L,2)-NODE(J,2)
                        RD3=DABS(DX*NORM(ID,2)-DY*NORM(ID,1))
                        TAU=DSQRT(EPS**2+RD3**2)
                        SUB2=SUB2+DLOG(TAU)*WNODE(L)
                    END DO
                    KER1(J,ID)=KER1(J,ID)+ADD1-SUB1
                    KER2(J,ID)=KER2(J,ID)+ADD2-SUB2
                END IF
            END IF
        D=D+ND(K)
        END DO
    END DO
    C=C+ND(I)
END DO


RETURN
END
!**********************************************************************
       SUBROUTINE SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!      USING GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION
!======================================================================
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
       REAL*8    PHI(NNODE),PPHI(NNODE)
       REAL*8    KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8    H1(NNODE,NNODE),Q1(NNODE),TEMP(NNODE)
       REAL*8  SUM,A,SIG,G1(NNODE,NNODE),P1(NNODE)

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
!-------------------
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       DO I=1,NNODE
          TEMP(I)=0.D0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO
!*************GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION*********
      DO I=1,NNODE
        G1(I,NNODE+1)=TEMP(I)
      END DO
      DO I=1,NNODE
         SUM=0.D0
         DO K=I,NNODE
            IF (G1(K,I) .NE. 0) THEN
               IF (K .NE. I) THEN
               IF (G1(I,I) .EQ. 0.D0) THEN
                 WRITE(*,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 WRITE(22,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 STOP
               END IF
               A=G1(K,I)/G1(I,I)
               DO J=I,NNODE+1
                  G1(K,J)=G1(K,J)-A*G1(I,J)
               END DO
               END IF
            END IF
            SUM=SUM+G1(K,I)
          END DO
          IF (SUM .EQ. 0.D0) THEN
          WRITE(*,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          WRITE(22,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          STOP
          END IF
      END DO

      P1(NNODE)=G1(NNODE,NNODE+1)/G1(NNODE,NNODE)
      DO I=NNODE-1,1,-1
         SIG=0.D0
         DO J=I+1,NNODE
            SIG=G1(I,J)*P1(J)+SIG
          END DO
         P1(I)=(G1(I,NNODE+1)-SIG)/G1(I,I)
      END DO
!=================================
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,ND,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!======================================================================
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,ND(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)

!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+ND(I)
       END DO
       
!-----MOVE KER1 AND KER2---- 
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,ND(L)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+ND(L)
          END DO
        END DO

       P1=MATMUL(H1,Q1)

!*************SOLVE BY CALLING LAPACK*********
CALL DGESV(NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+ND(I)
       END DO
       
      RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK2_1(NPL,PHI,PPHI,KER1,KER2,G1,H1,NNODE,ND,BUNTYP,IPIV)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,ND(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+ND(I)
       END DO
       
!-----MOVE KER1 AND KER2---- 
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,ND(L)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+ND(L)
          END DO
        END DO

       P1=MATMUL(H1,Q1)
       
!*************SOLVE BY CALLING LAPACK*********
CALL DGETRF(NNODE,NNODE,G1,NNODE,IPIV,INFO)
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)
       
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+ND(I)
       END DO
                     
      RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK2_2(NPL,PHI,PPHI,G1,H1,NNODE,ND,BUNTYP,IPIV)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,ND(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+ND(I)
       END DO

       P1=MATMUL(H1,Q1)

!*************SOLVE BY CALLING LAPACK*********
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+ND(I)
       END DO

      RETURN
    END
!********************************************************************
SUBROUTINE TANGENT_DERIV(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE,PHI,DPDS,DPDSS)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NARC,NPL,NNODE,ND(NPL),NS(NPL)
	REAL*8 AWT(NARC),ART(NARC)
	REAL*8 NODE(NNODE,2),PHI(NNODE),DPDS(NNODE),DPDSS(NNODE)
             
	DPDS=0.D0
!*********************ALL PLANES EXCEPT INLET AND OUTLET*********************
    K=0
    DO I=1,NPL
      IF (I.NE.2 .AND. I.NE.NPL)THEN ! avoid vertical wall
        CALL FS_CSDIFF(ND(I),NODE(K+1:NS(I),1),NODE(K+1:NS(I),2),PHI(K+1:NS(I)),AWT,ART,NARC,DPDS(K+1:NS(I)),DPDSS(K+1:NS(I)))
      END IF
      K=K+ND(I)
    END DO

!*********************ON RIGHT WALL*********************
    CALL FS_CSDIFF(ND(2),NODE(NS(1)+1:NS(2),2),NODE(NS(1)+1:NS(2),1),PHI(NS(1)+1:NS(2)),AWT,ART,NARC,DPDS(NS(1)+1:NS(2)),DPDSS(NS(1)+1:NS(2)))
    
!*********************ON LEFT WALL*********************
    CALL FS_CSDIFF(ND(NPL),NODE(NS(NPL-1)+1:NS(NPL),2),NODE(NS(NPL-1)+1:NS(NPL),1),PHI(NS(NPL-1)+1:NS(NPL)),AWT,ART,NARC,DPDS(NS(NPL-1)+1:NS(NPL)),DPDSS(NS(NPL-1)+1:NS(NPL)))
    DPDS(NS(NPL-1)+1)=0.D0

      RETURN
    END
!********************************************************************
SUBROUTINE TANGENT_DERIV_N(NARC,NPL,NNODE,ND,NS,AWT,ART,NODE,PPHI,DPNDS,DPNDSS)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NARC,NPL,NNODE,ND(NPL),NS(NPL)
	REAL*8 AWT(NARC),ART(NARC)
	REAL*8 NODE(NNODE,2),PPHI(NNODE),DPNDS(NNODE),DPNDSS(NNODE)
             
	DPNDS=0.D0
!*********************ALL PLANES EXCEPT INLET AND OUTLET*********************
    K=0
    DO I=1,NPL
      IF (I.NE.2 .AND. I.NE.NPL)THEN ! avoid vertical wall
        CALL FS_CSDIFF(ND(I),NODE(K+1:NS(I),1),NODE(K+1:NS(I),2),PPHI(K+1:NS(I)),AWT,ART,NARC,DPNDS(K+1:NS(I)),DPNDSS(K+1:NS(I)))
      END IF
      K=K+ND(I)
    END DO

!*********************ON RIGHT WALL*********************
    CALL FS_CSDIFF(ND(2),NODE(NS(1)+1:NS(2),2),NODE(NS(1)+1:NS(2),1),PPHI(NS(1)+1:NS(2)),AWT,ART,NARC,DPNDS(NS(1)+1:NS(2)),DPNDSS(NS(1)+1:NS(2)))

!*********************ON LEFT WALL*********************
    CALL FS_CSDIFF(ND(NPL),NODE(NS(NPL-1)+1:NS(NPL),2),NODE(NS(NPL-1)+1:NS(NPL),1),PPHI(NS(NPL-1)+1:NS(NPL)),AWT,ART,NARC,DPNDS(NS(NPL-1)+1:NS(NPL)),DPNDSS(NS(NPL-1)+1:NS(NPL)))

      RETURN
    END
!********************************************************************
SUBROUTINE FIRST_ORDER_TERMS(NARC,NPL,ND,NS,NNODE,AWT,ART,DEP,GRAV,MU,CLEN,ARDZONE,NODE,NORM,PHI,PPHI,DPDS,PHIT,DP,DPDT)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NARC,NPL,NNODE,ND(NPL),NS(NPL)
    REAL*8 DEP,GRAV,MU,CLEN,ARDZONE
	REAL*8 NODE(NNODE,2),NORM(NNODE,2),PHI(NNODE),PPHI(NNODE),DP(NNODE,2)
	REAL*8 DPDS(NNODE),PHIT(NNODE),DPDT(NNODE),AWT(NARC),ART(NARC)

!*********************ON FREE SURFACE*********************
    DO I=1,NS(1)
        IF(I.EQ.1) THEN
            DP(I,1)=-PPHI(NNODE)
            DPDS(I)=(DP(I,1)-PPHI(I)*NORM(I,1))/NORM(I,2)
            DP(I,2)=-DPDS(I)*NORM(I,1)+PPHI(I)*NORM(I,2)
		!ELSE IF(I.EQ.NS(1)) THEN
		!    DP(I,1)=PPHI(I+1)
        !    DPDS(I)=(DP(I,1)-PPHI(I)*NORM(I,1))/NORM(I,2)
        !    DP(I,2)=-DPDS(I)*NORM(I,1)+PPHI(I)*NORM(I,2)
        ELSE
			DP(I,1)=DPDS(I)*NORM(I,2)+PPHI(I)*NORM(I,1)
			DP(I,2)=-DPDS(I)*NORM(I,1)+PPHI(I)*NORM(I,2)
		END IF
    END DO
        
    DO I=1,NS(1)
        IF (NODE(I,1)>=ARDZONE)THEN
          DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*(NODE(I,1)-ARDZONE)/(CLEN-ARDZONE)*PHI(I)
        ELSE
          DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)
        END IF
        PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
    END DO

!*********************ON LEFT AND RIGHT WALL AND BOTTOM*********************
    DO I=NS(1)+1,NNODE
        DP(I,1)=DPDS(I)*NORM(I,2)+PPHI(I)*NORM(I,1)
        DP(I,2)=-DPDS(I)*NORM(I,1)+PPHI(I)*NORM(I,2)
    END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE FS_CSDIFF(N,X,Y,PHI,AWT,ART,NARC,DPDS,DPDSS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,N,NE,NARC
      REAL*8 X(N),Y(N),PHI(N),DPDS(N),DPDSS(N),AWT(NARC),ART(NARC),S(N),DARC,U
      REAL*8 B(N),C(N),D(N),Y0,Y1,Y2

!---CALCULATE THE ARC LENGTH AND TANGENTIAL COORDINATE
      CALL SPLINE(N,X,Y,B,C,D)
      S(1)=0.D0
      DO I=2,N
          DARC=0.D0
          DO J=1,NARC
              U=0.5D0*(1-ART(J))*X(I-1)+0.5D0*(1+ART(J))*X(I)
              CALL SEVAL(N,U,X,Y,B,C,D,Y0,Y1,Y2)
              DARC=DARC+DSQRT(1.D0+Y1**2)*0.5D0*(X(I)-X(I-1))*AWT(J)
          END DO
          S(I)=S(I-1)+DARC
      END DO

!---calculate dPHI/dS
      CALL SPLINE(N,S,PHI,B,C,D)
      DO I=1,N
          CALL SEVAL(N,S(I),S,PHI,B,C,D,Y0,DPDS(I),DPDSS(I))
      END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC(NNODE,NORM,DP,DPDSS,DPNDS,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,NNODE
    REAL*8 DPNDX,DPNDY,DPDNN
	REAL*8 NORM(NNODE,2),DP(NNODE,2),DPDSS(NNODE),DPNDS(NNODE),ACCMO(NNODE)
    
    DO I=1,NNODE
      DPDNN=-DPDSS(I)
      DPNDX=DPDNN*NORM(I,1)+DPNDS(I)*NORM(I,2)
      DPNDY=DPDNN*NORM(I,2)-DPNDS(I)*NORM(I,1)
      ACCMO(I)=DP(I,1)*DPNDX+DP(I,2)*DPNDY
    END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUNDT(NPL,NNODE,ND,NS,BUNTYP,OUTYP,PHIT,PPHIT,DPDS,DPDSS,WC,ACC,ACCMO)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8 (A-H,O-Z)
       INTEGER NNODE,ND(NPL),NS(NPL),BUNTYP(NPL),OUTYP
	   REAL*8 ACC,WC
       REAL*8 PHIT(NNODE),PPHIT(NNODE),DPDS(NNODE),DPDSS(NNODE),ACCMO(NNODE)

       PPHIT=0.D0
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,ND(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			 PHIT(J)=PHIT(J)
			 PPHIT(J)=0.D0
            ELSE
			  IF (I==2)THEN
                  IF (OUTYP==0)THEN
                      PPHIT(J)=0.D0
                      PHIT(J)=0.D0 
                  ELSE
                      PPHIT(J)=WC*DPDSS(J)
                      PHIT(J)=0.D0
                  END IF
			  ELSE IF (I==NPL)THEN
			    PPHIT(J)=-ACC-ACCMO(J)
			    PHIT(J)=0.D0
			  ELSE
			    PPHIT(J)=0.D0-ACCMO(J)
			    PHIT(J)=0.D0
			  END IF
            END IF
          END DO
          N=N+ND(I)
       END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE ALE_VEL(NPL,NNODE,NS,NORM,DP,DPDT,DPM,DPDTM,CV)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,NPL,NNODE,NS(NPL)
      REAL*8 NORM(NNODE,2),DP(NNODE,2),DPDT(NNODE),DPM(NNODE,2),DPDTM(NNODE),CV(NNODE,2)
      
!---CALCULATE MESH VELOCITY      
      DO I=1,NS(1)
          DPM(I,1)=0.D0
          DPM(I,2)=DP(I,2)+(DP(I,1)-DPM(I,1))*NORM(I,1)/NORM(I,2)
      END DO
      
!---CALCULATE CONVECTIVE VELOCITY AND DPDT IN ALE FRAME
      DO I=1,NS(1)
          CV(I,1)=DP(I,1)-DPM(I,1)
          CV(I,2)=DP(I,2)-DPM(I,2)
          DPDTM(I)=DPDT(I)-CV(I,1)*DP(I,1)-CV(I,2)*DP(I,2)
      END DO
      
      RETURN
    END
!**********************************************************************
SUBROUTINE ALE_BFD(DELTTIME,NPL,NNODE,NS,NODE,PHI,DPM,DPDTM,NEW_NODE,NEW_PHI)
!**********************************************************************
      IMPLICIT NONE 
      INTEGER I,NPL,NNODE,NS(NPL)
      REAL*8 DELTTIME
      REAL*8 NODE(NNODE,2),PHI(NNODE),NEW_NODE(NNODE,2),NEW_PHI(NNODE)
	  REAL*8 DPM(NNODE,2),DPDTM(NNODE) !,D2PM(NS(1),2),D2PDTM(NS(1))

      DO I=1,NS(1)
		NEW_PHI(I)=PHI(I)+DELTTIME*DPDTM(I) !+0.5D0*DELTTIME**2*D2PDTM(I)
		NEW_NODE(I,1)=NODE(I,1)+DELTTIME*DPM(I,1) !+0.5D0*DELTTIME**2*D2PM(I,1)
		NEW_NODE(I,2)=NODE(I,2)+DELTTIME*DPM(I,2) !+0.5D0*DELTTIME**2*D2PM(I,2)
      END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE FS_SMOOTH(NNRL,MDEG,N,Y)
!**********************************************************************
      IMPLICIT NONE
      INTEGER N,NNRL,MDEG,FLAG
      REAL*8 Y(N)
      
      CALL savgol_filter(NNRL,NNRL,0,MDEG,N,Y,flag)
      
      RETURN
    END
!********************************************************************
SUBROUTINE PRESSURE(ICON,TIME,THO,GRAV,DEP,NPL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
!********************************************************************
    IMPLICIT NONE
    INTEGER  I,ICON,NPL,NNODE,NS(NPL)
    REAL*8 TIME,DEP,THO,GRAV,P1,P2,P3,P_ATM
    REAL*8 NODE(NNODE,2),PHIT(NNODE),DP(NNODE,2),PR(NNODE)
    REAL*8 CP1(NS(1)),CP2(NS(1)),CP3(NS(1)),CP(NS(1))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON ! NS(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP)
	CP(I)=CP1(I)+CP2(I)+CP3(I)
	ENDDO
	P_ATM=CP(ICON)

!----PRESSURE ON ZONE 1 BOUNDARY
	DO I=NS(1)+2,NNODE-1
	P1=THO*PHIT(I)
	P2=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP)
	PR(I)=P_ATM-(P1+P2+P3)
	END DO

      RETURN
      END
!********************************************************************
SUBROUTINE DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
				 &SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,M,IL,IR,NPL,NGA,NFIELD,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
	  REAL*8 THO,GRAV,DEP,P_ATM
	  REAL*8 HB,DX,DY,PI2,TEMP,P1,P2,P3
	  REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
	  REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 DNODE(NFIELD,2),DVX(NFIELD),DVY(NFIELD),DPHIT(NFIELD),DPR(NFIELD)
	  REAL*8 KER1(NFIELD,NNODE),KER2(NFIELD,NNODE)
      REAL*8 H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
	  REAL*8 WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
	
	PI2=2.D0*DACOS(-1.D0)

!----CREATE DOMAIN POINT
	L=1
	DO I=2,NS(1)-1
	  CALL BWLOC(NODE(I,1),NS(NPL-1)-NS(2),NODE(NS(2)+1:NS(NPL-1),1),NS(2),IL,IR)
	  HB=NODE(IL,2)+(NODE(I,1)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))*(NODE(IR,2)-NODE(IL,2))+0.01D0 ! keep it a little far away from the boundary
	  DY=-NODE(I,2)/NELEM(NPL)
		DO J=2,NELEM(NPL)
		  TEMP=NODE(I,2)+DY*(J-1)
			IF(TEMP>HB)THEN
			DNODE(L,1)=NODE(I,1)
			DNODE(L,2)=TEMP !NODE(I,2)+DY*(J-1)
			L=L+1
			END IF
		END DO
	END DO

!--- SET A DUMMY NODE TO USELESS DNODE
	DO I=L,NFIELD
	DNODE(I,:)=DNODE(1,:)
	END DO

	KER1=0.D0
	KER2=0.D0
	DVX=0.D0
!----CALCULATE X VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((YFUNC(M)-DNODE(I,2))**2-(XFUNC(M)-DNODE(I,1))**2)*NORM(J,1)-&
					&2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,2))/&
					&((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(XFUNC(M)-DNODE(I,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVX=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2

	KER1=0.D0
	KER2=0.D0
	DVY=0.D0
!----CALCULATE Y VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((XFUNC(M)-DNODE(I,1))**2-(YFUNC(M)-DNODE(I,2))**2)*NORM(J,2)-&
					&2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,1))/&
					&((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(YFUNC(M)-DNODE(I,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVY=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2
	
	KER1=0.D0
	KER2=0.D0
	DPHIT=0.D0
!----CALCULATE PARTIAL POTENTIAL OVER TIME BY BIE
      DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(-1.D0)/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*&
					&((XFUNC(M)-DNODE(I,1))*NORM(J,1)+(YFUNC(M)-DNODE(I,2))*NORM(J,2))*TEMP
            G(K)=G(K)+DLOG(1.D0/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**0.5D0)*TEMP	 
		 END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DPHIT=(MATMUL(KER2,PPHIT)-MATMUL(KER1,PHIT))/PI2

!----CALCULATE PRESSURE DISTRIBUTION IN DOMAIN
	DO I=1,NFIELD
	P1=THO*DPHIT(I)
	P2=THO*0.5D0*(DVX(I)**2+DVY(I)**2)
	P3=THO*GRAV*(DNODE(I,2)-DEP)
	DPR(I)=P_ATM-(P1+P2+P3)
	END DO

	WRITE(11,'(6000(1X,F15.7))') NODE(:,1),DNODE(:,1)
	WRITE(11,'(6000(1X,F15.7))') NODE(:,2),DNODE(:,2)
	WRITE(11,'(6000(1X,F15.7))') DP(:,1),DVX
	WRITE(11,'(6000(1X,F15.7))') DP(:,2),DVY
	WRITE(11,'(6000(1X,F15.7))') PR,DPR

      RETURN
      END
!********************************************************************
      SUBROUTINE GAUSS(WT,RT,NGA)
!********************************************************************
      INTEGER NGA
      REAL*8   WT(NGA),RT(NGA)

      SELECT CASE(NGA)
       CASE(3)
        WT(1)=0.55555555
        WT(2)=0.88888889
        WT(3)=0.55555555
        RT(1)=0.77459667
        RT(2)=0.D0
        RT(3)=-0.77459667
       CASE(4)
        WT(1)=0.65214515
        WT(2)=0.34785484
        WT(3)=0.34785484
        WT(4)=0.65214515
        RT(1)=0.33998104
        RT(2)=0.86113631
        RT(3)=-0.86113631
        RT(4)=-0.33998104
       CASE(5)
        WT(1)=0.23692689
        WT(2)=0.47862867
        WT(3)=0.56888889
        WT(4)=0.47862867
        WT(5)=0.23692689
        RT(1)=0.90617985
        RT(2)=0.53846931
        RT(3)=0.D0
        RT(4)=-0.53846931
        RT(5)=-0.90617985
	 CASE(6)
	  WT(1)=0.17132449
	  WT(2)=0.36076157
	  WT(3)=0.46791393
	  WT(4)=0.46791393
	  WT(5)=0.36076157
	  WT(6)=0.17132449
	  RT(1)=0.93246951
	  RT(2)=0.66120938
	  RT(3)=0.23861918
	  RT(4)=-0.23861918
	  RT(5)=-0.66120938
	  RT(6)=-0.9346951
       CASE(8)
        WT(1)=0.1012285362903763D0
        WT(2)=0.2223810344533745D0
        WT(3)=0.3137066458778873D0
        WT(4)=0.3626837833783620D0
        WT(8)=0.1012285362903763D0
        WT(7)=0.2223810344533745D0
        WT(6)=0.3137066458778873D0
        WT(5)=0.3626837833783620D0
        RT(1)=0.9602898564975363D0
        RT(2)=0.7966664774136267D0
        RT(3)=0.5255324099163290D0
        RT(4)=0.1834346424956498D0
        RT(8)=-0.9602898564975363D0
        RT(7)=-0.7966664774136267D0
        RT(6)=-0.5255324099163290D0
        RT(5)=-0.1834346424956498D0
       CASE(10)
        WT(1)=0.D06667134
        WT(2)=0.14945134
        WT(3)=0.21908636
        WT(4)=0.26926671
        WT(5)=0.29552422
        WT(10)=0.D06667134
        WT(9)=0.14945134
        WT(8)=0.21908636
        WT(7)=0.26926671
        WT(6)=0.29552422
        RT(1)=0.97390652
        RT(2)=0.86506336
        RT(3)=0.67940956
        RT(4)=0.43339539
        RT(5)=0.14887433
        RT(10)=-0.97390652
        RT(9)=-0.86506336
        RT(8)=-0.67940956
        RT(7)=-0.43339539
        RT(6)=-0.14887433
      END SELECT

      RETURN
      END
!********************************************************************
	SUBROUTINE CONVERGE(N,P1,P2,E1LOC,E1,E2)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,K(1),E1LOC
	REAL*8 P1(N),P2(N),E1,E2

	E1=MAXVAL(DABS(P1-P2))
	E1LOC=MAXLOC(DABS(P1-P2),1)

	E2=0.D0
	DO I=1,N
	E2=E2+(P1(I)-P2(I))**2
	END DO
	E2=DSQRT(E2/N)

	RETURN
	END
!********************************************************************
	SUBROUTINE BWLOC(PX,N,X,IST,IL,IR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,IST,IL,IR
	REAL*8 PX,X(N)
    
	DO I=1,N-1
		IF(X(I)>PX.AND.X(I+1)<=PX)THEN
		IR=IST+I
		IL=IST+I+1
		GOTO 777
		END IF
	END DO
777 CONTINUE

	RETURN
	END
!********************************************************************
	SUBROUTINE WAVE_SPD(GRAV,OMEGA,D,C)
!********************************************************************
	IMPLICIT NONE
	INTEGER I
	REAL*8 K,K2,GRAV,OMEGA,D,C,PI,F0,F1
	PI=DACOS(-1.D0)

	K=1.D0
	DO I=1,100
	F0=K*DTANH(K*D)-OMEGA**2/GRAV
	F1=DTANH(K*D)+K-K*(DTANH(K*D)**2)
	K2=K-(F0)/(F1)
		IF((K2-K)/K<=0.000001D0) THEN
		GOTO 717
		END IF
	K=K2
	END DO
	717 CONTINUE

	C=DSQRT(GRAV*DTANH(K*D)/K)

	RETURN
	END
!********************************************************************
SUBROUTINE COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,CFLOC,NNODE,NELM,LN(NELM,2)
    REAL*8 TIME,DELTTIME,U,V,VE
	REAL*8 NODE(NNODE,2),DP(NNODE,2),JCB(NELM),CN(NELM),CFL

  DO I=1,NELM
    U=DSQRT(DP(LN(I,1),1)**2+DP(LN(I,1),2)**2)
    V=DSQRT(DP(LN(I,2),1)**2+DP(LN(I,2),2)**2)
    VE=MAX(U,V)
    CN(I)=0.5D0*VE*DELTTIME/JCB(I)
  END DO
  CFL=MAXVAL(CN)
  CFLOC=MAXLOC(CN,1)
  
	WRITE(23,*) TIME,CFL

    IF (CFL>=250.D0)THEN
    WRITE(22,*) TIME,"CFL=",CFL,"@ ELEMENT",CFLOC
    STOP
    END IF    

	RETURN
	END
!**********************************************************************
SUBROUTINE LOBATTO(N,X,W,A,B)
!**********************************************************************
IMPLICIT NONE
INTEGER I,ITER,N
REAL*8 A,B,C,D,X(N),W(N)
REAL*8 X0,X1,E,P1,P2,P3,PI
REAL*8,PARAMETER::EMAX=0.00000000001D0
REAL*8,EXTERNAL::LEGENDRE

PI=DACOS(-1.D0)
X(1)=-1.D0
W(1)=2.D0/N/(N-1)
X(N)=1.D0
W(N)=W(1)

DO I=2,N-1
ITER=1
E=10.D0
X0=(1.D0-3.D0*(N-2)/8.D0/(N-1)**3)*DCOS((4.D0*I-3)*PI/(4.D0*N-3))
    DO WHILE (E>=EMAX.AND.ITER<=1000)
    P1=(LEGENDRE(N-2,X0)-X0*LEGENDRE(N-1,X0))*(N-1)/(1.0-X0**2)
    P2=(2.D0*X0*P1-N*(N-1)*LEGENDRE(N-1,X0))/(1.0-X0**2)
    P3=(2.D0*X0*P2-(N*(N-1)-2)*P1)/(1.0-X0**2)
    X1=X0-2.D0*P1*P2/(2.0*P2**2-P1*P3)
    E=ABS(X1-X0)
    ITER=ITER+1
    X0=X1
    END DO
X(N-I+1)=X1
END DO

DO I=2,N-1
W(I)=2.D0/N/(N-1)/LEGENDRE(N-1,X(I))**2
END DO

C=(B-A)/2.D0
D=(B+A)/2.D0
DO I=1,N
X(I)=C*X(I)+D
W(I)=C*W(I)
END DO

RETURN
END
!**********************************************************************
FUNCTION LEGENDRE(N,X)
!**********************************************************************
INTEGER N
REAL*8 X,FI,PI,PIM1,PIM2
REAL*8 LEGENDRE
INTEGER I
IF (N.EQ.0) THEN
    LEGENDRE=1
    ELSEIF (N.EQ.1) THEN
    LEGENDRE=X
    ELSE
    PIM1=1
    PI=X
        DO I=2,N
        FI=I
        PIM2=PIM1
        PIM1=PI
        PI=((I+I-1)*X*PIM1-(I-1)*PIM2)/FI
        END DO
    LEGENDRE=PI
ENDIF
END
      
subroutine savgol_filter(nl,nr,ld,m,n1,y,flag)
!-----------------------------------------------------------------------------------
! This routine is used to perform the Savitzky-Golay algorithm.
!-----------------------------------------------------------------------------------
!    nl:: input, integer, the number of leftward data points used.
!    nr:: input, integer, the number of rightward data points used.
!    ld:: input, integer, the order of the derivative desired.
!     m:: input, integer, the order of the smoothing polynomial.
!    n1:: input, integer, the number of data points.
! y(n1):: input/output, real values, the data to be smoothed.
!  flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.26.
!-----------------------------------------------------------------------------------
! Dependence:: subroutine savgol.
! -----------------------------------------------------------------------------------

    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m, n1
    real   (kind=8), intent(inout):: y(n1)
    integer(kind=4), intent(out):: flag
    ! Local variables.
    integer(kind=4):: i, j, xl(nl+nr+1)
    real   (kind=8):: y0(n1), coef(nl+nr+1)

    xl(1) = 0

    y0 = y

    do i=1, nl

        xl(i+1) = -i

    end do

    do i=1, nr

        xl(1+nl+i) = nr-i+1

    end do

    call savgol(nl,nr,ld,m,coef,flag)

    if (flag/=0) return

    do i=1, n1-nr

        y(i) = 0.0

        do j=1, nl+nr+1

            if (i+xl(j) .gt. 0) then

                y(i) = y(i) + coef(j)*y0(i+xl(j))

            end if

        end do

    end do

    if (ld==0) then

        y(1:nl) = y0(1:nl)

        y(n1-nr+1:n1) = y0(n1-nr+1:n1)

    else 

        y(1:nl) = y(nl+1)

        y(n1-nr+1:n1) = y(n1-nr)
 
    end if

    return

end subroutine savgol_filter

subroutine savgol(nl,nr,ld,m,coef,flag)
!-----------------------------------------------------------------------------------
! This routine is used to calculate a set of Savitzky-Golay filter coefficients.
!-----------------------------------------------------------------------------------
!            nl:: input, integer, the number of leftward data points used.
!            nr:: input, integer, the number of rightward data points used.
!            ld:: input, integer, the order of the derivative desired.
!             m:: input, integer, the order of the smoothing polynomial.
! coef(nl+nr+1):: output, real values, calculated coefficents in wrap-around order.
!          flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.20.
!-----------------------------------------------------------------------------------
! Dependence:: subroutine ludcmp;
!              subroutine lubksb.
!-----------------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.646 IN Press et al.
! -----------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m
    real   (kind=8), intent(inout):: coef(nl+nr+1)
    integer(kind=4), intent(out):: flag
    ! Local variables.
    integer(kind=4):: imj, ipj, k, kk, mm, indx(m+1)

    real   (kind=8):: d, fac, summ, a(m+1,m+1), b(m+1)

    flag = 0

    if (nl < 0 .or. nr < 0 .or. ld > m .or. nl+nr < m) then

        flag = 1

        return

    end if

    do ipj=0, 2*m

        summ = 0.0

        if (ipj .eq. 0) summ = 1.0

        do k=1, nr

            summ = summ + (float(k))**ipj

        end do

        do k=1, nl

            summ = summ + (float(-k))**ipj

        end do

        mm = min(ipj, 2*m-ipj)

        do imj=-mm, mm, 2

            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = summ

        end do

    end do

    call ludcmp(a,m+1,indx,d,flag)

    if (flag .ne. 0) return

    b = 0.0

    b(ld+1) = 1.0

    call lubksb(a,m+1,indx,b)

    coef = 0.0

    do k=-nl, nr

        summ = b(1)

        fac = 1.0

        do mm=1, m

            fac = fac * k

            summ = summ + b(mm+1) * fac

        end do

        kk = mod(nl+nr+1-k, nl+nr+1) + 1

        coef(kk) = summ

    end do

    return

end subroutine savgol

subroutine lubksb(a,n,indx,b)
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, the LU decomposition of a matrix.
!       n:: input, integer, the dimenstion of the matrix.
! indx(n):: input,  integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!    b(n):: output, real values, the solution vector X for 
!                   linear equations A*X=B.
!-------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.18.
!-------------------------------------------------------------------------
! Dependence:: No.--------------------------------------------------------
!-------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.39 IN Press et al.
! -------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n, indx(n)
    real   (kind=8), intent(in):: a(n,n)
    real   (kind=8), intent(inout):: b(n)
   ! Local variables.
    integer(kind=4):: i, ii, j, ll
    real   (kind=8):: summ

    ii = 0

    do i=1, n

        ll = indx(i)

        summ = b(ll)

        b(ll) = b(i)

        if (ii .ne. 0) then

            do j=ii, i-1

                summ = summ - a(i,j) * b(j)

            end do

        else if (summ .ne. 0.0) then
            
            ii = i

        end if

        b(i) = summ

    end do

    do i=n, 1, -1

        summ = b(i)

        do j=i+1, n

            summ = summ - a(i,j) * b(j)

        end do

        b(i) = summ / a(i,i)

    end do

    return

end subroutine lubksb

subroutine ludcmp(a,n,indx,d,flag)
!-------------------------------------------------------------------------
!This routine is used in combination with lubksb to solve 
!linear equations or invert a matrix.
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, a matrix to be decomposed.
!       n:: input, integer, the dimension of the matrix.
! indx(n):: output, integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!       d:: output, integer, output as 1 or -1 depending on whether 
!           the number of row interchanges was even or odd.
!    flag:: output, integer, error message, 0=success, 1=singular matrix.
!-------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.20.
!-------------------------------------------------------------------------
! Dependence:: No.--------------------------------------------------------
!-------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.38 IN Press et al.
! ------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n
    integer(kind=4), intent(out):: indx(n), flag
    real   (kind=8), intent(inout):: a(n,n)
    real   (kind=8), intent(out):: d
    ! Local variables.
    integer(kind=4):: i, j, k, imax
    real   (kind=8):: aamax, dum, summ, vv(n)
   
    indx = 0

    flag = 0

    d = 1.0

    do i=1, n

        aamax = 0.0

        do j=1, n

            if (abs(a(i,j)) .gt. aamax)  aamax = abs(a(i,j))

        end do

        if (aamax .eq. 0.0) then 
 
            flag = 1

            return

        end if

        vv(i) = 1.0/aamax

    end do

    do j=1, n

        do i=1, j-1

            summ = a(i,j)

            do k=1, i-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

        end do

        aamax = 0.0

        do i=j, n

            summ = a(i,j)

            do k=1, j-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

            dum = vv(i) * abs(summ)

            if (dum .ge. aamax) then
       
                imax = i
          
                aamax = dum

            end if

        end do

        if (j .ne. imax) then

            do k=1, n

                dum = a(imax,k)

                a(imax,k) = a(j,k)

                a(j,k) = dum
  
            end do

            d = -d

            vv(imax) = vv(j)

        end if

        indx(j) = imax

        if (a(j,j) .eq. 0.0) a(j,j) = tiny(0.0D+00)

        if (j .ne. n) then

            dum = 1.0 / a(j,j)

            do i=j+1, n

                a(i,j) = a(i,j) * dum

            end do

        end if

    end do

    return

    end subroutine ludcmp    
      subroutine SEVAL (N,U,X,Y,B,C,D,V_0,V_1,V_2)
!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     SEVAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------
      REAL *8 B(N),C(N),D(N),X(N),Y(N),U,DX
      real *8 V_0,V_1,V_2
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I)) GO TO 101
      IF (U.LE.X(I+1)) GO TO 301
101 I = 1
      J = N+1
201 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 201

!     SPLINE EVALUATION

   301 DX = U-X(I)
!      SEVAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      V_0 = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      V_1 = B(I)+2.D0*C(I)*DX+3.D0*D(I)*DX**2
      V_2 = 2.D0*C(I)+6.D0*D(I)*DX
      RETURN
      END

      SUBROUTINE SPLINE (N,X,Y,B,C,D)
!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION B(N),C(N),D(N),X(N),Y(N)
      NM1 = N-1
      IF (N.LT.2) RETURN
      IF (N.LT.3) GO TO 501

!     BUILD THE TRIDIAGONAL SYSTEM
!     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO I = 2,NM1
      D(I) = X(I+1)-X(I)
      B(I) = 2.D0*(D(I-1)+D(I))
      C(I+1) = (Y(I+1)-Y(I))/D(I)
      C(I) = C(I+1)-C(I)
      END DO

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.D0
      C(N) = 0.D0
      IF (N.EQ.3) GO TO 151
      C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

151 DO I = 2,N
      T = D(I-1)/B(I-1)
      B(I) = B(I)-T*D(I-1)
      C(I) = C(I)-T*C(I-1)
      END DO

!     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO L = 1,NM1
      I = N-L
      C(I) = (C(I)-D(I)*C(I+1))/B(I)
      END DO

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
      DO I = 1,NM1
      B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
      D(I) = (C(I+1)-C(I))/D(I)
      C(I) = 3.D0*C(I)
      END DO
      C(N) = 3.D0*C(N)
      D(N) = D(NM1)
      RETURN

!     CAS N = 2

501 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.D0
      D(1) = 0.D0
      B(2) = B(1)
      C(2) = 0.D0
      D(2) = 0.D0
      RETURN
    END
!********************************************************************
 subroutine periodic(h,wave_period,wave_height,wave_amplitude)
!********************************************************************
! the subroutine is used to give the stroke of the wavemaker
 implicit none
 integer i
 real*8 pi,k0,k1,hs,s,e
 real*8 h,wave_period,wave_height,wave_amplitude
 pi=dacos(-1.d0)
    k0 = 1.d0
    e = 1.d0
    do while (e > 1.e-7)
        k1 = k0-(k0*dtanh(k0*h)-(2.d0*pi/wave_period)**2/9.81)/(dtanh(k0*h)+k0-k0*(dtanh(k0*h))**2)
        e = dabs((k1-k0)/k0)
        hs = 2.d0*(dcosh(2.d0*k1*h)-1.d0)/(dsinh(2.d0*k1*h)+2.d0*k1*h)
        s = wave_height/hs
        wave_amplitude = 0.5d0*s
        k0 = k1
    end do
    
    return
    end
!********************************************************************
    subroutine solitary(h0,wave_H,g,Npp,total_Nc,total_step,dt,dis,vb_x,acc)
!********************************************************************
    implicit none
    integer i,Npp,step,total_Nc,total_step,No(Npp)
    real*8 h0,wave_H,wave_K,wave_C,g,dt,x0
    real*8 vb_x(total_Nc,0:total_step-1),vb_z(total_Nc,0:total_step-1),x(total_Nc,0:total_step-1)
    real*8 dis(total_Nc,0:total_step-1),acc(total_Nc,0:total_step-1)
    
    call fenton_parameters(h0,wave_H,wave_K,wave_C,g)
    
    x0=-5.d0*wave_c
    x=0.d0
    do i=1,Npp
        No(i)=i
    end do

    do step = 0,total_step-1
        call wmbc(total_Nc,total_step,vb_x,vb_z,step,dt,wave_H,wave_K,wave_C,h0,x0,x,g,No,Npp)
    end do
    
    call get_wmk_dis_acc(total_Nc,total_step,vb_x,dt,dis,acc)

    return
    end
!********************************************************************
    subroutine get_wmk_dis_acc(total_Nc,total_step,vb_x,dt,dis,acc)
!********************************************************************
    implicit none
    integer i,j,total_Nc,total_step
    real*8 vb_x(total_Nc,0:total_step-1),dt,dis(total_Nc,0:total_step-1),acc(total_Nc,0:total_step-1)
    
    dis=0.d0
    acc=0.d0
    do i=1,total_Nc
        do j=1,total_step-1
            dis(i,j)=dis(i,j-1)+0.5d0*dt*(vb_x(i,j-1)+vb_x(i,j))
            if (j==total_step-1)then
                acc(i,j)=acc(i,j-1)
            else if (j==0)then
                acc(i,j)=(vb_x(i,j+1)-vb_x(i,j))/dt
            else
                acc(i,j)=0.5d0*(vb_x(i,j+1)-vb_x(i,j-1))/dt
            end if
        end do
    end do

    return
    end
!********************************************************************
    subroutine wmbc(total_Nc,total_step,vb_x,vb_z,step,dt,wave_H,wave_K,wave_C,h0,x0,x,g,No,Npp)
!********************************************************************
    ! the subroutine is used to give the velocity of the wavemaker
    ! x(total_Nc,total_step) is the coordinate of the observation point at time t
    ! x0 = initial position of wave crest
    implicit none
   integer j,total_Nc,total_step,step
   integer::zz,Npp,No(Npp)
   real*8 vb_x(total_Nc,0:total_step-1),vb_z(total_Nc,0:total_step-1),dt,wave_H,wave_K,wave_C,h0,x0
   real*8 x(total_Nc,0:total_step-1),zeta,kx,t,alpha,kk,cc,g,XX,S,xxx
   
   t=step*dt
   alpha=wave_H/h0
   kk=wave_K*h0
   cc=wave_C/(g*h0)**0.5
   
   do j=1,total_Nc
      do zz=1,Npp
          
         if(j.eq.No(zz)) then
             xxx=x(j,step) ! xxx: xi = position of the wave paddle at time t (step)
             XX=(xxx-x0-wave_C*t)/h0 ! XX: X = xi - Ct -x0
             S=1.d0/cosh(kk*XX) ! sech(K*X)
              zeta=S**2*alpha
              zeta=zeta+(-.75*S**2+.75*S**4)*alpha**2
              zeta=zeta+(.625*S**2-1.8875*S**4+1.2625*S**6)*alpha**3
           zeta=zeta+(-1.36817*S**2+3.88033*S**4-4.68304*S**6+2.17088*S**8)*alpha**4
        zeta=zeta+(1.86057*S**2-7.45136*S**4+12.7637*S**6-11.4199*S**8+4.24687*S**10)*alpha**5
    zeta=zeta+(-2.57413*S**2+13.2856*S**4-31.1191*S**6+40.1068*S**8-28.4272*S**10+8.728*S**12)*alpha**6
zeta=zeta+(3.4572*S**2-22.782*S**4+68.258*S**6-116.974*S**8+120.49*S**10-71.057*S**12+18.608*S**14)*alpha**7
zeta=zeta+(-4.6849*S**2+37.67*S**4-139.28*S**6+301.442*S**8-411.416*S**10+355.069*S**12-180.212*S**14+41.412*S**16)*alpha**8
zeta=zeta+(6.191*S**2-60.57*S**4+269.84*S**6-712.125*S**8+1217.98*S**10-1384.37*S**12+1023.07*S**14-450.29*S**16+90.279*S**18)*alpha**9
              zeta=zeta*h0
              vb_x(j,step)=wave_C*zeta/(h0+zeta) ! the velocity of the wavemaker
              goto 111
              
         else
            vb_x(j,step)=0.d0
         end if
         
      end do
     111  vb_z(j,step)=0.d0
   end do
   
    end
!********************************************************************
 subroutine fenton_parameters(h0,wave_H,wave_K,wave_C,g)
!********************************************************************
! the subroutine is used to give the parameters of the solitary wave by Fenton theory
! wave_H = wave height
! h0 = still water depth
! g = gravity
      real*8 wave_H,wave_K,wave_C,h0,g,alpha
      alpha=wave_H/h0
      wave_K=sqrt(3./4.d0*alpha)
      wave_K=wave_K*(1-.625*alpha+.554688*alpha**2-.561535*alpha**3+.567095*alpha**4-.602969*alpha**5 &
                  +.624914*alpha**6-.670850*alpha**7+.700371*alpha**8)
      wave_C=1+alpha
      wave_C=wave_C-.05*alpha**2-.0428571*alpha**3-.0342857*alpha**4-.0315195*alpha**5-.0292784*alpha**6 &
          -.0268451*alpha**7-.0302634*alpha**8-.0219347*alpha**9
      wave_K=wave_K/h0
      wave_C=wave_C**0.5*(g*h0)**0.5
    end
!********************************************************************
	SUBROUTINE SIMP(N,X,W,A,B)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N
	REAL*8 A,B,X(N),W(N)

	DO I=1,N
	X(I)=A+(B-A)*(I-1)/(N-1)
	W(I)=(B-A)/(N-1)/3.D0
	END DO

	DO I=2,N-1,2
	W(I)=4.D0*W(I)
	END DO

	DO I=3,N-2,2
	W(I)=2.D0*W(I)
	END DO

	RETURN
    END
!**********************************************************************
SUBROUTINE TRAP(N,X,W,A,B)
!**********************************************************************
IMPLICIT NONE
INTEGER I,N
REAL*8 A,B,X(N),W(N)

DO I=1,N
X(I)=A+(B-A)*(I-1)/(N-1)
W(I)=(B-A)/(N-1)
	IF(I==1.OR.I==N)THEN
	W(I)=W(I)/2.D0
	END IF
END DO

RETURN
END
!********************************************************************
	SUBROUTINE INCLU_ANGLE(X1,X2,P,THETA)
!********************************************************************
	IMPLICIT NONE
	REAL*8 A,B,C,X1(2),X2(2),P(2),THETA
    
    A = DSQRT((X1(1)-P(1))**2+(X1(2)-P(2))**2)
    B = DSQRT((X2(1)-P(1))**2+(X2(2)-P(2))**2)
    C = DSQRT((X1(1)-X2(1))**2+(X1(2)-X2(2))**2)
    THETA = DACOS(0.5D0*(C**2-A**2-B**2)/A/B)
    
	RETURN
    END
!********************************************************************
	SUBROUTINE search_closest(P,N,Q,NI,ID,DST)
!********************************************************************
	IMPLICIT NONE
    INTEGER I,N,NI,J(1),ID
	REAL*8 P(2),Q(N,2),DX,DY,RD(N),DST
    
    DO I=1,N
        DX=Q(I,1)-P(1)
        DY=Q(I,2)-P(2)
        RD(I)=DSQRT(DX*DX+DY*DY)
    END DO
    DST=MINVAL(RD)
    J=MINLOC(RD)
    ID=NI+J(1)-1

	RETURN
	END