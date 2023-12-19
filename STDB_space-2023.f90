!     Spatial distribution of b values
!
!     A versatile b-value calculator. Used for calculating time and space variations of b-values, 
!              it includes eight Mc calculation methods and four b-value calculation methods. 
!              Additional functions include synthetic seismic catalogs for seismic events, and more.
!
!     CSTDB is an extremely simple b-value calculator that includes most algorithms on the market, and the 
!             grouping between different algorithms is flexible and efficient. Just import a seismic catalog.
!
!     Refer to:
!     Weicheng Gong, Huayuan Cheng, Yajing Gao, Qing Li, and Yunqiang Sun (2023). Spatial-Temporal Distribution of b-value in the Eastern Tibetan Plateau. (in revision)
!      
!     INPUT
!     Catalog Format (the following numbers represent column numbers):
!     1.Year 2.Month 3.Day 4.Hour 5.Minute 6.Latitude 7.Longitude 8.Magnitude
!
!     OUTPUT
!     b-value as a function of spaace: b_value.txt
!     Mc as a function of spaace: Mc.txt
    
    program bzhi

    implicit none
    DIMENSION AM(10000,2),TIME_FW1(5),TIME_FW2(5),NUMER2(1000000,4)
    DIMENSION INDEXX(1),CATLOG(1000000,3),FREQUENCY(100,2),JD(2),WD(2)
    DIMENSION CUMULATIVE(100,4),BZHU(10000,4),TIME(1000000,5),AZHU(10000)
    DIMENSION NUMER(1000000,4),BZHUR(10000,4),AMZ(10000,2),BZHUR1(10000,4)
    DIMENSION AZ(1,2),BZ(1,2),CZ(1,2),DZ(1,2),X(100),Y(100),A(3),&
              AZHUR(10000),Z(100),Z2(100),MMC(10000),YX(4),RR3(10),RR4(10),&
              LO1(10,2),LO2(10,2),LO3(10,2),LO4(10,2),LI2(10),ACL(10000),&
              NUMER22(1000000,4),LA1(4,2),MMC2(10000)
    DOUBLE PRECISION D,D1,D2,D3,MAXN,MINN,PP,FREQUENCY,B_DOUBLE,MAG,RE,R1,R2
    DOUBLE PRECISION MC,M_AVERAGE,AA,CATLOG,CUMULATIVE,BAM,BAM2,BZHUR1
    DOUBLE PRECISION BZHU,B_VALE,NUMER,BZHUR,AM,AMZ,NUMER2,JD,WD,AZ,BZ,CZ,DZ,&
                     A,X,Y,DT1,DT2,DT3,AZHU,AZHUR,MI,A8,Z,Z2,MMC,MCC,BEE,HMM,&
                     DI1,DI2,YX,XZJ,R3,R4,RR3,RR4,LO1,LO2,LO3,LO4,CAA,CAA1,CAA2,&
                     NIU,ACL,LGN,MEAX,MEAY,X1,X2,Y1,Y2,PROJ_X,PROJ_Y,MK,BK,DISTAN,&
                     NUMER22,LA1,MMC2,LID,LIU,MAX_BB,BTZ,MIN_BB,BY,JG,P
    INTEGER I1,I2,I,AD,BD,COUNT,J,K,K1,Q,INDEXX,INDEXX2,J11,TEMP,CCON,AKB
    INTEGER A1,A2,II,QQ,QQ1,II3,II4,PP2,KK,JJ,WQ,EW,WQQ,EWE,QWEE,QQW,KKK
    INTEGER KII,TIME,CCN,TIME_FW1,TIME_FW2,I22,KT,IOP,TI,RG,RG2,IL,LIMI,&
            N,M,PO,BE,RF,MH,CG,GG,M1,YZZ,CC,GC,RRG,RR,LI1,LI2,LO,&
            NLA,PPL,BN,QWS,GRT,DIE,PQ,BIN,R
    CHARACTER*40 CATA
    CHARACTER*132 stamp
        
    WRITE (*,*) 'Graphic type:  1-Flat surface (recommended) 2-Vertical surface'
    READ (*,*) GRT
    
    IF (GRT.EQ.1) GOTO 11 
    WRITE (*,*) 'Section length (km)  '
    READ (*,*) JD(2)
    JD(1)=0
    
    WRITE (*,*) 'Depth range (km)-- Dep1:  Dep2: '
    READ (*,*) WD(1),WD(2)
    
    WRITE (*,*) 'Surface section (the first point)-- lat:    lon:  (From big to small)'
    READ (*,*) Y1,X1
    
    WRITE (*,*) 'Surface section (the second point)-- lat:    lon:  '
    READ (*,*) Y2,X2
    GOTO 889
    
11  WRITE (*,*) 'longitude range (100,108):'
    READ (*,*) JD(1),JD(2)

    
889 IF (GRT.EQ.2) GOTO 362
    WRITE (*,*) 'latitude range (28,35):'
    READ (*,*) WD(1),WD(2)
    
    WRITE (*,*) 'Whether to filter the output grid range: 1-YES  2-NO'
    READ (*,*) RG
    IF (RG.EQ.1) GOTO 363
    GOTO 362
363 WRITE(*,*) 'Range:   1-ellipse  2-rectangle'
    READ (*,*) RRG
    
    WRITE (*,*) 'LMSF-DATA:   1-YES   2-NO'
    READ (*,*) RG2
    
    IF (RG2.EQ.1) GOTO 362
    IF (RRG.EQ.1) GOTO 366
    WRITE (*,*) 'Lower-left point: example-(100,28)'
    READ (*,*) AZ(1,1),AZ(1,2)
    WRITE (*,*) 'Top-left point: '
    READ (*,*) BZ(1,1),BZ(1,2)
    WRITE (*,*) 'Top-rigth point: '
    READ (*,*) DZ(1,1),DZ(1,2)
    WRITE (*,*) 'Lower-right point: '
    READ (*,*) CZ(1,1),CZ(1,2) 
    GOTO 362
    
366 WRITE (*,*) 'Coordinates of the center of an ellipse (Latitude / longitude)'
    READ (*,*) YX(1),YX(2)
    WRITE (*,*) ' Length:   Height: '
    READ (*,*) YX(3),YX(4)
    WRITE (*,*) ' Rotation Angle (starting with parallel latitude/counterclockwise): '
    READ (*,*) XZJ
    
362 WRITE (*,*) 'Time range:  1-YES   2-NO'
    READ (*,*) TI
    IF (TI.EQ.1) GOTO 116
    GOTO 117
    
116 WRITE (*,*) 'Start time:  YEAR  MOUTH  DAY  HOUR  MINUTE'
    READ (*,*) TIME_FW1(1),TIME_FW1(2),TIME_FW1(3),TIME_FW1(4),TIME_FW1(5)
    WRITE (*,*) 'Deadline:  YEAR  MOUTH  DAY  HOUR  MINUTE'
    READ (*,*) TIME_FW2(1),TIME_FW2(2),TIME_FW2(3),TIME_FW2(4),TIME_FW2(5)
    
117 WRITE (*,*) 'Whether to filter the input directory range:  1-YES   2-NO'
    READ (*,*) DIE
    IF (DIE.EQ.1) GOTO 126
    GOTO 127
    
126 WRITE (*,*) 'Lower-left point--lat1:   lon1:  '
    READ (*,*) LA1(1,1),LA1(1,2)
    WRITE (*,*) 'Top-left point--lat2:   lon2: '
    READ (*,*) LA1(2,1),LA1(2,2)
    WRITE (*,*) 'Top-rigth point--lat3:   lon3: '
    READ (*,*) LA1(3,1),LA1(3,2)
    WRITE (*,*) 'Lower-right point--lat4:   lon4: '
    READ (*,*) LA1(4,1),LA1(4,2)    
    
127 WRITE (*,*) 'Each grid minimum seismic data  '
    READ (*,*) LIMI
    WRITE (*,*) 'Minimum magnitude  '
    READ (*,*) MAG
    WRITE (*,*) 'Mesh accuracy--X   unit    /km    '
    READ (*,*) MEAX
    WRITE (*,*) 'Mesh accuracy--Y   unit    /km    '
    READ (*,*) MEAY
    WRITE (*,*) 'Initial search radius  unit    /km    '
    READ (*,*) R1
    WRITE (*,*) 'Maximum search radius  unit    /km    '
    READ (*,*) R2
    WRITE (*,*) 'Number of search range segments  '
    READ (*,*) RR
    WRITE (*,*) 'Whether localised treatment is applied to some areas  1-YES 2-NO'
    READ (*,*) LO
    IF (LO.EQ.2) GOTO 118
    WRITE (*,*) 'Number of local areas: '
    READ (*,*) NLA
    
    DO 664 I=1,NLA
    WRITE (*,'("local areas",1X,I2)') I
    WRITE (*,*) 'Lower-left point: example-(100,28)'
    READ (*,*) LO1(I,1),LO1(I,2)
    WRITE (*,*) 'Top-left point: '
    READ (*,*) LO2(I,1),LO2(I,2)
    WRITE (*,*) 'Top-rigth point: '
    READ (*,*) LO3(I,1),LO3(I,2)
    WRITE (*,*) 'Lower-right point: '
    READ (*,*) LO4(I,1),LO4(I,2)
    WRITE (*,*) 'Each grid minimum seismic data (local areas',I,'):'
    READ (*,*) LI2(I)
    WRITE (*,*) 'Initial search radius (local areas',I,'):'
    READ (*,*) RR3(I)
    WRITE (*,*) 'Maximum search radius (local areas',I,'):'
    READ (*,*) RR4(I)
664 CONTINUE
    
118 WRITE (*,*) 'Calculation rule  1- b value  2- n value'
    READ (*,*) BN
    WRITE (*,*) 'b value calculation rule  1-Maximum likelihood method  2-Minimum quadratic fit'
    READ (*,*) BE
    WRITE (*,*) 'Mc Calculation method  1-EMR  2-MAXC  3-GFT  4-MBS  5-MBS(JS) 6-MBASS'
    READ (*,*) MH
    
    IF (MH.EQ.5) THEN
    WRITE (*,*) 'Box width (0.5¡ú5)£º'
    READ (*,*) BIN
    ELSE IF (MH.EQ.4) THEN
    WRITE (*,*) 'Stable platform (0.03-0.05)£º'
    READ (*,*) BY
    WRITE (*,*) 'MAGNITUDE INTERVAL (0.1)£º'
    READ (*,*) JG
    ELSE IF (MH.EQ.1) THEN
    WRITE (*,*) 'Significance of KS test (default:0.05)£º'
    READ (*,*) P
    ELSE IF (MH.EQ.3) THEN
    WRITE (*,*) 'Goodness of fit (90%¡ú90)£º'
    READ (*,*) R
    END IF
    
    WRITE (*,*) 'b-value limit (0.2,2.5)  '
    READ (*,*) LID,LIU
    WRITE (*,*) 'Enter catalog filename (year/month/day/hour/miute/lat/lon/magnitude/depth):'
    READ (*,'(a)') CATA
    
!   read file
888 COUNT=0
    CCN=0 
    OPEN(101,FILE=CATA)
    DO WHILE(NOT(EOF(101)))
    READ(101,*)
        COUNT=COUNT+1
    END DO
    CLOSE(101)
    WRITE(*,'("Total data:",1X,I8)') COUNT
    
    OPEN(101,FILE=CATA)
    DO K=1,COUNT
        READ(101,*) (TIME(K,KII),KII=1,5),NUMER2(K,1),NUMER2(K,2),NUMER2(K,3),NUMER2(K,4)
    END DO
    CLOSE(101)
    
    IF (TI.EQ.2) GOTO 115
    
!   Time Series Screening
    DO 100 K=1,COUNT
    IF (NUMER2(K,3).GE.MAG) GOTO 120
    GOTO 100
120 IF ((TIME(K,1).LE.TIME_FW2(1)).AND.(TIME(K,1).GE.TIME_FW1(1))) GOTO 10
    GOTO 100
10  IF (TIME(K,1).EQ.TIME_FW2(1)) GOTO 20
    IF (TIME(K,1).EQ.TIME_FW1(1)) GOTO 70
    GOTO 60
70  IF (TIME(K,2).GE.TIME_FW1(2)) GOTO 80 
    GOTO 100
80  IF (TIME(K,3).GE.TIME_FW1(3)) GOTO 90
    GOTO 100
90  IF (TIME(K,4).GE.TIME_FW1(4)) GOTO 110
    GOTO 100
110 IF (TIME(K,5).GE.TIME_FW1(5)) GOTO 60
    GOTO 100
    
20  IF (TIME(K,2).LE.TIME_FW2(2)) GOTO 30    
    GOTO 100
30  IF (TIME(K,3).LE.TIME_FW2(3)) GOTO 40
    GOTO 100
40  IF (TIME(K,4).LE.TIME_FW2(4)) GOTO 50
    GOTO 100
50  IF (TIME(K,5).LE.TIME_FW2(5)) GOTO 60
    GOTO 100
    
60  CCN=CCN+1
    NUMER(CCN,1)=NUMER2(K,1)
    NUMER(CCN,2)=NUMER2(K,2)
    NUMER(CCN,3)=NUMER2(K,3)
    NUMER(CCN,4)=NUMER2(K,4)
    
100 CONTINUE
    GOTO 54
    
115 DO 55 I=1,COUNT
    NUMER(I,1)=NUMER2(I,1)
    NUMER(I,2)=NUMER2(I,2)
    NUMER(I,3)=NUMER2(I,3)
    NUMER(I,4)=NUMER2(I,4)
55  CONTINUE
    CCN=COUNT
    
54  CONTINUE  
    
    IF (DIE.EQ.2) GOTO 45
    I1=0
    DO 87 I=1,CCN
    IF ((LA1(2,1)-LA1(1,1))*(NUMER(I,1)-LA1(1,2))-(LA1(2,2)-LA1(1,2))*(NUMER(I,2)-LA1(1,1)).LE.0) GOTO 181
    GOTO 87
181 IF ((LA1(3,1)-LA1(2,1))*(NUMER(I,1)-LA1(2,2))-(LA1(3,2)-LA1(2,2))*(NUMER(I,2)-LA1(2,1)).LE.0) GOTO 182
    GOTO 87
182 IF ((LA1(4,1)-LA1(3,1))*(NUMER(I,1)-LA1(3,2))-(LA1(4,2)-LA1(3,2))*(NUMER(I,2)-LA1(3,1)).LE.0) GOTO 183
    GOTO 87
183 IF ((LA1(1,1)-LA1(4,1))*(NUMER(I,1)-LA1(4,2))-(LA1(1,2)-LA1(4,2))*(NUMER(I,2)-LA1(4,1)).LE.0) GOTO 184
    GOTO 87
184 I1=I1+1
    NUMER22(I1,1)=NUMER(I,1)
    NUMER22(I1,2)=NUMER(I,2)
    NUMER22(I1,3)=NUMER(I,3)
    NUMER22(I1,4)=NUMER(I,4)
87  CONTINUE
    
    
    NUMER(1:I1,1:4)=NUMER22(1:I1,1:4)
    CCN=I1
    
    OPEN(113,FILE='SYSJ.txt')
    DO I=1,CCN
    WRITE(113,'(F10.2,F10.2,F10.2)') NUMER(I,1),NUMER(I,2),NUMER(I,3)
    END DO
    
45  CONTINUE  
    
!   Section processing
    IF (GRT.EQ.1) GOTO 52
    MK=(Y2-Y1)/(X2-X1)
    BK=Y1-MK*X1
    
    DO 53 I=1,CCN
    PROJ_X=(MK*(NUMER(I,1)-BK)+NUMER(I,2))/(MK**2+1)
    PROJ_Y=MK*PROJ_X+BK
    DISTAN=SQRT((((X1-PROJ_X)*111)**2)+(((Y1-PROJ_Y)*111)**2))
    NUMER(I,1)=NUMER(I,4)
    NUMER(I,2)=DISTAN
53  CONTINUE

!   Plane treatment
52  J=0
    IOP=0
    J11=0
    QQ=0
    CCON=0
    TEMP=0
    BD=NINT(ABS(WD(2)-WD(1))/MEAY)
    AD=NINT(ABS(JD(2)-JD(1))/MEAX)

!   Grid creation
!   Screening criteria for different regions
    LI1=LIMI
    R3=R1
    R4=R2

    RE=(R2-R1)*1.0/(RR*1.0)  
    RR=RR+1
    DO 552 I1=1,BD
    DO 551 I2=1,AD
    KT=1
    IOP=0
    
250 DO 251 AKB=1,RR
    J=0
    Q=0
    CATLOG=0.0
    
    DO 550 I=1,CCN
    D1=WD(1)+MEAY*I1
    D=ABS(NUMER(I,1)-D1-MEAY/2)
    
    IF (D.LE.(R1+(AKB-1)*RE)) GOTO 226   ! Determine whether it is within the search radius
    GOTO 550
226 D3=JD(1)+MEAX*I2
    D2=ABS(NUMER(I,2)-D3-MEAX/2)
    
    IF (D2.LE.(R1+(AKB-1)*RE)) GOTO 227
    GOTO 550
227 J=J+1
    Q=Q+1
    CATLOG(Q,1)=NUMER(I,1)
    CATLOG(Q,2)=NUMER(I,2)
    CATLOG(Q,3)=NUMER(I,3)
550 CONTINUE
    
!   special area processing
    IF (LO.EQ.2) GOTO 38
    DO 37 PPL=1,NLA
    IF ((LO2(PPL,1)-LO1(PPL,1))*(D1-0.1-LO1(PPL,2))-(LO2(PPL,2)-LO1(PPL,2))*(D3-0.1-LO1(PPL,1)).LE.0) GOTO 97
    GOTO 37
97  IF ((LO3(PPL,1)-LO2(PPL,1))*(D1-0.1-LO2(PPL,2))-(LO3(PPL,2)-LO2(PPL,2))*(D3-0.1-LO2(PPL,1)).LE.0) GOTO 98
    GOTO 37
98  IF ((LO4(PPL,1)-LO3(PPL,1))*(D1-0.1-LO3(PPL,2))-(LO4(PPL,2)-LO3(PPL,2))*(D3-0.1-LO3(PPL,1)).LE.0) GOTO 96
    GOTO 37
96  IF ((LO1(PPL,1)-LO4(PPL,1))*(D1-0.1-LO4(PPL,2))-(LO1(PPL,2)-LO4(PPL,2))*(D3-0.1-LO4(PPL,1)).LE.0) GOTO 36
    GOTO 37
36  LIMI=LI2(PPL)
    R1=RR3(PPL)
    R2=RR4(PPL)
    GOTO 38
37  CONTINUE
    LIMI=LI1
    R1=R3
    R2=R4
    
38  CONTINUE  
    GOTO 300
    
!   Expand the search radius
500 IOP=IOP+1
    IF (IOP.GE.RR) GOTO 350 
    FREQUENCY=0
    CUMULATIVE=0
    KT=KT+1
    J11=J11-1
    GOTO 250
    
300 IF (J.GE.LIMI) GOTO 301
    GOTO 251
301 IF (AKB.EQ.KT) EXIT
251 CONTINUE
    
!   n value calculation
    IF (BN.EQ.1) GOTO 69
    IF (J.LT.LIMI) GOTO 56
    J11=J11+1
    CAA1=0.
    CAA2=0.
    
    MINN=MINval(CATLOG(1:Q,3))
    DO 78 I=1,J
    CAA=CATLOG(I,3)-MINN
    CAA1=CAA1+CAA*CAA/J
78  CAA2=CAA2+CAA/J
    
    NIU=CAA1/(CAA2*CAA2)
    B_VALE=NIU
    B_DOUBLE=0
    MC=0
    IF ((NIU.GE.0.5).AND.(NIU.LE.5)) GOTO 918
    GOTO 500
    
############---b value main program---############
##################################################
69  IF (J.GE.LIMI) GOTO 66  ! Judging whether the minimum number of grids is met
    GOTO 56
66  J11=J11+1
    MAXN=MAXval(CATLOG(1:Q,3))
    MINN=MINval(CATLOG(1:Q,3))

!   data sorting   
    DO 99 A1=1,Q-1
    DO 89 A2=A1+1,Q
    IF(CATLOG(A1,3).GT.CATLOG(A2,3)) GOTO 216
    GOTO 89
216 AA=CATLOG(A1,3)
    CATLOG(A1,3)=CATLOG(A2,3)
    CATLOG(A2,3)=AA
89  CONTINUE
99  CONTINUE
    
!   Array remove duplicates
    QQ=0
    DO 886 II=1,Q-1
    IF(ABS(CATLOG(II,3)-CATLOG(II+1,3)).LE.0.001) GOTO 46
    GOTO 47
46  CYCLE
47  QQ=QQ+1
    FREQUENCY(QQ,1)=CATLOG(II,3)
886 CONTINUE
    
    QQ1=QQ+1
    FREQUENCY(QQ1,1)=CATLOG(Q,3)
    
!   statistics array  
    DO 778 II3=1,QQ1
    DO 777 II4=1,Q
    IF(ABS(FREQUENCY(II3,1)-CATLOG(II4,3)).LE.0.001) GOTO 887
    GOTO 777
887 CCON=CCON+1
777 CONTINUE
    FREQUENCY(II3,2)=CCON
    CCON=0
778 CONTINUE
    TEMP=0     

!   b value calculation
    DO 666 PP=MAXN,MINN-0.1,-0.1
    TEMP=TEMP+1
    CUMULATIVE(TEMP,1)=PP+0.00001
    
    DO 665 JJ=1,QQ1
    IF (ABS(PP-FREQUENCY(JJ,1)).LE.0.001) GOTO 48
    GOTO 665
48  CUMULATIVE(TEMP,2)=FREQUENCY(JJ,2)
665 CONTINUE    
    
    CUMULATIVE(TEMP,3)=CUMULATIVE(TEMP,1)*CUMULATIVE(TEMP,2)
    CUMULATIVE(TEMP,4)=SUM(CUMULATIVE(1:TEMP,2))
    
666 CONTINUE
    
    DO 456 JJ=TEMP,1,-1
    IF (CUMULATIVE(JJ,2).EQ.0) GOTO 456
    HMM=CUMULATIVE(JJ,1)+0.1
    EXIT
456 CONTINUE
    
!   Find the maximum event corresponding to Mc
    INDEXX=MAXLOC(CUMULATIVE(:,2))
    INDEXX2=INDEXX(1)
    
    M=2
    M1=M+1
    
!   Judgment b method    
    IF (BE.EQ.2) GOTO 555
    GOTO 556
    
!   Minimum quadratic fit
555 DO 554 PO=1,INDEXX2
    X(PO)=CUMULATIVE(PO,1)
    Y(PO)=LOG10(CUMULATIVE(PO,4))*1.0
554 CONTINUE
    N=INDEXX2
    GOTO 557
    
!   Maximum likelihood method  
556 CONTINUE
    DO 558 PO=1,TEMP
    X(PO)=CUMULATIVE(PO,1)
    Y(PO)=CUMULATIVE(PO,4)
    Z(PO)=CUMULATIVE(PO,2)
    Z2(PO)=CUMULATIVE(PO,3)
558 CONTINUE
    N=TEMP
    QWS=INT(Y(TEMP))
    
!   Judgment Mc method
    IF (MH.EQ.3) GOTO 469
    IF (MH.EQ.2) GOTO 912
    IF ((MH.EQ.4).OR.(MH.EQ.5)) GOTO 916
    IF (MH.EQ.1) GOTO 917
    IF (MH.EQ.6) GOTO 415
    
!   best fit method
469 CALL GFT(X,Y,Z2,A,N,M,R,MC,GG,B_VALE,B_DOUBLE,HMM,LIMI)
    IF (GG.NE.1) GOTO 500
    GOTO 918

!   (MBASS)
415 CALL MBASS(X,Z,N,MC,INDEXX2,CG)
    IF (CG.NE.1) GOTO 500
    M_AVERAGE=SUM(Z2(1:INDEXX2))*1.0/Y(INDEXX2)
    GOTO 913    
    
!   Mc by b-value stability 
916 CALL MBS(X,Y,Z2,N,B_VALE,B_DOUBLE,CC,MH,MC,HMM,LIMI,BIN,BY,JG)
    IF (CC.EQ.0) GOTO 500
    GOTO 918
    
!   Entire-magnitude-range method
917 CALL EMR(X,Y,Z,Z2,A,N,M,CC,B_VALE,B_DOUBLE,MC,HMM,QWS,P) 
    IF (CC.EQ.0) GOTO 500
    GOTO 918
    
!   Least squares fit to calculate b value  
557 CALL NH_QBXF(X,Y,N,A,M,M1)
    
    IF (BE.EQ.2) GOTO 911 
    GOTO 912
911 IF ((-A(2).LE.0.0).OR.(-A(2).GE.3.0)) GOTO 500

!   Maximum curvature-method
912 MC=CUMULATIVE(INDEXX2,1)
    LGN=CUMULATIVE(INDEXX2,4)
    M_AVERAGE=SUM(CUMULATIVE(1:INDEXX2,3))*1.0/CUMULATIVE(INDEXX2,4)*1.0
    
913 B_VALE=(M_AVERAGE-MC)*2.30258509
    B_VALE=1.0/B_VALE*1.0
    B_DOUBLE=B_VALE*1.0/SQRT(CUMULATIVE(INDEXX2,4))*1.0
    
918 BZHU(J11,1)=D3-MEAX
    BZHU(J11,2)=D1-MEAY
    IF (B_VALE >= MIN_BB .AND. B_VALE <= MAX_BB) THEN
    B_VALE = B_VALE + BTZ
    END IF
    BZHU(J11,3)=B_VALE
    BZHU(J11,4)=B_DOUBLE
    AZHU(J11)=-A(2)
    MMC(J11)=MC
    AM(J11,1)=J    
    AM(J11,2)=MAXN-MINN
    ACL(J11)=LGN+B_DOUBLE*MC

                
    IF (BN.EQ.2) GOTO 67
    IF ((B_VALE.LE.LID).OR.(B_VALE.GE.LIU)) GOTO 500 ! Determine whether the b-value limit is met
67  CALL fdate (stamp)
    IF (BE.EQ.2) WRITE (*,'(2X,A,2X,I5,2X,2(1X,F6.2),2X,F6.2)') stamp(11:20),&
                            J11,(BZHU(J11,I22),I22=1,2),AZHU(J11)
    IF (BE.EQ.1) WRITE (*,'(2X,A,2X,I5,2X,2(1X,F6.2),1X,F6.2,2(1X,F6.2))') stamp(11:20),&
                            J11,(BZHU(J11,I22),I22=1,2),MC,(BZHU(J11,I22),I22=3,4)
56  CONTINUE
    
!   clear data group
    CATLOG=0
    FREQUENCY=0
    CUMULATIVE=0
    J=0
    Q=0
    GOTO 360
350 J11=J11-1
    
360 CONTINUE
551 CONTINUE
552 CONTINUE
    
!   b value main program calculation ends 
!####################################################
    
!####################################################
!   Determine whether to output range filtering
    IF (RG.EQ.1) GOTO 337
    GOTO 334
337 IF (RG2.EQ.1) GOTO 338
    GOTO 335
    
!   out rang (LMSF)
338 IF (RRG.EQ.1) GOTO 383
    AZ(1,1)=102
    AZ(1,2)=29
    BZ(1,1)=101
    BZ(1,2)=30
    CZ(1,1)=105.8
    CZ(1,2)=33.6
    DZ(1,1)=106.8
    DZ(1,2)=32.6

!   Rectange range output
335 IF (RRG.EQ.1) GOTO 384
    CONTINUE
    WQ=0
    DO 996 EW=1,J11
    IF ((BZ(1,1)-AZ(1,1))*(BZHU(EW,2)-AZ(1,2))-(BZ(1,2)-AZ(1,2))*(BZHU(EW,1)-AZ(1,1)).LE.0) GOTO 997
    GOTO 996
997 IF ((CZ(1,1)-BZ(1,1))*(BZHU(EW,2)-BZ(1,2))-(CZ(1,2)-BZ(1,2))*(BZHU(EW,1)-BZ(1,1)).LE.0) GOTO 998
    GOTO 996
998 IF ((DZ(1,1)-CZ(1,1))*(BZHU(EW,2)-CZ(1,2))-(DZ(1,2)-CZ(1,2))*(BZHU(EW,1)-CZ(1,1)).LE.0) GOTO 999
    GOTO 996
999 IF ((AZ(1,1)-DZ(1,1))*(BZHU(EW,2)-DZ(1,2))-(AZ(1,2)-DZ(1,2))*(BZHU(EW,1)-DZ(1,1)).LE.0) GOTO 995
    GOTO 996
995 WQ=WQ+1
    BZHUR(WQ,1)=BZHU(EW,1)
    BZHUR(WQ,2)=BZHU(EW,2)
    BZHUR(WQ,3)=BZHU(EW,3)
    BZHUR(WQ,4)=BZHU(EW,4)
    AZHUR(WQ)=AZHU(EW)
    MMC2(WQ)=MMC(EW)
996 CONTINUE
    GOTO 895
    
!   out rang (LMSF)
383 YX(1)=104
    YX(2)=31.3
    YX(3)=3
    YX(4)=0.9
    XZJ=38

!   Ellipse range output
384 WQ=0
    DO 970 EW=1,J11
    DI1=(((BZHU(EW,1)-YX(1))*COSD(XZJ)+(BZHU(EW,2)-YX(2))*SIND(XZJ))**2)/(YX(3)**2)
    DI2=(((BZHU(EW,1)-YX(1))*SIND(XZJ)-(BZHU(EW,2)-YX(2))*COSD(XZJ))**2)/(YX(4)**2)
    IF ((DI1+DI2).GT.1) GOTO 970
    WQ=WQ+1
    BZHUR(WQ,1)=BZHU(EW,1)
    BZHUR(WQ,2)=BZHU(EW,2)
    BZHUR(WQ,3)=BZHU(EW,3)
    BZHUR(WQ,4)=BZHU(EW,4)
    AZHUR(WQ)=AZHU(EW)
    MMC2(WQ)=MMC(EW)
970 CONTINUE
    GOTO 895
    
334 DO 991 IL=1,J11
    BZHUR(IL,1)=BZHU(IL,1)
    BZHUR(IL,2)=BZHU(IL,2)
    BZHUR(IL,3)=BZHU(IL,3)
    BZHUR(IL,4)=BZHU(IL,4)
    AZHUR(IL)=AZHU(IL)
    MMC2(IL)=MMC(IL)
991 CONTINUE
      
    WQ=J11
    
#############---output file program---##############
####################################################
895 WRITE(*,'("Total seismic data:",1X,I8)') COUNT
    WRITE(*,'("Number of seismic data in the time range:",1X,I8)') CCN
    WRITE(*,'("The total number of grid:",1X,I8)') J11
    WRITE(*,'("Within the scope of effective number of grid:",1X,I8)') WQ
    
    OPEN(103,FILE='b_value.txt')
    DO 612 KK=1,WQ
    IF (BE.EQ.1) GOTO 611
    WRITE(103,'(F10.2,F10.2,F10.2)') BZHUR(KK,1),BZHUR(KK,2),AZHUR(KK)
    GOTO 612
611 IF (GRT.EQ.2) THEN
        BZHUR(KK,2)=-BZHUR(KK,2)
    END IF
    WRITE(103,'(F10.2,F10.2,F10.2)') BZHUR(KK,1),BZHUR(KK,2),BZHUR(KK,3)
612 CONTINUE
    CLOSE(103)
    
    BEE=0.
    OPEN (104,FILE='b_error.txt')
    DO 401 KK=1,WQ
    WRITE(104,'(F10.2,F10.2,F10.2)') BZHUR(KK,1),BZHUR(KK,2),BZHUR(KK,4)
    BEE=BEE+BZHUR(KK,4)
401 CONTINUE
    CLOSE (104)
    
    MCC=0.
    OPEN (185,FILE='Mc.txt')
    DO 402 KK=1,WQ
    WRITE(185,'(3(F10.2,2X))') BZHUR(KK,1),BZHUR(KK,2),MMC2(KK)
    MCC=MCC+MMC(KK)
402 CONTINUE
    CLOSE (185)
    
    WRITE(*,'("Average Mc magnitude:",1X,F6.2)') MCC/WQ*1.0
    WRITE(*,'("Average B-value uncertainty:",1X,F6.2)') BEE/WQ*1.0
    
    OPEN (195,FILE='Average.txt')
    WRITE (195,'(F6.2,2X,F6.2)') MCC/WQ*1.0,BEE/WQ*1.0
    CLOSE (195)
    
    PAUSE
    end program bzhi

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCC----CDFX-----CCCCC   
      SUBROUTINE CDFX(A,B,N,XIG,MU,S)
      DOUBLE PRECISION A,B,H,S,XIG,MU,X
      INTEGER N,I
      real*8, parameter :: pi = 3.1415926
      S=0.
      
      H=ABS(B-A)/REAL(N)
      DO I=1,N
      X=A+(I-1)*H
      S=S+1/(SQRT(2*PI)*XIG)*EXP((-(X-MU)**2)/(2*XIG**2))
      END DO
      
      RETURN
      END     
    
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
#CCCCC---------------Create rank---------------CCCCCCCCCC
      SUBROUTINE CEF(X,Y,N,TOP,IT)
      INTEGER N,I,TOP,IT,J
      DOUBLE PRECISION X(N),Y(IT),STEP
      
      DO I=1,IT
      NUMER=0
      DO J=1,N
      STEP=TOP*1.0/IT*I*1.0
      IF (STEP.GE.X(J)) NUMER=NUMER+1
      END DO
      Y(I)=NUMER*1.0/N*1.0
      END DO
      
      RETURN
      END SUBROUTINE CEF    
    
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
#CCCCC---------------Create rank---------------CCCCCCCCCC
    SUBROUTINE CRANK(X,Y,N)
    INTEGER N,I,J,NX(N),NX2(N),NUMER,I2,AVER,NU
    DOUBLE PRECISION X(N),Y(N),TEMP
    NX=0
    J=0
    
    DO I=1,N
    DO J=I+1,N
    IF (X(I)>X(J)) then
    TEMP=X(I)
    X(I)=X(J)
    X(J)=TEMP
    END IF
    END DO
    END DO
    
    DO I=1,N-1
    IF (X(I).EQ.X(I+1)) THEN
    NX(I)=I
    NX(I+1)=I+1
    END IF
    END DO
    
    J=0
    DO I=1,N-1
    IF (NX(I).NE.0) THEN
    J=J+1
    NX2(J)=NX(I)
    END IF
    END DO
    
    I=0
    DO WHILE (I<N)
10  I=I+1
    NU=0
    AVER=0
    DO I2=1,J
    IF (I+NU.EQ.NX2(I2)) THEN
    NU=NU+1
    END IF
    END DO
    IF (NU.NE.0) THEN
    DO I2=1,NU
    AVER=AVER+I-1+I2
    END DO
    DO I2=1,NU
    Y(I-1+I2)=AVER*1.0/NU*1.0
    END DO
    I=I+NU-1
    GOTO 10
    END IF
    Y(I)=I
    END DO
    
    RETURN
    END SUBROUTINE CRANK    
    
##########################################
#######----Least squares fit----##########
    SUBROUTINE NH_QBXF(X,Y,N,A,M,M1)
    DIMENSION X(N),Y(N),A(M1),IX(20),H(20)
    DOUBLE PRECISION X,Y,A,H,HA,HH,Y1,Y2,H1,H2,D,HM,IX
    INTEGER I,J,L,M,N,M1,II
    DO 75 I=1,M1
75   A(I)=0.0
    IF (M.GE.N) M=N-1
    IF (M.GE.20) M=19
    M1=M+1
    HA=0.0
    IX(1)=1

    IX(M1)=N
    
    L=(N-1)/M
    J=L
    DO 18 I=2,M
        IX(I)=J+1
        J=J+L
18  CONTINUE
28  HH=1.0
    DO 38 I=1,M1
        A(I)=Y(IX(I))
        H(I)=-HH
        HH=-HH
38  CONTINUE
    DO 59 J=1,M
    II=M1
    Y2=A(II)
    H2=H(II)
    DO 49 I=J,M
        D=X(IX(II))-X(IX(M1-I))
        Y1=A(M-I+J)
        H1=H(M-I+J)
        A(II)=(Y2-Y1)/D
        H(II)=(H2-H1)/D
        II=M-I+J
        Y2=Y1
        H2=H1
49  CONTINUE
59  CONTINUE
    HH=-A(M1)/H(M1)
    DO 68 I=1,M1
68  A(I)=A(I)+H(I)*HH
    DO 88 J=1,M-1
        II=M-J
        D=X(IX(II))
        Y2=A(II)
        DO 78 K=M1-J,M
            Y1=A(K)
            A(II)=Y2-D*Y1
            Y2=Y1
            II=K
78      CONTINUE
88  CONTINUE
    HM=ABS(HH)
    IF (HM.LE.HA) THEN
        A(M1)=-HM
        RETURN
    END IF
    A(M1)=HM
    HA=HM
    IM=IX(1)
    H1=HH
    J=1
    DO 105 I=1,N
        IF(I.EQ.IX(J)) THEN
            IF (J.LT.M1)J=J+1
        ELSE
        H2=A(M)
        DO 91 K=M-1,1,-1
91      H2=H2*X(I)+A(K)
        H2=H2-Y(I)
        IF (ABS(H2).GT.HM) THEN
            HM=ABS(H2)
            H1=H2
            IM=I
        END IF
    END IF
105 CONTINUE
    IF(IM.EQ.IX(1)) RETURN
    I=1
119 IF(IM.GE.IX(I)) THEN
        I=I+1
        IF (I.LE.M1) GOTO 119
    END IF
    IF (I.GT.M1) I=M1
    IF (I.EQ.(I/2)*2) THEN
        H2=HH
    ELSE
        H2=-HH
    END IF 
    IF(H1*H2.GE.0.0)THEN
        IX(I)=IM
        GOTO 28
    END IF
    IF(IM.LT.IX(1)) THEN
        DO 129 J=M,1,-1
129     IX(J+1)=IX(J)
        IX(1)=IM
        GOTO 28
    END IF
    IF(IM.GT.IX(M1)) THEN
        DO 139 J=2,M1
139     IX(J-1)=IX(J)
        IX(M1)=IM
        GOTO 28
    END IF
    IX(I-1)=IM
    GOTO 28
    END

#CCCCCCCCCC--median-based analysis of the segment slope--CCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCCCCCCCCCCCCCC----MBASS-----CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE MBASS(X,Y2,N,MC,INDEX,CG)
    INTEGER N,I,CG,J,J2,SA_SITE(1),PM1(1),N1,N2,INDEX,I1,U(400),U1(20)
    DOUBLE PRECISION X(N),Y2(N),MC,SLOPE(N,2),Y22(N),P(37,2),&
                     SA(N),SR(N),Wcrit,SW,XIG,Z(N),&
                     PM(20,2),PM2(1),Y(N)
    
    P(:,1)=(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
            21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 /)
    P(:,2)=(/ 12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,&
            2.228,2.201,2.179,2.160,2.145,2.131,2.120,2.110,2.101,&
            2.093,2.086,2.08,2.074,2.069,2.063,2.06,2.056,2.052,2.048,&
            2.045,2.042,2.021,2.009,2.003,1.994,1.990,1.987,1.984 /)
    
    U(:)=(/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
              0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,&
              0,0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,&
              0,0,0,0,1,2,3,4,4,5,6,7,8,9,10,11,11,12,13,13,&
              0,0,0,1,2,3,5,6,7,8,9,11,12,13,14,15,17,18,19,20,&
              0,0,1,2,3,5,6,8,10,11,13,14,16,17,19,21,22,24,25,27,&
              0,0,1,3,5,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,&
              0,0,2,4,6,8,10,13,15,17,19,22,24,26,29,31,34,36,38,41,&
              0,0,2,4,7,10,12,15,17,20,23,26,28,31,34,37,39,42,45,48,&
              0,0,3,5,8,11,14,17,20,23,26,29,33,36,39,42,45,48,52,55,&
              0,0,3,6,9,13,16,19,23,26,30,33,37,40,44,47,51,55,58,62,&
              0,1,4,7,11,14,18,22,26,29,33,37,41,45,49,53,57,61,65,69,&
              0,1,4,8,12,16,20,24,28,33,37,41,45,50,54,59,63,67,72,76,&
              0,1,5,9,13,17,22,26,31,36,40,45,50,55,59,64,67,74,78,83,&
             0,1,5,10,14,19,24,29,34,39,44,49,54,59,64,70,75,80,85,90,&
             0,1,6,11,15,21,26,31,37,42,47,53,59,61,70,75,81,86,92,98,&
            0,2,6,11,17,22,28,34,39,45,51,57,63,67,75,81,87,93,99,105,&
           0,2,7,12,18,24,30,36,42,48,55,61,67,74,80,86,93,99,106,112,&
          0,2,7,13,19,25,32,38,45,52,58,65,72,78,85,92,99,106,113,119,&
         0,2,8,13,20,27,34,41,48,55,62,69,76,83,90,98,105,112,119,127/)
    
    DO I=1,20
    U1(I)=I
    END DO
    
    CG=1
    
    DO I=1,N
    Y22(I)=LOG10(Y2(I))
    IF (Y2(I).LE.1) Y22(I)=0
    END DO
    
    DO I=1,N-1
    SLOPE(I,1)=(Y22(I+1)-Y22(I))/ABS(X(I+1)-X(I))
    END DO
    
    DO 40 I=1,N-2
    CALL CRANK(SLOPE(1:N-1+1-I,1),SLOPE(1:N-1+1-I,2),N-1+1-I)
    J2=0
    DO 10 J=1,N-1+1-I
    J2=J2+1
    SR(J2)=SUM(SLOPE(1:J,2))
    SA(J2)=ABS(2*SR(J2)-J*(N-1+1-I+1))
10  CONTINUE
    SA_SITE=MAXLOC(SA(1:J2))
    
    N1=SA_SITE(1)
    N2=N-1+1-I-N1
    Wcrit=N1*(N-1+1-I+1)/2
    SW=(N1*N2*(N-1+1-I+1)/12)**0.5
    
    IF (SR(SA_SITE(1)).LT.Wcrit) THEN
    XIG=0.5
    ELSE IF (SR(SA_SITE(1)).GT.Wcrit) THEN
    XIG=-0.5
    ELSE
    XIG=0.
    END IF
    
    Z(I)=(SR(SA_SITE(1))-Wcrit+XIG)/SW
    IF (Z(I).LE.-100) THEN
    CONTINUE
    END IF
    
    DO I1=1,20
    PM(I1,1)=ABS(U1(I1)-N1)
    PM(I1,2)=ABS(U1(I1)-N2)
    END DO
    PM1=MINLOC(PM(:,1))
    PM2=MINLOC(PM(:,2))
    IF (ABS(Z(I)).GE.U((PM1(1)-1)*20+PM2(1))) GOTO 20
    
40  CONTINUE
    CG=0
    GOTO 30
20  MC=X(N-I+1)
    INDEX=N-I+1
30  RETURN
    END          

#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCC----GFT-CCCCCCCCC
      SUBROUTINE GFT(X,Y,W,AQ,NQ,MQ,R,MC,GG,BV,BD,HMM,NN)
      DIMENSION X(NQ),Y(NQ),AQ(MQ+1),X2(NQ),Y2(NQ),YY22(NQ),&
               MII(100,3),WZ(1),RR(100),W(NQ),W2(NQ)
      INTEGER NQ,MQ,I,J,RF,GG,M1,NN,KI,KII,WZ,IK,RR,AYQ,R
      DOUBLE PRECISION X,Y,AQ,MC,X2,SUME,YY,Y2,HMM,YY22,MII,YY2,&
                      W,W2,AV,MA,BV,YZ,YYZ(100),BD
      
      
      J=0
      KI=0
      IK=0
      KII=0
      SUME=0.0
      YY=0.0
      MC=HMM
      GG=1
      M1Q=MQ+1
      
      GOTO 41
39    MC=MC+0.1
      J=0
      YY=0.0
      SUME=0.0
      IF (MC.GT.8.0) GOTO 12
      GOTO 41
12    IF (KI.EQ.0) GOTO 13
      WZ=MAXLOC(MII(1:KII,3))
      MC=MII(WZ(1),1)
      BV=MII(WZ(1),2)
      BD=MII(WZ(1),3)
      GOTO 36
13    GG=0
      GOTO 36
      
41    DO 31 I=1,NQ
      IF ((X(I)-MC).GE.-0.0001) GOTO 32
      GOTO 31
32    J=J+1
      X2(J)=X(I)
      Y2(J)=Y(I)
      W2(J)=W(I)
31    CONTINUE
      YM=J
      
      IF (J.EQ.0) GOTO 12
      IF (NN.GT.Y2(J)) GOTO 12
      
301   MA=SUM(W2(1:J))*1.0/Y2(J)*1.0
      BV=1.0/((MA-MC)*2.30258509)
      BD=BV/SQRT(Y2(J))
      
      YYZ(1:J)=EXP(-BV*X2(1:J)*LOG(10.0))
      YYZ(1:J)=YYZ(1:J)/MAXVAL(YYZ(1:J))*Y2(J)
      
      YY=0.
      YY2=0.
      DO I=1,J
      YY=YY+ABS(Y2(I)-YYZ(I))
      YY2=YY2+Y2(I)
      END DO
      
      RF=100-INT(YY/YY2*100)
      IK=IK+1
      RR(IK)=RF
      IF (RF.GE.R) GOTO 36
      IF (RF.GE.80) GOTO 49
      GOTO 39
49    KI=1
      KII=KII+1
      MII(KII,1)=MC
      MII(KII,2)=BV
      MII(KII,3)=BD
      
      GOTO 39
36    RETURN
      END
    
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCC----MBS----CCCCCCCCCC
      SUBROUTINE MBS(X,Y,Z,N,BV,BD,CC,MH,MC,HMM,NN,BIN,BY,JG)
      DIMENSION X(N),Y(N),X2(N),Y2(N),Z(N),Z2(N),&
               X3(N),Y3(N),Z3(N),BB(100,4),BB2(100,4),&
               WZ(1)
      DOUBLE PRECISION X,Y,MB,MC,Z,X2,Y2,Z2,MA,JG,&
         MC2,BV,BD,X3,Y3,Z3,BV2,BD2,MA2,DB,DB2,MM,BAVE,&
         HMM,BB,BB2,MC3,BY,CHUC(100,4)
      INTEGER(kind=4) N,I,K,K2,CC,NN,MH,KI,IIK,IIK2,WZ,MM2,&
                     BIN,I1,J,AZ(1)
      MC=HMM
      KI=0
      MC2=MC+JG
      K=0
      IIK=0
      IIK2=0
      K2=0
      BAVE=0
      CC=1
      MM=0.0
      J=1
      GOTO 70
      
60    IF (MC.GT.8.0) GOTO 80
      GOTO 90
80    IF ((MH.EQ.4).AND.(J.GT.1)) THEN
      AZ=MINLOC(CHUC(1:J-1,4))
      BV=CHUC(AZ(1),1)
      BD=CHUC(AZ(1),2)
      MC=CHUC(AZ(1),3)
      GOTO 50
      END IF
      CC=0
      GOTO 50
90    MC=MC+0.1
      MC2=MC+JG
      K=0
      K2=0
      BAVE=0
      MM=0
      
70    DO 10 I=1,N
      IF ((X(I)-MC).GE.-0.0001) GOTO 20
      GOTO 10
20    K=K+1
      X2(K)=X(I)
      Y2(K)=Y(I)
      Z2(K)=Z(I)
10    CONTINUE
      
      IF (K.EQ.0) GOTO 80
      IF (NN.GT.Y2(K)) GOTO 80
      
      MA=SUM(Z2(1:K))/Y2(K)*1.0
      BV=1.0/((MA-MC)*2.30258509)
      BD=BV/SQRT(Y2(K))
      
      DO 30 I=1,N
      IF ((X(I)-MC2).GE.-0.0001) GOTO 40
      GOTO 30
40    K2=K2+1
      X3(K2)=X(I)
      Y3(K2)=Y(I)
      Z3(K2)=Z(I)
30    CONTINUE
      
      IF (K2.EQ.0) GOTO 80
      IF (NN.GT.Y3(K2)) GOTO 80
      IF (MH.EQ.5) GOTO 121
      GOTO 110
      
      
121   DO I1=1,BIN
      MC3 = MC+I1*0.1
      K2 = 0
      
      DO 21 I=1,N
      IF ((X(I)-MC3).GE.-0.0001) GOTO 11
      GOTO 21
11    K2=K2+1
      X3(K2)=X(I)
      Y3(K2)=Y(I)
      Z3(K2)=Z(I)
21    CONTINUE
      
      IF (K2.EQ.0) GOTO 80
      IF (NN.GT.Y3(K2)) GOTO 80
      
      MA2=SUM(Z3(1:K2))/Y3(K2)
      BV2=1/((MA2-MC3)*2.30258509)
      BD2=BV2/SQRT(Y3(K2))
      BAVE=BAVE+BV2
      
      END DO
      GOTO 150
      
110   MA2=SUM(Z3(1:K2))/Y3(K2)
      BV2=(MA2-MC2)*2.30258509
      BV2=1.0/BV2
      BD2=BV2/SQRT(Y3(K2))
      
      IF (MH.EQ.4) GOTO 130
      
150   DO I=1,K
      MM=MM+ABS(X2(I)-MA)
      END DO
      
      DB=ABS(BAVE/BIN-BV)
      DB2=2.3*(BV**2)*SQRT(MM/(Y2(K)*(Y2(K)-1)))*&
         SQRT(Y2(K)*(Y2(K)-1))**0.6
      
      IF (DB2.GE.DB) GOTO 50
      GOTO 60
      
130   IF (ABS(BV2-BV).LE.BY) GOTO 50 
      IF (ABS(BV2-BV).LE.BY+0.08) THEN
      CHUC(J,1)=BV
      CHUC(J,2)=BD
      CHUC(J,3)=MC
      CHUC(J,4)=ABS(BV2-BV)
      J=J+1
      END IF
      GOTO 60
   
50    CONTINUE
      RETURN
      END
    
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCC-------EMR-------CCCCCCCCC
      SUBROUTINE EMR(X,Y,Z,W,A,N,M,CC,BV2,BD2,MC,HMM,NUM,P)
      DIMENSION X(N),Y(N),A(M+2),X2(N),Y2(N),Z(N),Z2(N),&
               DJJ(100),W(N),W2(N),MJJ(100),WPP(1),&
               BB(100,2),WZ(1),MCB(100),&
               EMRA(100,3),YZ(100),Z22(100),X22(100),QW(1)
      DOUBLE PRECISION X,Y,A,X2,MC,Y2,Z,Z2,DJ,DJJ,MC2,&
         MA2,BV2,BD2,W,W2,MJJ,ABZ,HMM,BB,MCB,&
         AV,MIU,XGM1,XGM,YZ,Z22,X22,MGASS,ZDZ,XS,&
         SZ2(100,4),D_REAL(100),D_CRIT,S,A_X,B_X,ZX2,CA(50,2),&
         MAX_IN,GL(100),P,KS_P
      INTEGER N,M,I,I2,CC,M1,WPP,WZ,I1,QW,&
             NUM,IJ,SZ1(100),TOP
      EXTERNAL MGASS
      
      MC=HMM+0.1
      I2=0
      I1=0
      CC=1
      M1=M+1
      DJJ=0.0
      J=1
      IJ=0
      GOTO 300
      
200   MC=MC+0.1
      I2=0
      I1=0
      A=0
      CA=0
      J=J+1
      X2=0
      Y2=0
      Z2=0
      W2=0
      
300   DO 10 I=1,N
      IF ((X(I)-MC).GE.0.0001) GOTO 20
      I1=I1+1
      X22(I1)=X(I)
      
      GOTO 10
20    I2=I2+1
      X2(I2)=X(I)
      Y2(I2)=Y(I)
      Z2(I2)=Z(I)
      W2(I2)=W(I)
      IF (Y(I).LT.1) Y2(I2)=0
10    CONTINUE 
      
      IF (I2.EQ.0) GOTO 280
      GOTO 500

#C     
500   MC2=MC
      MA2=SUM(W2(1:I2))*1.0/Y2(I2)*1.0
      BV2=(MA2-MC2)*2.30258509
      BV2=1.0/BV2*1.0
      BD2=BV2*1.0/SQRT(Y(I2))*1.0
      MCB(J)=MC
      BB(J,1)=BV2
      BB(J,2)=BD2
      
      DO I=1,I2
      Z22(I)=EXP(-BV2*X2(I)*2.3)
      IF (Z22(I).LT.0) Z22(I)=0.
      END DO
      Z22(1:I2)=Z22(1:I2)/MAXVAL(Z22(1:I2))*Z(I2)
      
#C   
      MIU=0.
      DO 55 I=1,N
55    MIU=MIU+X(I)/N
      
      XGM1=0.
      DO 56 I=1,N
56    XGM1=XGM1+(X(I)-MIU)**2/N
      XGM=SQRT(XGM1)*1.0
      
      A_X=-100
      DO I=1,I1
      B_X=X22(I)
      CALL CDFX(A_X,B_X,10000,XGM,MIU,S)
      YZ(I)=S
      END DO
      
      CALL CDFX(A_X,X2(I2),10000,XGM,MIU,S)
      
      YZ(1:I1)=YZ(1:I1)*Z(I2)/S

      DO 59 I=1,I1
59    Z22(I2+I)=YZ(I)
      
#C     
      IF (J.EQ.5) THEN
      OPEN (533,FILE='533.txt')
      DO I=1,N
      WRITE (533,'(F5.2,1X,F10.3)') X(I),Z22(I)
      END DO
      CLOSE(533)
      OPEN (534,FILE='534.txt')
      DO I=1,N
      WRITE (534,'(F5.2,1X,F10.3)') X(I),Z(I)
      END DO
      CLOSE(534)
      END IF
      
#C-----
      TOP=INT(MAX(MAXVAL(Z22(1:N)),MAXVAL(Z(1:N))))+1
      CALL CEF(Z(1:N),CA(1:50,1),N,TOP,50)
      CALL CEF(Z22(1:N),CA(1:50,2),N,TOP,50)
      MAX_IN=MAXVAL(ABS(CA(:,1)-CA(:,2)))
      
      GL(J)=MAX_IN
      KS_P=SQRT(-LOG(P/2.0)*0.5)*SQRT(2.0/N)
#C      IF (MAX_IN.LT.KS_P) GOTO 100 
      GOTO 200
      
280   IF (J.NE.1) THEN
      WZ=MINLOC(GL(1:J))
      BV2=BB(WZ(1),1)
      BD2=BB(WZ(1),2)
      MC=MCB(WZ(1))
      ELSE
      CC=0    
      END IF
100   RETURN 
      END