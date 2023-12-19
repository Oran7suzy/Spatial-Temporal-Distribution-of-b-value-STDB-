C     Spatial Distribution of b-value (STDB)
C
C     A versatile b-value calculator. Used for calculating time and space variations of b-values, 
C              it includes eight Mc calculation methods and four b-value calculation methods. 
C              Additional functions include synthetic seismic catalogs for seismic events, and more.
C
C     CSTDB is an extremely simple b-value calculator that includes most algorithms on the market, and the 
C             grouping between different algorithms is flexible and efficient. Just import a seismic catalog.
C
C     Refer to:
C     Weicheng Gong, Huayuan Cheng, Yajing Gao, Qing Li, and Yunqiang Sun (2023). Spatial-Temporal Distribution of b-value in the Eastern Tibetan Plateau. (in revision)
C      
C     INPUT
C     Catalog Format (the following numbers represent column numbers):
C     1.Year 2.Month 3.Day 4.Hour 5.Minute 6.Latitude 7.Longitude 8.Magnitude
C
C     OUTPUT
C     b-value as a function of time: TIME_B.txt
C     Mc as a function of time: TIME_Mc.txt
C     Test the results of the artificial earthquake catalog: CYC_B.TXT; CYC_MC.TXT
      
      DIMENSION TIME(1000000,5),NUMER11(1000000,3),
     1    TIME1(1000000,5),NUMER2(1000000,3),
     2    CATLOG(1000000,3),TIME2(5),FREQUENCY(100,2),
     3    CUMULATIVE(100,4),BZHU(5000,7,6),INDEXX(1),BV(5000,2,6),
     4    X(5000),Y(5000),Z(5000),Z2(5000),MMC(5000,6),TIME_FW1(5),
     5    TIME_FW2(5),TIME22(1000000,5),NUMER22(1000000,3),BBT(2),SS1(2)
     6    ,SS(5000),YA(5000),AZD(5000,6),XX1(5000),YYQ(5000),ZZQ(5000),
     7    ZZQ2(5000),A(10)
      INTEGER COUNT,TIME,TIME1,TIME2,WQ,WWQ,YB,IPP,OP,QQ,QQ1,CCON,
     1        INDEXX2,BZHU,EW,A1,A2,II,II3,II4,TEMP,JJ,YBB,SZ,
     2        SJ,OPP,CM,IPP2,OP2,LP,AP,user_variable,CM2
     3        ,FW,I,FW2,PO,MH,M,M1,CG,GG,R,CC,YZZ,GC,N,OOP(6),KK2
     4        ,TIME_FW1,TIME_FW2,CCN,TIME22,IOP,JJK,KK,KK1,J,SUI,ML,QWS
     5        ,II2,WZ1,TEMP1,BIN,num_vertices,crossings,CHUSS,CHUS,WZ(1)
     6        ,I2,I3,CYC,CYC_T,CYC_D,CYC_W,RMS,TEMPP,NUME(10000,6),
     7        LIMI,APN(6),CYC2,JK,XS1,CS,J1,J4,TIME3(1000000,5),B_MO,CS2
      DOUBLE PRECISION NUMER11,BV,B_VALE,B_DOUBLE,NUMER2,
     1        CATLOG,FREQUENCY,CUMULATIVE,MAXN1,MINN1,PP,M_AVERAGE,
     2        MC,BBTIME,X,Y,Z,Z2,A,MCC,BEE,HMM,MMC,NUMER22
     3        ,AB1,AB2,BBT,ABT,SS,SS1,S,AAW,AA1,BB1,YA,AZD,EMA,MGASS,
     4        MERRF,MGAM1,MGAM2,MIU,XIGMA,MAG_MIN,MAG_MAX,PP2,BB2,AA2
     5        ,QWN,QW,QQW,XX1,YYQ,ZZQ,ZZQ2,A_X,B_X,S_CDF,Y2(5000),
     6        vertices_x(100),vertices_y(100),MI_MAG,RM,RP,XX2(500000,2)
     7        ,XX3(500,2),GETMAX,X22(5000),Y22(5000),SUM_NU,SUM_XX2,BY,
     8        MJ,XS,XS2,JG,P,INTVAL,F3,F1,F2,P2,AA
      EXTERNAL MGASS,MERRF,MGAM1,MGAM2
      CHARACTER*132 stamp
      CHARACTER*40 CATA,CATA2
      COMMON/kuser/user_variable(1)
      LOGICAL isUnique


      COUNT=0
      CCON=0
      IOP=0
      OPP=0
      CCN=0
      LP=0
      OOP=0
      user_variable(1)=0
      
     
      
      WRITE (*,*) 'Seismic catalogue: 1-Artificial(OK1993) 
     1 2-Artificial(WW2005) 3-Real'
      READ (*,*) ML
      
      WRITE (*,*) 'b_value--Algorithm to choose: 1-MLE 2-LSM 3-BP'
      READ (*,*) B_MO
      
      WRITE (*,*) 'Mc--Algorithm to choose:  1-EMR  2-MAXC  3-GFT  4-MBS
     1  5-MBS-JS  6-MBASS  7-MAH  8-NDT'
      READ (*,*) MH
      
      IF (ML.EQ.3) THEN
      WRITE (*,*) 'Enter catalog filename (DATA.TXT):'
      READ (*,'(a)') CATA
      GOTO 546
      END IF
      
      WRITE (*,*) 'μ=   σ=   b=    !μ<=M_max'
      READ (*,*) MIU,XIGMA,BB2
      WRITE (*,*) 'Earthquake magnitude range:  Min=  Max='
      READ (*,*) MAG_MIN,MAG_MAX
      WRITE (*,*) 'Minimum seismic data：'
      READ (*,*) LIMI
      WRITE (*,*) 'Whether to cycle test: 1-yes 2-no'
      READ (*,*) CYC
      IF (CYC.EQ.1) THEN
      WRITE (*,*) 'MAXIMUM TOTAL:'
      READ (*,*) CYC_T
      WRITE (*,*) 'MINIMUM TOTAL:'
      READ (*,*) CYC_D
      WRITE (*,*) 'STEP:'
      READ (*,*) CYC_W
      ELSE
      WRITE (*,*) 'Total number of catalog: '
      READ (*,*) QWS
      END IF
      WRITE (*,*) 'REJECTION NUMBER: '
      READ (*,*) RMS
      WRITE (*,*) 'Whether to test all Mc methods: 1-yes 2-no '
      READ (*,*) CYC2
      IF (CYC2.EQ.1) THEN
      WRITE (*,*) 'MBS: Box width (0.5→5)：'
      READ (*,*) BIN
      WRITE (*,*) 'MBS: MAGNITUDE INTERVAL (0.1)：'
      READ (*,*) JG
      WRITE (*,*) 'MBS-JS: Stable platform (0.03-0.05)：'
      READ (*,*) BY
      END IF
      GOTO 963
      
546   CONTINUE
      
65    IF (MH.EQ.5) THEN
      WRITE (*,*) 'Box width (0.5→5)：'
      READ (*,*) BIN
      ELSE IF (MH.EQ.4) THEN
      WRITE (*,*) 'Stable platform (0.03-0.05)：'
      READ (*,*) BY
      WRITE (*,*) 'MAGNITUDE INTERVAL (0.1)：'
      READ (*,*) JG
      ELSE IF (MH.EQ.1) THEN
      WRITE (*,*) 'Significance of KS test (default:0.05)：'
      READ (*,*) P
      ELSE IF (MH.EQ.3) THEN
      WRITE (*,*) 'Goodness of fit (90%→90)：'
      READ (*,*) R
      ELSE IF (MH.EQ.7) THEN
      WRITE (*,*) 'Differential Convergence Criteria (0.005-0.0001)：'
      READ (*,*) P
      ELSE IF (MH.EQ.8) THEN
      WRITE (*,*) 'Number of bootstrap samples (default:50)：'
      READ (*,*) CS2
      END IF
      
      
      
      WRITE (*,*) 'Each Window minimum seismic data：'
      READ (*,*) LIMI

244   WRITE (*,*) 'Calculation model:  1-TIME  2-NUMER  3-ENTIRE AREA'
      READ (*,*) CM
      
      IF (CM.EQ.1) GOTO 223
      IF (CM.EQ.3) GOTO 545
      WRITE (*,*) 'The number of earthquakes within one time step:'
      READ (*,*) SZ
      IF (SZ.GE.LIMI) GOTO 535
      WRITE (*,*) 'Input the number does not meet the minimum number 
     1 within the grid'
      PAUSE
      STOP
535   GOTO 545
223   WRITE (*,*) 'Time step (s):  '
      READ (*,*) CM2
      
      WRITE (*,*) 'MIN MAGNITUDE:'
      READ (*,*) MI_MAG
      
545   WRITE (*,*) 'Input earthquake catalog range:   1-YES   2-NO'
      READ (*,*) FW
      IF (FW.EQ.1) GOTO 665
      GOTO 667
665   WRITE (*,*) '1-LMSF 2-XSHF 3-KLF 4-Tibetan Plateau South
     1    5-LMSF AND XSHF 6-OTHER '
      READ (*,*) FW2
      IF (FW2.NE.6) GOTO 667
      
      WRITE (*,*) 'Input area range file (points clockwise):'
      READ (*,'(a)') CATA2
      
      num_vertices=0
      OPEN(115,FILE=CATA2)
      DO WHILE(NOT(EOF(115)))
      num_vertices=num_vertices+1
      READ(115,*) vertices_x(num_vertices),vertices_y(num_vertices)
      END DO
      CLOSE(115)
      
667   IF (CM.EQ.3) GOTO 666
      GOTO 450
666   WRITE (*,*) 'Start time:  YEAR；MOUTH；DAY；HOUR；MINUTE'
      READ (*,*) TIME_FW1(1),TIME_FW1(2),TIME_FW1(3),TIME_FW1(4)
     1,TIME_FW1(5)
      WRITE (*,*) 'Deadline:  YEAR；MOUTH；DAY；HOUR；MINUTE'
      READ (*,*) TIME_FW2(1),TIME_FW2(2),TIME_FW2(3),TIME_FW2(4)
     1,TIME_FW2(5)
      
450   OPEN(101,FILE=CATA)
      DO WHILE(NOT(EOF(101)))
      COUNT=COUNT+1
      READ(101,*) (TIME(COUNT,I),I=1,5),NUMER11(COUNT,2),
     1                NUMER11(COUNT,1),NUMER11(COUNT,3)
      END DO
      CLOSE(101)
      
      JK=1
      
      IF (CM.EQ.3) GOTO 898
      GOTO 896
898   DO 103 K=1,COUNT
      IF (TIME_FW2(1).EQ.TIME_FW1(1)) GOTO 104
      IF ((TIME(K,1).LT.TIME_FW2(1)).AND.
     1(TIME(K,1).GT.TIME_FW1(1))) GOTO 60
      IF (TIME(K,1).EQ.TIME_FW2(1)) GOTO 20
      IF (TIME(K,1).EQ.TIME_FW1(1)) GOTO 70
      GOTO 103
      
104   IF (TIME(K,1).EQ.TIME_FW1(1)) GOTO 114
      GOTO 103
114   IF ((TIME(K,2).LT.TIME_FW2(2)).AND.
     1(TIME(K,2).GT.TIME_FW1(2))) GOTO 60
      IF (TIME(K,2).EQ.TIME_FW2(2)) GOTO 105
      IF (TIME(K,2).EQ.TIME_FW1(2)) GOTO 106
      
      GOTO 103
      
105   IF (TIME(K,3).LE.TIME_FW2(3)) GOTO 107    
      GOTO 103
107   IF (TIME(K,4).LE.TIME_FW2(4)) GOTO 108
      GOTO 103
108   IF (TIME(K,5).LE.TIME_FW2(5)) GOTO 60
      GOTO 103

106   IF (TIME(K,3).GE.TIME_FW1(3)) GOTO 109    
      GOTO 103
109   IF (TIME(K,4).GE.TIME_FW1(4)) GOTO 119
      GOTO 103
119   IF (TIME(K,5).GE.TIME_FW1(5)) GOTO 60
      GOTO 103
      
70    IF (TIME(K,2).GE.TIME_FW1(2)) GOTO 80 
      GOTO 103
80    IF (TIME(K,3).GE.TIME_FW1(3)) GOTO 90
      GOTO 103
90    IF (TIME(K,4).GE.TIME_FW1(4)) GOTO 110
      GOTO 103
110   IF (TIME(K,5).GE.TIME_FW1(5)) GOTO 60
      GOTO 103
      
20    IF (TIME(K,2).LE.TIME_FW2(2)) GOTO 30    
      GOTO 103
30    IF (TIME(K,3).LE.TIME_FW2(3)) GOTO 40
      GOTO 103
40    IF (TIME(K,4).LE.TIME_FW2(4)) GOTO 50
      GOTO 103
50    IF (TIME(K,5).LE.TIME_FW2(5)) GOTO 60
      GOTO 103
      
60    CCN=CCN+1
      NUMER22(CCN,1)=NUMER11(K,1)
      NUMER22(CCN,2)=NUMER11(K,2)
      NUMER22(CCN,3)=NUMER11(K,3)
      TIME22(CCN,1)=TIME(K,1)
      TIME22(CCN,2)=TIME(K,2)
      TIME22(CCN,3)=TIME(K,3)
      TIME22(CCN,4)=TIME(K,4)
      TIME22(CCN,5)=TIME(K,5)
103   CONTINUE
      
      IF (CM.EQ.3) GOTO 897
      GOTO 896
897   COUNT=CCN
      DO 895 K=1,COUNT
      NUMER11(K,1)=NUMER22(K,1)
      NUMER11(K,2)=NUMER22(K,2)
      NUMER11(K,3)=NUMER22(K,3)
      TIME(K,1)=TIME22(K,1)
      TIME(K,2)=TIME22(K,2)
      TIME(K,3)=TIME22(K,3)
      TIME(K,4)=TIME22(K,4)
      TIME(K,5)=TIME22(K,5)
895   CONTINUE   
      
      
896   IF (FW.EQ.1) GOTO 337
      GOTO 335
337   IF (FW2.EQ.1) GOTO 338
      IF (FW2.EQ.2) GOTO 339
      IF (FW2.EQ.3) GOTO 340
      IF (FW2.EQ.4) GOTO 341
      IF (FW2.EQ.5) GOTO 342
      GOTO 336
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCC----
338   vertices_x(1)=106.728
      vertices_y(1)=33.383
      vertices_x(2)=105.8
      vertices_y(2)=33.2
      vertices_x(3)=103.8
      vertices_y(3)=32.05
      vertices_x(4)=102.6
      vertices_y(4)=30.8
      vertices_x(5)=102.2
      vertices_y(5)=30.011
      vertices_x(6)=102.181
      vertices_y(6)=29.45958
      vertices_x(7)=102.7
      vertices_y(7)=29.4
      vertices_x(8)=103.325
      vertices_y(8)=30.212
      vertices_x(9)=104.05
      vertices_y(9)=30.95
      vertices_x(10)=105.083
      vertices_y(10)=31.728
      vertices_x(11)=106.291
      vertices_y(11)=32.166
      vertices_x(12)=107.046
      vertices_y(12)=32.856
      
      num_vertices=12
      GOTO 336
      
CCCCCCCCCCCCCCC------  
339   vertices_x(1)=99.01 
      vertices_y(1)=32.65
      vertices_x(2)=100.0254 
      vertices_y(2)=32.16046
      vertices_x(3)=101.65 
      vertices_y(3)=31.05608
      vertices_x(4)=102.4004 
      vertices_y(4)=29.84923
      vertices_x(5)=102.6827 
      vertices_y(5)=29.15188
      vertices_x(6)=102.572 
      vertices_y(6)=28.36127
      vertices_x(7)=102.5542 
      vertices_y(7)=28.02815
      vertices_x(8)=101.9042 
      vertices_y(8)=28.02815
      vertices_x(9)=101.822 
      vertices_y(9)=28.36127
      vertices_x(10)=101.8327 
      vertices_y(10)=29.15188
      vertices_x(11)=101.6504 
      vertices_y(11)=29.84923
      vertices_x(12)=100.72 
      vertices_y(12)=30.92608
      vertices_x(13)=100.0254 
      vertices_y(13)=31.45046
      vertices_x(14)=99.01 
      vertices_y(14)=31.88
     
      num_vertices=14
      GOTO 336
      
CCCCCCCCCCCCCCC---   
340   vertices_x(1)=99.21422 
      vertices_y(1)=35.29692
      vertices_x(2)=99.37691 
      vertices_y(2)=35.17615
      vertices_x(3)=99.50872 
      vertices_y(3)=35.10263
      vertices_x(4)=100.502 
      vertices_y(4)=34.64884
      vertices_x(5)=101.8791 
      vertices_y(5)=34.44755
      vertices_x(6)=102.0393 
      vertices_y(6)=34.37779
      vertices_x(7)=102.4625 
      vertices_y(7)=34.22925
      vertices_x(8)=103.11 
      vertices_y(8)=34.16
      vertices_x(9)=104.41 
      vertices_y(9)=33.41
      vertices_x(10)=105.41 
      vertices_y(10)=33.51
      vertices_x(11)=105.41 
      vertices_y(11)=32.81
      vertices_x(12)=104.41 
      vertices_y(12)=32.61
      vertices_x(13)=103.11 
      vertices_y(13)=33.46
      vertices_x(14)=102.4625 
      vertices_y(14)=33.52925
      vertices_x(15)=102.0393 
      vertices_y(15)=33.67779
      vertices_x(16)=101.8791 
      vertices_y(16)=33.74755
      vertices_x(17)=100.502 
      vertices_y(17)=33.94884
      vertices_x(18)=99.50872 
      vertices_y(18)=34.40263
      vertices_x(19)=99.37691 
      vertices_y(19)=34.47615
      vertices_x(20)=99.21422 
      vertices_y(20)=34.59692
      
      num_vertices=20
      GOTO 336
CCCCCCCCCCCCCCC------
341   vertices_x(1)=100.01
      vertices_y(1)=28.01
      vertices_x(2)=108.01
      vertices_y(2)=28.01
      vertices_x(3)=108.01
      vertices_y(3)=35.01
      vertices_x(4)=100.01
      vertices_y(4)=35.01
      
      num_vertices=4
      GOTO 336
CCCCCCCCCCCCCCC---
342   vertices_x(1)=102.01
      vertices_y(1)=28.801
      vertices_x(2)=101.01
      vertices_y(2)=30.101
      vertices_x(3)=105.701
      vertices_y(3)=33.801
      vertices_x(4)=106.701
      vertices_y(4)=32.601
      num_vertices=4
      
336   WQ=0

      DO EW=1,COUNT
      CHUSS=0
      DO I=1,num_vertices
      J=I+1
      IF (I.EQ.num_vertices) J=1
      CALL INTERSECTION(vertices_x(I),vertices_y(I),
     1                    vertices_x(J),vertices_y(J),
     2                    NUMER11(EW,1),NUMER11(EW,2),CHUS)
      IF (CHUS.EQ.1) CHUSS = CHUSS+ 1
      END DO

      IF ((mod(CHUSS, 2)==1).AND.(NUMER11(EW,3).GE.MI_MAG)) then
      WQ=WQ+1
      NUMER2(WQ,1)=NUMER11(EW,1)
      NUMER2(WQ,2)=NUMER11(EW,2)
      NUMER2(WQ,3)=NUMER11(EW,3)
      TIME1(WQ,1)=TIME(EW,1)
      TIME1(WQ,2)=TIME(EW,2)
      TIME1(WQ,3)=TIME(EW,3)
      TIME1(WQ,4)=TIME(EW,4)
      TIME1(WQ,5)=TIME(EW,5)
      END IF
      END DO
      
      GOTO 556
      
C------
335   DO 995 I=1,COUNT
      NUMER2(I,1)=NUMER11(I,1)
      NUMER2(I,2)=NUMER11(I,2)
      NUMER2(I,3)=NUMER11(I,3)
      TIME1(I,1)=TIME(I,1)
      TIME1(I,2)=TIME(I,2)
      TIME1(I,3)=TIME(I,3)
      TIME1(I,4)=TIME(I,4)
      TIME1(I,5)=TIME(I,5)
995   CONTINUE
      WQ=COUNT

C------
556   OPEN (124,FILE='SJ.txt')
      DO 134 EW=1,WQ
      WRITE (124,'(2(F6.2,1X),F3.1)') (NUMER2(EW,I),I=1,3)
134   CONTINUE
      CLOSE (124)
      
      OPEN (125,FILE='SJ2.txt')
      DO 135 EW=1,WQ
      BBTIME=(TIME1(EW,1)+((TIME1(EW,2)-1)*43200+TIME1(EW,3)*1440+TIME1
     1        (EW,4)*60+TIME1(EW,5))*1.0/525600.0)*1.0
      WRITE (125,'(F8.3,3X,F3.1)') BBTIME,NUMER2(EW,3)
135   CONTINUE
      CLOSE (125)
            
C------
      IF (CM.EQ.1) GOTO 990
      IF (CM.EQ.2) GOTO 880
      IF (CM.EQ.3) GOTO 870
      
CCCCCCCCCC---
870   DO 860 K=1,WQ
      IOP=IOP+1
      CATLOG(IOP,1)=NUMER2(IOP,1)
      CATLOG(IOP,2)=NUMER2(IOP,2)    
      CATLOG(IOP,3)=NUMER2(IOP,3)
      TIME3(IOP,1)=TIME1(IOP,1)
      TIME3(IOP,2)=TIME1(IOP,2)
      TIME3(IOP,3)=TIME1(IOP,3)
      TIME3(IOP,4)=TIME1(IOP,4)
      TIME3(IOP,5)=TIME1(IOP,5)
860   CONTINUE
      MAXN1=MAXval(CATLOG(1:IOP,3))
      MINN1=MINval(CATLOG(1:IOP,3))
      SZ=IOP
      AP=1
      TIME2(1)=TIME1(IOP,1)
      TIME2(2)=TIME1(IOP,2)
      TIME2(3)=TIME1(IOP,3)
      TIME2(4)=TIME1(IOP,4)
      TIME2(5)=TIME1(IOP,5)
      GOTO 360
      
CCCCCCCCCC----
880   WWQ=FLOOR(WQ*1.0/SZ*1.0)
      YYB=0
      IF ((WWQ*SZ-WQ).NE.0) THEN
      WWQ=WWQ+1
      YYB=ABS(WWQ*SZ-WQ)
      END IF
      
      CALL TIME_sorting(TIME1(:WQ,:5),NUMER2(:WQ,:3),WQ,5,3) 
      SZ2=SZ
      
      DO 100 IPP=1,WWQ
      IF ((YYB.NE.0).AND.(IPP.EQ.WWQ)) SZ=YYB
      DO 200 OP=1,SZ
      CATLOG(OP,1)=NUMER2(OP+(IPP-1)*SZ2,1)
      CATLOG(OP,2)=NUMER2(OP+(IPP-1)*SZ2,2)    
      CATLOG(OP,3)=NUMER2(OP+(IPP-1)*SZ2,3)    
      TIME3(OP,1)=TIME1(OP+(IPP-1)*SZ2,1)
      TIME3(OP,2)=TIME1(OP+(IPP-1)*SZ2,2)
      TIME3(OP,3)=TIME1(OP+(IPP-1)*SZ2,3)
      TIME3(OP,4)=TIME1(OP+(IPP-1)*SZ2,4)
      TIME3(OP,5)=TIME1(OP+(IPP-1)*SZ2,5)
200   CONTINUE
      MAXN1=MAXval(CATLOG(1:SZ,3))
      MINN1=MINval(CATLOG(1:SZ,3))
      CALL LOOKF(TIME3(1:SZ,:5),SZ,5,1,J4)
      IF (IPP.EQ.WWQ) THEN
      CONTINUE
      END IF
      TIME2(1)=TIME3(J4,1)
      TIME2(2)=TIME3(J4,2)
      TIME2(3)=TIME3(J4,3)
      TIME2(4)=TIME3(J4,4)
      TIME2(5)=TIME3(J4,5)
      GOTO 363
      
C-------
990   CONTINUE
      J3=0

      CALL TIME_sorting(TIME1(:WQ,:5),NUMER2(:WQ,:3),WQ,5,3) 
      CALL LOOKF(TIME1(1:WQ,:5),WQ,5,1,J)
      CALL LOOKF(TIME1(1:WQ,:5),WQ,5,2,J1)
      
      SJ=FLOOR(((TIME1(J,1)-TIME1(J1,1))*31536000+
     1    (TIME1(J,2)-TIME1(J1,2))*2628000+
     2    (TIME1(J,3)-TIME1(J1,3))*86400+
     3    (TIME1(J,4)-TIME1(J1,4))*3600+
     4    (TIME1(J,5)-TIME1(J1,5))*60)*1.0/CM2*1.0)+1
      
889   DO 100 IPP2=1,SJ
      DO 120 OP2=1,WQ
      F1=TIME1(OP2,1)+TIME1(OP2,2)*1.0/12.0+
     2   TIME1(OP2,3)*1.0/12.0/30.0+
     3   TIME1(OP2,4)*1.0/12.0/30.0/24.0+
     4   TIME1(OP2,5)*1.0/12.0/30.0/24.0/60.0
      F2=TIME1(J1,1)+TIME1(J1,2)*1.0/12.0+
     2   TIME1(J1,3)*1.0/12.0/30.0+
     3   TIME1(J1,4)*1.0/12.0/30.0/24.0+
     4   TIME1(J1,5)*1.0/12.0/30.0/24.0/60.0+IPP2*CM2*1.0/31536000.0
      F3=TIME1(J1,1)+TIME1(J1,2)*1.0/12.0+
     2   TIME1(J1,3)*1.0/12.0/30.0+
     3   TIME1(J1,4)*1.0/12.0/30.0/24.0+
     4   TIME1(J1,5)*1.0/12.0/30.0/24.0/60.0+(IPP2-1)*CM2*1.0/31536000.0
      IF ((F1.LE.F2).AND.(F1.GE.F3-0.0000001)) GOTO 130
      GOTO 120
130   OPP=OPP+1
      CATLOG(OPP,1)=NUMER2(OP2,1)
      CATLOG(OPP,2)=NUMER2(OP2,2)    
      CATLOG(OPP,3)=NUMER2(OP2,3) 
      TIME3(OPP,1)=TIME1(OP2,1)
      TIME3(OPP,2)=TIME1(OP2,2)
      TIME3(OPP,3)=TIME1(OP2,3)
      TIME3(OPP,4)=TIME1(OP2,4)
      TIME3(OPP,5)=TIME1(OP2,5)
120   CONTINUE
      
      IF (OPP.LT.LIMI) GOTO 100
      J3=J3+1
      
      user_variable(1)=user_variable(1)+OPP
      LP=OPP
      MAXN1=MAXval(CATLOG(1:LP,3))
      MINN1=MINval(CATLOG(1:LP,3))
      CALL LOOKF(TIME3(1:OPP,:5),OPP,5,1,J4)
      TIME2(1)=TIME3(J4,1)
      TIME2(2)=TIME3(J4,2)
      TIME2(3)=TIME3(J4,3)
      TIME2(4)=TIME3(J4,4)
      TIME2(5)=TIME3(J4,5)
      OPP=0

      SZ=LP
      AP=J3-OOP(JK)
      GOTO 360
363   AP=IPP-OOP(JK)

360   IF (B_MO.EQ.3) GOTO 412
      IF (SZ.EQ.0) GOTO 350
C--
      CALL AR_AR(CATLOG(:SZ,3),SZ,CUMULATIVE,TEMP,MAXN1)
      
      OPEN (156,FILE='FMD_single.txt')
      OPEN (155,FILE='FMD_cumulate.txt')
      DO 558 PO=1,TEMP
      X(PO)=CUMULATIVE(PO,1)  
      Y(PO)=CUMULATIVE(PO,4)  
      Z(PO)=CUMULATIVE(PO,2)  
      Z2(PO)=CUMULATIVE(PO,3) 
      
      AA1=LOG10(CUMULATIVE(PO,2)+1)
      BB1=LOG10(CUMULATIVE(PO,4)+1) 
      WRITE (155,'(F5.1,1X,F10.4,1X)') CUMULATIVE(PO,1),BB1
      WRITE (156,'(F5.1,1X,F10.4,1X)') CUMULATIVE(PO,1),AA1
558   CONTINUE
      N=TEMP
      CLOSE (155)
      CLOSE (156)
      
      QWS=INT(Y(TEMP))
      IF (ML.EQ.3) GOTO 559
      

963   TEMP=0
      DO PP2=MAG_MAX,MAG_MIN-0.1,-0.1
      TEMP=TEMP+1
      X(TEMP)=PP2
      END DO
      
      CALL ART_EC(X,Y,Y2,TEMP,MAG_MAX,MAG_MIN,XIGMA,MIU,BB2,ML)
      
C     Prepare for multiple cycle sampling use
      TEMPP=TEMP
      Y22(1:TEMPP)=Y2(1:TEMPP)
      X22(1:TEMPP)=X(1:TEMPP)
      AP=1
      
C     Cycle sampling star
      IF (CYC.EQ.2) GOTO 551
      DO 552 II=1,CYC_W
      QWS=INT(CYC_D+(II-1)*(CYC_T-CYC_D)/CYC_W)

551   I2=0
      SUM_XX2=0.
      CALL RANDOM_SEED()  ! Ensure that the number of data is not repeated
      DO 517 I=1,RMS
      IF (SUM_XX2.GE.QWS) EXIT
      CALL RANDOM_NUMBER(RM)  ! Generate random magnitude
      CALL RANDOM_NUMBER(RP)  ! Generate random probability
      RM = RM*(MAG_MAX-MAG_MIN)
      RP = RP*(MAXVAL(Y22(1:TEMPP))-MINVAL(Y22(1:TEMPP)))
      WZ = MINLOC(ABS(X22(1:TEMPP)-RM))
      IF (Y22(WZ(1)).GE.RP) THEN
      I2=I2+1
      XX2(I2,1)=X22(WZ(1))
      IF (RP.LE.0.00001) RP=0.
      XX2(I2,2)=RP
      SUM_XX2=SUM_XX2+RP*QWS
      END IF
517   CONTINUE

      I3=0
      XX3=0.0 
      DO I=1,I2
      isUnique=ALL(XX2(I,1)/=XX3(:,1))
      IF(isUnique) THEN
      I3=I3+1
      XX3(I3,1)=XX2(I,1)
      CALL SORTSUM(XX2(1:I2,:),I2,2,1,2,XX3(I3,1),SUM_NU)
      XX3(I3,2)=SUM_NU
      END IF
      END DO
      CALL SORTNUM(XX3(1:I3,:),I3,2,1)
      
      X(1:I3)=XX3(1:I3,1)
      Y2(1:I3)=XX3(1:I3,2)*QWS
      DO I=1,I3
      Y(I)=SUM(Y2(1:I))
      END DO 
      Z(1:I3)=XX3(1:I3,2)*QWS
      Z2(1:I3)=Z(1:I3)*X(1:I3)
      
      TEMP=I3
      N=TEMP
      
C       
      TEMP1=0
      DO PP2=MAG_MAX,MAG_MIN-0.1,-0.01
      TEMP1=TEMP1+1
      XX1(TEMP1)=PP2
      YYQ(TEMP1)=1/(SQRT(2*3.1415)*XIGMA)
     1*EXP((-(XX1(TEMP1)-MIU)**2)/(2*XIGMA**2))
      END DO
      
      DO I=1,TEMP1
      IF (I.EQ.1) THEN
      ZZQ(I)=YYQ(I)
      ELSE
      ZZQ(I)=ABS(YYQ(I)+YYQ(I-1))/2*0.01+ZZQ(I-1)
      END IF
      END DO
      
      DO I=1,TEMP1
      IF (I.EQ.1) THEN
      ZZQ2(I)=ZZQ(I)
      ELSE
      ZZQ2(I)=ABS(ZZQ(I)+ZZQ(I-1))/2*0.01+ZZQ2(I-1)
      END IF
      END DO
      
C    
      QW=0
      QQW=1/REAL(QWS)
      DO I=1,TEMP1
      QWN=XX1(I)
      IF (ZZQ2(I).GE.QQW) EXIT
      END DO

C     
      TIME2(1:5)=1.0
      OPEN (55,FILE='NIHE_single.txt')
      OPEN (56,FILE='NIHE_cumulate.txt')
      OPEN (57,FILE='CUMULATIVE1.txt')
      OPEN (54,FILE='CUMULATIVE2.txt')
      DO 557 I=1,TEMP
      AA1=LOG10(Z(I)+1)
      BB1=LOG10(Y(I)+1) 
      WRITE (55,'(F5.1,1X,F10.4,1X)') X(I),AA1
      WRITE (56,'(F5.1,1X,F10.4,1X)') X(I),BB1
      WRITE (57,'(F5.1,1X,F15.8)') X(I),Y(I)
      WRITE (54,'(F5.1,1X,F15.8)') X(I),Y2(I)
557   CONTINUE
      CLOSE (54)
      CLOSE (55)
      CLOSE (56)
      CLOSE (57)
      
      
      INDEXX=MAXLOC(Z(1:TEMP))
      INDEXX2=INT(INDEXX(1))
      N=TEMP
      
      JK=1
      IF (CYC.EQ.1) THEN
      AP=II-OOP(JK)
      ELSE
      AP=1
      END IF
C----
      IF (CYC2.EQ.2) GOTO 559
      DO 501 JK=1,6
      MH=JK
      AP=II-OOP(JK)
      APN(JK)=II-OOP(JK)
      
559   M=2
      M1=M+1
      IF (MH.EQ.3) GOTO 411
      IF (MH.EQ.2) GOTO 410
      IF ((MH.EQ.4).OR.(MH.EQ.5)) GOTO 413
      IF (MH.EQ.1) GOTO 414
      IF (MH.EQ.6) GOTO 415
      IF (MH.EQ.7) GOTO 417
      IF (MH.EQ.9) GOTO 416
      IF (MH.EQ.8) GOTO 418
      
C--------(GFT)
411   CALL GFT(X,Y,Z2,A,N,M,R,MC,GG,LIMI,INDEXX2)
      IF ((GG.EQ.0).AND.(CYC2.NE.1)) GOTO 350
      IF ((GG.EQ.0).AND.(CYC2.EQ.1)) GOTO 502
      GOTO 412
      
C-------(MAXC)      
410   CALL MAXC(X,Z,N,MC,INDEXX2,LIMI,CG)
      IF (CG.EQ.1) GOTO 412
      IF (CYC2.EQ.1) GOTO 502
      GOTO 350
      
C--------(NDT)        
418   CALL NDT(X,Z,N,MC,INDEXX2,LIMI,CC,CS2)
      IF (CC.EQ.1) GOTO 412
      IF (CYC2.EQ.1) GOTO 502
      GOTO 350
      
C--------(MBASS)
415   CALL MBASS(X,Z,N,MC,INDEXX2,CG)
      IF (CG.EQ.1) GOTO 412
      IF (CYC2.EQ.1) GOTO 502
      GOTO 350     

C---------USELESS---(MBASS2)
416   CALL KSM(X,Y,Z,Z2,N,CC,B_VALE,B_DOUBLE,MC,HMM,INTVAL,LIMI,P)
      IF (CC.EQ.1) GOTO 918
      IF (CYC2.EQ.1) GOTO 502
      GOTO 350       

C----------(MBS)    
413   CALL MBS(X,Y,Z2,N,CC,MH,MC,LIMI,BIN,BY,JG,INDEXX2)
      IF (((YZZ.EQ.0).OR.(CC.EQ.0)).AND.(CYC2.NE.1)) GOTO 350
      IF (((YZZ.EQ.0).OR.(CC.EQ.0)).AND.(CYC2.EQ.1)) GOTO 502
      GOTO 412
  
C----------(EMR)    
414   CALL EMR(X,Y,Z,Z2,A,N,M,CC,MC,QWS,P,INDEXX2) 
      IF ((CC.EQ.0).AND.(CYC2.NE.1)) GOTO 350
      IF ((CC.EQ.0).AND.(CYC2.EQ.1)) GOTO 502
      GOTO 412     

C----------(HAM) 
417   CALL HAM(X,Y,Z,N,MC,INDEXX2,P2,CC)
      IF ((CC.EQ.0).AND.(CYC2.NE.1)) GOTO 350
      IF ((CC.EQ.0).AND.(CYC2.EQ.1)) GOTO 502
      
C    
412   IF (B_MO.EQ.1) THEN
      CALL MLE(X,Y,Z2,N,INDEXX2,MC,B_VALE,B_DOUBLE)
      ELSE IF (B_MO.EQ.2) THEN
      CALL LSM(X,Y,N,INDEXX2,B_VALE,B_DOUBLE,CC)  
      IF (CC.EQ.0) GOTO 350
      ELSE IF (B_MO.EQ.3) THEN
      CALL BP(TIME3(:SZ,:5),CATLOG(:SZ,3),SZ,B_VALE,B_DOUBLE)
      MC=0
      END IF
      
C-------
918   BZHU(AP,1,JK)=TIME2(1)
      BZHU(AP,2,JK)=TIME2(2)
      BZHU(AP,3,JK)=TIME2(3)
      BZHU(AP,4,JK)=TIME2(4)
      BZHU(AP,5,JK)=TIME2(5)
      BV(AP,1,JK)=B_VALE
      BV(AP,2,JK)=B_DOUBLE
      MMC(AP,JK)=MC
      IF (CYC.EQ.1) NUME(AP,JK)=QWS
      
C-----
      DO 518 PO=1,TEMP
      IF (ABS(MC-X(PO)).GT.0.0001) GOTO 518
      AZD(AP,JK)=LOG10(Y(PO))+B_VALE*MC
      EXIT
518   CONTINUE
      
C-------     
      IF (((B_VALE.LE.0.1).OR.(B_VALE.GT.5.0)).AND.(CYC2.NE.1)) GOTO 350
      IF (((B_VALE.LE.0.1).OR.(B_VALE.GT.5.0)).AND.(CYC2.EQ.1)) GOTO 502
      CALL fdate (stamp)
      IF (CYC.NE.1) GOTO 400
      IF (CYC2.EQ.1) GOTO 401
      WRITE (*,'(A,2X,I5,1X,F5.2,1X,F5.2,1X,F5.2)') 
     1    stamp(11:20),AP,B_VALE,B_DOUBLE,MC
      GOTO 555

401   WRITE (*,'(A,2X,I5,I5,1X,F5.2,1X,F5.2,1X,F5.2)') 
     1    stamp(11:20),AP,JK,B_VALE,B_DOUBLE,MC
      GOTO 501
502   OOP(JK)=OOP(JK)+1      
501   CONTINUE
      GOTO 555
      
400   WRITE (*,'(A,2X,I5,1X,I4,2(1X,2I2.2),
     1            1X,F5.2,1X,F5.2,1X,F5.2)') 
     2 stamp(11:20),AP,(TIME2(I22),I22=1,5),B_VALE,B_DOUBLE,MC
      GOTO 555
      
C--------
350   OOP(JK)=OOP(JK)+1
555   CATLOG=0.0
      FREQUENCY=0.0
      CUMULATIVE=0.0
      X=0.0
      Y=0.0
      Y2=0.0
      Z=0.0
      Z2=0.0
100   CONTINUE
552   CONTINUE
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC--------
      IF (CYC.EQ.1) GOTO 503
      KK2=0
      BEE=0.0
      MCC=0.0
      OPEN(102,FILE='TIME_B.txt')
      OPEN(112,FILE='TIME_Mc.txt')
      DO 300 KK=1,AP
      KK2=KK2+1
      BBTIME=(BZHU(KK,1,1)+((BZHU(KK,2,1)-1)*43200+BZHU(KK,3,1)*1440+
     1        BZHU(KK,4,1)*60+BZHU(KK,5,1))*1.0/525600.0)*1.0
      WRITE(102,'(F8.3,2X,2(2X,F8.5))') BBTIME,BV(KK,1,1),BV(KK,2,1)
      WRITE(112,'(F8.3,1X,F5.2))') BBTIME,MMC(KK,1)
      BEE=BEE+BV(KK,2,1)
300   CONTINUE
      CLOSE(112)
      CLOSE(102)
      
      JJK=0
      DO 302 KK=1,AP-1
      KK1=KK+JJK
      AB1=BV(KK,1,1)-1.0
      AB2=BV(KK+1,1,1)-1.0
      BBT(1)=(BZHU(KK,1,1)+((BZHU(KK,2,1)-1)*43200+BZHU(KK,3,1)*1440+
     1        BZHU(KK,4,1)*60+BZHU(KK,5,1))*1.0/525600.0)*1.0
      BBT(2)=(BZHU(KK+1,1,1)+((BZHU(KK+1,2,1)-1)*43200+BZHU(KK+1,3,1)*
     1        1440+BZHU(KK+1,4,1)*60+BZHU(KK+1,5,1))*1.0/525600.0)*1.0
      ABT=ABS(BBT(1)-BBT(2))
      BBT(1)=BBT(1)-2000
      BBT(2)=BBT(2)-2000
      X(1)=BBT(1)
      X(2)=BBT(2)
      Y(1)=AB1
      Y(2)=AB2
      N=2
      M=2
      M1=M+1
      S=0.0
      CALL NH_QBXF(X,Y,N,A,M,M1)
      IF ((AB1*AB2).LT.0) GOTO 303
      GOTO 305
      
303   DO 307 J=1,2
      IF (J.NE.1) GOTO 308
      AAW=A(1)
      ABT=ABS(BBT(J)-AAW)
      S=ABT*AB1/2.0
      GOTO 309
308   AAW=-A(1)/A(2)
      ABT=ABS(BBT(J)-AAW)
      S=ABT*AB2/2.0
309   SS1(J)=S
307   CONTINUE
      
      SS(KK1)=SS1(1)
      SS(KK1+1)=SS1(2)
      JJK=JJK+1
      GOTO 302
      
305   S=(AB1+AB2)*ABT*1.0/2.0
      SS(KK1)=S
302   CONTINUE
      
      DO 301 KK=1,KK2
      MCC=MCC+MMC(KK,1)
301   CONTINUE
      
      WRITE (*,*) 'The b-value calculation is finished'
      WRITE(*,'("Total data:",I8,2X,"Effective number:",
     1                I8,2X,"Data points:",I6)') COUNT,WQ,KK2
      
      IF (CM.EQ.3) GOTO 506
      WRITE(*,'("Average Mc:",1X,F6.2)') MCC/KK2*1.0
      WRITE(*,'("Average b-value uncertainty:",1X,F6.2)') BEE/KK2*1.0
      GOTO 507
      
506   WRITE(*,'("Mc:",1X,F6.2)') MC
      WRITE(*,'("b-value:",1X,F6.4)') B_VALE
      WRITE(*,'("b-value uncertainty:",1X,F6.4)') B_DOUBLE
      
507   CONTINUE
      OPEN (195,FILE='Average.txt')
      WRITE (195,'(F6.2,2X,F6.2)') MCC/KK2*1.0,BEE/KK2*1.0
      CLOSE (195)
      
      OPEN (196,FILE='JF.txt')
      DO 331 I=1,KK1
      WRITE (196,*) SS(I)*100
331   CONTINUE
      CLOSE (196)
      
      WRITE(*,'("Expectation:",1X,F6.2)') QWN
      
      IF (CYC.EQ.2) GOTO 615
503   OPEN (508,FILE='CYC_B.TXT')
      OPEN (509,FILE='CYC_MC.TXT')
      
      IF (CYC2.EQ.1) GOTO 504
      DO I=1,AP
      IF (NUME(I,1).EQ.0) CYCLE
      WRITE(508,'(I10,1X,F8.5)') NUME(I,1),BV(I,1,1)
      WRITE(509,'(I10,1X,F8.5)') NUME(I,1),MMC(I,1)
      END DO
      GOTO 505
      
504   DO JK=1,6
      IF (JK.LT.4) THEN
      WRITE(508,'(A4,I2)') '> -Z',JK-4
      WRITE(509,'(A4,I2)') '> -Z',JK-4
      ELSE
      WRITE(508,'(A4,I1)') '> -Z',JK-4
      WRITE(509,'(A4,I1)') '> -Z',JK-4
      END IF
      DO I=1,APN(JK)
      MJ=MMC(I,JK)
      IF (NUME(I,JK).EQ.0) EXIT
      WRITE(508,'(I10,1X,F8.5)') NUME(I,JK),BV(I,1,JK)
      WRITE(509,'(I10,1X,F8.5)') NUME(I,JK),MJ
      END DO
      END DO
      
505   CLOSE (508)
      CLOSE (509)

615   WRITE(*,*) 'CLOSE ALL'
      PAUSE
      STOP
      END
      
      
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCC---------Sorting time array---------------CCC      
      SUBROUTINE TIME_sorting(X,Y,N,N2,N3)
      INTEGER I,J,N,N2,P(5),X(N,N2),N3
      DOUBLE PRECISION Y(N,N3),P2(3)

      DO 10 I=1,N-1
      DO 20 J=I+1,N
          
      IF (X(I,1).GT.X(J,1)) THEN
      P(1:5)=X(J,1:5)
      X(J,1:5)=X(I,1:5)
      X(I,1:5)=P(1:5)
      P2(1:3)=Y(J,1:3)
      Y(J,1:3)=Y(I,1:3)
      Y(I,1:3)=P2(1:3)
      ELSE IF (X(I,1).EQ.X(J,1)) THEN
      IF (X(I,2).GT.X(J,2)) THEN
      P(1:5)=X(J,1:5)
      X(J,1:5)=X(I,1:5)
      X(I,1:5)=P(1:5)
      P2(1:3)=Y(J,1:3)
      Y(J,1:3)=Y(I,1:3)
      Y(I,1:3)=P2(1:3)
      ELSE IF (X(I,2).EQ.X(J,2)) THEN  
      IF (X(I,3).GT.X(J,3)) THEN
      P(1:5)=X(J,1:5)
      X(J,1:5)=X(I,1:5)
      X(I,1:5)=P(1:5)
      P2(1:3)=Y(J,1:3)
      Y(J,1:3)=Y(I,1:3)
      Y(I,1:3)=P2(1:3)
      ELSE IF (X(I,3).EQ.X(J,3)) THEN    
      IF (X(I,4).GT.X(J,4)) THEN
      P(1:5)=X(J,1:5)
      X(J,1:5)=X(I,1:5)
      X(I,1:5)=P(1:5)
      P2(1:3)=Y(J,1:3)
      Y(J,1:3)=Y(I,1:3)
      Y(I,1:3)=P2(1:3)
      ELSE IF (X(I,4).EQ.X(J,4)) THEN   
      IF (X(I,5).GE.X(J,5)) THEN    
      P(1:5)=X(J,1:5)
      X(J,1:5)=X(I,1:5)
      X(I,1:5)=P(1:5)
      P2(1:3)=Y(J,1:3)
      Y(J,1:3)=Y(I,1:3)
      Y(I,1:3)=P2(1:3)
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF

20    CONTINUE
10    CONTINUE   
      

      
      RETURN
      END SUBROUTINE TIME_sorting
          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCC---------Find maximum/minimum value of time array---------------CCC      
C     P: 1-MAX  2-MIN
      SUBROUTINE LOOKF(X,N,N2,P,J)
      INTEGER I,J,N,N2,P,X(N,N2)
      
      IF (P.EQ.1) THEN
      J=1
      DO 10 I=1,N
      IF (X(J,1).GT.X(I,1)) THEN
      CONTINUE
      ELSE IF (X(J,1).EQ.X(I,1)) THEN
      IF (X(J,2).GT.X(I,2)) THEN
      CONTINUE
      ELSE IF (X(J,2).EQ.X(I,2)) THEN  
      IF (X(J,3).GT.X(I,3)) THEN
      CONTINUE
      ELSE IF (X(J,3).EQ.X(I,3)) THEN    
      IF (X(J,4).GT.X(I,4)) THEN
      CONTINUE
      ELSE IF (X(J,4).EQ.X(I,4)) THEN   
      IF (X(J,5).GE.X(I,5)) THEN    
      CONTINUE
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      ELSE
      GOTO 20
      END IF
      GOTO 10
20    J=I
10    CONTINUE   
      
      ELSE IF (P.EQ.2) THEN
      J=1
      DO 30 I=1,N
      IF (X(J,1).LT.X(I,1)) THEN
      CONTINUE
      ELSE IF (X(J,1).EQ.X(I,1)) THEN
      IF (X(J,2).LT.X(I,2)) THEN
      CONTINUE
      ELSE IF (X(J,2).EQ.X(I,2)) THEN  
      IF (X(J,3).LT.X(I,3)) THEN
      CONTINUE
      ELSE IF (X(J,3).EQ.X(I,3)) THEN    
      IF (X(J,4).LT.X(I,4)) THEN
      CONTINUE
      ELSE IF (X(J,4).EQ.X(I,4)) THEN   
      IF (X(J,5).LE.X(I,5)) THEN    
      CONTINUE
      ELSE
      GOTO 40
      END IF
      ELSE
      GOTO 40
      END IF
      ELSE
      GOTO 40
      END IF
      ELSE
      GOTO 40
      END IF
      ELSE
      GOTO 40
      END IF
      GOTO 30
40    J=I
30    CONTINUE 
      
      ELSE
      WRITE (6,*) 'ERROR: P'
      END IF
      
      RETURN
      END SUBROUTINE LOOKF
          
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCCC---------------KS-test---------------CCCCCCCCCC
      SUBROUTINE CEF(X,Y,N,TOP,IT,CO)
      INTEGER N,I,IT,J
      DOUBLE PRECISION X(N),Y(IT),STEP,CO(IT),TOP
      
      DO I=1,IT
      NUMER=0
      DO J=1,N
      STEP=TOP*1.0/IT*I*1.0
      IF (STEP+0.0001.GE.X(J)) NUMER=NUMER+1
      END DO
      Y(I)=NUMER*1.0/N*1.0
      CO(I)=STEP
      END DO
      
      RETURN
      END SUBROUTINE CEF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCCC---------------Create rank---------------CCCCCCCCCC
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
10    I=I+1
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
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC---Add the same parameters in the found group---CCCC  
C     JEMP(N,M)
C     NUM=key
      SUBROUTINE SORTSUM(JEMP,N,M,NUM1,NUM2,VALUE,SUM_NU)
      INTEGER N,M,NUM1,NUM2
      DOUBLE PRECISION JEMP(N,M),TEMP,VALUE,SUM_NU
      LOGICAL MASK(N)
      MASK = (JEMP(:,NUM1)==VALUE)
      SUM_NU = SUM(PACK(JEMP(:,NUM2),MASK))
      
      RETURN
      END SUBROUTINE SORTSUM         
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC----Find the maximum value for a given value in an array----CCCCC   
C     ARR1-2(ROW)
C     NUM1=key
C     NUM2=value
      FUNCTION GETMAX(ARR1,ARR2,VALUE,ROW) 
      INTEGER ROW
      DOUBLE PRECISION ARR1(ROW),ARR2(ROW),VALUE
      DOUBLE PRECISION MAXNUM,GETMAX
      LOGICAL MASK(ROW)
      MASK = (ARR1(:)==VALUE)
      MAXNUM = MAXVAL(PACK(ARR2(:),MASK))
      
      GETMAX=MAXNUM
      END FUNCTION GETMAX
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC----INTERSECTION-Determine whether two line segments intersect subroutine----CCCCC           
      SUBROUTINE INTERSECTION(X1,Y1,X2,Y2,X3,Y3,CHUS)
      INTEGER CHUS
      DOUBLE PRECISION X1,Y1,X2,Y2,X3,Y3,X4,Y4,T1,T2,
     1                 XL1(2),XL2(2),XL3(2),T3,T4,XL4(2),
     2                 XL5(2),XL6(2)
      
      X4 = 180
      Y4 = Y3
      CHUS = 0

      XL1(1)=X4-X3
      XL1(2)=Y4-Y3
      XL2(1)=X1-X3
      XL2(2)=Y1-Y3
      XL3(1)=X2-X3
      XL3(2)=Y2-Y3
      
      XL4(1)=X2-X1
      XL4(2)=Y2-Y1
      XL5(1)=X3-X1
      XL5(2)=Y3-Y1
      XL6(1)=X4-X1
      XL6(2)=Y4-Y1
      
      T1 = XL1(1)*XL2(2)-XL1(2)*XL2(1)
      T2 = XL1(1)*XL3(2)-XL1(2)*XL3(1)
      
      T3 = XL4(1)*XL5(2)-XL4(2)*XL5(1)
      T4 = XL4(1)*XL6(2)-XL4(2)*XL6(1)
      
      IF (T1*T2<0.AND.T3*T4<0) then
          CHUS=1
      ELSE
          CHUS=0
      END IF
      
      RETURN
      end subroutine INTERSECTION         


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC---Sort subroutines (largest to smallest)---CCCCCC  
C     JEMP(N,M)
C     NUM=key
      SUBROUTINE SORTNUM(JEMP,N,M,NUM)
      INTEGER I,J,N,M,NUM,j1
      DOUBLE PRECISION JEMP(N,M),TEMP
      
      do i = 1, N
      do j = i + 1, N
      if (JEMP(i, 1) < JEMP(j, 1)) then
      temp = JEMP(i, 1)
      JEMP(i, 1) = JEMP(j, 1)
      JEMP(j, 1) = temp
      
      do j1 = 2, M
      temp = JEMP(i, j1)
      JEMP(i, j1) = JEMP(j, j1)
      JEMP(j, j1) = temp
      end do
      end if
      end do
      end do
      
      RETURN
      END SUBROUTINE SORTNUM
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC----CDFX-Calculate the cumulative normal distribution function----CCCCC   
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC-------KSM--------CCCCCCCCCCCCCCCC
      SUBROUTINE KSM(X,Y,Z,W,N,CC,BV,BD,MC,HMM,INTVAL,LIMI,P)
      DIMENSION X(N),Y(N),X1(N),Y1(N),Z(N),Z1(N),W(N),W1(N),
     2          BB(100,2),WZ(1),MCB(100),X2(N),Y2(N),Z2(N),W2(N)
      DOUBLE PRECISION X,Y,X2,MC,Y2,Z,Z2,Y1,X1,Z1,W1,
     1    MA,BV,BD,W,W2,HMM,BB,MCB,CA(50,2),P,
     4    MAX_IN,GL(100),INTVAL,KS_P,SLOP(100),SLOP_MID
      INTEGER N,I,I2,CC,WZ,I1,TOP,J,LIMI

      I1=0
      CC=1
 
      DO 10 I=1,N
      IF (I.GE.2) THEN
      I1=I1+1
      SLOP(I1)=(LOG(Y(I)+1)-LOG(Y(I-1)+1))/0.1
      IF (I.GE.5) THEN
      SLOP_MID=SUM(SLOP(1:I1-1))/(I1-1)*0.15
      IF ((I.GE.3).AND.(SLOP(I1).LT.SLOP_MID)) GOTO 20
C      IF (SLOP(I1).LT.0.4) GOTO 20
      END IF
      END IF
10    CONTINUE 
      CC=0
      GOTO 100
      
20    MC=X(I)
      MA=SUM(W(1:I))/Y(I)
      BV=1.0/((MA-MC)*2.30258509)
      BD=BV/SQRT(Y(I))*1.0
      
100   RETURN 
      END    
      
CCCCCCCCCC--median-based analysis of the segment slope--CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCC----MBASS-----CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MBASS(X,Y2,N,MC,INDEX,CG)
      INTEGER N,I,CG,J,J2,SA_SITE(1),PM1(1),N1,N2,INDEX,I1,U(400),U1(20)
      DOUBLE PRECISION X(N),Y2(N),MC,SLOPE(N,2),Y22(N),P(37,2),
     1                 SA(N),SR(N),Wcrit,SW,XIG,Z(N),
     2                 PM(20,2),PM2(1),Y(N)
      
      P(:,1)=(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     1        21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 /)
      P(:,2)=(/ 12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,
     1        2.228,2.201,2.179,2.160,2.145,2.131,2.120,2.110,2.101,
     2        2.093,2.086,2.08,2.074,2.069,2.063,2.06,2.056,2.052,2.048,
     3        2.045,2.042,2.021,2.009,2.003,1.994,1.990,1.987,1.984 /)
      
      U(:)=(/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     1          0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,
     2          0,0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,
     3          0,0,0,0,1,2,3,4,4,5,6,7,8,9,10,11,11,12,13,13,
     4          0,0,0,1,2,3,5,6,7,8,9,11,12,13,14,15,17,18,19,20,
     5          0,0,1,2,3,5,6,8,10,11,13,14,16,17,19,21,22,24,25,27,
     6          0,0,1,3,5,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,
     7          0,0,2,4,6,8,10,13,15,17,19,22,24,26,29,31,34,36,38,41,
     8          0,0,2,4,7,10,12,15,17,20,23,26,28,31,34,37,39,42,45,48,
     9          0,0,3,5,8,11,14,17,20,23,26,29,33,36,39,42,45,48,52,55,
     9          0,0,3,6,9,13,16,19,23,26,30,33,37,40,44,47,51,55,58,62,
     9          0,1,4,7,11,14,18,22,26,29,33,37,41,45,49,53,57,61,65,69,
     9          0,1,4,8,12,16,20,24,28,33,37,41,45,50,54,59,63,67,72,76,
     9          0,1,5,9,13,17,22,26,31,36,40,45,50,55,59,64,67,74,78,83,
     9         0,1,5,10,14,19,24,29,34,39,44,49,54,59,64,70,75,80,85,90,
     9         0,1,6,11,15,21,26,31,37,42,47,53,59,61,70,75,81,86,92,98,
     9        0,2,6,11,17,22,28,34,39,45,51,57,63,67,75,81,87,93,99,105,
     9       0,2,7,12,18,24,30,36,42,48,55,61,67,74,80,86,93,99,106,112,
     9      0,2,7,13,19,25,32,38,45,52,58,65,72,78,85,92,99,106,113,119,
     9     0,2,8,13,20,27,34,41,48,55,62,69,76,83,90,98,105,112,119,127
     9 /)
      
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
10    CONTINUE
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
      
40    CONTINUE
      CG=0
      GOTO 30
20    MC=X(N-I+1)
      INDEX=N-I+1
30    RETURN
      END      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC----GFT-----CCCCCCCCCCCCCC
      SUBROUTINE GFT(X,Y,W,AQ,NQ,MQ,R,MC,GG,NN,TIM)
      DIMENSION X(NQ),Y(NQ),AQ(MQ+1),X2(NQ),Y2(NQ),YY22(NQ),
     1          MII(100,2),WZ(1),RR(100),W(NQ),W2(NQ)
      INTEGER NQ,MQ,I,J,RF,GG,M1,NN,KI,KII,WZ,IK,RR,AYQ,R,TIM
      DOUBLE PRECISION X,Y,AQ,MC,X2,SUME,YY,Y2,YY22,MII,YY2,
     1                 W,W2,AV,MA,BV,YZ,YYZ(100),BD
      
      J=0
      KI=0
      IK=0
      KII=0
      SUME=0.0
      YY=0.0
      MC=X(NQ)
      GG=1
      M1Q=MQ+1
      
      GOTO 41
      
39    IF (MC.GE.X(1)) GOTO 12
      MC=X(NQ-IK)
      J=0
      YY=0.0
      SUME=0.0
      GOTO 41
12    IF (KI.EQ.0) GOTO 13
      WZ=MAXLOC(MII(1:KII,2))
      MC=MII(WZ(1),1)

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
      MII(KII,2)=RF
      
      GOTO 39
      
36    DO I=1,NQ
      IF (ABS(X(I)-MC).LT.0.001) THEN
      TIM=I
      GOTO 33
      END IF
      END DO
      
33    RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC----MBS-----CCCCCCCCCCCCCC
      SUBROUTINE MBS(X,Y,Z,N,CC,MH,MC,NN,BIN,BY,JG,TIM)
      DIMENSION X(N),Y(N),X2(N),Y2(N),Z(N),Z2(N),
     1          X3(N),Y3(N),Z3(N),WZ(1)
      DOUBLE PRECISION X,Y,MB,MC,Z,X2,Y2,Z2,MA,JG,
     1    MC2,BV,BD,X3,Y3,Z3,BV2,BD2,MA2,DB,DB2,MM,BAVE
     2    ,MC3,BY,CHUC(100,2)
      INTEGER(kind=4) N,I,K,K2,CC,NN,MH,KI,IIK,IIK2,WZ,MM2,
     1                BIN,I1,J,AZ(1),IJJ,TIM
      
      MC=X(N)
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
      IJJ=0
      GOTO 70
      
60    IF (MC.GE.X(1)) GOTO 80
      IJJ=IJJ+1
      GOTO 90
80    IF ((MH.EQ.4).AND.(J.GT.1)) THEN
      AZ=MINLOC(CHUC(1:J-1,2))
      MC=CHUC(AZ(1),1)
      GOTO 50
      END IF
      CC=0
      GOTO 50
90    MC=X(N-IJJ)
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
      DB2=2.3*(BV**2)*SQRT(MM/(Y2(K)*(Y2(K)-1)))*
     1    SQRT(Y2(K)*(Y2(K)-1))**0.6
      
      IF (DB2.GE.DB) GOTO 50
      GOTO 60
      
130   IF (ABS(BV2-BV).LE.BY) GOTO 50 
      IF (ABS(BV2-BV).LE.BY+0.001) THEN
      CHUC(J,1)=MC
      CHUC(J,2)=ABS(BV2-BV)
      J=J+1
      END IF
      GOTO 60
   
50    DO I=1,N
      IF (ABS(X(I)-MC).LT.0.001) THEN
      TIM=I
      GOTO 33
      END IF
      END DO
      
33    CONTINUE
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC-------EMR-CCCCCCCCC
      SUBROUTINE EMR(X,Y,Z,W,A,N,M,CC,MC,NUM,P,TIM)
      DIMENSION X(N),Y(N),A(M+2),X2(N),Y2(N),Z(N),Z2(N),
     1          DJJ(100),W(N),W2(N),MJJ(100),WPP(1),
     2          BB(100),WZ(1),MCB(100),CO(100,2),
     3          EMRA(100,3),YZ(100),Z22(100),X22(100),QW(1)
      DOUBLE PRECISION X,Y,A,X2,MC,Y2,Z,Z2,DJ,DJJ,MC2,
     1    MA2,BV2,BD2,W,W2,MJJ,ABZ,BB,MCB,CO,
     2    AV,MIU,XGM1,XGM,YZ,Z22,X22,MGASS,ZDZ,
     3    SZ2(100,4),D_REAL(100),D_CRIT,S,A_X,B_X,ZX2,CA(50,2),
     4    MAX_IN,GL(100),P,KS_P,TOP,EN(10000,2)
      INTEGER N,M,I,I2,CC,M1,WPP,WZ,I1,QW,TIM,SUM_EN,SUM_EN2,
     1        NUM,IJ,SZ1(100),XS
      EXTERNAL MGASS
      
      MC=X(N)
      I2=0
      I1=0
      CC=1
      M1=M+1
      DJJ=0.0
      J=1
      IJ=0
      GOTO 300
      
200   IF (MC.EQ.X(1)) GOTO 280
      MC=X(N-J)
      J=J+1
      I2=0
      I1=0
      A=0
      CA=0
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

C     （G-R law）
500   MC2=MC
      MA2=SUM(W2(1:I2))*1.0/Y2(I2)*1.0
      BV2=(MA2-MC2)*2.30258509
      BV2=1.0/BV2*1.0
      BD2=BV2*1.0/SQRT(Y(I2))*1.0
      MCB(J)=MC
      BB(J)=I2
      
      DO I=1,I2
      Z22(I)=EXP(-BV2*X2(I)*2.3)
      IF (Z22(I).LT.0) Z22(I)=0.
      END DO
      Z22(1:I2)=Z22(1:I2)/MAXVAL(Z22(1:I2))*Z(I2)
      
C     (CDF)
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
      
C-----K-S
      TOP=MAXVAL(X(1:N))
      XS=INT((TOP+0.1)*10)
      
      SUM_EN=0
      SUM_EN2=0
      DO I=1,N
      EN(1+SUM_EN:INT(Z(I)+SUM_EN),1)=X(I)
      EN(1+SUM_EN2:INT(Z22(I)+SUM_EN2),2)=X(I)
      SUM_EN=SUM_EN+INT(Z(I))
      SUM_EN2=SUM_EN2+INT(Z22(I))
      END DO
      
      CALL CEF(EN(1:SUM_EN,1),CA(:XS,1),SUM_EN,TOP,XS,CO(:XS,1))
      IF (SUM_EN.EQ.0) CA(:XS,1)=0
      CALL CEF(EN(1:SUM_EN2,2),CA(:XS,2),SUM_EN2,TOP,XS,CO(:XS,2))
      IF (SUM_EN2.EQ.0) CA(:XS,2)=0
      MAX_IN=MAXVAL(ABS(CA(:,1)-CA(:,2)))
      
      GL(J)=MAX_IN
      KS_P=SQRT(-LOG(P/2.0)*0.5)*SQRT(2.0/N)
      IF (MAX_IN.LE.KS_P) GOTO 100 
      GOTO 200
      
280   IF (J.NE.1) THEN
      WZ=MINLOC(GL(1:J-1))
      TIM=BB(WZ(1))
      MC=MCB(WZ(1))
      ELSE
      CC=0    
      END IF
100   RETURN 
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC-------HAM--------CCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HAM(X,Y,Z,N,MC,TIM,P,CC)
      DIMENSION X(N),Y(N),Z(N),X2(N),Y2(N),Z2(N),TIMM(N),V(N-1,2),
     1          NV(N,2),NH(N,2),NH2(N,2),A(3),DT(N),INDX(1)
      DOUBLE PRECISION X,Y,Z,X2,Y2,Z2,MC,P,LS,LS1,V,NV,NH,A,NH2,
     1                 DT
      INTEGER N,TIM,TIMM,I,I2,NHP,M,M1,CC,INDX
      
      CC=0
      M=2
      M1=3
      MC=X(N)
      DO J=1,N-1
          
      I2=0
      DO I=1,N
      IF ((X(I)-MC).GT.0) THEN
      I2=I2+1
      X2(I2)=X(I)
      Y2(I2)=Y(I)
      Z2(I2)=Z(I)
      TIMM(J)=I2
      IF (Z(I).LE.0) THEN
      Z2(I2)=1
      IF (I2.EQ.1) THEN
      Y2(I2)=1
      ELSE
      Y2(I2)=Y2(I2-1)+Z2(I2)   
      END IF
      END IF
      END IF

      END DO 
      
      LS=0
      LS1=0
      DO I=1,TIMM(J)
      LS=LS+1.0/Y2(I)
      LS1=LS1+1.0/Z2(I)
      END DO
      V(J,1)=I2/LS
      V(J,2)=I2/LS1
      MC=X(N-J)
      END DO
      
      DO I=1,N-1
      NV(I,1)=X(I)
      NV(I,2)=V(I,2)
      END DO
      
      NHP=INT((N-1)/3)+1
      DO I=1,NHP
      NH(I,1)=NV(I+NHP,1)
      NH(I,2)=NV(I+NHP,2)
      END DO
      
      CALL NH_QBXF(NH(:NHP,1),NH(:NHP,2),NHP,A,M,M1)
      DO I=1,N-1
      NH2(I,1)=NV(N-I,1)
      NH2(I,2)=NV(N-I,1)*A(2)+A(1)
      DT(I)=ABS(NH2(I,2)-NV(N-I,2))
      IF (DT(I).LT.P) THEN
      MC=NH2(I,1)
      TIM=N-I
      CC=1
      GOTO 10
      END IF
      END DO
      
      INDX=MINLOC(DT(:N-1))
      MC=NH2(INDX(1),1)
      TIM=N-INDX(1)
      CC=1
10    RETURN 
      END      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC-------MAXC--------CCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MAXC(X,Z,N,MC,TIM,LIM,CC)
      DIMENSION X(N),Z(N),A(1)
      DOUBLE PRECISION X,Z,MC
      INTEGER N,TIM,A,LIM,CC
      
      CC=1
      IF ((SUM(Z).LT.LIM).OR.(N.LT.3)) THEN
      CC=0
      GOTO 10
      END IF
      A=MAXLOC(Z)
      TIM=A(1)
      MC=X(TIM)
      
10    RETURN 
      END       
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCC-------NDT-----CCCCCCCCCCCCCCCCCCC
      SUBROUTINE NDT(X,Z,N,MC,TIM,LIM,CC,CS)
      DIMENSION X(N),Z(N),WZ(1),LS(1000,2),EN(10000),ZD(100),CO(100),
     1          MCB(1000)
      DOUBLE PRECISION X,Z,MC,RM,RP,LS,EN,ZD,CO,MCB,MAG_MIN,
     1                 MAG_MAX
      INTEGER N,TIM,LIM,CC,WZ,COUNT,ND,I1,J,I,SUM_EN,XS,CS,SUM_Z
      
      IF (SUM(Z(:N)).LT.LIM) GOTO 10
      
      CC=1
      EN=0.
      X2=0.
      Z2=0.
      
      MAG_MAX=MAXVAL(X(:N))
      MAG_MIN=MINVAL(X(:N))
      CALL RANDOM_SEED()
      
C     自举抽样法
      DO I=1,CS
      
      SUM_Z=0
      I1=0
      LS=0.
      EN=0.
      ZD=0.
      CO=0.
      DO J=1,1000
      CALL RANDOM_NUMBER(RM)
      CALL RANDOM_NUMBER(RP)
      RM = RM*(MAG_MAX-MAG_MIN)
      WZ = MINLOC(ABS(X(:)-RM))
      RP = RP*Z(WZ(1))
      IF (RP.LT.1) CYCLE
      SUM_Z=SUM_Z+INT(RP)
      IF (SUM_Z.GE.SUM(Z)*0.8) EXIT
      I1=I1+1
      LS(I1,1)=X(WZ(1))
      LS(I1,2)=INT(RP)
      END DO
      
      SUM_EN=0
      DO J=1,I1
      EN(1+SUM_EN:INT(LS(J,2)+SUM_EN))=LS(J,1)
      SUM_EN=SUM_EN+INT(LS(J,2))
      END DO
      
C     判断数组中的不重复点
C      CALL QSORT(LS(:I1,2),I1,2,1)
C      ND=0
C      DO I=1,I1
C      ND=ND+1
C      I1=0
C      DO J=I,I1
C      IF (ABS(LS(I,1)-LS(J,1)).LE.0.0001) THEN
C      X2(ND)=LS(I,1)
C      Z2(ND)=Z2(ND)+LS(I,2)
C      I1=I1+1
C      END IF
C      EXIT
C      END DO
C      I=I1+I
C      END DO
      
      XS=INT((MAXVAL(EN(:SUM_EN))+0.1)*10)
      CALL CEF(EN(:SUM_EN),ZD(:XS),SUM_EN,MAXVAL(EN(:SUM_EN)),
     1         XS,CO(:XS))
      
C     找到单次样本中符合显著度的震级
      DO J=1,XS
      IF (ZD(J).GE.0.95) THEN
      MCB(I)=CO(J)
      EXIT
      END IF
      END DO
      
      END DO
      
      XS=INT(MAXVAL(MCB(:CS)))*10
      CALL CEF(MCB(:CS),ZD(:XS),CS,MAXVAL(MCB(:CS)),XS,CO(:XS))
      
C     找到全部样本中符合显著度的震级
      DO J=1,XS
      IF (ZD(J).GE.0.95) THEN
      MC=CO(J)
      EXIT
      END IF
      END DO
      
      DO J=1,N
      IF (ABS(X(J)-MC).LE.0.0001) THEN
      TIM=J
      EXIT
      END IF
      END DO
      
10    RETURN 
      END     
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC-------b-positive------CCCCCCCCCCCCCCCCCCC
      SUBROUTINE BP(X,MA,N,BV,BD)
      DOUBLE PRECISION MA,BV,BD,TIME,MAG,TEMP,M1,M2,MC
      INTEGER X,N,I,J,I1,I2
      DIMENSION X(N,5),MA(N),TIME(N),MAG(N),TEMP(2),M1(N),M2(N)
      PARAMETER (LN10=2.302585093)
      
      TIME(:)=X(:,1)*1.0+X(:,2)/12.0+X(:,3)/365.0+
     1        X(:,4)/8760.0+X(:,5)/525600.0
      MAG(:)=MA(:)
      
      DO I=1,N-1
      DO J=I,N
      IF (TIME(I).GT.TIME(J)) THEN
      TEMP(1) = TIME(I)
      TEMP(2) = MAG(I)
      TIME(I) = TIME(J)
      TIME(J) = TEMP(1)
      MAG(I) = MAG(J)
      MAG(J) = TEMP(2)
      END IF
      END DO
      END DO
      
      I1=0
      I2=0
      OPEN (881,FILE='BP_EARTHQUAKE.TXT')
      DO I=1,N-1
      TEMP(1)=MAG(I+1)-MAG(I)
      IF (TEMP(1).GT.0) THEN
      I1=I1+1
      M1(I1)=TEMP(1)
      ELSE
      I2=I2+1
      M2(I2)=TEMP(1)
      END IF
      WRITE(881,*) TIME(I),TEMP(1)
      END DO
      CLOSE(881)
      
      MC=MINVAL(M1(1:I1))
      IF (MC.LT.0.2) MC=0.2
      BV=1.0/(LN10*(SUM(M1(1:I1))/I1-MC))
      BD=BV/SQRT(I1*1.0)
      
      RETURN 
      END      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC-------MLE-------CCCCCCCCCCCCCCCCCCC
      SUBROUTINE MLE(X,Y,Z2,N,INDEX,MC,BV,BD)
      DOUBLE PRECISION X,BV,BD,MC,Y,Z2,MC_AV
      INTEGER N,I,INDEX
      DIMENSION X(N),Y(N),Z2(N)
      PARAMETER (LN10=2.30258509)
      
      MC=X(INDEX)
      MC_AV=SUM(Z2(1:INDEX))*1.0/Y(INDEX)
      BV=1.0/((MC_AV-MC+0.05)*LN10)
      BD=BV/SQRT(Y(INDEX))
      
      RETURN 
      END   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC-------LSM------CCCCCCCCCCCCCCCCCCC
      SUBROUTINE LSM(X,Y,N,INDEX,BV,BD,CC)
      DOUBLE PRECISION X,Y,BV,BD,X2,Y2,A
      INTEGER N,I,INDEX,M,M1,CC
      DIMENSION X(N),Y(N),X2(INDEX),Y2(INDEX),A(3)
      PARAMETER (LN10=2.30258509)
      
      M=2
      M1=3
      CC=1
      
      DO I=1,INDEX
      X2(I)=X(I)
      Y2(I)=LOG10(Y(I))
      IF (Y(I).LE.1.0) Y2(I)=0
      END DO
      
      CALL NH_QBXF(X2,Y2,INDEX,A,M,M1)
      
      BV=-A(2)
      BD=BV/SQRT(Y(INDEX))
      
      IF (ABS(A(1)-100).LT.0.001) CC=0
      
      RETURN 
      END 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCC-----Chebyshev fitting------CCCCCCCCCCCCC
      SUBROUTINE NH_QBXF(X,Y,N,A,M,M1)
      DIMENSION X(N),Y(N),A(M1),IX(20),H(20)
      DOUBLE PRECISION X,Y,A,H,HA,HH,Y1,Y2,H1,H2,D,HM
      INTEGER I,J,L,M,N,M1,II,IX
      
      IF (N.LT.3) THEN
      A(1)=-100
      GOTO 140
      END IF
      
      DO 5 I=1,M1
5     A(I)=0.0
      IF (M.GE.N) M=N-1
      IF (M.GE.20) M=19
      M1=M+1
      HA=0.0
      IX(1)=1
      IX(M1)=N
      L=(N-1)/M
      J=L
      DO 10 I=2,M
          IX(I)=J+1
          J=J+L
10    CONTINUE
20    HH=1.0
      DO 30 I=1,M1
          A(I)=Y(IX(I))
          H(I)=-HH
          HH=-HH
30    CONTINUE
      DO 50 J=1,M
      II=M1
      Y2=A(II)
      H2=H(II)
      DO 40 I=J,M
          D=X(IX(II))-X(IX(M1-I))
          Y1=A(M-I+J)
          H1=H(M-I+J)
          A(II)=(Y2-Y1)/D
          H(II)=(H2-H1)/D
          II=M-I+J
          Y2=Y1
          H2=H1
40    CONTINUE
50    CONTINUE
      HH=-A(M1)/H(M1)
      DO 60 I=1,M1
60    A(I)=A(I)+H(I)*HH
      DO 80 J=1,M-1
          II=M-J
          D=X(IX(II))
          Y2=A(II)
          DO 70 K=M1-J,M
              Y1=A(K)
              A(II)=Y2-D*Y1
              Y2=Y1
              II=K
70        CONTINUE
80    CONTINUE
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
      DO 100 I=1,N
          IF(I.EQ.IX(J)) THEN
              IF (J.LT.M1)J=J+1
          ELSE
          H2=A(M)
          DO 90 K=M-1,1,-1
90        H2=H2*X(I)+A(K)
          H2=H2-Y(I)
          IF (ABS(H2).GT.HM) THEN
              HM=ABS(H2)
              H1=H2
              IM=I
          END IF
      END IF
100   CONTINUE
      IF(IM.EQ.IX(1)) RETURN
      I=1
110   IF(IM.GE.IX(I)) THEN
          I=I+1
          IF (I.LE.M1) GOTO 110
      END IF
      IF (I.GT.M1) I=M1
      IF (I.EQ.(I/2)*2) THEN
          H2=HH
      ELSE
          H2=-HH
      END IF 
      IF(H1*H2.GE.0.0)THEN
          IX(I)=IM
          GOTO 20
      END IF
      IF(IM.LT.IX(1)) THEN
          DO 120 J=M,1,-1
120       IX(J+1)=IX(J)
          IX(1)=IM
          GOTO 20
      END IF
      IF(IM.GT.IX(M1)) THEN
          DO 130 J=2,M1
130       IX(J-1)=IX(J)
          IX(M1)=IM
          GOTO 20
      END IF
      IX(I-1)=IM
      GOTO 20
      
140   RETURN
      END
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC---Gamma function---CCCCCCCCCCCCCCCCCCCCCC
      FUNCTION MGAM1(X)
      DOUBLE PRECISION MGAM1,X
      DOUBLE PRECISION Y,T,S,U,A(11)
      DATA A/0.0000677106,-0.0003442342,0.0015397681,
     *    -0.0024467480,0.0109736958,-0.0002109075,
     *    0.0742379071,0.0815782188,0.4118402518,
     *    0.4227843370,1.0/
      IF (X.LE.0.) THEN
      WRITE( * ,* )'ERR **x<0!'
      MGAM1=-1.0
      RETURN
      END IF
      Y=X
      IF (Y.LE.1.0) THEN
      T=1.0/(Y * (Y+1.0))
      Y=Y+2.0
      ELSE IF (Y.LE.2.0) THEN
      T=1.0/Y
      Y=Y+1.0
      ELSE IF (Y.LE.3.0)THEN
      T=1.0
      ELSE
      T=1.0
10    IF(Y.GT.3.0)THEN
      Y=Y-1.0
      T=T*Y
      GOTO 10
      END IF
      END IF
      S=A(1)
      U=Y-2.0
      DO 20 I=1,10
20    S=S*U+A(I+1)
      S=S*T
      MGAM1=S
      RETURN
      END
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC---incomplete gamma function---CCCCCCCCCCCCCCCCCC
      FUNCTION MGAM2(A,X)
      DOUBLE PRECISION MGAM2,A,X
      DOUBLE PRECISION MGAM1,P,Q,D,S,S1,P0,Q0,P1,Q1,QQ
      
      IF ((A.LE.0.0).OR.(X.LT.0.0)) THEN
      IF (A.LE.0.0)THEN
      WRITE(*,*) 'ERR * * A<=01'
      END IF
      IF (X.LT.0.0) THEN
      WRITE(*,*) 'ERR * * X<0!'
      END IF
      MGAM2=-1.0
      END IF
      IF(X+1.0.EQ.1.0) THEN
      MGAM2=0.0
      RETURN
      END IF
      IF(X.GT.1.0D+35) THEN
      MGAM2=1.0
      RETURN
      END IF
      Q=LOG(X)
      Q=A*Q
      QQ=EXP(Q)
      IF (X.LT.1.0+A)THEN
      P=A
      D=1.0/A
      S=D
      DO 10 N=1,100
      P=1.0+P
      D=D*X/P
      S=S+D
      IF (ABS(D).LT.ABS(S)*1.0D-07) THEN
      S=S*EXP(-X)*QQ/MGAM1(A)
      MGAM2=S
      RETURN
      END IF
10    CONTINUE
      ELSE
      S=1.0/X
      P0=0.0
      P1=1.0
      Q0=1.0
      Q1=X
      DO 20 N=1,100
      P0=P1+(N-A)*P0
      Q0=Q1+(N-A)*Q0
      P=X*PO+N*P1
      Q=X*Q0+N*Q1
      IF (ABS(Q)+1.0.NE.1.0)THEN
      S1=P/Q
      P1=P
      Q1=Q
      IF (ABS((S1-S)/S1).LT.1.0D-07)THEN
      S=S1*EXP(-X)*QQ/MGAM1(A)
      MGAM2=1.0-S
      RETURN
      END IF
      S=S1
      END IF
      P1=P
      Q1=Q
20    CONTINUE
      END IF
      WRITE (*,*) 'A TOO LARGE !'
      S=1.0-S*EXP(-X)*QQ/MGAM1(A)
      MGAM2=S
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC---error function---CCCCCCCCCCCCCCCCCCCCCC
      FUNCTION MERRF(X)
      DOUBLE PRECISION MERRF,X,MGAM2
      IF (X.GE.0) THEN
      MERRF=MGAM2(0.5D0,X*X)
      ELSE
      MERRF=-MGAM2(0.5D0,X*X)
      END IF
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC---normal distribution function---CCCCCCCC
      FUNCTION MGASS(A,D,X)
      DOUBLE PRECISION MGASS,A,D,X,MERRF
      IF (D.LE.0) D=1.0D-10
      MGASS=0.5+0.5*MERRF((X-A)/(SQRT(2.0)*D))
      RETURN
      END
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC---Generate artificial earthquake catalog---CCCCC
      SUBROUTINE ART_EC(X,Y,Y2,TEMP,MAG_MAX,MAG_MIN,XIGMA,MIU,BB2,ML)
      INTEGER TEMP,ML,I,CS
      DIMENSION X(TEMP),Y(TEMP),Y2(TEMP)
      DOUBLE PRECISION PP2,X,Y,Y2,XIGMA,MIU,BB2,MAG_MAX,MAG_MIN,A_X,
     1                 S_CDF,B_X,XS2,XS
      
C     计算观测的地震概率密度函数
      A_X=-100
      IF (ML.EQ.1) THEN
      DO I=1,TEMP
      B_X=X(I)
      CALL CDFX(A_X,B_X,10000,XIGMA,MIU,S_CDF)
      Y2(I)=S_CDF*EXP(-BB2*X(I)*2.3)
      END DO
      
      DO I=1,TEMP
      Y(I)=SUM(Y2(1:I))
      END DO    
      Y(1:TEMP)=Y(1:TEMP)/MAXVAL(Y(1:TEMP))   ! Cumulative probability density
      
      DO 516 I=1,TEMP
      IF (I.EQ.1) THEN
      Y2(I)=Y(I)
      GOTO 516
      END IF
      Y2(I)=Y(I)-Y(I-1)   ! Probability density
516   CONTINUE
      
      ELSE
      CS=0
      DO I=1,TEMP    
      B_X=X(I)
      CALL CDFX(A_X,B_X,10000,XIGMA,MIU,S_CDF)
      IF (X(I).GT.MIU) THEN
      Y2(I)=EXP(-BB2*X(I)*2.3)
      XS1=I
      ELSE
      IF (CS.EQ.0) THEN
      XS2=EXP(-BB2*X(I)*2.3)
      END IF
      CS=CS+1
      Y2(I)=S_CDF
      END IF
      END DO
      
      XS=Y2(XS1+1)/XS2
      Y2(1:XS1)=Y2(1:XS1)*XS

      DO I=1,TEMP
      Y(I)=SUM(Y2(1:I))
      END DO    
      Y(1:TEMP)=Y(1:TEMP)/MAXVAL(Y(1:TEMP))   ! Cumulative probability density
      
      DO 519 I=1,TEMP
      IF (I.EQ.1) THEN
      Y2(I)=Y(I)
      GOTO 519
      END IF
      Y2(I)=Y(I)-Y(I-1)   ! Probability density
519   CONTINUE
      END IF
      
      RETURN
      END 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--Arrange array--CCCCCCCCCCCCCCCCC
      SUBROUTINE AR_AR(CATLOG,SZ,CUMULATIVE,TEMP,MAXN1)
      INTEGER A1,A2,SZ,QQ,II,II3,II4,JJ,QQ1,CCON,TEMP
      DIMENSION CATLOG(SZ),FREQUENCY(100,2),CUMULATIVE(100,4)
      DOUBLE PRECISION CATLOG,FREQUENCY,AA,PP,CUMULATIVE,MAXN1
      
      FREQUENCY=0
      CUMULATIVE=0
      
      DO A1=1,SZ-1
      DO A2=A1+1,SZ
      IF(CATLOG(A1).GT.CATLOG(A2)) THEN
          AA=CATLOG(A1)
          CATLOG(A1)=CATLOG(A2)
          CATLOG(A2)=AA
      END IF
      END DO
      END DO
      
C--------
      QQ=0
      DO II=1,SZ-1
      IF(ABS(CATLOG(II)-CATLOG(II+1)).LT.0.01) THEN
          CYCLE
      ELSE
          QQ=QQ+1
          FREQUENCY(QQ,1)=CATLOG(II)
      END IF
      END DO
      QQ1=QQ+1
      FREQUENCY(QQ1,1)=CATLOG(SZ)
      
C--
      CCON=0
      DO II3=1,QQ1
      DO II4=1,SZ
      IF(ABS(FREQUENCY(II3,1)-CATLOG(II4)).LT.0.01) THEN
              CCON=CCON+1
      END IF
      END DO
      FREQUENCY(II3,2)=CCON
      CCON=0
      END DO
      
C    
      OPEN (214,FILE='FREQUENCY.txt')
      DO 74 II3=1,QQ1
      WRITE (214,'(F10.2,1X,F10.1)') FREQUENCY(II3,1),FREQUENCY(II3,2)
74    CONTINUE
      CLOSE (214)
      
C-----
      TEMP=0 
      DO PP=MAXN1,0.0,-0.1
      TEMP=TEMP+1
      CUMULATIVE(TEMP,1)=PP
      DO JJ=1,QQ1
      IF (ABS(PP-FREQUENCY(JJ,1)).LT.0.1) THEN
          CUMULATIVE(TEMP,2)=FREQUENCY(JJ,2)
      END IF
      END DO
      CUMULATIVE(TEMP,3)=CUMULATIVE(TEMP,1)*CUMULATIVE(TEMP,2)
      CUMULATIVE(TEMP,4)=SUM(CUMULATIVE(1:TEMP,2))
      END DO
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCC---QSORT---CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE QSORT(X,N,M,BZ)
      INTEGER N,M,BZ,I,I1
      DOUBLE PRECISION X,TEMP
      DIMENSION X(N,M),TEMP(M)
      
      DO I=1,N-1
      DO I1=I,N
      IF(X(I,BZ).GT.X(I1,BZ)) THEN
      DO I2=1,M
      TEMP(I2)=X(I,I2)
      X(I,I2)=X(I1,I2)
      X(I1,I2)=TEMP(I2)
      END DO
      END IF
      END DO   
      END DO
      
      RETURN
      END
