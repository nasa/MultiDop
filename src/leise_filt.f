      SUBROUTINE T5FLTR(Y,N1,N2,N3,NSTEP)
C                                               JIM LEISE 8/80
C    *****************************************************************
C    HELLO,
C    I AM A MULTIDIMENSIONAL LOW-PASS FILTER WHICH NEEDS NO EXTRA
C    ARRAY SPACE.  THUS, THE FILTERED ANSWER IS RETURNED IN THE SAME
C    ARRAY Y(N1,N2,N3) THAT THE DATA IS INPUT.  THE CENTRAL FILTER
C    IS A LINEAR 5-PT FILTER AND THE BOUNDARY FILTER IS COMPUTED
C    USING A MIRROR EXTENSION.  THUS, THE TOTAL FILTER IS LINEAR.
C
C          ********** NSTEP CONTROL FOR 1-DIM **********
C        STEP RESTRICTION:  5*2**(NSTEP-1) .LE. MAX(N1,N2,N3)
C         PASSBAND .LE. 2**(NSTEP+2)  POINTS/CYCLE
C         STOPBAND .GE. 2**(NSTEP)  POINTS/CYCLE.
C
C          ********** MULTIDIMENSIONAL USE **********
C    PARAMETER CONTROL FOR THE THREE DIMENSIONS CAN BE REALIZED
C    VIA COMMON/FLTRPL/ WHERE NS CORRESPONDS TO NSTEP.  IF THIS
C    COMMON IS NOT USED, THE VALUES OF NS ARE DEFAULTED TO NSTEP
C    -I.E. NSTEP IS USED IN PLACE OF ANY ZEROS.
C    ******************************************************************
C
         DIMENSION Y(1),KORD(5),NET(5),NNS(3)
         COMMON/FLTRPL/NS(3)
C    INITIALIZATION OF NS FOR CSD APPLICATIONS (4/15/82)
         DATA NS/0,0,0/
C
C    Definition of input variables
C    Y      - input variable to filter
C    N1     - Number of points in 1st dimension
C    N2     - Number of points in 2nd dimension
C    N3     - Number of points in 3rd dimension
C    NSTEP  - Number of steps to perform
C    NDIM   - Number of dimensions
C
C    INITIALIZE THE 3-D ARITHMETIC.
C    Determines how many spatial dimensions there are, defaults to 1 for N2 and N3 do not exist
         NDIM=1
         IF(N2.GT.1)NDIM=2
         IF(N3.GT.1)NDIM=3
C    Sets KORD of up to 5 dimensions, the first three are spatial 
         KORD(1)=MAX0(1,N1)
         KORD(2)=MAX0(1,N2)
         KORD(3)=MAX0(1,N3)
         KORD(4)=KORD(1)
         KORD(5)=KORD(2)
         NET(1)=1
         NET(2)=KORD(1)
         NET(3)=KORD(1)*KORD(2)
         NET(4)=NET(1)
         NET(5)=NET(2)
C
C    DEFAULT PARAMETER TRANSFER.
         MPYRMD=0
C        Loop until line 10 over each dimension
         DO 10 N=1,NDIM
         NNS(N)=NS(N)
         IF(NS(N).EQ.0)NNS(N)=NSTEP
         IF(KORD(N).LT.5)NNS(N)=0
 10      MPYRMD=MAX0(MPYRMD,NNS(N)+NNS(N)-1)
         IF(MPYRMD.LE.0)RETURN
         MSTEP=(MPYRMD+1)/2
C
C    ***** START THE MAIN LOOP *****
         K1=1
         DO 50 MAIN=1,MPYRMD
         DO 40 N=1,NDIM
C    SAMPLING CHECKS.
         IF(10*K1.GT.KORD(N))NNS(N)=MIN0(NNS(N),MAIN)
         IF((MAIN.GE.NNS(N)).AND.(MPYRMD-MAIN.GE.NNS(N)))GO TO 40
C
C    THE 3-D ARITHMETIC.
         M1=K1*NET(N)
         M2=M1+M1
         M3=M2+M1
         ISTOP=KORD(N+1)
         JSTOP=KORD(N+2)
         DO 30 I=1,ISTOP
         DO 30 J=1,JSTOP
         KSTRT=1+(I-1)*NET(N+1)+(J-1)*NET(N+2)
         KSTOP=KSTRT+(KORD(N)-1)*NET(N)
         KN=KSTRT-NET(N)
         DO 30 K=1,K1
         KN=KN+NET(N)
C        WRITE(*,*) ' I, J, K K1', I, J, K, K1
         LN=KN+((KSTOP-KN)/M1)*M1
C
C    FILTER THE ENDS USING A MIRROR EXTENSION.
         YKN=.875*Y(KN)+.1875*Y(KN+M1)-.0625*Y(KN+M2)
         YLN=.875*Y(LN)+.1875*Y(LN-M1)-.0625*Y(LN-M2)
         YKN1=.1875*Y(KN)+.625*Y(KN+M1)+.25*Y(KN+M2)-.0625*Y(KN+M3)
         YLN1=.1875*Y(LN)+.625*Y(LN-M1)+.25*Y(LN-M2)-.0625*Y(LN-M3)
C
C    DO THE CENTRAL 5-PT FILTER.
         YM2=Y(KN)
         YM1=Y(KN+M1)
         MSTRT=KN+M2
         MSTOP=LN-M2
C
         DO 20 M=MSTRT,MSTOP,M1
         YSAVE=Y(M)
         Y(M)=.625*Y(M)+.25*(YM1+Y(M+M1))-.0625*(YM2+Y(M+M2))
         YM2=YM1
 20      YM1=YSAVE
C
         Y(KN+M1)=YKN1
         Y(LN-M1)=YLN1
         Y(KN)=YKN
         Y(LN)=YLN
C
 30      CONTINUE
 40      CONTINUE
C    UPDATE THE SAMPLING INCREMENT.
         K1=K1+K1
         IF(MAIN.GE.MSTEP)K1=K1/4
 50      CONTINUE
C
         RETURN
      END
