! FROM UPDATE LIBRARY         LMDBIB                CY=2       02/08/85

!************************************************************************
!*                                                                      *
!* C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                      *
!*                                                               FFT99  *
!*                                                               FFT991 *
!*                                                                      *
!*                                                                      *
!* SUBPROGRAM       SUBROUTINE    FFT99                                 *
!*                                FFT991                                *
!*                                                                      *
!* PURPOSE          PERFORM MULTIPLE FAST FOURIER TRANSFORMS            *
!*                                                                      *
!*                                                                      *
!* VERSION          CYBER                         CRAY-1                *
!*                                                                      *
!*                  JAN 1979 ORIGINAL             JAN 1979 ORIGINAL     *
!*                                                                      *
!* USAGE                                                                *
!*                  CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                  CALL FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                                                                      *
!* ARGUMENTS        1.DIMENSION                                         *
!*                       A(IDIM),WORK((N+1)*M),TRIGS(3*N/2),IFAX(10)    *
!*                       WORK IS A WORK ARRAY                           *
!*             ! >>>     TRIGS(DIMENSION 3*N/2 - ADD 1 IF N/2 IS ODD)   *
!*                                                                      *
!*                  2.INPUT                                             *
!*                      A - AN ARRAY CONTAINING THE INPUT DATA OR       *
!*                          COEFFICIENT VECTORS.                        *
!*                         THIS ARRAY IS OVERWRITTEN BY THE RESULTS.    *
!*                      TRIGS AND IFAX - ARRAYS SET UP BY FFTRIG AND FAX*
!*                                     - SEE WRITEUP OF FFTRIG AND FAX  *
!*                      INC - THE WORD INCREMENT BETWEEN SUCCESSIVE     *
!*                           ELEMENTS OF EACH DATA OR COEFFICIENT VECTOR*
!*                           E.G. INC=1 FOR CONSECUTIVELY STORED DATA.  *
!*                      JUMP - THE WORD INCREMENT BETWEEN THE FIRST     *
!*                            ELEMENTS OF SUCCESSIVE DATA OR COEFFICIENT*
!*                            VECTORS.                                  *
!*                      N - THE LENGTH OF EACH TRANSFORM. (SEE NOTE X)  *
!*                      M - THE NUMBER OF TRANSFORMS TO BE DONE         *
!*                          SIMULTANEOUSLY.                             *
!*                      ISIGN - +1 FOR A TRANSFORM FROM FOURIER         *
!*                              COEFFICIENTS TO DATA VALUES.            *
!*                              -1 FOR A TRANSFORM FROM DATA VALUES     *
!*                              TO FOURIER COEFFICIENTS.                *
!*                                                                      *
!*                  3.OUTPUT                                            *
!*                      A - CONTAINS EITHER THE COEFFICIENTS OR THE     *
!*                          DATA VALUES,DEPENDING ON ISIGN.             *
!*                          IN EACH CASE N INDEPENDENT QUANTITIES       *
!*                          OCCUPY N+2 WORDS.   THE COEFFICIENTS ARE    *
!*                          STORED AS SUCCESSIVE PAIRS OF REAL AND      *
!*                          IMAGINARY PARTS -                           *
!*                          A(K),B(K) , K=0,1,...N/2                    *
!*                          B(0) AND B(N/2) ARE STORED ALTHOUGH THEY    *
!*                          MUST BE 0.                                  *
!*                      FOR FFT99 THE DATA IS STORED WITH EXPLICIT      *
!*                          PERIODICITY -                               *
!*                          X(N-1),X(0),X(1),....X(N-1),X(0)            *
!*                      FOR FFT991 THE DATA APPEARS AS -                *
!*                          X(0),X(1),X(2),......X(N-1),0,0             *
!*                                                                      *
!* NOTES            1. ON CRAY-1, ARRANGE DATA SO THAT JUMP IS NOT A    *
!*                     MULTIPLE OF 8 (TO AVOID MEMORY BANK CONFLICTS)   *
!*                                                                      *
!* WRITE UP         COMPUTER BULLETIN B6.6/1                            *
!*                                                                      *
!* ENTRY POINTS        FFT99,FFT991                                     *
!*                                                                      *
!* COMMON BLOCKS    NONE                                                *
!*                                                                      *
!* I/O              NONE                                                *
!*                                                                      *
!* PRECISION        SINGLE                                              *
!*                                                                      *
!* OTHER ROUTINES   FFT99A,FFT99B,VPASSM          (CY)                  *
!*       REQUIRED   CAL99,CPASS                   (CR)                  *
!*                                                                      *
!*                                                                      *
!* 7/80                      FFT99-1                                    *
!*                                                                      *
!************************************************************************
!*                                                                      *
!* C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                      *
!*                                                               FFT99  *
!*                                                               FFT991 *
!*                                                                      *
!* ACCESS (OBJECT)  CYBER:                                              *
!*                           ATTACH,ECLIB.                              *
!*                           LDSET(LIB=ECLIB)                           *
!*                  CRAY 1:                                             *
!*                           LDR(LIB=ECLIB...)                          *
!*                                                                      *
!* ACCESS (SOURCE)           ATTACH,OLDPL,ECLIBPL                       *
!*                                                                      *
!*                  CYBER :         %DEFINE CYBER                       *
!*                  CRAY:           %DEFINE CRAY                        *
!*                                  %C    FFT99,FFT991                  *
!*                                                                      *
!* LANGUAGE         FORTRAN                                             *
!*                  BUT CRAY IMPLEMENTATION OF PASS IS IN CAL           *
!*                                                                      *
!* SPECIALIST       CLIVE TEMPERTON                                     *
!*                                                                      *
!* HISTORY          WRITTEN BY C.TEMPERTON      JAN     1979            *
!*                                                                      *
!* ALGORITHM        THE ALGORITHM IS THE SELF-SORTING (TEMPERTON)       *
!*                  VERSION OF THE FAST FOURIER TRANSFORM               *
!*                                                                      *
!* REFERENCES       ECMWF TECHNICAL REPORT NO.3                         *
!*                  ECMWF INTERNAL REPORT NO.21 -   C.TEMPERTON         *
!*                                                                      *
!* OBJECT SIZE               FFT991  FFT99  (OCTAL WORDS)               *
!*                  CYBER:    2665    2676                              *
!*                  CRAY :    1250    1260                              *
!*                                                                      *
!*                                                                      *
!* ACCURACY                                                             *
!*                                                                      *
!* TIMING           VECTORIZATION IS ON VECTORS OF LENGTH M.      (CR)  *
!*                  HENCE TIMING IS STRONGLY DEPENDENT ON M.            *
!*                  TIME PER TRANSFORM ON CRAY-1 (MICROSECONDS)         *
!*                  N    M=4    M=16    M=64                            *
!*                 64     46      17      10                            *
!*                128     81      33      21                            *
!*                180    150      58      37                            *
!*                192    149      58      36                            *
!*                240    192      76      49                            *
!*                256    191      76      49                            *
!*                288    219      89      58                            *
!*                300    253     102      68                            *
!*                320    248     101      66                            *
!*                360    286     118      79                            *
!*               1024    898     359     238                            *
!*                                                                      *
!* PORTABILITY      STANDARD FORTRAN                                    *
!*                  STANDARD CAL  (CR)                                  *
!*                                                                      *
!* SYSTEM ROUTINES  NONE                                                *
!*        REQUIRED                                                      *
!*                                                                      *
!* 7/80                      FFT99-1                                    *
!*                                                                      *
!************************************************************************
!
!     SUBROUTINE 'FFT99' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS ' WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
!
!     A IS THE ARRAY CONTAINING INPUT ' OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA "VECTOR"
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
!         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
!
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
SUBROUTINE FFT99 (A, WORK, TRIGS, IFAX, INC, JUMP, N, LOT, ISIGN)
  USE SHARE_VARS
  IMPLICIT NONE
  INTEGER N
  REAL (KIND=PR), DIMENSION (N) :: A, WORK, TRIGS
  INTEGER, DIMENSION (1) :: IFAX
  INTEGER :: INC, ISIGN, JUMP, LOT
  INTEGER :: NFAX, MX, NH, INK, IGO, IBASE, JBASE, I, J, L, M
  INTEGER :: K, IA, IB, LA
      NFAX=IFAX(1)
      MX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+MX
   20 CONTINUE
!
      IGO=60
      GO TO 40
!
!     PREPROCESSING (ISIGN=+1)
!     ------------------------
!
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!
!     COMPLEX TRANSFORM
!     -----------------
!
   40 CONTINUE
      IA=INC+1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS, &
        INK,2,JUMP,MX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS, &
         2,INK,MX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!
      IF (ISIGN.EQ.-1) GO TO 130
!
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+MX
      JBASE=JBASE+JUMP
  100 CONTINUE
!
!     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1
      DO 120 L=1,LOT
      A(IA)=A(IB)
      A(IB+INC)=A(IA+INC)
      IA=IA+JUMP
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!
!     POSTPROCESSING (ISIGN=-1):
!     --------------------------
!
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
  140 CONTINUE
END SUBROUTINE FFT99

SUBROUTINE FFT99A (A, WORK, TRIGS, INC, JUMP, N, LOT)
!     subroutine fft99a - preprocessing step for fft99, isign=+1
!     (spectral to gridpoint transform)
  USE SHARE_VARS
  IMPLICIT NONE
  INTEGER N
  REAL (KIND=PR), DIMENSION (N) :: A, WORK, TRIGS
  REAL (KIND=PR) :: C, S
  INTEGER :: INC, JUMP, LOT
  INTEGER :: NH, MX, INK, IA, IB, JA, JB, K, L
  INTEGER :: JBASE, IABASE, IBBASE, JABASE, JBBASE
      NH=N/2
      MX=N+1
      INK=INC+INC
!
!     A(0) ' A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+MX
      JB=JB+MX
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
!
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))- &
         (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+ &
         (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+ &
         (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))- &
         (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+MX
      JB=JB+MX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      DO 40 L=1,LOT
      WORK(JA)=2.0*A(IA)
      WORK(JA+1)=-2.0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+MX
   40 CONTINUE
!
   50 CONTINUE
end subroutine fft99a

SUBROUTINE FFT99B (WORK, A, TRIGS, INC, JUMP, N, LOT)
!     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN = -1
!     (GRIDPOINT TO SPECTRAL TRANSFORM)
  USE SHARE_VARS
  IMPLICIT NONE
  INTEGER N
  REAL (KIND=PR), DIMENSION (N) :: WORK, A, TRIGS
  REAL (KIND=PR) :: C,S, SCALE
  INTEGER :: INC, JUMP, LOT
  INTEGER :: NH, MX, INK, IA, IB, JA, JB, IABASE, IBBASE, JABASE, JBBASE
  INTEGER :: K, L

      NH=N/2
      MX=N+1
      INK=INC+INC
!
!     A(0) ' A(N/2)
      SCALE=1.0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.0
      A(JB+INC)=0.0
      IA=IA+MX
      IB=IB+MX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      SCALE=0.5*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
!
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB)) &
        +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB)) &
        -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1))) &
         +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1))) &
         -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+MX
      IB=IB+MX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0*SCALE
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+MX
      JA=JA+JUMP
   40 CONTINUE
!
   50 CONTINUE
END SUBROUTINE FFT99B

!     SUBROUTINE 'FFT991' - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM
!
!************************************************************************
!*                                                                     *
!*C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                     *
!*                                                              FFT99  *
!*                                                              FFT991 *
!*                                                                     *
!*                                                                     *
!*SUBPROGRAM       SUBROUTINE    FFT99                                 *
!*                               FFT991                                *
!*                                                                     *
!*PURPOSE          PERFORM MULTIPLE FAST FOURIER TRANSFORMS            *
!*                                                                     *
!*                                                                     *
!*VERSION          CYBER                         CRAY-1                *
!*                                                                     *
!*                 JAN 1979 ORIGINAL             JAN 1979 ORIGINAL     *
!*                                                                     *
!*USAGE                                                                *
!*                 CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                 CALL FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   *
!*                                                                     *
!*ARGUMENTS        1.DIMENSION                                         *
!*                      A(IDIM),WORK((N+1)*M),TRIGS(3*N/2),IFAX(10)    *
!*                      WORK IS A WORK ARRAY                           *
!*                                                                     *
!*                 2.INPUT                                             *
!*                     A - AN ARRAY CONTAINING THE INPUT DATA OR       *
!*                         COEFFICIENT VECTORS.                        *
!*                        THIS ARRAY IS OVERWRITTEN BY THE RESULTS.    *
!*                     TRIGS AND IFAX - ARRAYS SET UP BY FFTRIG AND FAX*
!*                                    - SEE WRITEUP OF FFTRIG AND FAX  *
!*                     INC - THE WORD INCREMENT BETWEEN SUCCESSIVE     *
!*                          ELEMENTS OF EACH DATA OR COEFFICIENT VECTOR*
!*                          E.G. INC = 1 FOR CONSECUTIVELY STORED DATA.  *
!*                     JUMP - THE WORD INCREMENT BETWEEN THE FIRST     *
!*                           ELEMENTS OF SUCCESSIVE DATA OR COEFFICIENT*
!*                           VECTORS.                                  *
!*                     N - THE LENGTH OF EACH TRANSFORM. (SEE NOTE X)  *
!*                     M - THE NUMBER OF TRANSFORMS TO BE DONE         *
!*                         SIMULTANEOUSLY.                             *
!*                     ISIGN - +1 FOR A TRANSFORM FROM FOURIER         *
!*                             COEFFICIENTS TO DATA VALUES.            *
!*                             -1 FOR A TRANSFORM FROM DATA VALUES     *
!*                             TO FOURIER COEFFICIENTS.                *
!*                                                                     *
!*                 3.OUTPUT                                            *
!*                     A - CONTAINS EITHER THE COEFFICIENTS OR THE     *
!*                         DATA VALUES,DEPENDING ON ISIGN.             *
!*                         IN EACH CASE N INDEPENDENT QUANTITIES       *
!*                         OCCUPY N+2 WORDS.   THE COEFFICIENTS ARE    *
!*                         STORED AS SUCCESSIVE PAIRS OF REAL AND      *
!*                         IMAGINARY PARTS -                           *
!*                         A(K),B(K) , K = 0,1,...N/2                    *
!*                         B(0) AND B(N/2) ARE STORED ALTHOUGH THEY    *
!*                         MUST BE 0.                                  *
!*                     FOR FFT99 THE DATA IS STORED WITH EXPLICIT      *
!*                         PERIODICITY -                               *
!*                         X(N-1),X(0),X(1),....X(N-1),X(0)            *
!*                     FOR FFT991 THE DATA APPEARS AS -                *
!*                         X(0),X(1),X(2),......X(N-1),0,0             *
!*                                                                     *
!*NOTES            1. ON CRAY-1, ARRANGE DATA SO THAT JUMP IS NOT A    *
!*                    MULTIPLE OF 8 (TO AVOID MEMORY BANK CONFLICTS)   *
!*                                                                     *
!*WRITE UP         COMPUTER BULLETIN B6.6/1                            *
!*                                                                     *
!*ENTRY POINTS        FFT99,FFT991                                     *
!*                                                                     *
!*COMMON BLOCKS    NONE                                                *
!*                                                                     *
!*I/O              NONE                                                *
!*                                                                     *
!*PRECISION        SINGLE                                              *
!*                                                                     *
!*OTHER ROUTINES   FFT99A,FFT99B,VPASSM          (CY)                  *
!*      REQUIRED   CAL99,CPASS                   (CR)                  *
!*                                                                     *
!*                                                                     *
!*7/80                      FFT99-1                                    *
!*                                                                     *
!************************************************************************
!*                                                                     *
!*C06-SUMMATION OF SERIES                                       B6.1/3 *
!*                                                                     *
!*                                                              FFT99  *
!*                                                              FFT991 *
!*                                                                     *
!*ACCESS (OBJECT)  CYBER:                                              *
!*                          ATTACH,ECLIB.                              *
!*                          LDSET(LIB = ECLIB)                           *
!*                 CRAY 1:                                             *
!*                          LDR(LIB = ECLIB...)                          *
!*                                                                     *
!*ACCESS (SOURCE)           ATTACH,OLDPL,ECLIBPL                       *
!*                                                                     *
!*                 CYBER :         %DEFINE CYBER                       *
!*                 CRAY:           %DEFINE CRAY                        *
!*                                 %C    FFT99,FFT991                  *
!*                                                                     *
!*LANGUAGE         FORTRAN                                             *
!*                 BUT CRAY IMPLEMENTATION OF PASS IS IN CAL           *
!*                                                                     *
!*SPECIALIST       CLIVE TEMPERTON                                     *
!*                                                                     *
!*HISTORY          WRITTEN BY C.TEMPERTON      JAN     1979            *
!*                                                                     *
!*ALGORITHM        THE ALGORITHM IS THE SELF-SORTING (TEMPERTON)       *
!*                 VERSION OF THE FAST FOURIER TRANSFORM               *
!*                                                                     *
!*REFERENCES       ECMWF TECHNICAL REPORT NO.3                         *
!*                 ECMWF INTERNAL REPORT NO.21 -   C.TEMPERTON         *
!*                                                                     *
!*OBJECT SIZE               FFT991  FFT99  (OCTAL WORDS)               *
!*                 CYBER:    2665    2676                              *
!*                 CRAY :    1250    1260                              *
!*                                                                     *
!*                                                                     *
!*ACCURACY                                                             *
!*                                                                     *
!*TIMING           VECTORIZATION IS ON VECTORS OF LENGTH M.      (CR)  *
!*                 HENCE TIMING IS STRONGLY DEPENDENT ON M.            *
!*                 TIME PER TRANSFORM ON CRAY-1 (MICROSECONDS)         *
!*                 N    M=4    M=16    M=64                            *
!*                64     46      17      10                            *
!*               128     81      33      21                            *
!*               180    150      58      37                            *
!*               192    149      58      36                            *
!*               240    192      76      49                            *
!*               256    191      76      49                            *
!*               288    219      89      58                            *
!*               300    253     102      68                            *
!*               320    248     101      66                            *
!*               360    286     118      79                            *
!*              1024    898     359     238                            *
!*                                                                     *
!*PORTABILITY      STANDARD FORTRAN                                    *
!*                 STANDARD CAL  (CR)                                  *
!*                                                                     *
!*SYSTEM ROUTINES  NONE                                                *
!*       REQUIRED                                                      *
!*                                                                     *
!*7/80                      FFT99-1                                    *
!*                                                                     *
!************************************************************************
!
!     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!     THAT IN MRFFT2
!
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS ' WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
!
!     A IS THE ARRAY CONTAINING INPUT ' OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA "VECTOR"
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
!
!     **!*N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
SUBROUTINE FFT991 (A, WORK, TRIGS, IFAX, INC, JUMP, N, LOT, ISIGN)
  USE SHARE_VARS
  IMPLICIT NONE
  INTEGER N
  REAL (KIND=PR), DIMENSION (N) :: A, WORK, TRIGS
  INTEGER, DIMENSION (1) :: IFAX
  INTEGER :: INC, JUMP, LOT, ISIGN
  INTEGER :: NFAX, MX, NH, INK, IGO, IBASE, JBASE, I, J, K, L, M, IA, IB, LA 
      NFAX=IFAX(1)
      MX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+MX
   20 CONTINUE
!
      IGO=60
      GO TO 40
!
!     PREPROCESSING (ISIGN=+1)
!     ------------------------
!
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!
!     COMPLEX TRANSFORM
!     -----------------
!
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS, &
        INK,2,JUMP,MX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS, &
         2,INK,MX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!
      IF (ISIGN.EQ.-1) GO TO 130
!
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+MX
      JBASE=JBASE+JUMP
  100 CONTINUE
!
!     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
      DO 120 L=1,LOT
      A(IB)=0.0
      A(IB+INC)=0.0
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!
!     POSTPROCESSING (ISIGN=-1):
!     --------------------------
!
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
  140 CONTINUE
END SUBROUTINE FFT991

!     subroutine 'vpassm' - multiple version of 'vpassa'
!     performs one pass through data
!     as part of multiple complex fft routine
!     a is first real input vector
!     b is first imaginary input vector
!     c is first real output vector
!     d is first imaginary output vector
!     trigs is precalculated table of sines ' cosines
!     inc1 is addressing increment for a and b
!     inc2 is addressing increment for c and d
!     inc3 is addressing increment between a's & b's
!     inc4 is addressing increment between c's & d's
!     lot is the number of vectors
!     n is length of vectors
!     ifac is current factor of n
!     la is product of previous factors
SUBROUTINE VPASSM (A, B, C, D, TRIGS, INC1, INC2, INC3, INC4, LOT, N, IFAC, LA)
  USE SHARE_VARS
  IMPLICIT NONE
  INTEGER N
  REAL (KIND=PR), DIMENSION (N) :: A, B, C, D, TRIGS
  REAL (KIND=PR), PARAMETER  :: SIN36 = 0.587785252292473, &
       COS36 = 0.809016994374947, &
       SIN72 = 0.951056516295154, COS72 = 0.309016994374947, &
       SIN60 = 0.866025403784437
  REAL (KIND=PR) :: C1, C2, C3, C4, S1, S2, S3, S4
  INTEGER :: INC1, INC2, INC3, INC4, LOT, IFAC, LA
  INTEGER :: M, IINK, JINK, JUMP, IGO, IA, IB, JA, JB, IBASE, JBASE, I, J, K, L, IJK
  INTEGER :: IC, JC, ID, JD, IE, JE, KE, LA1, KB, KC, KD
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
!     CHECK FACTORS ARE CORRECT - ENSURE NON-NEGATIVE
      GO TO (10,50,90,130),IGO
!
!     CODING FOR FACTOR 2
!
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 3
!
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)= &
         C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))) &
        -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)= &
         S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))) &
        +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)= &
         C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))) &
        -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)= &
         S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))) &
        +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 4
!
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)= &
         C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
        -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)= &
         S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
        +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)= &
         C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))) &
        -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)= &
         S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))) &
        +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)= &
         C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))) &
        -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)= &
         S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))) &
        +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 5
!
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
       -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
       +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
       +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
       -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
       -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
       +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
       +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
       -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)= &
         C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
           -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
        -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
           +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)= &
         S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
           -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
        +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
           +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)= &
         C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
           +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
        -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
           -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)= &
         S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
           +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
        +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
           -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))) 
      C(JC+J)= &
         C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
           -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
        -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
           +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)= &
         S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
           -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
        +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
           +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)= &
         C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
           +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
        -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
           -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)= &
         S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
           +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
        +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
           -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
END SUBROUTINE VPASSM


!***********************************************************************
!                                                                      *
! C06-SUMMATIOM
! C06-SUMMATION OF SERIES                                      B6.1/3  *
!                                                                      *
!                                                              FFTRIG  *
!                                                              FAX     *
!                                                                      *
!                                                                      *
! SUP
! SUBPROGRAM       SUBROUTINE   FFTRIG                                 *
!                               FAX                                    *
!                                                                      *
! PURPOSE          SETUP ROUTINES FOR FFT PACKAGES                     *
!                                                                      *
!                                                                      *
! VERSION          CYBER                         CRAY-1                *
!                                                                      *
!                  JAN 1979 ORIGINAL             JAN 1979 ORIGINAL     *
!                                                                      *
! USAGE                                                                *
!                  CALL FFTRIG(TRIGS,N,3)                              *
!                  CALL FAX   (IFAX ,N,3)                              *
!                                                                      *
! ARGUMENTS        1.DIMENSION                                         *
!                       TRIGS(DIMENSION 3*N/2 - ADD 1 IF N/2 IS ODD)   *
!                       IFAX(10)                                       *
!                                                                      *
!                  2.INPUT                                             *
!                      N - THE LENGHT OF THE TRANSFORMS TO BE PERFORMED*
!                          N MUST BE EVEN.                             *
!                          THE NUMBER OF WORDS OF IFAX USED INCREASES  *
!                          LOGARITHMICALLY WITH N.                     *
!                          IFAX(10) SUFFICES FOR PRACTICAL PURPOSES.   *
!                          (TRANSFORMS OF LENGHT AT LEAST 10000)       *
!                                                                      *
!                  3.OUTPUT                                            *
!                      TRIGS - FFTRIG RETURNS AN ARRAY OF TRIGONOMETRIC*
!                              FUNCTION VALUES SUBSEQUENTLY USED BY    *
!                              FFT ROUTINES.                           *
!                      IFAX  - FAX FACTORIZES N/2 INTO A PRODUCT OF    *
!                              4"S AND 2"S AND HIGHER PRIME NUMBERS.   *
!                              IFAX(1) CONTAINS THE NUMBER OF FACTORS. *
!                              AND THE FACTORS THEMSELVES ARE STORED   *
!                              IN ASCENDING ORDER IN IFAX(2),IFAX(3).. *
!                              IF FAX IS CALLED WITH N ODD ,IFAX(1)    *
!                              IS SET TO -99(ERROR CONDITION) AND NO   *
!                              FACTORIZATION IS DONE.                  *
!                                                                      *
! WRITE UP         NONE                                                *
!                                                                      *
! ENTRY POINTS           FFTRIG,  FAX                                  *
!                                                                      *
! COMMON BLOCKS    NONE                                                *
! I/O              NONE                                                *
! PRECISION        SINGLE                                              *
! OTHER ROUTINES   NONE                                                *
!       REQUIRED                                                       *
! 7/80                     FFTRIG-1                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
! CO6-SUMMATION OF SERIES                                       B6.1/3 *
!                                                                      *
!                                                              FFTRIG  *
!                                                              FAX     *
!                                                                      *
! ACSSES (OBJECT)  CYBER:                                              *
!                           ATTACH,ECLIB.                              *
!                           LDSET(LIB=ECLIB)                           *
!                  CRAY 1:                                             *
!                           LDR(LIB=ECLIB...)                          *
!                                                                      *
! ACCESS (SOURCE)           ATTACH,OLDPL,ECLIBPL                       *
!                                                                      *
!                  CYBER :         %DEFINE CYBER                       *
!                  CRAY:           %DEFINE CRAY                        *
!                                  %C   FFTRIG,   FAX                  *
!                                                                      *
! LANGUAGE         FORTRAN                                             *
!                                                                      *
! SPECIALIST       CLIVE TEMPERTON                                     *
!                                                                      *
! HISTORY          WRITTEN BY C.TEMPERTON      JAN     1979            *
!                                                                      *
! ALGORITHM                                                            *
! REFERENCES                                                           *
!                                                                      *
! OBJECT SIZE               FFTRIG  FAX  (OCTAL WORDS)                 *
!                  CYBER:     145   127                                *
!                  CRAY :     221   157                                *
!                                                                      *
!                                                                      *
! ACCURACY                                                             *
!                                                                      *
! TIMING                                                               *
!                                                                      *
! PORTABILITY      STANDARD FORTRAN                                    *
!                                                                      *
! SYSTEM ROUTINES  NONE                                                *
!        REQUIRED                                                      *
!
! 7/80                      FFTRIG-2                                   *
!
!***********************************************************************
SUBROUTINE FAX (IFAX, N, MODE)
  IMPLICIT NONE
  INTEGER, DIMENSION (10) :: IFAX
  INTEGER :: N, MODE
  INTEGER :: NN, IAB, NN2, I, K, L, INC, ITEM, II, ISTOP, NFAX
      NN=N
      IAB=IABS(MODE)
      IF (IAB.EQ.1) GO TO 10
      IF (IAB.EQ.8) GO TO 10
      NN=N/2
      NN2=NN+NN
      IF (NN2.EQ.N)  GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
!     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40

!     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
!     INC ALTERNATIVELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
!     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
!     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
END SUBROUTINE FAX

SUBROUTINE FFTRIG (TRIGS, N, MODE)
!    FFTRIG RETURNS AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES
!           SUBSEQUENTLY USED BY    F F T    ROUTINES
!    SEE COMMENTS IN ROUTINE    F A X
  USE SHARE_VARS
  IMPLICIT NONE
  INTEGER :: N
  REAL (KIND=PR), DIMENSION (3*N+2) :: TRIGS
  REAL (KIND=PR) :: ANGLE, DEL, PI1
  INTEGER*4 MODE
  INTEGER :: IMODE, NN, L, I, NH, LA
      PI1=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI1+PI1)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2

      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(LA+I)=COS(ANGLE)
      TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) THEN
         RETURN
      END IF
      DEL=0.5*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=2.0*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5*DEL
      DO 50 I=2,N
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
 END SUBROUTINE FFTRIG

 subroutine fftcc (a, work, trigs, ifax, inc, jump, n, nft, isign)
! fft complex-->complex using temperton's real-->complex fft
! same arguments as fft991
! n must be even
   use share_vars
   implicit none
   integer n
   real (kind=pr), dimension (n) :: a, work, trigs
   integer :: ifax, inc, jump, nft, isign
   if(isign == 1) then
      call trccc (a, work, inc, jump, n, nft, isign)
      call fft991 (a, work, trigs, ifax, inc, jump, n, nft, isign)
      return
   else
      call fft991 (a, work, trigs, ifax, inc, jump, n, nft, isign)
      call trccc (a, work, inc, jump, n, nft, isign)
      return
   endif
 end subroutine fftcc

!***********************************************************************
!
! this routine is used to compute complex to complex transforms
!       using temperton's routine fft991
!
!
!
! after ft along x1, one has
!
!    sum   f(x)cos(k1x1)      sum   f(x)sin(k1x1)
!     x1                       x1
!
!--------------------------------------------------------------
! after ft along x2 using fft991 one has
!
!     sum  f(x)cos(k1x1)sin(k2x2)   sum f(x)sin(k1x1)sin(k2x2)
!      x                             x
!
!     sum  f(x)cos(k1x1)cos(k2x2)   sum f(x)sin(k1x1)cos(k2x2)
!      x                             x
!
!---------------------------------------------------------------
! what we want is
!
!     sum  f(x)( cos(k1x1)cos(k2x2)-sin(k1x1)sin(k2x2) )
!      x
!
!     sum  f(x)( cos(k1x1)sin(k2x2)+sin(k1x1)cos(k2x2) )
!      x
!
!---------------------------------------------------------------
 SUBROUTINE TRCCC (F, W, INC, JUMP, N, NFT, ISIGN)
   USE SHARE_VARS
   IMPLICIT NONE
   REAL (KIND=PR), DIMENSION (1) :: F
   REAL (KIND=PR), DIMENSION (2,0:1) :: W
   INTEGER :: INC, JUMP, N, NFT, ISIGN
   INTEGER :: NH, NFTH, I1, INCD, IFT, I2, I1P, I2P, J, I1N, I2N 
      NH=N/2
      NFTH=NFT/2
      I1=1
      INCD=INC*2
      IF(ISIGN.EQ.-1)THEN
        DO 1 IFT=1,NFTH
        I2=I1+1
        I1P=I1+INC
        I2P=I2+INC
        DO 3 J=0,NH
        W(1,J)=F(I1+INCD*J)-F(I2P+INCD*J)
        W(2,J)=F(I2+INCD*J)+F(I1P+INCD*J)
        W(1,N-J)=F(I1+INCD*J)+F(I2P+INCD*J)
        W(2,N-J)=F(I2+INCD*J)-F(I1P+INCD*J)
  3     CONTINUE
        DO 5 J=0,N-1
        F(I1+INC*J)=W(1,J)
        F(I2+INC*J)=W(2,J)
  5     CONTINUE
        I1=I1+JUMP*2
  1     CONTINUE
      ELSE
        DO 2 IFT=1,NFTH
        I2=I1+1
        I1N=I1+N*INC
        I2N=I2+N*INC
        DO 4 J=1,NH-1
        W(1,2*J)=.5*(F(I1+INC*J)+F(I1N-INC*J))
        W(2,2*J)=.5*(F(I2+INC*J)+F(I2N-INC*J))
        W(1,2*J+1)=.5*(F(I2+INC*J)-F(I2N-INC*J))
        W(2,2*J+1)=.5*(F(I1N-INC*J)-F(I1+INC*J))
  4     CONTINUE
        W(1,0)=F(I1)
        W(2,0)=F(I2)
        W(1,1)=0.
        W(2,1)=0.
        W(1,N)=F(I1N)
        W(2,N)=F(I2N)
        W(1,N+1)=0.
        W(2,N+1)=0.
        DO 6 J=0,N+1
        F(I1+INC*J)=W(1,J)
        F(I2+INC*J)=W(2,J)
  6     CONTINUE
        I1=I1+JUMP*2
  2     CONTINUE
      ENDIF
 END SUBROUTINE TRCCC
 
 subroutine set99 (trigs, ifax, n)
   use share_vars
   implicit none
   integer :: n
   integer, dimension (1) :: ifax
   real (kind=pr), dimension (n) :: trigs
   integer :: lmode
   lmode = 3
   call fftrig (trigs, n, lmode)
   call fax (ifax, n, lmode)
 end subroutine set99
