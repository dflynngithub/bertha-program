      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*16 FILE1,FILE2,FILE3
C
c      WRITE(*,*) 'DIMENSION?'
c      READ(*,*) NDIM
       NDIM = 1804
c      WRITE(*,*) 'FILE 1 NAME?'
c      READ(*,*) FILE1
c      WRITE(*,*) 'FILE 2 NAME?'
c      READ(*,*) FILE2
c      WRITE(*,*) 'FILE 3 NAME? (OUTPUT)'
c      READ(*,*) FILE3
       FILE1 = 'NUCGRAD_r'
       FILE2 = 'NUCSCHF_r'
       FILE3 = 'VOLEFFECT'
C      
       CALL DIFF(FILE1,FILE2,FILE3,NDIM)
c
c      CALL DIAG(FILE3,NDIM)

      
      CALL GNUMTRX(FILE1,NDIM)
      CALL GNUMTRX(FILE2,NDIM)
      CALL GNUMTRX(FILE3,NDIM)

      stop

      
10    CONTINUE
      WRITE(*,*) 'Value of an element... I,J?'
      READ(*,*) I,J
C
      CALL ELEMENT(FILE3,NDIM,I,J,VAL)
      
      WRITE(*,*) 'FOCK3(I,J) = ',VAL
      GOTO 10
C
      STOP
      END
C
      SUBROUTINE DIFF(FILE1,FILE2,FILE3,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FOCK1(NDIM,NDIM),FOCK2(NDIM,NDIM),FOCK3(NDIM,NDIM)
      CHARACTER*16 FILE1,FILE2,FILE3
C
      OPEN(UNIT=11,FILE=TRIM(FILE1)//'.dat',STATUS='OLD')
      REWIND(UNIT=11)
      DO I=1,NDIM
        READ(11,*) (FOCK1(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=11)
C
      OPEN(UNIT=12,FILE=TRIM(FILE2)//'.dat',STATUS='OLD')
      REWIND(UNIT=12)
      DO I=1,NDIM
        READ(12,*) (FOCK2(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=12)
      
      DO I=1,NDIM
        DO J=1,NDIM
          FOCK3(I,J) = FOCK1(I,J)-FOCK2(I,J)
        ENDDO
      ENDDO
C
      OPEN(UNIT=13,FILE=TRIM(FILE3)//'.dat',STATUS='UNKNOWN')
      REWIND(UNIT=13)
      DO I=1,NDIM
        WRITE(13,*) (FOCK3(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=13)
C
      RETURN
      END    
C
C
      SUBROUTINE ELEMENT(FILE3,NDIM,I,J,VAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FOCK3(NDIM,NDIM)
      CHARACTER*16 FILE3
C
      OPEN(UNIT=13,FILE=FILE3,STATUS='OLD')
      REWIND(UNIT=13)
      DO K=1,NDIM
        READ(13,*) (FOCK3(K,L),L=1,NDIM)
      ENDDO
      CLOSE(UNIT=13)
C
      VAL = FOCK3(I,J)
C
      RETURN
      END      
C
C
      SUBROUTINE DIAG(FILE3,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FOCK3(NDIM,NDIM)
      CHARACTER*16 FILE3
C
      OPEN(UNIT=13,FILE=FILE3,STATUS='OLD')
      REWIND(UNIT=13)
      DO K=1,NDIM
        READ(13,*) (FOCK3(K,L),L=1,NDIM)
      ENDDO
      CLOSE(UNIT=13)
C
      DO K=1,NDIM
        WRITE(*,*) K,FOCK3(K,K)
      ENDDO
C
      RETURN
      END   
C
C
      SUBROUTINE GNUMTRX(TITLE,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   GGGGGG  NN    NN UU    UU MM       MM TTTTTTTT RRRRRRR  XX     XX  C
C  GG    GG NNN   NN UU    UU MMM     MMM    TT    RR    RR  XX   XX   C
C  GG       NNNN  NN UU    UU MMMM   MMMM    TT    RR    RR   XX XX    C
C  GG       NN NN NN UU    UU MM MM MM MM    TT    RR    RR    XXX     C
C  GG   GGG NN  NNNN UU    UU MM  MMM  MM    TT    RRRRRRR    XX XX    C
C  GG    GG NN   NNN UU    UU MM   M   MM    TT    RR    RR  XX   XX   C
C   GGGGGG  NN    NN  UUUUUU  MM       MM    TT    RR    RR XX     XX  C
C                                                                      C
C -------------------------------------------------------------------- C
C  GNUMTRX EXPORTS AN ARRAY TO AN EXTERNAL DATA FILE AND PLOTS IT.     C
C**********************************************************************C
C
      CHARACTER*16 TITLE
C
      DIMENSION ARRAY(NDIM,NDIM)
C
C     ACCESS EXTERNAL DATA FILE
      OPEN(UNIT=8,FILE=TRIM(TITLE)//".dat",STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO I=1,NDIM
        READ(8, *) (ARRAY(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=8)
C
      XEND = DFLOAT(NDIM)-0.5D0
      YEND = DFLOAT(NDIM)-0.5D0
C
C     WRITE GNUPLOT MAKE FILE
      OPEN(UNIT=9,FILE=TRIM(TITLE)//'.gnuplot',
     &                                                 STATUS='REPLACE')
      WRITE(9,'(A)') '#'//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#  Usage:'
      WRITE(9,'(A)') '#  gnuplot < '//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Terminal output specs'
      WRITE(9,'(A)') 'set terminal pdf size 4,4'
      WRITE(9,'(A)') 'set output "'//TRIM(TITLE)//'.pdf"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') 'load "pals/jet2.pal"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A,F6.1,A)') 'set xrange [-0.5:',XEND,']'
      WRITE(9,'(A,F6.1,A)') 'set yrange [',YEND,':-0.5] reverse'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot data to file'
      WRITE(9,'(A,I2,A,I2,A)') 'plot "'//TRIM(TITLE)//'.dat"'
     &                                    //' matrix with image notitle'
      CLOSE(UNIT=9)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot '//TRIM(TITLE)//'.gnuplot')
      CALL SYSTEM('xdg-open '//TRIM(TITLE)//'.pdf')
C
      RETURN
      END
           
