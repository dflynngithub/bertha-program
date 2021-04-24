C
C     CALCULATION TREE AND HAMILTONIAN
      LOGICAL     READIN
      CHARACTER*5 TREE
      CHARACTER*4 HMLT
      COMMON/OPTN/READIN,TREE,HMLT
C
C     MOLECULAR GEOMETRY
      CHARACTER*6 SHAPE
      COMMON/GEOM/SHAPE
C
C     SELF-CONSISTENT FIELD CALCULATION CHOICES
      LOGICAL     RACAH1,SSSSI3,SSSSI4,GAUNT2,VSDNSB,VSDNSC,
     &            DIISQC,DAMPFC,SCHWRZ,EQFILE,RCFILE,INTSYM,
     &            PRM1CT,BGAUNT,TRIANG,MQNSYM,CDPRM,PRM1KL
      COMMON/TGGL/RACAH1,SSSSI3,SSSSI4,GAUNT2,VSDNSB,VSDNSC,
     &            DIISQC,DAMPFC,SCHWRZ,EQFILE,RCFILE,INTSYM,
     &            PRM1CT,BGAUNT,TRIANG,MQNSYM,CDPRM,PRM1KL
C
C     PARALLEL COMPUTATION OPTIONS
      LOGICAL OPENMP
C     NPRCSR = OMP_GET_NUM_PROCS()
      COMMON/PRLL/NPRCSR,OPENMP