!  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
!  Licensed under the NCSA open source license

      module rimp2_shared
         use omp_lib
         double precision,allocatable,dimension(:,:) :: eij,eab,B32
         integer:: NAUXBASD,NACT,NVIR
      end module



      program mp2CorrEng
      use rimp2_shared
      implicit double precision(a-h,o-z)

      ! energy var
      double precision:: E2

      ! var dec
      character(80) :: filename,tmp
      double precision,allocatable,dimension(:) :: EIG

      ! timing
      double precision:: dt

      ! Generate input data
      call Initialization(E2_ref)

      ! Warming up
      E2=0.0D0
      CALL RIMP2_ENERGY_WHOLE(E2)

      ! Actual computation
      ! tic
      st=omp_get_wtime()

      ! corr energy accumulation
      E2=0.0D00
      CALL RIMP2_ENERGY_WHOLE(E2)

      ! toc
      et=omp_get_wtime()
      dt = et - st

      write(*,'(5X,A25,I5)') "Number of OMP threads:", omp_get_max_threads()
      write(*,'(5X,A25,E12.5)') "Computed mp2 corr energy:", E2
      write(*,'(5X,A25,E12.5)') "Ref. mp2 corr energy:", E2_ref
      write(*,'(5X,A25,F8.3,A)') "Wall time:", dt," sec"
      ediff=E2-E2_ref
      if(abs(ediff/E2_ref).le.1.0D-6) then
         write(*,'(10X,A)') "Passed :-) "
      else
         write(*,'(10X,A)') "Failed :-("
      endif

      END !*************************************************************



      SUBROUTINE RIMP2_ENERGY_WHOLE(E2)
      use rimp2_shared
      implicit double precision(a-h,o-z)

      ! output
      double precision :: E2

      ! local data
      double precision,save :: E2_omp
      double precision,allocatable,dimension(:,:,:),save :: QVV

      E2_omp = 0.0D00
      ALLOCATE(QVV(NVIR,NACT,NVIR))

! Hint: You may add target data map here
      IACT=1
      DO JACT=1,NACT
         CALL RIMP2_ENERGYIJ( E2_omp, B32(1,IACT),B32(1,JACT),      &
                 QVV,IACT,JACT )
      ENDDO
! Hint: You may add target data map here

      E2 = E2 + E2_omp
      DEALLOCATE(QVV)

      END !*************************************************************




      SUBROUTINE RIMP2_ENERGYIJ( E2,BI,BJ,QVV,IACT,JACT )
      use rimp2_shared
      implicit double precision(a-h,o-z)

      double precision :: E2
      double precision :: BI(NAUXBASD*NVIR*NACT),BJ(NAUXBASD*NVIR)
      double precision :: QVV(NVIR,JACT,NVIR)

! Hint: You may add target data map with use_device_ptr here
      CALL DGEMM &
           ('T','N',  &
                     NVIR*JACT,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BJ,NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*JACT)

! Hint: target region starts
      DO IC=1,JACT
         E2_t = 0.0D00
          DO IB=1,NVIR
            DO IA=1,NVIR
               Tijab=QVV(IA,IC,IB)/(eij(IACT+IC-1,JACT)-eab(IA,IB))
               Q_t=QVV(IA,IC,IB)+QVV(IA,IC,IB)
               E2_t=E2_t + Tijab*(Q_t-QVV(IB,IC,IA))
            ENDDO
         ENDDO

         FAC=2.0D00
         IF(IACT+IC-1.EQ.JACT) FAC=1.0D00
         E2 = E2 + FAC*E2_t
      ENDDO
! Hint: target region ends

      END




      SUBROUTINE Initialization(E2_ref)
      use rimp2_shared
      implicit double precision(a-h,o-z)

      double precision, intent(inout) :: E2_ref
     
      ! var dec
      character(80) :: filename
      double precision,allocatable,dimension(:) :: EIG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generate an arbitrary input that has the same structure of GAMESS data !!!

      ! input file name
      call get_command_argument(1,filename)

      if (filename.eq.'benz') then
         write(*,'(5x,A)') 'Using the structure of benz.kern'
         NAUXBASD=420
         NCOR=6
         NACT=15
         NVIR=93
         NBF=120
      else if (filename.eq.'cor'.or.filename.eq.'') then
         write(*,'(5x,A)') 'Using the structure of cor.kern'
         NAUXBASD=1512
         NCOR=24
         NACT=54
         NVIR=282
         NBF=384
      else if (filename.eq.'c60') then
         write(*,'(5x,A)') 'Using the structure of c60.kern'
         NAUXBASD=3960
         NCOR=60
         NACT=120
         NVIR=360
         NBF=540
      else if (filename.eq.'w30') then
         write(*,'(5x,A)') 'Using the structure of w30.kern'
         NAUXBASD=2520
         NCOR=30
         NACT=120
         NVIR=570
         NBF=750
      else if (filename.eq.'w60') then
         write(*,'(5x,A)') 'Using the structure of w60.kern'
         NAUXBASD=5040
         NCOR=60
         NACT=240
         NVIR=1140
         NBF=1500
      else
         write(*,'(5x,A,/,5x,A)') 'Error!','One of the following should be used as an input:'
         write(*,'(10x,A)') 'benz, cor, c60, w30, w60.'
         stop
      endif

      ! Write the variables used
      write(*,'(5x,A,5I5)') 'NAUXBASD,NCOR,NACT,NVIR,NBF = ',NAUXBASD,NCOR,NACT,NVIR,NBF

      ! Generate MO energy
      ALLOCATE(EIG(NBF))
      EIG=1.0d0
      do ii=1,NACT
         EIG(ii+NCOR) = 2.0d0
      enddo

      ! Generate B32
      ALLOCATE(B32(NAUXBASD*NVIR,NACT))
      B32=0.0d0
      do iact=1,NACT
         B32(1,iact) = 1.0d1/NACT
      enddo

      ! Compute the corresponding mp2 corr energy
      E2_ref=0.5d4/NACT/NACT

!!!!!!!!!!!!!!!!!!!!!!!! Finish the input generation !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,'(5x,A)') 'Memory Footprint:'
      write(*,'(10x,A10,I7,A,I5,A,F11.4,A)') 'B32(',NAUXBASD*NVIR,',',NACT,') = ',NAUXBASD*NVIR*NACT*8.D-6,' MB'
      write(*,'(10x,A10,I7,A,I5,A,F11.4,A)') 'eij(',NACT,',',NACT,') = ',NACT*NACT*8.D-6,' MB'
      write(*,'(10x,A10,I7,A,I5,A,F11.4,A)') 'eab(',NVIR,',',NVIR,') = ',NVIR*NVIR*8.D-6,' MB'
      write(*,'(10x,A10,I4,A,I3,A,I4,A,F11.4,A)') 'QVV(',NVIR,',',NACT,',',NVIR,') = ',NVIR*NACT*NVIR*8.D-6,' MB'

      ! some parameters
      NOCC=NCOR+NACT

      ! virt-virt MO energy pairs
      ALLOCATE(eab(NVIR,NVIR))
      DO IB=1,NVIR
         DO IA=1,IB
            eab(IA,IB) = EIG(IA+NOCC) + EIG(IB+NOCC)
            eab(IB,IA) = eab(IA,IB)
         ENDDO
      ENDDO

      ! occ-occ MO energy pairs
      ALLOCATE(eij(NACT,NACT))
      DO JJ=1,NACT
         DO II=1,JJ
            eij(II,JJ)=EIG(II+NCOR) + EIG(JJ+NCOR)
         ENDDO
      ENDDO


      END



