      module sharedData
         ! mpi stuff
         integer:: ME, NPROC
         logical:: MASWRK
      end !*************************************************************


      program mp2CorrEng
      use omp_lib
      use sharedData
      implicit double precision(a-h,o-z)

      ! energy var
      double precision:: E2, E2_mpi

      ! var dec
      character(80) :: filename
      double precision,allocatable,dimension(:) :: EIG
      double precision,allocatable,dimension(:,:) :: eij,eab,B32
      integer,allocatable :: TRI(:,:)



      ! mpi init
      include 'mpif.h'
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, ME, ierror)
      MASWRK=ME.eq.0

#ifdef CPU
write(*,*) 'You are running the code with CPU OpenMP'
#else
write(*,*) 'Yar are running the code serially'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! read input for gpu kernel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! input file name
      call get_command_argument(1,filename)

      ! open input file
      open(unit=500, file=filename,status='old', & 
                        iostat=ierr,action='read')

      ! read parameters
      read(500,*) NAUXBASD,NCOR,NACT,NVIR,NBF

      ! read MO energy
      ALLOCATE(EIG(NBF))
      do ii=1,NBF
         read(500,*) EIG(ii)
      enddo        

      ! read B32
      ALLOCATE(B32(NAUXBASD*NVIR,NACT))
      do iact=1,NACT
         do ixvrt=1,NAUXBASD*NVIR
            read(500,*) B32(ixvrt,iact)
         enddo
      enddo

      ! read mp2 corr energy
      read(500,*) E2_ref
!!!!!!!!!!! finish reading input for gpu kernel !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  

      st=omp_get_wtime()

      ! some parameters
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR
      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*NORB



      ! eliminate extra ranks when NPROC > NACT
      IF(ME.GT.NACT-1) then
         write(*,*) "rank skipped wwww", ME
         GOTO 120
      endif


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


      ! trapezoidal decomposition occ-occ pairs
      CALL RIMP2_TRAPE_DEC(LddiActStart,LddiActEnd,NACT)


      ! buffer array for nested loop fusion
      ALLOCATE(TRI(NACT*(NACT+1)/2,2))
      LddiActStart_t=LddiActStart
      IF(MASWRK) LddiActStart_t=LddiActStart+1
      ISIZE=0
      DO JACT=LddiActStart_t,LddiActEnd
         DO IACT=1,JACT-1
            ISIZE=ISIZE+1
            TRI(ISIZE,1)=JACT
            TRI(ISIZE,2)=IACT
         ENDDO
      ENDDO

      ! init E2_mpi
      E2_mpi=0.0D00

      ! corr energy accumulation
      CALL RIMP2_ENERGY_WHOLE                            &
          (E2_mpi,TRI,ISIZE,B32, eij,eab,                &
           LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR)


      et=omp_get_wtime()

      write(*,'(A20,F10.3)') "timing: ", et-st



  120 CONTINUE


      E2=0.0D00
      call MPI_REDUCE(E2_mpi,E2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) 


      if(MASWRK) THEN
         write(*,'(A25,I15)') "mpi processes:", NPROC
         write(*,'(A25,I15)') "omp threads:", omp_get_max_threads()
         write(*,'(A25,F15.10)') "kernel mp2 corr energy:", E2
         write(*,'(A25,F15.10)') "ref mp2 corr energy:", E2_ref
         ediff=E2-E2_ref

         if(abs(ediff).le.1.0D-6) then
            write(*,*) "passed"
         else
            write(*,*) "failed"
         endif
      endif


      call MPI_FINALIZE(ierr)

      END !*************************************************************



      SUBROUTINE RIMP2_TRAPE_DEC(Istart,Iend,N)
      use sharedData
      implicit double precision(a-h,o-z)

      IF(N .LE. NPROC) THEN
         Istart = ME + 1
         Iend = ME + 1
         RETURN
      ENDIF
   

      TOT = (N*(N+1)) / NPROC 
      Start = 1.0
      Istart=nint(Start)
      DO II = 0,ME
         TMP = TOT + Start*Start - Start
         End = (SQRT(4*TMP + 1.0) - 1.0)/2.0
         Iend=nint(End)
         IF(II.LT.ME) Start = End + 1
         Istart=nint(Start)
      ENDDO

      END !*************************************************************





      SUBROUTINE RIMP2_ENERGY_WHOLE                         &
                (E2,TRI,ISIZE, B32, eij,eab,                &
                 LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR)
      use omp_lib

      implicit double precision(a-h,o-z)

      ! output
      double precision :: E2

      ! input
      integer :: TRI(NACT*(NACT+1)/2,2)
      double precision :: eij(NACT,NACT)
      double precision :: eab(NVIR,NVIR)
      double precision :: B32(NAUXBASD*NVIR,NACT)

      ! local data
      double precision,save :: E2_omp
#ifdef CPU
!$omp threadprivate(E2_omp)
#endif
!       double precision,allocatable,save :: QVV(:,:)
! !$omp threadprivate(QVV)
      double precision,allocatable,dimension(:,:),save :: QVV, QVV_t
#ifdef CPU
!$omp threadprivate(QVV,QVV_t)
#endif



      ! turn off dynamics threads
      CALL OMP_SET_DYNAMIC(.FALSE.)

      ! env num threads
      Nthreads=omp_get_max_threads()

#ifdef CPU
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiActStart,LddiActEnd, Istart,Iend,                      &   
!$omp        NACT,NVIR,NAUXBASD,                                        &
!$omp        B32,eij,eab,E2,TRI,ISIZE)                                  &
!$omp private(IACT,JACT,w43,w44)
#endif

      E2_omp = 0.0D00

      ALLOCATE(QVV(NVIR,NVIR))
      ALLOCATE(QVV_t(NVIR,NVIR))
         
!!!!! IACT .ne. JACT
#ifdef CPU
!$omp do schedule(DYNAMIC)
#endif
      DO KK=1,ISIZE

         JACT=TRI(KK,1)
         IACT=TRI(KK,2)       

         CALL RIMP2_ENERGYIJ                                            &
             (E2_omp, B32(1,IACT),B32(1,JACT),                          &
              eij(IACT,JACT),eab, QVV,QVV_t,                                  &
              IACT,JACT,2.0D00,NVIR,NAUXBASD)

      ENDDO !KK
#ifdef CPU
!$omp end do nowait
#endif

!!!!! IACT .eq. JACT
#ifdef CPU
!$omp do schedule(DYNAMIC)
#endif
      DO JACT=LddiActStart,LddiActEnd           
         CALL RIMP2_ENERGYIJ                                            &
             (E2_omp, B32(1,JACT),B32(1,JACT),                          &
              eij(JACT,JACT),eab, QVV,QVV_t,                                  &
              JACT,JACT,1.0D00,NVIR,NAUXBASD)
      ENDDO !JACT
#ifdef CPU
!$omp end do
#endif


#ifdef CPU
!$omp atomic
#endif
      E2 = E2 + E2_omp

      DEALLOCATE(QVV,QVV_t)

#ifdef CPU
!$OMP END PARALLEL
#endif

      END !*************************************************************






      SUBROUTINE RIMP2_ENERGYIJ                                         &
                (E2, BI,BJ,eij,eab,                                     &
                 QVV,QVV_t, IACT,JACT,FAC,NVIR,NAUXBASD)

      implicit double precision(a-h,o-z)

      ! output
      double precision :: E2

      ! input
      double precision :: BI(NAUXBASD*NVIR)
      double precision :: BJ(NAUXBASD*NVIR)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij

      ! buffer
      double precision :: QVV(NVIR,NVIR)
      double precision :: QVV_t(NVIR,NVIR)

      ! double precision,allocatable,dimension(:,:) :: QVV_t


      ! ALLOCATE(QVV_t(NVIR,NVIR))


      CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BJ,NAUXBASD,                                  &
                  0.0D00, QVV_t,NVIR)

      QVV = transpose(QVV_t)
      QVV = QVV_t + QVV_t - QVV
      QVV_t = QVV_t/(eij - eab)
      E2_t = DDOT(NVIR*NVIR,QVV,1,QVV_t,1)
      E2 = E2 + FAC*E2_t
      
      ! DEALLOCATE(QVV_t)

    
      END
