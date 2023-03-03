!  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
!  Licensed under the NCSA open source license

#if defined(INTEL_OFFLOAD)
  include "mkl_omp_offload.f90"
#endif
      module rimp2_shared
      use omp_lib

#if defined(HIPBLAS)
      use iso_c_binding
      use hipfort
      use hipfort_check
      use hipfort_hipblas
      type(c_ptr) :: HIPBLAS_handle = c_null_ptr
#endif

#if defined(CUBLAS) || defined(CUBLASXT)
        use cublasf
#endif
      ! mpi stuff
      integer:: ME, NPROC
      logical:: MASWRK

#if defined(CUBLAS) || defined(CUBLASXT)
        integer(c_int)   :: cublas_return
        type(c_ptr)      :: cublas_handle
#if defined(CUBLASXT)
          integer(c_int),dimension(1) :: cublasXt_deviceId
#endif
#endif

      end module

      module rimp2_input
        double precision,allocatable,dimension(:,:) :: eij,eab
#if defined(HIPBLAS)        
        double precision,allocatable,target,dimension(:,:) :: B32
#else
        double precision,allocatable,dimension(:,:) :: B32
#endif

        double precision,allocatable,dimension(:) :: EIG
        integer:: NAUXBASD,NCOR,NACT,NVIR,NBF,NQVV
        double precision:: E2_ref
      end module




      program mp2CorrEng
      use rimp2_shared
      use rimp2_input
      implicit double precision(a-h,o-z)

      ! energy var
      double precision:: E2, E2_mpi

      ! timing
      double precision:: dt_mpi, dt_min, dt_max, dt_mean


#if !defined(NOMPI)
        ! mpi init
        include 'mpif.h'
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, ME, ierror)
#else
        NPROC=1
        ME=0
#endif
      MASWRK=ME.eq.0


      if(MASWRK) THEN
#ifdef CPU
        write(*,*) 'You are running the code with CPU OpenMP'
#elif defined(NVBLAS)
        write(*,*) 'You are running the code with nvblas on GPU'
#elif defined(INTEL_OFFLOAD)
        write(*,*) 'You are running the code with mkl on GPU'
#elif defined(CUBLAS)
        write(*,*) 'You are running the code with cublas on GPU'
#elif defined(CUBLASXT)
        write(*,*) 'You are running the code with cublasxt on GPU'
#elif defined(HIPBLAS)
        write(*,*) 'You are running the code with HIPBLAS on GPU' 
#else
        write(*,*) 'You are running the code serially'
#endif
      endif

#if defined(CUBLAS)
        cublas_return = cublascreate_v2(cublas_handle)
#elif defined(CUBLASXT)
        cublas_return = cublasXtcreate(cublas_handle)
        cublasXt_deviceId(1) = 0      ! only for one GPU
        cublas_return = cublasXtDeviceSelect(cublas_handle, 1, cublasXt_deviceId)  ! only for one GPU
        cublas_return = cublasXtSetBlockDim(cublas_handle, 2048)
#elif defined(HIPBLAS)
        call hipblasCheck(hipblasCreate(HIPBLAS_handle))
#endif


      ! Read or generate input data
      call Initialization(MASWRK)

!!!!!!!!!!! finish reading input for gpu kernel !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      ! eliminate extra ranks when NPROC > NACT
      IF(ME.GT.NACT-1) then
        write(*,*) "rank skipped wwww", ME
        GOTO 120
      endif

      ! trapezoidal decomposition occ-occ pairs
      CALL RIMP2_TRAPE_DEC(LddiActStart,LddiActEnd,NACT)

      ! Warming up
      E2_mpi=0.0D0
      st=omp_get_wtime()
      ! corr energy accumulation
#if defined(COMBINED_VERSION)
      CALL RIMP2_ENERGY_WHOLE_COMBINED(E2_mpi,         &
           LddiActStart,LddiActEnd)
#else
      CALL RIMP2_ENERGY_WHOLE(E2_mpi,B32,eij,eab,         &
           LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR,NQVV)
#endif
      ! Actual computation
      ! init E2_mpi
      E2_mpi=0.0D00
      print *, "time: ", omp_get_wtime()-st
      ! tic
      st=omp_get_wtime()

      ! corr energy accumulation
#if defined(COMBINED_VERSION)
      CALL RIMP2_ENERGY_WHOLE_COMBINED(E2_mpi,         &
           LddiActStart,LddiActEnd)
#else
      CALL RIMP2_ENERGY_WHOLE(E2_mpi,B32,eij,eab,         &
           LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR,NQVV)
#endif

      ! toc
      et=omp_get_wtime()
      dt_mpi = et - st

#if defined(CUBLAS)
        cublas_return = cublasdestroy_v2(cublas_handle)
#elif defined(CUBLASXT)
        cublas_return = cublasXtdestroy(cublas_handle)
#elif defined(HIPBLAS)
       call hipblasCheck(hipblasDestroy(HIPBLAS_handle))
#endif

  120 CONTINUE


#if !defined(NOMPI)
        E2=0.0D00
        call MPI_REDUCE(E2_mpi,E2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) 
        dt_mean = 0.0D0
        dt_min = 0.0D0
        dt_max = 0.0D0
        call MPI_REDUCE(dt_mpi,dt_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(dt_mpi,dt_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(dt_mpi,dt_mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        dt_mean = dt_mean/NPROC
#else
        E2=E2_mpi
        dt_min = dt_mpi
        dt_max = dt_mpi
        dt_mean = dt_mpi
#endif

      if(MASWRK) THEN
        ediff=E2-E2_ref
        Rel_E2_error=abs(ediff/E2_ref)
        write(*,'(/,5X,A)') "Results:"
        write(*,'(10X,A45,I5)') "Number of MPI ranks   = ", NPROC
        write(*,'(10X,A45,I5)') "Number of OMP threads = ", omp_get_max_threads()
        ! write(*,'(10X,A45,E12.5)') "Computed  MP2 corr. energy = ", E2
        ! write(*,'(10X,A45,E12.5)') "Reference MP2 corr. energy = ", E2_ref
        write(*,'(10X,A45,E12.5)') "Rel. error of computed MP2 corr. energy = ", Rel_E2_error
        write(*,'(10X,A45,F8.3,A)') "Wall time (minimum)   = ", dt_min," sec"
        write(*,'(10X,A45,F8.3,A)') "Wall time (mean)      = ", dt_mean," sec"
        write(*,'(10X,A45,F8.3,A)') "Wall time (maximum)   = ", dt_max," sec"
        if(Rel_E2_error.le.1.0D-6) then
           write(*,'(10X,A)') "Passed :-) "
        else
           write(*,'(10X,A)') "Failed :-("
        endif
      endif

      ! deallocate arrays
      call finalization


#if !defined(NOMPI)
        call MPI_FINALIZE(ierr)
#endif

      END !*************************************************************



      SUBROUTINE RIMP2_TRAPE_DEC(Istart,Iend,N)
      use rimp2_shared
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




      SUBROUTINE RIMP2_ENERGY_WHOLE_COMBINED(E2,             &
                 LddiActStart,LddiActEnd)

      use rimp2_shared
      use rimp2_input
#if defined(INTEL_OFFLOAD)
#if defined(MKL_ILP64)
use onemkl_blas_omp_offload_ilp64
#else
use onemkl_blas_omp_offload_lp64
#endif
#endif
      implicit double precision(a-h,o-z)

      ! output
      double precision :: E2
#if defined(HIPBLAS)
      integer :: m ,n,k
      integer :: lda, ldb, ldc
      double precision,parameter :: alpha=1.0D00, beta=0.0D0
      type(c_ptr) :: da = c_null_ptr, db = c_null_ptr, dc = c_null_ptr

      !local data
      double precision,allocatable,target,dimension(:,:,:),save :: QVV
#else
      !local data
      double precision,allocatable,dimension(:,:,:),save :: QVV
#endif
      ! local data
      double precision,save :: E2_omp

#ifdef CPU
        !$omp threadprivate(E2_omp,QVV)
#endif

      ! turn off dynamics threadss
      CALL OMP_SET_DYNAMIC(.FALSE.)

#ifdef CPU
        ! env num threads
        Nthreads=omp_get_max_threads()

        !$OMP PARALLEL NUM_THREADS(Nthreads)                              &
        !$omp default(none)                                               &
        !$omp shared(LddiActStart,LddiActEnd,                             &   
        !$omp        NACT,NVIR,NAUXBASD,B32,eij,eab,E2,NQVV)              &
        !$omp private(IACT,JACT,iQVV,IACTmod,num,num2,fac,ib_n,ic_n,ia_n2,ic_n2,tijab,q_t)
#endif

      E2_omp = 0.0D00
      ALLOCATE(QVV(NVIR,NQVV,NVIR))

#if defined(INTEL_OFFLOAD) || defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target enter data map(alloc: QVV) 
        !$omp target enter data map(to: eij,eab,B32) 
#endif

#ifdef CPU
        !$omp do schedule(DYNAMIC)
#endif
      DO JACT=LddiActStart,LddiActEnd
        DO IACTmod=1,(JACT-1)/NQVV+1
          IACT = (IACTmod-1)*NQVV+1
          IF(IACTmod*NQVV>JACT) THEN
            iQVV = JACT-(IACTmod-1)*NQVV
          ELSE
            iQVV = NQVV
          ENDIF



#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target data use_device_ptr(B32,QVV)
#endif
#if defined(CUBLAS) || defined(CUBLASXT)
#if defined(CUBLAS)
          cublas_return =  CUBLASDGEMM_v2 &
                   (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N,              &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                     1.0D00, B32(1,IACT),NAUXBASD,                                  &
                     B32(1,JACT),NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)
        cublas_return = cudaDeviceSynchronize()

#elif defined(CUBLASXT)
          cublas_return =  cublasXtDgemm  &
                   (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N,              &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                     1.0D00, B32(1,IACT),NAUXBASD,                                  &
                     B32(1,JACT),NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)
        cublas_return = cudaDeviceSynchronize()
#endif
#elif defined(HIPBLAS)
        m=NVIR*iQVV
        n=NVIR
        k=NAUXBASD
        da=c_loc(B32(1,IACT))
        lda=NAUXBASD
        db=c_loc(B32(1,JACT))
        ldb=NAUXBASD
        dc=c_loc(QVV(1,1,1))
        ldc=NVIR*iQVV

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle,HIPBLAS_OP_T,HIPBLAS_OP_N, &
                 m,n,k,alpha,da,lda,db,ldb,beta,dc,ldc))
        call hipCheck(hipDeviceSynchronize())

#elif defined(NVBLAS) || defined(CPU)
        CALL DGEMM &
           ('T','N',  &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                  1.0D00, B32(1,IACT),NAUXBASD,                                  &
                          B32(1,JACT),NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)
#elif defined(INTEL_OFFLOAD)
         !$omp target variant dispatch use_device_ptr(B32,QVV)

        CALL DGEMM &
           ('T','N',  &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                  1.0D00, B32(1,IACT),NAUXBASD,                                  &
                          B32(1,JACT),NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)


        !$omp end target variant dispatch

#else
        print *, "incorrect macro choices"
        stop 1
#endif

#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp end target data
#endif



#if defined(INTEL_OFFLOAD) || defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
#if defined(THREADLIMIT)
        !$omp target teams distribute parallel do collapse(3) reduction(+:E2_omp) private(ib,ia,ic,num,ib_n,ic_n,num2,ia_n2,ic_n2,tijab,q_t,fac) thread_limit(TLIMIT) num_teams(NTEAMS)
#else
        !$omp target teams distribute parallel do collapse(3) reduction(+:E2_omp)
#endif
#endif
        DO IC=1,iQVV
           DO IB=1,NVIR
              DO IA=1,NVIR
                 num = IC+(IB-1)*iQVV 
                 IB_n = (num -1)/ NQVV + 1
                 IC_n = mod(num-1, NQVV) + 1
                 num2 = IC+(IA-1)*iQVV 
                 IA_n2 = (num2 - 1) / NQVV + 1
                 IC_n2 = mod(num2-1, NQVV) + 1

                 Tijab=QVV(IA,IC_n,IB_n)/(eij(IACT+IC-1,JACT)-eab(IA,IB))
                 Q_t=QVV(IA,IC_n,IB_n)+QVV(IA,IC_n,IB_n)
                 FAC=2.0D00
                 IF(IACT+IC-1.EQ.JACT) FAC=1.0D00
                 E2_omp = E2_omp + FAC*Tijab*(Q_t-QVV(IB,IC_n2,IA_n2))

              ENDDO
           ENDDO
        ENDDO
#if defined(INTEL_OFFLOAD) || defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
!$omp end target teams distribute parallel do
#endif
        ENDDO
      ENDDO
#ifdef CPU
        !$omp end do
#endif

#if defined(INTEL_OFFLOAD) || defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target exit data map(release: QVV,eij,eab,B32) 
#endif

#ifdef CPU
        !$omp atomic
#endif
      E2 = E2 + E2_omp

      DEALLOCATE(QVV)


#ifdef CPU
        !$OMP END PARALLEL
#endif

      END !*************************************************************


      SUBROUTINE RIMP2_ENERGY_WHOLE(E2,B32,eij,eab,             &
                 LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR,NQVV)
      use rimp2_shared
      implicit double precision(a-h,o-z)

      ! output
      double precision :: E2

      ! input
      double precision :: eij(NACT,NACT)
      double precision :: eab(NVIR,NVIR)
      double precision :: B32(NAUXBASD*NVIR,NACT)

      ! local data
      double precision,save :: E2_omp
      double precision,allocatable,dimension(:,:,:),save :: QVV
#ifdef CPU
        !$omp threadprivate(E2_omp,QVV)
#endif

      ! turn off dynamics threadss
      CALL OMP_SET_DYNAMIC(.FALSE.)

#ifdef CPU
        ! env num threads
        Nthreads=omp_get_max_threads()

        !$OMP PARALLEL NUM_THREADS(Nthreads)                              &
        !$omp default(none)                                               &
        !$omp shared(LddiActStart,LddiActEnd,                             &   
        !$omp        NACT,NVIR,NAUXBASD,B32,eij,eab,E2,NQVV)              &
        !$omp private(IACT,JACT,iQVV,IACTmod)
#endif

      E2_omp = 0.0D00
      ALLOCATE(QVV(NVIR,NQVV,NVIR))

#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target enter data map(alloc: QVV) 
        !$omp target enter data map(to: eij,eab,B32) 
#endif

#ifdef CPU
        !$omp do schedule(DYNAMIC)
#endif
      DO JACT=LddiActStart,LddiActEnd
        DO IACTmod=1,(JACT-1)/NQVV+1
          IACT = (IACTmod-1)*NQVV+1
          IF(IACTmod*NQVV>JACT) THEN
            iQVV = JACT-(IACTmod-1)*NQVV
          ELSE
            iQVV = NQVV
          ENDIF

          CALL RIMP2_ENERGYIJ(E2_omp, B32(1,IACT),B32(1,JACT),           &
                 eij,eab, QVV,IACT,JACT,NACT, NVIR,NAUXBASD,iQVV)

        ENDDO
      ENDDO
#ifdef CPU
        !$omp end do
#endif

#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target exit data map(release: QVV,eij,eab,B32) 
#endif

#ifdef CPU
        !$omp atomic
#endif
      E2 = E2 + E2_omp

      DEALLOCATE(QVV)

#ifdef CPU
        !$OMP END PARALLEL
#endif

      END !*************************************************************




      SUBROUTINE RIMP2_ENERGYIJ(E2,BI,BJ,eij,eab,                     &
                 QVV,IACT,JACT,NACT,NVIR,NAUXBASD,iQVV)
      use rimp2_shared
      implicit double precision(a-h,o-z)

      ! output
      double precision :: E2

      ! input
      double precision :: BI(NAUXBASD*NVIR*iQVV)
      double precision :: BJ(NAUXBASD*NVIR)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij(NACT,NACT)

      ! buffer
      double precision :: QVV(NVIR,iQVV,NVIR)
#if defined(HIPBLAS)
      integer ::  m, n, k
      integer :: lda, ldb, ldc
      double precision,parameter ::  alpha=1.0D00, beta=0.0D0
      type(c_ptr) :: da = c_null_ptr, db = c_null_ptr, dc = c_null_ptr
#endif



#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target data use_device_ptr(BI,BJ,QVV)
#endif
#if defined(CUBLAS) || defined(CUBLASXT)
#if defined(CUBLAS)
          cublas_return =  CUBLASDGEMM_v2 &
                   (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N,              &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BJ,NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)
        cublas_return = cudaDeviceSynchronize()

#elif defined(CUBLASXT)
          cublas_return =  cublasXtDgemm  &
                   (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N,              &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BJ,NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)
        cublas_return = cudaDeviceSynchronize()
#endif
#elif defined(HIPBLAS)
        da=c_loc(BI(1))
        db=c_loc(BJ(1))
        dc=c_loc(QVV(1,1,1))
        m=NVIR*iQVV
        n=NVIR
        k=NAUXBASD
        lda=NAUXBASD
        ldb=NAUXBASD
        ldc=NVIR*iQVV

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
         HIPBLAS_OP_T, & ! integer(kind(hipblas_op_n)), value
         HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
         m, & ! integer(c_int), value
         n, & ! integer(c_int), value
         k, & ! integer(c_int), value
         alpha, & ! real(c_double)
         da, & ! type(c_ptr), value 
         lda, & ! integer(c_int), value
         db, & ! type(c_ptr), value
         ldb, & ! integer(c_int), value
         beta, & ! real(c_double)
         dc, & ! type(c_ptr), value
         ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#else
        CALL DGEMM &
           ('T','N',  &
                     NVIR*iQVV,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BJ,NAUXBASD,                                  &
                  0.0D00, QVV,NVIR*iQVV)
#endif

#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp end target data
#endif

#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp target map(tofrom:E2)
        !$omp teams distribute reduction(+:E2) 
#endif
      DO IC=1,iQVV
        E2_t = 0.0D00
#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
          !$omp parallel do reduction(+:E2_t) collapse(2)
#endif
        DO IB=1,NVIR
          DO IA=1,NVIR
            Tijab=QVV(IA,IC,IB)/(eij(IACT+IC-1,JACT)-eab(IA,IB))
            Q_t=QVV(IA,IC,IB)+QVV(IA,IC,IB)
            E2_t=E2_t + Tijab*(Q_t-QVV(IB,IC,IA))
          ENDDO
        ENDDO
#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
          !$omp end parallel do
#endif
         FAC=2.0D00
         IF(IACT+IC-1.EQ.JACT) FAC=1.0D00
         E2 = E2 + FAC*E2_t
      ENDDO
#if defined(NVBLAS) || defined(CUBLAS) || defined(CUBLASXT) || defined(HIPBLAS)
        !$omp end teams distribute
        !$omp end target
#endif

      END



      SUBROUTINE Initialization(MASWRK)
      use rimp2_input
      implicit double precision(a-h,o-z)

      logical, intent(in):: MASWRK

      ! var dec
      character(80) :: filename,tmp
      integer :: nwater,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Generate an arbitrary input that has the same structure of GAMESS data !!!

      ! input file name
      call get_command_argument(1,filename)

      if (filename.eq.'benz.kern'.or.filename.eq.'cor.kern'        &
      .or.filename.eq.'c60.kern' .or.filename.eq.'w30.kern'        &
      .or.filename.eq.'w60.kern' ) then
        if (MASWRK) then
          write(*,'(5x,A,A)') 'Reading data from ',filename
        endif
        call Read_Input_File(filename)
      else
        if (filename.eq.'benz.rand') then
          if (MASWRK) then
            write(*,'(5x,A)') 'Generating arbitrary input data with the structure of benz.kern'
          endif
          NAUXBASD=420
          NCOR=6
          NACT=15
          NVIR=93
          NBF=120
        else if (filename.eq.'cor.rand') then
          if (MASWRK) then
            write(*,'(5x,A)') 'Generating arbitrary input data with the structure of cor.kern'
          endif
          NAUXBASD=1512
          NCOR=24
          NACT=54
          NVIR=282
          NBF=384
        else if (filename.eq.'c60.rand') then
          if (MASWRK) then
            write(*,'(5x,A)') 'Generating arbitrary input data with the structure of c60.kern'
          endif
          NAUXBASD=3960
          NCOR=60
          NACT=120
          NVIR=360
          NBF=540
        else if (filename.eq.'w30.rand') then
          if (MASWRK) then
            write(*,'(5x,A)') 'Generating arbitrary input data with the structure of w30.kern'
          endif
          NAUXBASD=2520
          NCOR=30
          NACT=120
          NVIR=570
          NBF=750
        else if (filename.eq.'w60.rand') then
          if (MASWRK) then
            write(*,'(5x,A)') 'Generating arbitrary input data with the structure of w60.kern'
          endif
          NAUXBASD=5040
          NCOR=60
          NACT=240
          NVIR=1140
          NBF=1500
        else
          read(filename,*,iostat=ierr) nwater
          if ( (ierr.eq.0) .and. (nwater .gt. 0) ) then
            if (MASWRK) then
              write(*,'(5x,A,I3,A)') 'Generating arbitrary input data with the structure of',nwater,' water clusters'
            endif
            NAUXBASD=84*nwater
            NCOR=nwater
            NACT=4*nwater
            NVIR=19*nwater
            NBF=25*nwater
          else
            if (MASWRK) then
              write(*,'(5x,A,/,5x,A)') 'Error!','One of the followings should be used as an input:'
              write(*,'(10x,A)') 'benz.kern, cor.kern, c60.kern, w30.kern, or w60.kern for actual data sets, or'
              write(*,'(10x,A)') 'benz.rand, cor.rand, c60.rand, w30.rand, or w60.rand for arbitrary data sets.'
            endif
            stop
          endif
        endif
    
        ! Generate MO energy
        ALLOCATE(EIG(NBF))
        EIG=1.0d0
        do ii=1,NACT
           EIG(ii+NCOR) = 2.0d0
        enddo
  
        ! Generate B32
        ALLOCATE(B32(NAUXBASD*NVIR,NACT))
        B32=1.0d0/NAUXBASD
  
        ! Compute the corresponding mp2 corr energy
        E2_ref=0.5d0*(NVIR**2)*(NACT**2)/(NAUXBASD**2)

      endif

!!!!!!!!!!!!!!!!!!!!!!!! Finish the input generation !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      IF (command_argument_count().gt.1) then
        call get_command_argument(2,tmp)
        read(tmp,*) NQVV
      ELSE
        NQVV=NACT
      ENDIF

      if(MASWRK) THEN
        write(*,'(5x,A,5I6)') 'NAUXBASD,NCOR,NACT,NVIR,NBF = ',NAUXBASD,NCOR,NACT,NVIR,NBF
        write(*,'(5x,A,I10)') 'NQVV =',NQVV
        write(*,'(5x,A)') 'Memory Footprint:'
        write(*,'(10x,A4,I8,A,I7,A,F11.4,A)') 'B32(',NAUXBASD*NVIR,',',NACT,') = ',NAUXBASD*NVIR*NACT*8.D-6,' MB'
        write(*,'(10x,A4,I8,A,I7,A,F11.4,A)') 'eij(',NACT,',',NACT,') = ',NACT*NACT*8.D-6,' MB'
        write(*,'(10x,A4,I8,A,I7,A,F11.4,A)') 'eab(',NVIR,',',NVIR,') = ',NVIR*NVIR*8.D-6,' MB'
        write(*,'(10x,A4,I5,A,I4,A,I5,A,F11.4,A)') 'QVV(',NVIR,',',NACT,',',NVIR,') = ',NVIR*NACT*NVIR*8.D-6,' MB'
      endif

      END

      SUBROUTINE Finalization
        use rimp2_input
        implicit double precision(a-h,o-z)
        integer status
        deallocate( B32, stat=status)
        deallocate( EIG, stat=status)
        deallocate( eij, stat=status)
        deallocate( eab, stat=status)
      end SUBROUTINE Finalization

      SUBROUTINE Read_Input_File(filename)
      use rimp2_input
      implicit double precision(a-h,o-z)

      character(80), intent(in) :: filename
     
      ! open input file
      open(unit=500, file=filename,status='old',form="unformatted", &
#if defined(INTEL)
             access='direct',recl=10,iostat=ierr,action='read')
#else
             access='direct',recl=40,iostat=ierr,action='read')
#endif

      ! read parameters
      read(500,iostat=ierr,rec=1) NAUXBASD
      read(500,iostat=ierr,rec=2) NCOR
      read(500,iostat=ierr,rec=3) NACT
      read(500,iostat=ierr,rec=4) NVIR
      read(500,iostat=ierr,rec=5) NBF

      ! write MO energy
      ALLOCATE(EIG(NBF))
      do ii=1,NBF
        irec=ii+5
        read(500,iostat=ierr,rec=irec) EIG(ii)
      enddo

      ! read B32
      ALLOCATE(B32(NAUXBASD*NVIR,NACT))
      jrec=irec
      do iact=1,NACT
        do ixvrt=1,NAUXBASD*NVIR
          jrec=jrec+1
          read(500,iostat=ierr,rec=jrec) B32(ixvrt,iact)
        enddo
      enddo
      ! read mp2 corr energy
      read(500,iostat=ierr,rec=jrec+1) E2_ref

      END
