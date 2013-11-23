!INCLUDE 'MATHIO.F90'
!INCLUDE 'TENSOR.F90'
! ############### CONST ###################
MODULE CONST
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
END MODULE CONST
! ############### MODEL ###################
MODULE MODEL
	USE CONST
	INTEGER :: LEN = 6
	INTEGER :: MAX_CUT = 12  ! 16
	REAL    :: MAX_ERR = 0.
	INTEGER :: SWEEPS = 3
END MODULE MODEL
! ################ MATH ################### 
MODULE MATH
	USE TENSORIAL
	! definition of huge tensor
	TYPE HUGE_TENSOR
		TYPE(TENSOR) :: TEN ! tensor content
		REAL :: LEV ! scale by EXP(LEV)
	END TYPE HUGE_TENSOR
CONTAINS
! -------- huge tensor -----------
! rescale huge tensor
SUBROUTINE RESCALE(A)
! inout: A - huge tensor to be rescaled (in place)
	TYPE(HUGE_TENSOR), INTENT(INOUT) :: A
	! local variables
	REAL :: R ! rescale ratio
	
	R = SCALE(A%TEN%VALS) ! get the scale
	IF (R > 0.) THEN
		A%TEN%VALS = A%TEN%VALS/R ! rescale the vals
		A%LEV = A%LEV + LOG(R) ! update the lev
	END IF
END SUBROUTINE RESCALE
! -------- combinatorial ---------
! cal binomial
FUNCTION BINOMIAL(L, K) RESULT(N)
! N = L!/(L-K)!K!
	INTEGER, INTENT(IN) :: L, K
	INTEGER :: N
	! local variable
	INTEGER :: I, M
	
	! check validity
	IF (L < K) THEN ! if L < K, ill defined
		N = 0 ! set N = 0 and return
		RETURN
	END IF
	! now L >= K
	N = 1 ! initialize
	M = MIN(K, L-K) ! see which is smaller 
	DO I = 1, M ! run over M
		N = N*(L-I+1)/I ! cal binomial product
	END DO
END FUNCTION BINOMIAL
! end of module MATH
END MODULE MATH
! ############## DATAPOOL #################
MODULE DATAPOOL
! tensor conventions
! MPO (by default B-type):
!   A-type       B-type
!     ╭            ╮
!     3            3
! ─ 2 A 1 ─    ─ 1 B 2 ─
!     4            4
!     ╰            ╯
! MPS (by default A-type):
!   A-type       B-type
!     │            │
!     2            2
! ─ 1 A 3 ─    ─ 3 B 1 ─
	USE MODEL
	USE MATH
	! indices
	INTEGER :: L
	! tensors
	TYPE(TENSOR), POINTER :: TA, TB ! site-MPO
	TYPE(TENSOR), POINTER :: WA(:), WB(:) ! MPS
	TYPE(HUGE_TENSOR), POINTER :: DA(:), DB(:) ! block-MPO
	TYPE(TENSOR) :: P ! Schmidt tensor
	! variables
	REAL :: F       ! free energy
	REAL, POINTER :: FS(:)   ! free energy lattice dependence
	REAL, POINTER :: SS(:,:) ! entanglement entropy
	! private data storages (only assessable through pointers)
	TYPE(TENSOR), TARGET, PRIVATE :: T(2)
	TYPE(TENSOR), TARGET, ALLOCATABLE, PRIVATE :: W(:)
	TYPE(HUGE_TENSOR), TARGET, ALLOCATABLE, PRIVATE :: D(:)
	REAL, TARGET, ALLOCATABLE, PRIVATE :: RWORK(:)
CONTAINS
! allocate space for datapool
SUBROUTINE DAT_ALLOCATE()
	! allocation control
	INTEGER, SAVE :: LEN0 = 0
	
	! check validity of system size LEN
	IF (MODULO(LEN,2) == 1 .OR. LEN < 4) THEN
		WRITE (*,'(A)') 'DAT_ALLOCATE::xlen: LEN must be even and greater than 4.'
		STOP
	END IF
	! now LEN is valid, compare with LEN0
	IF (LEN /= LEN0) THEN ! if size changed
		IF (ALLOCATED(W)) DEALLOCATE (W) ! try deallocate
		ALLOCATE(W(2*LEN)) ! reallocate new size
		IF (ALLOCATED(D)) DEALLOCATE (D) ! try deallocate
		ALLOCATE(D(2*LEN)) ! reallocate new size
		IF (ALLOCATED(RWORK)) DEALLOCATE (RWORK) ! try deallocate
		ALLOCATE(RWORK(2*LEN+LEN/2)) ! reallocate new size
		LEN0 = LEN ! keep the current size
	END IF
	! make initial association
	CALL DAT_ASSOCIATE('I')
END SUBROUTINE DAT_ALLOCATE
! set pointers to data blocks
SUBROUTINE DAT_ASSOCIATE(MODE)
! input: MODE - DMRG mode: 'I' - initial, 'F' - forward, 'B' - backward
	CHARACTER, INTENT(IN) :: MODE
	
	SELECT CASE(MODE)
		CASE ('I','F') ! initial or forward
			TA => T(1)
			TB => T(2)
			WA => W(:LEN)
			WB => W(LEN+1:)
			DA => D(:LEN)
			DB => D(LEN+1:)
		CASE ('B') ! backward
			TA => T(2)
			TB => T(1)
			WA => W(LEN+1:)
			WB => W(:LEN)
			DA => D(LEN+1:)			
			DB => D(:LEN)
	END SELECT
	SS(1:LEN,1:2) => RWORK(1:2*LEN)
	FS => RWORK(2*LEN+1:2*LEN+LEN/2)
END SUBROUTINE DAT_ASSOCIATE
! put MPO to DATAPOOL
SUBROUTINE SET_MPO(MPO)
! * DATAPOOL must be initiated before calling me
	TYPE(TENSOR), INTENT(IN) :: MPO
	
	IF (ASSOCIATED(TB) .AND. ASSOCIATED(TA)) THEN
		TB = MPO ! put MPO to TB
		! set TA from TB MPO
		TA = TEN_TRANS(TB,[2,1,3,4]) ! given by 1 <-> 2
	ELSE
		WRITE(*,'(A)')'SET_MPO::nptr: pointer TA or TB not associated.'
	END IF
END SUBROUTINE SET_MPO
! end of module DATAPOOL
END MODULE DATAPOOL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
CONTAINS
! ----------- MPO tensors -----------
! make MPO tensor B given the quantum number M
FUNCTION MAKE_B(M) RESULT(B)
! input: M - physical quantum number
! output: B - MPO tensor
	USE CONST
	USE MATH
	INTEGER, INTENT(IN) :: M
	TYPE(TENSOR) :: B
	! local variables
	COMPLEX :: A(0:1,0:1,0:1,0:1), Z
	INTEGER :: I1, I2, I3, I4
	
	A = Z0 ! clear
	Z = Z1/SQRT(REAL(BINOMIAL(4,M+2))) ! save pre-factor
	FORALL (I1=0:1,I2=0:1,I3=0:1,I4=0:1,-I1-I2+I3+I4==M)
		A(I4,I2,I1,I3) = Z*(2*I1-1)*(2*I2-1) ! assign
	END FORALL
	! make B tensor
	B%DIMS = [2,2,2,2] ! set dimension
	CALL TEN_ABSTRACT(A,B)
END FUNCTION MAKE_B
! ----------- DMRG kernel -----------
! DMRG controlling routine
SUBROUTINE DMRG(MPO)
! data transfered by DATAPOOL
! DATAPOOL: L. MODEL: LEN, SWEEPS
	USE DATAPOOL
	TYPE(TENSOR), INTENT(IN) :: MPO 
	! local variables
	INTEGER :: ITER
	
	CALL DAT_ALLOCATE() ! allocate data pool
	CALL SET_MPO(MPO) ! set MPO tensors
	! iDMRG (warm up)
	DO L = 1, LEN/2
		CALL DMRG_STEP('I')
	END DO
	! fDMRG (first sweep)
	DO L = LEN/2+1, LEN
		CALL DMRG_STEP('F')
	END DO
	! fDMRG sweeps
	DO ITER = 1, SWEEPS
		! fFDMRG (backward sweep)
		DO L = 1, LEN
			CALL DMRG_STEP('B')
		END DO
		! fDMRG (forward sweep)
		DO L = 1, LEN
			CALL DMRG_STEP('F')
		END DO
	END DO
END SUBROUTINE DMRG
! make one DMRG step
SUBROUTINE DMRG_STEP(MODE)
! input: MODE - DMRG mode: 'I' - initial, 'F' - forward, 'B' - backward
	USE DATAPOOL ! for DAT_ASSOCIATE
	CHARACTER, INTENT(IN) :: MODE
	
	CALL DAT_ASSOCIATE(MODE)
	SELECT CASE(MODE)
		CASE ('I')
			CALL IDMRG()
		CASE ('F','B')
			CALL FDMRG()
	END SELECT
	CALL ISSUE(MODE)
END SUBROUTINE DMRG_STEP
! infinite-size DMRG step
SUBROUTINE IDMRG()
! perform one step of iDMRG
! data transfered by DATAPOOL
! DATAPOOL: TA, TB, WA, WB, DA, DB, P. MODEL: MAX_CUT, MAX_ERR
! on exit DA, DB, WA, WB, P0, P will be updated
	USE DATAPOOL
	! local tensors
	TYPE(TENSOR) :: W, TS
	TYPE(TENSOR), SAVE :: P0
	! local variables
	INTEGER :: DPHY
	COMPLEX :: TVAL
	
	IF (L == 1) THEN
		! initialize tensors
		DPHY = TB%DIMS(3) ! get physical dim
		! set initial Schmidt spectrum
		P0 = EYE_TEN([1,1])
		P = EYE_TEN([DPHY,DPHY],SQRT(Z1/DPHY))
		! set initial MPS (boundary)
		!   2       2
		! 1 WA 3  3 WB 1
		! initial MPS transfers d.o.f. from leg 2 to 3, with leg 1 dummy
		WA(1) = TEN_PROD(TENSOR([1],[0],[Z1]),EYE_TEN([DPHY,DPHY]))
		WB(1) = WA(1) ! WB is the same as WA
		F = -LOG(ABS(EVALUATE(TB,TA))) ! set an initial F
		! set initial block MPO (and rescale)
		DA(1)%TEN = TA
		DA(1)%LEV = 0.
		DB(1)%TEN = TB
		DB(1)%LEV = 0.
		CALL RESCALE(DA(1))
		CALL RESCALE(DB(1))
	ELSE ! L >= 2
		! estimate trial W
		P0%VALS = Z1/P0%VALS ! cal P0^(-1) -> P0
		W = MAKE_W5(WA(L-1), WB(L-1), P, P0)
		P0 = P ! P0 is used, update to P
		TS = MAKE_TS(TA, TB, DA(L-1)%TEN, DB(L-1)%TEN) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = -(DA(L-1)%LEV+DB(L-1)%LEV+LOG(ABS(TVAL))) ! calculate F
		! SVD split W, and update WA, WB, P
		CALL SVD(W,[1,2],[3,4],WA(L),WB(L),P,MAX_CUT,MAX_ERR)
		! update DA, DB and rescale
		DA(L)%TEN = NEW_DX(TA, WA(L), DA(L-1)%TEN)
		DA(L)%LEV = DA(L-1)%LEV
		DB(L)%TEN = NEW_DX(TB, WB(L), DB(L-1)%TEN)
		DB(L)%LEV = DB(L-1)%LEV
		CALL RESCALE(DA(L))
		CALL RESCALE(DB(L))
	END IF
END SUBROUTINE IDMRG
! finite-size DMRG step (forward convention)
SUBROUTINE FDMRG()
! perform one step of fDMRG (forward convention)
! data transfered by DATAPOOL
! DATAPOOL: TA, TB, WA, WB, DA, P. MODEL: LEN, MAX_CUT, MAX_ERR
! on exit DA, WA, P will be updated, DB not touched, WB destroyed
	USE DATAPOOL
	! local tensors
	TYPE(TENSOR) :: W, TS
	COMPLEX :: TVAL
	REAL, SAVE :: F0

	IF (L == 1) THEN
		! WA(1), no update, obtained from the last sweep
		! update DA, to restart
		! SSB treatment to be implemented here
		DA(1)%TEN = NEW_DX(TA, WA(1))
		DA(1)%LEV = 0.
		CALL RESCALE(DA(1)) ! rescale
		F = F0 ! retrieve F from last save
	ELSEIF (L == LEN-1) THEN
		W = MAKE_W3(P, WB(2), WB(1)) ! estimate trial W
		TS = MAKE_TS(TA, TB, DA(LEN-2)%TEN) ! construct TS (boundary)
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = -(DA(LEN-2)%LEV+LOG(ABS(TVAL))) ! calculate F
		F0 = F ! save F at the last site
		! SVD split W, update WA, P
		CALL SVD(W,[1,2],[3,4],WA(LEN-1),WB(1),P,MAX_CUT,MAX_ERR)
		! update DA and rescale
		DA(LEN-1)%TEN = NEW_DX(TA, WA(LEN-1), DA(LEN-2)%TEN)
		DA(LEN-1)%LEV = DA(LEN-2)%LEV
		CALL RESCALE(DA(LEN-1))
	ELSEIF (L == LEN) THEN
		! update the ending WA by P*WB
		WA(LEN) = TEN_TRANS(TEN_PROD(P,WB(1),[2],[3]),[1,3,2])
		! update DA and rescale
		DA(LEN)%TEN = NEW_DX(TA, WA(LEN), DA(LEN-1)%TEN)
		DA(LEN)%LEV = DA(LEN-1)%LEV
		CALL RESCALE(DA(LEN))
		! return F for the whole lattice
		F = -(DA(LEN)%LEV + LOG(ABS(EVALUATE(DA(LEN)%TEN))))
	ELSE ! 2 <= L <= LEN-2
		W = MAKE_W3(P, WB(LEN-L+1), WB(LEN-L)) ! estimate trial W
		TS = MAKE_TS(TA, TB, DA(L-1)%TEN, DB(LEN-L-1)%TEN) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = -(DA(L-1)%LEV+DB(LEN-L-1)%LEV+LOG(ABS(TVAL))) ! calculate F
		! SVD split W, update WA, P
		CALL SVD(W,[1,2],[3,4],WA(L),WB(LEN-L),P,MAX_CUT,MAX_ERR)
		! update DA and rescale
		DA(L)%TEN = NEW_DX(TA, WA(L), DA(L-1)%TEN)
		DA(L)%LEV = DA(L-1)%LEV
		CALL RESCALE(DA(L))
	END IF
END SUBROUTINE FDMRG
! --------- tensor makers ---------
! estimate iDMG trial state
FUNCTION MAKE_W5(WA, WB, P, PI) RESULT(W)
! call by IDMRG
	TYPE(TENSOR), INTENT(IN) :: WA, WB, P, PI
	TYPE(TENSOR) :: W
	
	! W = 
	!               2                 4
	!               │                 │
	!               2                 2
	! 1 ─ 1 P 2 ─ 3 WB 1 ─ 1 PI 2 ─ 1 WA 3 ─ 2 P 1 ─ 3
	W = TEN_PROD(TEN_PROD(TEN_PROD(P,WB,[2],[3]),PI,[2],[1]),TEN_PROD(P,WA,[2],[3]),[3],[2])
END FUNCTION MAKE_W5
! estimate fDMRG trial state
FUNCTION MAKE_W3(P, W1, W2) RESULT(W)
! called by IDMRG
! W1, W2 tensors can be either A or B type
	TYPE(TENSOR), INTENT(IN) :: P, W1, W2
	TYPE(TENSOR) :: W
	
	! W = 
	!               2        4
	!               │        │
	!               2        2
	! 1 ─ 1 P 2 ─ 3 W1 1 ─ 3 W2 1 ─ 3
	W = TEN_PROD(TEN_PROD(P,W1,[2],[3]),W2,[2],[3])
END FUNCTION MAKE_W3
! construct system tensor
FUNCTION MAKE_TS(TA, TB, DA, DB) RESULT (TS)
! input: TA,TB - site MPO, DA - A-block MPO, DB - B-block MPO
! output: TS - system MPO
! DB is optional: present => bulk, missing => boundary 
	TYPE(TENSOR), INTENT(IN) :: TA, TB, DA
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: DB
	TYPE(TENSOR) :: TS
	
	IF (PRESENT(DB)) THEN ! bulk algorithm
		! TS (bulk) = 
		!     1        3        7        5
		!     │        │        │        │
		!     3        3        3        3
		! ╭ 2 DA 1 ─ 2 TA 1 ─ 1 TB 2 ─ 1 DB 2 ╮
		! │   4        4        4        4    │
		! │   │        │        │        │    │
		! │   2        4        8        6    │
		! ╰───────────────────────────────────╯
		TS = TEN_PROD(TEN_PROD(DA,TA,[1],[2]),TEN_PROD(DB,TB,[1],[2]),[1,4],[1,4])
	ELSE ! boundary algorithm
		! TS (boundary) = 
		! 1        3        7             5
		! │        │        │             │
		! 3        3        3             1
		! DA 1 ─ 2 TA 1 ─ 1 TB 2 ─ 2 DA ⊗ I
		! 4        4        4             2
		! │        │        │             │
		! 2        4        8             6
		TS = TEN_PROD(TEN_PROD(TEN_PROD(DA,TA,[1],[2]),EYE_TEN([1,1])),TB,[1,4],[2,1])
	END IF
END FUNCTION MAKE_TS
! update DX given TX and WX
FUNCTION NEW_DX(TX, WX, DX) RESULT (DX1)
! input: TX - MPO, WX - MPS, DX - block MPO
! output: DX1 - new block MPO by packing TX and WX into the old
! DX optional: if present => bulk, missing => boundary
! zipper-order contraction algorithm is used
	TYPE(TENSOR), INTENT(IN) :: TX, WX
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: DX
	TYPE(TENSOR) :: DX1
	
	IF (PRESENT(DX)) THEN ! bulk algorithm
		! DX1 (bulk) = 
		!       ╭──── 1 WX* 3 ─ 3
		!       │        2
		!       │        │
		!       3        3
		! 2 ─ 2 DX 1 ─ 2 TX 1 ─ 1
		!       4        4
		!       │        │
		!       │        2 
		!       ╰───── 1 WX 3 ─ 4
		DX1 = TEN_PROD(TX,TEN_PROD(TEN_PROD(DX,TEN_CONJG(WX),[3],[1]),WX,[3],[1]),[2,3,4],[1,3,5])
	ELSE ! boundary algorithm
		! DX1 (boundary) = 
		!   ╭ 1 WX* 3 ─ 3
		!   │    2
		!   │    │
		!   ╯    3
		! 2 ── 2 TX 1 ─ 1
		!   ╮    4
		!   │    │
		!   │    2 
		!   ╰─ 1 WX 3 ─ 4
		DX1 = TEN_FLATTEN(TEN_PROD(TEN_PROD(TX,TEN_CONJG(WX),[3],[2]),WX,[3],[2]), [1,0,2,3,5,0,4,0,6])
	END IF
END FUNCTION NEW_DX
! ----------- Solvers ------------
! anneal the state W to the ground state of TS
FUNCTION ANNEAL(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
! return the corresponding eigen value TVAL
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	
! +++++++++ choose a solver here +++++++++++
!	TVAL = ANNEAL0(TS, W) ! self-made (fas, unstable)
	TVAL = ANNEAL1(TS, W) ! package (slower, robust)
! ++++++++++++++++++++++++++++++++++++++++++
END FUNCTION ANNEAL
! self-made Arnoldi algorithm (fast, unstable)
FUNCTION ANNEAL0(TS, W) RESULT (TVAL)
! called by ANNEAL
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! parameters
	INTEGER, PARAMETER :: N = 16 ! Krylov space dimension
	INTEGER, PARAMETER :: MAX_ITER = 500 ! max interation
	REAL, PARAMETER :: TOL = 1.E-12 ! allowed error of Tval
	! local variables
	INTEGER :: DIM, I, J, K, ITER, INFO
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:)
	COMPLEX, ALLOCATABLE :: V(:,:)
	COMPLEX :: TVAL0, A(N,N), D(N),Z(N,N),WORK(65*N)
	REAL :: RWORK(2*N)
	LOGICAL :: EXH
	
	! unpack data from tensor
	! collect leg-combined inds in TS and W (remember to +1)
	LINDS = COLLECT_INDS(TS,[1,3,5,7])+1
	RINDS = COLLECT_INDS(TS,[2,4,6,8])+1
	WINDS = COLLECT_INDS(W,[1,2,3,4])+1
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(V(DIM,N))
	V = Z0 ! clear
	! dump W to the 1st col of V
	FORALL (I = 1:SIZE(W%INDS))
		V(WINDS(I),1) = W%VALS(I)
	END FORALL
	V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
	! prepare to start Arnoldi iteration
	TVAL0 = Z0
	EXH = .FALSE. ! space exhausted flag
	! use Arnoldi iteration algorithm
	DO ITER = 1, MAX_ITER
		A = Z0 ! initialize Heisenberg matrix
		! construct Krylov space
		DO K = 2, N
			! apply TS to V(:,K-1) -> V(:,K)
			V(:,K) = Z0
			DO I = 1,SIZE(TS%INDS)
				V(LINDS(I),K) = V(LINDS(I),K) + TS%VALS(I)*V(RINDS(I),K-1)
			END DO
			! orthogonalization by stabilized Gram–Schmidt process
			DO J = 1, K-1
				A(J,K-1) = DOT_PRODUCT(V(:,J),V(:,K))
				V(:,K) = V(:,K) - A(J,K-1)*V(:,J)
			END DO
			! cal the norm of residual vector
			A(K,K-1) = SQRT(DOT_PRODUCT(V(:,K),V(:,K)))
			! if it is vanishing, the Arnoldi iteration has broken down
			IF (ABS(A(K,K-1))<TOL) THEN
				EXH = .TRUE.
				EXIT ! exit the iteration, stop construction of Krylov space
			END IF
			! otherwise, normalize the residual vect to a new basis vect
			V(:,K) = V(:,K)/A(K,K-1)
		END DO !K
		! now the Heisenberg matrix has been constructed
		! the action of TS is represented on the basis V as A
		! call LAPACK to diagonalize A
		CALL ZGEEV('N','V',N,A,N,D,Z,1,Z,N,WORK,65*N,RWORK,INFO)
		! now D holds the eigen vals, and Z the eigen vects
		! find the max abs eigen val
		I = LOC_MAX_MODE(D)
		TVAL = D(I) ! save it in Tval
		! reorganize the eigen vector from Z(:,I)
		V(:,1) = MATMUL(V,Z(:,I)) ! save to 1st col of V
		V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
		! check convergence
		! if space exhausted, or relative error < tol
		IF (EXH .OR. ABS((TVAL-TVAL0)/TVAL) < TOL) THEN
			! Arnoldi iteration has converge
			EXIT ! exit Arnoldi interation
		ELSE ! not converge, next iteration
			TVAL0 = TVAL ! save TVAL to TVAL0
		END IF
	END DO ! next Arnoldi interation
	! if exceed max iteration
	IF (ITER > MAX_ITER) THEN !then power iteration has not converge
		WRITE (*,'(A)') 'ANNEAL::fcov: Arnoldi iteration failed to converge.'
	END IF
	! reconstruct W tensor for output
	W%INDS = [(I,I=0,DIM-1)]
	W%VALS = V(:,1)
END FUNCTION ANNEAL0
! by calling to LAPACK and ARPACK (slower, robust)
FUNCTION ANNEAL1(TS, W) RESULT (TVAL)
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! parameters
    CHARACTER(*),PARAMETER :: BMAT = 'I'    ! regular eigen problem
    CHARACTER(*),PARAMETER :: WHICH = 'LM'  ! largest magnitude mode
    CHARACTER(*),PARAMETER :: HOWMNY = 'A'  ! get Ritz eigen vecs
    LOGICAL, PARAMETER     :: RVEC = .TRUE. ! calculate eigen vecs
    INTEGER, PARAMETER :: MAXITER = 300 ! max Arnoldi iterations
    REAL, PARAMETER    :: TOL = 1.D-12 ! tolerance
    INTEGER, PARAMETER :: NCV = 16 ! Krylov vector space dim
    INTEGER, PARAMETER :: NEV = 6  ! num of eigenvalues to be compute
    INTEGER, PARAMETER :: LWORKL = 3*NCV**2 + 5*NCV
	! local variables
	INTEGER :: NCALL ! num of calls of TS*W by ARPACK
	CHARACTER, PARAMETER :: BACK(20) = CHAR(8)
	! tensor related
	INTEGER :: I, P ! iterator/index
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:) ! inds to collect
	INTEGER :: N ! Hilbert space dim
	COMPLEX, ALLOCATABLE :: D(:)   ! eigen values
	COMPLEX, ALLOCATABLE :: Z(:,:) ! eigen vectors
	COMPLEX, ALLOCATABLE :: A(:,:) ! mat rep of TS
	! ARPACK related
	INTEGER :: IDO, INFO ! communication/information flag
	INTEGER :: NCONV  ! num of converged eigen vals
	INTEGER :: IPARAM(11), IPNTR(14) ! used by ARPACK
	COMPLEX :: SIGMA  ! shift, not referenced
	LOGICAL, ALLOCATABLE :: SEL(:) ! eigen vect selection
	REAL,    ALLOCATABLE :: RWORK(:) ! work space
	COMPLEX, ALLOCATABLE :: RESID(:) ! residual vect (to put initial vect)
	COMPLEX, ALLOCATABLE :: V(:,:)   ! Arnoldi basis vectors
	COMPLEX, ALLOCATABLE :: WORKL(:), WORKEV(:) ! work spaces
	COMPLEX, ALLOCATABLE, TARGET :: WORKD(:) ! work space
	COMPLEX, POINTER :: W0(:), W1(:) ! pointer to WORKD
	
	! cal total dim of W
	N = PRODUCT(W%DIMS) ! this will be the Hilbert space dim
	IF (N <= NCV) THEN ! for small scale problem
		! call LAPACK directly
		! unpack data from TS
		A = TEN2MAT(TS,[1,3,5,7],[2,4,6,8])
		IF (.NOT.(SIZE(A,1) == N .AND. SIZE(A,2) == N)) THEN
			WRITE (*,'(A)') 'ANNEAL::xdim: incompatible dimension between TS and W.'
			STOP
		END IF
		! allocate for work space
		ALLOCATE(D(N), Z(N,N), WORKD(65*N), RWORK(2*N))
		! call LAPACK to diagonalize A
		CALL ZGEEV('N','V', N, A, N, D, Z, 1, Z, N, WORKD,&
			65*N, RWORK, INFO)
		IF (INFO /= 0) THEN ! error handling
			IF (INFO > 0) THEN
				WRITE (*,'(A)') 'ANNEAL::fcnv: the QR algorithm failed.'
			ELSE
				WRITE (*,'(A,I3)') 'ANNEAL::xinf: ZGEEV error ', INFO
			END IF
		END IF
		! now D holds the eigen vals, and Z the eigen vects
		! find the max abs eigen val
		P = LOC_MAX_MODE(D)
	ELSE ! for large scale problem
		! call ARPACK
		! allocate for work space
		ALLOCATE(D(NCV), Z(N,NEV), RESID(N), V(N,NCV), WORKD(3*N),&
			WORKL(LWORKL), WORKEV(3*NCV), RWORK(NCV), SEL(NCV), STAT = INFO)
		! collect leg-combined inds in TS and W (remember to +1)
		LINDS = COLLECT_INDS(TS,[1,3,5,7])+1
		RINDS = COLLECT_INDS(TS,[2,4,6,8])+1
		WINDS = COLLECT_INDS(W,[1,2,3,4])+1
		! unpack initial vect from W tensor input
		RESID = Z0 ! clear
		FORALL (I = 1:SIZE(W%INDS))
			RESID(WINDS(I)) = W%VALS(I) ! dump
		END FORALL
		RESID = RESID/SQRT(DOT_PRODUCT(RESID,RESID)) ! normalize
		! set ARPACK parameters
		IPARAM(1) = 1       ! ishift
		IPARAM(3) = MAXITER ! maxitr
		IPARAM(7) = 1       ! model 1: A*x = lambda*x
		! now it is ready to start Arnoldi iteration
		IDO  = 0 ! initialize communication flag
		INFO = 1 ! RESID contains the initial residual vector
		! call ARPACK ----------------------
		! main loop (reverse communication)
		NCALL = 0
		WRITE (6,'(A,I5,A)',ADVANCE = 'NO') 'dim:',N,', iter:'
		DO WHILE (.TRUE.) ! repeatedly call ZNAUPD
			! until converge or exceed MAXITER
			CALL ZNAUPD(IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, &
				V, N, IPARAM, IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO)
			! take action according to communication flag IDO
			SELECT CASE (IDO)
				CASE (-1,1) ! compute W1 = TS*W0
					! monitor each call
					NCALL = NCALL +1 ! NCALL inc
					WRITE (6,'(I4)',ADVANCE = 'NO') NCALL
					! perform action of TS to W0 -> W1
					! point W0, W1 to the correct part of WORKD
					W0 => WORKD(IPNTR(1):IPNTR(1)+N-1)
					W1 => WORKD(IPNTR(2):IPNTR(2)+N-1)
					! carry out sparse matrix multiplication
					W1 = Z0 ! clear
					DO I = 1,SIZE(TS%INDS)
						W1(LINDS(I)) = W1(LINDS(I)) + TS%VALS(I)*W0(RINDS(I))
					END DO
					! monitor back space
					WRITE (6,'(4A1)',ADVANCE = 'NO') BACK(1:4)
				CASE (99) ! done
					EXIT ! exit Arnoldi iteration
			END SELECT
		END DO
		! clear the monitor
		WRITE (6,'(20A1)',ADVANCE = 'NO') BACK
		! Either we have convergence or there is an error    
		IF (INFO /= 0) THEN ! error handling
			SELECT CASE (INFO)
				CASE (1)
					WRITE (*,'(A)') 'ANNEAL::fcnv: Arnoldi iteration failed to converge.'
				CASE (3)
					WRITE (*,'(A)') 'ANNEAL::xsft: no shift can be applied, try larger NCV.'        	
				CASE DEFAULT
					WRITE (*,'(A,I3)') 'ANNEAL::xinf: ZNAUPD error ', INFO
					STOP
			END SELECT
		ELSE
			! No fatal errors occurred.
			! Post-Process using ZNEUPD.
			CALL ZNEUPD(RVEC, HOWMNY, SEL, D, Z, N, SIGMA, WORKEV, &
				BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, N, IPARAM, &
				IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO)
			IF (INFO /=0 ) THEN ! error handling
				WRITE (*,'(A,I3)') 'ANNEAL::xinf: ZNEUPD error ', INFO
				STOP
			END IF
			NCONV = IPARAM(5) ! get the num of convergent modes
		END IF
		! now D(:NCONV) holds the eigen vals
		! and Z(:,:NCONV) the eigen vects
		! find the max abs eigen val
		P = LOC_MAX_MODE(D(1:NCONV)) ! get position
	END IF
	! prepare output
	TVAL = D(P) ! output eigen val as TVAL
	! reconstruct W tensor from the corresponding mode
	W%INDS = [(I,I=0,N-1)] ! make index
	W%VALS = Z(:,P) ! dump data
END FUNCTION ANNEAL1
! locate max mode from eigen values
FUNCTION LOC_MAX_MODE(D) RESULT (P)
! input: D - eigen values of TS
! output: P - position where the max mode is located
	COMPLEX, INTENT(IN) :: D(:)
	INTEGER :: P
	! local variables
	INTEGER :: I
	COMPLEX :: LGDP, LGDI
	COMPLEX, ALLOCATABLE :: LGD(:)
	REAL, PARAMETER :: TOL = 1.E-7 ! allowed relative error
	
	! take log of D
	LGD = LOG(D)
	! assuming  P at the beginning
	P = 1
	LGDP = LGD(P) ! record the in-hand lgD
	! start searching for better choice
	DO I = 2, SIZE(D)
		LGDI = LGD(I) ! get the current lgD
		! compare with the one in hand
		IF (REAL(LGDI) > REAL(LGDP)+TOL) THEN
			! lgD(I) > lgD(P)
			P = I ! catch I
			LGDP = LGD(P) ! record in-hand
		ELSEIF (REAL(LGDI) > REAL(LGDP)-TOL) THEN
			! lgD(I) ~ lgD(P) (degenerated in magnitude)
			IF (ABS(IMAG(LGDI)) < ABS(IMAG(LGDP)) - TOL) THEN
				! |phase of I| < |phase of P|, D(I) is more real
				P = I ! catch I
				LGDP = LGD(P) ! record in-hand
			ELSEIF (ABS(IMAG(LGDI)) < ABS(IMAG(LGDP))+TOL) THEN
				! |phase of I| ~ |phase of P| (phase degenerated)
				IF (IMAG(LGDI) > IMAG(LGDP)) THEN
					! choose the positive one
					P = I ! catch I
					LGDP = LGD(P) ! record in-hand
				END IF
			END IF
		END IF
	END DO
END FUNCTION LOC_MAX_MODE
! ----------- Measure -----------
! measure operators O's on MPS X
FUNCTION MEASURE(OS) RESULT (M)
! input: OS - MPOs
! output: M - correlation function
! M is output in terms of a tensor:
! leg <-> operator, index on leg <-> position of operator 
	USE MODEL
	USE MATH ! BINOMIAL
	TYPE(TENSOR), INTENT(IN) :: OS(:)
	TYPE(TENSOR) :: M
	! local tensors
	TYPE(TENSOR) :: E
	! local variables
	INTEGER :: I, K, N
	
	K = SIZE(OS) ! get the num of operators
	! deal with some special cases
	IF (K == 0) THEN ! if no operator
		M = ZERO_TEN([INTEGER::]) ! return null tensor
		RETURN	
	END IF
	IF (LEN == 0) THEN ! if no lattice
		M = ZERO_TEN([(0,I=1,K)]) ! return null tensor
		RETURN
	END IF
	! now L and K are both non-zero
 	N = BINOMIAL(LEN, K) ! cal the size of correlation func
 	IF (N == 0) THEN ! if not able to put down operators
 		M = ZERO_TEN([(LEN,I=1,K)]) ! return null tensor
 		RETURN
 	END IF
 	! now N is also non-zero, the normal case
	! allocate space for M
	ALLOCATE(M%DIMS(K),M%INDS(N),M%VALS(N))
	M%DIMS = LEN ! each leg can go through the lattice LEN
	M%INDS = 0  ! clean up
	M%VALS = Z0 ! clean up
	! carry out recursive measurement
	CALL HSET_M(LEN, M%INDS, M%VALS, OS) ! enter without environment
	! on exit, M has been filled with measurement data, return
END FUNCTION MEASURE
! measurement kernel (simple)
RECURSIVE SUBROUTINE SSET_M(LA, INDS, VALS, OS, E0)
! input: LA - num of remaining sites
!        INDS - indices, VALS - measurement values
!        OS - MPOs, E0 - environment
! output: INDS, VALS - modified by new measurement data
! E0 optional: if missing, make it initialized
	USE MODEL ! LEN
	USE DATAPOOL ! WA(:)
	INTEGER, INTENT(IN)    :: LA
	INTEGER, INTENT(INOUT) :: INDS(:)
	COMPLEX, INTENT(INOUT) :: VALS(:)
	TYPE(TENSOR), INTENT(IN) :: OS(:)
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: E0
	! local tensor
	TYPE(TENSOR) :: E
	! local variables
	INTEGER :: K, I, N1, N2
		
	K = SIZE(OS) ! get the num of operators
	! check E0 input
	IF (PRESENT(E0)) THEN ! if given
		E = E0 ! use it as the environment
	ELSE ! if not given, make it from OS(K)
		SELECT CASE (SIZE(OS(K)%DIMS)) ! branch by # of legs
			CASE (2) ! 2-leg case
				E = TEN_PROD(EYE_TEN([1,1]),EYE_TEN([1,1]))
			CASE (4) ! 4-leg case
				E = TEN_PROD(EYE_TEN(OS(K)%DIMS([1,2])),EYE_TEN([1,1]))
			CASE DEFAULT
				WRITE (*,'(A)') 'SSET_M::ornk: rank of observable tensors must be 2 or 4.'
				STOP
		END SELECT
	END IF
	! now E has been set, lay down operators
	IF (K == 1) THEN ! for the last operator OS(1)
		! make direct measurement on the current lattice
		DO I = LA, 1, -1
			! record the index of measurement
			INDS(I) = I - 1
			! grow E with OS(1) at WA(I) and evaluate
			VALS(I) = EVALUATE(GROW(E, WA(I), OS(1))) ! rec the data
			! now E can be updated to I-1 lattice
			E = GROW(E, WA(I)) ! by absorbing WA(I)
		END DO
	ELSE ! K > 1 more than one operator remaining
		N1 = SIZE(VALS) ! get the size of correlation function
		DO I = LA, K, -1 ! only down to K, otherwise can't fit operators
			! cal the workspace range
			N2 = N1 ! upper range from previous lower range
			N1 = N1*(I - K)/I ! lower range update
			! grow E with the last OS(K) at WS(I)
			! passing the remaining OS to the lower level SET_M
			! with the workspace bounded by N1+1:N2 
			CALL SSET_M(I-1,INDS(N1+1:N2),VALS(N1+1:N2),OS(:K-1),GROW(E,WA(I),OS(K)))
			! on return, VALS contains the measurement values
			! and INDS contains the indices in the lower level
			! for this level, inds must be lifted by LEN and shifted by (I-1)
			INDS(N1+1:N2) = INDS(N1+1:N2) + (I-1)*LEN**(K-1)
			! the measurement has complete
			! now E can be updated to I-1 lattice
			E = GROW(E, WA(I)) ! by absorbing WA(I)
		END DO
	END IF
END SUBROUTINE SSET_M
! measurement kernel (huge)
RECURSIVE SUBROUTINE HSET_M(LA, INDS, VALS, OS, E0)
! input: LA - num of remaining sites
!        INDS - indices, VALS - measurement values
!        OS - MPOs, E0 - environment
! output: INDS, VALS - modified by new measurement data
! E0 optional: if missing, make it initialized
	USE MODEL ! LEN
	USE MATH  ! HUGE_TENSOR
	USE DATAPOOL ! WA(:), TB
	INTEGER, INTENT(IN)    :: LA
	INTEGER, INTENT(INOUT) :: INDS(:)
	COMPLEX, INTENT(INOUT) :: VALS(:)
	TYPE(TENSOR), INTENT(IN) :: OS(:)
	TYPE(HUGE_TENSOR), OPTIONAL, INTENT(IN) :: E0
	! local tensor
	TYPE(HUGE_TENSOR) :: E, E1
	! local variables
	INTEGER :: K, I, N1, N2
		
	K = SIZE(OS) ! get the num of operators
	! check E0 input
	IF (PRESENT(E0)) THEN ! if given
		E = E0 ! use it as the environment
	ELSE ! if not given, make it from OS(K)
		SELECT CASE (SIZE(OS(K)%DIMS)) ! branch by # of legs
			CASE (4) ! 4-leg case
				E%TEN = TEN_PROD(EYE_TEN(OS(K)%DIMS([1,2])),EYE_TEN([1,1]))
				E%LEV = 0.
				CALL RESCALE(E)
			CASE DEFAULT
				WRITE (*,'(A)') 'HSET_M::ornk: rank of observable tensors must be 4.'
				STOP
		END SELECT
	END IF
	! now E has been set, lay down operators
	IF (K == 1) THEN ! for the last operator OS(1)
		! make direct measurement on the current lattice
		DO I = LA, 1, -1
			! record the index of measurement
			INDS(I) = I - 1
			IF (I > 1) THEN
				! grow E with OS(1) at WA(I) and evaluate
				VALS(I) = EVALUATE(GROW(E%TEN, WA(I), OS(1)), DA(I-1)%TEN)
				VALS(I) = VALS(I)*EXP(E%LEV + DA(I-1)%LEV+FS(LEN/2))
				! now E can be updated to I-1 lattice
				E%TEN = GROW(E%TEN, WA(I), TB) ! by absorbing WA(I) and TB
				CALL RESCALE(E)
			ELSE ! for the last site I = 1
				! grow E with OS(1) at WA(1) and evaluate
				VALS(1) = EVALUATE(GROW(E%TEN, WA(1), OS(1)))
				VALS(1) = VALS(1)*EXP(E%LEV+FS(LEN/2))
			END IF
		END DO
	ELSE ! K > 1 more than one operator remaining
		N1 = SIZE(VALS) ! get the size of correlation function
		DO I = LA, K, -1 ! only down to K, otherwise can't fit operators
			! cal the workspace range
			N2 = N1 ! upper range from previous lower range
			N1 = N1*(I - K)/I ! lower range update
			! grow E with the last OS(K) at WS(I)
			E1%TEN = GROW(E%TEN,WA(I),OS(K))
			E1%LEV = E%LEV
			! passing the remaining OS to the lower level SET_M
			! with the workspace bounded by N1+1:N2 
			CALL HSET_M(I-1,INDS(N1+1:N2),VALS(N1+1:N2),OS(:K-1),E1)
			! on return, VALS contains the measurement values
			! and INDS contains the indices in the lower level
			! for this level, inds must be lifted by LEN and shifted by (I-1)
			INDS(N1+1:N2) = INDS(N1+1:N2) + (I-1)*LEN**(K-1)
			! the measurement has complete
			! now E can be updated to I-1 lattice
			E%TEN = GROW(E%TEN, WA(I), TB) ! by absorbing WA(I) and TB
			CALL RESCALE(E)
		END DO
	END IF
END SUBROUTINE HSET_M
! grow environment tensor
FUNCTION GROW(E, W, O) RESULT(E1)
! input: E - environment tensor, W - MPS tensor
! output: E1 - new environment tensor
! if O is absent, treat O = 1
	TYPE(TENSOR), INTENT(IN) :: E, W
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: O
	TYPE(TENSOR) :: E1
	
	! to prevent accessing undefined W, check here
	IF (.NOT. ALLOCATED(W%DIMS)) THEN ! if tensor not defined
		WRITE (*,'(A)') 'GROW::xdef: MPS tensor W not defined.'
		STOP
	END IF
	IF (PRESENT(O)) THEN ! O is given
		SELECT CASE (SIZE(O%DIMS)) ! branch by # of legs of O
			CASE (2) ! 2-leg case
				! E1 (2-leg O) = 
				! 3 ─ 1 W* 3 ──────╮
				!       2          │
				!       │          │
				!       1          3
				!       O    1 ─ 1 E 2 ─ 2
				!       2          4
				!       │          │
				!       2          │
				! 4 ─ 1 W 3 ───────╯
				E1 = TEN_PROD(O,TEN_PROD(TEN_PROD(E,TEN_CONJG(W),[3],[3]),W,[3],[3]),[1,2],[4,6])			
			CASE (4) ! 4-leg case
				! E1 (4-leg O) = 
				! 3 ─ 1 W* 3 ────╮
				!       2        │
				!       │        │
				!       3        3
				! 1 ─ 1 O 2 ── 1 E 2 ─ 2
				!       4        4
				!       │        │
				!       2        │
				! 4 ─ 1 W 3 ─────╯
				E1 = TEN_PROD(O,TEN_PROD(TEN_PROD(E,TEN_CONJG(W),[3],[3]),W,[3],[3]),[2,3,4],[1,4,6])			
		END SELECT
	ELSE ! O is absent
		! E1 (missing O) = 
		! 3 ─ 1 W* 3 ──────╮
		!       2          3
		!       │    1 ─ 1 E 2 ─ 2
		!       2          4
		! 4 ─ 1 W 3 ───────╯
		E1 = TEN_PROD(TEN_PROD(E,TEN_CONJG(W),[3],[3]),W,[3,5],[3,2])
	END IF
END FUNCTION GROW
! evaluate environment tensor
FUNCTION EVALUATE(E, D) RESULT(Z)
! input: E - environment tensor, D - block MPO, optional
! output: Z - measurement value of E, by tensor trace
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: E
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: D
	COMPLEX :: Z
	
	IF  (PRESENT(D)) THEN
		! contract E with D
		!     ╭───────╮
		!     3       3
		! ╭ 2 D 1 ─ 1 E 2 ╮
		! │   4       4   │
		! │   ╰───────╯   │
		! ╰───────────────╯
		Z = ZVAL(TEN_PROD(E,D,[1,2,3,4],[1,2,3,4]))
	ELSE
		! trace out the legs of E
		!  ╭───╮
		!  │   3
		! ╭┼ 1 E 2 ╮
		! ││   4   │
		! │╰───╯   │
		! ╰────────╯
		Z = ZVAL(TEN_TRACE(E,[1,3],[2,4]))
	END IF
END FUNCTION EVALUATE
! --------- physical -----------
! issue physical data after DMRG step
SUBROUTINE ISSUE(MODE)
	USE MODEL
	USE DATAPOOL
	CHARACTER, INTENT(IN) :: MODE
	! local variables
	INTEGER :: I1, I2
	CHARACTER(2) :: JN
	
	! site labeling system
	SELECT CASE(MODE)
		CASE('I')
			I1 = L
			I2 = L+1
			JN = '<>'
		CASE('F')
			I1 = L
			I2 = L+1
			JN = '->'
			IF (I2 == LEN+1) I2 = 1
		CASE('B')
			I1 = LEN-L
			I2 = LEN-L+1
			JN = '<-'
			IF (I1 == 0) I1 = LEN
	END SELECT
	! free energy recording system
	IF (MODE == 'I') THEN
		FS(L) = F
	ELSE ! for fDMRG, only improve the last free energy
		FS(LEN/2) = (FS(LEN/2)*(LEN-1) + F)/LEN
	END IF
	! entanglement entropy recording system
	IF (L == 0 .OR. L == LEN) THEN
		SS(L,1) = 0.
		SS(L,2) = 0.
	ELSE
		SS(L,1) = ENTROPY(P, 1)
		SS(L,2) = ENTROPY(P, 2)
	END IF
	! broadcast system
	WRITE (*,'(I3,A,I3,A,F10.6,A,F6.3,A,F5.2,A)') I1,JN,I2,': F = ', F, ', S = ', SS(L,1)/LOG(2.), 'bit ', SVD_ERR*100,'%'
END SUBROUTINE ISSUE
! cal entanglement entropy
FUNCTION ENTROPY(P, N) RESULT (S)
! input: P - Schmidt spectrum, N - order of entropy
! update: S - entanglement entropy
	TYPE(TENSOR), INTENT(IN) :: P
	INTEGER, INTENT(IN) :: N
	REAL :: S
	! local variables
	REAL, ALLOCATABLE :: ES(:)
	
	ES = REALPART(P%VALS)**2 ! entanglement weight
	SELECT CASE(N)
		CASE (1) ! Shannon entropy
			S = -SUM(ES*LOG(ES))
		CASE (2:) ! Rényi entropy
			S = LOG(SUM(ES**N))/(1-N)
	END SELECT
END FUNCTION ENTROPY
! ----------- Debug ------------
! reconstruct wave function
FUNCTION MPS(WA) RESULT (W)
! input: WA - MPS tensors
! output: W - MPS state vector
	TYPE(TENSOR), INTENT(IN) :: WA(:)
	COMPLEX, ALLOCATABLE :: W(:)
	! local variables
	TYPE(TENSOR) :: WT
	INTEGER :: NW, IW, I
	
	NW = SIZE(WA) ! get size of WA
	IF (NW == 0) THEN
		W = [COMPLEX::]
	END IF
	WT = WA(1)
	DO IW = 2, NW
		WT = TEN_FLATTEN(TEN_PROD(WT, WA(IW),[3],[1]),[1,0,2,3,0,4])
	END DO
	WT = TEN_TRACE(WT,[1],[3])
	ALLOCATE(W(WT%DIMS(1)))
	FORALL (I = 1:SIZE(WT%INDS))
		W(WT%INDS(I)+1) = WT%VALS(I)
	END FORALL
END FUNCTION MPS
! end of module PHYSICS
END MODULE PHYSICS
! ################ TASK ####################
MODULE TASK
	USE PHYSICS
CONTAINS
! ------------ Data --------------
SUBROUTINE COLLECT()
	USE MODEL
	USE DATAPOOL
	TYPE(TENSOR) :: MPO, OS(2), DAT(3)
	INTEGER :: I
	CHARACTER(9) :: FILE
	
	WRITE (*,'(A,I3)') 'cut = ', MAX_CUT
	MPO = MAKE_B(0)
	CALL DMRG(MPO)
	OS(1) = MAKE_B(1)
	OS(2) = MAKE_B(-1)
	DAT(1) = MEASURE(OS)
	DAT(2)%DIMS = SHAPE(FS)
	CALL TEN_ABSTRACT(CMPLX(FS), DAT(2))
	DAT(3)%DIMS = SHAPE(SS)
	CALL TEN_ABSTRACT(CMPLX(SS), DAT(3))
	WRITE(FILE,'(A,I0.3,A,I0.2,A,I0.1)')'L',LEN,'D',MAX_CUT,'S',SWEEPS
	CALL TEN_SAVE('data center (AKLT)/'//FILE,DAT)
END SUBROUTINE COLLECT
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	USE DATAPOOL
	
END SUBROUTINE TEST
! test MAKE_B
SUBROUTINE TEST_B()
	TYPE(TENSOR) :: B
	
	B = MAKE_B(0)
	CALL TEN_PRINT(B)
END SUBROUTINE TEST_B
! test DMRG
SUBROUTINE TEST_DMRG()
	USE MODEL
	TYPE(TENSOR) :: MPO
	REAL :: T0, T1
	
	WRITE (*,'(A,I3)') 'cut = ', MAX_CUT
	MPO = MAKE_B(0)
	CALL CPU_TIME(T0)
	CALL DMRG(MPO)
	CALL CPU_TIME(T1)
	PRINT '(A,F7.3,A)', 'timing:', T1-T0, 's'
END SUBROUTINE TEST_DMRG
! test MEASURE
SUBROUTINE TEST_MEASURE()
	USE MATHIO
	USE MODEL
	USE DATAPOOL
	TYPE(TENSOR) :: MPO, M, OS(2)
	INTEGER :: I
	
	WRITE (*,'(A,I3)') 'cut = ', MAX_CUT
	MPO = MAKE_B(0)
	CALL DMRG(MPO)
	OS(1) = MAKE_B(1)
	OS(2) = MAKE_B(-1)
	M = MEASURE(OS)
	CALL TEN_SAVE('M',M)
	CALL EXPORT('FS',FS)
	CALL EXPORT('SS',SS)
END SUBROUTINE TEST_MEASURE
! test transfer matrix
SUBROUTINE TEST_TRANSF()
	USE MATHIO
	INTEGER, PARAMETER :: N = 2, L = 100
	TYPE(TENSOR) :: T
	INTEGER :: I, J, D
	INTEGER, ALLOCATABLE :: LBS(:)
	COMPLEX, ALLOCATABLE :: M(:,:), S(:,:), A(:,:), B(:,:)
	COMPLEX :: Z(2,L), A0, F
	
	T = MAKE_B(0)
	DO I = 1, N
		!           3
		!       ╭───┴───╮
		!       3       3
		! 1 ─ 1 T 2 ─ 1 T 2 ─ 2
		!       4       4
		!       ╰───┬───╯
		!           4
		T = TEN_FLATTEN(TEN_PROD(T,T,[2],[1]),[1,0,4,0,2,5,0,3,6])
	END DO
	T = TEN_TRACE(T,[1],[2])
	M = TEN2MAT(T,[1],[2])
	D = SIZE(M,1)
	CALL EXPORT('M',M)
!	ALLOCATE(LBS(2**N))
!	LBS = 0
!	LBS(1) = 3
!	LBS(2**N) = 3 
!	S = TEN2MAT(PAULI_MAT(LBS),[1],[2])
!	A = M
!	A0 = SUM([(A(I,I),I=1,D)])
!	A = A/A0
!	F = LOG(A0)/LOG(2.)
!	DO J = 1, L
!		B = MATMUL(S,A)
!		Z(1,J) = SUM([(B(I,I),I=1,D)])
!		Z(2,J) = F
!		A = MATMUL(M,A)
!		A0 = SUM([(A(I,I),I=1,D)])
!		A = A/A0
!		F = F + LOG(A0)/LOG(2.)
!	END DO
!	CALL EXPORT('Z',Z)
END SUBROUTINE TEST_TRANSF
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	INTEGER :: I
	PRINT *, '------------ DMRG -------------'

!	CALL TEST()
!	CALL TEST_B()
!	CALL TEST_DMRG()
!	CALL TEST_MEASURE()
!	CALL TEST_TRANSF()
	CALL COLLECT()
END PROGRAM MAIN