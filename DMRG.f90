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
!	REAL    :: BETA = 0.440687
	REAL    :: BETA = 0.
	REAL    :: THETA = 1.*PI
	INTEGER :: LEN = 20 ! must be even and larger than 4
	INTEGER :: MAX_CUT = 18
	REAL    :: MAX_ERR = 0.
END MODULE MODEL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
CONTAINS
! ------------ set MPO tensor ---------------
! square lattice MPO
SUBROUTINE SET_MPO(T)
	USE MODEL
	TYPE(TENSOR), INTENT(OUT) :: T ! MPO tensor output
	! local variables
	TYPE(TENSOR) :: X, Y, U, S
	COMPLEX, ALLOCATABLE :: A(:,:)
	COMPLEX :: Q, B
	
! ++++++++ set the vertex tensor here ++++++++
	Q = THETA * ZI
	B = BETA * Z1
	X =  TENSOR([2,2,2,2],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[EXP(2*B),EXP(Q/4.),EXP(Q/4.),EXP(Z0),EXP(Q/4.),EXP(-2*B - Q/2.),EXP(Z0),EXP(-Q/4.),EXP(Q/4.),EXP(Z0),EXP(-2*B + Q/2.),EXP(-Q/4.),EXP(Z0),EXP(-Q/4.),EXP(-Q/4.),EXP(2*B)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! symm SVD of X to unitary U and diagonal S
	CALL SYSVD(X,[1,4],[3,2],U,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	U = TEN_PROD(U,S,[3],[1]) ! attach the half S to U
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	U = TEN_PROD(Y,U,[1],[2])
	T = TEN_PROD(U,U,[1,3],[3,1])
!	T%VALS = T%VALS/2/EXP(ABS(BETA))
END SUBROUTINE SET_MPO
! ----------- DMRG -----------
! DMRG Kernel
SUBROUTINE DMRG(T, A, B)
! input: T - MPO site tensor
! output: A,B - MPS tensor, S - entanglement matrix
	USE MODEL
	TYPE(TENSOR), INTENT(IN)  :: T
	TYPE(TENSOR), INTENT(OUT) :: A(LEN), B(LEN)
	! local tensor
	TYPE(TENSOR) :: TA(LEN), TB(LEN), S
	! local variables
	INTEGER :: L, ITER
	COMPLEX :: TVAL
	! parameters
	INTEGER, PARAMETER :: SWEEPS = 1
	
	! check validity of system size LEN
	IF (MODULO(LEN,2) == 1 .OR. LEN <= 4) THEN
		WRITE (*,'(A)') 'DMRG::xlen: LEN must be even and greater than 4.'
		STOP
	END IF
	! iDMRG (warm up)
	DO L = 1, LEN/2
		CALL IDMRG(L, T, TA, TB, A, B, S, TVAL)
		CALL SHOW_DMRG_STATUS('I',L,TVAL,S%VALS)
	END DO
	! fDMRG (first sweep)
	DO L = LEN/2+1, LEN
		CALL FDMRG(L, T, TA, TB, A, B, S, TVAL)
		CALL SHOW_DMRG_STATUS('F',L,TVAL,S%VALS)
	END DO
	! fDMRG sweeps
	DO ITER = 1, SWEEPS
		! fFDMRG (backward sweep)
		DO L = 1, LEN
			CALL FDMRG(L, T, TB, TA, B, A, S, TVAL)
			CALL SHOW_DMRG_STATUS('B',L,TVAL,S%VALS)
		END DO
		! fFDMRG (forward sweep)
		DO L = 1, LEN
			CALL FDMRG(L, T, TA, TB, A, B, S, TVAL)
			CALL SHOW_DMRG_STATUS('F',L,TVAL,S%VALS)
		END DO
	END DO
END SUBROUTINE DMRG
! infinite-size DMRG step
SUBROUTINE IDMRG(L, T, TA, TB, A, B, S, TVAL)
! perform one step of iDMRG at L
! input: L - lattice position, T - MPO
! inout: TA, TB, A, B, S0, S will be updated
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: T
	TYPE(TENSOR), INTENT(INOUT) :: TA(:), TB(:), A(:), B(:), S
	COMPLEX, INTENT(INOUT) :: TVAL
	! local tensors
	TYPE(TENSOR) :: W, TS
	TYPE(TENSOR), SAVE :: S0
	! local variables
	INTEGER :: DPHY
	
	IF (L == 1) THEN
		! initialize tensors
		! boundary MPO
		TA(1) = T
		TB(1) = T
		! boundary MPS
		DPHY = T%DIMS(1) ! get physical dim
		A(1) = TEN_PROD(TENSOR([1],[0],[Z1]),EYE_TEN([DPHY,DPHY]))
		B(1) = A(1)
		! initial entanglement mat
		S0 = EYE_TEN([1,1])
		S = EYE_TEN([DPHY,DPHY],SQRT(Z1/DPHY))
		TVAL = Z1 ! return initial Tval
	ELSE ! L >= 2
		! estimate trial W
		S0%VALS = Z1/S0%VALS ! cal S0^(-1) -> S0
		W = MAKE_W5(A(L-1), B(L-1), S, S0)
		S0 = S ! S0 is used, update to S
		TS = GET_TS(T, TA(L-1), TB(L-1)) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		! SVD split W, and update A, B, S
		CALL SVD(W,[1,2],[3,4],A(L),B(L),S,MAX_CUT,MAX_ERR)
		! update TA, TB
		TA(L) = GET_TX(TA(L-1), T, A(L))
		TB(L) = GET_TX(TB(L-1), T, B(L))
	END IF
END SUBROUTINE IDMRG
! finite-size DMRG step
SUBROUTINE FDMRG(L, T, TA, TB, A, B, S, TVAL)
! perform one step of fDMRG at L (forward dirction)
! input: L - lattice position, T - MPO
! inout: TA, A, S will be updated
! on output: TB is not touched, but B is updated
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: T
	TYPE(TENSOR), INTENT(INOUT) :: TA(:), TB(:), A(:), B(:), S
	COMPLEX, INTENT(INOUT) :: TVAL
	! local tensors
	TYPE(TENSOR) :: W, TS
	COMPLEX, SAVE :: TVAL0
	
	IF (L == 1) THEN
		! A(1), no update, obtained from the last sweep
		! update TA, to restart
		! SSB treatment to be implemented here
		TA(1) = TEN_FLATTEN(TEN_PROD(TEN_PROD(TEN_CONJG(A(1)),T,[2],[1]),A(1), [4],[2]), [2,0,1,4,5,0,6,0,3])
		TVAL = TVAL0 ! retrieve Tval from last save
	ELSEIF (L == LEN-1) THEN
		! estimate trial W
		W = MAKE_W3(S, B(2), B(1))
		TS = GET_TS1(T, TA(LEN-2)) ! construct TS (boundary)
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		TVAL0 = TVAL ! save Tval at the last site
		! SVD split W, update A, S
		CALL SVD(W,[1,2],[3,4],A(LEN-1),B(1),S,MAX_CUT,MAX_ERR)
		! update TA
		TA(LEN-1) = GET_TX(TA(LEN-2), T, A(LEN-1))
	ELSEIF (L == LEN) THEN
		! update the ending A by S*B
		A(LEN) = TEN_TRANS(TEN_PROD(S,B(1),[2],[3]),[1,3,2])
		! update TA
		TA(LEN) = GET_TX(TA(LEN-1), T, A(LEN))
		! cal Tval by Tr (for the whole system)
		TS = TEN_TRACE(TA(LEN),[1,2],[3,4])
		TVAL = TS%VALS(1)
	ELSE ! 2 <= L <= LEN-2
		! estimate trial W
		W = MAKE_W3(S, B(LEN-L+1), B(LEN-L))
		TS = GET_TS(T, TA(L-1), TB(LEN-L-1)) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		! SVD split W, update A, S
		CALL SVD(W,[1,2],[3,4],A(L),B(LEN-L),S,MAX_CUT,MAX_ERR)
		! update TA
		TA(L) = GET_TX(TA(L-1), T, A(L))
	END IF
END SUBROUTINE FDMRG
! estimate 2-site trial state
FUNCTION MAKE_W5(A, B, S, SI) RESULT(W)
	TYPE(TENSOR), INTENT(IN) :: A, B, S, SI
	TYPE(TENSOR) :: W
	
	W = TEN_PROD(TEN_PROD(TEN_PROD(S,B,[2],[3]),SI,[2],[1]),TEN_PROD(S,A,[2],[3]),[3],[2])
END FUNCTION MAKE_W5
! estimate 1-site trial state
FUNCTION MAKE_W3(S, X1, X2) RESULT(W)
	TYPE(TENSOR), INTENT(IN) :: S, X1, X2
	TYPE(TENSOR) :: W
	
	W = TEN_PROD(TEN_PROD(S,X1,[2],[3]),X2,[2],[3])
END FUNCTION MAKE_W3
! construct 2-block-2-site system tensor
FUNCTION GET_TS(T, TA, TB) RESULT (TS)
	TYPE(TENSOR), INTENT(IN) :: T, TA, TB
	TYPE(TENSOR) :: TS
	
	TS = TEN_PROD(TEN_PROD(TA,T,[4],[2]),TEN_PROD(TB,T,[4],[2]),[2,6],[2,6])
END FUNCTION GET_TS
! construct 1-block-2-site system tensor (at boundary)
FUNCTION GET_TS1(T, TA) RESULT (TS)
	TYPE(TENSOR), INTENT(IN) :: T, TA
	TYPE(TENSOR) :: TS
	
	TS = TEN_PROD(TEN_PROD(TEN_PROD(TA,T,[4],[2]),EYE_TEN([1,1])),T,[2,6],[2,4])
END FUNCTION GET_TS1
! update TX given T and X
FUNCTION GET_TX(TX, T, X) RESULT (TX1)
! input: TX - block tensor, T - site tensor (MPO), X - projector
! return the enlarged block tensor TX1
	TYPE(TENSOR), INTENT(IN) :: TX, T, X
	TYPE(TENSOR) :: TX1
	
	! zipper-order contraction algorithm
	TX1 = TEN_PROD(TEN_PROD(TEN_PROD(TEN_CONJG(X),TX,[1],[1]),X,[4],[1]),T,[1,4,5],[1,2,3])
END FUNCTION GET_TX
! ----------- Solvers ------------
! anneal the state W to fix point of TS
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
	! collect leg-combined inds in TS and W
	LINDS = COLLECT_INDS(TS,[1,3,5,7])
	RINDS = COLLECT_INDS(TS,[2,4,6,8])
	WINDS = COLLECT_INDS(W,[1,2,3,4])
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(V(0:DIM-1,N)) ! note the rank matching 0:DIM-1
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
			WORKL(LWORKL), WORKEV(3*NCV), RWORK(NCV), SEL(NCV))
		! collect leg-combined inds in TS and W
		LINDS = COLLECT_INDS(TS,[1,3,5,7])
		RINDS = COLLECT_INDS(TS,[2,4,6,8])
		WINDS = COLLECT_INDS(W,[1,2,3,4])
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
					W0(0:N-1) => WORKD(IPNTR(1):IPNTR(1)+N-1)
					W1(0:N-1) => WORKD(IPNTR(2):IPNTR(2)+N-1)
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
			END SELECT
			STOP
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
! 2-site correlation of MPS
FUNCTION CORR2(X, O1, O2) RESULT (CORR)
! input: X - MPS tensor, O1, O2 -  observables
! output: CORR - correlation function
	USE MODEL
	TYPE(TENSOR), INTENT(IN) :: X(LEN), O1, O2
	COMPLEX :: CORR(LEN)
	! local tensors
	TYPE(TENSOR) :: D2
	! local variables
	INTEGER :: L
	
	! assuming X(:) is left-canonical, go backward
	D2 = MAKE_D(X(LEN), O2) ! make initial D2
	! now D2 carries the measurement of O2 at the boundary
	CORR = Z0 ! clear CORR
	DO L = LEN-1, 1, -1
		! measure O1 on site L with D2
		CORR(L) = MEASURE(D2, X(L), O1)
		! grow D2
		IF (L > 1) D2 = GROW_D(D2, X(L))
	END DO
END FUNCTION CORR2
! make double tensor
FUNCTION MAKE_D(X, O) RESULT(D)
! input: X - MPS tensor, O - observable matrix
! output: D - double tensor, D = XH*O*X
! if O is absent, treat O = 1
	TYPE(TENSOR), INTENT(IN) :: X
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: O
	TYPE(TENSOR) :: D
	
	IF (PRESENT(O)) THEN
		D = TEN_PROD(TEN_PROD(TEN_CONJG(X),O,[2],[1]),X,[2,3],[3,2])
	ELSE ! O is absent
		D = TEN_PROD(TEN_CONJG(X),X,[2,3],[2,3])
	END IF
END FUNCTION MAKE_D
! grow double tensor
FUNCTION GROW_D(D, X, O) RESULT(D1)
! input: D - double tensor, X - MPS tensor
! output: D1 - new double tensor, D1 = (XH*O*X)*D
! if O is absent, treat O = 1
	TYPE(TENSOR), INTENT(IN) :: D, X
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: O
	TYPE(TENSOR) :: D1
	
	IF (PRESENT(O)) THEN
		D1 = TEN_PROD(TEN_PROD(TEN_PROD(TEN_CONJG(X),O,[2],[1]),X,[3],[2]),D,[2,4],[1,2])
	ELSE ! O is absent
		D1 = TEN_PROD(TEN_PROD(TEN_CONJG(X),D,[3],[1]),X,[2,3],[2,3])
	END IF
END FUNCTION GROW_D
! take measurement
FUNCTION MEASURE(D, X, O) RESULT (R)
! input: D - double tensor, X - MPS tensor, O - observable mat
! output: R - result of measurement R = Tr (XH*O*X)*D
! if O is absent, treat O = 1
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: D, X
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: O
	COMPLEX :: R
	! local tensor
	TYPE(TENSOR) :: TR
	
	IF (PRESENT(O)) THEN
		TR = TEN_PROD(TEN_PROD(TEN_PROD(TEN_CONJG(X),O,[2],[1]),X,[1,3],[1,2]),D,[1,2],[1,2])
	ELSE ! O is absent
		TR = TEN_PROD(TEN_PROD(TEN_CONJG(X),X,[1,2],[1,2]),D,[1,2],[1,2])
	END IF
	! take the measurement value
	IF (SIZE(TR%VALS) == 0) THEN ! if TR is empty
		R = Z0 ! result is zero
	ELSE ! if TR is finite
		R = TR%VALS(1) ! put the result
	END IF
END FUNCTION MEASURE
! ----------- Debug ------------
! show status after DMRG step
SUBROUTINE SHOW_DMRG_STATUS(MODE,L,TVAL,SVALS)
	USE MODEL
	CHARACTER, INTENT(IN) :: MODE
	INTEGER, INTENT(IN) :: L
	COMPLEX, INTENT(IN) :: TVAL
	COMPLEX, INTENT(IN) :: SVALS(:)
	! local variables
	REAL, ALLOCATABLE :: S(:)
	INTEGER :: I1, I2, NS
	CHARACTER(2) :: JN
	COMPLEX :: LGT
	
	SELECT CASE (MODE)
		CASE ('I')
			I1 = L
			I2 = L+1
			JN = '<>'
		CASE ('F')
			I1 = L
			I2 = L+1
			JN = '->'
			IF (I2 == LEN+1) I2 = 1
		CASE ('B')
			I1 = LEN-L
			I2 = LEN-L+1
			JN = '<-'
			IF (I1 == 0) I1 = LEN
	END SELECT
	IF (L == LEN) THEN
		S = [1.]
	ELSE
		S = REALPART(SVALS)
	END IF
	NS = SIZE(SVALS)
	LGT = LOG(TVAL)
	IF (NS <= 6) THEN
		WRITE (*,'(I3,A,I3,A,F10.6,A,F6.3,A,6F6.3)') I1,JN,I2,': (', REALPART(LGT),',', IMAGPART(LGT), ') S:', S
	ELSE
		WRITE (*,'(I3,A,I3,A,F10.6,A,F6.3,A,4F6.3,A,F6.3)') I1,JN,I2,': (', REALPART(LGT),',', IMAGPART(LGT), ') S:', S(1:4), '  ... ', S(NS)
	END IF
END SUBROUTINE SHOW_DMRG_STATUS
! end of module PHYSICS
END MODULE PHYSICS
! ################ TASK ####################
MODULE TASK
	USE PHYSICS
CONTAINS
! ------------ Data --------------
! collect data
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	INTEGER :: I,J
	CHARACTER, PARAMETER :: BACK(21) = [(CHAR(8),I=1,21)], BLANK(21) = [(' ',I=1, 21)]
	
	WRITE (*,'(A)',ADVANCE='NO') 'Arnoldi iteration:'
	DO I = 1, 5
		WRITE (*,'(I3)',ADVANCE='NO') I
		CALL SLEEP(1)
		WRITE (*,'(3A1)',ADVANCE='NO') BACK(1:3)
	END DO
	WRITE (*,'(21A1)',ADVANCE='NO') BACK
	WRITE (*,'(21A1)',ADVANCE='NO') BLANK
	WRITE (*,'(21A1)',ADVANCE='NO') BACK
	WRITE (*,'(A)') 'ABC'
END SUBROUTINE TEST
! test MPO
SUBROUTINE TEST_MPO()
	TYPE(TENSOR) :: T
	
	CALL SET_MPO(T)
!	CALL TEN_PRINT(T)
	CALL TEN_SAVE('T',TEN_TRANS(T,[2,4,1,3]))
END SUBROUTINE TEST_MPO
! test DMRG
SUBROUTINE TEST_DMRG()
	USE MATHIO
	USE MODEL
	TYPE(TENSOR) :: T, A(LEN), B(LEN), S
	COMPLEX :: C(LEN)
	INTEGER :: ITER
	REAL :: T0, T1
	
	CALL SET_MPO(T)
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA
	CALL CPU_TIME(T0)
	CALL DMRG(T, A, B)
	CALL CPU_TIME(T1)
	PRINT *, T1-T0
!	C = CORR2(A, PAULI_MAT([3]), PAULI_MAT([3]))
!	CALL EXPORT('C',C)
END SUBROUTINE TEST_DMRG
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	INTEGER :: I
	PRINT *, '------------ DMRG -------------'
		
!	CALL TEST()
!	CALL TEST_MPO()
	CALL TEST_DMRG()
END PROGRAM MAIN