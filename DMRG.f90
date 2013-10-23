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
	REAL    :: BETA = 0.
	REAL    :: THETA = 1.*PI
	INTEGER :: LEN = 8 ! must be even and larger than 4
	INTEGER :: MAX_CUT = 8
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
	COMPLEX :: Q, B
	
! ++++++++ set the vertex tensor here ++++++++
	Q = THETA*ZI
	B = BETA*Z1
	X = TENSOR([2,2,2,2],[1,2,3,4,6,7,8,9,11,12,13,14],[EXP(Q/4.),EXP(Q/4.),Z1,EXP(Q/4.),Z1,EXP(-Q/4.),EXP(Q/4.),Z1,EXP(-Q/4.),Z1,EXP(-Q/4.),EXP(-Q/4.)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! symm SVD of X to unitary U and diagonal S
	CALL SYSVD(X,[1,4],[2,3],U,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	U = TEN_PROD(U,S,[3],[1]) ! attach the half S to U
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	U = TEN_PROD(Y,U,[1],[2])
	T = TEN_PROD(U,U,[1,3],[3,1])
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
	INTEGER, PARAMETER :: MAX_SWEEP = 0
	
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
	DO ITER = 1, MAX_SWEEP
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
! anneal the state W to fix point of TS
FUNCTION ANNEAL(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
!	USE MATHIO
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! parameters
	INTEGER, PARAMETER :: N = 16 ! Krylov space dimension
	INTEGER, PARAMETER :: MAX_ITER = 50 ! max interation
	REAL, PARAMETER :: TOL = 1.E-12 ! allowed error of Tval
	! local variables
	INTEGER :: DIM, I, J, K, ITER, INFO
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:)
	COMPLEX, ALLOCATABLE :: Q(:,:)
	COMPLEX :: TVAL0, H(N,N), E(N),VL(0,0),VR(N,N),WORK(65*N)
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
	ALLOCATE(Q(0:DIM-1,N))
	Q = Z0
	FORALL (I = 1:SIZE(W%INDS))
		Q(WINDS(I),1) = W%VALS(I)
	END FORALL
	Q(:,1) = Q(:,1)/SQRT(DOT_PRODUCT(Q(:,1),Q(:,1))) ! normalize
	! prepare to start Arnoldi iteration
	TVAL0 = Z0
	EXH = .FALSE. ! space exhausted flag
	! use Arnoldi iteration algorithm
	DO ITER = 1, MAX_ITER
		H = Z0 ! initialize Heisenberg matrix
		! construct Krylov space
		DO K = 2, N
			! apply TS to Q(:,K-1) -> Q(:,K)
			Q(:,K) = Z0
			DO I = 1,SIZE(TS%INDS)
				Q(LINDS(I),K) = Q(LINDS(I),K) + TS%VALS(I)*Q(RINDS(I),K-1)
			END DO
			! orthogonalization by stabilized Gram–Schmidt process
			DO J = 1, K-1
				H(J,K-1) = DOT_PRODUCT(Q(:,J),Q(:,K))
				Q(:,K) = Q(:,K) - H(J,K-1)*Q(:,J)
			END DO
			! cal the norm of residual vector
			H(K,K-1) = SQRT(DOT_PRODUCT(Q(:,K),Q(:,K)))
			! if it is vanishing, the Arnoldi iteration has broken down
			IF (ABS(H(K,K-1))<TOL) THEN
				EXH = .TRUE.
				EXIT ! exit the iteration, stop construction of Krylov space
			END IF
			! otherwise, normalize the residual vect to a new basis vect
			Q(:,K) = Q(:,K)/H(K,K-1)
		END DO !K
		! now the Heisenberg matrix has been constructed
		! the action of TS is represented on the basis Q as H
		! call LAPACK to diagonalize H
		CALL ZGEEV('N','V',N,H,N,E,VL,1,VR,N,WORK,65*N,RWORK,INFO)
		! now E holds the eigen vals, and VR the eigen vects
		! find the max abs eigen val
		I = MAXLOC(ABS(E),1,IMAG(E)>-TOL/2)
		TVAL = E(I) ! save it in Tval
		! reorganize the eigen vector
		Q(:,1) = MATMUL(Q,VR(:,I)) ! save to 1st col of Q
		Q(:,1) = Q(:,1)/SQRT(DOT_PRODUCT(Q(:,1),Q(:,1))) ! normalize
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
	W%VALS = Q(:,1)
END FUNCTION ANNEAL
! update TX given T and X
FUNCTION GET_TX(TX, T, X) RESULT (TX1)
! input: TX - block tensor, T - site tensor (MPO), X - projector
! return the enlarged block tensor TX1
	TYPE(TENSOR), INTENT(IN) :: TX, T, X
	TYPE(TENSOR) :: TX1
	
	! zipper-order contraction algorithm
	TX1 = TEN_PROD(TEN_PROD(TEN_PROD(TEN_CONJG(X),TX,[1],[1]),X,[4],[1]),T,[1,4,5],[1,2,3])
END FUNCTION GET_TX
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
	REAL, ALLOCATABLE :: S(:)
	
	S = [2.,2.,2.,2.,2.,1.,1.,1.,0.5,0.5,0.,0.,0.]
	SVD_CUT = 9
	SVD_ERR = 0.
	CALL FIND_DCUT(S, SVD_CUT, SVD_ERR)
	PRINT *, SVD_CUT, SVD_ERR
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
	
	CALL SET_MPO(T)
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA
	CALL DMRG(T, A, B)
	C = CORR2(A, PAULI_MAT([3]), PAULI_MAT([3]))
	CALL EXPORT('C',C)
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