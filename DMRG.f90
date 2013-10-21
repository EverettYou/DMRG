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
!	REAL :: BETA = 0.440687
	REAL :: BETA = 0.5
	REAL :: THETA = 0.*PI
	INTEGER :: DCUT = 8
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
SUBROUTINE DMRG(T, M)
! perform infinite-size DMRG to find the max MPS of MPO
! input: T - MPO tensor
! output: M - MPS tensor
	USE MODEL
	TYPE(TENSOR), INTENT(IN)  :: T ! input MPO tensor
	TYPE(TENSOR), INTENT(OUT) :: M ! output MPS tensors
	! local tensors
	TYPE(TENSOR) :: TE, A, LP, LF, W, TS
	! local variables
	COMPLEX :: TVAL
	INTEGER :: ITER
	! parameters
	INTEGER, PARAMETER :: MAX_ITER = 10
	
	! initialize tensors
	CALL DMRG_INITIALIZATION(T, TE, A, LP, LF)
	DO ITER = 1, MAX_ITER
		! estimate starting state
		W = STATE_ESTIMATE(A, LP, LF)
		LP = LF ! LP has been used, update to LF
		! construct system T tensor
!		CALL TEN_SAVE('TE',TE)
		TS = MAKE_SYSTEM(TE, T)
		CALL TEN_SAVE('TS',TEN_FLATTEN(TS,[1,3,5,7,0,2,4,6,8]))
		! anneal the state W to the fixed point of TS
		TVAL = ANNEAL(TS, W)
		WRITE (*,'(A,2G12.4)') 'Tval = ', REAL(TVAL), IMAG(TVAL)
		! SYSVD decompose W to A and LF
		CALL SYSVD(W,[1,2],[3,4],A,LF,DCUT)
		! update environment tensor TE
		CALL UPDATE_ENVIRONMENT(TE, T, A)
		! check convergence
		! now LP and LF holds the successive Schmidt spectrums
		WRITE (*,'(8F7.3)') ABS(LF%VALS)**2
	END DO
	! suppose the entanglement spectrum converges
	! construct the MPS tensor
	LP%VALS = Z1/LP%VALS
	A = TEN_PROD(LP,A,[2],[1])
	LF%VALS = SQRT(LF%VALS)
	M = TEN_PROD(TEN_PROD(LF,A,[2],[1]),LF,[3],[1])
END SUBROUTINE DMRG
! initialization routine
SUBROUTINE DMRG_INITIALIZATION(T, TE, A, LP, LF)
! called by DMRG
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: T
	TYPE(TENSOR), INTENT(OUT) :: TE, A, LP, LF
	! local variables
	INTEGER :: DP, DI, I
	
	DP = T%DIMS(1)
	DI = T%DIMS(2)
	TE = TENSOR([1,DI,1,DI],[(I,I=0,DI-1)]*(DI+1),[(Z1,I=1,DI)])
!	A  = TENSOR([1,DP,1],[(I,I=0,DP-1)],[(Z1/SQRT(2.),I=1,DP)])
	A  = TENSOR([1,DP,1],[0],[Z1])
	LP = TENSOR([1,1],[0],[Z1])
	LF = LP
END SUBROUTINE DMRG_INITIALIZATION
! estimate starting state
FUNCTION STATE_ESTIMATE(A, LP, LF) RESULT(W)
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: A, LP, LF
	TYPE(TENSOR) :: W
	! local tensors
	TYPE(TENSOR) :: LPI, A1
	
	! cal LP^(-1/2)
	LPI = LP ! move data to LPI
	LPI%VALS = Z1/SQRT(LPI%VALS) ! 1/sqrt
	! contract LPI and LF into A
	A1 = TEN_PROD(TEN_PROD(LF,A,[2],[3]),LPI,[2],[1])
	! construct W tensor
	W = TEN_PROD(A1,A1,[3],[3])
END FUNCTION STATE_ESTIMATE
! construct system T
FUNCTION MAKE_SYSTEM(TE, T) RESULT (TS)
	TYPE(TENSOR), INTENT(IN) :: TE, T
	TYPE(TENSOR) :: TS
	! local tensors
	TYPE(TENSOR) :: TB
	
	! construct block T tensor
	TB = TEN_PROD(TE,T,[4],[2])
	! contract block tensor to system tensor
	TS = TEN_PROD(TB,TB,[2,6],[2,6])
END FUNCTION MAKE_SYSTEM
! anneal the state W to fix point of TS
FUNCTION ANNEAL(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
	USE MATHIO
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
	LINDS = COLLECT_INDS(TS,[2,4,6,8])
	RINDS = COLLECT_INDS(TS,[1,3,5,7])
	WINDS = COLLECT_INDS(W,[1,2,3,4])
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(Q(0:DIM-1,N))
	! dump data from tensor W to vector Q(:,1)
	Q(:,1) = Z0
	FORALL (I = 1:SIZE(W%INDS))
		Q(WINDS(I),1) = W%VALS(I)
	END FORALL
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
				Q(RINDS(I),K) = Q(RINDS(I),K) + TS%VALS(I)*Q(LINDS(I),K-1)
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
!		CALL EXPORT('Q',Q)
!		CALL EXPORT('H',H)
		! now the Heisenberg matrix has been constructed
		! the action of TS is represented on the basis Q as H
		! call LAPACK to diagonalize H
		CALL ZGEEV('N','V',N,H,N,E,VL,1,VR,N,WORK,65*N,RWORK,INFO)
!		PRINT *, E
!		CALL EXPORT('VR',VR)
		! now E holds the eigen vals, and VR the eigen vects
		! find the max abs eigen val
		I = MAXLOC(ABS(E),1,IMAG(E)>-TOL/2)
		TVAL = E(I) ! save it in Tval
		! reorganize the eigen vector
		Q(:,1) = MATMUL(Q,VR(:,I)) ! save to 1st col of Q
		Q(:,1) = Q(:,1)/SQRT(DOT_PRODUCT(Q(:,1),Q(:,1))) ! normalize
!		CALL EXPORT('Q1',Q(:,1))
		! check convergence
		! if space exhausted, or relative error < tol
		IF (EXH .OR. ABS((TVAL-TVAL0)/TVAL) < TOL) THEN
			! Arnoldi iteration has converge
			EXIT ! exit Arnoldi interation
		ELSE ! not converge, next iteration
			TVAL0 = TVAL ! save TVAL to TVAL0
			PRINT *, TVAL
		END IF
	END DO ! next Arnoldi interation
	! if exceed max iteration
	IF (ITER > MAX_ITER) THEN !then power iteration has not converge
		WRITE (*,'(A)') 'ANNEAL::fcov: Arnoldi iteration failed to converge, eigen state has not been found.'
	END IF
	! reconstruct W tensor for output
	W%INDS = [(I,I=0,DIM-1)]
	W%VALS = Q(:,1)
!	CALL TEN_PRINT(W)
END FUNCTION ANNEAL
! anneal the state W to fix point of TS
FUNCTION ANNEAL1(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! local variables
	INTEGER :: DIM, I, ITER
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:)
	COMPLEX, ALLOCATABLE :: W0(:), W1(:)
	COMPLEX :: TVAL0
	! parameters
	INTEGER, PARAMETER :: MAX_ITER = 500 ! max interation
	REAL, PARAMETER :: TOL = 1.E-12 ! allowed error of Tval
	
	! unpack data from tensor
	! collect leg-combined inds in TS and W
	LINDS = COLLECT_INDS(TS,[2,4,6,8])
	RINDS = COLLECT_INDS(TS,[1,3,5,7])
	WINDS = COLLECT_INDS(W,[1,2,3,4])
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate vector for W
	ALLOCATE(W0(0:DIM-1), W1(0:DIM-1))
	! dump data from tensor W to vector W0
	W0 = Z0
	FORALL (I = 1:SIZE(W%INDS))
		W0(WINDS(I)) = W%VALS(I)
	END FORALL
	! use power iteration algorithm
	TVAL0 = Z0
	DO ITER = 1, MAX_ITER
		! apply TS to W0 -> W1
		W1 = Z0
		DO I = 1,SIZE(TS%INDS)
			W1(RINDS(I)) = W1(RINDS(I)) + TS%VALS(I)*W0(LINDS(I))
		END DO
		! cal the diagonal val
!		TVAL = (W1(RINDS(1))/ABS(W1(RINDS(1))))*SQRT(SUM(ABS(W1)**2))
		TVAL = SQRT(SUM(ABS(W1)**2))
		! perform normalization
		W1 = W1/TVAL
		! stabilize update and renormalize
		W0 = W0 + W1
		W0 = W0/SQRT(SUM(ABS(W0)**2))
		! check convergence
		IF (ABS((TVAL-TVAL0)/TVAL) < TOL) THEN ! if relative error < tol
			! power iteration has converge
			EXIT ! exit power interation
		ELSE ! not converge, next iteration
			TVAL0 = TVAL ! save TVAL to TVAL0
		END IF
	END DO
	! if exceed max iteration
	IF (ITER > MAX_ITER) THEN !then power iteration has not converge
		WRITE (*,'(A)') 'ANNEAL::fcov: Power iteration failed to converge, eigen state has not been found.'
	END IF
	! reconstruct W tensor for output
	W%INDS = [(I,I=0,DIM-1)]
	W%VALS = W0
!	CALL TEN_PRINT(W)
END FUNCTION ANNEAL1
! anneal the state W to fix point of TS
FUNCTION ANNEAL0(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! local tensors
	TYPE(TENSOR) :: W1
	! local variables
	COMPLEX :: TVAL0
	INTEGER :: ITER
	! parameters
	INTEGER, PARAMETER :: MAX_ITER = 500 ! max interation
	REAL, PARAMETER :: TOL = 1.E-12 ! allowed error of Tval
	
	! use power iteration algorithm
	TVAL0 = Z0
	DO ITER = 1, MAX_ITER
		! apply TS to W -> W1
		W1 = TEN_PROD(TS,W,[2,4,6,8],[1,2,3,4])
		! cal the diagonal val
		TVAL = (W1%VALS(1)/ABS(W1%VALS(1)))*SQRT(SUM(ABS(W1%VALS)**2))
		! perform gauge-fixing and normalization
		W1%VALS = W1%VALS/TVAL
		! stabilize update and renormalize
		W = TEN_ADD(W1, W)
		W%VALS = W%VALS/SQRT(SUM(ABS(W%VALS)**2))
		! check convergence
		IF (ABS((TVAL-TVAL0)/TVAL) < TOL) THEN ! if relative error < tol
			! power iteration has converge
			RETURN ! return with TVAL, and annealed W
		ELSE ! not converge, next iteration
			TVAL0 = TVAL ! save TVAL to TVAL0
		END IF
	END DO
	! if exceed max iteration, then power iteration has not converge
	WRITE (*,'(A)') 'ANNEAL::fcov: Power iteration failed to converge, eigen state has not been found.'
END FUNCTION ANNEAL0
! update TE given T and A
SUBROUTINE UPDATE_ENVIRONMENT(TE, T, A)
! input: TE - environment tensor, T - site tensor (MPO), A - projector
! on output, TE is updated
	TYPE(TENSOR), INTENT(INOUT) :: TE
	TYPE(TENSOR), INTENT(IN) :: T, A
	! local tensor
	TYPE(TENSOR) :: TE1
	
	! zipper-order contraction algorithm
	TE1 = TEN_PROD(TEN_PROD(TEN_PROD(TEN_CONJG(A),TE,[1],[1]),A,[4],[1]),T,[1,4,5],[1,2,3])
!	IF (ALL(TE1%DIMS == TE%DIMS)) THEN
!		TE = TEN_ADD(TE,TE1)
!		TE%VALS = TE%VALS/2
!	ELSE
!		TE = TE1
!	END IF
	TE = TE1
END SUBROUTINE UPDATE_ENVIRONMENT
! ----------- MEASURE -----------
! correlation of MPS
SUBROUTINE CORR(M,O1,O2)
! input: M - MPS tensor, O1, O2 -  observables
! output:
	TYPE(TENSOR), INTENT(IN) :: M, O1, O2
	! local tensors
	TYPE(TENSOR) :: MC, M0, M1, M2
	
	MC = TEN_CONJG(M)
	M0 = TEN_FLATTEN(TEN_PROD(MC,M,[2],[2]),[1,3,0,2,4])
	M1 = TEN_FLATTEN(TEN_PROD(MC,TEN_PROD(O1,M,[2],[2]),[2],[1]),[1,3,0,2,4])
	M2 = TEN_FLATTEN(TEN_PROD(MC,TEN_PROD(O2,M,[2],[2]),[2],[1]),[1,3,0,2,4])
	CALL TEN_SAVE('M0',M0)
	CALL TEN_SAVE('M1',M1)
	CALL TEN_SAVE('M2',M2)
END SUBROUTINE CORR
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
	TYPE(TENSOR) :: T, TE, A, LP, LF, TS
	
	CALL SET_MPO(T)
	CALL DMRG_INITIALIZATION(T,TE,A,LP,LF)
	TS = MAKE_SYSTEM(TE, T)
	TS = TEN_FLATTEN(TS,[1,3,7,5,0,2,4,8,6])
	CALL TEN_SAVE('TS',TS)
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
	USE MODEL
	TYPE(TENSOR) :: T, M
	
	CALL SET_MPO(T)
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A)') 'DCUT = ', DCUT, ', THETA = ', THETA/PI, '*PI, BETA = ', BETA
	CALL DMRG(T, M)
	CALL CORR(M, PAULI_MAT([3]), PAULI_MAT([3]))
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