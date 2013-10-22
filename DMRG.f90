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
	REAL :: BETA = 0.3
	REAL :: THETA = 0.*PI
	INTEGER :: DCUT = 6
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
! iDMRG Kernel
SUBROUTINE IDMRG(T, L, TL, AL, S)
! perform infinite-size DMRG to find the max MPS of MPO
! input: T - MPO tensor, L - system size
! output: TL - compressed transfer matrix at each level
!         AL - projection operator in TL
!         S  - final entanglement tensor 
	USE MODEL
	TYPE(TENSOR), INTENT(IN)  :: T ! input MPO tensor
	INTEGER, INTENT(IN)       :: L ! system size
	TYPE(TENSOR), INTENT(OUT) :: TL(0:L), AL(0:L), S
	! local tensors
	TYPE(TENSOR) :: S0, W, TS
	! local variables
	COMPLEX :: TVAL
	INTEGER :: DP, DI, I
	
	! initialize tensors
	DP = T%DIMS(1) ! get physical dim
	DI = T%DIMS(2) ! get MPO internal dim
	TL(0) = TENSOR([1,DI,1,DI],[(I,I=0,DI-1)]*(DI+1),[(Z1,I=1,DI)])
	AL(0) = TENSOR([1,DP,1],[(I,I=0,DP-1)],[(Z1/SQRT(2.),I=1,DP)])
	S0 = TENSOR([1,1],[0],[Z1])
	S  = S0
	! start iDMRG iteration, collect tensors
	DO I = 1, L ! loop over all levels
		! estimate the trial state
		W = STATE2(AL(I-1), S0, S)
		S0 = S ! S0 has been used, update to S
		! construct system tensor TS
		TS = TS4(TL(I-1), T)
		! anneal the state W to the fixed point of TS
		TVAL = ANNEAL(TS, W, [1,3,5,7],[2,4,6,8],[1,2,3,4])
		! SYSVD decompose W to A and S
		CALL SYSVD(W,[1,2],[3,4],AL(I),S,DCUT)
		! construct new environment tensor
		TL(I) = GET_TE(TL(I-1), T, AL(I))
		WRITE (*,'(I3,A,G10.4,A,G10.4,A,100F6.3)') I, ': (', REALPART(TVAL),',', IMAGPART(TVAL), ') S:', REALPART(S%VALS)
	END DO ! I
END SUBROUTINE IDMRG
! fDMRG Kernel
SUBROUTINE FDMRG(T, L, TL, AL, S)
! input: T - MPO tensor, L - system size
! inout: initial TL, AL, S are obtained from IDMRG
!        on return, TL, AL, S are updated
	TYPE(TENSOR), INTENT(IN) :: T ! MPO tensor
	INTEGER, INTENT(IN)      :: L ! system size
	TYPE(TENSOR), INTENT(INOUT) :: TL(0:L), AL(0:L), S
	! local tensors
	TYPE(TENSOR) :: M, TS, B
	! local variables
	INTEGER :: DIMO, DIMS, I
	COMPLEX :: TVAL
	
	! initialize trial state
	M = STATE1(AL(L),S)
	! construct initial environment tensor
	DIMO = T%DIMS(2) ! get MPO internal dim
	DIMS = S%DIMS(1) ! get MPS internal dim
	TL(0) = TEN_TRANS(TEN_PROD(EYE_TEN([DIMO,DIMO]),EYE_TEN([DIMS,DIMS])),[3,1,4,2])
	! start fDMRG iteration
	DO I = 1, L ! sweep over the lattice
		! make system tensor TS
		TS = TS3(TL(L-I), T, TL(I-1))
		! anneal M to the ground sate of TS
		TVAL = ANNEAL(TS, M, [1,3,5],[2,4,6],[1,2,3])
		! split M to A and B, also return S
		CALL M2AB(M, AL(I), B, S)
		! construct new environment tensor
		TL(I) = GET_TE(TL(I-1), T, AL(I))
		! estimate next trial state
		IF (I < L) M = STATE1(AL(L-I), B)
		WRITE (*,'(I3,A,G10.4,A,G10.4,A,100F6.3)') I, ': (', REALPART(TVAL),',', IMAGPART(TVAL), ') S:', REALPART(S%VALS)
	END DO
END SUBROUTINE FDMRG
! fDMRG Kernel
SUBROUTINE FDMRG1(T, L, TL, AL, S)
! input: T - MPO tensor, L - system size
! inout: initial TL, AL, S are obtained from IDMRG
!        on return, TL, AL, S are updated
	TYPE(TENSOR), INTENT(IN) :: T ! MPO tensor
	INTEGER, INTENT(IN)      :: L ! system size
	TYPE(TENSOR), INTENT(INOUT) :: TL(0:L), AL(0:L), S
	! local tensors
	TYPE(TENSOR) :: M, TS, B, TN(0:L), AN(0:L)
	! local variables
	INTEGER :: DIMO, DIMS, I
	COMPLEX :: TVAL
	
	! initialize trial state
	M = STATE1(AL(L),S)
	! construct initial environment tensor
	DIMO = T%DIMS(2) ! get MPO internal dim
	DIMS = S%DIMS(1) ! get MPS internal dim
	TN(0) = TEN_TRANS(TEN_PROD(EYE_TEN([DIMO,DIMO]),EYE_TEN([DIMS,DIMS])),[3,1,4,2])
	! start fDMRG iteration
	DO I = 1, L ! sweep over the lattice
		! make system tensor TS
		TS = TS3(TL(L-I), T, TN(I-1))
		! anneal M to the ground sate of TS
		TVAL = ANNEAL(TS, M, [1,3,5],[2,4,6],[1,2,3])
		! split M to A and B, also return S
		CALL M2AB(M, AN(I), B, S)
		! construct new environment tensor
		TN(I) = GET_TE(TN(I-1), T, AN(I))
		! estimate next trial state
		IF (I < L) M = STATE1(AL(L-I), B)
		WRITE (*,'(I3,A,G10.4,A,G10.4,A,100F6.3)') I, ': (', REALPART(TVAL),',', IMAGPART(TVAL), ') S:', REALPART(S%VALS)
	END DO
	! move data
	TL = TN
	AL = AN
END SUBROUTINE FDMRG1
! estimate 2-site trial state
FUNCTION STATE2(A, S0, S) RESULT(W)
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: A, S0, S
	TYPE(TENSOR) :: W
	! local tensors
	TYPE(TENSOR) :: SI, A1
	
	! cal S0^(-1/2)
	SI = S0 ! move data to SI
	SI%VALS = Z1/SQRT(SI%VALS) ! 1/sqrt
	! contract SI and S into A
	A1 = TEN_PROD(TEN_PROD(S,A,[2],[3]),SI,[2],[1])
	! construct W tensor
	W = TEN_PROD(A1,A1,[3],[3])
END FUNCTION STATE2
! estimate 1-site trial state
FUNCTION STATE1(A, B) RESULT(M)
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: A, B
	TYPE(TENSOR) :: M
	
	M = TEN_PROD(A,B,[3],[1])
END FUNCTION STATE1
! construct 2-block-2-site system tensor
FUNCTION TS4(TE, T) RESULT (TS)
	TYPE(TENSOR), INTENT(IN) :: TE, T
	TYPE(TENSOR) :: TS
	! local tensors
	TYPE(TENSOR) :: TB
	
	! construct block T tensor
	TB = TEN_PROD(TE,T,[4],[2])
	! contract block tensor to system tensor
	TS = TEN_PROD(TB,TB,[2,6],[2,6])
END FUNCTION TS4
! construct 2-block-1-site system tensor
FUNCTION TS3(TL, T, TR) RESULT (TS)
	TYPE(TENSOR), INTENT(IN) :: TL, T, TR
	TYPE(TENSOR) :: TS
	
	TS = TEN_PROD(TEN_PROD(TL,T,[4],[2]),TR,[2,6],[4,2])
END FUNCTION TS3
! anneal the state W to fix point of TS
FUNCTION ANNEAL(TS, W, LLEGS, RLEGS, WLEGS) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
	USE MATHIO
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	INTEGER, INTENT(IN) :: LLEGS(:), RLEGS(:), WLEGS(:) 
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
	LINDS = COLLECT_INDS(TS,LLEGS)
	RINDS = COLLECT_INDS(TS,RLEGS)
	WINDS = COLLECT_INDS(W,WLEGS)
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(Q(0:DIM-1,N))
	! dump data from tensor W to vector Q(:,1)
	Q(:,1) = Z0
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
		WRITE (*,'(A)') 'ANNEAL::fcov: Arnoldi iteration failed to converge, eigen state has not been found.'
	END IF
	! reconstruct W tensor for output
	W%INDS = [(I,I=0,DIM-1)]
	WINDS = COLLECT_INDS(W,WLEGS)
	W%VALS = Q(WINDS,1)
END FUNCTION ANNEAL
! update TE given T and A
FUNCTION GET_TE(TE, T, A) RESULT (TE1)
! input: TE - environment tensor, T - site tensor (MPO), A - projector
! on output, TE is updated
	TYPE(TENSOR), INTENT(IN) :: TE, T, A
	TYPE(TENSOR) :: TE1
	
	! zipper-order contraction algorithm
	TE1 = TEN_PROD(TEN_PROD(TEN_PROD(TEN_CONJG(A),TE,[1],[1]),A,[4],[1]),T,[1,4,5],[1,2,3])
END FUNCTION GET_TE
! split M to A and B
SUBROUTINE M2AB(M, A, B, S)
! input: M - state tensor
! output: A - projection tensor, B - link tensor, S - singular value matrix (of B)
	TYPE(TENSOR), INTENT(IN)  :: M
	TYPE(TENSOR), INTENT(OUT) :: A, B, S
	! local variables
	REAL, ALLOCATABLE :: SVAL(:), RWORK(:)
	COMPLEX, ALLOCATABLE :: BA(:,:), U(:,:), VH(:,:), WORK(:)
	INTEGER :: N1, N2, L, LWORK, INFO, I
	
	! dump M tensor to the matrix BA
	BA = TEN2MAT(M,[1],[2,3])
	! prepare to perform SVD to BA
	! get actual size of BA
	N1 = SIZE(BA,1)
	N2 = SIZE(BA,2)
	L = MIN(N1,N2)
	! calculate the size of the workspace
	LWORK = MAX(2*L+MAX(N1,N2),MIN(L*(L+65),2*L+32*(N1+N2)))
	! allocate for workspace
	ALLOCATE(SVAL(L), U(N1,L), VH(L,N2), WORK(LWORK), RWORK(5*L))
	! now everything is ready for SVD, call LAPACK
	CALL ZGESVD('S', 'S', N1, N2, BA, N1, SVAL, U, N1, VH, L, WORK, LWORK, RWORK, INFO)
	DEALLOCATE(BA, WORK, RWORK) ! make room for later
	! now BA = U.diag_mat(SVAL).VH
	! we should have A = VH, S = diag_mat(SVAL), B = U.S
	! construct A from VH
	A = MAT2TEN(VH,M%DIMS([3,2,1]),[3],[2,1])
	S = DIAG_MAT(CMPLX(SVAL))
	B = TEN_PROD(MAT2TEN(U,[N1,L],[1],[2]),S,[2],[1])
END SUBROUTINE M2AB
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
	TYPE(TENSOR) :: T
	
	CALL SET_MPO(T)
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
	INTEGER, PARAMETER :: L = 8
	TYPE(TENSOR) :: T, TL(0:L), AL(0:L), S
	INTEGER :: ITER
	
	CALL SET_MPO(T)
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A)') 'DCUT = ', DCUT, ', THETA = ', THETA/PI, '*PI, BETA = ', BETA
	WRITE (*, '(A)') 'iDMRG'
	CALL IDMRG(T, L, TL, AL, S)
	DO ITER = 1,10
		WRITE (*, '(A)') 'fDMRG'
		CALL FDMRG(T, L ,TL, AL, S)
	END DO
!	CALL CORR(M, PAULI_MAT([3]), PAULI_MAT([3]))
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