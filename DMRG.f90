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
	COMPLEX :: Q = EXP(ZI*PI)
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
	INTEGER :: DCUT
	
! ++++++++ set the vertex tensor here ++++++++
	X =  TENSOR([2,2,2,2],[1,2,3,4,6,7,8,9,11,12,13,14],[Q**0.25,Q**0.25,Z1,Q**0.25,Z1,Q**(-0.25),Q**0.25,Z1,Q**(-0.25),Z1,Q**(-0.25),Q**(-0.25)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! symm SVD of X to unitary U and diagonal S
	CALL SYSVD(X,[1,4],[2,3],U,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	U = TEN_PROD(U,S,[3],[1]) ! attach the half S to U
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	T = TEN_PROD(TEN_PROD(Y,U,[2],[1]),TEN_PROD(Y,U,[2],[2]),[2,3],[3,2])
END SUBROUTINE SET_MPO
! ----------- DMRG -----------
! DMRG Kernel
SUBROUTINE DMRG(T,G,L,DCUT)
! perform infinite-size DMRG to find the max MPS of MPO
! input: T - MPO tensor
! output: G - entanglement tensor, L - Schmidt spectrum
! parameter: DCUT - internal dim cutoff
	USE CONST
	TYPE(TENSOR), INTENT(IN)  :: T ! input MPO tensor
	TYPE(TENSOR), INTENT(OUT) :: G, L ! output MPS tensors
	INTEGER, INTENT(IN) :: DCUT
	! local tensors
	TYPE(TENSOR) :: TE, A, LP, LF, W, TS
	! local variables
	COMPLEX :: TVAL
	INTEGER :: ITER
	! parameters
	INTEGER, PARAMETER :: MAX_ITER = 7
	
	! initialize tensors
	CALL DMRG_INITIALIZATION(TE, A, LP, LF)
	DO ITER = 1, MAX_ITER
		! estimate starting state
		W = STATE_ESTIMATE(A, LP, LF)
		LP = LF ! LP has been used, update to LF
		! construct system T tensor
		TS = MAKE_SYSTEM(TE, T)
		! anneal the state W to the fixed point of TS
!		IF (ITER == 2) THEN
!			CALL TEN_PRINT(TS)
!			CALL TEN_PRINT(W)
!		END IF
		TVAL = ANNEAL(TS, W)
		PRINT *, TVAL
		! SYSVD decompose W to A and LF
		CALL SYSVD(W,[1,2],[3,4],A,LF,DCUT)
		! update environment tensor TE
		CALL UPDATE_ENVIRONMENT(TE, T, A)
		! check convergence
		! now LP and LF holds the successive Schmidt spectrums
		CALL TEN_PRINT(TEN_PROD(LF,LF,[1],[1]))
	END DO
END SUBROUTINE DMRG
! initialization routine
SUBROUTINE DMRG_INITIALIZATION(TE,A,LP,LF)
! called by DMRG
	USE CONST
	TYPE(TENSOR), INTENT(OUT) :: TE, A, LP, LF
	
	TE = TENSOR([1,1,4],[0],[Z1])
	A  = TENSOR([1,2,1],[0,1],[Z1,Z1]/SQRT(2.))
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
	TB = TEN_PROD(TE,T,[3],[2])
	! contract block tensor to system tensor
	TS = TEN_PROD(TB,TB,[5],[5])
END FUNCTION MAKE_SYSTEM
! anneal the state W to fix point of TS
FUNCTION ANNEAL(TS, W) RESULT (TVAL)
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
END FUNCTION ANNEAL
! update TE given T and A
SUBROUTINE UPDATE_ENVIRONMENT(TE, T, A)
! input: TE - environment tensor, T - site tensor (MPO), A - projector
! on output, TE is updated
	TYPE(TENSOR), INTENT(INOUT) :: TE
	TYPE(TENSOR), INTENT(IN) :: T, A
	! local tensor
	TYPE(TENSOR) :: TE1
	
	! zipper-order contraction algorithm
	TE1 = TEN_PROD(TEN_PROD(TEN_PROD(TE,TEN_CONJG(A),[1],[1]),A,[1],[1]),T,[2,1,4],[1,2,3])
!	IF (ALL(TE1%DIMS == TE%DIMS)) THEN
!		TE = TEN_ADD(TE,TE1)
!		TE%VALS = TE%VALS/2
!	ELSE
!		TE = TE1
!	END IF
	TE = TE1
END SUBROUTINE UPDATE_ENVIRONMENT
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
	TYPE(TENSOR) :: T, G, L
	
	CALL SET_MPO(T)
	CALL DMRG(T, G, L, 6)
END SUBROUTINE TEST
! test MPO
SUBROUTINE TEST_MPO()
	TYPE(TENSOR) :: T
	
	CALL SET_MPO(T)
	CALL TEN_PRINT(T)
	CALL TEN_SAVE('T',T)
END SUBROUTINE TEST_MPO
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	INTEGER :: I
	PRINT *, '------------ DMRG -------------'
		
	CALL TEST()
!	CALL TEST_MPO()
END PROGRAM MAIN