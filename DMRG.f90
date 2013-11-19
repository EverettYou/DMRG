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
	! physical
	REAL    :: THETA = 1.*PI ! theta-term
	REAL    :: CROSS = 1.    ! crossing: 1. = allow, 0. = avoid
	REAL    :: BETA = 0.440687    ! inverse temperature 0.440687
	! lattice
	INTEGER :: LEN = 20
	INTEGER :: CIR = 1024
	! algorithm
	INTEGER :: MAX_CUT = 8  ! 16
	REAL    :: MAX_ERR = 0.
	INTEGER :: SWEEPS = 1
END MODULE MODEL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
<<<<<<< HEAD
	! definition of huge tensor
	TYPE HUGE_TENSOR
		TYPE(TENSOR) :: TEN ! tensor content
		REAL :: LEV ! scale by EXP(LEV)
	END TYPE HUGE_TENSOR
	REAL :: SCHUR_TOL
CONTAINS
! ----------- tensor -------------
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
! ------- linear algebra ---------
! orthogonalization
SUBROUTINE ORTHOGONALIZE(V, D)
! in-place orthogonalization
! D - orthogonalized space dimension
	COMPLEX, INTENT(INOUT) :: V(:,:)
	INTEGER, INTENT(INOUT) :: D
	! local variables
	INTEGER :: N, I, J
	COMPLEX :: B
	REAL, PARAMETER :: TOL = 1.E-12
	
	N = SIZE(V,2)
	IF (D > N) D = N
	DO I = 1, D
		DO J = 1, I-1
			B = DOT_PRODUCT(V(:,J),V(:,I))
			V(:,I) = V(:,I) - B * V(:,J)
		END DO
		B = DOT_PRODUCT(V(:,I),V(:,I))
		IF (ABS(B) <= TOL) THEN
			D = I - 1
			RETURN
		ELSE
			V(:,I) = V(:,I)/SQRT(B)
		END IF
	END DO
	D = N
END SUBROUTINE ORTHOGONALIZE
! random orthogonal basis
SUBROUTINE RAND_BASIS(V, D)
! inout: V - hold the basis, in: D - num of basis
	USE CONST
	COMPLEX, INTENT(INOUT) :: V(:,:)
	INTEGER, INTENT(INOUT) :: D
	REAL, ALLOCATABLE :: R(:,:)
	
	ALLOCATE(R(SIZE(V,1),D))
	CALL RANDOM_NUMBER(R)
	R = R - 0.5
	V = Z0
	V(:,:D) = CMPLX(R)
	CALL ORTHOGONALIZE(V, D)
END SUBROUTINE RAND_BASIS
! Arnoldi iteration
SUBROUTINE ARNOLDI(A, N, V, D, VALS, LINDS, RINDS)
! construct Krylov basis V, and represent the action of sparse matrix in A
! A - representation matrix
! N - dim of Krylov basis
! V - Krylov space basis
! D - untouched subspace dim
! VALS, TLINS, TRINDS - vals, l and r inds of T mat
	USE MATHIO
	USE CONST
	COMPLEX, INTENT(INOUT) :: A(:,:)
	COMPLEX, TARGET, INTENT(INOUT) :: V(:,:)
	COMPLEX, INTENT(IN) :: VALS(:)
	INTEGER, INTENT(IN) :: N, LINDS(:), RINDS(:)
	INTEGER, INTENT(INOUT) :: D
	! local variables
	INTEGER :: DIM, I, J, K
	COMPLEX, TARGET, ALLOCATABLE :: W(:)
	COMPLEX, POINTER :: MV(:)
	REAL, PARAMETER :: TOL = 1.E-12
	
	DIM = SIZE(V,1) ! get Hilbert space dimension
	! represent T on the first D basis of V
	A = Z0 ! clear
	I = 1  ! leading vector pointer
	IF (D >= N) THEN
		WRITE(*,'(A)')'ARNOLDI::full: Krylov space is full.'
		K = N + 1 ! no new vector to be added
		D = N ! Krylov space cut off at N
		IF (.NOT.ALLOCATED(W)) THEN
			ALLOCATE(W(DIM)) ! get additional space
			MV => W ! point to that
		END IF
	ELSE ! set new vector pointer
		K = D + 1
		MV => V(:,K) ! point work space to V(K)
	END IF ! here K <= N
	DO WHILE(I <= N .AND. I <= D)
		! apply T: V(K) <- T*V(I)
		CALL MV_ACTION(V(:,I), MV, VALS, LINDS, RINDS)
		DO J = 1, D
			! cal overlap with previous
			A(J,I) = DOT_PRODUCT(V(:,J),MV)
			! project out previous components
			MV = MV - A(J,I)*V(:,J)
		END DO
		IF (K <= N) THEN ! if K is not outside
			! cal norm of residual V(K)
			A(K,I) = SQRT(DOT_PRODUCT(MV, MV))
			! if normalizable
			IF (ABS(A(K,I)) > TOL) THEN
				MV = MV/A(K,I) ! normalize
				D = K ! Krylov space dim enlarged to K
				! if K+1 available, K inc
				IF (K < N) THEN
					K = K + 1 ! K will not exceed N
					MV => V(:,K) ! point to the new V(K)
				ELSE ! V has used up allocate additional work space
					K = N + 1 ! K is out
					IF (.NOT.ALLOCATED(W)) THEN
						ALLOCATE(W(DIM)) ! get additional space
						MV => W ! point to that
					END IF
				END IF
			! if V(K) null, D will not increase
			END IF
		END IF
		I = I + 1 ! next vector
	END DO
END SUBROUTINE ARNOLDI
! apply M to V1 resulting in V2
SUBROUTINE MV_ACTION(V1, V2, VALS, LINDS, RINDS)
! M is given in sparse form,
! M(LINDS,RINDS) = VALS
	USE CONST
	COMPLEX, INTENT(IN) :: V1(:), VALS(:)
	INTEGER, INTENT(IN) :: LINDS(:), RINDS(:)
	COMPLEX, INTENT(INOUT) :: V2(:)
	! local variable
	INTEGER :: NREC, IREC
	
	NREC = SIZE(VALS) ! get num of records
	V2 = Z0 ! clear
	DO IREC = 1,NREC ! for all records
		V2(LINDS(IREC)) = V2(LINDS(IREC)) + VALS(IREC)*V1(RINDS(IREC))
	END DO
END SUBROUTINE MV_ACTION
! Schur decomposition
SUBROUTINE SCHUR(A, V, TOL)
! in-place Schur for A
! returning A in Schur form and unitary V
! such that A -> V*A*V^H
	COMPLEX, INTENT(INOUT) :: A(:,:)
	COMPLEX, ALLOCATABLE, INTENT(OUT) :: V(:,:)
	REAL, OPTIONAL :: TOL
	! local variables
	COMPLEX, ALLOCATABLE :: W(:), WORK(:)
	INTEGER :: N, SDIM, INFO
	REAL, ALLOCATABLE :: RWORK(:)
	LOGICAL, ALLOCATABLE :: BWORK(:)
	
	! set tolerance
	IF (PRESENT(TOL)) THEN
		SCHUR_TOL = TOL
	ELSE
		SCHUR_TOL = 1.E-12
	END IF
	! get size
	N = SIZE(A,1)
	! allocate space
	ALLOCATE(W(N),V(N,N),WORK(2*N),RWORK(N),BWORK(N))
	! call LAPACK for Schur decomposition
	CALL ZGEES('V','S',SEL,N,A,N,SDIM,W,V,N,WORK,2*N,RWORK,BWORK,INFO)
	! truncate and output
	IF (SDIM /= 0 .AND. SDIM /= N) THEN
		V = V(:,1:SDIM)
	END IF
END SUBROUTINE SCHUR
! Schur decomposition selection function
FUNCTION SEL(Z) RESULT(Q)
	COMPLEX, INTENT(IN) :: Z
	LOGICAL :: Q
	
	Q = (ABS(Z)>SCHUR_TOL)
END FUNCTION SEL
! matrix power
RECURSIVE FUNCTION MATPOWER(M, N) RESULT (A)
! return A = M^N
	USE CONST
	COMPLEX, INTENT(IN) :: M(:,:)
	INTEGER, INTENT(IN) :: N
	COMPLEX, ALLOCATABLE :: A(:,:)
	! local variables
	INTEGER :: N1
	COMPLEX, ALLOCATABLE :: B(:,:)
	
	SELECT CASE(N)
		CASE (1)
			A = M
		CASE (2)
			A = MATMUL(M,M)
		CASE (3:)
			IF (MODULO(N,2) == 0) THEN
				N1 = N/2
				B = MATPOWER(M,N1)
				A = MATMUL(B,B)
			ELSE
				N1 = (N-1)/2
				B = MATPOWER(M,N1)
				A = MATMUL(MATMUL(B,B),M)
			END IF
		CASE DEFAULT
			A = M
			WRITE(*,'(A)')'MATPOWER::ipow: zero and negative power is not supported yet.'
	END SELECT
END FUNCTION MATPOWER
! ----------- sort -------------
! ordering (by merge sort)
FUNCTION D_ORDERING(F) RESULT(X)
! return the ordering of F in X
! such that F(X) follows the descending order
	REAL, INTENT(IN) :: F(:) ! input data array
	INTEGER, TARGET  :: X(SIZE(F)) ! output ordering
	! local variables
	INTEGER, TARGET  :: Y(SIZE(F)) ! work space
	INTEGER, POINTER :: A(:), B(:)  ! pointer to X,Y
 	INTEGER :: S(SIZE(F) + 1) ! sector initial marker
	INTEGER :: IS, NS
	INTEGER :: LA0, LA1, LA2, IA1, IA2, IB
	LOGICAL :: AB
	REAL :: F0
	
	! prepare the initial X as the index for F
	X = [(IS, IS = 1, SIZE(F))]
	IF (SIZE(F) <= 1) RETURN ! if size(F) <= 1, no need to sort 
	! the algorithm first partition X into sectors
	! the initial position of each sector in S
	! last element in S = end position + 1
	! initialize sector pointer
	NS = SIZE(F) + 1
	S = [(IS, IS = 1, NS)]! each element is a sector on entrance
	! start merge sort iteration
	AB = .TRUE. ! set AB to T, s.t. A => X, B => Y first
	DO WHILE (NS > 2) !  if there are more than 1 sector to merge
		! point A, B to X, Y alternatively
		IF (AB) THEN
			A => X
			B => Y
		ELSE
			A => Y
			B => X
		END IF
		AB = .NOT. AB ! flip AB
		! start to merge sectors in A (2 by 2) to B
		IB = 1 ! rewind IB to the begining
		! IB++, as data move from A -> B
		DO IS = 2, NS, 2 ! IS proceed 2 by 2
			! merge A(S(IS-1:IS)), A(S(IS:IS+1)) to B(MS:MS+1)
			! prepare splits for A from S table
			LA0 = S(IS-1)
			LA1 = S(IS)
			! the first sector is always well-defined
			! but the second sector may be empty
			IF (IS < NS) THEN ! if next sector available
				LA2 = S(IS+1) ! get split from S table
			ELSE ! when IS = NS, the second sector is empty
				LA2 = LA1 ! set split to LA1, nothing will be taken
			END IF
			! merge A(LA0:LA1-1) and A(LA1:LA2-1)
			! ----[-----)[-----)-------
			!    LA0    LA1   LA2
			IA1 = LA0 ! start IA1 from LA0
			IA2 = LA1 ! start IA2 from LA1
			DO WHILE (IA1 < LA1 .AND. IA2 < LA2)
				IF (F(A(IA1)) > F(A(IA2))) THEN ! F(A(IA1)) is greater
					B(IB) = A(IA1) ! greater first
					IB = IB + 1    ! IB inc
					IA1 = IA1 + 1  ! next
				ELSEIF (F(A(IA1)) < F(A(IA2))) THEN ! F(A(IA2)) is greater
					B(IB) = A(IA2) ! greater first
					IB = IB + 1	   ! IB inc
					IA2 = IA2 + 1  ! next
				ELSE ! F(A(IA1)) == F(A(IA2)), duplicated
					F0 = F(A(IA1)) ! rec F
					! rec all duplicates with order preserved					
					! wind IA1 to next distinct entries
					DO WHILE (F(A(IA1))==F0)
						B(IB) = A(IA1) ! rec duplicated entry
						IB = IB + 1    ! IB inc
						IA1 = IA1 + 1  ! next IA1
						IF (IA1 == LA1) EXIT ! exit if reach the end
					END DO
					! wind IA2 to next distinct entries
					DO WHILE (F(A(IA2))==F0)
						B(IB) = A(IA2) ! rec duplicated entry
						IB = IB + 1    ! IB inc
						IA2 = IA2 + 1  ! next IA2
						IF (IA2 == LA2) EXIT ! exit if reach the end
					END DO
				END IF
			END DO ! IA1 == LA1 or IA2 == LA2
			! now one of the sector has been exausted
			! check which
			IF (IA1 == LA1) THEN ! sector 1 exausted
				! dump the rest of sector 2 -> B
				DO WHILE (IA2 < LA2)
					B(IB) = A(IA2)
					IA2 = IA2 + 1
					IB = IB + 1
				END DO
			ELSE ! IA2 == LA2 , sector 2 exausted
				! dump the rest of sector 1 -> B
				DO WHILE (IA1 < LA1)
					B(IB) = A(IA1)
					IA1 = IA1 + 1
					IB = IB + 1
				END DO
			END IF
			! now the sublists are merged into B
			! the next B sector start with the crrent IB
			! record it into the stack
			S(IS/2 + 1) = IB
		END DO
		! all the sublists in A are merged into B
		NS = IS/2 ! set the new stack top
	END DO
	! now B is sorted
	! copy data to X
	X = B
!	PRINT '(100F5.1)',F(X)
END FUNCTION D_ORDERING
! end of module MATH
END MODULE MATH
! ############## DATAPOOL #################
MODULE DATAPOOL
=======
>>>>>>> parent of 1e75352... Use HUGE_TENSOR and DATAPOOL
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
<<<<<<< HEAD
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
	REAL :: F ! free energy
	! private data storages (only assessable through pointers)
	TYPE(TENSOR), TARGET, PRIVATE :: T(2)
	TYPE(TENSOR), TARGET, ALLOCATABLE, PRIVATE :: W(:)
	TYPE(HUGE_TENSOR), TARGET, ALLOCATABLE, PRIVATE :: D(:)
	! temp
	TYPE(TENSOR) :: TS_OUT
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
END SUBROUTINE DAT_ASSOCIATE
! end of module DATAPOOL
END MODULE DATAPOOL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
=======
>>>>>>> parent of 1e75352... Use HUGE_TENSOR and DATAPOOL
CONTAINS
! ------------ set MPO tensor ---------------
! square lattice MPO
SUBROUTINE SET_MPO(TB)
! output: TB - MPO tensor (by default B-type)
	USE MODEL
	TYPE(TENSOR), INTENT(OUT) :: TB ! MPO tensor output
	! local variables
	TYPE(TENSOR) :: X, Y, U, S
	COMPLEX, ALLOCATABLE :: WA(:,:)
	COMPLEX :: Q, B
	
! ++++++++ set the vertex tensor here ++++++++
	Q = THETA * ZI
	B = BETA * Z1
	X =  TENSOR([2,2,2,2],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[EXP(2*B),EXP(Q/4.),EXP(Q/4.),EXP(Z0),EXP(Q/4.),CROSS*EXP(-2*B - Q/2.),EXP(Z0),EXP(-Q/4.),EXP(Q/4.),EXP(Z0),CROSS*EXP(-2*B + Q/2.),EXP(-Q/4.),EXP(Z0),EXP(-Q/4.),EXP(-Q/4.),EXP(2*B)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! symm SVD of X to unitary U and diagonal S
	CALL SYSVD(X,[1,4],[3,2],U,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	U = TEN_PROD(S,U,[2],[3]) ! attach the half S to U
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	U = TEN_PROD(U,Y,[3],[3])
	TB = TEN_TRANS(TEN_PROD(U,U,[2,4],[4,2]),[1,3,2,4])
END SUBROUTINE SET_MPO
! ----------- DMRG -----------
! DMRG controlling routine
SUBROUTINE DMRG(TB, WA)
! input: T - MPO site tensor
! output: WA - MPS tensor
	USE MODEL
	TYPE(TENSOR), INTENT(IN)  :: TB
	TYPE(TENSOR), INTENT(OUT) :: WA(LEN)
	! local tensor
	TYPE(TENSOR) :: TA, WB(LEN), DA(LEN), DB(LEN), P
	! local variables
	INTEGER :: L, ITER
	COMPLEX :: F, FA(LEN), FB(LEN)
	
	! check validity of system size LEN
	IF (MODULO(LEN,2) == 1 .OR. LEN < 4) THEN
		WRITE (*,'(A)') 'DMRG::xlen: LEN must be even and greater than 4.'
		STOP
	END IF
	! set TA from TB MPO
	TA = TEN_TRANS(TB,[2,1,3,4]) ! given by 1 <-> 2
	! iDMRG (warm up)
	DO L = 1, LEN/2
		CALL IDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
		CALL SHOW_DMRG_STATUS('I',L,F,P)
	END DO
	! fDMRG (first sweep)
	DO L = LEN/2+1, LEN
		CALL FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
		CALL SHOW_DMRG_STATUS('F',L,F,P)
	END DO
	! fDMRG sweeps
	DO ITER = 1, SWEEPS
		! fFDMRG (backward sweep)
		DO L = 1, LEN
			CALL FDMRG(L, TB, TA, DB, DA, FB, FA, WB, WA, P, F)
			CALL SHOW_DMRG_STATUS('B',L,F,P)
		END DO
		! fFDMRG (forward sweep)
		DO L = 1, LEN
			CALL FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
			CALL SHOW_DMRG_STATUS('F',L,F,P)
		END DO
	END DO
END SUBROUTINE DMRG
! infinite-size DMRG step
SUBROUTINE IDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
! perform one step of iDMRG at L
! input: L - lattice position
!        TA,TB - MPO
!        DA, DB - block MPO, FA, FB - level of block MPO
!        WA, WB - MPS, P - Schmidt spectrum 
!        F - free energy
! inout: DA, DB, WA, WB, S0, S will be updated
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: TA, TB
	TYPE(TENSOR), INTENT(INOUT) :: DA(:), DB(:), WA(:), WB(:), P
	COMPLEX, INTENT(INOUT) :: F, FA(:), FB(:)
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
		F = Z0 ! set an initial F
		! set initial block MPO (and rescale)
		DA(1) = TA
		DB(1) = TB
		CALL RESCALE(DA(1), Z0, FA(1))
		CALL RESCALE(DB(1), Z0, FB(1))
	ELSE ! L >= 2
		! estimate trial W
		P0%VALS = Z1/P0%VALS ! cal P0^(-1) -> P0
		W = MAKE_W5(WA(L-1), WB(L-1), P, P0)
		P0 = P ! P0 is used, update to P
		TS = MAKE_TS(TA, TB, DA(L-1), DB(L-1)) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = FA(L-1)+FB(L-1)+LOG(TVAL) ! calculate F
		! SVD split W, and update WA, WB, P
		CALL SVD(W,[1,2],[3,4],WA(L),WB(L),P,MAX_CUT,MAX_ERR)
		! update DA, DB and rescale
		DA(L) = NEW_DX(TA, WA(L), DA(L-1))
		DB(L) = NEW_DX(TB, WB(L), DB(L-1))
		CALL RESCALE(DA(L), FA(L-1), FA(L))
		CALL RESCALE(DB(L), FB(L-1), FB(L))
	END IF
	
	TS_OUT = TS
END SUBROUTINE IDMRG
! finite-size DMRG step (forward convention)
SUBROUTINE FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
! perform one step of fDMRG at L (forward dirction)
! input: L - lattice position
!        TA,TB - MPO
!        DA, DB - block MPO, LA, LB - level of block MPO
!        WA, WB - MPS, P - Schmidt spectrum 
!        F - free energy
! inout: DA, WA, S will be updated
! on output: DB is not touched, WB is destroyed
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: TA, TB
	TYPE(TENSOR), INTENT(INOUT) :: DA(:), DB(:), WA(:), WB(:), P
	COMPLEX, INTENT(INOUT) :: F, FA(:), FB(:)
	! local tensors
	TYPE(TENSOR) :: W, TS
	COMPLEX :: TVAL
	COMPLEX, SAVE :: F0

	IF (L == 1) THEN
		! WA(1), no update, obtained from the last sweep
		! update DA, to restart
		! SSB treatment to be implemented here
		DA(1) = NEW_DX(TA, WA(1))
		CALL RESCALE(DA(1), Z0, FA(1)) ! rescale
		F = F0 ! retrieve F from last save
	ELSEIF (L == LEN-1) THEN
		! estimate trial W
		W = MAKE_W3(P, WB(2), WB(1))
		TS = MAKE_TS(TA, TB, DA(LEN-2)) ! construct TS (boundary)
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = FA(LEN-2)+LOG(TVAL) ! calculate F
		F0 = F ! save F at the last site
		! SVD split W, update WA, P
		CALL SVD(W,[1,2],[3,4],WA(LEN-1),WB(1),P,MAX_CUT,MAX_ERR)
		! update DA and rescale
		DA(LEN-1) = NEW_DX(TA, WA(LEN-1), DA(LEN-2))
		CALL RESCALE(DA(LEN-1), FA(LEN-2), FA(LEN-1))
	ELSEIF (L == LEN) THEN
		! update the ending WA by P*WB
		WA(LEN) = TEN_TRANS(TEN_PROD(P,WB(1),[2],[3]),[1,3,2])
		! update DA and rescale
		DA(LEN) = NEW_DX(TA, WA(LEN), DA(LEN-1))
		CALL RESCALE(DA(LEN), FA(LEN-1), FA(LEN))
		F = FA(LEN) ! return F for the whole lattice
	ELSE ! 2 <= L <= LEN-2
		! estimate trial W
		W = MAKE_W3(P, WB(LEN-L+1), WB(LEN-L))
		TS = MAKE_TS(TA, TB, DA(L-1), DB(LEN-L-1)) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = FA(L-1)+FB(LEN-L-1)+LOG(TVAL) ! calculate F
		! SVD split W, update WA, P
		CALL SVD(W,[1,2],[3,4],WA(L),WB(LEN-L),P,MAX_CUT,MAX_ERR)
		! update DA and rescale
		DA(L) = NEW_DX(TA, WA(L), DA(L-1))
		CALL RESCALE(DA(L), FA(L-1), FA(L))
	END IF
END SUBROUTINE FDMRG
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
! rescale DX
SUBROUTINE RESCALE(DX, FX0, FX)
! input: DX - unscaled block MPO, FX0 - last level
! output: DX - scaled block MPO, FX - this scaling level
	USE CONST
	TYPE(TENSOR), INTENT(INOUT) :: DX
	COMPLEX, INTENT(IN)  :: FX0
	COMPLEX, INTENT(OUT) :: FX
	! local variables
	COMPLEX :: Z
	
	! get the scale of D by evaluation
	Z = EVALUATE(DX)
	IF (ABS(Z) /= 0.) THEN
		! rescale by a factor of Z
		DX%VALS = DX%VALS/Z
		! use it to calculate new level
		FX = FX0 + LOG(Z)
	END IF
END SUBROUTINE RESCALE
! cal entanglement entropy
FUNCTION ENTROPY(P) RESULT (S)
! input: P - Schmidt spectrum
! update: S - entanglement entropy
	TYPE(TENSOR), INTENT(IN) :: P
	REAL :: S
	! local variables
	REAL, ALLOCATABLE :: EVALS(:)
	
	EVALS = REALPART(P%VALS)**2
	S = -SUM(EVALS*LOG(EVALS))
END FUNCTION ENTROPY
! ----------- Solvers ------------
! anneal the state W to the ground state of TS
FUNCTION ANNEAL(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
! return the corresponding eigen value TVAL
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	
<<<<<<< HEAD
	TVAL = ANNEAL0(TS,W)
=======
! +++++++++ choose a solver here +++++++++++
	TVAL = ANNEAL0(TS, W) ! self-made (fas, unstable)
!	TVAL = ANNEAL1(TS, W) ! package (slower, robust)
! ++++++++++++++++++++++++++++++++++++++++++
>>>>>>> parent of 1e75352... Use HUGE_TENSOR and DATAPOOL
END FUNCTION ANNEAL
! anneal the state W to the ground state of TS
FUNCTION ANNEAL0(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
! return the corresponding eigen value TVAL
	USE CONST
	USE MATH
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! parameters
	INTEGER, PARAMETER :: N = 16 ! Krylov space dimension
	INTEGER, PARAMETER :: MAX_ITER = 500 ! max interation
	! local variables
	INTEGER :: DIM, D, DS, I, ITER
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), ORDS(:)
	COMPLEX, ALLOCATABLE :: V(:,:), VS(:,:), V1(:)
	COMPLEX :: A(N,N), FID
	REAL, PARAMETER :: TOL = 1.E-12
	
	! leg indexing (remember to +1)
	LINDS = COLLECT_INDS(TS,[1,3,5,7]) + 1
	RINDS = COLLECT_INDS(TS,[2,4,6,8]) + 1
	! get dim of Hilbert space
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(V(DIM,N),V1(DIM))
	V = Z0 ! clear
	! dump W to the 1st col of V
	CALL TEN_REPRESENT(V(:,1), W)
	V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
	D = 1 ! set subspace dim
	! start Arnoldi iteration
	DO ITER = 1, MAX_ITER
		! Arnoldi iteration to get effective action of T
		IF (D > N/2) D = N/2 ! if D too full, cut to half
		CALL ARNOLDI(A,N,V,D,TS%VALS,LINDS,RINDS)
		! now TS is represented as A(:D,:D) by V(:,:D)
		! Schur decomposition
		CALL SCHUR(A(:D,:D),VS)
		DS = SIZE(VS,2) ! get Schur cutoff dimension
		! update Krylov subspace
		V(:,:DS) = MATMUL(V(:,:D),VS)
		! test convergence by fidelity
		CALL MV_ACTION(V(:,1),V1,TS%VALS,LINDS,RINDS) ! one more action on ground state
		FID = DOT_PRODUCT(V(:,1),V1)/A(1,1)
		PRINT '(I3,G10.3,6F8.3)', ITER, ABS(FID - Z1), A(1,1),A(2,2),A(3,3)
		! get mode ordering according to large power of A
		ORDS = GET_ORDS(A(:DS,:DS)) ! ordering without zeros
		! now A(ORDS,ORDS) is the density matrix kernel
		! and V(:,ORDS) holds the basis
		D = SIZE(ORDS)      ! rec num of modes
		V(:,:D) = V(:,ORDS) ! move these modes to the front
		IF (ABS(FID - Z1)< TOL .AND. ITER > 3) THEN
			! error less than tol, converged
			TVAL = A(1,1) ! return max eigen val
			EXIT ! exit do
		END IF
	END DO ! next Arnoldi iteration
	IF (ITER > MAX_ITER) WRITE(*,'(A)')'ANNEAL2::fcnv: Arnoldi iteration fails to converge, the returned subspace may not be accurate.'
	! now V(:,:D) contains the converged subspace
	! and A(ORDS,ORDS) is the density matrix kernel
	PRINT '(A,I2)', 'D = ', D
	! reconstruct W tensor for output
	W%INDS = [(I,I=0,DIM-1)]
	W%VALS = V(:,1)
END FUNCTION ANNEAL0
! get orderings of modes, according to large power of A
FUNCTION GET_ORDS(A) RESULT(ORDS)
! input: A - Schur matrix
! output: ORDS - orderings (without zeros)
	USE CONST ! Z0
	USE MODEL ! CIR
	USE MATH
	COMPLEX, INTENT(IN) :: A(:,:)
	INTEGER, ALLOCATABLE :: ORDS(:)
	! local variables
	COMPLEX, ALLOCATABLE :: B(:,:)
	REAL, ALLOCATABLE :: C(:), C0
	INTEGER :: I, N
	REAL, PARAMETER :: TOL = 1.E-12
	
	! put down Schur matrix, with normalization
	B = A/A(1,1)
	! power by CIR
	B = MATPOWER(B, CIR)
	! get abs diagonal of B 
	N = SIZE(B,1)    ! first get dim of B
	ALLOCATE(C(N))   ! then prepare space to hold diagonal
	FORALL (I = 1:N) ! fill C with abs diagonal
		C(I) = ABS(B(I,I))
	END FORALL
	! scaled by the max val
	C0 = MAXVAL(C)
	IF (C0 > 1) C = C/C0
	! now ready to order
	ORDS = D_ORDERING(C) ! get ordering
	ORDS = PACK(ORDS, C(ORDS) > TOL) ! exclude zeros from the pick
	! now C(ORDS) returns in descending order without zeros
END FUNCTION GET_ORDS
! ----------- Measure -----------
! measure operators O's on MPS X
FUNCTION MEASURE(WS, OS) RESULT (M)
! input: WS - MPSs, OS - MPOs
! output: M - correlation function
! M is output in terms of a tensor:
! leg <-> operator, index on leg <-> position of operator 
	USE MODEL
	TYPE(TENSOR), INTENT(IN) :: WS(:), OS(:)
	TYPE(TENSOR) :: M
	! local tensors
	TYPE(TENSOR) :: E
	! local variables
	INTEGER :: L, K, N
	
	L = SIZE(WS) ! get the num of sites
	K = SIZE(OS) ! get the num of operators
	! deal with some special cases
	IF (K == 0) THEN ! if no operator
		M = ZERO_TEN([INTEGER::]) ! return null tensor
		RETURN	
	END IF
	IF (L == 0) THEN ! if no lattice
		M = ZERO_TEN([(L,L=1,K)]) ! return null tensor
		RETURN
	END IF
	! now L and K are both non-zero
 	N = BINOMIAL(L, K) ! cal the size of correlation func
 	IF (N == 0) THEN ! if not able to put down operators
 		M = ZERO_TEN([(L,N=1,K)]) ! return null tensor
 		RETURN
 	END IF
 	! now N is also non-zero, the normal case
	! allocate space for M
	ALLOCATE(M%DIMS(K),M%INDS(N),M%VALS(N))
	M%DIMS = L  ! each leg can go through the lattice L
	M%INDS = 0  ! clean up
	M%VALS = Z0 ! clean up
	! carry out recursive measurement
	CALL SET_M(L, M%INDS, M%VALS, WS, OS) ! enter without environment
	! on exit, M has been filled with measurement data, return
END FUNCTION MEASURE
! measurement kernel
RECURSIVE SUBROUTINE SET_M(L0, INDS, VALS, WS, OS, E0)
! input: L0 - total lattice size
!        INDS - indices, VALS - measurement values
!        WS - MPSs, OS - MPOs, E0 - environment
! output: INDS, VALS - modified by new measurement data
! E0 optional: if missing, make it initialized
	INTEGER, INTENT(IN) :: L0
	INTEGER, INTENT(INOUT) :: INDS(:)
	COMPLEX, INTENT(INOUT) :: VALS(:)
	TYPE(TENSOR), INTENT(IN) :: WS(:), OS(:)
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: E0
	! local tensor
	TYPE(TENSOR) :: E
	! local variables
	INTEGER :: L, K, I, N1, N2
		
	L = SIZE(WS) ! get the num of sites
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
				WRITE (*,'(A)') 'SET_M::ornk: rank of observable tensors must be 2 or 4.'
				STOP
		END SELECT
	END IF
	! now E has been set, lay down operators
	IF (K == 1) THEN ! for the last operator OS(1)
		! make direct measurement on the current lattice
		DO I = L, 1, -1
			! record the index of measurement
			INDS(I) = I - 1
			! grow E with OS(1) at WS(I) and evaluate
			VALS(I) = EVALUATE(GROW(E, WS(I), OS(1))) ! rec the data
			! now E can be updated to I-1 lattice
			E = GROW(E, WS(I)) ! by absorbing WS(I)
		END DO
	ELSE
		N1 = SIZE(VALS) ! get the size of correlation function
		DO I = L, 1, -1
			! cal the workspace range
			N2 = N1 ! upper range from previous lower range
			N1 = N1*(I-K)/I ! lower range update
			! grow E with the last OS(K) at WS(I)
			! passing the remaining WS and OS to the lower level SET_M
			! with the workspace bounded by N1+1:N2 
			CALL SET_M(L0,INDS(N1+1:N2),VALS(N1+1:N2),WS(:I-1),OS(:K-1),GROW(E,WS(I),OS(K)))
			! on return, VALS contains the measurement values
			! and INDS contains the indices in the lower level
			! for this level, inds must be lifted by L0 and shifted by (I-1)
			INDS(N1+1:N2) = INDS(N1+1:N2) + (I-1)*L0**(K-1)
			! the measurement has complete
			! now E can be updated to I-1 lattice
			E = GROW(E, WS(I)) ! by absorbing WS(I)
		END DO
	END IF
END SUBROUTINE SET_M
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
FUNCTION EVALUATE(E) RESULT(Z)
! input: E - environment tensor
! output: Z - measurement value of E, by tensor trace
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: E
	COMPLEX :: Z
	! local tensor
	TYPE(TENSOR) :: E0
	
	! trace out the legs of E
	!  ╭───╮
	!  │   3
	! ╭┼ 1 E 2 ╮
	! ││   4   │
	! │╰───╯   │
	! ╰────────╯
	E0 = TEN_TRACE(E,[1,3],[2,4])
	! check that E0 is a number
	IF (SIZE(E0%DIMS) /= 0) THEN
		WRITE (*,'(A)') 'EVALUATE::xleg: not evaluated to a number.'
		STOP
	END IF
	! prepare output
	IF (SIZE(E0%VALS) == 0) THEN ! if E0 empty
		Z = Z0 ! return zero
	ELSE ! if not empty
		Z = E0%VALS(1) ! return the value
	END IF
END FUNCTION EVALUATE
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
! ----------- Debug ------------
! show status after DMRG step
SUBROUTINE SHOW_DMRG_STATUS(MODE,L,F,P)
	USE MODEL
	CHARACTER, INTENT(IN) :: MODE
	INTEGER, INTENT(IN) :: L
	COMPLEX, INTENT(IN) :: F
	TYPE(TENSOR), INTENT(IN) :: P
	! local variables
	REAL :: S
	INTEGER :: I1, I2
	CHARACTER(2) :: JN
	
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
	IF (L == 0 .OR. L == LEN) THEN
		S = 0.
	ELSE
		S = ENTROPY(P)
	END IF
	WRITE (*,'(I3,A,I3,A,F10.6,A,F6.3,A,F10.6,A,F5.2,A)') I1,JN,I2,': F = (', REALPART(F),',', IMAGPART(F), '), S = ', S/LOG(2.), 'bit ', SVD_ERR*100,'%'
END SUBROUTINE SHOW_DMRG_STATUS
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
! collect data
SUBROUTINE COLLECT(BETAS)
	USE MODEL
	REAL, INTENT(IN) :: BETAS(:) ! a list of beta to be sampled
	! local tensors
	TYPE(TENSOR) :: T, WS(LEN)
	TYPE(TENSOR) :: M, OS(2)
	! local variables
	INTEGER :: I, N
	
	! prepare observables to be measured
	DO I = 1,2 ! two operators are both sigma_3
		OS(I) = PAULI_MAT([3])
	END DO
	N = SIZE(BETAS) ! get size of beta list
	DO I = 1, N ! for each beta in the list
		BETA = BETAS(I) ! (remember to set beta!!!)
		CALL SET_MPO(T) ! update MPO, using new beta
		! prepare to launch DMRG
		WRITE (*,'(A,I3,A,F5.2,A,F5.2,A,F5.2)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA, ', crossing = ', CROSS
		! launch DMRG to find ground state MPS given MPO
		CALL DMRG(T, WS)
		! take measurements on the MPS
		M = MEASURE(WS, OS)
		! save result to disk
		CALL TEN_SAVE('./data center (IS)/M'//FILENAME(),M)
	END DO 
END SUBROUTINE COLLECT
! make filename
FUNCTION FILENAME()
! file naming system
! call by COLLECT, will use MODEL parameters
	USE MODEL
	CHARACTER(15) :: FILENAME
	
	WRITE (FILENAME,'(A,I2,A,I1,A,I1,A,I1,A,I5.4)') 'D',MAX_CUT,'S',SWEEPS,'Q',NINT(THETA/PI),'C',NINT(CROSS),'B',NINT(BETA*10000)
END FUNCTION FILENAME
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	TYPE(TENSOR) :: TB
	
<<<<<<< HEAD
	CALL DAT_ALLOCATE()
	CALL SET_MPO()
	DO L = 1, LEN/2
		CALL DMRG_STEP('I')
	END DO
!	CALL TEN_SAVE('TS',TS_OUT)

!	TYPE(TENSOR) :: TS	
	
=======
	CALL SET_MPO(TB)
	CALL TEN_SAVE('TB',TB)
>>>>>>> parent of 1e75352... Use HUGE_TENSOR and DATAPOOL
END SUBROUTINE TEST
! test order
SUBROUTINE TEST_ORD()
	USE MATH
	REAL :: F(5) = [0.,2.,3.,1.,1.01]
	INTEGER, ALLOCATABLE :: X(:)
	
	X = D_ORDERING(F)
	X = PACK(X,F(X)>0.01)
	PRINT '(100I6)', X
	PRINT '(100F6.2)', F(X)
END SUBROUTINE TEST_ORD
! test DMRG
SUBROUTINE TEST_DMRG()
	USE MATHIO
	USE MODEL
	TYPE(TENSOR) :: TB, WA(LEN)
	REAL :: T0, T1
	
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A,F5.2)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA, ', crossing = ', CROSS
	CALL CPU_TIME(T0)
	CALL SET_MPO(TB)
	CALL DMRG(TB, WA)
	CALL CPU_TIME(T1)
	PRINT *, T1-T0
END SUBROUTINE TEST_DMRG
! test MEASURE
SUBROUTINE TEST_MEASURE()
	USE MATHIO
	USE MODEL
	TYPE(TENSOR) :: TB, WA(LEN)
	TYPE(TENSOR) :: MA, OS(2)
	COMPLEX, ALLOCATABLE :: W(:)
	INTEGER :: I
	
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A,F5.2)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA, ', crossing = ', CROSS
	CALL SET_MPO(TB)
	CALL DMRG(TB, WA)
	DO I = 1,2
		OS(I) = PAULI_MAT([3])
	END DO
	CALL TEN_SAVE('WA',WA)
!	W = MPS(WA)
!	CALL EXPORT('W',W)
	MA = MEASURE(WA, OS)
	CALL TEN_SAVE('MA',MA)
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
	
	CALL SET_MPO(T)
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
	ALLOCATE(LBS(2**N))
	LBS = 0
	LBS(1) = 3
	LBS(2**N) = 3 
	S = TEN2MAT(PAULI_MAT(LBS),[1],[2])
	A = M
	A0 = SUM([(A(I,I),I=1,D)])
	A = A/A0
	F = LOG(A0)/LOG(2.)
	DO J = 1, L
		B = MATMUL(S,A)
		Z(1,J) = SUM([(B(I,I),I=1,D)])
		Z(2,J) = F
		A = MATMUL(M,A)
		A0 = SUM([(A(I,I),I=1,D)])
		A = A/A0
		F = F + LOG(A0)/LOG(2.)
	END DO
	CALL EXPORT('Z',Z)
END SUBROUTINE TEST_TRANSF
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	INTEGER :: I
	PRINT *, '------------ DMRG -------------'

	CALL TEST()
!	CALL TEST_ORD()
!	CALL TEST_DMRG()
!	CALL TEST_MEASURE()
!	CALL TEST_TRANSF()
!	CALL COLLECT([0.4406])
END PROGRAM MAIN