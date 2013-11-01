! ###########################################
! TENSORIAL module
! defines a derived type TENSOR to represent sparse tensor.
! The type TENSOR has three components:
! DIMS - dimensions, integer array (/DIM1, DIM2, .../)
! INDS - indices, integer array (/IND1, IND2, .../)
! VALS - values, complex array (/VAL1, VAL2, .../)
! A tensor T with its K th non-zero elements at T(I1, I2, I3, ...)
! will be stored as
! INDS(K) = (I1-1) + (I2-1)*DIM1 + (I3-1)*DIM1*DIM2 + ...
! VALS(K) = T(I1, I2, I3, ...)
! and K runs over all the non-zero elements in T.
! The size of INDS given the num of non-zero elements.
! NOTE:
! * The first entry T(1,1,...) is indexed by IND = 0.
! * INDS is always sorted in ascending order.
! * Scalar is represented as TENSOR([INTEGER::],[0],[val])
MODULE TENSORIAL
	! type definition of sparse TENSOR
	TYPE TENSOR
		INTEGER, ALLOCATABLE :: DIMS(:) ! dimensions
		INTEGER, ALLOCATABLE :: INDS(:) ! indices
		COMPLEX, ALLOCATABLE :: VALS(:) ! values
	END TYPE TENSOR
	REAL, PARAMETER :: TOL = 1.E-5 	! global tolerance level
	! recent SVD record
	INTEGER :: SVD_CUT = 0  ! SVD cut dimension
	REAL    :: SVD_ERR = 0. ! SVD truncation error
	INTERFACE TEN_PRINT
		MODULE PROCEDURE TEN_PRINT_0
		MODULE PROCEDURE TEN_PRINT_1
	END INTERFACE
	INTERFACE TEN_SAVE
		MODULE PROCEDURE TEN_SAVE_0
		MODULE PROCEDURE TEN_SAVE_1
	END INTERFACE
	INTERFACE TEN_LOAD
		MODULE PROCEDURE TEN_LOAD_0
		MODULE PROCEDURE TEN_LOAD_1
	END INTERFACE
CONTAINS
! Tensor Construction ----------------------
! return zero tensor of dimension DIMS
FUNCTION ZERO_TEN(DIMS) RESULT(TEN)
	INTEGER, INTENT(IN) :: DIMS(:) ! dimensions of the tensor
	TYPE(TENSOR) :: TEN ! tensor to return
	
	TEN = TENSOR(DIMS,[INTEGER::],[COMPLEX::])
END FUNCTION ZERO_TEN
! return diagonal matrix from DIAG
FUNCTION DIAG_MAT(DIAG) RESULT(MAT)
	COMPLEX, INTENT(IN) :: DIAG(:) ! diagonal elements input
	TYPE(TENSOR) :: MAT ! diagonal matrix to return
	! local variable
	INTEGER :: N, M, I
	REAL    :: LB
	INTEGER :: LOCS(SIZE(DIAG))
	
	N = SIZE(DIAG) ! get dimension
	MAT%DIMS = [N,N] ! set dimension
	LB = TOL * NORM(DIAG) ! set truncation
	! elements in DIAG with abs below LB will be considered as zero
	! get locs of non-zero elements in DIAG
	M = 0
	DO I = 1, N
		IF (ABS(DIAG(I))>LB) THEN
			M = M + 1
			LOCS(M) = I
		END IF
	END DO
	! now LOCS(1:M) holds the locs of non-zeros in DIAG
	MAT%INDS = (LOCS(1:M)-1)*(N+1) ! set inds
	MAT%VALS = DIAG(LOCS(1:M)) ! move vals
END FUNCTION DIAG_MAT
! return Pauli matrix according to the labels LBS
FUNCTION PAULI_MAT(LBS) RESULT(MAT)
! direct product of Pauli matrix
! example: if LBS = (/1,2,3/)
! will return sigma_1 ox sigma_2 ox sigma_3 as 8*8 matrix
	INTEGER, INTENT(IN) :: LBS(:) ! labels to input
	TYPE(TENSOR) :: MAT ! Pauli matrix as TENSOR to return
	! local variables
	INTEGER :: NLB, N, L, LL, I
	COMPLEX, PARAMETER :: Z1 = (1.,0.), ZI = (0.,1.)
	INTEGER, PARAMETER :: IJS(2,0:1,0:3) = &
		RESHAPE([[[0,0],[1,1]], &
				 [[1,0],[0,1]], &
				 [[1,0],[0,1]], &
				 [[0,0],[1,1]]],[2,2,4])
	! IJS(r/c, entry, label) returns the row or column index (in 0, 1)
	! of non-zero entries in elementary Pauli matrix
	! r/c: 1 - row index, 2 - column index
	! entry: 0 - first entry, 1 - second entry
	! label: 0,1,2,3 for 4 elementary Pauli matrices
	COMPLEX, PARAMETER :: ZS(0:1,0:3) = &
		RESHAPE([[Z1,Z1],[Z1,Z1],[ZI,-ZI],[Z1,-Z1]],[2,4])
	! ZS(entry, label) returns the value of non-zero entries
	! in elementary Pauli matrix
	! entry: 0 - first entry, 1 - second entry
	! label: 0,1,2,3 for 4 elementary Pauli matrices

	NLB = SIZE(LBS) ! num of labels
	N = 2**NLB      ! size of Pauli matrix
	MAT%DIMS = (/N,N/) ! set dimensions
	! num of non-zeros = size = N
	! initialize with auto reallocation
	MAT%INDS = SPREAD(0,1,N)  ! INDS = 0
	MAT%VALS = SPREAD(Z1,1,N) ! VALS = Z1
	! loop over each non-zero entry
	! each entry is labeled by an integer L as binary number
	! e.g. L = 01101 in a 5-labeled Pauli matrix
	! the I th bit labels the entry in the I th elementary Pauli matrix
	! 0 - first entry of elementary Pauli matrix
	! 1 - second entry of elementary Pauli matrix
	DO L = 0, N-1
		LL = L ! LL is temp, use to get each bit of L
		! I runs over elementary Pauli matrices (by labels)
		DO I = 0, NLB-1
			! ind is accumulated for each I by
			! + 2^I * (row_ind + size N * col_ind)
			MAT%INDS(L+1) = MAT%INDS(L+1) + &
				2**I * SUM((/1,N/) * IJS(:, MOD(LL,2), LBS(NLB-I)))
			! val is multiplied for each I by
			MAT%VALS(L+1) = MAT%VALS(L+1) * ZS(MOD(LL,2), LBS(NLB-I))
			! for every next label, bit ++ by dividing LL by 2
			LL = LL/2
		END DO !I = 0, NLB-1
	END DO !L = 0, N-1
END FUNCTION PAULI_MAT
! return identity tensor of dimension DIMS
FUNCTION EYE_TEN(DIMS, VAL) RESULT (TEN)
	INTEGER, INTENT(IN) :: DIMS(:) ! dimensions of the tensor
	COMPLEX, OPTIONAL :: VAL
	TYPE(TENSOR) :: TEN ! tensor to return
	! local variables
	INTEGER :: ILEG, PD, SKP, I, N
	
	! calculate the skip
	SKP = 1 ! initialize with 1
	! iteratively get all dimension product
	PD = 1 ! starting from the first dim
	DO ILEG = 1, SIZE(DIMS) - 1
		PD = PD * DIMS(ILEG)
		SKP = SKP + PD
	END DO ! ILEG
	! construct identity tensor
	TEN%DIMS = DIMS ! set dims
	N = MINVAL(DIMS) ! num of diagonal elements
	TEN%INDS = [(I-1,I=1,N)]*SKP
	IF (PRESENT(VAL)) THEN
		TEN%VALS = [(VAL,I=1,N)]	
	ELSE
		TEN%VALS = [((1.,0.),I=1,N)]
	END IF
END FUNCTION EYE_TEN
! Representation ---------------------------
! represent sparse tensor TEN by dense tensor TAR
SUBROUTINE TEN_REPRESENT(TAR, TEN)
! * caller must make sure TEN%DIMS consistent with the shape of TAR
	TYPE(TENSOR), INTENT(IN) :: TEN ! sparse tensor input
	COMPLEX, INTENT(INOUT)   :: TAR(0:PRODUCT(TEN%DIMS)-1) ! dense tensor to be added to
	
	TAR = (0.,0.) ! clear TAR
	TAR(TEN%INDS) = TEN%VALS ! simply put into TAR
END SUBROUTINE TEN_REPRESENT
! reconstruct sparse tensor TEN from dense tensor TAR
SUBROUTINE TEN_ABSTRACT(TAR, TEN)
! * caller must make sure TEN%DIMS consistent with the shape of TAR
! on input, TEN must contain the dimension info TEN%DIMS
	TYPE(TENSOR), INTENT(INOUT) :: TEN ! sparse tensor output
	COMPLEX, INTENT(IN) :: TAR(0:PRODUCT(TEN%DIMS)-1) ! dense tensor for abstract
	! local variables
	INTEGER :: I
	
	! first dump all data from TAR to TEN
	TEN%INDS = [(I,I=0,PRODUCT(TEN%DIMS)-1)] ! making the index list
	TEN%VALS = TAR(0:PRODUCT(TEN%DIMS)-1) ! move data
	! then call TEN_CHOP to eliminate zero elements in TEN
	CALL TEN_CHOP(TEN)
END SUBROUTINE TEN_ABSTRACT
! extract component of sparse tensor TEN in dense tensor TAR
FUNCTION TEN_EXTRACT(TAR, TEN) RESULT(VAL)
! * caller must make sure TEN%DIMS consistent with the shape of TAR
	TYPE(TENSOR), INTENT(IN) :: TEN
	COMPLEX, INTENT(IN)      :: TAR(0:PRODUCT(TEN%DIMS)-1)
	COMPLEX :: VAL
	! local variable
	INTEGER :: I
	
	! cal inner product
	VAL = (0.,0.)
	DO I = 1, SIZE(TEN%INDS)
		VAL = VAL + TAR(TEN%INDS(I))*CONJG(TEN%VALS(I))
	END DO
	! normalization
	VAL = VAL / SUM(ABS(TEN%VALS)**2)
END FUNCTION TEN_EXTRACT
! add sparse tensor TEN to dense tensor TAR with coefficient VAL
SUBROUTINE TEN_DEPOSIT(TAR, TEN, VAL)
! * caller must make sure TEN%DIMS consistent with the shape of TAR
	TYPE(TENSOR), INTENT(IN) :: TEN ! sparse tensor input
	COMPLEX, INTENT(INOUT)   :: TAR(0:PRODUCT(TEN%DIMS)-1) ! dense tensor to be added to
	COMPLEX, OPTIONAL, INTENT(IN) :: VAL ! coefficient
	
	IF (PRESENT(VAL)) THEN ! if coefficient presented
		TAR(TEN%INDS) = TAR(TEN%INDS) + VAL*TEN%VALS ! add with coefficient
	ELSE ! if coefficient is not specified
		TAR(TEN%INDS) = TAR(TEN%INDS) + TEN%VALS ! simply add to TAR
	END IF
END SUBROUTINE TEN_DEPOSIT
! Arithmetics ------------------------------
! conjugate tensor
FUNCTION TEN_CONJG(TEN) RESULT(TEN0)
! elementary complex conjugate of a tensor
	TYPE(TENSOR), INTENT(IN) :: TEN
	TYPE(TENSOR) :: TEN0
	
	TEN0 = TEN ! copy tensor to result
	TEN0%VALS = CONJG(TEN0%VALS) ! take conjugate
END FUNCTION TEN_CONJG
! addition of two sparse tensors
FUNCTION TEN_ADD(TEN1, TEN2) RESULT(TEN0)
	TYPE(TENSOR), INTENT(IN) :: TEN1, TEN2 ! input tensors to be added
	TYPE(TENSOR) :: TEN0 ! resulting tensor for output
	! local variables
	INTEGER :: NREC1, NREC2, NREC0, IREC1, IREC2, IREC0
	INTEGER, ALLOCATABLE :: INDS(:)
	COMPLEX, ALLOCATABLE :: VALS(:)
	REAL :: EPS
	
	! check input tensor
	CALL TEN_CHECK('TEN_ADD', TEN1)
	CALL TEN_CHECK('TEN_ADD', TEN2)
	! check tensor dimension consistency
	IF (.NOT.ALL(TEN1%DIMS==TEN2%DIMS)) THEN
		WRITE (*,'(A)') 'TEN_ADD::xdim: Attempt to add tensors of different dimensions.'
		STOP
	END IF
	! set dimension
	TEN0%DIMS = TEN1%DIMS
	! get num of records
	NREC1 = SIZE(TEN1%INDS)
	NREC2 = SIZE(TEN2%INDS)
	NREC0 = NREC1 + NREC2 
	! set eps, val < eps will be treated as 0
	EPS = MAX(NORM(TEN1%VALS),NORM(TEN2%VALS))*TOL
	! prepare temp storage for resulting INDS and VALS
	ALLOCATE(INDS(NREC0),VALS(NREC0))
	! start merge iteration
	IREC1 = 1
	IREC2 = 1
	IREC0 = 0
	! IREC1, IREC2 runs over records in TEN1 and TEN2
	DO WHILE (IREC1 <= NREC1 .AND. IREC2 <= NREC2)
		IF (IREC0 == 0 .OR. ABS(VALS(IREC0)) > EPS) THEN
		! IREC0 inc only if it is at a non-zero position
			IREC0 = IREC0 + 1
		END IF
		IF (TEN1%INDS(IREC1) < TEN2%INDS(IREC2)) THEN
			! if IREC1 is currently the smaller indexed record
			! move TEN1 data to temp storage
			INDS(IREC0) = TEN1%INDS(IREC1)
			VALS(IREC0) = TEN1%VALS(IREC1)
			IREC1 = IREC1 + 1
		ELSEIF (TEN1%INDS(IREC1) > TEN2%INDS(IREC2)) THEN
			! if IREC2 is currently the smaller indexed record
			! move TEN2 data to temp storage
			INDS(IREC0) = TEN2%INDS(IREC2)
			VALS(IREC0) = TEN2%VALS(IREC2)
			IREC2 = IREC2 + 1
		ELSE
			! if IREC1 and IREC2 is equally indexed
			! TEN1 + TEN2 data to temp storage
			INDS(IREC0) = TEN1%INDS(IREC1) ! here TEN1%INDS(IREC1)==TEN2%INDS(IREC2)
			VALS(IREC0) = TEN1%VALS(IREC1) + TEN2%VALS(IREC2)
			IREC1 = IREC1 + 1
			IREC2 = IREC2 + 1
		END IF
	END DO ! upon exit, one of TEN1, TEN2 has records exhausted
	! dump the remaining data to temp storage
	! but before dumping, should check if current rec is zero
	! if so, the record pointer -1
	IF (ABS(VALS(IREC0)) < EPS) IREC0 = IREC0 - 1
	DO WHILE (IREC1 <= NREC1) ! will do if TEN1 still remains
		IREC0 = IREC0 + 1
		INDS(IREC0) = TEN1%INDS(IREC1)
		VALS(IREC0) = TEN1%VALS(IREC1)
		IREC1 = IREC1 + 1
	END DO
	DO WHILE (IREC2 <= NREC2) ! will do if TEN2 still remains
		IREC0 = IREC0 + 1
		INDS(IREC0) = TEN2%INDS(IREC2)
		VALS(IREC0) = TEN2%VALS(IREC2)
		IREC2 = IREC2 + 1
	END DO
	! now all data added into temp storage
	! num of record = IREC0
	! auto allocate, and move data to TEN0
	TEN0%INDS = INDS(1:IREC0)
	TEN0%VALS = VALS(1:IREC0)
	! free temp storage
	DEALLOCATE(INDS,VALS)
END FUNCTION TEN_ADD
! summation of sparse tensor array
RECURSIVE FUNCTION TEN_SUM(TENS) RESULT(TEN)
! caller must make sure every tensor in TENS has the same DIMS
! use divide and conquer algorithm 
	TYPE(TENSOR), INTENT(IN) :: TENS(1:) ! an array of tensors to sum up
	TYPE(TENSOR) :: TEN ! the resulting tensor of the sum
	! local variables
	INTEGER :: NTEN
	
	NTEN = SIZE(TENS) ! num of tensors to sum
	IF (NTEN > 1) THEN
		! if there are more than one tensors in TENS
		! recursively partition TENS for summation
		IF (MOD(NTEN,2)==0) THEN ! even partition
			TEN = TEN_ADD(TEN_SUM(TENS(1:NTEN/2)),TEN_SUM(TENS(NTEN/2+1:NTEN)))
		ELSE ! odd partition
			TEN = TEN_ADD(TEN_SUM(TENS(1:(NTEN-1)/2)),TEN_SUM(TENS((NTEN+1)/2:NTEN)))		
		END IF
	ELSE ! if only one tensor in TENS
		TEN = TENS(1) ! return that tensor
	END IF
END FUNCTION TEN_SUM
! tensor linear combination
FUNCTION TEN_COMBINE(TENS, COEFS) RESULT(TEN)
! linear combine tensors by coefs
! user must make sure that COEFS has the same size as TENS
	COMPLEX, INTENT(IN) :: COEFS(:)
	TYPE(TENSOR), INTENT(IN) :: TENS(:)
	TYPE(TENSOR) :: TEN
	! local variables
	TYPE(TENSOR) :: TS(SIZE(TENS))
	INTEGER :: NTEN, ITEN, NT
	REAL :: EPS
	
	! get num of tensors
	NTEN = SIZE(TENS)
	! set eps, coef below eps will be treated as 0
	EPS = NORM(COEFS)*TOL
	! num of non-zero tensors
	NT = 0
	! multiply coefs to tens
	DO ITEN = 1, NTEN ! go through each tensor
		IF (ABS(COEFS(ITEN)) > EPS) THEN
		! if coef not 0
			NT = NT+1
			TS(NT) = TENSOR(TENS(ITEN)%DIMS,TENS(ITEN)%INDS,COEFS(ITEN)*TENS(ITEN)%VALS)
		END IF
	END DO ! ITEN
	! sum up
	TEN = TEN_SUM(TS(1:NT))
END FUNCTION TEN_COMBINE
! tensor trace (contraction)
FUNCTION TEN_TRACE(TEN, LEGS1, LEGS2) RESULT(TEN0)
! contract LEGS1 with LEGS2 of tensor TEN
! retun the resulting tensor as TEN0
! NOTE:
! * LEGS1 and LEGS2 must have the same size
! EXAMPLE:
! * Trace of a matrix: TEN_TRACE(mat,[1],[2])
	TYPE(TENSOR), INTENT(IN) :: TEN ! tensor input
	INTEGER, INTENT(IN) :: LEGS1(:), LEGS2(:) ! legs to contract
	TYPE(TENSOR) :: TEN0 ! tensor output
	! local variables
	INTEGER :: NREC
	INTEGER, ALLOCATABLE :: RLEGS(:) ! remaining legs
	INTEGER, ALLOCATABLE :: LINDS1(:), LINDS2(:), RINDS(:) ! inds
	LOGICAL, ALLOCATABLE :: CRASH(:)
	INTEGER, ALLOCATABLE :: INDS(:), ORD(:) ! temp INDS and its ordering
	COMPLEX, ALLOCATABLE :: VALS(:) ! temp VALS
	
	! check input tensor
	CALL TEN_CHECK('TEN_TRACE', TEN)
	! prepare indices and ordering
	! get remaining legs by excluding leading legs 1 and 2
	RLEGS = REMAINING_LEGS(SIZE(TEN%DIMS),[LEGS1,LEGS2])
	! split index list for leading and remaining legs
	LINDS1 = COLLECT_INDS(TEN, LEGS1)
	LINDS2 = COLLECT_INDS(TEN, LEGS2)
	RINDS  = COLLECT_INDS(TEN, RLEGS)
	! put the dimensions (taken from remaining legs)
	TEN0%DIMS = TEN%DIMS(RLEGS)
	! get mask of index crashing
	CRASH = (LINDS1 == LINDS2) ! crash if two leading ind match
	DEALLOCATE(RLEGS,LINDS1,LINDS2) ! not ref later
	! pack the inds and vals to temp storage
	INDS = PACK(RINDS, CRASH)
	VALS = PACK(TEN%VALS, CRASH)
	DEALLOCATE(CRASH) ! not ref later
	NREC = SIZE(INDS) ! get num of records
	! del duplicates in INDS while adding associated vals together
	CALL MERGE_DUPLICATES(INDS, VALS, NREC)
	! move data to TEN0 by auto-realloc
	TEN0%INDS = INDS(:NREC)
	TEN0%VALS = VALS(:NREC)
END FUNCTION TEN_TRACE
! tensor product (including contraction)
FUNCTION TEN_PROD(TEN1, TEN2, LEGS1, LEGS2) RESULT(TEN0)
! tensor product of TEN1 and TEN2
! with LEGS1 of TEN1 contracted with LEGS2 of TEN2
! return the resulting tensor as TEN0
! * caller must make sure TEN1%DIMS(LEG1) and TEN2%DIMS(LEG2) match
! EXAMPLES: (for matrices mat1, mat2)
! 1 Tensor product : TEN_PROD(mat1,mat2) -> rank-4 tensor
! 2 Matrix multiplication: TEN_PROD(mat1,mat2,[2],[1]) -> rank-2 (matrix)
! 3 Trace contraction: TEN_PROD(mat1,mat2,[2,1],[1,2]) -> rank-0 (scalar)
	TYPE(TENSOR), INTENT(IN) :: TEN1, TEN2 ! input tensor for contraction
	INTEGER, OPTIONAL, INTENT(IN) :: LEGS1(:), LEGS2(:) ! legs to be connected
	TYPE(TENSOR) :: TEN0 ! resulting tensor
	! local variable
	INTEGER :: NREC1, NREC2, IREC1, IREC2, NREC0, IREC0, IREC20
	INTEGER :: RDIM1, RDIM2, LDIM, PDIM, LIND, I, INFO
	INTEGER, ALLOCATABLE :: RLEGS1(:), RLEGS2(:) ! remaining legs
	INTEGER, ALLOCATABLE :: LINDS1(:), LINDS2(:) ! leading inds
	INTEGER, ALLOCATABLE :: RINDS1(:), RINDS2(:) ! remaining inds
	INTEGER, ALLOCATABLE :: ORD1(:), ORD2(:) ! orderings of leading inds
	INTEGER, ALLOCATABLE :: INDS(:) ! temp INDS
	COMPLEX, ALLOCATABLE :: VALS(:) ! temp VALS
	LOGICAL, ALLOCATABLE :: MASK(:) ! mask for non-zeros
	
!	PRINT *, 'TEN_PROD'
	! check input tensor
	CALL TEN_CHECK('TEN_PROD', TEN1)
	CALL TEN_CHECK('TEN_PROD', TEN2)
	! check leg-rank and dimension consistency
	IF (PRESENT(LEGS1)) THEN
		! given LEGS1, check its legs not exceeding the rank
		IF (.NOT.(ALL(LEGS1 >= 1) .AND. ALL(LEGS1 <= SIZE(TEN1%DIMS)))) THEN
			WRITE (*,'(A)') 'TEN_PROD::xleg: leg label lying outside the tensor rank.'
			STOP
		END IF
		IF (PRESENT(LEGS2)) THEN
			! given LEGS2, check its legs not exceeding the rank
			IF (.NOT.(ALL(LEGS2 >= 1) .AND. ALL(LEGS2 <= SIZE(TEN2%DIMS)))) THEN
				WRITE (*,'(A)') 'TEN_PROD::xleg: leg label lying outside the tensor rank.'
				STOP
			END IF
			! given both LEGS1, LEGS2, check dimension match
			IF (.NOT. ALL(TEN1%DIMS(LEGS1) == TEN2%DIMS(LEGS2))) THEN
				WRITE (*,'(A)') 'TEN_PROD::xdim: contracting dimensions mismatch.'
				STOP
			END IF
		ELSE ! given LEGS1, missing LEGS2, impossible to contract
			WRITE (*,'(A)') 'TEN_PROD::xdim: contracting dimensions mismatch.'
			STOP
		END IF
	ELSE
		IF (PRESENT(LEGS2)) THEN ! given LEGS2, missing LEGS1, impossible to contract
			WRITE (*,'(A)') 'TEN_PROD::xdim: contracting dimensions mismatch.'
			STOP
		END IF
		! missing both LEGS1, LEGS2 is ok
		! and the dimensions always match in this case
	END IF
	! prepare indices and ordering
	IF (PRESENT(LEGS1) .AND. SIZE(LEGS1)>0) THEN
		! get the remaining legs
		RLEGS1 = REMAINING_LEGS(SIZE(TEN1%DIMS), LEGS1)
		! split index list for leading and remaining legs
		LINDS1 = COLLECT_INDS(TEN1, LEGS1)
		RINDS1 = COLLECT_INDS(TEN1, RLEGS1)
		! get ordering of leading indices by merge sort
		ORD1 = ORDERING(LINDS1)
	ELSE ! if leading leg not presented
		RLEGS1 = [(I,I=1,SIZE(TEN1%DIMS))] ! all legs are remaining
		LINDS1 = SPREAD(0,1,SIZE(TEN1%INDS)) ! leading inds = all 0
		RINDS1 = TEN1%INDS ! remaining inds = inds
		ORD1 = [(I,I=1,SIZE(TEN1%INDS))] ! default ordering 
	END IF
	IF (PRESENT(LEGS2) .AND. SIZE(LEGS2)>0) THEN
		! get the remaining legs
		RLEGS2 = REMAINING_LEGS(SIZE(TEN2%DIMS), LEGS2)
		! spit index list for leading and remaining legs
		LINDS2 = COLLECT_INDS(TEN2, LEGS2)
		RINDS2 = COLLECT_INDS(TEN2, RLEGS2)
		! get ordering of leading indices by merge sort
		ORD2 = ORDERING(LINDS2)
	ELSE ! if leading leg not presented
		RLEGS2 = [(I,I=1,SIZE(TEN2%DIMS))] ! all legs are remaining
		LINDS2 = SPREAD(0,1,SIZE(TEN2%INDS)) ! leading inds = all 0
		RINDS2 = TEN2%INDS ! remaining inds = inds
		ORD2 = [(I,I=1,SIZE(TEN2%INDS))] ! default ordering 
	END IF
	! put the dimensions
	TEN0%DIMS = [TEN1%DIMS(RLEGS1),TEN2%DIMS(RLEGS2)]
	! total dim of remaining legs in TEN1 and TEN2
	RDIM1 = PRODUCT(TEN1%DIMS(RLEGS1)) ! will be used later in IND calculation
	RDIM2 = PRODUCT(TEN2%DIMS(RLEGS2))
	LDIM = PRODUCT(TEN1%DIMS)/RDIM1 ! total dim of leading legs
	DEALLOCATE (RLEGS1,RLEGS2) ! free space
	! get num of records
	NREC1 = SIZE(TEN1%INDS)
	NREC2 = SIZE(TEN2%INDS)
	! allocate space for temp storage
	! the worst case is that all LINDS1 = LINDS2 = the same index
	! in this case there will be totally NREC1*NREC2 crashes
	! however dense storage requires RDIM1*RDIM2 space
	! if NREC1*NREC2 is worse than RDIM1*RDIM2, switch to the dense method
	IF (1.*NREC1*NREC2 < 1.*RDIM1*RDIM2) THEN ! use REAL to prevent overflow
		! if record pairs < dense storage, use sparse method
		PDIM = NREC1*NREC2 ! max possible total dims
!		PRINT *, 'spa', NREC1, NREC2, PDIM
		ALLOCATE(INDS(PDIM),VALS(PDIM),STAT=INFO) ! allocate for space by PDIM
		IF (INFO /= 0) THEN ! if failed to allocate
			WRITE (*,'(A)') 'TEN_PROD::fail: fail to allocate space for new tensor.'
			STOP
		END IF
		! search for index crash
		IREC1 = 1
		IREC2 = 1
		IREC0 = 0
		DO WHILE (IREC1 <= NREC1 .AND. IREC2 <= NREC2)
			! run through both records
			IF (LINDS1(ORD1(IREC1)) < LINDS2(ORD2(IREC2))) THEN
				! if LINDS1 is smaller
				IREC1 = IREC1 + 1 ! IREC1 inc
			ELSEIF (LINDS1(ORD1(IREC1)) > LINDS2(ORD2(IREC2))) THEN
				! if LINDS2 is smaller
				IREC2 = IREC2 + 1 ! IREC2 inc
			ELSE ! now LINDS1 == LINDS2, a crash happens
				LIND = LINDS1(ORD1(IREC1)) ! rec the leading index
				IREC20 = IREC2 ! rec initial IREC2
				! perform double-iteration among degenerated leading indices
				DO WHILE (IREC1 <= NREC1 .AND. LINDS1(ORD1(IREC1)) == LIND)
					IREC2 = IREC20 ! rewind IREC2
					DO WHILE (IREC2 <= NREC2 .AND. LINDS2(ORD2(IREC2)) == LIND)
						IREC0 = IREC0 + 1
						INDS(IREC0) = RINDS1(ORD1(IREC1)) + RDIM1 * RINDS2(ORD2(IREC2))
						VALS(IREC0) = TEN1%VALS(ORD1(IREC1)) * TEN2%VALS(ORD2(IREC2))
						IREC2 = IREC2 + 1
					END DO ! IREC2 iteration
					IREC1 = IREC1 + 1
				END DO ! IREC1 iteration
			END IF
		END DO ! upon return, one LINDS must have been exhausted
		! the remaining LINDS has nothing to crash with, discard
		DEALLOCATE(LINDS1,LINDS2,RINDS1,RINDS2,ORD1,ORD2) ! free the space
		! now INDS(1:IREC0) and VALS(1:IREC0) holds the indices and values
		NREC0 = IREC0 ! update num of records 
!		PRINT *, NREC0, '<', PDIM
!		PRINT *, INDS(:NREC0)
!		PRINT *, VALS(:NREC0)
		! del duplicates in INDS while adding associated vals together
		CALL MERGE_DUPLICATES(INDS, VALS, NREC0)
		! move data to TEN0 by auto allocation
		TEN0%INDS = INDS(:NREC0)
		TEN0%VALS = VALS(:NREC0)
	ELSE ! if record pairs > dense storage, switch to dense method
		PDIM = RDIM1*RDIM2 ! dense storage
!		PRINT *, 'den', RDIM1, RDIM2, PDIM
		! initialization
		ALLOCATE(TEN0%INDS(PDIM),TEN0%VALS(PDIM),STAT=INFO)
		IF (INFO /= 0) THEN ! if failed to allocate
			WRITE (*,'(A)') 'TEN_PROD::fail: fail to allocate space for new tensor.'
			STOP
		END IF
		TEN0%INDS = [(I,I=0,PDIM-1)]
		TEN0%VALS = (0.,0.)
		! search for index crash
		IREC1 = 1
		IREC2 = 1
!		IREC0 = 0
		DO WHILE (IREC1 <= NREC1 .AND. IREC2 <= NREC2)
			! run through both records
			IF (LINDS1(ORD1(IREC1)) < LINDS2(ORD2(IREC2))) THEN
				! if LINDS1 is smaller
				IREC1 = IREC1 + 1 ! IREC1 inc
			ELSEIF (LINDS1(ORD1(IREC1)) > LINDS2(ORD2(IREC2))) THEN
				! if LINDS2 is smaller
				IREC2 = IREC2 + 1 ! IREC2 inc
			ELSE ! now LINDS1 == LINDS2, a crash happens
				LIND = LINDS1(ORD1(IREC1)) ! rec the leading index
				IREC20 = IREC2 ! rec initial IREC2
				! perform double-iteration among degenerated leading indices
				DO WHILE (IREC1 <= NREC1 .AND. LINDS1(ORD1(IREC1)) == LIND)
					IREC2 = IREC20 ! rewind IREC2
					DO WHILE (IREC2 <= NREC2 .AND. LINDS2(ORD2(IREC2)) == LIND)
!						IREC0 = IREC0 + 1
						I = 1 + RINDS1(ORD1(IREC1)) + RDIM1 * RINDS2(ORD2(IREC2))
						TEN0%VALS(I) = TEN0%VALS(I) + TEN1%VALS(ORD1(IREC1)) * TEN2%VALS(ORD2(IREC2))
						IREC2 = IREC2 + 1
					END DO ! IREC2 iteration
					IREC1 = IREC1 + 1
				END DO ! IREC1 iteration
			END IF
		END DO ! upon return, one LINDS must have been exhausted
		! the remaining LINDS has nothing to crash with, discard
		! chop of zero elements in TEN0
		CALL TEN_CHOP(TEN0)
	END IF
!	PRINT *, LDIM, 'crash', IREC0
END FUNCTION TEN_PROD
! get the remaining legs
FUNCTION REMAINING_LEGS(NLEG, LLEGS) RESULT(RLEGS)
! called by TEN_TR, TEN_PROD
! given tensor rank and leading legs, return remaining legs
	INTEGER, INTENT(IN)  :: NLEG ! total num of legs (= tensor rank)
	INTEGER, INTENT(IN)  :: LLEGS(:) ! leading legs
	INTEGER, ALLOCATABLE :: RLEGS(:) ! remaining legs
	! local variables
	LOGICAL :: RQ(NLEG) ! RQ: T - remaining, F - leading 
	INTEGER :: LEG
	
	! check leg labels
	IF (.NOT.(ALL(LLEGS >= 1) .AND. ALL(LLEGS <= NLEG))) THEN
		WRITE (*,'(A)') 'REMAINING_LEGS::xleg: leg label lying outside the tensor rank.'
		STOP
	END IF
	RQ = .TRUE. ! assuming all legs are remaining
	! turn off the flag for leading legs
	RQ(LLEGS) = .FALSE.
	! now only the remaining legs have RQ = .T.
	! collect the remaining legs
	RLEGS = PACK([(LEG,LEG = 1, NLEG)],RQ)
END FUNCTION REMAINING_LEGS
! Transformations ---------------------------
! chop zeros in the tensor
SUBROUTINE TEN_CHOP(TEN, EPS0)
	TYPE(TENSOR), INTENT(INOUT) :: TEN
	REAL, OPTIONAL, INTENT(IN) :: EPS0 ! the chopping level
	! local variables
	REAL :: EPS
	LOGICAL, ALLOCATABLE :: MASK(:)
	
	IF (PRESENT(EPS0)) THEN ! if the chopping level provided
		EPS = EPS0 ! use it
	ELSE ! if not, set appropriate chopping level
		EPS = NORM(TEN%VALS)*TOL
	END IF
	! get the mask for non-zero elements
	MASK = ABS(TEN%VALS)>EPS
	! pack data according to the mask
	! INDS and VALS are auto reallocated to the reduced size (F95)
	TEN%INDS = PACK(TEN%INDS,MASK)
	TEN%VALS = PACK(TEN%VALS,MASK)
END SUBROUTINE TEN_CHOP
! flatten tensor to matrix
FUNCTION TEN2MAT(TEN, LLEGS, RLEGS0) RESULT(MAT)
! matrix: row <- leading legs, col <- remaining legs
! this is more efficient than TEN_FLATTEN
	TYPE(TENSOR), INTENT(IN) :: TEN ! tensor input
	INTEGER, INTENT(IN) :: LLEGS(:) ! leading legs
	INTEGER, OPTIONAL, INTENT(IN) :: RLEGS0(:)
	COMPLEX, ALLOCATABLE :: MAT(:,:)
	! local variables
	INTEGER, ALLOCATABLE :: RLEGS(:), LINDS(:), RINDS(:)
	INTEGER :: LDIM, RDIM, I
	
	! check input tensor
	CALL TEN_CHECK('TEN2MAT', TEN)
	IF (PRESENT(RLEGS0)) THEN ! if the remaining leg is specified
		RLEGS = RLEGS0 ! use the specification
	ELSE ! if not specified
		! get remaining legs according to the leading legs
		RLEGS = REMAINING_LEGS(SIZE(TEN%DIMS),LLEGS)
	END IF
	! collect row and col inds
	! +1 to shift the inds to start with 1
	LINDS = COLLECT_INDS(TEN, LLEGS)+1
	RINDS = COLLECT_INDS(TEN, RLEGS)+1
	! cal dims of the matrix
	LDIM = PRODUCT(TEN%DIMS(LLEGS))
	RDIM = PRODUCT(TEN%DIMS(RLEGS))
	! allocate the matrix
	ALLOCATE(MAT(LDIM,RDIM))
	MAT = (0.,0.) ! initialize with zeros
	! fill in the matrix
	FORALL (I = 1:SIZE(TEN%INDS))
		MAT(LINDS(I),RINDS(I)) = TEN%VALS(I)
	END FORALL
END FUNCTION TEN2MAT
! restore tensor from matrix
FUNCTION MAT2TEN(MAT, DIMS, LLEGS, RLEGS) RESULT(TEN)
! construct TEN from MAT by row -> leading legs, col -> remaining legs
! the resulting tensor dims is specified as DIMS
! * caller must make sure the total dim is consistent with the matrix dim
	COMPLEX, INTENT(IN) :: MAT(:,:)
	INTEGER, INTENT(IN) :: LLEGS(:), DIMS(:)
	INTEGER, OPTIONAL, INTENT(IN) :: RLEGS(:)
	TYPE(TENSOR) :: TEN
	! local variables
	INTEGER :: DIM, I
	INTEGER, ALLOCATABLE :: LEGS(:), INDS(:), ORD(:)
	
	DIM = PRODUCT(DIMS) ! get total dimension
	! contruct TEN by dumping data from MAT
	! the tensor legs will be arraged as [LLEGS, RLEGS]
	IF (PRESENT(RLEGS)) THEN ! if the remaining leg is specified
		LEGS = [LLEGS, RLEGS]! use the specification
	ELSE ! if not specified
		! get remaining legs according to the leading legs
		LEGS = [LLEGS, REMAINING_LEGS(SIZE(DIMS),LLEGS)]
	END IF
	TEN%DIMS = DIMS(LEGS)
	TEN%INDS = [(I,I=0,DIM-1)]
	TEN%VALS = RESHAPE(MAT,[DIM])
	! chop off zero elements in TEN
	CALL TEN_CHOP(TEN)
	! perform tensor transpose to restore the leg ordering
	TEN = TEN_TRANS(TEN,ORDERING(LEGS))
END FUNCTION MAT2TEN
! tensor transpose
FUNCTION TEN_TRANS(TEN, LEGS) RESULT(TEN0)
! perform in-place rearrangement of legs in the order of LEGS
! such that T0(l1,l2,...) = T(LEGS(l1),LEGS(l2),...)
! NOTES:
! * LEGS must have same size as TEN%DIMS
! * To transpose a matrix, use TEN_TRANS(mat,[2,1])
	TYPE(TENSOR), INTENT(IN) :: TEN ! tensor input
	INTEGER, INTENT(IN) :: LEGS(:) ! new legs
	TYPE(TENSOR) :: TEN0 ! tensor output
	! local variables
	INTEGER, ALLOCATABLE :: INDS(:), ORD(:)
	
	! collect new inds and get ordering
	INDS = COLLECT_INDS(TEN, LEGS)
	ORD = ORDERING(INDS)
	! make transposed tensor
	TEN0 = TENSOR(TEN%DIMS(LEGS), INDS(ORD), TEN%VALS(ORD))
END FUNCTION TEN_TRANS
! tensor flatten (leg combination)
FUNCTION TEN_FLATTEN(TEN, LEGS) RESULT(TEN0)
! groups of legs are separated by 0 
! LEGS = [l11,l12,...,0,l21,l22,...,0,...]
! given i, combining all legs lij (for all j)
! to make the leg i in the result
	TYPE(TENSOR), INTENT(IN) :: TEN ! tensor input
	INTEGER, INTENT(IN) :: LEGS(:) ! leg grouping instruction
	TYPE(TENSOR) :: TEN0 ! resulting tensor
	! local variable
	INTEGER :: NLEG, ILEG1, ILEG2, NGRP, DIM
	INTEGER, ALLOCATABLE :: INDS(:), ORD(:)
	
	NLEG = SIZE(LEGS) ! get instruction length
	TEN0%DIMS = [INTEGER::] ! initialize DIMS by empty
	INDS = SPREAD(0,1,SIZE(TEN%INDS)) ! initialize INDS to 0
	! go through the instructions
	! ILEG1 and ILEG2 marks the head and end of each group of legs
	ILEG2 = NLEG ! assuming the last group ends at NLEG 
	DO WHILE (ILEG2 >= 1) ! back-to-front reading
		IF (LEGS(ILEG2) == 0) THEN ! if end is at 0 position
			ILEG2 = ILEG2 - 1 ! end dec
			CYCLE ! cycle to the prev non-zero
		END IF ! if the end is at non-0
		! now end is at the end of a group
		ILEG1 = ILEG2 - 1 ! set the head right before the end
		! wind the head to the next 0
		DO WHILE (ILEG1 >= 1 .AND. LEGS(ILEG1)/=0)
			! if the head is not at a 0
			ILEG1 = ILEG1 - 1 ! head dec
			CYCLE ! cycle the head to 0
		END DO
		! now head is at the first 0 before the group begins
		! or at position 0 (outside LEGS)
		! the group of legs is scoped by ILEG1+1:ILEG2
		! update dim to the dim of current group of legs
		DIM = PRODUCT(TEN%DIMS(LEGS(ILEG1+1:ILEG2)))
		! prepend to new tensor dims
		TEN0%DIMS = [DIM, TEN0%DIMS]
		! shift the prev inds by dim, add with new inds
		INDS = DIM*INDS + COLLECT_INDS(TEN, LEGS(ILEG1+1:ILEG2))
		ILEG2 = ILEG1 - 1 ! set the end before current head
		! initiate another round of searching for group ending
	END DO ! now INDS holds the flattened inds
	! get its ordering
	ORD = ORDERING(INDS)
	! scatter the data to new tensor
	TEN0%INDS = INDS(ORD)
	TEN0%VALS = TEN%VALS(ORD)
END FUNCTION TEN_FLATTEN
! Tensor SVDs -------------------------------
! Standard SVD
SUBROUTINE SVD(TEN, LEGS1, LEGS2, TEN1, TEN2, SMAT, MAX_DCUT, MAX_ERR)
! split TEN into TEN1*SMAT*TEN2, where SMAT is singular value matrix
! with original legs grouped by LEGS1 and LEGS2
! the extra leg is the last leg in TEN1 and TEN2
	TYPE(TENSOR), INTENT(IN) :: TEN ! input tensor
	INTEGER, INTENT(IN) :: LEGS1(:), LEGS2(:)
	TYPE(TENSOR), INTENT(OUT) :: TEN1, TEN2, SMAT
	INTEGER, OPTIONAL, INTENT(IN) :: MAX_DCUT ! max Dcut of extra leg
	REAL, OPTIONAL, INTENT(IN) :: MAX_ERR ! max truncation error of extra leg
	! local variables
	COMPLEX, ALLOCATABLE :: A(:,:) ! A(M,N) matrix
	INTEGER :: M, N, L, I
	! SVD related
	REAL, ALLOCATABLE :: S(:), RWORK(:)
	COMPLEX, ALLOCATABLE :: U(:,:), VH(:,:), WORK(:)
	INTEGER :: LWORK, LRWORK, INFO
		
	! flatten the tensor to matrix by grouping LEGS1 and LEGS2 respectively
	A = TEN2MAT(TEN,LEGS1,LEGS2)
	! get actual size of the matrix
	M = SIZE(A,1)
	N = SIZE(A,2)
	L = MIN(M,N)
	! calculate the size of the workspace
	LWORK = MAX(2*L+MAX(M,N),MIN(L*(L+65),2*L+32*(M+N)))
	LRWORK = 5*L
	! allocate for workspace
	ALLOCATE(S(L), U(M,L), VH(L,N), WORK(LWORK), RWORK(LRWORK))
	! now everything is ready for SVD, call LAPACK
	CALL ZGESVD('S', 'S', M, N, A, M, S, U, M, VH, L, WORK, LWORK, RWORK, INFO)
	! now A = U.diag_mat(S).VH
	DEALLOCATE(A, WORK, RWORK) ! make room for TEN1 and TEN2
	IF (PRESENT(MAX_DCUT)) THEN ! if max cut specified 
		SVD_CUT = MAX_DCUT ! use it to get a suitable cut
	ELSE ! if max cut not specified
		SVD_CUT = SIZE(S) ! include all
	END IF
	IF (PRESENT(MAX_ERR)) THEN ! if max error specified 
		SVD_ERR = MAX_ERR ! use it to get a suitable Dcut
	ELSE ! if max error not specified
		SVD_ERR = 0. ! require exact
	END IF
	! set actual Dcut to L
	CALL FIND_DCUT(S, SVD_CUT, SVD_ERR)
	L = SVD_CUT
	! construct TEN1 and TEN2
	! rec tensor dimensions
	TEN1%DIMS = [TEN%DIMS(LEGS1),L]
	TEN2%DIMS = [TEN%DIMS(LEGS2),L]
	! make inds
	TEN1%INDS = [(I,I=0,M*L-1)]
	TEN2%INDS = [(I,I=0,N*L-1)]
	! fill in the vals
	TEN1%VALS = RESHAPE(U(:,1:L),[M*L])
	TEN2%VALS = RESHAPE(TRANSPOSE(VH(1:L,:)),[N*L])
	! chop off zero elements in TEN1, TEN2
	CALL TEN_CHOP(TEN1)
	CALL TEN_CHOP(TEN2)
	SMAT = DIAG_MAT(CMPLX(S(1:L)))
END SUBROUTINE SVD
! Kronecker Product SVD
SUBROUTINE KPSVD(TEN, LEGS1, LEGS2, TEN1, TEN2, MAX_DCUT, MAX_ERR)
! split TEN into the TEN1 and TEN2 (with one extra leg connecting them)
! with original legs grouped by LEGS1 and LEGS2
! the extra leg is the last leg in TEN1 and TEN2
	TYPE(TENSOR), INTENT(IN) :: TEN ! input tensor
	INTEGER, INTENT(IN) :: LEGS1(:), LEGS2(:)
	TYPE(TENSOR), INTENT(OUT) :: TEN1, TEN2
	INTEGER, OPTIONAL, INTENT(IN) :: MAX_DCUT ! max Dcut of extra leg
	REAL, OPTIONAL, INTENT(IN) :: MAX_ERR ! max truncation error of extra leg
	! local variables
	COMPLEX, ALLOCATABLE :: A(:,:) ! A(M,N) matrix
	COMPLEX :: GU, GV, G
	INTEGER :: M, N, L, I
	INTEGER, ALLOCATABLE :: INDS1(:), INDS2(:)
	! SVD related
	REAL, ALLOCATABLE :: S(:), RWORK(:)
	COMPLEX, ALLOCATABLE :: U(:,:), VT(:,:), WORK(:)
	INTEGER :: LWORK, LRWORK, INFO
		
	! flatten the tensor to matrix by grouping LEGS1 and LEGS2 respectively
	A = TEN2MAT(TEN,LEGS1,LEGS2)
	! get actual size of the matrix
	M = SIZE(A,1)
	N = SIZE(A,2)
	L = MIN(M,N)
	! calculate the size of the workspace
	LWORK = MAX(2*L+MAX(M,N),MIN(L*(L+65),2*L+32*(M+N)))
	LRWORK = 5*L
	! allocate for workspace
	ALLOCATE(S(L), U(M,L), VT(L,N), WORK(LWORK), RWORK(LRWORK))
	! now everything is ready for SVD, call LAPACK
	CALL ZGESVD('S', 'S', M, N, A, M, S, U, M, VT, L, WORK, LWORK, RWORK, INFO)
	! now A = U.diag_mat(S).VT
	DEALLOCATE(A, WORK, RWORK) ! make room for TEN1 and TEN2
	IF (PRESENT(MAX_DCUT)) THEN ! if max Dcut specified 
		SVD_CUT = MAX_DCUT ! use it to get a suitable Dcut
	ELSE ! if max Dcut not specified
		SVD_CUT = SIZE(S) ! Dcut at zeros
	END IF
	IF (PRESENT(MAX_ERR)) THEN ! if max error specified 
		SVD_ERR = MAX_ERR ! use it to get a suitable Dcut
	ELSE ! if max error not specified
		SVD_ERR = 0. ! require exact
	END IF
	! set actual Dcut to L
	CALL FIND_DCUT(S, SVD_CUT, SVD_ERR)
	L = SVD_CUT
	! construct TEN1 and TEN2
	! rec tensor dimensions
	TEN1%DIMS = [TEN%DIMS(LEGS1),L]
	TEN2%DIMS = [TEN%DIMS(LEGS2),L]
	! make inds
	TEN1%INDS = [(I,I=0,M*L-1)]
	TEN2%INDS = [(I,I=0,N*L-1)]
	! distribute singular values
	DO I = 1, L
		! prepare for phase rectification
		GU = SUM(U(:,I))
		GV = SUM(VT(I,:))
		G = SQRT(CONJG(GU)*GV) ! get phase diff between U and V
		G = REAL(G*GU)*G ! make sure Re U positive
		IF (ABS(G) > 0.) THEN ! normalize G to a phase factor
			G = G/ABS(G)
		ELSE
			G = (1.,0.)
		END IF
		U(:,I) = G*SQRT(S(I))*U(:,I)
		VT(I,:) = CONJG(G)*SQRT(S(I))*VT(I,:)
	END DO
	! fill in the vals
	TEN1%VALS = RESHAPE(U,[M*L])
	TEN2%VALS = RESHAPE(TRANSPOSE(VT),[N*L])
	! chop off zero elements in TEN1, TEN2
	CALL TEN_CHOP(TEN1)
	CALL TEN_CHOP(TEN2)
END SUBROUTINE KPSVD
! symmetric KPSVD (for transpose symmetric tensor)
SUBROUTINE SYSVD(TEN, LEGS1, LEGS2, TEN1, TEN0, MAX_DCUT, MAX_ERR)
! * caller must make sure that TEN(LEGS1,LEGS2) = TEN(LEGS2,LEGS1)
! split TEN into the TEN1 - TEN0 - TEN1^T (with one extra leg connecting them)
! with original legs grouped by LEGS1 and LEGS2
! the extra leg is the last leg in TEN1
	TYPE(TENSOR), INTENT(IN) :: TEN ! input tensor
	INTEGER, INTENT(IN) :: LEGS1(:), LEGS2(:)
	TYPE(TENSOR), INTENT(OUT) :: TEN1, TEN0
	INTEGER, OPTIONAL, INTENT(IN) :: MAX_DCUT ! max Dcut of extra leg
	REAL, OPTIONAL, INTENT(IN) :: MAX_ERR ! max truncation error of extra leg
	! local variables
	COMPLEX, ALLOCATABLE :: A(:,:) ! A(N,N) matrix, suppose to be symmetric
	COMPLEX :: GU, GV, G
	INTEGER :: M, N, L, I
	INTEGER, ALLOCATABLE :: INDS1(:), INDS2(:)
	! SVD related
	REAL, ALLOCATABLE :: S(:), RWORK(:)
	COMPLEX, ALLOCATABLE :: U(:,:), VH(:,:), WORK(:)
	INTEGER :: LWORK, LRWORK, INFO
		
	! flatten the tensor to matrix by grouping LEGS1 and LEGS2 respectively
	A = TEN2MAT(TEN,LEGS1,LEGS2)
	! get actual size of the matrix
	M = SIZE(A,1)
	N = SIZE(A,2)
	L = MIN(M,N)
	! calculate the size of the workspace
	LWORK = MAX(2*L+MAX(M,N),MIN(L*(L+65),2*L+32*(M+N)))
	LRWORK = 5*L
	! allocate for workspace
	ALLOCATE(S(L), U(M,L), VH(L,N), WORK(LWORK), RWORK(LRWORK))
	! now everything is ready for SVD, call LAPACK
	CALL ZGESVD('S', 'S', M, N, A, M, S, U, M, VH, L, WORK, LWORK, RWORK, INFO)
	! now A = U.diag_mat(S).VH
	DEALLOCATE(A, WORK, RWORK) ! make room for TEN1
	IF (PRESENT(MAX_DCUT)) THEN ! if max Dcut specified 
		SVD_CUT = MAX_DCUT ! use it to get a suitable Dcut
	ELSE ! if max Dcut not specified
		SVD_CUT = SIZE(S) ! Dcut at zeros
	END IF
	IF (PRESENT(MAX_ERR)) THEN ! if max error specified 
		SVD_ERR = MAX_ERR ! use it to get a suitable Dcut
	ELSE ! if max error not specified
		SVD_ERR = 0. ! require exact
	END IF
	! set actual Dcut to L
	CALL FIND_DCUT(S, SVD_CUT, SVD_ERR)
	L = SVD_CUT
	! construct TEN0 as diag mat of singular values
	TEN0 = DIAG_MAT(CMPLX(S(1:L)))
	! gauge fixing
	U(:,1:L) = MATMUL(U(:,1:L),INV_SQRT(MATMUL(CONJG(VH(1:L,:)),U(:,1:L))))
	! construct TEN1
	! rec tensor dimensions
	TEN1%DIMS = [TEN%DIMS(LEGS1),L]
	! make inds
	TEN1%INDS = [(I,I=0,M*L-1)]
	! fill in the vals
	TEN1%VALS = RESHAPE(U(:,1:L),[M*L])
	! chop off zero elements in TEN1
	CALL TEN_CHOP(TEN1)
END SUBROUTINE SYSVD
! High-Order SVD (in-place)
SUBROUTINE HOSVD(TEN, US, LEGS0, DCUTS)
! perform in-place HOSVD on TEN and return the unitary U's
! A_{i} = S_{j} * {U_ij} 
! input tensor A, on exit TEN is replaced by the core tensor S
! LEGS0 specifies the legs that should be diagonalized
! by default, all legs will be diagonalized
! DCUTS specifies the dimension cuts for each leg that is specified in LEGS0
! * US must be allocatable on calling
	TYPE(TENSOR), INTENT(INOUT) :: TEN
	TYPE(TENSOR), ALLOCATABLE, INTENT(OUT) :: US(:)
	INTEGER, OPTIONAL, INTENT(IN) :: LEGS0(:), DCUTS(:)
	! local variables
	INTEGER :: NLEG, ILEG, RANK, LLEG
	INTEGER, ALLOCATABLE :: LEGS(:), RLEGS(:)
	! SVD related
	REAL, ALLOCATABLE :: S(:), RWORK(:)
	COMPLEX, ALLOCATABLE :: A(:,:), U(:,:), VT(:,:), WORK(:)
	INTEGER :: LD, RD, M, N, L, LWORK, LRWORK, INFO
	
	! the legs to perform SVD
	IF (PRESENT(LEGS0)) THEN
		NLEG = SIZE(LEGS0)
		LEGS = LEGS0
	ELSE
		NLEG = SIZE(TEN%DIMS)
		LEGS = [(ILEG,ILEG=1,NLEG)]
	END IF
	ALLOCATE(US(NLEG)) ! allocate for US
	! prepare space for SDV
	LD = MAXVAL(TEN%DIMS(LEGS)) ! max leading dim 
	RD = MAX(PRODUCT(TEN%DIMS)/MINVAL(TEN%DIMS(LEGS)),LD) ! max remaining dim
	! estimate size of work space, note that RD >= LD
	LWORK = MAX(2*LD+RD,MIN(LD*(LD+65),2*LD+32*SUM(TEN%DIMS)))
	LRWORK = 5*LD
	! allocate for work space
	ALLOCATE(S(LD), U(LD,LD), VT(LD,RD), WORK(LWORK), RWORK(LRWORK))
	DO ILEG = 1, NLEG ! loop over all 
		! flatten tensor to matrix, leading leg: LEGS(ILEG)
		A = TEN2MAT(TEN,[LEGS(ILEG)])
		! get actual size of the matrix
		M = SIZE(A,1)
		N = SIZE(A,2)
		L = MIN(M,N)
		IF (PRESENT(DCUTS)) L = MIN(L,DCUTS(ILEG))
		! call LAPACK to perform SVD on A
		CALL ZGESVD('S', 'S', M, N, A, M, S, U, LD, VT, LD, WORK, LWORK, RWORK, INFO)
		! now transformation matrix is stored in U(:M,:L)
		! organize U matrix into a tensor and record
		US(ILEG) = MAT2TEN(U(:M,:L),[M,L],[1])
	END DO
	DEALLOCATE(A, S, U, VT, WORK, RWORK) ! make room for TEN construction
	! now U for all legs are obtained in US in the form of 2-tensor
	! perform Tucker product to get the core tensor
	! A_{i} * {U^*_ij} = S_{j}
	RANK = SIZE(TEN%DIMS) ! keep the tensor rank for later use
	DO ILEG = 1, NLEG
		LLEG = LEGS(ILEG) ! rec the leading leg
		! perform tensor product to apply U^H to the leading leg
		TEN = TEN_PROD(TEN,TEN_CONJG(US(ILEG)),[LLEG],[1])
		! after the product, the leading leg is arranged to the end
		! like: [RLEGs..., LLEG], which must be reordered using ORDERING
		! rectify the leg arrangement by tensor transpose
		TEN = TEN_TRANS(TEN,ORDERING([REMAINING_LEGS(RANK,[LLEG]),LLEG]))
	END DO
END SUBROUTINE HOSVD
! Auxiliary procedures ----------------------
! find suitable Dcut from singular val list
SUBROUTINE FIND_DCUT(S, CUT, ERR)
! called by KPSVD, SYSVD
! S(:) - list of singular vals, assumed to have ordered
! on input:
! CUT - upper bound of dim cut
! ERR  - upper bound of relative error
! the suitable cut must:
! 1. not exceeding max cut
! 2. S(1:CUT) contains no zeros
! 3. degenerate singular vals should not be separated
! 4. under above, relative err not exceeding max err
! under the above conditions, find cut as small as possible
! on output:
! CUT and ERR are modified to the actual value
	REAL, INTENT(IN) :: S(:)
	INTEGER, INTENT(INOUT) :: CUT
	REAL, INTENT(INOUT) :: ERR
	! local variables
	REAL :: EPS, S0, S1, SE, SA
	INTEGER :: NS, GAP
	
	IF (SIZE(S) == 0) THEN ! if S empty 
		CUT = 0  ! return with cut = 0
		ERR = 0. ! which is exact
		RETURN
	END IF
	! if S not empty
	S0 = S(1) ! get the leading S val
	IF (S0 == 0.) THEN ! if leading S val = 0
		CUT = 0  ! return with cut = 0 
		ERR = 0. ! which is exact
		RETURN
	END IF
	! use S(1) to estimate the resolution
	EPS = S0*TOL ! estimate EPS
	! singular vals within EPS will be considered the same
	NS = MIN(SIZE(S),CUT) ! set upper bound to the max cut
	! move the upper bound down to the edge of zero zone
	DO WHILE (ABS(S(NS)) < EPS) ! if NS is in the zero zone
		NS = NS - 1 ! decrease NS
	END DO
	! now NS rest in a position that is within the max cut
	! and out of the zero zone (1. 2. satisfied)
	! suitable cut must be find <= NS
	! first, find the first GAP within NS
	GAP = NS ! initialize with NS
	IF (NS < SIZE(S)) THEN ! take care of degeneracy if not at the end
		! as now NS < SIZE(S), S(NS+1) is safe to reference 
		DO WHILE (ABS(S(GAP)-S(GAP+1)) < EPS) ! if degenerate
			GAP = GAP - 1 ! decrease GAP
			IF (GAP == 0) THEN
				! this means all are degenerated up to NS
				! have to give up 3. and set GAP back to NS
				GAP = NS
				EXIT
			END IF
		END DO ! until GAP is at a gap-left
	END IF ! now NS is either at the list end or at a gap-left
	! searching for suitable cut (using 3. 4.)
	SA = SUM(S**2) ! prepare total S^2
	SE = SUM(S(GAP+1:)**2) ! initialize SE
	IF (GAP == 1) THEN ! if gap at the begining
		CUT = GAP ! must at least return the non-zero first mode
		ERR = SE/SA
		RETURN
	ELSE ! otherwise, return with the gap before exceeding max err
		! if already exceed max err
		IF (SE/SA >= ERR) THEN ! return with the last gap
			CUT = GAP
			ERR = SUM(S(GAP+1:)**2)/SA
			RETURN
		END IF
		! if within max err, can try lower cut
		DO NS = GAP, 1, -1 ! going downward
			! accumulate this S^2 to SE
			SE = SE + S(NS)**2
			! if exceed max err
			IF (SE/SA >= ERR) THEN ! return with the last gap
				CUT = GAP
				ERR = SUM(S(GAP+1:)**2)/SA
				RETURN
			END IF
			! check the next gap between:
			! ... NS-1 | NS ...
			! if not degenerated
			IF (NS > 1 .AND. ABS(S(NS-1)-S(NS)) > EPS) THEN
				! discover a new gap in the S vals
				GAP = NS-1 ! update the gap-left position
			END IF
		END DO ! I <- I-1
	END IF
	! if the subroutine has not return yet
	! it means the max err is too large that even if NS = 1
	! the max err has not been exceeded
	CUT = GAP ! in this case, set CUT to the latest gap
	ERR = SUM(S(GAP+1:)**2)/SA ! and return with actual err
END SUBROUTINE FIND_DCUT
! collect index for legs
FUNCTION COLLECT_INDS(TEN, LEGS) RESULT(LINDS)
! first get dims and dim-prods for legs
! then for each IND, the leg-collected ind is
! LIND = mod(IND/P1,D1) + D1*(mod(IND/P2,D2) + D2*(mod(IND/P3,D3) + ...))
	TYPE(TENSOR), INTENT(IN) :: TEN
	INTEGER, INTENT(IN) :: LEGS(:)
	! input indices, legs, dimensions and product of dimensions
	INTEGER, ALLOCATABLE :: LINDS(:) ! resulting indices of given legs
	! local variable
	INTEGER :: PRDS(SIZE(TEN%DIMS))
	INTEGER :: DS(SIZE(LEGS)), PS(SIZE(LEGS))
	INTEGER :: NLEG, ILEG, IBLOCK
	
	NLEG = SIZE(LEGS) ! get num of leading legs
	IF (NLEG == 0) THEN ! if no leg to collect
		ALLOCATE(LINDS(SIZE(TEN%INDS)))
		LINDS = 0 ! return a list of 0
		RETURN
	END IF ! otherwise, collect leading index
	! prepare dimension products
	! PRDS(l) = prod(DIMS(1:l-1)) 
	PRDS(1) = 1 ! prod of DIMS(null) = 1
	DO ILEG = 1, SIZE(TEN%DIMS) - 1
		PRDS(ILEG+1) = PRDS(ILEG)*TEN%DIMS(ILEG)
	END DO
	! initialize the first block
	IBLOCK = 1	
	DS(1) = TEN%DIMS(LEGS(1))
	PS(1) = PRDS(LEGS(1))
	DO ILEG = 2, NLEG ! for the rest of the legs (if any)
		! seek the possibility of blocking
		IF (LEGS(ILEG) == LEGS(ILEG-1)+1) THEN
			! if current leg is continued from the previous
			! it belongs to the same block, do dimension extension
			DS(IBLOCK) = DS(IBLOCK)*TEN%DIMS(LEGS(ILEG))
		ELSE ! if current leg is not continued from the previous
			! start a new block
			IBLOCK = IBLOCK + 1 ! block marker inc
			DS(IBLOCK) = TEN%DIMS(LEGS(ILEG))
			PS(IBLOCK) = PRDS(LEGS(ILEG))
		END IF
	END DO
	! now DS(1:IBLOCK) holds the D's for each block
	! and PS(1:IBLOCK) holds the P's for each block
	! cal LINDS by popping out block-by-block
	LINDS = MOD(TEN%INDS/PS(IBLOCK),DS(IBLOCK)) ! top block
	DO WHILE (IBLOCK > 1)
		IBLOCK = IBLOCK - 1 ! pop top block, block id dec
		LINDS(:) = MOD(TEN%INDS/PS(IBLOCK),DS(IBLOCK)) + DS(IBLOCK) * LINDS
	END DO
	! now LINDS holds the leg-collected indices to return
END FUNCTION COLLECT_INDS
! merge duplicates
SUBROUTINE MERGE_DUPLICATES(INDS, VALS, NREC)
! called by TEN_TR, TEN_PROD
! merge duplicated inds, while adding associated vals together
	INTEGER, INTENT(INOUT) :: NREC ! num of record
	! on entrance: records of all inds, on exit: the reduced records
	INTEGER, INTENT(INOUT) :: INDS(NREC) ! inds to del duplicates
	COMPLEX, INTENT(INOUT) :: VALS(NREC) ! vals to merge
	! local variables
	INTEGER :: IREC0, IREC1, IREC2
	INTEGER :: ORD(NREC)
	REAL :: EPS
	
	! get ordering of INDS by merge sort
	ORD = ORDERING(INDS)
	! set eps, val < eps will be treated as 0
	EPS = NORM(VALS)*TOL
	! for duplicated ind add associated vals together
	IREC1 = 1
	IREC0 = 0
	! use IREC1 and IREC2 to mark a sector of duplicates
	! use IREC0 to mark the position for sorted recs
	DO WHILE (IREC1 <= NREC) ! within the records
		IREC2 = IREC1 + 1 ! set sector end
		DO WHILE (IREC2 <= NREC .AND. INDS(ORD(IREC2)) == INDS(ORD(IREC1)))
			! if index duplicated, val added to IREC1 position
			VALS(ORD(IREC1)) = VALS(ORD(IREC1)) + VALS(ORD(IREC2))
			IREC2 = IREC2 + 1
		END DO
		! now IREC2 gets out of the duplicated sector
		! and the gathered val-sum is pointed by IREC1
		! if the val-sum is not 0, then move data to the front
		IF (ABS(VALS(ORD(IREC1))) > EPS) THEN
			! instead of moving both INDS and VALS
			! transfer the ORD index to the front (IREC1 -> IREC0)
			IREC0 = IREC0 + 1 ! IREC0 inc
			! because IREC0 <= IREC1 always true, data will not interfere
			ORD(IREC0) = ORD(IREC1)
		END IF
		IREC1 = IREC2 ! sector beginning moving to next
	END DO
	! now INDS(ORD(1:IREC0)) holds ordered non-duplicated index
	! and VALS(ORD(1:IREC0)) holds the corresponding values
	NREC = IREC0 ! update record num
	! scatter the data to the front
	INDS(:NREC) = INDS(ORD(:NREC))
	VALS(:NREC) = VALS(ORD(:NREC))
END SUBROUTINE MERGE_DUPLICATES
! ordering (by merge sort)
FUNCTION ORDERING(F) RESULT(X)
! return the ordering of F in X
	INTEGER, INTENT(IN) :: F(:) ! input data array
	INTEGER, TARGET  :: X(SIZE(F)) ! output ordering
	! local variables
	INTEGER, TARGET  :: Y(SIZE(F)) ! work space
	INTEGER, POINTER :: A(:), B(:)  ! pointer to ORE, ORE
 	INTEGER :: S(SIZE(F) + 1) ! sector initial marker
	INTEGER :: IS, NS
	INTEGER :: LA0, LA1, LA2, IA1, IA2, IB
	LOGICAL :: AB
	INTEGER :: F0
	
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
				IF (F(A(IA1)) < F(A(IA2))) THEN ! F(A(IA1)) is smaller
					B(IB) = A(IA1) ! smaller first
					IB = IB + 1    ! IB inc
					IA1 = IA1 + 1  ! next
				ELSEIF (F(A(IA1)) > F(A(IA2))) THEN ! F(A(IA2)) is smaller
					B(IB) = A(IA2) ! smaller first
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
END FUNCTION ORDERING
! estimate the norm of dim-1 complex array
FUNCTION NORM(ZS) RESULT(X)
	COMPLEX, INTENT(IN) :: ZS(:) ! array input
	REAL :: X ! norm output
	! local variables
	INTEGER, PARAMETER :: MAX_SIZE = 100
	INTEGER :: NZ
	REAL, ALLOCATABLE :: H(:) ! space to harvest random reals
	
	NZ = SIZE(ZS)
	IF (NZ < MAX_SIZE) THEN
		X = SQRT(SUM(ABS(ZS)**2)/NZ)
	ELSE
		ALLOCATE(H(MAX_SIZE)) ! allocate space of rand reals
		! H is of the size MAX_SIZE
		CALL RANDOM_NUMBER(H) ! 0 <= H < 1
		! to map H in [0,1) to {1..NZ}, scaled by (NZ-1) and shift by 1 
		X = SQRT(SUM(ABS(ZS(1+FLOOR(H*(NZ-1))))**2)/MAX_SIZE)
		DEALLOCATE(H) ! free H
	END IF
END FUNCTION NORM
! matrix inverse square root by Schur decomposition
FUNCTION INV_SQRT(A) RESULT(B)
! input: A must be a square matrix
! output: INV_SQRT(A) = A^(-1/2)
	COMPLEX, INTENT(IN) :: A(:,:)
	COMPLEX, ALLOCATABLE :: B(:,:)
	! local variable
	LOGICAL :: SEL
	LOGICAL :: BWORK(0)
	INTEGER :: I, N, SDIM, LDVS, LWORK, INFO
	COMPLEX, ALLOCATABLE :: W(:), VS(:,:), WORK(:)
	REAL, ALLOCATABLE :: RWORK(:)
	
	N = SIZE(A,1) ! dim of A
	B = A ! move data to B
	! allocate for space
	LDVS = N
	LWORK = 65*N
	ALLOCATE(W(N),VS(LDVS,N),WORK(LWORK),RWORK(N))
	! call LAPACK to perform Schur decomposition
	CALL ZGEES('V','N',SEL,N,B,N,SDIM,W,VS,LDVS,WORK,LWORK,RWORK,BWORK,INFO)
	! now W stores the eigen values
	! construct inverse square root matrix
	FORALL (I=1:N)
		B(I,I) = (1.,0.)/SQRT(B(I,I))
	END FORALL
	! restore the original basis
	B = MATMUL(MATMUL(VS,B),CONJG(TRANSPOSE(VS)))
END FUNCTION INV_SQRT
! I/O System -------------------------------
! tensor print
SUBROUTINE TEN_PRINT_0(TEN)
	TYPE(TENSOR), INTENT(IN) :: TEN ! input tensor to print out
	! local variable
	INTEGER :: L, NDIM, IDIM, IND
	
	NDIM = SIZE(TEN%DIMS)
	IF (NDIM == 0) THEN ! null tensor
		WRITE(*,'("TENSOR null")')
		DO L = 1, SIZE(TEN%VALS)
			WRITE(*,'("("I3") -> "2F7.3)') TEN%INDS(L), TEN%VALS(L)
		END DO
	ELSE ! normal tensor
		WRITE(*,'("TENSOR")',ADVANCE='NO') ! print 'TENSOR'
		DO IDIM = 1, NDIM-1
			WRITE (*,'(I3" x")',ADVANCE='NO') TEN%DIMS(IDIM)
		END DO
		WRITE (*,'(I3)') TEN%DIMS(NDIM)
		! for each record of non-zeros (labeled by L)
		DO L = 1, SIZE(TEN%VALS)
			IND = TEN%INDS(L) ! get IND
			! decompose IND into each dimension and print
			WRITE(*,'("(")',ADVANCE='NO') ! print a '(' first
			DO IDIM = 1, NDIM
				! cal the index by MOD+1, and print
				WRITE (*,'(I3)',ADVANCE='NO') MOD(IND,TEN%DIMS(IDIM)) + 1
				! print ',' between
				IF (IDIM < NDIM) WRITE(*,'(",")',ADVANCE='NO')
				! next dimension by dividing IND
				IND = IND/TEN%DIMS(IDIM)
			END DO
			WRITE(*,'(") -> "2F7.3)') TEN%VALS(L) ! end by printing the value
		END DO
	END IF
END SUBROUTINE TEN_PRINT_0
! tensor print (array)
SUBROUTINE TEN_PRINT_1(TENS)
	TYPE(TENSOR), INTENT(IN) :: TENS(:) ! input tensor array
	INTEGER :: NTEN, ITEN
	
	NTEN = SIZE(TENS) ! get size of array
	DO ITEN = 1, NTEN ! print each tensor in the array
		CALL TEN_PRINT_0(TENS(ITEN))
	END DO
END SUBROUTINE TEN_PRINT_1
! save tensor to disk
SUBROUTINE TEN_SAVE_0(FILENAME, TEN)
! save TEN to the file named FILENAME
! can be loaded by Mathematica by TensorLoad
	CHARACTER(*), INTENT(IN) :: FILENAME
	TYPE(TENSOR), INTENT(IN) :: TEN
	INTEGER :: NTEN, NDIM, NREC
	
	PRINT *, '>>', FILENAME, '.ten'
    OPEN (UNIT = 99, FILE = FILENAME//'.ten', &
        STATUS = 'REPLACE', ACCESS = 'STREAM')
    NTEN = 1 ! get num of tensors
    NDIM = SIZE(TEN%DIMS) ! get num of legs
    NREC = SIZE(TEN%INDS) ! get num of recs
    WRITE(99) NTEN ! put NTEN
    WRITE(99) NDIM ! put NDIM
    WRITE(99) NREC ! put NREC
    ! dump data from TEN
    WRITE(99) TEN%DIMS
    WRITE(99) TEN%INDS
    WRITE(99) TEN%VALS
    CLOSE(99) ! close stream
END SUBROUTINE TEN_SAVE_0
! load tensor from disk
SUBROUTINE TEN_LOAD_0(FILENAME, TEN)
! save TEN to the file named FILENAME
! can be loaded by Mathematica by TensorLoad
	CHARACTER(*), INTENT(IN) :: FILENAME
	TYPE(TENSOR), INTENT(OUT) :: TEN
	INTEGER :: NTEN, NDIM, NREC
	
	PRINT *, '<<', FILENAME, '.ten'
    OPEN (UNIT = 99, FILE = FILENAME//'.ten', &
        STATUS = 'UNKNOWN', ACCESS = 'STREAM')
    READ(99) NTEN ! get num of tensors
    IF (NTEN /= 1) THEN ! check consistency (to prevent reading old tensor file)
    	WRITE (*,'(A)') 'TEN_LOAD::nten: expect only one tensor in the data file.'
    	STOP 
    END IF
    READ(99) NDIM, NREC ! get num of legs and recs
    ! declare TEN by allocation
	ALLOCATE(TEN%DIMS(NDIM), TEN%INDS(NREC), TEN%VALS(NREC))
	! load data into TEN
    READ(99) TEN%DIMS
    READ(99) TEN%INDS
    READ(99) TEN%VALS
    CLOSE(99) ! close stream
END SUBROUTINE TEN_LOAD_0
! save tensor to disk (array)
SUBROUTINE TEN_SAVE_1(FILENAME, TENS)
! save TEN to the file named FILENAME
! can be loaded by Mathematica by TensorLoad
	CHARACTER(*), INTENT(IN) :: FILENAME
	TYPE(TENSOR), INTENT(IN) :: TENS(:)
	INTEGER :: NTEN, ITEN, NDIM, NREC
	
	PRINT *, '>>', FILENAME, '.ten'
    OPEN (UNIT = 99, FILE = FILENAME//'.ten', &
        STATUS = 'REPLACE', ACCESS = 'STREAM')
    NTEN = SIZE(TENS) ! get num of tens
    WRITE(99) NTEN ! put NTEN
    DO ITEN = 1, NTEN ! for each tensor in the array
		NDIM = SIZE(TENS(ITEN)%DIMS) ! get num of legs
		NREC = SIZE(TENS(ITEN)%INDS) ! get num of recs
		WRITE(99) NDIM ! put NDIM
		WRITE(99) NREC ! put NREC
		! dump data from TEN
		WRITE(99) TENS(ITEN)%DIMS
		WRITE(99) TENS(ITEN)%INDS
		WRITE(99) TENS(ITEN)%VALS
	END DO
    CLOSE(99) ! close stream
END SUBROUTINE TEN_SAVE_1
! load tensor from disk (array)
SUBROUTINE TEN_LOAD_1(FILENAME, TENS)
! save TEN to the file named FILENAME
! can be loaded by Mathematica by TensorLoad
	CHARACTER(*), INTENT(IN) :: FILENAME
	TYPE(TENSOR), ALLOCATABLE, INTENT(OUT) :: TENS(:)
	INTEGER :: NTEN, ITEN, NDIM, NREC
	
	PRINT *, '<<', FILENAME, '.ten'
    OPEN (UNIT = 99, FILE = FILENAME//'.ten', &
        STATUS = 'UNKNOWN', ACCESS = 'STREAM')
    READ(99) NTEN ! get num of tens
    ALLOCATE(TENS(NTEN)) ! allocate tensor array
    DO ITEN = 1, NTEN ! for each tensor
		READ(99) NDIM, NREC ! get num of legs and recs
		! declare TEN by allocation
		ALLOCATE(TENS(ITEN)%DIMS(NDIM), TENS(ITEN)%INDS(NREC), TENS(ITEN)%VALS(NREC))
		! load data into TEN
		READ(99) TENS(ITEN)%DIMS
		READ(99) TENS(ITEN)%INDS
		READ(99) TENS(ITEN)%VALS
    END DO
    CLOSE(99) ! close stream
END SUBROUTINE TEN_LOAD_1
! Accessibility check ----------------------
SUBROUTINE TEN_CHECK(CALLER, TEN)
	CHARACTER(*), INTENT(IN) :: CALLER
	TYPE(TENSOR), INTENT(IN) :: TEN
	
	IF (.NOT. (ALLOCATED(TEN%DIMS) .AND. ALLOCATED(TEN%INDS) .AND. ALLOCATED(TEN%VALS))) THEN
		WRITE (*,'(A)') CALLER//'::xten: input tensor not defined.'
		STOP
	END IF
	IF (SIZE(TEN%INDS) /= SIZE(TEN%VALS)) THEN
		WRITE (*,'(A)') CALLER//'::xlen: input tensor INDS and VALS have different sizes.'
		STOP
	END IF
END SUBROUTINE TEN_CHECK
! end of the module
END MODULE TENSORIAL