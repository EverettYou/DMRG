! square lattice MPO
SUBROUTINE SET_MPO()
! data transfered by DATAPOOL
! DATAPOOL: TA, TB. MODEL: THETA, BETA, CROSS
! * DATAPOOL must be initiated before calling me
	USE DATAPOOL
	! local variables
	TYPE(TENSOR) :: X, Y, U, UA, UB, S
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
	! set TA from TB MPO
	TA = TEN_TRANS(TB,[2,1,3,4]) ! given by 1 <-> 2
END SUBROUTINE SET_MPO
