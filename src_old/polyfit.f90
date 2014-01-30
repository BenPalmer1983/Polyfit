Program polyfitprog

! University of Birmingham
! Ben Palmer
!


!Setup Modules
Use kinds				!data kinds
Use initialise			!initialise program
Use maths				!maths functions


!force declaration of all variables
  Implicit None

!declare variables
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPoints
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsA
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsB
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: dataPointsX
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: dataPointsY
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsA
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsB
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: dataPointsS
  Real(kind=SingleReal), Dimension( : ), Allocatable :: coefficientsS
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: testSquare
  Integer(kind=StandardInteger) :: i,j,k,polyOrder
  Real(kind=DoubleReal) :: rss
  
!initialise
  Call runInitialise()
  
  polyOrder=5
  Allocate(dataPointsX(1:21))
  Allocate(dataPointsY(1:21))
  Allocate(coefficients(0:polyOrder))
  Allocate(dataPoints(1:21,1:2))
  Allocate(dataPointsA(1:6,1:2))
  
  dataPointsX(1) = 7970.0132
  dataPointsY(1) = -1680.3405
  dataPointsX(2) = 7994.188
  dataPointsY(2) = -1680.6499
  dataPointsX(3) = 8018.4116
  dataPointsY(3) = -1680.9008
  dataPointsX(4) = 8042.6865
  dataPointsY(4) = -1681.0769
  dataPointsX(5) = 8067.0107
  dataPointsY(5) = -1681.01
  dataPointsX(6) = 8091.3813
  dataPointsY(6) = -1680.9685
  dataPointsX(7) = 8115.8008
  dataPointsY(7) = -1680.9734
  dataPointsX(8) = 8140.272
  dataPointsY(8) = -1681.0387
  dataPointsX(9) = 8164.7896
  dataPointsY(9) = -1680.9819
  dataPointsX(10) = 8189.3589
  dataPointsY(10) = -1680.8303
  dataPointsX(11) = 8213.9756
  dataPointsY(11) = -1680.7107
  dataPointsX(12) = 8238.6436
  dataPointsY(12) = -1680.7112
  dataPointsX(13) = 8263.3584
  dataPointsY(13) = -1680.4535
  dataPointsX(14) = 8288.1221
  dataPointsY(14) = -1680.1946
  dataPointsX(15) = 8312.9365
  dataPointsY(15) = -1679.8531
  dataPointsX(16) = 8337.8018
  dataPointsY(16) = -1679.5502
  dataPointsX(17) = 8362.7168
  dataPointsY(17) = -1679.3109
  dataPointsX(18) = 8387.6797
  dataPointsY(18) = -1678.856
  dataPointsX(19) = 8412.6934
  dataPointsY(19) = -1678.2948
  dataPointsX(20) = 8437.7529
  dataPointsY(20) = -1678.0349
  dataPointsX(21) = 8462.8672
  dataPointsY(21) = -1677.4668
  
  dataPoints(1,1) = 7970.0132
  dataPoints(1,2) = -1680.3405
  dataPoints(2,1) = 7994.188
  dataPoints(2,2) = -1680.6499
  dataPoints(3,1) = 8018.4116
  dataPoints(3,2) = -1680.9008
  dataPoints(4,1) = 8042.6865
  dataPoints(4,2) = -1681.0769
  dataPoints(5,1) = 8067.0107
  dataPoints(5,2) = -1681.01
  dataPoints(6,1) = 8091.3813
  dataPoints(6,2) = -1680.9685
  dataPoints(7,1) = 8115.8008
  dataPoints(7,2) = -1680.9734
  dataPoints(8,1) = 8140.272
  dataPoints(8,2) = -1681.0387
  dataPoints(9,1) = 8164.7896
  dataPoints(9,2) = -1680.9819
  dataPoints(10,1) = 8189.3589
  dataPoints(10,2) = -1680.8303
  dataPoints(11,1) = 8213.9756
  dataPoints(11,2) = -1680.7107
  dataPoints(12,1) = 8238.6436
  dataPoints(12,2) = -1680.7112
  dataPoints(13,1) = 8263.3584
  dataPoints(13,2) = -1680.4535
  dataPoints(14,1) = 8288.1221
  dataPoints(14,2) = -1680.1946
  dataPoints(15,1) = 8312.9365
  dataPoints(15,2) = -1679.8531
  dataPoints(16,1) = 8337.8018
  dataPoints(16,2) = -1679.5502
  dataPoints(17,1) = 8362.7168
  dataPoints(17,2) = -1679.3109
  dataPoints(18,1) = 8387.6797
  dataPoints(18,2) = -1678.856
  dataPoints(19,1) = 8412.6934
  dataPoints(19,2) = -1678.2948
  dataPoints(20,1) = 8437.7529
  dataPoints(20,2) = -1678.0349
  dataPoints(21,1) = 8462.8672
  dataPoints(21,2) = -1677.4668
  
  
  dataPointsA(1,1) = 7970.0132
  dataPointsA(1,2) = -1680.3405
  dataPointsA(2,1) = 8115.8008
  dataPointsA(2,2) = -1680.9734
  dataPointsA(3,1) = 8213.9756
  dataPointsA(3,2) = -1680.7107
  dataPointsA(4,1) = 8288.1221
  dataPointsA(4,2) = -1680.1946
  dataPointsA(5,1) = 8387.6797
  dataPointsA(5,2) = -1678.856
  dataPointsA(6,1) = 8462.8672
  dataPointsA(6,2) = -1677.4668
  
  !coefficients = polyfita(dataPointsX,dataPointsY,polyOrder)
  !rss = calcResidualSquareSum(dataPoints,coefficients)
	
  !print *,""
  !print *,""
  !do i=0,(size(coefficients)-1)
  !  print *,i,coefficients(i)
  !enddo
  
  !print *,rss
  !print *,""
  !print *,""
  
  
  !coefficients = polynomialFit(dataPointsA,polyOrder)
  !rss = calcResidualSquareSum(dataPointsA,coefficients)
	
  !print *,""
  !print *,""
  !do i=0,(size(coefficients)-1)
   ! print *,i,coefficients(i)
  !enddo
  
  !print *,rss
  !print *,""
  !print *,""
  
  
  Allocate(testSquare(1:3,1:3))
  Allocate(coefficientsA(0:5))
  Allocate(coefficientsB(0:8))
  
  testSquare(1,1) = 4.0D0
  testSquare(1,2) = 7.0D0
  testSquare(1,3) = 3.5D0
  testSquare(2,1) = 1.0D0
  testSquare(2,2) = 17.0D0
  testSquare(2,3) = 5.0D0
  testSquare(3,1) = 8.3D0
  testSquare(3,2) = 3.1D0
  testSquare(3,3) = 5.7D0
  
  !coefficientsA = PolynomialRegression(dataPoints,8)
  
  !print *,coefficientsA(0),coefficientsA(1),coefficientsA(2),&
  !coefficientsA(3),coefficientsA(4),coefficientsA(5)
  
  print *,"Time: ",ProgramTime ()
  coefficientsB = PolyFit(dataPoints,5)
  print *,"Time: ",ProgramTime ()
  
  print *,coefficientsB(0),coefficientsB(1),coefficientsB(2),&
  coefficientsB(3),coefficientsB(4),coefficientsB(5),coefficientsB(6),&
  coefficientsB(7),coefficientsB(8)
  
  
  print *,"Rss: ",CalcResidualSquareSum(dataPoints,coefficientsB)	
  
    
  !coefficientsA = PolyFit(dataPoints,5)
  
  !print *,coefficientsA(0),coefficientsA(1),coefficientsA(2),&
  !coefficientsA(3),coefficientsA(4),coefficientsA(5)
	

End