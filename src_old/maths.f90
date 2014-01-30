Module maths

!----------------------------------------!
! Maths functions                        !
! Ben Palmer, University of Birmingham   !
!----------------------------------------!


! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
!Make private
  Private  
  
! Polynomial Related Functions
  Public :: SolvePolynomial  
  Public :: CalcPolynomial  
! Interpolation and Regression Functions 
  Public :: PolyFit  
  Public :: PolynomialInterpolation  
  Public :: QuadraticInterpolation  
  Public :: QuadraticInterpolationCalc  
  Public :: CubicInterpolationCalc  
  Public :: CalcResidualSquareSum  
  
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 
! -
! 
! Polynomial Related Functions
! - SolvePolynomial
! - CalcPolynomial
!  
!  
! Interpolation and Regression Functions 
! - PolyFit
! - PolynomialInterpolation 
! - PolynomialRegression
! - PolynomialRegressionVandermonde
! - QuadraticInterpolation
! - QuadraticInterpolationCalc
! - CubicInterpolationCalc
! - CalcResidualSquareSum
!
! Matrix Functions
! - InvertSquareMatrix
! 
! Random Number Related Functions
! - RandomSeed
! - RandomDataPoints
! 
  

  
  
!------------------------------------------------------------------------!
! Polynomial Related Functions
!------------------------------------------------------------------------!     
  
  Function SolvePolynomial (coefficients, lower, upper) RESULT (output)    	
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: upper, lower, output
	Real(kind=DoubleReal) :: x,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget, factor, difference
	Integer(kind=StandardInteger) i,j,k,maxLoops	
!Set values
	convergenceTarget = 0
	convergence = 1000
	convergenceThreshold = 0.00001
	maxLoops = 0
!set start value for x
	factor = 0.5
	difference = upper - lower
	x = lower + factor * difference	
	do while(convergence.gt.convergenceThreshold.and.maxLoops.le.10000)
	  maxLoops = maxLoops + 1
	  difference = factor * difference
	  y = 0
	  do i=0,size(coefficients)-1
	    y = y + x**(i) * coefficients(i)			
	  enddo
	  dydx = 0
	  do i=1,size(coefficients)-1
	    dydx = dydx + i * x**(i-1) * coefficients(i)			
	  enddo		  
	  convergence = abs(convergenceTarget - y)
	  if(convergence.gt.convergenceThreshold)then
	    if((dydx.lt.0.and.y.ge.0).or.(dydx.ge.0.and.y.lt.0))then
	      x = x + difference	
	    else
	      x = x - difference
	    endif
	  endif
	enddo	
	output = x
  End Function solvePolynomial   
  
   
  Function CalcPolynomial (polyCoefficients, x, derivative) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k, derivative
	Integer(kind=StandardInteger) :: coeffCount 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyWorking
    Real(kind=DoubleReal) :: x, y
!count coefficients
	coeffCount = size(polyCoefficients)	
!force no or positive derivative
	if(derivative.lt.0)then
	  derivative = 0
	endif
!run calc
	if(derivative.gt.coeffCount)then
!if order of derivative higher than coeffs/order polynomial, then answer is zero	
	  y = 0.0D0
	  Allocate(polyWorking(0:0))
!transfer to polyWorking
	  polyWorking(0) = 0.0D0
	  coeffCount = 1
	else
	  if(derivative.eq.0)then
!if no derivative, use input coeffs	
        Allocate(polyWorking(0:(coeffCount-1)))
!transfer to polyWorking
	    do i=0,(coeffCount-1)
		  polyWorking(i) = polyCoefficients(i)
		enddo
	  else
!loop through each derivative
	    do i=1,derivative
		  coeffCount = coeffCount - 1
		  do j=0,(coeffCount-1)
		    polyCoefficients(j) = (j+1)*polyCoefficients(j+1)
		  enddo
		enddo	
        Allocate(polyWorking(0:(coeffCount-1)))
!transfer to polyWorking
	    do i=0,(coeffCount-1)
		  polyWorking(i) = polyCoefficients(i)
		enddo
	  endif
	  y = 0.0D0
	  do i=0,(coeffCount-1)
	    y = y+1.0D0*polyWorking(i)*x**(1.0D0*i)
	  enddo 
	endif 
!returns y
  End Function CalcPolynomial
  
  
  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
! Interpolation and Regression Functions 
!------------------------------------------------------------------------! 
  
  Function PolyFit(inputPoints, order) RESULT (coefficients) 
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,order
	Real(kind=DoubleReal) :: rssThreshold, rssA, rssB, rssC, rssBest
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsA
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsB
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsTemp
    Integer(kind=StandardInteger), Dimension(1:1000) :: randomSeed
!set random seed
    randomSeed = SetRandomSeedArray()
!allocate arrays
    Allocate(coefficients(0:order))
    Allocate(coefficientsA(0:order))
    Allocate(coefficientsB(0:order))
    Allocate(coefficientsTemp(0:order))
!Use vandermonde method and calculate rss
	coefficientsA = PolynomialRegressionVandermonde(inputPoints, order)
	rssA = CalcResidualSquareSum(inputPoints,coefficientsA)
!calculate and compare rss for each
    rssThreshold = rssA
	coefficientsB = PolynomialRegression(inputPoints, order, rssThreshold)
	rssB = CalcResidualSquareSum(inputPoints,coefficientsB)
	if(rssA.lt.rssB)then
!use vandermonde result
      do i=0,order
	    coefficients(i) = coefficientsA(i) 
	  enddo
    else
!improve regression attempt if better than vandermonde method	  
	  rssBest = rssB
	  print *,rssB
	  do i=1,4
	    rssThreshold = rssThreshold*0.2**i
		coefficientsTemp = PolynomialRegression(inputPoints, order, rssThreshold)
		rssB = CalcResidualSquareSum(inputPoints,coefficientsTemp)
		if(rssB.lt.rssBest)then
		  do j=0,order
		    coefficientsB(j) = coefficientsTemp(j)
			rssBest = rssB
		  enddo	
		else
          exit
        endif		  
	  enddo
	  do i=0,order
	    coefficients(i) = coefficientsB(i) 
	  enddo
	endif
  End Function PolyFit      
  
  
  
  Function PolynomialInterpolation(inputPoints) RESULT (coefficients)   
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i, col, row
	Integer(kind=StandardInteger) :: order, coefficientCount
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsTemp
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix
!set variables
    coefficientCount = size(inputPoints,1)
	order = size(inputPoints,1)-1
!Allocate arrays
    Allocate(xMatrix(1:coefficientCount,1:coefficientCount))
    Allocate(xMatrixInverse(1:coefficientCount,1:coefficientCount))
    Allocate(yMatrix(1:coefficientCount))
    Allocate(coefficientsTemp(1:coefficientCount))
    Allocate(coefficients(0:order))
!Fill arrays
	do row=1,coefficientCount
	  do col=1,coefficientCount
        xMatrix(row,col) = 1.0D0*inputPoints(row,1)**(1.0D0*(col-1))
      enddo
	  yMatrix(row) = 1.0D0*inputPoints(row,2)
    enddo
!invert xMatrix
    xMatrixInverse = InvertSquareMatrix(xMatrix)
!find coefficients
    coefficientsTemp = matMul(xMatrixInverse,yMatrix)
!move coefficients
    do i=1,coefficientCount
      coefficients(i-1) = coefficientsTemp(i)
	enddo
  End Function PolynomialInterpolation  
    
	
	
  Function PolynomialRegression(inputPoints, order, rssThresholdIn) RESULT &
  (coefficients)   
!force declaration of all variables
	Implicit None
!declare variables
    Real(kind=DoubleReal),optional :: rssThresholdIn
    Integer(kind=StandardInteger) :: i,j,k,col,row
	Integer(kind=StandardInteger) :: order, coefficientCount, polyCount, loops
	Real(kind=DoubleReal) :: randNumber, rss, bestRss, decayFactor, rssThreshold
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: samplePoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: tempCoefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: testCoefficients
!set variables	
	coefficientCount = order + 1
	if(present(rssThresholdIn))then
	  rssThreshold = rssThresholdIn
	else
	  rssThreshold = 1000
	endif
	if(rssThreshold.lt.1)then
	  !rssThreshold = 1
	endif
!allocate arrays
    Allocate(coefficients(0:order))
    Allocate(tempCoefficients(0:order))
    Allocate(testCoefficients(0:order))
    Allocate(samplePoints(1:coefficientCount,1:2))
!fill coefficients array with zeros
    do j=0,order
	  testCoefficients(j) = 0.0D0
	  coefficients(j) = 0.0D0
	  tempCoefficients(j) = 0.0D0
	enddo
!find best fit by interpolation between random data points
    polyCount = 0
	if(order.lt.3)then
	  loops = 25
	else  
	  loops = 50*((order-2)**2)
	endif
	do i=1,(48+order**3)
	  samplePoints = RandomDataPoints(inputPoints,coefficientCount,1)
	  tempCoefficients = PolynomialInterpolation(samplePoints)
	  rss = CalcResidualSquareSum(inputPoints,tempCoefficients)
	  if(rss.lt.rssThreshold)then 
	    if(polyCount.eq.0)then
	      do j=0,order
		    coefficients(j) = tempCoefficients(j)
		  enddo
		  bestRss = rss		
		  polyCount = polyCount + 1  
		else 
	      do j=0,order		
		    testCoefficients(j) = &
			   (polyCount * coefficients(j) + tempCoefficients(j)) /&
			   (polyCount + 1)
		  enddo
		  rss = CalcResidualSquareSum(inputPoints,testCoefficients)
		  if(rss.lt.bestRss)then
		    bestRss = rss	
		    polyCount = polyCount + 1
	        do j=0,order
		      coefficients(j) = 1.0D0*testCoefficients(j)
		    enddo
		  endif		  		  
		endif	  
	  endif	
	enddo
  End Function PolynomialRegression  
  

  
  Function PolynomialRegressionVandermonde(inputPoints, order)&
  RESULT (coefficients) 
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,k,col,row
	Integer(kind=StandardInteger) :: order, dataPoints, matrixSize, exponentValue
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
!set variables
    dataPoints = size(inputPoints,1)
	matrixSize = order+1
!allocate arrays
    Allocate(xMatrix(1:matrixSize,1:matrixSize))
    Allocate(yMatrix(1:matrixSize))
    Allocate(coefficients(0:order))
!Build Least Squares Fitting Vandermonde matrix
    do row=1,matrixSize
	  do col=1,matrixSize
	    exponentValue = 1.0D0*row+1.0D0*col-2.0D0
		xMatrix(row,col) = 0.0D0
		do k=1,dataPoints
	      xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*inputPoints(k,1)&
		  **exponentValue
		enddo
	  enddo
	enddo
	do row=1,matrixSize
	  exponentValue = 1.0D0*row-1.0D0
	  yMatrix(row) = 0.0D0
	  do k=1,dataPoints
	    yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*inputPoints(k,2)*&
		inputPoints(k,1)**exponentValue
	  enddo
	enddo
!invert xMatrix
	xMatrix = InvertSquareMatrix(xMatrix)
!multiply inverse by y to get coefficients
    yMatrix = matMul(xMatrix,yMatrix)
!save coefficients
    do i=0,order
	  coefficients(i) = yMatrix(i+1)
	enddo  
  End Function PolynomialRegressionVandermonde  
  
  
  
  Function QuadraticInterpolation(xA,yA,xB,yB,xC,yC) RESULT (coefficients)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC
	Real(kind=DoubleReal) :: aA,aB,aC
	Allocate(coefficients(0:2))  	
	aA = 1.0D0*(yA/((xA-xB)*(xA-xC)))
    aB = 1.0D0*(yB/((xB-xA)*(xB-xC)))
	aC = 1.0D0*(yC/((xC-xA)*(xC-xB)))
	coefficients(2) = 1.0D0*(aA+aB+aC)
	coefficients(1) = -1.0D0*(aA*(xB+xC)+aB*(xA+xC)+aC*(xA+xB))
	coefficients(0) = 1.0D0*(aA*xB*xC+aB*xA*xC+aC*xA*xB)
  End Function QuadraticInterpolation
  !------------------------------
  Function QuadraticInterpolationCalc(xA,yA,xB,yB,xC,yC,x) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC,x,y
	Allocate(coefficients(0:2)) 
	coefficients = QuadraticInterpolation(xA,yA,xB,yB,xC,yC)
	y = 1.0D0*coefficients(0)+1.0D0*coefficients(1)*x+1.0D0*coefficients(2)*x**2
  End Function QuadraticInterpolationCalc
  
  
  
  Function CubicInterpolationCalc(xA,yA,xB,yB,xC,yC,xD,yD,x) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC,xD,yD,x,y
	Real(kind=DoubleReal) :: aA,aB,aC,aD
	
	aA = 1.0D0*yA*(((x-xB)*(x-xC)*(x-xD))/((xA-xB)*(xA-xC)*(xA-xD)))
	aB = 1.0D0*yB*(((x-xA)*(x-xC)*(x-xD))/((xB-xA)*(xB-xC)*(xB-xD)))
	aC = 1.0D0*yC*(((x-xA)*(x-xB)*(x-xD))/((xC-xA)*(xC-xB)*(xC-xD)))
	aD = 1.0D0*yD*(((x-xA)*(x-xB)*(x-xC))/((xD-xA)*(xD-xB)*(xD-xC)))
		
	y = 1.0D0*(aA+aB+aC+aD)
  End Function CubicInterpolationCalc
  
  
  
  Function CalcResidualSquareSum(inputPoints,polyCoefficients) RESULT (rss) 
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: order
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
	Real(kind=DoubleReal) :: rss, x, y	
	order = size(polyCoefficients)-1	
	rss = 0.0D0
	do i=1,size(inputPoints,1)
	  x = 1.0D0*inputPoints(i,1)
	  y = calcPolynomial(polyCoefficients,x,0)
	  rss = rss + (y-inputPoints(i,2))**2  
	enddo  
  End Function CalcResidualSquareSum   
  
  
  
  
  
  
!------------------------------------------------------------------------!
! Matrix Functions
!------------------------------------------------------------------------! 
  
  Function InvertSquareMatrix(xMatrix) RESULT (xMatrixInverse)  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,row,col,rowb
	Integer(kind=StandardInteger) :: matrixSize,optimiseSum,optimiseExponent
	Real(kind=DoubleReal) :: xA, xB
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixWorking
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: xMatrixRow
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: optimiseArray	
	Logical :: optimising	
!if a square matrix
	if(size(xMatrix,1).eq.size(xMatrix,2))then
	  matrixSize = size(xMatrix,1)
!Allocate arrays
	  Allocate(xMatrixInverse(1:matrixSize,1:matrixSize))
	  Allocate(xMatrixWorking(1:matrixSize,1:2*matrixSize))	 
	  Allocate(xMatrixRow(1:2*matrixSize)) 
	  Allocate(optimiseArray(1:matrixSize)) 
!Fill working array
      do row=1,matrixSize
	    do col=1,(2*matrixSize)
	      xMatrixWorking(row,col) = 0.0D0
	    enddo
	  enddo
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
	    enddo
	  enddo
!Fill inverse array
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixInverse(row,col) = 0.0D0
	    enddo
	  enddo
!optimise matrix order
      do row=1,matrixSize		!row at a time
	    optimiseSum = 0
	    do col=1,matrixSize 
	      if(xMatrixWorking(row,col).ne.0.0D0)then
	        optimiseExponent = matrixSize - col
		    optimiseSum = optimiseSum + 2**optimiseExponent
		  endif
	    enddo
        optimiseArray(row) = optimiseSum 
	  enddo	
	  optimising = .true.
      do while(optimising)
        optimising = .false.
        do row=1,(matrixSize-1)        !loop through rows
          if(optimiseArray(row).lt.optimiseArray(row+1))then
		    optimising = .true.
!reorder optimising array
            i = optimiseArray(row)
            j = optimiseArray(row+1)
		    optimiseArray(row) = j
		    optimiseArray(row+1) = i
!reorder xMatrixWorking
            do col=1,(2*matrixSize)   
		      xA = 1.0D0*xMatrixWorking(row,col)
              xB = 1.0D0*xMatrixWorking(row+1,col)
		      xMatrixWorking(row,col) = 1.0D0*xB
		      xMatrixWorking(row+1,col) = 1.0D0*xA
		    enddo
		  endif
	    enddo
      enddo	
!Make identity in rhs matrix
      do row=1,matrixSize
	    do col=1,matrixSize
		  if(row.eq.col)then
	        xMatrixWorking(row,col+matrixSize) = 1.0D0
		  endif
	    enddo
	  enddo
!make lower triangle of zeros	  
	  do row=1,matrixSize-1
	    do rowb=row+1,matrixSize	  
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixRow(col) = 1.0D0*&
		    ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
			xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
		  enddo
!replace row values
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
		  enddo
	    enddo
!force zeros in the lower triangle
        do rowb=row+1,matrixSize
	      xMatrixWorking(rowb,row) = 0.0D0
	    enddo
	  enddo
!re-force zeros in the lower triangle
      do row=1,matrixSize
	    do col=1,matrixSize
		  if(row.gt.col)then
		    xMatrixWorking(row,col) = 0.0D0
		  endif
		enddo
	  enddo
!make upper triangle of zeros	
	  do row=matrixSize,2,-1
	    do rowb=row-1,1,-1	  
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixRow(col) = 1.0D0*&
		    ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
			xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
		  enddo
!replace row values
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
		  enddo
	    enddo
!force zeros in the upper triangle
        do rowb=row-1,1,-1
	      xMatrixWorking(rowb,row) = 0.0D0
	    enddo
	  enddo
!Divide rhs by diagonal on lhs and store in inverse
	  do row=1,matrixSize
	    do col=1,matrixSize
		  xMatrixInverse(row,col) = 1.0D0*&
		  xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
		enddo
	  enddo
	endif
  End Function InvertSquareMatrix  
  
  
  
  
  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
! Random Number Related Functions
!------------------------------------------------------------------------!    
  
  Function SetRandomSeedArray() RESULT (randomSeed)  
    Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: clockReturn, gap, seedTemp, multiple
    Integer(kind=StandardInteger), Dimension(1:1000) :: randomSeed
	Integer(kind=StandardInteger), Dimension(0:9) :: randomIntegers
  !random number seed from cpu
    randomIntegers(0) = 25441
    randomIntegers(1) = 37261
    randomIntegers(2) = 1622261
    randomIntegers(3) = 162982	
    randomIntegers(4) = 72635	
    randomIntegers(5) = 9927151
    randomIntegers(6) = 91
    randomIntegers(7) = 6452
    randomIntegers(8) = 448327
    randomIntegers(9) = 9253411
	j = 0
	Call SYSTEM_CLOCK(clockReturn)
	gap = clockReturn - 10*floor(1.0D0*clockReturn/10)
    do i=1,1000	  
	  k = j + gap
	  if(k.gt.9)then
	    k = k - 10
	  endif
	  seedTemp = i * (randomIntegers(j)-randomIntegers(k))
	  seedTemp = abs(seedTemp)
	  multiple = floor(1.0D0 * seedTemp / 1.0D0 * clockReturn)
	  seedTemp = seedTemp - multiple * clockReturn
	  randomSeed(i) = abs(clockReturn - seedTemp)
	  if(j.eq.9)then
	    j = 0
	  endif
	  j = j + 1
	enddo	
	Call RANDOM_SEED(put=randomSeed)
  End Function SetRandomSeedArray
  
  
  
  Function RandomDataPoints(inputPoints,noPoints, tether) RESULT (outputPoints)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,noPoints,point,noInputPoints,tether
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: selectedPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputPoints
	Logical :: pointOk
!set variables	
	noInputPoints = size(inputPoints,1) 
!Allocate arrays
    Allocate(selectedPoints(1:noPoints))
    Allocate(outputPoints(1:noPoints,1:size(inputPoints,2)))
!fill selected points array with 0s
	do i=1,noPoints
	  selectedPoints(i) = 0  
	enddo
	if(tether.eq.1)then
	  !tether start and end points
	  selectedPoints(1) = 1
	  selectedPoints(noPoints) = noInputPoints
	endif
	if(tether.eq.2)then
	  !tether start and end points
	  Call RANDOM_NUMBER(randNumber)
	  if(randNumber.lt.0.5)then
	    selectedPoints(1) = 1
	  else
	    selectedPoints(1) = 2
	  endif
	  Call RANDOM_NUMBER(randNumber)
	  if(randNumber.lt.0.5)then
	    selectedPoints(noPoints) = noInputPoints
	  else
	    selectedPoints(noPoints) = noInputPoints-1
	  endif
	endif
!pick random points
    i = 1
	do while(i.le.noPoints)
	  Call RANDOM_NUMBER(randNumber)  
	  point = ceiling(randNumber*noInputPoints)
	  pointOk = .true.
	  do j=1,noPoints
	    if(selectedPoints(j).eq.point)then
		  pointOk = .false.
		endif		
	  enddo
	  if(pointOk)then
		selectedPoints(i) = point
	    i = i + 1
	  endif
	enddo
!Transfer data to output array
    do i=1,noPoints
	  do j=1,size(inputPoints,2)
        outputPoints(i,j) = inputPoints(selectedPoints(i),j)
	  enddo
	enddo
  End Function RandomDataPoints  
  
  
  
End Module maths  