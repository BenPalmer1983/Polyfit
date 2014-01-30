Program polyfitprog

! University of Birmingham
! Ben Palmer
!


!Setup Modules
Use kinds				!data kinds
Use initialise			!initialise program
Use maths				!maths functions
Use input			!input


!force declaration of all variables
  Implicit None

!declare variables
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
  Real(kind=DoubleReal) :: rss, fittingTime
  Integer(kind=StandardInteger) :: i
  
!initialise and run input
  Call runInitialise()  
  Call runInput
  
!fit data
  fittingTime = ProgramTime()
  Allocate(coefficients(0:inputOrder))
  coefficients = PolyFit(inputDataPoints,inputOrder)
  rss = CalcResidualSquareSum(inputDataPoints,coefficients)  
  fittingTime = ProgramTime() - fittingTime
  
!Output coefficients
  do i=0,inputOrder
    print *,"x^",i,"*",coefficients(i)
  enddo 
  print *,"Rss: ",rss
  print *,"Time:",fittingTime
End