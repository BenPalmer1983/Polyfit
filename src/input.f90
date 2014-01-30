Module input

! Setup Modules
  Use kinds
  Use maths
  Use initialise


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Character(len=255)  :: inputFileName
  Integer(kind=StandardInteger) :: inputOrder
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputDataPoints
  
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runInput				    !Subroutine  
  Public :: inputFileName			    !Variable
  Public :: inputOrder		            !Variable
  Public :: inputDataPoints		        !Variable

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runInput()	
	!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
	Call readUserInput()
	Call readInputFile()	
  End Subroutine runInput  
  
!read in user input data, input from the command line
  Subroutine readUserInput()
  !force declaration of all variables
	Implicit None  
!Read in command line arguments
    call get_command_argument(1,inputFileName)
  End Subroutine readUserInput
  
!read in input file
  Subroutine readInputFile()      
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, dataCounter
	Character(len=32) :: buffera, bufferb, bufferc, bufferd
	Character(len=255) :: bufferLongA
	Character(len=255) :: fileRow
	Real(kind=DoubleReal) :: dataX, dataY
!open & read in file	
    dataCounter = 0
  	Open(UNIT=1,FILE=inputFileName) 
    do i=1,maxFileRows 
!Read in line
	  Read(1,'(A255)',IOSTAT=ios) fileRow
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
	  if(fileRow(1:16).eq."#polynomialorder")then
	    Read(1,'(A255)',IOSTAT=ios) fileRow
		Read(fileRow,*) inputOrder	  
	  endif
	  if(fileRow(1:5).eq."#data")then
!read data lines
        do j=1,maxFileRows 
	      Read(1,'(A255)',IOSTAT=ios) fileRow
		  if (ios/=0) then
	        EXIT 
	      end if
		  if (fileRow(1:1).eq."#") then
	        EXIT 
	      end if
		  dataCounter = dataCounter + 1
		enddo
	  endif 
    enddo
!close file	
	CLOSE(1)	
!Allocate array
    Allocate(inputDataPoints(1:dataCounter,1:2))	
!open file and read in data
	    dataCounter = 0
  	Open(UNIT=1,FILE=inputFileName) 
    do i=1,maxFileRows 
!Read in line
	  Read(1,'(A255)',IOSTAT=ios) fileRow
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
	  if(fileRow(1:5).eq."#data")then
!read data lines
        do j=1,maxFileRows 
	      Read(1,'(A255)',IOSTAT=ios) fileRow
		  if (ios/=0) then
	        EXIT 
	      end if
		  if (fileRow(1:1).eq."#") then
	        EXIT 
	      end if
		  dataCounter = dataCounter + 1
		  read(fileRow,*) dataX, dataY
		  inputDataPoints(dataCounter,1) = dataX
		  inputDataPoints(dataCounter,2) = dataY
		enddo
	  endif 
    enddo
!close file	
	CLOSE(1)    
  End Subroutine readInputFile   
End Module input