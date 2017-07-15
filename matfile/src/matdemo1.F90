#include "fintrf.h"
!
#if 0
!
!     matdemo1.F
!     .F file need to be preprocessed to generate .for equivalent
!
#endif
!
!     matdemo1.f
!
!     This is a simple program that illustrates how to call the MATLAB
!     MAT-file functions from a FORTRAN program.  This demonstration
!     focuses on writing MAT-files.
!
! Copyright 1984-2011 The MathWorks, Inc.
!=====================================================================
! 

!     matdemo1 - create a new MAT-file from scratch.
!
program matdemo1

   !     Declarations
   implicit none

   mwPointer matOpen, mxCreateDoubleMatrix, mxCreateString 
   mwPointer matGetVariable, mxGetPr
   mwPointer mp, pa1, pa2, pa3, pa0
   integer stat, matClose, mxIsFromGlobalWS
   integer matPutVariable, matPutVariableAsGlobal, matDeleteVariable
   integer*4 mxIsNumeric, mxIsChar

   double precision dat(12)
   data dat / 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 /

   mwSize M, N
   parameter(M=3) 
   parameter(N=4) 


   !  Open MAT-file for writing

   write(6,*) 'Creating MAT-file matdemo.mat ...'
   mp = matOpen('matdemo.mat', 'w')
   if (mp .eq. 0) then
      write(6,*) 'Can''t open ''matdemo.mat'' for writing.'
      write(6,*) '(Do you have write permission in this directory?)'
      stop
   end if

   ! Create variables

   pa0 = mxCreateDoubleMatrix(M,N,0)  ! Complex flag .false. = 0, Empty matrix
   call mxCopyReal8ToPtr(dat, mxGetPr(pa0), M*N)

   pa1 = mxCreateDoubleMatrix(2,6,0)  ! Create Empty matrix

   pa2 = mxCreateString('MATLAB: The language of computing')

   pa3 = mxCreateString('MATLAB: The language of computing')

   stat = matPutVariableAsGlobal(mp, 'NumericGlobal', pa0) 
   if (stat .ne. 0) then
      write(6,*) 'matPutVariableAsGlobal ''Numeric Global'' failed'
      stop
   end if
   write(*,*) 'matPutVariableAsGlobal ''Numeric Global'' succeeded'

!   stat = matPutVariable(mp, 'Numeric', pa1)
!   if (stat .ne. 0) then
!      write(6,*) 'matPutVariable ''Numeric'' failed'
!      stop
!   end if
!      write(6,*) 'matPutVariable ''Numeric'' succeeded'

   stat = matPutVariable(mp, 'String', pa2)
   if (stat .ne. 0) then
      write(6,*) 'matPutVariable ''String'' failed'
      stop
   end if
      write(6,*) 'matPutVariable ''String'' succeeded'

   stat = matPutVariable(mp, 'String2', pa3)
   if (stat .ne. 0) then
      write(6,*) 'matPutVariable ''String2'' failed'
      stop
   end if
      write(6,*) 'matPutVariable ''String2'' succeeded'

   !
   !     Whoops! Forgot to copy the data into the first matrix -- 
   !     it is now blank.  Well, ok, this was deliberate.  This 
   !     demonstrates that matPutVariable will overwrite existing 
   !     matrices.
   !
   call mxCopyReal8ToPtr(dat, mxGetPr(pa1), 2*6)
   stat = matPutVariable(mp, 'Numeric', pa1)
   if (stat .ne. 0) then
      write(6,*) 'matPutVariable ''Numeric'' failed 2nd time'
      stop
   end if
      write(6,*) 'matPutVariable ''Numeric'' succeeded 2nd time'

   !
   !     We will delete String2 from the MAT-file.
   !
   stat = matDeleteVariable(mp, 'String2')
   if (stat .ne. 0) then
      write(6,*) 'matDeleteVariable ''String2'' failed'
      stop
   end if
      write(6,*) 'matDeleteVariable ''String2'' succeeded'
   !     
   !     Finally, read back in MAT-file to make sure we know what we put
   !     in it.
   !
   stat = matClose(mp)
   if (stat .ne. 0) then
      write(6,*) 'Error closing MAT-file'
      stop
   end if
      write(6,*) 'Succeeded closing MAT-file'
   !
   mp = matOpen('matdemo.mat', 'r')
   if (mp .eq. 0) then
      write(6,*) 'Can''t open ''matdemo.mat'' for reading.'
      stop
   end if
   !
   pa0 = matGetVariable(mp, 'NumericGlobal')
   if (mxIsFromGlobalWS(pa0) .eq. 0) then
      write(6,*) 'Invalid non-global matrix written to MAT-file'
      stop
   end if
   
   pa1 = matGetVariable(mp, 'Numeric')
   if (mxIsNumeric(pa1) .eq. 0) then
      write(6,*) 'Invalid non-numeric matrix written to MAT-file'
      stop
   end if
   

   pa2 = matGetVariable(mp, 'String')

   if (mxIsChar(pa2) .eq. 0) then
      write(6,*) 'Invalid non-string matrix written to MAT-file'
      stop
   end if
   
   pa3 = matGetVariable(mp, 'String2')
   if (pa3 .ne. 0) then
      write(6,*) 'String2 not deleted from MAT-file'
      stop
   end if
   
   !     clean up memory
   call mxDestroyArray(pa0)
   call mxDestroyArray(pa1)
   call mxDestroyArray(pa2)
   call mxDestroyArray(pa3)

   stat = matClose(mp)
   if (stat .ne. 0) then
      write(6,*) 'Error closing MAT-file'
      stop
   end if
   
   write(6,*) 'Done creating MAT-file'
   stop
   end

