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
      integer status, matClose, mxIsFromGlobalWS
      integer matPutVariable, matPutVariableAsGlobal, matDeleteVariable
      integer*4 mxIsNumeric, mxIsChar

      double precision dat(9)
      data dat / 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 /
      
      mwSize M, N
      parameter(M=3) 
      parameter(N=3) 

!
!     Open MAT-file for writing
!
      write(6,*) 'Creating MAT-file matdemo.mat ...'
      mp = matOpen('matdemo.mat', 'w7.3')
      if (mp .eq. 0) then
         write(6,*) 'Can''t open ''matdemo.mat'' for writing.'
         write(6,*) '(Do you have write permission in this directory?)'
         stop
      end if
!
!     Create variables
!
      pa0 = mxCreateDoubleMatrix(M,N,0)  ! Complex flag .false. = 0
      call mxCopyReal8ToPtr(dat, mxGetPr(pa0), M*N)
!
      pa1 = mxCreateDoubleMatrix(3,3,0)
!
      pa2 = mxCreateString('MATLAB: The language of computing')
!
      pa3 = mxCreateString('MATLAB: The language of computing')
!
      status = matPutVariableAsGlobal(mp, 'NumericGlobal', pa0)
      if (status .ne. 0) then
         write(6,*) 'matPutVariableAsGlobal ''Numeric Global'' failed'
         stop
      end if
      status = matPutVariable(mp, 'Numeric', pa1)
      if (status .ne. 0) then
         write(6,*) 'matPutVariable ''Numeric'' failed'
         stop
      end if
      status = matPutVariable(mp, 'String', pa2)
      if (status .ne. 0) then
         write(6,*) 'matPutVariable ''String'' failed'
         stop
      end if
      status = matPutVariable(mp, 'String2', pa3)
      if (status .ne. 0) then
         write(6,*) 'matPutVariable ''String2'' failed'
         stop
      end if

!
!     Whoops! Forgot to copy the data into the first matrix -- 
!     it is now blank.  Well, ok, this was deliberate.  This 
!     demonstrates that matPutVariable will overwrite existing 
!     matrices.
!
      call mxCopyReal8ToPtr(dat, mxGetPr(pa1), M*N)
      status = matPutVariable(mp, 'Numeric', pa1)
      if (status .ne. 0) then
         write(6,*) 'matPutVariable ''Numeric'' failed 2nd time'
         stop
      end if

!
!     We will delete String2 from the MAT-file.
!
      status = matDeleteVariable(mp, 'String2')
      if (status .ne. 0) then
         write(6,*) 'matDeleteVariable ''String2'' failed'
         stop
      end if
!     
!     Finally, read back in MAT-file to make sure we know what we put
!     in it.
!
      status = matClose(mp)
      if (status .ne. 0) then
         write(6,*) 'Error closing MAT-file'
         stop
      end if
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
!
      pa1 = matGetVariable(mp, 'Numeric')
      if (mxIsNumeric(pa1) .eq. 0) then
         write(6,*) 'Invalid non-numeric matrix written to MAT-file'
         stop
      end if
!

      pa2 = matGetVariable(mp, 'String')

      if (mxIsChar(pa2) .eq. 0) then
         write(6,*) 'Invalid non-string matrix written to MAT-file'
         stop
      end if
!
      pa3 = matGetVariable(mp, 'String2')
      if (pa3 .ne. 0) then
         write(6,*) 'String2 not deleted from MAT-file'
         stop
      end if
!
!     clean up memory
      call mxDestroyArray(pa0)
      call mxDestroyArray(pa1)
      call mxDestroyArray(pa2)
      call mxDestroyArray(pa3)

      status = matClose(mp)
      if (status .ne. 0) then
         write(6,*) 'Error closing MAT-file'
         stop
      end if
!
      write(6,*) 'Done creating MAT-file'
      stop
      end
