! This is a Fortran program that read FDS 1st combined binary data 
! and write them into MATLAB ".mat" format file. 
!
#include "fintrf.h"

Program fds2mat

   implicit none

   ! Variable annocement
   integer, parameter :: fb = selected_real_kind(6)
   character(255) :: infilename, outfilename
   integer, parameter :: time_dim = 1001
   integer ::  infileid, outfileid
   integer :: nmlstat, infilestat, outfilestat
   integer :: nx, ny, nz, nt
   logical :: lwtime, lwquan

   character(255) :: varname, varunit, regularname
   integer:: tdim1,tdim2
   integer:: xdim1,xdim2,ydim1,ydim2,zdim1,zdim2
   integer:: i,j,k,l
   integer*8:: ii
   real(fb), allocatable, dimension(:):: time
   real(fb), allocatable, dimension(:):: xcor,ycor,zcor
   real(fb), allocatable, dimension(:,:,:,:):: quan
   real(fb), allocatable, dimension(:):: quanline

   character(255) :: frmt

   ! Annoncement of Matlab Fortran Interface Pointers
   mwPointer matOpen
   mwPointer mxCreateNumericArray, mxCreateString
   mwPointer mxGetNumberOfElements, mxGetDimensions, mxGetData
   mwPointer matGetVariable

   mwPointer mpfile  ! file handle
   mwPointer ptime, pquan  !var handle

   mwSize mndim,mlength,mtdim       ! integer*8 type
   parameter (mndim=4)  ! dimension of quantity
   parameter(mtdim=1)
   mwSize mdims(mndim)
   mwSize mtlen(mtdim)

   integer mstat, matClose, mxIsFromGlobalWS
   integer*4 matPutVariable, matPutVariableAsGlobal, matDeleteVariable
   integer*4 classid, complexflag
   integer*4 mxClassIDFromClassName
   character*(32) :: classname='single'


   ! Open fds2slcf output binary file
   !
   ! ================================================
   ! File structure: 
   !
   !------------- File header ----------- 
   ! character(255):: Variable name
   ! character(255):: Variable Unit

   !------------- Dimensions ------------
   ! integer::  tdim1, tdim2  (Time dimension)
   ! integer:: xdim1,xdim2,ydim1,ydim2,zdim1,zdim2 (Space dimension)
   !
   !------------- Coordinates -----------
   ! real,dimension(:) time(i), i=tdim1,tdim2   (Time serials)
   ! real,dimension(:) xcor(i), i=xdim1,xdim2   (X coordinates)
   ! real,dimension(:) ycor(i), i=ydim1,ydim2   (Y coordinates)
   ! real,dimension(:) zcor(i), i=zdim1,zdim2   (Z coordinates)
   !
   !------------- Quantity Data Block ------------
   !! real, dimension(:,:,:,:) 
   ! ((((quan(i,j,k,nt),i=nqx1,nqx2),j=nqy1,nqy2),k=nqz1,nqz2),nt=t1,t2)
   ! 
   ! End of file ================================================


   ! Some initial values
   infileid=20
   outfileid=30
   print*, "===================================================================="
   print*
   print*, "Start to executing the program fds2mat"

   ! Read namelist file
   write(*,*) "Enter the input file name: "
   read(*,'(A)')  infilename 
   print *
   write(*,*) "Enter the output file name:"
   read(*,'(A)') outfilename

   ! Read input data file
   open(unit=infileid, file=infilename, form="unformatted", status='unknown',iostat=infilestat)
   if ( infilestat /= 0 ) then
      write(*,*) "*** Error open ",trim(infilename),", Please check the access of inputfile"
      stop
   endif

   write(*,*) 'Reading binary data from file:'
   write(*,'(4X,A)') trim(infilename)
   ! Read file header
   read(infileid,end=99) varname  
   read(infileid,end=99) varunit 

   ! Read dimensions 
   read(infileid,end=99) tdim1, tdim2
   read(infileid,end=99) xdim1,xdim2,ydim1,ydim2,zdim1,zdim2

   nx=xdim2-xdim1+1
   ny=ydim2-ydim1+1
   nz=zdim2-zdim1+1
   nt=tdim2-tdim1+1

   allocate(time(1:nt))
   allocate(xcor(1:nx))
   allocate(ycor(1:ny))
   allocate(zcor(1:nz))
   allocate(quan(1:nx,1:ny,1:nz,1:nt))

   ! Read coordinates in time and space, Sequence should not change
   read(infileid,end=99) (time(i),i=1, nt)
   read(infileid,end=99) (xcor(i),i=1, nx)
   read(infileid,end=99) (ycor(i),i=1, ny)
   read(infileid,end=99) (zcor(i),i=1, nz)

   ! Read quantity data block
   read(infileid,end=99) ((((quan(i,j,k,l),i=1,nx),j=1,ny),k=1,nz),l=1,nt)

   99 close(infileid)

   ! Reshape 4-D quan to 1-d quanline
   allocate(quanline(1:size(quan)))
   ii=1
   do l=1,nt
      do k=1,nz
         do j=1,ny
            do i=1,nx
               quanline(ii)=quan(i,j,k,l)
               ii=ii+1
            end do
         end do
      end do
   end do


   ! ---------- Data check -------------
   ! Dimension output
   print *
   write(*,*) "Dimensions  is:"
   frmt='(3X,A,3X,I4)'
   write(*,frmt) "X-axis dimension:",nx
   write(*,frmt) "Y-axis dimension:",ny
   write(*,frmt) "Z-axis dimension:",nz
   write(*,frmt) "Time dimension:",nt

   !
   ! ===========================================================================
   ! *                      Write values to .mat file
   ! ============================================================================
   !

   ! Open MAT-file for writing
   write(*,*) "Creating MAT-file: ",trim(outfilename)
   mpfile = matOpen(outfilename, 'u')
   if (mpfile == 0) then
      print*
      write(*,*) "Can't open: ",trim(outfilename)
      write(*,*) "  file may not exist, Try creating one with write purpose"
      mpfile = matOpen(outfilename, 'w')
      if (mpfile == 0 ) then
         print*
         write(*,*) "Can''t open ",trim(outfilename), "directory may not accessable"
         write(*,*) "Program will exit"
         stop
      end if
   end if

   ! Check time field
   ptime = matGetVariable(mpfile, 'time')
   if (ptime == 0 ) then
      print*
      write(*,*) "'Time' field doesn't exist, writing..."
      lwtime = .true.
   else 
      write(*,*) "'Time' field exists, skip..."
      lwtime = .false.
   endif

   call namecheck(varname, regularname)  ! varname check 

   pquan = matGetVariable(mpfile, trim(regularname))
   if (pquan == 0 ) then
      print *
      write(*,*) "'",trim(regularname),"' field doesn't exist, writing ..."
      lwquan = .true.
   else
      write(*,*) "'",trim(regularname),"' field exists, skip ..."
      lwquan = .false.
   end if

   ! Prepare dimension for MATLAB variable creation 
   classid=mxClassIDFromClassName(classname)
   complexflag=0

   ! define some dimensions
   mdims(1)=nx
   mdims(2)=ny
   mdims(3)=nz
   mdims(4)=nt
   mlength=nx*ny*nz*nt

   mtlen=nt

   ! Create MATLAB object variables with point and transfer fortran values to them
   print*


   if (lwtime) then
      ptime = mxCreateNumericArray(mtdim,mtlen,classid,complexflag)
      call mxCopyReal4ToPtr(time, mxGetData(ptime), mtlen)
      mstat = matPutVariable(mpfile, 'time', ptime)  ! Variable name must not have '-', '='
      write(*,*) "Writing NumericArray 'time'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable 'time'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable 'time' succeeded"
   endif

   ! Write Values to matfile pointed by mpfile
   if (lwquan) then
      pquan = mxCreateNumericArray(mndim,mdims,classid,complexflag)
      call mxCopyReal4ToPtr(quanline, mxGetData(pquan), mlength)
      mstat = matPutVariable(mpfile, trim(regularname), pquan)  
      write(*,*) "Writing NumericArray '",trim(regularname),"'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable '",trim(regularname),"'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable '",trim(regularname), "' succeeded"
   end if


   1000 write(*,*) "--- Closing MAT-file, Program Ends"
   print*
   mstat = matClose(mpfile)
   if (mstat /= 0) then
      write(*,*) "***Error closing MAT-file, Program Stop"
      print*
      stop
   endif
   stop

end program fds2mat


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
!                      SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine namecheck(str_in, str_out)
   ! This subroutine eliminate informal characters in Strings act as
   ! MATLAB variable name.
   character(255) :: str_in, str_out, temp
   character(1) :: space = ' '
   character(1) :: minus = '-'
   character(1) :: equal = '='
   character(1) :: plus = '+'

   character(1) :: dash = '_'

   integer :: ilen,i
   logical :: iffind

   temp=adjustl(str_in)
   do  i=1,len_trim(temp)
      if(temp(i:i) .eq. space ) then
         temp(i:i) = dash
      else if (temp(i:i) .eq. minus) then
         temp(i:i) = dash
      else if (temp(i:i) .eq. equal) then
         temp(i:i) = dash
      else if (temp(i:i) .eq. plus) then
         temp(i:i) = dash
      end if
   end do
   str_out=temp
   return 

end subroutine namecheck

