! This is a Fortran program that read FDS 1st combined binary data 
! and write them into MATLAB ".mat" format file. 
!
!#include "fintrf.h"

Program fds2mat

   ! Variable annocement
   character(255) :: infilename, outfilename, nmlname
   integer, parameter :: fb = selected_real_kind(6)
   integer, parameter :: time_dim = 1001
   integer :: nmlid, infileid, outfileid
   integer :: nmlstat, infilestat, outfilestat
   character(255) ::  quantname
   namelist /quantity/ infilename, outfilename

   character(255) :: varname, varunit
   integer:: tdim1,tdim2
   integer:: xdim1,xdim2,ydim1,ydim2,zdim1,zdim2
   integer:: i,j,k,l
   real, allocatable, dimension(:):: time
   real, allocatable, dimension(:):: xcor,ycor,zcor
   real, allocatable, dimension(:,:,:,:):: quan

   character(255) :: frmt

   ! Some initial values

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
   nmlid=10
   infileid=20
   outfileid=30


   ! Read namelist file
   write(*,*) " Enter the namelist file name: "
   read(*,'(A)')  nmlname 
   open(unit=nmlid, file=nmlname, form="formatted", status='unknown',iostat=nmlstat)
   if ( nmlstat /= 0 ) then
      write(*,*) "*** Error open ",trim(nmlname), ", File may not exist"
      stop
   endif

   read(nmlid, nml=quantity)
   write(*,*) "Read namelist quantity:"
   write(*,nml=quantity)
   close(nmlid)

   ! Read input data file
   open(unit=infileid, file=infilename, form="unformatted", status='unknown',iostat=infilestat)
   if ( infilestat /= 0 ) then
      write(*,*) "*** Error open ",trim(infilename),", Please check the access of inputfile"
      stop
   endif

   write(*,*) ' Reading binary data from file :',trim(infilename)
   ! Read file header
   read(infileid,end=99) varname  
   read(infileid,end=99) varunit 

   ! Read dimensions 
   read(infileid,end=99) tdim1, tdim2
   read(infileid,end=99) xdim1,xdim2,ydim1,ydim2,zdim1,zdim2

   allocate(time(tdim1:tdim2))
   allocate(xcor(xdim1:xdim2))
   allocate(ycor(ydim1:ydim2))
   allocate(zcor(zdim1:zdim2))
   allocate(quan(xdim1:xdim2, ydim1:ydim2, zdim1:zdim2, tdim1:tdim2))

   ! Read coordinates in time and space
   read(infileid,end=99) (time(i),i=tdim1, tdim2)

   read(infileid,end=99) (xcor(i),i=xdim1, xdim2)
   read(infileid,end=99) (ycor(i),i=ydim1, ydim2)
   read(infileid,end=99) (zcor(i),i=zdim1, zdim2)

   ! Write quantity data block
   read(infileid,end=99) ((((quan(i,j,k,l),i=xdim1,xdim2),j=ydim1,ydim2),k=zdim1,zdim2),l=tdim1,tdim2)


   99 close(infileid)

   ! ---------- Data check -------------
   ! 1. Dimension Check
   print *
   frmt="1X,A,1X,I1,':',I4"
   write(*,*) "Dimension Read from ",trim(infilename)," is:"
   write(*,frmt) "Time dimension:",tdim1,tdim2
   write(*,frmt) "X-axis dimension:",xdim1,xdim2
   write(*,frmt) "Y-axis dimension:",ydim1,ydim2
   write(*,frmt) "Z-axis dimension:",zdim1,zdim2




end program fds2mat
