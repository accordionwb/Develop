! This file implemented MATLAB matrix API
! Environmental Variables LD_LIBRARY_PATH must contains $MATLAB_ROOT/bin/glnxa64
!
#include "fintrf.h"

Program SLCF2MAT

   implicit none

   ! Variable annocement
   integer, parameter :: fb = selected_real_kind(6)
   character(255) :: infilename, matfilename
   integer, parameter :: time_dim = 1002
   integer, parameter :: file_dim = 500
   integer ::  infileid, outfileid
   integer :: nmlstat, infilestat, outfilestat
   integer :: nx, ny, nz, nt
   integer :: ierr, NMESHES, nm, noc
   integer :: i, j, k, l, ix, iy, iza
   integer :: kk
   integer:: tdim1,tdim2
   integer:: xdim1,xdim2,ydim1,ydim2,zdim1,zdim2
   integer :: idum, ifile, nsam, nv, mv
   integer :: imin,imax,jmin,jmax,kmin,kmax
   integer :: i1, i2, j1, j2, k1, k2, t1, t2, i3, j3, k3  ! dimension
   integer :: iix,iiy,iiz,nix,niy,niz
   integer :: i10, i20, j10, j20, k10, k20
   integer :: ncount, ior_input, npatch, ijbar, jkbar
   integer :: ii, nxp, nyp, nzp, n, i_sample,jj
   integer :: allerr,iuser_label,idummy
   integer :: nqx1,nqy1,nqz1,nqx2,nqy2,nqz2
   logical :: lwtime, lwquan, lwdims, lwxor, lwyor, lwzor

   character(255) ::   regularname
   integer*4, dimension(4) :: dims
   real(fb), allocatable, dimension(:):: time
   real(fb), allocatable, dimension(:) :: cqx,cqy,cqz
   real(fb), allocatable, dimension(:) :: xcor,ycor, zcor
   real(fb), allocatable, dimension(:,:,:,:):: quantity
   real(fb), allocatable, dimension(:):: quanline

   ! Annoncement of Matlab Fortran Interface Pointers
   mwPointer matOpen
   mwPointer mxCreateNumericArray, mxCreateString
   mwPointer mxGetNumberOfElements, mxGetDimensions, mxGetData
   mwPointer matGetVariable

   mwPointer mpfile  ! file handle
   mwPointer ptime, pquan, pdims, pxor, pyor, pzor !var handle

   mwSize mqdim,mlength,mtdim, mddim, mxdim, mydim, mzdim      ! integer*8 type
   parameter (mqdim=4)  ! dimension of quantity
   parameter (mtdim=1)
   parameter (mddim=1)
   parameter (mxdim=1)
   parameter (mydim=1)
   parameter (mzdim=1)
   mwSize mdims(mqdim)
   mwSize mtlen(mtdim)
   mwSize mdlen(mddim)
   mwSize mxlen(mxdim)
   mwSize mylen(mydim)
   mwSize mzlen(mzdim)

   integer mstat, matClose, mxIsFromGlobalWS
   integer*4 matPutVariable, matPutVariableAsGlobal, matDeleteVariable
   integer*4 classid1, classid2, complexflag
   integer*4 mxclassidFromclassname
   character*(8) :: realname='single'
   character*(8) :: intname='int32'



   interface
      subroutine parse(buffer,sb_toks,se_toks,n_toks)
         implicit none
         character(*), intent(inout) :: buffer
         integer, dimension(*), intent(out) :: sb_toks, se_toks
         integer, intent(out) :: n_toks
      end subroutine parse
   end interface

   real(fb) :: xs, xf, ys, yf, zs, zf  ! user defined bounds
   real(fb) :: d1, d2, d3, d4, tbeg, tend
   real(fb), dimension(file_dim) :: xtemp,ytemp,ztemp
   character(256) :: auto_slice_label,auto_slice_unit
   integer :: auto_slice_flag, n_auto_slices, iauto, iztemp
   integer, dimension(1) :: izmin1,ixmin1,iymin1
   integer :: izmin,ixmin,iymin
   real(fb):: ixoff,iyoff,izoff,Dspace
   integer, dimension(file_dim) :: auto_slice_lists,temp_slice_lists
   real(fb), dimension(file_dim) :: auto_slice_z,auto_slice_y,auto_slice_x
   real(fb) :: ax1, ay1, az1, az2
   real(fb) :: eps, zmin


   type mesh_type
      real(fb), pointer, dimension(:) :: x,y,z
      real(fb) :: d1,d2,d3,d4
      integer :: ibar,jbar,kbar,ierr,nxp,nyp,nzp
   end type mesh_type

   type (mesh_type), dimension(:), allocatable, target :: mesh
   type(mesh_type), pointer :: m

   real(fb), dimension(60) :: sum
   integer, allocatable, dimension(:) :: ix1,ix2,iy1,iy2,iz1,iz2
   integer, allocatable, dimension(:) :: auto_slices
   real(fb), allocatable, dimension(:,:,:,:) :: q
   real(fb), allocatable, dimension(:,:,:) :: f
   logical, allocatable, dimension(:,:,:) :: already_used
   logical :: new_plot3d=.true., iflist=.false., iffind
   integer ior_slcf, stationID
   character(1) :: ans,seq
   character(4) :: ext1, ext2
   character(4) :: choice
   character(30) :: unitjunk
   character(256) gridfile,qname,chid,qfile,junk,frmt,slcf_label_dummy
   character(256) outfile, outfile1, outfile2, nonspace
   character(256), dimension(file_dim) :: pl3d_file,slcf_file,bndf_file,slcf_text,bndf_text,slcf_unit,bndf_unit,slcf_label
   character(256), dimension(file_dim) :: slice_labels, slice_units
   character(20), dimension(file_dim) :: bndf_type
   character(20) :: bndf_type_chosen
   integer :: nslice_labels, slice_exist, new_slice_seq
   real(fb), dimension(file_dim) :: x1, x2, y1, y2, z1, z2
   integer,  dimension(60) :: ib,is
   integer,  dimension(file_dim) :: slcf_mesh
   logical, dimension(file_dim) :: file_exists
   integer :: rcode
   integer :: nfiles_exist
   real(fb) :: zoffset
   integer :: zoffset_flag
   logical :: exists
   integer :: lu_in, nargs, iarg, lenstring 
   character(256) :: buffer, filein
   integer, dimension(256) :: sb_toks, se_toks
   integer :: n_toks
   integer :: Error_status
   character(256) :: arg

   zoffset=0.0
   zoffset_flag=0
   eps=0.001
   filein='stdin'

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

   write(*,*) "======================= Program Start =========================="
   print *
   write(*,*) ' Enter job id string (chid):'

   read(*,'(a)') chid

   ! check to see if the .smv file exists

   gridfile = trim(chid)//'.smv'
   inquire(file=gridfile,exist=exists)
   if (.not.exists) then
      write(6,*)"*** fatal Error: the file: ",trim(gridfile)," does not exist"
      stop
   end if

   ! open the .smv file

   open(11,file=gridfile,status='old',form='formatted')
   print*
   write(*,*) "Open SMV file: ",trim(gridfile)


   ! determine the number of meshes

   rewind(11)

   call search('NMESHES',7,11,ierr)
   if (ierr.eq.1) then
      write(6,*) ' warning: assuming 1 mesh'
      NMESHES = 1
   else
      read(11,*) NMESHES
      write(*,'(A,I3)') "Number of mesh is: ",NMESHES
   end if

   allocate(mesh(NMESHES))

   ! get the coordinates of the meshes

   rewind(11)
   read_smv: do nm=1,NMESHES

      m=>mesh(nm)
      call search('GRID',4,11,ierr)
      read(11,*) m%ibar,m%jbar,m%kbar
      allocate(m%x(0:m%ibar))
      allocate(m%y(0:m%jbar))
      allocate(m%z(0:m%kbar))

      call search('TRNX',4,11,ierr)
      read(11,*) noc
      do i=1,noc
         read(11,*)
      end do
      do i=0,m%ibar
         read(11,*) idum,m%x(i)
      end do

      call search('TRNY',4,11,ierr)
      read(11,*) noc
      do i=1,noc
         read(11,*)
      end do
      do j=0,m%jbar
         read(11,*) idum,m%y(j)
      end do

      call search('TRNZ',4,11,ierr)
      read(11,*) noc
      do i=1,noc
         read(11,*)
      end do
      do k=0,m%kbar
         read(11,*) idum,m%z(k)
      end do

   end do read_smv

   if (iflist) then
      write(6,*) ' domain selection:'
      write(6,*) '   y - domain size is limited'
      write(6,*) '   n - domain size is not limited'
      write(6,*) '   z - domain size is not limited and z levels are offset'
      read(*,'(a)') ans
      print *
   else
      ans = 'n'
   end if

   call toupper(ans,ans)
   if (ans(1:1).eq.'y') then
      write(6,*) ' Enter min/max x, y and z'
      read(*,*) xs,xf,ys,yf,zs,zf
   else
      xs = -100000.
      xf =  100000.
      ys = -100000.
      yf =  100000.
      zs = -100000.
      zf =  100000.
   end if


   ! if ans is z or z then subtract zoffset from z data (for multi-level cases)

   if (ans(1:1).eq.'z') then
      zoffset_flag=1
   else
      zoffset_flag=0
   end if


   !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   !
   !  Key ieteration search slcf (i), list slcf file and variables
   !		1. Print all exist slice file and their contains.
   !  	2. Generate unique slice variable list.
   !
   !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   ! Set some initial values
   slcf_label = 'null'
   slcf_mesh = 1
   rewind(11)
   nfiles_exist=0  ! record variables cover multiple slice files
   nslice_labels=0

   Search_slcf: do i=1,file_dim    

      ! General check and read SLCF headers
      call search2('SLCF',4,'SLCC',4,11,ierr,choice)
      if (ierr.eq.1) exit search_slcf
      backspace(11)
      read(11,*) junk,slcf_mesh(i)  ! Slcf dimensions
      read(11,'(a)') slcf_file(i)  ! slcf file name
      read(11,'(a)') slcf_text(i)  ! slcf variable defination

      ! Check if slice file exist
      read(11,*)    ! skip dummy line
      read(11,'(a)') slcf_unit(i)
      open(12,file=trim(slcf_file(i)),form='unformatted',status='old', iostat=rcode)
      if (rcode.ne.0) then
         close(12)
         cycle ! jump to search_slcf for another run
      end if

      ! record of exist slice files
      nfiles_exist=nfiles_exist+1  

      ! TASK 1: Read and Print what individual slice file contains.

      read(12) unitjunk
      read(12) unitjunk
      read(12) unitjunk
      read(12) i1,i2,j1,j2,k1,k2
      close(12)

      nm=slcf_mesh(i)
      m=>mesh(nm)
      x1(i)=m%x(i1)
      x2(i)=m%x(i2)
      y1(i)=m%y(j1)
      y2(i)=m%y(j2)
      z1(i)=m%z(k1)
      z2(i)=m%z(k2)

      if (iflist) then
         write(6,'(i3,1x,a,1x,a)')i,trim(slcf_text(i)),trim(slcf_file(i))
         write(6,'(3x,a,6(1x,f8.2))')'slice bounds:',m%x(i1),m%x(i2),m%y(j1),m%y(j2),m%z(k1),m%z(k2)
         write(6,*) 
      end if

      ! TASK 2: Create unique list of slice types      

      slice_exist=0

      do ii=1, nslice_labels ! loop over labels kind
         if(trim(slcf_text(i)).eq.slice_labels(ii))then
            slice_exist=1 ! set flag slice already recorded
            exit  
         end if
      end do

      if (slice_exist.eq.0) then ! new slice
         nslice_labels=nslice_labels+1  ! add record 
         slice_labels(nslice_labels)=trim(slcf_text(i)) ! Record new slice voariable
         slice_units(nslice_labels)=trim(slcf_unit(i))
      end if

   end do Search_slcf

   ! Error code: slice file not found
   write(*,'(A,I3)') "Number of slcf files: ",nfiles_exist
   IF (nfiles_exist.EQ.0)THEN
      print *
      WRITE(6,*)"There are no slice files to convert"
      STOP
   end if

   close(11)   ! close .svm file

   !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   !
   !  1. Print unique slice variable list 
   !  2. Select their stationary dimension 
   !  3. Adjust the mesh sequence to insure dimmension increases in all direction
   !
   !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   ! Task 1: Print the unique slice variable list.
   N_AUTO_SLICES=0
   print *
   write(6,*)' Enter slice file type index'
   DO II=1, nslice_labels
      write(6,'(2x,I3,1x,A)')II,TRIM(slice_labels(II))
   end do
   READ(*,*) iuser_label
   AUTO_SLICE_LABEL=TRIM(slice_labels(iuser_label))
   auto_slice_unit=trim(slice_units(iuser_label))

   ! Task 2: Select stationary dimension.
   print *
   write(*,*) "Choose the 2D slice plate (1/2/3):"
   write(*,'(3X,A)')  "1   Y-Z Plate"
   write(*,'(3X,A)')  "2   X-Z Plate"
   write(*,'(3X,A)')  "3   X-Y Plate"
   read(*,*) StationID
   print *
   write(*,*) "----The original Mesh sequence is: ----"

   n_auto_slices = 0

   Select case(stationID)

   Case(1) ! Y-Z plate
      boundloop1: DO I=1,nfiles_exist ! loop over slice files, equal to (I) in Search_slcf
         IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(I))) CYCLE boundloop1
         IF(x2(I)-x1(I).GT.EPS)CYCLE boundloop1

         iix=0
         do kk=1,nix
            if (abs(xtemp(kk) - x1(I)) .lt. eps) then
               iix=1
               exit
            end if
         end do
         if (iix .eq. 0) then
            nix=nix+1
            xtemp(nix)=x1(I)
         end if

         n_auto_slices = n_auto_slices + 1
         nm=slcf_mesh(I)
         m=>mesh(nm)
         write(6,'(A,I3,5(a,f8.3))')"Mesh:",nm,", ymin=",y1(I),", ymax=",y2(I),", zmin=",z1(I),", zmax=",z2(I),", x=",x1(I)
         write(6,'(3(3X,A,I3))') "The mesh size on each axis is: NX =",m%ibar,", NY=",m%jbar,", NZ=",m%kbar
         auto_slice_lists(n_auto_slices) = I
      end do boundloop1
      if (n_auto_slices .eq. 0) then
         write(*,*) "*** There is no slice found on Y-Z plate, Program Stop! "
         stop
      end if

   case(2) ! X-Z plate
      boundloop2: DO I=1,nfiles_exist ! loop over slice files, equal to (I) in Search_slcf
         IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(I))) CYCLE boundloop2
         IF(y2(I)-y1(I).GT.EPS)CYCLE boundloop2

         iiy=0
         do kk=1,niy
            if (abs(ytemp(kk) - y1(I)) .lt. eps) then
               iiz=1
               exit
            end if
         end do
         if (iiy .eq. 0) then
            niy=niy+1
            ytemp(niy)=x1(I)
         end if

         n_auto_slices = n_auto_slices + 1
         nm=slcf_mesh(I)
         m=>mesh(nm)
         write(6,'(A,I3,5(a,f8.3))')"Mesh:",nm,", xmin=",x1(I),", xmax=",x2(I),", zmin=",z1(I),", zmax=",z2(I),", y=",y1(I)
         write(6,'(3(3X,A,I3))') "The mesh size on each axis is: IX =",m%ibar,", IY=",m%jbar,", IZ=",m%kbar
         auto_slice_lists(n_auto_slices) = I
      end do boundloop2
      if (n_auto_slices .eq. 0) then
         write(*,*) "*** There is no slice found on X-Z plate, Program Stop! "
         stop
      end if

   case(3) ! X-Y plate
      ztemp=0.0
      niz=0
      boundloop3: DO I=1,nfiles_exist ! loop over slice files, equal to (I) in Search_slcf
         IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(I))) CYCLE boundloop3
         IF(z2(I)-z1(I).GT.EPS)CYCLE boundloop3

         iiz=0
         do kk=1,niz
            if (abs(ztemp(kk) - z1(I)) .lt. eps) then
               iiz=1
               exit
            end if
         end do
         if (iiz .eq. 0) then
            niz=niz+1
            ztemp(niz)=z1(I)
         end if

         n_auto_slices = n_auto_slices + 1
         nm=slcf_mesh(I)
         m=>mesh(nm)
         write(6,'(A,I3,5(a,f8.3))')"Mesh:",nm,", xmin=",x1(I),", xmax=",x2(I),", ymin=",y1(I),", ymax=",y2(I),", z=",z1(I)
         write(6,'(3(3X,A,I3))') "The mesh size on each axis is: IX =",m%ibar,", IY=",m%jbar,", IZ=",m%kbar
         auto_slice_lists(n_auto_slices) = I
      end do boundloop3
      if (n_auto_slices .eq. 0) then
         write(*,*) "*** There is no slice found on X-Y plate, Program Stop! "
         stop
      end if

   end Select

   ! Task 3:  Adjust mesh list sequence

   allocate(ix1(n_auto_slices))   ! array storing x y z coordinate
   allocate(ix2(n_auto_slices))
   allocate(iy1(n_auto_slices))
   allocate(iy2(n_auto_slices))
   allocate(iz1(n_auto_slices))
   allocate(iz2(n_auto_slices))

   print *
   write(*,*) "Need to adjust the coordinates offset? (y/n)"
   read(*,*) seq

   if (seq .eq. 'y' .or. seq .eq. 'Y') then
      write(*,*) "Provide approprate coordinates offset (xoff,yoff,zoff) and Dspace"
      read(*,*) ixoff,iyoff,izoff,Dspace
      print *
      write(*,*) "----New mesh coordinates after adjustment----"

      ! Select case

      select case(stationID)

      case(1)
         do J=1,n_auto_slices
            I=auto_slice_lists(J)
            do kk=1,nix
               if (abs(xtemp(kk)-x1(I)).lt. eps) then
                  ix1(J)=kk;
                  ix2(J)=kk;
                  exit
               end if
            end do
            iy1(J)=int((y1(I)-iyoff)/Dspace)
            iy2(J)=int((y2(I)-iyoff)/Dspace)
            iz1(J)=int((z1(I)-izoff)/Dspace)
            iz2(J)=int((z2(I)-izoff)/Dspace)
            write(6,'(A,I3,A,I3,5(a,f8.3))')"Seq:",J,", Mesh:",slcf_mesh(I),", ymin=",y1(I),", ymax=",y2(I),", zmin=",z1(I),", zmax=",z2(I),", x=",x1(I)
            write(6,'(3(3X,A,I3,2X,I3))') "The mesh size on each axis is: IX =",ix1(j),ix2(j),", IY=",iy1(j),iy2(j),", IZ=",iz1(j),iz2(J)
            if (J .eq. 1) then
               nqx1=ix1(J)
               nqy1=iy1(J)
               nqz1=iz1(J)
            end if
            if (J .eq. n_auto_slices) then
               nqx2=ix2(J)
               nqy2=iy2(J)
               nqz2=iz2(J)
               nx=nqx2-nqx1+1
               ny=nqy2-nqy1+1
               nz=nqz2-nqz1+1
               allocate(quantity(nqx1:nqx2,nqy1:nqy2,nqz1:nqz2,0:time_dim),stat=allerr)
               allocate(cqx(nqx1:nqx2))
               allocate(cqy(nqy1:nqy2))
               allocate(cqz(nqz1:nqz2))
               if (allerr .eq. 0) then
                  print * 
                  write(*,"(A,3(I3,':',I3),I3,':',I6)") " Quantity allocated by dimension: ", nqx1, nqx2, nqy1,nqy2,nqz1,nqz2,0,time_dim
               end if
            end if
         end do

      case(2)
         do J=1,n_auto_slices
            I=auto_slice_lists(J)
            ix1(J)=int((x1(I)-ixoff)/Dspace)
            ix2(J)=int((x2(I)-ixoff)/Dspace)
            do kk=1,niy
               if (abs(ytemp(kk)-y1(I)) .lt. eps) then
                  iy1(J)=kk;
                  iy2(J)=kk;
                  exit
               end if
            end do
            iz1(J)=int((z1(I)-izoff)/Dspace)
            iz2(J)=int((z2(I)-izoff)/Dspace)
            write(6,'(A,I3,A,I3,5(a,f8.3))')"Seq:",J,", Mesh:",slcf_mesh(I),", xmin=",x1(I),", xmax=",x2(I),", zmin=",z1(I),", zmax=",z2(I),", y=",y1(I)
            write(6,'(3(3X,A,I3,2X,I3))') "The mesh size on each axis is: IX =",ix1(j),ix2(j),", IY=",iy1(j),iy2(j),", IZ=",iz1(j),iz2(J)
            if (J .eq. 1) then
               nqx1=ix1(J)
               nqy1=iy1(J)
               nqz1=iz1(J)
            end if
            if (J .eq. n_auto_slices) then
               nqx2=ix2(J)
               nqy2=iy2(J)
               nqz2=iz2(J)
               nx=nqx2-nqx1+1
               ny=nqy2-nqy1+1
               nz=nqz2-nqz1+1
               allocate(quantity(nqx1:nqx2,nqy1:nqy2,nqz1:nqz2,0:time_dim),stat=allerr)
               allocate(cqx(nqx1:nqx2))
               allocate(cqy(nqy1:nqy2))
               allocate(cqz(nqz1:nqz2))
               if (allerr .eq. 0) then
                  print * 
                  write(*,"(A,3(I3,':',I3),I3,':',I6)") " Quantity allocated by dimension: ", nqx1, nqx2, nqy1,nqy2,nqz1,nqz2,0,time_dim
               end if
            end if
         end do

      case(3)

         do J=1,n_auto_slices
            I=auto_slice_lists(J)
            ix1(J)=int((x1(I)-ixoff)/Dspace)
            ix2(J)=int((x2(I)-ixoff)/Dspace)
            iy1(J)=int((y1(I)-iyoff)/Dspace)
            iy2(J)=int((y2(I)-iyoff)/Dspace)
            do kk=1,niz
               if (abs(ztemp(kk)-z1(I)) .lt. eps) then
                  iz1(J)=kk;
                  iz2(J)=kk;
                  exit
               end if
            end do
            write(6,'(A,I3,A,I3,5(a,f8.3))')"Seq:",J,", Mesh:",slcf_mesh(I),", xmin=",x1(I),", xmax=",x2(I),", ymin=",y1(I),", ymax=",y2(I),", z=",z1(I)
            write(6,'(3(3X,A,I3,2X,I3))') "The mesh size on each axis is: IX =",ix1(j),ix2(j),", IY=",iy1(j),iy2(j),", IZ=",iz1(j),iz2(J)
            if (J .eq. 1) then
               nqx1=ix1(J)
               nqy1=iy1(J)
               nqz1=iz1(J)
            end if
            if (J .eq. n_auto_slices) then
               nqx2=ix2(J)
               nqy2=iy2(J)
               nqz2=iz2(J)
               nx=nqx2-nqx1+1
               ny=nqy2-nqy1+1
               nz=nqz2-nqz1+1
               allocate(quantity(nqx1:nqx2,nqy1:nqy2,nqz1:nqz2,0:time_dim),stat=allerr)
               allocate(cqx(nqx1:nqx2))
               allocate(cqy(nqy1:nqy2))
               allocate(cqz(nqz1:nqz2))
               if (allerr .eq. 0) then
                  print * 
                  write(*,"(A,3(I3,':',I3),I3,':',I6)") " Quantity allocated by dimension: ", nqx1, nqx2, nqy1,nqy2,nqz1,nqz2,0,time_dim
               end if
            end if
         end do

      end select
   else
      ixoff=0
      iyoff=0
      izoff=0
   end if


   !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   !
   ! LOOP:  Read through N_Auto_Slices and write to output file
   !   	1. Read and combine all meshes of specific quantity into one piece
   !		2. Write the specific quantity into file
   !
   !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   ! Task 1: Input Dimension of IX,IY,IZ


   allocate(time(time_dim))

   Auto_List: DO JJ=1, n_auto_slices

      II=AUTO_SLICE_LISTS(JJ) ! I indicate file ID in slcf_file

      nm=slcf_mesh(II)
      m=>mesh(nm)
      allocate(f(0:m%ibar,0:m%jbar,0:m%kbar),stat=allerr)
      if (allerr .ne. 0) then
         write(*,*) " Error! variable f not allocated "
         write(*,*) "At iteration:", JJ
         stop
      end if

      f = 0.

      QFILE = SLCF_FILE(II)

      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
      READ(12)
      READ(12)
      READ(12)
      READ(12) I1,I2,J1,J2,K1,K2                     
      select case(stationID)
      case(1) ! Y-Z
         cqx(ix1(JJ):ix2(JJ))=x1(II)
         cqy(iy1(JJ):iy2(JJ))=m%y(J1:j2)
         cqz(iz1(JJ):iz2(JJ))=m%z(k1:k2)
      case(2) ! X-Z
         cqx(ix1(JJ):ix2(JJ))=m%x(I1:I2)
         cqy(iy1(JJ):iy2(JJ))=y1(II)
         cqz(iz1(JJ):iz2(JJ))=m%z(k1:k2)
      case(3) ! X-Y
         cqx(ix1(JJ):ix2(JJ))=m%x(I1:I2)
         cqy(iy1(JJ):iy2(JJ))=m%y(J1:j2)
         cqz(iz1(JJ):iz2(JJ))=z1(II)
      end select

      ! Read each slcf file with specific variable
      ncount = 0
      t1=ncount

      Read_loop: do
         read(12,end=99) time(ncount)
         read(12,end=99) (((f(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
         if (ncount .eq. 1) then
            print *
            write(*,*) " Reading file :", trim(qfile)
            write(*,'(3X,A,3(3X,I3,A,I3))') "Quantity index:",ix1(JJ),":",ix2(JJ),iy1(JJ),":",iy2(JJ),iz1(JJ),":",iz2(JJ)
            write(*,'(3X,A,3(3X,I3,A,I3))') "f index:",i1,":",i2,j1,":",j2,k1,":",k2
         end if

         quantity(ix1(JJ):ix2(JJ),iy1(JJ):iy2(JJ),iz1(JJ):iz2(JJ),ncount) = f(i1:i2,j1:j2,k1:k2)

         if (ncount .ge. time_dim) then
            print *
            write(*,'(A,I4,A)') " *** Fatal error: default time_dim = ",time_dim," is smaller than needed."
            write(*,'(A,f7.2)') " *** Try enlarge the parameter time_dim. Current time is: ",time(ncount)
            write(*,'(A,I4)') " *** Current iteration reaches: ",ncount
            stop
         end if
         ncount = ncount + 1
         t2=ncount
      end do Read_loop 
      99  close(12)
      deallocate(f)

   end do Auto_list
   nt=t2-t1+1

! #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
!
!                    Write individual quantity into .mat file 
!                       e.g. VELOCITY on all meshes
!
! #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#



   ! Output mat file name and path
   print *
   write(*,*) '##############################################################'
   write(*,*) '#'
   write(*,*) '#                Begin to write matfile'
   write(*,*) '#'
   write(*,*) '##############################################################'
   print *
   write(*,*) "Enter the output file name for current folder or full path :"
   read(*,'(A)') matfilename


   ! Reshape 4-D quan to 1-d quanline by colume (Fortran style)
   write(*,*) 'Reformating quantity into one line variables following the reshape rules...'
   allocate(quanline(1:size(quantity)))
   ii=1
   do l=t1,t2
      do k=nqz1,nqz2
         do j=nqy1,nqy2
            do i=nqx1,nqx2
               quanline(ii)=quantity(i,j,k,l)
               ii=ii+1
            end do
         end do
      end do
   end do

   allocate(xcor(1:nx))
   allocate(ycor(1:ny))
   allocate(zcor(1:nz))

   xcor(1:nx)=cqx
   ycor(1:ny)=cqy
   zcor(1:nz)=cqz

   !
   ! ===========================================================================
   ! *                      Write values to .mat file
   ! ============================================================================
   !

   ! Open MAT-file for writing
   write(*,*) "...Creating MAT-file: ",trim(matfilename)
   mpfile = matOpen(matfilename, 'u')
   if (mpfile == 0) then
      print*
      write(*,*) "Can't open: ",trim(matfilename)
      write(*,*) "  file may not exist, Try creating one with write access (v7.3 default)"
      mpfile = matOpen(matfilename, 'w7.3')
      if (mpfile == 0 ) then
         print*
         write(*,*) "Can''t open ",trim(matfilename), "directory may not accessable"
         write(*,*) "Program will exit"
         stop
      end if
   end if

   ! ---------- Data check -------------
   ! Dimension output
   print *
   write(*,*) "Dimensions  is:"
   frmt='(3X,A,3X,I4)'
   write(*,frmt) "X-axis dimension:",nx
   write(*,frmt) "Y-axis dimension:",ny
   write(*,frmt) "Z-axis dimension:",nz
   write(*,frmt) "Time dimension:",nt


   ! Prepare dimension for MATLAB variable creation 
   classid1=mxclassidFromclassname(realname)
   classid2=mxclassidFromclassname(intname)
   complexflag=0

   ! define some dimensions
   mdims(1)=nx
   mdims(2)=ny
   mdims(3)=nz
   mdims(4)=nt
   mlength=nx*ny*nz*nt

   dims(1)=nx
   dims(2)=ny
   dims(3)=nz
   dims(4)=nt

   mtlen=nt
   mdlen=4

   mxlen=nx
   mylen=ny
   mzlen=nz


   ! Check Time field
   ptime = matGetVariable(mpfile, 'time')
   if (ptime == 0 ) then
      print*
      write(*,*) "'Time' field doesn't exist, writing..."
      lwtime = .true.
   else 
      write(*,*) "'Time' field exists, skip..."
      lwtime = .false.
   endif

   ! Write time field
   if (lwtime) then
      ptime = mxCreateNumericArray(mtdim,mtlen,classid1,complexflag)
      call mxCopyReal4ToPtr(time, mxGetData(ptime), mtlen)
      mstat = matPutVariable(mpfile, 'time', ptime)  ! Variable name must not have '-', '='
      write(*,*) "Writing NumericArray 'time'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable 'time'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable 'time' succeeded"
   endif
   

   ! Check dimension field
   pdims = matGetVariable(mpfile, 'dims')
   if (pdims == 0 ) then
      print*
      write(*,*) "'Dims' field doesn't exist, writing..."
      lwdims = .true.
   else 
      write(*,*) "'Dims' field exists, skip..."
      lwdims = .false.
   endif

   ! Write dimension fields
   if (lwdims) then
      pdims = mxCreateNumericArray(mddim,mdlen,classid2,complexflag)
      call mxCopyInteger4ToPtr(dims, mxGetData(pdims), mdlen)
      mstat = matPutVariable(mpfile, 'dims', pdims)  ! Variable name must not have '-', '='
      write(*,*) "Writing NumericArray 'dims'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable 'dims'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable 'dims' succeeded"
   endif



   ! Check x-axis field
   pxor = matGetVariable(mpfile, 'xcor')
   if (pxor == 0 ) then
      print*
      write(*,*) "'xcor' field doesn't exist, writing..."
      lwxor = .true.
   else 
      write(*,*) "'xcor' field exists, skip..."
      lwxor = .false.
   endif

   ! Write X-coordinate fields
   if (lwxor) then
      pxor = mxCreateNumericArray(mxdim, mxlen, classid1, complexflag)
      call mxCopyReal4ToPtr(xcor, mxGetData(pxor), mxlen)
      mstat = matPutVariable(mpfile, 'xcor', pxor)
      write(*,*) "Writing NumericArray 'xcor'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable 'xcor'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable 'xcor' succeeded"
   endif



   ! Check Y-Axis field
   pyor = matGetVariable(mpfile, 'ycor')
   if (pyor == 0 ) then
      print*
      write(*,*) "'ycor' field doesn't exist, writing..."
      lwyor = .true.
   else 
      write(*,*) "'ycor' field exists, skip..."
      lwyor = .false.
   endif

   ! Write Y-Axis field
   if (lwyor) then
      pyor = mxCreateNumericArray(mydim, mylen, classid1, complexflag)
      call mxCopyReal4ToPtr(ycor, mxGetData(pyor), mylen)
      mstat = matPutVariable(mpfile, 'ycor', pyor)
      write(*,*) "Writing NumericArray 'ycor'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable 'ycor'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable 'ycor' succeeded"
   endif



   ! Check Z-Axis field
   pzor = matGetVariable(mpfile, 'zcor')
   if (pzor == 0 ) then
      print*
      write(*,*) "'zcor' field doesn't exist, writing..."
      lwzor = .true.
   else 
      write(*,*) "'zcor' field exists, skip..."
      lwzor = .false.
   endif

   ! Write Z-Axis field
   if (lwzor) then
      pzor = mxCreateNumericArray(mzdim, mzlen, classid1, complexflag)
      call mxCopyReal4ToPtr(zcor, mxGetData(pzor), mzlen)
      mstat = matPutVariable(mpfile, 'zcor', pzor)
      write(*,*) "Writing NumericArray 'zcor'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable 'zcor'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable 'zcor' succeeded"
   endif


   ! Check Main quantity field
   call namecheck(auto_slice_label, regularname)  ! varname check 

   pquan = matGetVariable(mpfile, trim(regularname))
   if (pquan == 0 ) then
      print *
      write(*,*) "'",trim(regularname),"' field doesn't exist, writing ..."
      lwquan = .true.
   else
      write(*,*) "'",trim(regularname),"' field exists, skip ..."
      lwquan = .false.
   end if


   ! Write Write Quantity field
   if (lwquan) then
      pquan = mxCreateNumericArray(mqdim,mdims,classid1,complexflag)
      call mxCopyReal4ToPtr(quanline, mxGetData(pquan), mlength)
      mstat = matPutVariable(mpfile, trim(regularname), pquan)  
      write(*,*) "Writing NumericArray '",trim(regularname),"'"
      if (mstat /= 0) then
         write(*,*) "***Error in matPutVariable '",trim(regularname),"'"
         goto 1000 ! Stop the program
      endif
      write(*,*)  "...Action matPutVariable '",trim(regularname), "' succeeded"
   end if



   ! Finish writing 
   1000 print*
   write(*,*) "==== Closing MAT-file, Program Ends ===="
   print*
   mstat = matClose(mpfile)
   if (mstat /= 0) then
      write(*,*) "***Error closing MAT-file, Program Stop"
      print*
      stop
   endif
   stop


End Program SLCF2MAT



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
!                      SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! *********************** search *******************************

subroutine search(string,length,lu,ierr)

   implicit none
   character(*), intent(in) :: string
   integer, intent(out) :: ierr
   integer, intent(in) :: lu, length
   character(20) :: junk

   search_loop: do 
      read(lu,'(a)',end=10) junk
      if (junk(1:length).eq.string(1:length)) exit search_loop
   end do search_loop

   ierr = 0
   return

   10 ierr = 1
   return

end subroutine search

subroutine search2(string,length,string2,length2,lu,ierr,choice)

   implicit none
   character(*), intent(in) :: string,string2
   character(*), intent(out) :: choice
   integer, intent(out) :: ierr
   integer, intent(in) :: lu, length, length2
   character(20) :: junk

   search_loop: do 
      read(lu,'(a)',end=10) junk
      if (junk(1:length).eq.string(1:length).or.junk(1:length2).eq.string2(1:length2)) then
         if (junk(1:length) .eq.string(1:length))   choice = junk(1:length)
         if (junk(1:length2).eq.string2(1:length2)) choice = junk(1:length2)
         exit search_loop
      end if
   end do search_loop

   ierr = 0
   return

   10 ierr = 1
   return

end subroutine search2

! *********************** parse *******************************

subroutine parse(buffer,sb_toks,se_toks,n_toks)

   ! parse buffer into a series of tokens each separated by blanks

   implicit none
   character(*), intent(inout) :: buffer
   integer, dimension(*), intent(out) :: sb_toks, se_toks
   integer, intent(out) :: n_toks
   integer :: i, intok, inquote, lenbuf
   character(len=1) :: c

   n_toks=0
   lenbuf = len(trim(buffer))
   if(lenbuf.eq.0)return ! buffer is blank so there are no tokens
   intok=0
   inquote=0
   do i = 1, lenbuf
      if(intok.eq.0)then  
         if(buffer(i:i).ne.' ')then  ! beginning of a new token since previous char
            intok=1                   ! was not in a token and this one is
            n_toks=n_toks + 1
            sb_toks(n_toks)=i
            if(buffer(i:i).eq."'")inquote=1
         end if
      else
         if(inquote.eq.1)then
            if(buffer(i:i).eq."'")then
               se_toks(n_toks)=i
               intok=0
               inquote=0
            end if
         end if
         if(buffer(i:i).eq.' ')then
            se_toks(n_toks)=i-1       ! previous char was in a token, this one is not
            intok=0                   ! so previous char is end of token
         end if
      end if
   end do
   if(buffer(lenbuf:lenbuf).ne.' ')se_toks(n_toks)=lenbuf ! last char in buffer is not blank
   ! so it is end of last token

   ! strip out single or double quotes if present
   do i = 1, n_toks
      c = buffer(sb_toks(i):sb_toks(i))
      if(c.eq."'")sb_toks(i)=sb_toks(i)+1
      c = buffer(se_toks(i):se_toks(i))
      if(c.eq."'")se_toks(i)=se_toks(i)-1
      if(se_toks(i).lt.sb_toks(i))then
         se_toks(i)=sb_toks(i)
         buffer(sb_toks(i):se_toks(i))=' '
      end if
   end do
end subroutine parse

! *********************** toupper *******************************

subroutine toupper(bufferin, bufferout)
   character(len=*), intent(in) :: bufferin
   character(len=*), intent(out) :: bufferout
   character(len=1) :: c

   integer :: lenbuf, i

   lenbuf=min(len(trim(bufferin)),len(bufferout))
   do i = 1, lenbuf
      c = bufferin(i:i)
      if(c.ge.'a'.and.c.le.'z')c=char(ichar(c)+ichar('a')-ichar('a'))
      bufferout(i:i)=c
   end do

end subroutine toupper

! ********************* nospace ********************************

subroutine nospace(str_in, str_out)
   character(256) :: str_in, str_out, temp
   character(1) :: space = ' '
   integer :: ilen,i
   logical :: iffind

   temp=adjustl(str_in)
   ilen=len_trim(temp)
   do  i=1,ilen
      if(temp(i:i) .eq. ' ') then
         temp(i:i)='-'
      end if
   end do
   str_out=temp
   return 

end subroutine nospace
! This is a Fortran program that read FDS 1st combined binary data 
! and write them into MATLAB ".mat" format file. 
!
#include "fintrf.h"





!******************** Namecheck *******************************************

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

