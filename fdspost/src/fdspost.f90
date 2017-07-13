program fdspost

  ! program to convert various fds output files to ascii

  implicit none

  interface
    subroutine parse(buffer,sb_toks,se_toks,n_toks)
      implicit none
      character(*), intent(inout) :: buffer
      integer, dimension(*), intent(out) :: sb_toks, se_toks
      integer, intent(out) :: n_toks
    end subroutine parse
  end interface

  character(255), parameter :: f2aversion='2.1.0'
  integer, parameter :: fb = selected_real_kind(6)
  integer, parameter :: file_dim = 500
  integer, parameter :: time_dim = 1010
  integer :: ierr, NMESHES, nm, noc, i, j, k, l
  integer :: idum, ifile, nsam, nv, mv
  integer :: i1, i2, j1, j2, k1, k2, t1, t2, i3, j3, k3  ! dimension
  integer :: i10, i20, j10, j20, k10, k20
  integer :: ncount, ior_input, npatch, ijbar, jkbar
  integer :: ii, nxp, nyp, nzp, n, i_sample
  real(fb) :: xs, xf, ys, yf, zs, zf  ! user defined bounds
  real(fb) :: d1, d2, d3, d4, tbeg, tend
  character(256) :: auto_slice_label
  integer :: auto_slice_flag, n_auto_slices, iauto, iztemp
  integer, dimension(1) :: izmin1
  integer :: izmin
  integer, dimension(file_dim) :: auto_slice_lists
  real(fb), dimension(file_dim) :: auto_slice_z
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
  integer, allocatable, dimension(:) :: ior,i1b,i2b,j1b,j2b,k1b,k2b
  integer, allocatable, dimension(:) :: auto_slices
  real(fb), allocatable, dimension(:,:,:,:,:) :: quantity
  real(fb), allocatable, dimension(:,:,:,:) :: q
  real(fb), allocatable, dimension(:,:,:) :: f
  real(fb), allocatable, dimension(:) :: time
  logical, allocatable, dimension(:,:,:) :: already_used
  logical :: new_plot3d=.true.
  integer ior_slcf
  character(1) :: ans
  character(4) :: ext1, ext2
  character(4) :: choice
  character(30) :: unitjunk
  character(256) gridfile,qname,chid,qfile,junk,frmt,outfile,outfile1,outfile2,slcf_label_dummy
  character(256), dimension(file_dim) :: pl3d_file,slcf_file,bndf_file,slcf_text,bndf_text,slcf_unit,bndf_unit,slcf_label
  character(256), dimension(file_dim) :: slice_labels
  character(20), dimension(file_dim) :: bndf_type
  character(20) :: bndf_type_chosen
  integer :: nslice_labels, slice_exist
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

  ! For debug usage
  !    write(*,*) "Selected_real_kind(6) = ",fb
  ! set a few default values

  zoffset=0.0
  zoffset_flag=0
  eps=0.001
  filein='stdin'

  ! parse command line arguments

  if (filein.eq.'stdin') then
    lu_in=5
  else
    lu_in=7
    inquire(file=filein,exist=exists)
    if (.not.exists) then
      write(6,*)"*** fatal Error: the file: ",trim(filein), " does not exist"
      stop
    endif
    open(lu_in,file=filein,status='old',form='formatted')
  endif

  write(6,*) ' Enter job id string (chid):'

  read(lu_in,'(a)') chid

  ! check to see if the .smv file exists

  gridfile = trim(chid)//'.smv'
  inquire(file=gridfile,exist=exists)
  if (.not.exists) then
    write(6,*)"*** fatal Error: the file: ",trim(gridfile)," does not exist"
    stop
  endif

  ! open the .smv file

  open(11,file=gridfile,status='old',form='formatted')
  write(*,*) "Open SMV file: ",trim(gridfile)
  print *


  ! determine the number of meshes

  rewind(11)

  call search('NMESHES',7,11,ierr)
  if (ierr.eq.1) then
    write(6,*) ' warning: assuming 1 mesh'
    NMESHES = 1
  else
    read(11,*) NMESHES
    write(*,*) "Number of mesh is: ",NMESHES
  endif

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
  enddo
  do i=0,m%ibar
  read(11,*) idum,m%x(i)
  enddo

  call search('TRNY',4,11,ierr)
  read(11,*) noc
  do i=1,noc
  read(11,*)
  enddo
  do j=0,m%jbar
  read(11,*) idum,m%y(j)
  enddo

  call search('TRNZ',4,11,ierr)
  read(11,*) noc
  do i=1,noc
  read(11,*)
  enddo
  do k=0,m%kbar
  read(11,*) idum,m%z(k)
  enddo

  enddo read_smv

  ! get sampling factor for data

  !    write(6,*) ' enter sampling factor for data?'
  !    write(6,*) ' (1 for all data, 2 for every other point, etc.)'

  !   read(lu_in,*) nsam
  nsam=1

  ! determine whether to limit domain size

  ! ans(1:1) may be 'y', 'n' or 'z' (or upper case equivalents)
  ! ans(2:2) may be 'a' or ' '

  ! y - domain size is limited
  ! n - domain size is not limited
  ! z - domain size is not limited z levels are offset by zoffset
  ! a - slice files are selected based on slice file type and location

  write(6,*) ' domain selection:'
  write(6,*) '   y - domain size is limited'
  write(6,*) '   n - domain size is not limited'
  write(6,*) '   z - domain size is not limited and z levels are offset'
  read(lu_in,'(a)') ans
  print *
  call toupper(ans,ans)
  if (ans(1:1).eq.'y') then
    write(6,*) ' Enter min/max x, y and z'
    read(lu_in,*) xs,xf,ys,yf,zs,zf
  else
    xs = -100000.
    xf =  100000.
    ys = -100000.
    yf =  100000.
    zs = -100000.
    zf =  100000.
  endif

  ! if ans is z or z then subtract zoffset from z data (for multi-level cases)

  if (ans(1:1).eq.'z') then
    zoffset_flag=1
  else
    zoffset_flag=0
  endif

  ! extract slcf data

  slcf_label = 'null'

  !    if (zoffset_flag.eq.0) then
  !        write(6,*) ' Enter starting and ending time for averaging (s)'
  !        read(lu_in,*) tbeg,tend
  !    else
  !        write(6,*) ' Enter starting and ending time for averaging (s) and zoffset (m)'
  !        read(lu_in,*) tbeg,tend,zoffset
  !    endif

  slcf_mesh = 1

  rewind(11)
  nfiles_exist=0

  ! ------------- key ieteration search slcf (i)
  search_slcf: do i=1,file_dim    !====================================================

  call search2('SLCF',4,'SLCC',4,11,ierr,choice)
  if (ierr.eq.1) exit search_slcf
  backspace(11)
  read(11,*) junk,slcf_mesh(i)  ! Slcf dimensions
  read(11,'(a)') slcf_file(i)  ! slcf file name
  read(11,'(a)') slcf_text(i)  ! slcf variable defination

  ! create unique list of slice types      

  nslice_labels=0
  slice_exist=0

  do ii=1, nslice_labels ! loop over labels kind
  if(trim(slcf_text(i)).eq.slice_labels(ii))then
    slice_exist=1 ! set flag slice already recorded
    exit  
  endif
  enddo

  if (slice_exist.eq.0) then ! new slice
    nslice_labels=nslice_labels+1  ! add record 
    slice_labels(nslice_labels)=trim(slcf_text(i)) ! Record new slice variable
  endif


  read(11,*) 
  read(11,'(a)') slcf_unit(i)
  open(12,file=trim(slcf_file(i)),form='unformatted',status='old', iostat=rcode)
  if (rcode.ne.0) then
    close(12)
    cycle ! jump to search_slcf for another run
  endif

  nfiles_exist=nfiles_exist+1

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

  write(6,'(i3,1x,a,1x,a)')i,trim(slcf_text(i)),trim(slcf_file(i))
  write(6,'(3x,a,6(1x,f8.2))')'slice bounds:',m%x(i1),m%x(i2),m%y(j1),m%y(j2),m%z(k1),m%z(k2)
  write(6,*) 

  enddo search_slcf

  if (nfiles_exist.eq.0)then
    write(6,*)"there are no slice files to convert"
    stop
  endif


  write(6,*)'How many variables to read:'
  read(lu_in,*) nv   ! ---- key iteration mv = 1:nv
  print *
  n_auto_slices=1


  sum = 0.

  varloop: do mv=1,nv  ! -----------------------------------

  slcf_label_dummy=' '
  write(6,'(a,i2)') ' Enter index for variable',mv
  read(lu_in,'(a)',iostat=Error_status) buffer
  print *

  if(Error_status.ne.0)then
    write(6,*)"*** fatal Error: read of variable index failed"
    stop
  endif

  call parse(buffer,sb_toks,se_toks,n_toks)

  if(n_toks.ge.1)then
    read(buffer(sb_toks(1):se_toks(1)),*)i
    if(n_toks.gt.1)slcf_label_dummy=buffer(sb_toks(2):se_toks(n_toks))
  else
    write(6,*)"*** fatal Error: index for variable ",mv," not entered"
    stop
  endif

  slcf_label(i) = slcf_label_dummy

  if (slcf_label(i)=='null' .or. slcf_label(i)==' ') slcf_label(i) = slcf_text(i)

  is(mv) = i
  qfile = slcf_file(i)

  if (mv.eq.1) then
    nm = slcf_mesh(i)
    m=>mesh(nm)
    allocate(q(0:m%ibar,0:m%jbar,0:m%kbar,nv))
    allocate(f(0:m%ibar,0:m%jbar,0:m%kbar))
    write(*,*) "This shows Q and F are successfully allocated"
    allocate(quantity(0:m%ibar,0:m%jbar,0:m%kbar,nv,time_dim))
    write(*,*) "To this point, quantity is also allocated"
    f = 0.
    q = 0.
  else
    if (slcf_mesh(i).ne.nm) then
      write(6,*) ' Error: all slices must have the same mesh'
      stop
    endif
  endif

  open(12,file=qfile,form='unformatted',status='old')
  write(*,*) "Open SLCF file: ",trim(qfile)," @line 357"
  read(12)
  read(12)
  read(12)
  read(12) i1,i2,j1,j2,k1,k2                    

  if (mv.eq.1) then
    i10=i1 ; i20=i2 ; j10=j1 ; j20=j2 ; k10=k1 ; k20=k2
    if (i1.eq.i2) ior_slcf = 1
    if (j1.eq.j2) ior_slcf = 2
    if (k1.eq.k2) ior_slcf = 3
  else
    if (i1.eq.i2 .and. i10.eq.i20) then
      i1=i10
      i2=i20
    endif
    if (j1.eq.j2 .and. j10.eq.j20) then
      j1=j10
      j2=j20
    endif
    if (k1.eq.k2 .and. k10.eq.k20) then
      k1=k10
      k2=k20
    endif
    if (((i1.ne.i10.or.i2.ne.i20).and.(i10.ne.i20)) .or. &
      ((j1.ne.j10.or.j2.ne.j20).and.(j10.ne.j20)) .or. &
      ((k1.ne.k10.or.k2.ne.k20).and.(k10.ne.k20))) then
      write(6,*) ' Error: slice files are incompatible'
      stop
    endif
  endif

  ncount = 1
  t1=ncount
  allocate(time(time_dim))

  ! read each slcf file with specific variable

  Read_loop: do
  read(12,end=99) time(ncount)
  !   write(*,*) "Reading file: ",trim(qfile),"at simulation time: ",time(ncount)

  read(12,end=99) (((f(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
  !        if (time.lt.tbeg) cycle read_loop
  !        if (time.gt.tend) exit read_loop
  t2=ncount
  quantity(i1:i2,j1:j2,k1:k2,mv,ncount) = f(i1:i2,j1:j2,k1:k2)
  q(i1:i2,j1:j2,k1:k2,mv) = q(i1:i2,j1:j2,k1:k2,mv)+ f(i1:i2,j1:j2,k1:k2)
  if (ncount .ge. time_dim) then
    write(*,"(A,I4,A)") " *** Fatal error: preset time_dim = ",time_dim," is smaller than needed."
    write(*,"(A,f7.2)") " *** Try enlarge the parameter time_dim. Current time is: ",time(ncount)
    write(*,"(A,I4)") " *** Current iteration reaches: ",ncount
    stop
  endif
  ncount = ncount + 1
  enddo Read_loop

  99  close(12)
  print *
  write(*,*) "File ",trim(qfile)," has been read with T_end = ",t2

  enddo varloop   ! varloop end -------------------------------------------------------




  !================ write out sample data to an ascii file =================================

  i_sample=1000
  ext1='.csv'
  ext2='.nc '

  print *
  write(6,*) 'Enter output file prefix: (the extension .csv and .nc will be added automatically)'
  read(lu_in,'(a)') outfile
  outfile1=trim(outfile)//trim(ext1)
  outfile2=trim(outfile)//trim(ext2)
  open(44,file=outfile1,form='formatted',status='unknown')

  i3 = i2 - i1 + 1
  j3 = j2 - j1 + 1
  k3 = k2 - k1 + 1

  ! one-dimensional section file

  if (i1.eq.i2 .and. j1.eq.j2 .and. k1.ne.k2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv,"(a,','),a)"
    ! write header
    write(44,frmt) 'z',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv,"(e12.5,','),e12.5)"
    write(6,*) ' writing to file z-axis      ',trim(outfile1)
    ! write data
    loop1: do k=k1,k2,nsam
    if (m%z(k).gt.zf .or. m%z(k).lt.zs) cycle loop1
    if (zoffset_flag.eq.1.and.m%z(k)-zoffset.lt.0.0) cycle loop1
    write(44,frmt) m%z(k)-zoffset,(q(i2,j2,k,l),l=1,nv)
    enddo loop1
  endif

  if(i1.eq.i2.and.j1.ne.j2.and.k1.eq.k2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv,"(a,','),a)"
    write(44,frmt) 'y',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv,"(e12.5,','),e12.5)"
    write(6,*) ' writing to file y-axis      ',trim(outfile1)
    loop2: do j=j1,j2,nsam
    if (m%y(j).gt.yf .or. m%y(j).lt.ys) cycle loop2
    write(44,frmt) m%y(j),(q(i2,j,k2,l),l=1,nv)
    enddo loop2
  endif

  if(i1.ne.i2.and.j1.eq.j2.and.k1.eq.k2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv,"(a,','),a)"
    write(44,frmt) 'x',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv,"(e12.5,','),e12.5)"
    write(6,*) ' writing to file x-axis      ',trim(outfile1)
    loop3: do i=i1,i2,nsam
    if (m%x(i).gt.xf .or. m%x(i).lt.xs) cycle loop3
    write(44,frmt) m%x(i),(q(i,j2,k2,l),l=1,nv)
    enddo loop3
  endif

  ! two-dimensional section file

  if(i1.eq.i2.and.j1.ne.j2.and.k1.ne.k2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv+1,"(a,','),a)"
    write(44,frmt) 'y','z',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm','m',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv+1,"(e12.5,','),e12.5)"
    write(6,*) ' writing y-z surface data to file  ',trim(outfile1)
    do k=k1,k2,nsam
    loop4: do j=j1,j2,nsam
    if (m%y(j).gt.yf .or. m%y(j).lt.ys) cycle loop4
    if (m%z(k).gt.zf .or. m%z(k).lt.zs) cycle loop4
    write(44,frmt) m%y(j),m%z(k)-zoffset,(q(i2,j,k,l),l=1,nv)
    enddo loop4
    enddo
  endif

  if (j1.eq.j2.and.i1.ne.i2.and.k1.ne.k2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv+1,"(a,','),a)"
    write(44,frmt) 'x','z',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm','m',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv+1,"(e12.5,','),e12.5)"
    write(6,*) ' writing x-z surface data to file  ',trim(outfile1)
    do k=k1,k2,nsam
    loop5: do i=i1,i2,nsam
    if (m%x(i).gt.xf .or. m%x(i).lt.xs) cycle loop5
    if (m%z(k).gt.zf .or. m%z(k).lt.zs) cycle loop5
    write(44,frmt) m%x(i),m%z(k)-zoffset,(q(i,j2,k,l),l=1,nv)
    enddo loop5
    enddo
  endif

  if(k1.eq.k2.and.i1.ne.i2.and.j1.ne.j2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv+1,"(a,','),a)"
    write(44,frmt) 'x','y',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm','m',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv+1,"(e12.5,','),e12.5)"
    write(6,*) ' writing x-y surface data to file  ',trim(outfile1)
    do j=j1,j2,nsam
    loop6: do i=i1,i2,nsam
    if (m%x(i).gt.xf .or. m%x(i).lt.xs) cycle loop6
    if (m%y(j).gt.yf .or. m%y(j).lt.ys) cycle loop6
    write(44,frmt) m%x(i),m%y(j),(q(i,j,k2,l),l=1,nv)
    enddo loop6
    enddo
  endif

  ! three-dimensional section file

  if (i1.ne.i2.and.j1.ne.j2.and.k1.ne.k2) then
    write(frmt,'(a,i2.2,a)') "(1x,",nv+2,"(a,','),a)"
    write(44,frmt) 'x','y','z',(trim(slcf_label(is(l))),l=1,nv)
    write(44,frmt) 'm','m','m',(trim(slcf_unit(is(l))),l=1,nv)
    write(frmt,'(a,i2.2,a)') "(",nv+2,"(e12.5,','),e12.5)"
    write(6,*) ' writing 3-d body data to file ',trim(outfile1)
    do k=k1,k2,nsam
    do j=j1,j2,nsam
    loop7: do i=i1,i2,nsam
    if (m%x(i).gt.xf .or. m%x(i).lt.xs) cycle loop7
    if (m%y(j).gt.yf .or. m%y(j).lt.ys) cycle loop7
    if (m%z(k).gt.zf .or. m%z(k).lt.zs) cycle loop7
    write(44,frmt) m%x(i),m%y(j),m%z(k),(q(i,j,k,l),l=1,nv)
    enddo loop7
    enddo
    enddo
  endif

  close(44)

  stop
end program fdspost

! ################### Subroutine Part #############################

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
  enddo search_loop

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
  endif
  enddo search_loop

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
    endif
  else
    if(inquote.eq.1)then
      if(buffer(i:i).eq."'")then
        se_toks(n_toks)=i
        intok=0
        inquote=0
      endif
    endif
    if(buffer(i:i).eq.' ')then
      se_toks(n_toks)=i-1       ! previous char was in a token, this one is not
      intok=0                   ! so previous char is end of token
    endif
  endif
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
  endif
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
