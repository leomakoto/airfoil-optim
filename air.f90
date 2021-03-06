! ******************************************************************
! ******************************************************************

program algencanma

  implicit none

  ! LOCAL SCALARS
  logical      :: checkder
  integer      :: allocerr,hnnzmax,hnnzmax1,hnnzmax2,hnnzmax3, &
       inform,jcnnzmax,m,n,nvparam,i
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm

  ! LOCAL ARRAYS
  character(len=15)     :: strtmp
  character(len=80)     :: specfnm,outputfnm,vparam(10)
  logical               :: coded(11)
  logical,      pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
              myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Number of variables
  ! This has to be manually set (for now)
  n = 8
  
  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  ! Set bounds for variables N1, N2, Au, and Al. They must all be in [0,1].
  l(1) = 0.45d0
  l(2) = 0.95d0
  u(1) = 0.50d0
  u(2) = 1.00d0
  x(1) = 0.5d0
  x(2) = 1.0d0
  do i = 3, n 
     l(i) =   0.1d-3
     u(i) =   2.0d0
     x(i) =   0.2d0
  end do
  
  ! Read initial guess. Comment this if no initial guess is provided.
  open(42, file = 'eval/feval.in', status = 'old', action = 'read')
  read(42, *) x
  close(42)
  
  ! write(*,*) x

  ! Constraints

  m = 0

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  equatn(1:m) = .false.
  lambda(1:m) = 0.0d0

  linear(1) = .false.
  linear(2) = .true.

  ! Coded subroutines

  coded(1:1)  = .true.  ! fsub
  coded(2:11) = .false. ! gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub,gjacpsub,hlsub,hlpsub

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = 4
  hnnzmax1 = 0
  hnnzmax2 = 1
  hnnzmax3 = 6
  hnnzmax  = hnnzmax1 + hnnzmax2 + hnnzmax3

  ! Checking derivatives?

  checkder = .false.

  ! Parameters setting

  epsfeas   = 1.0d-08
  epsopt    = 1.0d-08

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  efacc     = sqrt( epsfeas )
  eoacc     = sqrt( epsopt )

  outputfnm = ''
  specfnm   = ''

  nvparam = 0

  ! Optimize

  call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
       checkder,f,cnorm,snorm,nlpsupn,inform)

  deallocate(x,l,u,lambda,equatn,linear,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

  stop

end program algencanma

! ******************************************************************
! ******************************************************************

subroutine myevalf(n,x,f,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  
  ! AUXILIARY
  real(kind=8) :: coeffs(2)
  integer :: i

  flag = 0
  
  ! Write input for python script
  open(42, file = 'eval/feval.in', status = 'old', action = 'write')
  
  write(42, '(F0.8)', advance='no') x(1)
  do i = 2, n
  	write(42, '(" ", F0.8)', advance='no') x(i)
  end do
  
  close(42)
  
  ! Run python script
  call execute_command_line ("python airfoil.py")
  
  ! Read output of python script
  open(43, file = 'eval/feval.out', status = 'old', action = 'read')
  read(43, *) coeffs

  close(43)
  
  ! coeffs(1) is drag coeff and coeff2(2) is lift coeff
  f=coeffs(1)-coeffs(2)*coeffs(2)*coeffs(2)
  

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  flag = 0

  g(1) = 0.0d0
  g(2) = 1.0d0

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,n
  integer,      intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  flag = 0
  lmem = .false.

  hnnz = 0

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: ind,n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)

  flag = 0

  if ( ind .eq. 1 ) then
     c = x(1) ** 2 + 1 - x(n)
     return
  end if

  if ( ind .eq. 2 ) then
     c = 2.0d0 - x(1) - x(n)
     return
  end if

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in)  :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  flag = 0
  lmem = .false.

  if ( ind .eq. 1 ) then
     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = 1
     jcval(1) = 2.0d0 * x(1)
     jcvar(2) = 2
     jcval(2) = - 1.0d0

     return
  end if

  if ( ind .eq. 2 ) then
     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = 1
     jcval(1) = - 1.0d0
     jcvar(2) = 2
     jcval(2) = - 1.0d0

     return
  end if

end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in)  :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  flag = 0
  lmem = .false.

  if ( ind .eq. 1 ) then
     hcnnz = 1

     if ( hcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     hcrow(1) = 1
     hccol(1) = 1
     hcval(1) = 2.0d0

     return
  end if

  if ( ind .eq. 2 ) then
     hcnnz = 0
     return
  end if

end subroutine myevalhc

! ******************************************************************
! ******************************************************************

subroutine myevalfc(n,x,f,m,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: m,n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: c(m)

  flag = - 1

end subroutine myevalfc

! ******************************************************************
! ******************************************************************

subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,m,n
  integer,      intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: jcfun(lim),jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n),jcval(lim)

  flag = - 1

end subroutine myevalgjac

! ******************************************************************
! ******************************************************************

subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,   intent(inout) :: gotj
  integer,   intent(in)    :: m,n
  integer,   intent(out)   :: flag
  character, intent(in)    :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)    :: x(n)
  real(kind=8), intent(inout) :: p(m),q(n)
  real(kind=8), intent(out)   :: g(n)

  flag = - 1

end subroutine myevalgjacp

! ******************************************************************
! ******************************************************************

subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,m,n
  integer,      intent(out) :: flag,hlnnz
  real(kind=8), intent(in)  :: sf

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hlcol(lim),hlrow(lim)
  real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
  real(kind=8), intent(out) :: hlval(lim)

  flag = - 1

end subroutine myevalhl

! ******************************************************************
! ******************************************************************

subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(inout) :: goth
  integer,      intent(in)    :: m,n
  integer,      intent(out)   :: flag
  real(kind=8), intent(in)    :: sf

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
  real(kind=8), intent(out) :: hp(n)

  flag = - 1

end subroutine myevalhlp

