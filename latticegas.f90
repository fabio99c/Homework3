module latticegas

  implicit none

  logical, allocatable, public :: lattice(:,:) 
  integer, allocatable, public :: x(:),y(:)
  !integer, allocatable, public :: x_t(:,:),y_t(:,:)
  integer, dimension(33), public :: seed
  double precision,public, allocatable :: dx(:),dy(:)
  !real,public :: p_rho

  integer, public  :: Nsteps,Np,L,i_x 
  integer, public  :: istep,isubstep  
  integer, public  :: dir,i,j,nfail,njumps, sizer
  !integer, dimension(:), allocatable, public :: seed

  integer, public :: free(4),nfree   
  integer, public :: dxtrial(4),dytrial(4) 
  integer, public :: xnew(4),ynew(4)

  real, dimension(2), public :: rnd(2)
  !real,allocatable, public :: p_c(:)
  real,public :: rnd1
  double  precision,public :: dxsum,dysum,dxsqsum,dysqsum 
  double  precision,public :: t,deltat,drsqave,D,a,help
  
  contains

  subroutine init()
    call random_seed(put=seed)
	
    allocate(lattice(0:L-1,0:L-1))
    allocate(x(Np),y(Np)) 
    allocate(dx(Np),dy(Np))
    !allocate(p_c(Np))
    !allocate(x_t(Np,Np*Nsteps+1))
    !allocate(y_t(Np,Np*Nsteps+1))
	
	
    dxtrial(1)=+1; dytrial(1)= 0;   
    dxtrial(2)=-1; dytrial(2)= 0;   
    dxtrial(3)= 0; dytrial(3)=+1; 
    dxtrial(4)= 0; dytrial(4)=-1;
	
    nfail=0; njumps=0;
	
    i_x=2
    !p_rho=0
    
	if (Np >= L*L) then
		print *,'Number of particles > number of sites !!!' 
		STOP 'Too small lattice' 
	endif
	

    do i=0,L-1
		do j=0,L-1
			lattice(i,j) = .false. 
		end do
    end do
	

    ! Generate particles on lattice 
    do i=1,Np
		do ! Loop until empty position found
        !  To   be  on  safe  side,   check  that  upper   limit  not  reached
			call random_number(rnd)
			x(i)=int(rnd(1)*L);  if (x(i)>=L) x(i)=L-1;  
			y(i)=int(rnd(2)*L);  if (y(i)>=L) y(i)=L-1; 
			if (lattice(x(i),y(i))) then
			! Position already filled, loop to find new trial 
				cycle 
			else
				lattice(x(i),y(i))=.true. 
            !  Success, go  to next particle  
                exit 
            endif
        enddo
        !x_t(i,1)=x(i)
        !y_t(i,1)=y(i)
        dx(i)=0.0d0; dy(i)=0.0d0; 
    enddo
   

  end subroutine init  

  subroutine metropolis()
	do isubstep=1,Np ! Do all particles on average once every MC step
        ! Pick one particle at random 
        call random_number(rnd1)
        i=int(rnd1*(Np-1))+1;
        ! Find possible directions, store it in free() 
        nfree=0 
        do j=1,4
           xnew(j)=x(i)+dxtrial(j); 
           if  (xnew(j) >= L) xnew(j)=0; if (xnew(j)<0) xnew(j)=L-1; 
           ynew(j)=y(i)+dytrial(j); 
           if  (ynew(j) >= L) ynew(j)=0; if (ynew(j)<0) ynew(j)=L-1; 
           if (.not. lattice(xnew(j),ynew(j))) then
              ! Success: position free 
              nfree=nfree+1 
              free(nfree)=j 
           endif
        end do
        ! If no possible directions, get new particle 
        if (nfree == 0) then
           nfail=nfail+1 
           cycle 
        endif
        njumps=njumps+1
        !  Pick  one of  the  possible directions  randomly  
        !  Note that  the dir>nfree  check here  really is  needed!  
        call random_number(rnd1)
        dir=int(rnd1*nfree)+1; if (dir>nfree) dir=nfree 
        j=free(dir)
        ! Now x(i),y(i) is old position and xnew(j),ynew(j) new 
        ! Double check that new site really is free 
        if (lattice(xnew(j),ynew(j))) then
           print *,'ERROR:    THIS   SHOULD   BE   IMPOSSIBLE'   
           print *,i,j,dir,nfree  
           print *,free  
           print *,x(i),y(i),xnew(j),ynew(j) 
           STOP 'ERROR  new  site  bug'  
        endif
        !Empty  old  position  and  fill  new
        lattice(x(i),y(i))=.false. 
        lattice(xnew(j),ynew(j))=.true.
        x(i)=xnew(j); y(i)=ynew(j);          
        dx(i)=dx(i)+dxtrial(j); dy(i)=dy(i)+dytrial(j);
		!do j=1,Np
		!	x_t(j,i_x)=x(j)
		!	y_t(j,i_x)=y(j)
		!end do
        i_x=i_x+1
		
    end do
  end subroutine Metropolis
  
  !subroutine probability()
  !integer  :: k,j,i 
  !real, dimension(Np)::p_,p
  !p_=0;p=0;
  !do k=1,Np
  !   do i=2,(Np*Nsteps)-1
  !         if (x_t(k,i)/=x_t(k,i-1) .or. y_t(k,i)/=y_t(k,i-1)) then
  !              p_(k)=p_(k)+1.0
  !              if (x_t(k,i+1)==x_t(k,i-1) .and. y_t(k,i+1)==y_t(k,i-1)) then
  !                   p(k)=p(k)+1.0
  !              endif
  !              p_c(k)=p(k)/p_(k)
  !         endif
  !   end do
  !end do
  !p_rho = sum(p_c)/real(Np)
  !end subroutine probability

  subroutine metropolis_driver()
  t=0.0; 
  do istep=0,Nsteps-1 ! Loop over MC steps
	call metropolis()
	if (mod(istep*Np,1) == 0) then
        ! Calculate  and print intermediate results
        ! Get total displacement from dx,dy 
		dxsum=0.0d0; dysum=0.0d0; 
        dxsqsum=0.0d0; dysqsum=0.0d0; 
        do i=1,Np
           dxsum=dxsum+dx(i);   dysum=dysum+dy(i);   
           dxsqsum=dxsqsum+dx(i)*dx(i);
           dysqsum=dysqsum+dy(i)*dy(i); 
        enddo
        drsqave=(dxsqsum+dysqsum)/(1.0*Np)

        if (t>0.0) then
           !  Get  diffusion coefficient  by  proper scaling  
           D=drsqave*a*a/(4*t)
           !write(*,fmt='(3(a,1pe10.2))')'Density Np/L^2=',real(Np)/L**2,' :  t=',t,';   D=',D
        endif

     endif

     t=t+deltat  

  end do
  !print *,'Number of  failed jumps',nfail,' number of  successes',njumps
  
  end subroutine metropolis_driver
  
  subroutine deallocate_()
      deallocate(lattice)
      deallocate(x,y) 
      deallocate(dx,dy)
      !deallocate(p_c)
      !deallocate(x_t,y_t)
  end subroutine deallocate_

end module latticegas
