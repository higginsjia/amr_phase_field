module advance_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(mla,phi,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer i,dm,n,nlevs
    dm = mla%dim
    nlevs = mla%nlevel
    
    call update_phi(mla,phi,dx,dt,the_bc_tower)


  end subroutine advance


  subroutine update_phi(mla,phi,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    ! local variables
    integer :: n, nlevs,i,lo(mla%dim), hi(mla%dim),ng_p
    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    nlevs = mla%nlevel
    ng_p = phi(1)%ng

    do n=1,nlevs
				do i=1,nfabs(phi(n))
				   pp  => dataptr(phi(n),i)
				   lo = lwb(get_box(phi(n),i))
				   hi = upb(get_box(phi(n),i))
           !! psi, psi_nex
				   call update_phi_2d(pp(:,:,1,1),pp(:,:,1,2), ng_p,lo, hi, dx(n), dt)
				end do
    end do

    do n=1,nlevs
				do i=1,nfabs(phi(n))
				   pp  => dataptr(phi(n),i)
				   lo = lwb(get_box(phi(n),i))
				   hi = upb(get_box(phi(n),i))
				   call array_copy_2d(pp(:,:,1,1),pp(:,:,1,2), ng_p,lo, hi)
				end do
    end do
    
    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array,2,1,1)


  end subroutine update_phi



  subroutine update_phi_2d(PSI,PSI_NEX, ng_p, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision :: PSI(lo(1)-ng_p:,lo(2)-ng_p:),PSI_NEX(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: dx, dt
    ! local variables
    integer i,j
    double precision ::cipjp,cipjm,cimjp,cimjm,JR,JL,JT,JB,DERX,DERY,Atheta,Aptheta
    double precision ::lamda,dx2_in,dpsi,itemFir,itemSec,itemTrd
    lamda=3.19
    dx2_in=1.0/dx/dx
    

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! perform fvm calculation here
					!update stencil value use nine-point scheme
					cipjp=( PSI(i+1,j+1) + PSI(i,j+1) + PSI(i,j) + PSI(i+1,j) )/4.0d0
					cipjm=( PSI(i+1,j) + PSI(i,j) + PSI(i,j-1) + PSI(i+1,j-1) )/4.0d0
					cimjp=( PSI(i,j+1) + PSI(i-1,j+1) + PSI(i-1,j) + PSI(i,j) )/4.0d0
					cimjm=( PSI(i,j) + PSI(i-1,j) + PSI(i-1,j-1) + PSI(i,j-1) )/4.0d0

					!! Calculate right edge flux
					DERX= PSI(i+1,j)-PSI(i,j)
					DERY= cipjp - cipjm
					call set_aniso_parameters(DERX,DERY,Atheta,Aptheta)
					JR = Atheta * ( Atheta*DERX + Aptheta*DERY )

					!! Calculate left edge flux
					DERX= PSI(i,j)-PSI(i-1,j)
					DERY= cimjp - cimjm
					call set_aniso_parameters(DERX,DERY,Atheta,Aptheta)
					JL = Atheta * ( Atheta*DERX + Aptheta*DERY )

					!! Calculate top edge flux
					DERX= cipjp - cimjp
					DERY= PSI(i,j+1)-PSI(i,j)
					call set_aniso_parameters(DERX,DERY,Atheta,Aptheta)
					JT = Atheta * ( Atheta*DERY - Aptheta*DERX )

					!! Calculate bottom edge flux
					DERX= cipjm - cimjm
					DERY= PSI(i,j)-PSI(i,j-1)
					call set_aniso_parameters(DERX,DERY,Atheta,Aptheta)
					JB = Atheta * ( Atheta*DERY - Aptheta*DERX )
          
          !! Calculate At for time step
					DERX= PSI(i+1,j)-PSI(i-1,j)
					DERY= PSI(i,j+1)-PSI(i,j-1)
					call set_aniso_parameters(DERX,DERY,Atheta,Aptheta)
          !! Calculate items prepared for equation
					itemFir=dx2_in * (JR - JL + JT - JB)
					itemSec=PSI(i,j)-PSI(i,j)**3
					itemTrd=lamda*(-0.1)*(1.0d0 - PSI(i,j)**2 )**2
          !! dPsi 
					dpsi = (dt / Atheta**2) *(itemFir+itemSec-itemTrd)
          !! update psi_nex
					PSI_NEX(i,j)=PSI(i,j)+dpsi
       end do
    end do

  end subroutine update_phi_2d



	subroutine set_aniso_parameters(DERX,DERY,Atheta,Aptheta)
		implicit none
		double precision ::DERX,DERY,Atheta,Aptheta
		double precision ::theta,theta_diff,epsion,eps,pi_my
		epsion=0.04
    eps=1.0d-8
		pi_my= 3.1415926535898

		if(DERX*DERX.gt.eps)then
		  theta=tanh(DERY/DERX)
		else if(DERY*DERY.lt.eps)then
		  theta=pi_my*0.25
		else
		  theta=pi_my*0.5
		end if
		theta_diff=theta-0.0
		!change orientation here
		Atheta=1.0+epsion*cos(4.0*theta_diff)
		Aptheta=4.0*epsion*sin(4.0*theta_diff)
	end subroutine set_aniso_parameters


  ! array copy
  subroutine array_copy_2d(dest,src, ng_p, lo, hi)

    integer          :: lo(2), hi(2), ng_p
    double precision ::   dest(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision ::   src(lo(1)-ng_p:,lo(2)-ng_p:)
    ! local variables
    integer i,j
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          dest(i,j) = src(i,j)         
       end do
    end do
  end subroutine array_copy_2d


end module advance_module

