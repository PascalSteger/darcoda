!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::i,j,id,iu,iv,iw,ip
  real(dp)::rho1,rho2,T1,T2,xc,yc,zc,rad,radcl
  integer::ivar,irad,imet
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::zgas1,zgas2

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   
  id=1; iu=2; iv=3; iw=4; ip=ndim+2; imet=ndim+3;irad=ndim+4

  rho1=gravity_params(1)
  T1=gravity_params(2)
  zgas1=gravity_params(3)
  rho2=gravity_params(4)
  T2=gravity_params(5)
  zgas2=gravity_params(6)
  radcl=gravity_params(7)

  xc=boxlen/2.0
  yc=boxlen/2.0
  zc=boxlen/2.0

  radcl=radcl*3.08568e18/scale_l !into code units

  do i=1,nn
     q(i,id)=rho1       ! already in H/cc
     q(i,iu)=0.0
     q(i,iv)=0.0
     q(i,iw)=0.0
     q(i,ip)=q(i,id)*T1/scale_T2
     q(i,imet)=zgas1*0.02d0 
  enddo
  do i=1,nn
     rad=((x(i,1)-xc)**2+(x(i,2)-yc)**2+(x(i,3)-zc)**2)**0.5
     if(rad<=radcl) then
        q(i,id)=rho2 !H/cc
        q(i,iu)=0.0
        q(i,iv)=0.0
        q(i,iw)=0.0 
        q(i,ip)=q(i,id)*T2/scale_T2 !Is this true?
        q(i,imet)=zgas2*0.02d0  !*0.01
        !              write(*,*) q(i,id), q(i,ip), rho1, rho1*T1/scale_T2,T1,T2 
     endif
  enddo
  
  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
