!> Basic Lagrangian particle solver class:
!> Provides support for Lagrangian-transported objects in liquid phase
!> Add add-mass, buoyancy, and lubrication forces
module lptl_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use lpt_class,     only: lpt
   implicit none
   private

   public :: lptl

   type, extends(lpt) :: lptl

   contains

      procedure :: add_buoyancy
      procedure :: add_pressureGrad
      procedure :: collide=>lptl_collide

   end type lptl

   interface lptl
      procedure construct_from_lpt
   end interface lptl

contains

   function construct_from_lpt(cfg, name) result(self)
      use lpt_class, only: lpt
      implicit none
      type(lptl) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      self%lpt=lpt(cfg,name)
   end function construct_from_lpt

   ! NOTE: This should be call after collision. Cuz collision will zero the acceleration
   subroutine add_buoyancy(this, rho)
      implicit none
      class(lptl), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho

      integer :: i
      integer :: id, jd, kd
      real(WP), dimension(3) :: g

      g = this%gravity

      buoyancy: do i=1, this%np_

         if (this%p(i)%id.le.0) cycle buoyancy

         id = this%p(i)%ind(1); jd = this%p(i)%ind(2); kd = this%p(i)%ind(3)

         this%p(i)%Acol = this%p(i)%Acol - g*rho(id,jd,kd)/this%rho

      end do buoyancy

   end subroutine add_buoyancy

   subroutine add_pressureGrad(this, pgradx, pgrady, pgradz)
      implicit none
      class(lptl), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: pgradx, pgrady, pgradz

      integer :: i
      integer :: id, jd, kd

      pressureGrad: do i=1, this%np_

         if (this%p(i)%id.le.0) cycle pressureGrad

         id = this%p(i)%ind(1); jd = this%p(i)%ind(2); kd = this%p(i)%ind(3)

         this%p(i)%Acol = this%p(i)%Acol - (/pgradx(id,jd,kd), pgrady(id,jd,kd), pgradz(id,jd,kd)/)/this%rho

      end do pressureGrad

   end subroutine add_pressureGrad

   !> Resolve collisional interaction between particles, walls, and an optional IB level set
   !> Requires tau_col, e_n, e_w and mu_f to be set beforehand
   subroutine lptl_collide(this,dt,Gib,Nxib,Nyib,Nzib)
      implicit none
      class(lptl), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Gib  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Nxib !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Nyib !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Nzib !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell

      ! Check if all IB parameters are present
      check_G: block
         use messager, only: die
         if (present(Gib).and.(.not.present(Nxib).or..not.present(Nyib).or..not.present(Nzib))) &
            call die('[lpt collide] IB collisions need Gib, Nxib, Nyib, AND Nzib')
      end block check_G

      ! Start by zeroing out the collision force
      zero_force: block
         integer :: i
         do i=1,this%np_
            this%p(i)%Acol=0.0_WP
            this%p(i)%Tcol=0.0_WP
         end do
      end block zero_force

      ! Then share particles across overlap
      call this%share()

      ! We can now assemble particle-in-cell information
      pic_prep: block
         use mpi_f08
         integer :: i,ip,jp,kp,ierr
         integer :: mymax_npic,max_npic

         ! Allocate number of particle in cell
         allocate(npic(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); npic=0

         ! Count particles and ghosts per cell
         do i=1,this%np_
            ip=this%p(i)%ind(1); jp=this%p(i)%ind(2); kp=this%p(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
         end do

         ! Get maximum number of particle in cell
         mymax_npic=maxval(npic); call MPI_ALLREDUCE(mymax_npic,max_npic,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)

         ! Allocate pic map
         allocate(ipic(1:max_npic,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); ipic=0

         ! Assemble pic map
         npic=0
         do i=1,this%np_
            ip=this%p(i)%ind(1); jp=this%p(i)%ind(2); kp=this%p(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=i
         end do
         do i=1,this%ng_
            ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
            npic(ip,jp,kp)=npic(ip,jp,kp)+1
            ipic(npic(ip,jp,kp),ip,jp,kp)=-i
         end do

      end block pic_prep

      ! Finally, calculate collision force
      collision_force: block
         use mpi_f08
         use mathtools, only: Pi,normalize,cross_product
         integer :: i1,i2,ii,jj,kk,nn,ierr
         real(WP) :: d1,m1,d2,m2,d12,m12,buf
         real(WP), dimension(3) :: r1,v1,w1,r2,v2,w2,v12,n12,f_n,t12,f_t
         real(WP) :: k_n,eta_n,k_coeff,eta_coeff,k_coeff_w,eta_coeff_w,rnv,r_influ,delta_n,rtv
         real(WP), parameter :: aclipnorm=1.0e-6_WP
         real(WP), parameter :: acliptan=1.0e-9_WP
         real(WP), parameter :: rcliptan=0.05_WP

         ! Reset collision counter
         this%ncol=0

         ! Precompute coefficients for k and eta
         k_coeff=(Pi**2+log(this%e_n)**2)/this%tau_col**2
         eta_coeff=-2.0_WP*log(this%e_n)/this%tau_col
         k_coeff_w=(Pi**2+log(this%e_w)**2)/this%tau_col**2
         eta_coeff_w=-2.0_WP*log(this%e_w)/this%tau_col

         ! Loop over all local particles
         collision: do i1=1,this%np_

            ! Cycle if id<=0
            if (this%p(i1)%id.le.0) cycle collision

            ! Store particle data
            r1=this%p(i1)%pos
            v1=this%p(i1)%vel
            w1=this%p(i1)%angVel
            d1=this%p(i1)%d
            m1=this%rho*Pi/6.0_WP*d1**3

            ! First collide with walls
            d12=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=this%Wdist,bc='d')
            n12=this%Wnorm(:,this%p(i1)%ind(1),this%p(i1)%ind(2),this%p(i1)%ind(3))
            n12=-normalize(n12+[epsilon(1.0_WP),epsilon(1.0_WP),epsilon(1.0_WP)])
            rnv=dot_product(v1,n12)
            r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
            delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)

            ! Assess if there is collision
            if (delta_n.gt.0.0_WP) then
               ! Normal collision
               k_n=m1*k_coeff_w
               eta_n=m1*eta_coeff_w
               f_n=-k_n*delta_n*n12-eta_n*rnv*n12
               ! Tangential collision
               f_t=0.0_WP
               if (this%mu_f.gt.0.0_WP) then
                  t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
                  rtv = sqrt(sum(t12*t12))
                  if (rnv*dt/d1.gt.aclipnorm) then
                     if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                  else
                     if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
                  end if
                  if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
               end if
               ! Calculate collision force
               f_n=f_n/m1; f_t=f_t/m1
               this%p(i1)%Acol=this%p(i1)%Acol+f_n+f_t
               ! Calculate collision torque
               this%p(i1)%Tcol=this%p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t)
            end if

            ! Collide with IB
            if (present(Gib)) then
               d12=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=Gib,bc='n')
               n12(1)=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=Nxib,bc='n')
               n12(2)=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=Nyib,bc='n')
               n12(3)=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=Nzib,bc='n')
               buf = sqrt(sum(n12*n12))+epsilon(1.0_WP)
               n12 = -n12/buf
               rnv=dot_product(v1,n12)
               r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
               delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)

               ! Assess if there is collision
               if (delta_n.gt.0.0_WP) then
                  ! Normal collision
                  k_n=m1*k_coeff_w
                  eta_n=m1*eta_coeff_w
                  f_n=-k_n*delta_n*n12-eta_n*rnv*n12
                  ! Tangential collision
                  f_t=0.0_WP
                  if (this%mu_f.gt.0.0_WP) then
                     t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
                     rtv = sqrt(sum(t12*t12))
                     if (rnv*dt/d1.gt.aclipnorm) then
                        if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                     else
                        if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
                     end if
                     if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
                  end if
                  ! Calculate collision force
                  f_n=f_n/m1; f_t=f_t/m1
                  this%p(i1)%Acol=this%p(i1)%Acol+f_n+f_t
                  ! Calculate collision torque
                  this%p(i1)%Tcol=this%p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t)
               end if
            end if

            ! Loop over nearest cells
            do kk=this%p(i1)%ind(3)-1,this%p(i1)%ind(3)+1
               do jj=this%p(i1)%ind(2)-1,this%p(i1)%ind(2)+1
                  do ii=this%p(i1)%ind(1)-1,this%p(i1)%ind(1)+1

                     ! Loop over particles in that cell
                     do nn=1,npic(ii,jj,kk)

                        ! Get index of neighbor particle
                        i2=ipic(nn,ii,jj,kk)

                        ! Get relevant data from correct storage
                        if (i2.gt.0) then
                           r2=this%p(i2)%pos
                           v2=this%p(i2)%vel
                           w2=this%p(i2)%angVel
                           d2=this%p(i2)%d
                           m2=this%rho*Pi/6.0_WP*d2**3
                        else if (i2.lt.0) then
                           i2=-i2
                           r2=this%g(i2)%pos
                           v2=this%g(i2)%vel
                           w2=this%g(i2)%angVel
                           d2=this%g(i2)%d
                           m2=this%rho*Pi/6.0_WP*d2**3
                        end if

                        ! Compute relative information
                        d12=norm2(r1-r2)
                        if (d12.lt.10.0_WP*epsilon(d12)) cycle !< this should skip auto-collision
                        n12=(r2-r1)/d12
                        v12=v1-v2
                        rnv=dot_product(v12,n12)
                        r_influ=min(abs(rnv)*dt,0.1_WP*(d1+d2))
                        delta_n=min(0.5_WP*(d1+d2)+r_influ-d12,this%clip_col*0.5_WP*(d1+d2))

                        ! Assess if there is collision
                        if (delta_n.gt.0.0_WP) then
                           ! Normal collision
                           m12=m1*m2/(m1+m2)
                           k_n=m12*k_coeff
                           eta_n=m12*eta_coeff
                           f_n=-k_n*delta_n*n12-eta_n*rnv*n12
                           ! Tangential collision
                           f_t=0.0_WP
                           if (this%mu_f.gt.0.0_WP) then
                              t12 = v12-rnv*n12+cross_product(0.5_WP*(d1*w1+d2*w2),n12)
                              rtv = sqrt(sum(t12*t12))
                              if (rnv*dt*2.0_WP/(d1+d2).gt.aclipnorm) then
                                 if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                              else
                                 if (rtv*dt*2.0_WP/(d1+d2).lt.acliptan) rtv=0.0_WP
                              end if
                              if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
                           end if
                           ! Calculate collision force
                           f_n=f_n/m1; f_t=f_t/m1
                           this%p(i1)%Acol=this%p(i1)%Acol+f_n+f_t
                           ! Calculate collision torque
                           this%p(i1)%Tcol=this%p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t)
                           ! Add up the collisions
                           this%ncol=this%ncol+1
                        end if

                     end do

                  end do
               end do
            end do

            ! Deal with dimensionality
            if (this%cfg%nx.eq.1) then
               this%p(i1)%Acol(1)=0.0_WP
               this%p(i1)%Tcol(2)=0.0_WP
               this%p(i1)%Tcol(3)=0.0_WP
            end if
            if (this%cfg%ny.eq.1) then
               this%p(i1)%Tcol(1)=0.0_WP
               this%p(i1)%Acol(2)=0.0_WP
               this%p(i1)%Tcol(3)=0.0_WP
            end if
            if (this%cfg%nz.eq.1) then
               this%p(i1)%Tcol(1)=0.0_WP
               this%p(i1)%Tcol(2)=0.0_WP
               this%p(i1)%Acol(3)=0.0_WP
            end if

         end do collision

         ! Determine total number of collisions
         call MPI_ALLREDUCE(this%ncol,nn,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%ncol=nn/2

      end block collision_force

      ! Clean up
      if (allocated(npic)) deallocate(npic)
      if (allocated(ipic)) deallocate(ipic)

   end subroutine lptl_collide


end module lptl_class
