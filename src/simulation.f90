module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use fft3d_class,       only: fft3d
   use lowmach_class,     only: lowmach
   use lpt_class,         only: lpt
   use rotorDisk_class,      only: rotorDisk
   use blade_class,          only: blade
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use sgsmodel_class,    only: sgsmodel
   use string, only: str_medium
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   type(fft3d),     public :: ps
   type(ddadi),     public :: vs
   type(lowmach),       public :: fs
   type(sgsmodel),    public :: sgs
   type(timetracker),   public :: time
   type(rotorDisk),      public :: rd
   type(blade),          public :: bl
   type(lpt),         public :: lp

   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   type(partmesh) :: pmesh

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,lptfile

   !> Work arrays and fluid properties
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resRho
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:), allocatable :: rho
   real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp

   real(WP) :: int_RP

   !> Viscosity
   real(WP) :: visc

   !> cfl
   real(WP) :: cfl, cflc, lp_cfl

contains

   subroutine simulation_init
      use param, only: param_read
      use string, only: str_medium
      implicit none

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resRho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(rho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(srcUlp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcVlp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcWlp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays

      ! Initialize LPT solver
      initialize_lpt: block
         use random, only: random_uniform
         real(WP) :: dp, tankRadius, xStart, zStart,yStart, px, py, pz
         integer :: i,j,np,npPerAxis, tmp
         real(WP), dimension(3) :: diskPos
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Get number of particles
         call param_read('Number of particles',np)

         ! Get tank radius
         call param_read("Rotor Disk center", diskPos)
         call param_read('tank radius',tankRadius)
         npPerAxis = int(floor(tankRadius/(sqrt(2.0_WP)*cfg%min_meshsize)))*2
         xStart = -0.5_WP*npPerAxis*cfg%min_meshsize
         zStart = -0.5_WP*npPerAxis*cfg%min_meshsize
         yStart = diskPos(2)

         ! Root process initializes np particles randomly
         if (lp%cfg%amRoot) then
            call lp%resize(np)
            do i=1,np
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! assign position in the domain
               py = int(floor(real(i,WP)/(npPerAxis**2)))*cfg%min_meshsize+yStart
               tmp = mod(i,npPerAxis**2)
               px = int(floor(real(tmp,WP)/npPerAxis))*cfg%min_meshsize+xStart
               pz = mod(tmp,npPerAxis)*cfg%min_meshsize+zStart
               lp%p(i)%pos= [px,py,pz]
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               lp%p(i)%angVel=0.0_WP
               ! Zero out collision forces
               lp%p(i)%Acol=0.0_WP
               lp%p(i)%Tcol=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n)
         call param_read('Wall restitution',lp%e_w,default=lp%e_n)
         call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         if (lp%cfg%amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Number of particles', np
         end if
      end block initialize_lpt

      ! Create a low Mach flow solver with bconds
      create_flow_solver: block
         real(WP) :: rho0
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Assign density
         call param_read('Density',rho0) ; rho = rho0
         fs%rho=rho*(1.0_WP-lp%VF)
         ! Configure pressure solver
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)

         ! Initialize the flow field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! Form momentum
         call fs%rho_multiply()
         call fs%interp_vel(Ui,Vi,Wi)
         resRho = 0.0_WP
         call fs%get_div(drhodt=resRHO)

      end block create_flow_solver

      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs

      ! Create a rotor disk
      create_and_initialize_rotorDisk: block
         use param, only: param_read

         integer :: Nr, Na          ! Number of radial and azimuthal points

         call param_read('Number of radical points',Nr)
         call param_read('Number of azimuthal points',Na)

         ! create blade object
         bl=blade(Nr=Nr,Na=Na)
         ! read blade parameters
         call param_read('Blade radius',bl%radius)
         call param_read('Blade twists',bl%twist)
         call param_read('Blade Chords',bl%chord)
         call param_read('Blade AoA',bl%aoa)
         call param_read('Blade Cd',bl%cd)
         call param_read('Blade Cl',bl%cl)

         rd = rotorDisk(bl=bl, cfg=cfg)
         call param_read('Rotor Disk number of blades',rd%nblades)
         call param_read('Rotor Disk min radius',rd%minR)
         call param_read('Rotor Disk max radius',rd%maxR)
         call param_read('Rotor Disk center',rd%center)
         call param_read('Rotor Disk axis',rd%axis)
         call param_read('Rotor Disk reference direction',rd%ref_dir)
         call param_read('Rotor Disk angular velocity',rd%omega)

         ! prepare the rotor disk
         call rd%prepareRotorDisk()

         call rd%calculateForce(fs%rho,Ui,Vi,Wi)    ! Get volumetric force

      end block create_and_initialize_rotorDisk

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         pmesh=partmesh(nvar=0,nvec=0,name='lpt')
         call lp%update_partmesh(pmesh)
      end block create_pmesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='StirredTank')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('Gib',cfg%Gib)
         call ens_out%add_scalar('VF',cfg%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('diskArea',rd%area)
         call ens_out%add_vector('diskForce',rd%forceX,rd%forceY,rd%forceZ)
         call ens_out%add_scalar('Density',fs%rho)
         call ens_out%add_scalar("viscosity",fs%visc)
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp_cfl,'Particle CFL')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
      end block create_monitor

   end subroutine simulation_init

   subroutine simulation_run
      use sgsmodel_class, only: constant_smag
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call lp%get_cfl(time%dt,cflc=lp_cfl,cfl=lp_cfl)
         call fs%get_cfl(time%dt,cfl)
         time%cfl=max(lp_cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old density, velocity, and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Reset here fluid properties
         fs%visc=visc

         ! Turbulence modeling
         call fs%get_strainrate(SR=SR)
         call sgs%get_visc(type=constant_smag,dt=time%dtold,rho=fs%rho,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         where (sgs%visc.lt.-fs%visc)
            sgs%visc=-fs%visc
         end where
         fs%visc=fs%visc+sgs%visc

         call fs%get_div_stress(resU,resV,resW)

         ! =================== Particle Solver ===================
         call lp%collide(dt=time%dtmid)
         call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
            srcU=srcUlp,srcV=srcVlp,srcW=srcWlp)

         ! Update density based on particle volume fraction
         fs%rho=rho*(1.0_WP-lp%VF)
         resRho=(fs%RHO-fs%RHOold)/time%dtmid

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Update density based on particle volume fraction and multi-vd
            fs%rho=0.5_WP*(fs%rho+fs%rhoold)

            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)

            call rd%calculateForce(fs%rho,Ui,Vi,Wi)    ! Get volumetric force
            resU=resU+rd%forceX
            resV=resV+rd%forceY
            resW=resW+rd%forceZ

            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            ! Add momentum source term from lpt
            add_lpt_src: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                        resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                     end do
                  end do
               end do
            end block add_lpt_src

            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            call fs%rho_multiply()

            ! Apply IB forcing to enforce BC at the walls
            ibforcing: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%U(i,j,k)=fs%U(i,j,k)*sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                        fs%V(i,j,k)=fs%V(i,j,k)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                        fs%W(i,j,k)=fs%W(i,j,k)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
               call fs%rho_multiply()
            end block ibforcing

         
            ! Solve Poisson equation
            call fs%get_div(drhodt=resRho)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            call cfg%integrate(A=fs%psolv%rhs,integral=int_RP)
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide()
            ! ===================================================

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(drhodt=resRho)

         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call lp%get_max()
         call lptfile%write()

      end do

   end subroutine simulation_run

   subroutine simulation_final

      implicit none

      ! Deallocate work arrays
      deallocate_work_arrays: block
         deallocate(resU,resV,resW,resRho)
         deallocate(Ui,Vi,Wi)
         deallocate(SR)
         deallocate(srcUlp,srcVlp,srcWlp)
      end block deallocate_work_arrays

   end subroutine simulation_final

end module simulation
