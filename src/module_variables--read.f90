!------------------------------------------------------------------------------
! MODULE: variables
!------------------------------------------------------------------------------
!> \brief Post-treatment variables.
!!
!! This module contains all the _additional_ variables that need to be
!! defined during the post-treatment.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
module variables

  use parameters
  use input
  use parallel


  implicit none


  integer (ip) , parameter , private :: nstatmax = 10000 !< maximum number of statistical files


!> \brief Simulation derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type sim_type


     real (dp) , dimension (nstatmax)                     :: dtime_per_ite   !< non-dimensional time step of the simulation
     real (dp)                                            :: sumdtime        !< non-dimensional time of the simulation
     integer (ip)                                         :: nfiles          !< corrected total number of stat files to post-process
     integer (ip)                                         :: sfile , efile   !< start/end file for each MPI process (for subroutine transform)


  end type sim_type


!> \brief Instantaneous derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type inst_type


     !> special variable to reduce iterations
     real (dp) , dimension (:,:,:) , allocatable     :: T0

     !> vector of conserved variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: u

     !> timestep variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: dt

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS
!     real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: H


  end type inst_type


!> \brief Reynolds average derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type reyavg_type


     !> basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: rho
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS

     real (dp) , dimension (:,:,:,:,:) , allocatable :: tau_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: nu

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H

     !> other basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs

     !> custom variables
     real (dp) , dimension (:,:,:,:,:) , allocatable :: tau


  end type reyavg_type


!> \brief Favre average derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type favavg_type


     !> basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: tke_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: nu

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H

     !> other basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs

     !> custom variables


  end type favavg_type


!> \brief Reynolds fluctuation derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type reyfluc_type


     !> basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: rho
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et


     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS ! May be not need ?! (A.Techer)
     real (dp) , dimension (:,:,:) , allocatable     :: nu

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H

     !> other basic variables

     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs

     !> custom variables
     real (dp) , dimension (:,:,:,:,:) , allocatable :: tau



  end type reyfluc_type


!> \brief Favre fluctuation derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type favfluc_type


     !> basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS ! May be not need ?! (A.Techer)
     real (dp) , dimension (:,:,:) , allocatable     :: nu

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H

     !> other basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs

     !> custom variables


  end type favfluc_type


!> \brief Statistical derived type declaration.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  type stat_type


     !> _repeated_ reynolds averages for convergence plots
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_P
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_ux
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_vy
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_wz

     !> other averages for convergence plots
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_uxux
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_vyvy
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_wzwz
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_uxvy
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_uxwz
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_vywz
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_enstrophy


     !> variances (typical)
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_ux
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_vy
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_wz
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_uxvy
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_uxwz
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_vywz
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_uxvywz

     real (dp) , dimension (:,:,:) , allocatable     :: favvar_ux
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_vy
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_wz
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_uxvy
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_uxwz
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_vywz
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_uxvywz

     !> variances (other)
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_p
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_acu
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_ent
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_T
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_T_acu
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_T_ent
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_du1_dx1
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_p_du1_dx1
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_p


     !> energy transport equation
     real (dp) , dimension (:,:,:) , allocatable     :: tke_dissip
     real (dp) , dimension (:,:,:) , allocatable     :: tke_press_strain
     real (dp) , dimension (:,:,:,:) , allocatable   :: tke_transp
     real (dp) , dimension (:,:,:) , allocatable     :: tke_ra_ff_ux
     real (dp) , dimension (:,:,:) , allocatable     :: tke_ra_ff_vy
     real (dp) , dimension (:,:,:) , allocatable     :: tke_ra_ff_wz

     !> vorticity transport equation
     real (dp) , dimension (:,:,:) , allocatable     :: vor_conv
     real (dp) , dimension (:,:,:) , allocatable     :: ens_conv
     real (dp) , dimension (:,:,:) , allocatable     :: vor_stret
     real (dp) , dimension (:,:,:) , allocatable     :: ens_stret
     real (dp) , dimension (:,:,:) , allocatable     :: vor_dila
     real (dp) , dimension (:,:,:) , allocatable     :: ens_dila
     real (dp) , dimension (:,:,:) , allocatable     :: vor_baro
     real (dp) , dimension (:,:,:) , allocatable     :: ens_baro
     real (dp) , dimension (:,:,:) , allocatable     :: vor_visc
     real (dp) , dimension (:,:,:) , allocatable     :: ens_visc


     !> spatial correlations
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_ux
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_vy
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_wz
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_p

     !> Reynolds transport equation
     real (dp) , dimension (:,:,:) , allocatable     :: rey_dissip_11 , rey_dissip_22 , rey_dissip_33 , &
                                                        rey_dissip_12 , rey_dissip_13 , rey_dissip_23
     real (dp) , dimension (:,:,:) , allocatable     :: rey_press_strain_11 , rey_press_strain_22 , rey_press_strain_33 , &
                                                        rey_press_strain_12 , rey_press_strain_13 , rey_press_strain_23
     real (dp) , dimension (:,:,:,:) , allocatable   :: rey_transp_11 , rey_transp_22 , rey_transp_33 , &
                                                        rey_transp_12 , rey_transp_13 , rey_transp_23



     !> Fourier spectra (only the amplitude)
     real (dp) , dimension (:,:) , allocatable       :: spec_ux
     real (dp) , dimension (:,:) , allocatable       :: spec_vy
     real (dp) , dimension (:,:) , allocatable       :: spec_wz


  end type stat_type


contains


!> \brief Allocate instantaneous derived type.
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
  subroutine alloc_inst_type (inst)


    type (inst_type) , intent (inout)                            :: inst !< instantaneous derived type


    integer (ip) :: ok


    allocate ( inst % T0    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               inst % u     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv )            , &

               inst % dt    ( sx:ex , sy:ey , sz:ez , 1:5 )                               , & ! A.Techer
               
               inst % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               inst % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate inst type')


    if ( LES ) then

       allocate ( inst % mu_SGS( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! A.Techer
!                  inst % tau_iso_SGS( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) ,&
                  stat = ok )

       if ( ok > 0 ) call abort_mpi ('error allocate inst type 2')

    end if ! LES


    inst % T0 = 0.0_dp

    inst % u  = 0.0_dp  
    inst % T  = 0.0_dp
    inst % H  = 0.0_dp


  end subroutine alloc_inst_type


!> \brief Write instantaneous derived type declaration into a file. (Subroutine not used... ?!)
!!
!> \author
!! Institut Pprime\n
!! CNRS - Poitiers University - ENSMA\n
!! FRANCE, Chasseneuil du Poitou 86360
!  subroutine write_inst_type (inst)


!    type (inst_type) , intent (in)                               :: inst !< instantaneous derived type


!    integer (ip)                 :: ok , i , j , k , l , m
!    character (len_default)      :: number_ascii
!
!
!    write ( number_ascii , format_restart ) rank
!
!
!    open ( unit_stat , file   = trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) , &
!                       form   = 'unformatted'   ,                                                      &
!                       status = 'unknown'       ,                                                      &
!                       iostat = ok              )
!
!    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) )
!
!
!    write (unit_stat) (((( inst % u (i,j,k,l) , i = sx , ex ) , &
!                                                j = sy , ey ) , &
!                                                k = sz , ez ) , &
!                                                l = 1  , nv )
!
!    write (unit_stat) (((( inst % Ydot (i,j,k,l) , i = sx , ex  ) , &
!                                                   j = sy , ey  ) , &
!                                                   k = sz , ez  ) , &
!                                                   l = 1  , nrv )
!
!    write (unit_stat) ((( inst % omega (i,j,k) , i = sx , ex ) , &
!                                                 j = sy , ey ) , &
!                                                 k = sz , ez )
!
!    write (unit_stat) ((( inst % mu (i,j,k) , i = sx , ex ) , &
!                                              j = sy , ey ) , &
!                                              k = sz , ez )
!
!
!    write (unit_stat) ((( inst % kpa (i,j,k) , i = sx , ex ) , &
!                                               j = sy , ey ) , &
!                                               k = sz , ez )
!
!    write (unit_stat) ((( inst % ct (i,j,k) , i = sx , ex ) , &
!                                              j = sy , ey ) , &
!                                              k = sz , ez )
!
!    write (unit_stat) (((( inst % dm (i,j,k,l) , i = sx , ex  ) , &
!                                                 j = sy , ey  ) , &
!                                                 k = sz , ez  ) , &
!                                                 l = 1  , nrv )
!
!    write (unit_stat) (((( inst % tdr (i,j,k,l) , i = sx , ex  ) , &
!                                                  j = sy , ey  ) , &
!                                                  k = sz , ez  ) , &
!                                                  l = 1  , nrv )
!
!    write (unit_stat) ((((( inst % rd (i,j,k,l,m) , i = sx , ex  ) , &
!                                                    j = sy , ey  ) , &
!                                                    k = sz , ez  ) , &
!                                                    l = 1  , nrv ) , &
!                                                    m = 1  , nrv )
!
!    write (unit_stat) ((( inst % T (i,j,k) , i = sx , ex ) , &
!                                             j = sy , ey ) , &
!                                             k = sz , ez )
!
!    write (unit_stat) ((( inst % H (i,j,k) , i = sx , ex ) , &
!                                             j = sy , ey ) , &
!                                             k = sz , ez )
!
!    write (unit_stat) ((( inst % cp (i,j,k) , i = sx , ex ) , &
!                                              j = sy , ey ) , &
!                                              k = sz , ez )
!
!    write (unit_stat) ((( inst % W_i (i,j,k) , i = sx , ex ) , &
!                                               j = sy , ey ) , &
!                                               k = sz , ez )
!
!    write (unit_stat) (((( inst % W_scal_i (i,j,k,l) , i = sx , ex    ) , &
!                                                       j = sy , ey    ) , &
!                                                       k = sz , ez    ) , &
!                                                       l = 1  , nreac )
!
!    close (unit_stat)
!
!
!  end subroutine write_inst_type


!> \brief Read instantaneous derived type declaration from a file.

  subroutine read_inst_type ( inp , inst )

    type (inp_type)  , intent (in)                               :: inp  !< input derived type
    type (inst_type) , intent (inout)                            :: inst !< instantaneous derived type


    integer (ip)                    :: ok , irank , i , j , k , l , m

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank        = 0


    if ( nproc == 1 ) then ! sequential problem


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                      &
                                   status = 'unknown'       ,                                                      &
                                   iostat = ok              )

                if ( ok /= 0 ) &
                call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) )


                read (unit_stat) (((( inst % u (i,j,k,l) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz ) , &
                                                           l = 1  , nv )

                
             end do
          end do
       end do


    else ! parallel problem


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                      &
                          status = 'unknown'       ,                                                      &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) )


       read (unit_stat) (((( inst % u (i,j,k,l) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez ) , &
                                                  l = 1  , nv )


    end if


    close (unit_stat)


  end subroutine read_inst_type


!> \brief Allocate Reynolds average derived type.

  subroutine alloc_reyavg_type (reyavg)


    type (reyavg_type) , intent (out)                               :: reyavg !< Reynolds average derived type


    integer (ip) :: ok


    allocate ( reyavg % rho   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % tau ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate reyavg type')

    if ( LES ) then

       allocate ( reyavg % mu_SGS ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! A.Techer
                  reyavg % tau_iso_SGS( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                  reyavg % tau_SGS ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax )  , &

                  stat = ok )

       if ( ok > 0 ) call abort_mpi ('error allocate reyavg type 3')

    end if


    ! reyavg % Y     = 0.0_dp
    ! reyavg % rho   = 0.0_dp
    ! reyavg % ux    = 0.0_dp
    ! reyavg % vy    = 0.0_dp
    ! reyavg % wz    = 0.0_dp
    ! reyavg % et    = 0.0_dp

    ! reyavg % mu    = 0.0_dp
    ! reyavg % nu    = 0.0_dp

    ! reyavg % T     = 0.0_dp
    ! reyavg % H     = 0.0_dp


    ! reyavg % P     = 0.0_dp
    ! reyavg % Tt    = 0.0_dp
    ! reyavg % Ht    = 0.0_dp
    ! reyavg % cs    = 0.0_dp


    ! reyavg % tau   = 0.0_dp



  end subroutine alloc_reyavg_type


!> \brief Write Reynolds average derived type into a file.

  subroutine write_reyavg_type (reyavg)


    type (reyavg_type) , intent (in)                                :: reyavg !< Reynolds average derived type


    integer (ip)                 :: ok , i , j , k , l , m
    character (len_default)      :: number_ascii


    write ( number_ascii , format_restart ) rank


    open ( unit_stat , file   = trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) , &
                       form   = 'unformatted'   ,                                                        &
                       status = 'unknown'       ,                                                        &
                       iostat = ok              )

    if ( ok /= 0 ) then
       call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) )
    end if


    write (unit_stat) ((( reyavg % rho (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

    write (unit_stat) ((( reyavg % ux (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % vy (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % wz (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % et (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % mu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % nu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % T (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( reyavg % H (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( reyavg % P (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % Tt (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % Ht (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % M (i,j,k) , i = sx , ex ) , &
                                              j = sy , ey ) , &
                                              k = sz , ez )

   write (unit_stat) ((( reyavg % cs (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((((( reyavg % tau (i,j,k,l,m) , i = sx , ex )     , &
                                                      j = sy , ey )     , &
                                                      k = sz , ez )     , &
                                                      l = 1 , ndimmax ) , &
                                                      m = 1 , ndimmax )


   close (unit_stat)


  end subroutine write_reyavg_type


!> \brief Read Reynolds average derived type from a file.

  subroutine read_reyavg_type ( inp , reyavg )


    type (inp_type) , intent (in)                                 :: inp    !< input derived type
    type (reyavg_type) , intent (inout)                           :: reyavg !< Reynolds average derived type


    integer (ip)                    :: ok , irank , i , j , k , l , m

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank      = 0


    if ( nproc == 1 ) then ! sequential problem


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                        &
                                   status = 'unknown'       ,                                                        &
                                   iostat = ok              )

                if ( ok /= 0 ) &
                call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) )


                read (unit_stat) ((( reyavg % rho (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( reyavg % ux (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % vy (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % wz (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % et (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % mu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % nu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % T (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % H (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % P (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % Tt (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % Ht (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % M (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % cs (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((((( reyavg % tau (i,j,k,l,m) , i = ix , fx )     , &
                                                                  j = iy , fy )     , &
                                                                  k = iz , fz )     , &
                                                                  l = 1 , ndimmax ) , &
                                                                  m = 1 , ndimmax )


             end do
          end do
       end do


    else ! parallel problem


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                        &
                          status = 'unknown'       ,                                                        &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) )


       read (unit_stat) ((( reyavg % rho (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( reyavg % ux (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % vy (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % wz (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % et (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % mu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % nu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % T (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % H (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % P (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % Tt (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % Ht (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % M (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % cs (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((((( reyavg % tau (i,j,k,l,m) , i = sx , ex )     , &
                                                         j = sy , ey )     , &
                                                         k = sz , ez )     , &
                                                         l = 1 , ndimmax ) , &
                                                         m = 1 , ndimmax )


    end if


    close (unit_stat)


  end subroutine read_reyavg_type


!> \brief Allocate Favre average derived type.

  subroutine alloc_favavg_type (favavg)


    type (favavg_type) , intent (out)                               :: favavg !< Favre average derived type


    integer (ip) :: ok


    allocate ( favavg % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               
               favavg % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate favavg type')

    if ( LES ) then

       allocate ( favavg % mu_SGS   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! A.Techer
                  favavg % tke_SGS  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

                  stat = ok )

       if ( ok > 0 ) call abort_mpi ('error allocate favavg type 3')

    end if


    ! favavg % ux    = 0.0_dp
    ! favavg % vy    = 0.0_dp
    ! favavg % wz    = 0.0_dp
    ! favavg % et    = 0.0_dp

    ! favavg % mu    = 0.0_dp
    ! favavg % nu    = 0.0_dp

    ! favavg % T     = 0.0_dp
    ! favavg % H     = 0.0_dp

    ! favavg % P     = 0.0_dp
    ! favavg % Tt    = 0.0_dp
    ! favavg % Ht    = 0.0_dp
    ! favavg % cs    = 0.0_dp


  end subroutine alloc_favavg_type


!> \brief Write Favre average derived type to a file.

  subroutine write_favavg_type (favavg)


    type (favavg_type) , intent (in)                                 :: favavg !< Favre average derived type


    integer (ip)                 :: ok , i , j , k , l
    character (len_default)      :: number_ascii


    write ( number_ascii , format_restart ) rank


    open ( unit_stat , file   = trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) , &
                       form   = 'unformatted'   ,                                                        &
                       status = 'unknown'       ,                                                        &
                       iostat = ok              )

    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) )


    write (unit_stat) ((( favavg % ux (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % vy (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % wz (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % et (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % mu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % nu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % T (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( favavg % H (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( favavg % P (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % Tt (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % Ht (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % M (i,j,k) , i = sx , ex ) , &
                                              j = sy , ey ) , &
                                              k = sz , ez )

   write (unit_stat) ((( favavg % cs (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )


   close (unit_stat)


  end subroutine write_favavg_type


!> \brief Read Favre average derived from a file.

  subroutine read_favavg_type ( inp , favavg )


    type (inp_type) , intent (in)                                 :: inp    !< input derived type
    type (favavg_type) , intent (inout)                           :: favavg !< Favre average derived type


    integer (ip)                    :: ok , irank , i , j , k , l

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank      = 0


    if ( nproc == 1 ) then


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                        &
                                   status = 'unknown'       ,                                                        &
                                   iostat = ok              )

                if ( ok /= 0 ) &
                 call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) )


                read (unit_stat) ((( favavg % ux (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % vy (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % wz (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % et (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % mu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % nu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % T (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % H (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % P (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % Tt (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % Ht (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % M (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % cs (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

             end do
          end do
       end do


    else


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                        &
                          status = 'unknown'       ,                                                        &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) )


       read (unit_stat) ((( favavg % ux (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % vy (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % wz (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % et (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % mu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % nu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % T (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % H (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % P (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % Tt (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % Ht (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % M (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % cs (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

    end if


    close (unit_stat)


  end subroutine read_favavg_type


!> \brief Allocate Reynolds fluctuation derived type.

  subroutine alloc_reyfluc_type (reyfluc)


    type (reyfluc_type) , intent (out)                               :: reyfluc !< Reynolds fluctuation derived type


    integer (ip) :: ok


    allocate ( reyfluc % rho   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % tau ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate reyfluc type')


    ! reyfluc % rho   = 0.0_dp
    ! reyfluc % ux    = 0.0_dp
    ! reyfluc % vy    = 0.0_dp
    ! reyfluc % wz    = 0.0_dp
    ! reyfluc % et    = 0.0_dp

    ! reyfluc % mu    = 0.0_dp
    ! reyfluc % nu    = 0.0_dp

    ! reyfluc % T     = 0.0_dp
    ! reyfluc % H     = 0.0_dp

    ! reyfluc % P     = 0.0_dp
    ! reyfluc % Tt    = 0.0_dp
    ! reyfluc % Ht    = 0.0_dp
    ! reyfluc % cs    = 0.0_dp

    ! reyfluc % tau   = 0.0_dp


  end subroutine alloc_reyfluc_type


!> \brief Allocate Favre fluctuation derived type.

  subroutine alloc_favfluc_type (favfluc)


    type (favfluc_type) , intent (out)                               :: favfluc !< Favre fluctuation derived type


    integer (ip) :: ok


    allocate ( favfluc % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate favfluc type')

    ! favfluc % ux    = 0.0_dp
    ! favfluc % vy    = 0.0_dp
    ! favfluc % wz    = 0.0_dp
    ! favfluc % et    = 0.0_dp

    ! favfluc % mu    = 0.0_dp
    ! favfluc % nu    = 0.0_dp

    ! favfluc % T     = 0.0_dp
    ! favfluc % H     = 0.0_dp

    ! favfluc % P     = 0.0_dp
    ! favfluc % Tt    = 0.0_dp
    ! favfluc % Ht    = 0.0_dp
    ! favfluc % cs    = 0.0_dp


  end subroutine alloc_favfluc_type


!> \brief Allocate statistical derived type.

  subroutine alloc_stat_type (stat)


    type (stat_type) , intent (out)                                   :: stat !< statistical derived type


    integer (ip) :: ok


    allocate ( stat % reyavg_P      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % reyavg_uxux      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_vyvy      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_wzwz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_uxvy      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_uxwz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_vywz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_enstrophy ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % reyvar_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_uxvy   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_uxwz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_vywz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_uxvywz ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % favvar_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_uxvy   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_uxwz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_vywz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_uxvywz ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % reyvar_p         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_acu   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_ent   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_T         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_T_acu     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_T_ent     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favvar_du1_dx1   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_p_du1_dx1 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_p     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % tke_dissip       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_press_strain ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_transp       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )        , &
               stat % tke_ra_ff_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_ra_ff_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_ra_ff_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % vor_conv         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_conv         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_stret        ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_stret        ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_dila         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_dila         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_baro         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_baro         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_visc         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_visc         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % tmp_corr_ux      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_vy      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_wz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_p       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % rey_dissip_11    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_22    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_33    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_12    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_13    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_23    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % rey_press_strain_11 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_22 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_33 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_12 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_13 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_23 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &

               stat % rey_transp_11 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_22 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_33 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_12 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_13 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_23 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &

               stat % spec_ux  ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_vy  ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_wz  ( sz-ng:ez+ng , 1:ntimes )                                           , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate stat type')


    ! stat % reyavg_P         = 0.0_dp
    ! stat % reyavg_ux        = 0.0_dp
    ! stat % reyavg_vy        = 0.0_dp
    ! stat % reyavg_wz        = 0.0_dp

    ! stat % reyavg_uxux      = 0.0_dp
    ! stat % reyavg_vyvy      = 0.0_dp
    ! stat % reyavg_wzwz      = 0.0_dp
    ! stat % reyavg_uxvy      = 0.0_dp
    ! stat % reyavg_uxwz      = 0.0_dp
    ! stat % reyavg_vywz      = 0.0_dp
    ! stat % reyavg_enstrophy = 0.0_dp

    ! stat % reyvar_ux        = 0.0_dp
    ! stat % reyvar_vy        = 0.0_dp
    ! stat % reyvar_wz        = 0.0_dp
    ! stat % reyvar_uxvy      = 0.0_dp
    ! stat % reyvar_uxwz      = 0.0_dp
    ! stat % reyvar_vywz      = 0.0_dp
    ! stat % reyvar_uxvywz    = 0.0_dp

    ! stat % favvar_ux        = 0.0_dp
    ! stat % favvar_vy        = 0.0_dp
    ! stat % favvar_wz        = 0.0_dp
    ! stat % favvar_uxvy      = 0.0_dp
    ! stat % favvar_uxwz      = 0.0_dp
    ! stat % favvar_vywz      = 0.0_dp
    ! stat % favvar_uxvywz    = 0.0_dp

    ! stat % reyvar_p         = 0.0_dp
    ! stat % reyvar_rho       = 0.0_dp
    ! stat % reyvar_rho_acu   = 0.0_dp
    ! stat % reyvar_rho_ent   = 0.0_dp
    ! stat % reyvar_T         = 0.0_dp
    ! stat % reyvar_T_acu     = 0.0_dp
    ! stat % reyvar_T_ent     = 0.0_dp
    ! stat % favvar_du1_dx1   = 0.0_dp
    ! stat % reyvar_p_du1_dx1 = 0.0_dp
    ! stat % reyvar_rho_p     = 0.0_dp

    ! stat % tke_dissip       = 0.0_dp
    ! stat % tke_press_strain = 0.0_dp
    ! stat % tke_transp       = 0.0_dp
    ! stat % tke_ra_ff_ux     = 0.0_dp
    ! stat % tke_ra_ff_vy     = 0.0_dp
    ! stat % tke_ra_ff_wz     = 0.0_dp

    ! stat % vor_conv         = 0.0_dp
    ! stat % ens_conv         = 0.0_dp
    ! stat % vor_stret        = 0.0_dp
    ! stat % ens_stret        = 0.0_dp
    ! stat % vor_dila         = 0.0_dp
    ! stat % ens_dila         = 0.0_dp
    ! stat % vor_baro         = 0.0_dp
    ! stat % ens_baro         = 0.0_dp
    ! stat % vor_visc         = 0.0_dp
    ! stat % ens_visc         = 0.0_dp

    ! stat % tmp_corr_ux      = 0.0_dp
    ! stat % tmp_corr_vy      = 0.0_dp
    ! stat % tmp_corr_wz      = 0.0_dp
    ! stat % tmp_corr_p       = 0.0_dp

    ! stat % rey_dissip_11 = 0.0_dp       ; stat % rey_dissip_22 = 0.0_dp       ; stat % rey_dissip_33 = 0.0_dp
    ! stat % rey_dissip_12 = 0.0_dp       ; stat % rey_dissip_13 = 0.0_dp       ; stat % rey_dissip_23 = 0.0_dp
    ! stat % rey_press_strain_11 = 0.0_dp ; stat % rey_press_strain_22 = 0.0_dp ; stat % rey_press_strain_33 = 0.0_dp
    ! stat % rey_press_strain_12 = 0.0_dp ; stat % rey_press_strain_13 = 0.0_dp ; stat % rey_press_strain_23 = 0.0_dp
    ! stat % rey_transp_11 = 0.0_dp       ; stat % rey_transp_22 = 0.0_dp       ; stat % rey_transp_33 = 0.0_dp
    ! stat % rey_transp_12 = 0.0_dp       ; stat % rey_transp_13 = 0.0_dp       ; stat % rey_transp_23 = 0.0_dp

    ! stat % spec_ux          = 0.0_dp
    ! stat % spec_vy          = 0.0_dp
    ! stat % spec_wz          = 0.0_dp

  end subroutine alloc_stat_type


!> \brief Write statistical derived type.

  subroutine write_stat_type (stat)


    type (stat_type) , intent (in)                                :: stat !< statistical derived type


    integer (ip)                 :: ok , i , j , k , l
    character (len_default)      :: number_ascii


    write ( number_ascii , format_restart ) rank


    open ( unit_stat , file   = trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) , &
                       form   = 'unformatted'   ,                                                      &
                       status = 'unknown'       ,                                                      &
                       iostat = ok              )

    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) )


    write (unit_stat) ((( stat % reyavg_P (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyavg_ux (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_vy (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_wz (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_uxux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_vyvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_wzwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_uxvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_uxwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_vywz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_enstrophy (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % reyvar_ux (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_vy (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_wz (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyvar_vywz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % favvar_ux (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_vy (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_wz (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favvar_vywz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % reyvar_p (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_acu (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_ent (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

    write (unit_stat) ((( stat % reyvar_T (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyvar_T_acu (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % reyvar_T_ent (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % favvar_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

    write (unit_stat) ((( stat % reyvar_p_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_p (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % tke_dissip (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % tke_press_strain (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) (((( stat % tke_transp (i,j,k,l) , i = sx , ex      ) , &
                                                         j = sy , ey      ) , &
                                                         k = sz , ez      ) , &
                                                         l = 1  , ndimmax )

    write (unit_stat) ((( stat % tke_ra_ff_ux (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % tke_ra_ff_vy (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % tke_ra_ff_wz (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % vor_conv (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_conv (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % vor_stret (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % ens_stret (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % vor_dila (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_dila (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % vor_baro (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_baro (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % vor_visc (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_visc (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_p (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_11 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_22 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_33 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_12 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_13 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_23 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_11 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_22 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_33 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_12 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_13 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_23 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) (((( stat % rey_transp_11 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_22 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_33 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_12 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_13 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_23 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    ! pdfs and spectras do not need to be saved (not plotted in paraview)

    write (unit_stat) (( stat % spec_ux (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )

    write (unit_stat) (( stat % spec_vy (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )

    write (unit_stat) (( stat % spec_wz (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )


    close (unit_stat)


  end subroutine write_stat_type


!> \brief Read statistical derived type.

  subroutine read_stat_type ( inp , stat )


    type (inp_type) , intent (in)                                 :: inp  !< input derived type
    type (stat_type) , intent (inout)                             :: stat !< statistical derived type


    integer (ip)                    :: ok , irank , i , j , k , l

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank        = 0


    if ( nproc == 1 ) then


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                      &
                                   status = 'unknown'       ,                                                      &
                                   iostat = ok              )

                if ( ok /= 0 ) &
                call abort_mpi &
                ('error opening ' // trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) )


                read (unit_stat) ((( stat % reyavg_P (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyavg_ux (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_vy (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_wz (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_uxux (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_vyvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_wzwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_uxvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_uxwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_vywz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_enstrophy (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % reyvar_ux (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_vy (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_wz (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_uxvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyvar_uxwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyvar_vywz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyvar_uxvywz (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % favvar_ux (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_vy (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_wz (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_uxvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favvar_uxwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favvar_vywz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favvar_uxvywz (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % reyvar_p (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_acu (i,j,k) , i = ix , fx ) , &
                                                                     j = iy , fy ) , &
                                                                     k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_ent (i,j,k) , i = ix , fx ) , &
                                                                     j = iy , fy ) , &
                                                                     k = iz , fz )

                read (unit_stat) ((( stat % reyvar_T (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyvar_T_acu (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % reyvar_T_ent (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % favvar_du1_dx1 (i,j,k) , i = ix , fx ) , &
                                                                     j = iy , fy ) , &
                                                                     k = iz , fz )

                read (unit_stat) ((( stat % reyvar_p_du1_dx1 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_p (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )               

                read (unit_stat) ((( stat % tke_dissip (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % tke_press_strain (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) (((( stat % tke_transp (i,j,k,l) , i = ix , fx      ) , &
                                                                    j = iy , fy      ) , &
                                                                    k = iz , fz      ) , &
                                                                    l = 1  , ndimmax )

                read (unit_stat) ((( stat % tke_ra_ff_ux (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % tke_ra_ff_vy (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % tke_ra_ff_wz (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % vor_conv (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_conv (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % vor_stret (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % ens_stret (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % vor_dila (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_dila (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % vor_baro (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_baro (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % vor_visc (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_visc (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_ux (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_vy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_wz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_p (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_11 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_22 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_33 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_12 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_13 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_23 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_11 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_22 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_33 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_12 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_13 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_23 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) (((( stat % rey_transp_11 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_22 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_33 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_12 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_13 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_23 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (( stat % spec_ux (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )

                read (unit_stat) (( stat % spec_vy (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )

                read (unit_stat) (( stat % spec_wz (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )


             end do
          end do
       end do


    else ! parallel problem


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                      &
                          status = 'unknown'       ,                                                      &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) )


       read (unit_stat) ((( stat % reyavg_P (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyavg_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_uxux (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_vyvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_wzwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_uxvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_uxwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_vywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_enstrophy (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % reyvar_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyvar_vywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % favvar_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favvar_vywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % reyvar_p (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_acu (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_ent (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

       read (unit_stat) ((( stat % reyvar_T (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyvar_T_acu (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % reyvar_T_ent (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % favvar_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

       read (unit_stat) ((( stat % reyvar_p_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_p (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )


       read (unit_stat) ((( stat % tke_dissip (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % tke_press_strain (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) (((( stat % tke_transp (i,j,k,l) , i = sx , ex      ) , &
                                                           j = sy , ey      ) , &
                                                           k = sz , ez      ) , &
                                                           l = 1  , ndimmax )

       read (unit_stat) ((( stat % tke_ra_ff_ux (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % tke_ra_ff_vy (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % tke_ra_ff_wz (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % vor_conv (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_conv (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_stret (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % ens_stret (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_dila (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_dila (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_baro (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_baro (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_visc (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_visc (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_ux (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_vy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_wz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_p (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_11 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_22 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_33 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_12 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_13 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_23 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_11 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_22 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_33 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_12 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_13 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_23 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) (((( stat % rey_transp_11 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_22 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_33 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_12 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_13 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_23 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (( stat % spec_ux (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )

       read (unit_stat) (( stat % spec_vy (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )

       read (unit_stat) (( stat % spec_wz (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )


    end if


    close (unit_stat)


  end subroutine read_stat_type


end module variables
